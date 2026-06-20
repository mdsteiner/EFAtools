# Analytic expected-information standard errors for the unrotated ML solution.

# Independent reference: the unit expected information assembled densely as
# 1/2 tr(W dSigma_i W dSigma_j), the identification constraint off-diag(Lambda' Psi^-1
# Lambda) = 0 differentiated by central finite differences, bordered, inverted and scaled
# by 1 / (N - 1). This shares no code with the package's closed-form assembly.
.ref_se_information <- function(L, psi, N) {
  p <- nrow(L); k <- ncol(L); pk <- p * k; q <- pk + p
  Sigma <- tcrossprod(L); diag(Sigma) <- diag(Sigma) + psi
  W <- solve(Sigma)

  dS <- vector("list", q); idx <- 0
  for (b in seq_len(k)) for (a in seq_len(p)) {
    idx <- idx + 1; M <- matrix(0, p, p); M[a, ] <- L[, b]; M[, a] <- L[, b]; M[a, a] <- 2 * L[a, b]
    dS[[idx]] <- M
  }
  for (c in seq_len(p)) { idx <- idx + 1; M <- matrix(0, p, p); M[c, c] <- 1; dS[[idx]] <- M }
  WdS <- lapply(dS, function(M) W %*% M)
  A <- matrix(0, q, q)
  for (i in seq_len(q)) for (j in i:q) {
    v <- 0.5 * sum(WdS[[i]] * t(WdS[[j]])); A[i, j] <- v; A[j, i] <- v
  }

  g <- function(th) {
    Lm <- matrix(th[seq_len(pk)], p, k); ps <- th[pk + seq_len(p)]
    M <- crossprod(Lm, Lm / ps)
    M[upper.tri(M)]
  }
  th0 <- c(as.vector(L), psi); nc <- k * (k - 1L) / 2L
  Cmat <- matrix(0, nc, q); h <- 1e-6
  for (j in seq_len(q)) {
    tp <- th0; tp[j] <- tp[j] + h; tm <- th0; tm[j] <- tm[j] - h
    Cmat[, j] <- (g(tp) - g(tm)) / (2 * h)
  }

  Aug <- rbind(cbind(A, t(Cmat)), cbind(Cmat, matrix(0, nc, nc)))
  V <- solve(Aug)[seq_len(q), seq_len(q)] / (N - 1)
  se <- sqrt(diag(V))
  list(loadings_se = matrix(se[seq_len(pk)], p, k), uniquenesses_se = se[pk + seq_len(p)])
}


test_that("closed-form information SEs match an independent dense reference", {
  L <- matrix(c(0.70, 0.60, 0.50, 0.10, 0.00, 0.20,
                0.10, 0.05, 0.20, 0.70, 0.60, 0.50), 6, 2)
  psi <- 1 - rowSums(L^2)

  ours <- EFAtools:::.se_information_ml(L, psi, N = 250)
  ref <- .ref_se_information(L, psi, N = 250)

  expect_equal(ours$loadings_se, ref$loadings_se, tolerance = 1e-5)
  expect_equal(ours$uniquenesses_se, ref$uniquenesses_se, tolerance = 1e-5)
})


test_that("information uniqueness SEs match lavaan", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  R <- test_models$baseline$cormat
  N <- 500
  fit <- EFA(R, n_factors = 3, N = N, method = "ML", rotation = "none", se = "information")

  lav <- lavaan::efa(sample.cov = R, sample.nobs = N, nfactors = 3,
                     rotation = "none", estimator = "ML")
  lav_fit <- Filter(function(e) inherits(e, "lavaan"), lav)[[1]]
  pe <- lavaan::parameterEstimates(lav_fit)
  uq <- pe[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% rownames(R), ]
  uq <- uq[match(rownames(R), uq$lhs), ]

  # Uniquenesses are rotation-invariant, so they are directly comparable across the two
  # packages' identifications. The small gap is the N vs N - 1 scaling.
  expect_equal(unname(fit$boot.SE$uniquenesses), uq$se, tolerance = 0.02)
})


test_that("single-factor information loading and uniqueness SEs match lavaan", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  # With one factor there is no rotational indeterminacy, so the loadings are identified
  # (up to sign) the same way in both packages and their SEs are directly comparable.
  R <- stats::cor(GRiPS_raw)
  N <- nrow(GRiPS_raw)
  fit <- EFA(GRiPS_raw, n_factors = 1, method = "ML", rotation = "none", se = "information")

  lav <- lavaan::efa(sample.cov = R, sample.nobs = N, nfactors = 1,
                     rotation = "none", estimator = "ML")
  lav_fit <- Filter(function(e) inherits(e, "lavaan"), lav)[[1]]
  pe <- lavaan::parameterEstimates(lav_fit)
  ld <- pe[pe$op == "=~", ]; ld <- ld[match(colnames(R), ld$rhs), ]
  uq <- pe[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% colnames(R), ]
  uq <- uq[match(colnames(R), uq$lhs), ]

  expect_equal(as.vector(fit$boot.SE$unrot_loadings), ld$se, tolerance = 0.01)
  expect_equal(unname(fit$boot.SE$uniquenesses), uq$se, tolerance = 0.01)
})


test_that("information loading SEs track the bootstrap on an identifiable solution", {
  skip_on_cran()
  skip_if_not_slow()

  # Well-separated three-factor structure at a large N, so the unrotated solution is
  # cleanly identified and the bootstrap alignment is stable. The data are multivariate
  # normal and the model is correctly specified, the regime in which the expected-information
  # and the bootstrap SE are both consistent and coincide (Yuan & Hayashi, 2006); they are
  # different finite-sample estimators, so they agree in structure (strong correlation)
  # rather than to a tight per-element tolerance.
  set.seed(123)
  R0 <- test_models$baseline$cormat
  X <- matrix(rnorm(2500 * ncol(R0)), 2500) %*% chol(R0)
  colnames(X) <- colnames(R0)

  fa <- EFA(X, n_factors = 3, method = "ML", rotation = "none", se = "information")
  fb <- EFA(X, n_factors = 3, method = "ML", rotation = "none",
            se = "np-boot", b_boot = 250, seed = 1)

  expect_gt(cor(c(fa$boot.SE$unrot_loadings), c(fb$boot.SE$unrot_loadings)), 0.9)
  expect_lt(abs(median(fa$boot.SE$unrot_loadings / fb$boot.SE$unrot_loadings) - 1), 0.35)
})


test_that("information SEs populate the bootstrap SE/CI schema", {
  R <- test_models$baseline$cormat
  fit <- EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "none", se = "information")

  expect_setequal(names(fit$boot.SE), c("unrot_loadings", "uniquenesses"))
  expect_setequal(names(fit$boot.CI), c("unrot_loadings", "uniquenesses"))
  expect_null(fit$boot.arrays)

  expect_equal(dim(fit$boot.SE$unrot_loadings), dim(fit$unrot_loadings))
  expect_length(fit$boot.SE$uniquenesses, ncol(R))
  expect_false(anyNA(fit$boot.SE$unrot_loadings))

  # Wald intervals bracket the point estimates.
  expect_true(all(fit$boot.CI$unrot_loadings$lower <= fit$unrot_loadings))
  expect_true(all(fit$unrot_loadings <= fit$boot.CI$unrot_loadings$upper))
  expect_named(fit$boot.CI$uniquenesses, c("lower", "upper"))
})


test_that("information SEs cover only the unrotated solution under a (non-promax) rotation", {
  R <- test_models$baseline$cormat
  fit <- EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "oblimin", se = "information")

  expect_false(is.null(fit$rot_loadings))                 # the rotated solution is still returned
  expect_setequal(names(fit$boot.SE), c("unrot_loadings", "uniquenesses"))
  expect_null(fit$boot.SE$rot_loadings)                   # rotated SEs are not produced here

  # No unrotated-loading CI table is shown under a rotation, so neither a Wald CI table nor
  # its provenance note may appear (the note must not advertise hidden intervals).
  out <- testthat::capture_output(print(summary(fit)))
  expect_false(any(grepl("Wald", out, fixed = TRUE)))
})


test_that("information SEs are available from a correlation matrix and match the raw fit", {
  set.seed(5)
  R0 <- test_models$baseline$cormat
  X <- matrix(rnorm(800 * ncol(R0)), 800) %*% chol(R0)
  colnames(X) <- colnames(R0)

  raw <- EFA(X, n_factors = 3, method = "ML", rotation = "none", se = "information")
  cmat <- EFA(stats::cor(X), n_factors = 3, N = 800, method = "ML",
              rotation = "none", se = "information")

  expect_equal(raw$boot.SE$unrot_loadings, cmat$boot.SE$unrot_loadings, tolerance = 1e-5)
  expect_equal(raw$boot.SE$uniquenesses, cmat$boot.SE$uniquenesses, tolerance = 1e-5)
})


test_that("a Heywood case yields NA SEs with a classed warning", {
  # A communality above one drives a uniqueness below its zero boundary, where Psi^-1 and
  # the Lambda' Psi^-1 Lambda identification are no longer defined, so no analytic SE exists.
  L <- matrix(c(sqrt(1.11), 0.6, 0.5, 0.4), 4, 1)   # h2[1] = 1.11 > 1, so psi[1] < 0
  fit_out <- list(unrot_loadings = L)

  expect_warning(
    out <- EFAtools:::.se_information(fit_out, N = 200, ci = 0.95),
    class = "efa_se_unreliable"
  )
  expect_true(anyNA(out$SE$unrot_loadings))
})


test_that("unsupported se combinations abort early with a clear class", {
  R <- test_models$baseline$cormat

  expect_error(
    EFA(R, n_factors = 3, N = 500, method = "PAF", rotation = "none", se = "information"),
    class = "efa_se_unsupported"
  )
  expect_error(
    EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "promax", se = "information"),
    class = "efa_se_unsupported"
  )
  expect_error(
    EFA(R, n_factors = 3, N = 500, method = "PAF", rotation = "none", se = "sandwich"),
    class = "efa_se_unsupported"
  )
  expect_error(
    EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "none", se = "sandwich"),
    class = "efa_se_not_implemented"
  )
  # Sandwich is rejected up front for every non-PAF, non-promax estimator, not only ML.
  expect_error(
    EFA(R, n_factors = 3, N = 500, method = "ULS", rotation = "none", se = "sandwich"),
    class = "efa_se_not_implemented"
  )
  # A correlation matrix carries no sample size, which the analytic methods require.
  expect_error(
    EFA(R, n_factors = 3, method = "ML", rotation = "none", se = "information"),
    class = "efa_se_no_n"
  )
})


test_that("an analytic-SE fit prints and summarises without error", {
  R <- test_models$baseline$cormat
  fit <- EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "none", se = "information")

  expect_no_error(testthat::capture_output(print(fit)))
  expect_no_error(testthat::capture_output(print(summary(fit))))

  # The Wald provenance note (with the CI table) appears in the full summary, but not in the
  # brief print() view nor when intervals are suppressed with ci = "none".
  expect_true(any(grepl("Wald", testthat::capture_output(print(summary(fit))), fixed = TRUE)))
  expect_false(any(grepl("Wald", testthat::capture_output(print(fit)), fixed = TRUE)))
  expect_false(any(grepl("Wald", testthat::capture_output(print(summary(fit, ci = "none"))), fixed = TRUE)))
})
