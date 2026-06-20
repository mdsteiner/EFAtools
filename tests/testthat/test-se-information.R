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


# Match the columns of B to the columns of A by maximum absolute loading-column correlation,
# returning the permutation and signs that align B to A. Rotated solutions order and sign factors
# arbitrarily across packages and bootstrap replicates, so the rotated SEs must be aligned before
# they are compared.
.align_factor_cols <- function(A, B) {
  k <- ncol(A); perm <- integer(k); sgn <- numeric(k); used <- logical(k)
  for (j in seq_len(k)) {
    cors <- vapply(seq_len(k),
                   function(i) if (used[i]) NA_real_ else stats::cor(A[, j], B[, i]), 0)
    i <- which.max(abs(cors)); perm[j] <- i; sgn[j] <- sign(cors[i]); used[i] <- TRUE
  }
  list(perm = perm, sgn = sgn)
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


test_that("information SEs populate the rotated SE/CI schema under an oblique rotation", {
  R <- test_models$baseline$cormat
  fit <- EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "oblimin", se = "information")

  expect_setequal(names(fit$boot.SE),
                  c("unrot_loadings", "uniquenesses", "rot_loadings",
                    "communalities", "Phi", "Structure"))
  expect_setequal(names(fit$boot.CI),
                  c("unrot_loadings", "uniquenesses", "rot_loadings",
                    "communalities", "Phi", "Structure"))
  expect_null(fit$boot.arrays)

  expect_equal(dim(fit$boot.SE$rot_loadings), dim(fit$rot_loadings))
  expect_equal(dim(fit$boot.SE$Structure), dim(fit$rot_loadings))
  expect_equal(dim(fit$boot.SE$Phi), dim(fit$Phi))
  expect_length(fit$boot.SE$communalities, ncol(R))
  expect_false(anyNA(fit$boot.SE$rot_loadings))

  # The unit diagonal of Phi is fixed, so it carries no sampling variance.
  expect_true(all(diag(fit$boot.SE$Phi) == 0))

  # Wald intervals bracket every rotated point estimate.
  expect_true(all(fit$boot.CI$rot_loadings$lower <= fit$rot_loadings &
                    fit$rot_loadings <= fit$boot.CI$rot_loadings$upper))
  expect_true(all(fit$boot.CI$Structure$lower <= fit$Structure &
                    fit$Structure <= fit$boot.CI$Structure$upper))
  expect_named(fit$boot.CI$Phi, c("lower", "upper"))
})


test_that("information SEs under an orthogonal rotation omit Phi and the structure matrix", {
  R <- test_models$baseline$cormat
  fit <- EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "varimax", se = "information")

  expect_setequal(names(fit$boot.SE),
                  c("unrot_loadings", "uniquenesses", "rot_loadings", "communalities"))
  expect_null(fit$boot.SE$Phi)
  expect_null(fit$boot.SE$Structure)
  expect_equal(dim(fit$boot.SE$rot_loadings), dim(fit$rot_loadings))
  expect_false(anyNA(fit$boot.SE$rot_loadings))
  expect_false(anyNA(fit$boot.SE$communalities))
})


test_that("information SEs are produced for every supported native rotation", {
  # The lavaan oracle below validates the SE magnitude for a representative rotation; this guards
  # the remaining native criteria against a wrong criterion mapping in `.rotation_se_method` or a
  # warm-start reproduction failure that would silently degrade their rotated SEs to NA.
  R <- test_models$baseline$cormat
  orth <- c("quartimax", "equamax", "bentlerT", "geominT", "bifactorT")
  oblq <- c("quartimin", "bentlerQ", "geominQ", "bifactorQ")
  for (rot in c(orth, oblq)) {
    fit <- EFA(R, n_factors = 3, N = 500, method = "ML", rotation = rot, se = "information")
    expect_false(anyNA(fit$boot.SE$rot_loadings), info = rot)
    expect_equal(dim(fit$boot.SE$rot_loadings), dim(fit$rot_loadings), info = rot)
    expect_false(anyNA(fit$boot.SE$communalities), info = rot)
    if (rot %in% oblq) {
      expect_false(anyNA(fit$boot.SE$Phi), info = rot)
      expect_false(anyNA(fit$boot.SE$Structure), info = rot)
      expect_true(all(diag(fit$boot.SE$Phi) == 0), info = rot)
    } else {
      expect_null(fit$boot.SE$Phi, info = rot)
    }
  }
})


test_that("bifactor rotated SEs survive a general-factor reorder", {
  # When a group factor has a larger sum of squared loadings than the general factor,
  # `.reflect_and_order` moves the general factor out of the first column. The analytic SE
  # re-rotation must exempt the general factor wherever it landed (not blindly column 1), or it
  # cannot reproduce the reported loadings and the SEs collapse to NA. A weak general factor (0.30
  # on every variable) with two strong group factors (0.75) forces that reorder.
  set.seed(1)
  p <- 12
  Lp <- matrix(0, p, 3); Lp[, 1] <- 0.30; Lp[1:6, 2] <- 0.75; Lp[7:12, 3] <- 0.75
  Sig <- Lp %*% t(Lp); diag(Sig) <- 1
  X <- matrix(stats::rnorm(800 * p), 800) %*% chol(Sig)
  colnames(X) <- paste0("v", seq_len(p))

  for (rot in c("bifactorT", "bifactorQ")) {
    fit <- EFA(X, n_factors = 3, method = "ML", rotation = rot, se = "information")
    # the general factor (the column loading on every variable) is not in column 1 here
    expect_gt(which.max(apply(abs(unclass(fit$rot_loadings)), 2, min)), 1L)
    expect_false(anyNA(fit$boot.SE$rot_loadings), info = rot)
  }
})


test_that("rotated information loading and Phi SEs match lavaan's delta method", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  set.seed(42)
  R0 <- test_models$baseline$cormat
  X <- matrix(stats::rnorm(1500 * ncol(R0)), 1500) %*% chol(R0)
  colnames(X) <- colnames(R0)
  N <- 1500

  # normalize = FALSE matches lavaan's row.weights = "none"; oblimin (gam = 0) is quartimin, the
  # criterion lavaan rotates to under rotation = "oblimin".
  fit <- suppressWarnings(
    EFA(X, n_factors = 3, method = "ML", rotation = "oblimin", normalize = FALSE,
        se = "information"))
  ef <- lavaan::efa(sample.cov = stats::cor(X), sample.nobs = N, nfactors = 3,
                    rotation = "oblimin", rotation.se = "delta",
                    rotation.args = list(row.weights = "none"))
  lf <- Filter(function(e) inherits(e, "lavaan"), ef)[[1]]
  Ll <- lavaan::lavInspect(lf, "est")$lambda
  SEl <- lavaan::lavInspect(lf, "se")$lambda

  Lt <- unclass(fit$rot_loadings)
  al <- .align_factor_cols(Ll, Lt)
  Lt_a <- sweep(Lt[, al$perm, drop = FALSE], 2, al$sgn, "*")
  SEt_a <- fit$boot.SE$rot_loadings[, al$perm, drop = FALSE]

  # Both packages reach the same (identification-invariant) rotated solution.
  expect_equal(max(abs(Lt_a - Ll)), 0, tolerance = 0.01)
  # Their delta-method loading SEs agree to finite-sample noise: lavaan scales by N and EFAtools by
  # N - 1, and the optima differ at O(1 / sqrt(N)).
  expect_gt(stats::cor(c(SEt_a), c(SEl)), 0.97)
  expect_lt(mean(abs(SEt_a - SEl)), 0.004)

  # Factor-correlation SEs, matched by correlation value (the two packages order factors differently
  # so the off-diagonal positions do not line up).
  pe <- lavaan::parameterEstimates(lf)
  fc <- pe[pe$op == "~~" & pe$lhs != pe$rhs, ]
  ours <- data.frame(corr = fit$Phi[lower.tri(fit$Phi)],
                     se = fit$boot.SE$Phi[lower.tri(fit$boot.SE$Phi)])
  ours <- ours[order(ours$corr), ]
  fc <- fc[order(fc$est), ]
  expect_lt(max(abs(ours$se - fc$se)), 0.01)
})


test_that("rotated information SEs track the bootstrap", {
  skip_on_cran()
  skip_if_not_slow()

  set.seed(7)
  R0 <- test_models$baseline$cormat
  X <- matrix(stats::rnorm(2000 * ncol(R0)), 2000) %*% chol(R0)
  colnames(X) <- colnames(R0)

  fa <- EFA(X, n_factors = 3, method = "ML", rotation = "oblimin", se = "information")
  fb <- EFA(X, n_factors = 3, method = "ML", rotation = "oblimin",
            se = "np-boot", b_boot = 300, seed = 1)

  al <- .align_factor_cols(unclass(fa$rot_loadings), unclass(fb$rot_loadings))
  SEb <- fb$boot.SE$rot_loadings[, al$perm, drop = FALSE]

  # Bootstrap rotated-loading SEs carry target-rotation (Procrustes) alignment noise that the
  # analytic delta method does not, so the two are different finite-sample estimators (Yuan &
  # Hayashi, 2006) that agree in overall magnitude (median ratio near 1) rather than cell by cell.
  expect_lt(abs(stats::median(fa$boot.SE$rot_loadings / SEb) - 1), 0.3)

  SEs_b <- fb$boot.SE$Structure[, al$perm, drop = FALSE]
  expect_lt(abs(stats::median(fa$boot.SE$Structure / SEs_b) - 1), 0.3)

  # Communalities are rotation-invariant (no alignment), so they agree in pattern as well as
  # magnitude; the bootstrap has no communality slot, so they are recomputed from the replicate
  # loadings.
  comm_b <- apply(fb$boot.arrays$unrot_loadings, 3, function(L) rowSums(L^2))
  comm_b_se <- apply(comm_b, 1, stats::sd, na.rm = TRUE)
  expect_gt(stats::cor(fa$boot.SE$communalities, comm_b_se), 0.55)
  expect_lt(abs(stats::median(fa$boot.SE$communalities / comm_b_se) - 1), 0.4)
})


test_that("a Heywood case under a rotation yields NA rotated SEs with a classed warning", {
  # A communality above one drives a uniqueness below its zero boundary, where the expected
  # information is undefined, so no analytic SE -- rotated or unrotated -- exists.
  L <- matrix(c(sqrt(1.11), 0.6, 0.5, 0.4, 0.2, 0.7,
                0.1, 0.05, 0.2, 0.3, 0.15, 0.25), 6, 2)   # h2[1] > 1, so psi[1] < 0
  fit_out <- list(unrot_loadings = L)
  rot_info <- list(rotation = "oblimin", rotmat = diag(2), rot_loadings = L,
                   Phi = diag(2), normalize = FALSE, crit_args = list(gam = 0, delta = 0.01))

  expect_warning(
    out <- EFAtools:::.se_information_rotated(fit_out, rot_info, N = 200, ci = 0.95),
    class = "efa_se_unreliable"
  )
  expect_true(anyNA(out$SE$rot_loadings))
  expect_true(anyNA(out$SE$Phi))
  expect_true(anyNA(out$SE$communalities))
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
  # simplimax has a non-smooth (piecewise) criterion, so it has no usable analytic rotation Jacobian.
  expect_error(
    EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "simplimax", se = "information"),
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
