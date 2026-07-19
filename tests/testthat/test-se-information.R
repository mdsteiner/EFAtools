# Analytic expected-information standard errors for the ML solution.

# Independent reference for the correlation-structure information. theta = vec(Lambda) alone (the
# analysed diagonal is fixed at 1, so the uniquenesses are a function of Lambda, not free
# parameters): Gamma is assembled by an elementwise double loop over the pair index using the
# Olkin-Siotani formula, the model Jacobian Delta and the constraint Jacobian by central finite
# differences, and the information Delta' Gamma^-1 Delta is bordered, inverted and scaled by
# 1 / (N - 1). This shares no code with the package's vectorised assembly.
.ref_se_information <- function(L, psi, N) {
  p <- nrow(L); k <- ncol(L); pk <- p * k
  Sigma <- tcrossprod(L); diag(Sigma) <- 1
  prs <- utils::combn(p, 2L); n <- ncol(prs)

  # N Cov(r_ij, r_kl) under normality, one pair combination at a time.
  G <- matrix(0, n, n)
  for (s in seq_len(n)) for (t in seq_len(n)) {
    i <- prs[1, s]; j <- prs[2, s]; kk <- prs[1, t]; l <- prs[2, t]
    G[s, t] <- 0.5 * Sigma[i, j] * Sigma[kk, l] *
      (Sigma[i, kk]^2 + Sigma[i, l]^2 + Sigma[j, kk]^2 + Sigma[j, l]^2) +
      Sigma[i, kk] * Sigma[j, l] + Sigma[i, l] * Sigma[j, kk] -
      Sigma[i, j] * (Sigma[i, kk] * Sigma[i, l] + Sigma[j, kk] * Sigma[j, l]) -
      Sigma[kk, l] * (Sigma[i, kk] * Sigma[j, kk] + Sigma[i, l] * Sigma[j, l])
  }

  # sigma_offdiag(vec Lambda) and the gauge constraint off-diag(Lambda' Psi^-1 Lambda) = 0, both
  # differentiated by central finite differences.
  sig <- function(th) { Lm <- matrix(th, p, k); tcrossprod(Lm)[t(prs)] }
  g <- function(th) { Lm <- matrix(th, p, k); ps <- 1 - rowSums(Lm^2)
                      M <- crossprod(Lm, Lm / ps); M[upper.tri(M)] }
  th0 <- as.vector(L); nc <- k * (k - 1L) / 2L
  Delta <- matrix(0, n, pk); Cmat <- matrix(0, nc, pk); h <- 1e-6
  for (jj in seq_len(pk)) {
    tp <- th0; tp[jj] <- tp[jj] + h; tm <- th0; tm[jj] <- tm[jj] - h
    Delta[, jj] <- (sig(tp) - sig(tm)) / (2 * h)
    if (nc > 0) Cmat[, jj] <- (g(tp) - g(tm)) / (2 * h)
  }

  A <- crossprod(Delta, solve(G, Delta))
  Aug <- if (nc > 0) rbind(cbind(A, t(Cmat)), cbind(Cmat, matrix(0, nc, nc))) else A
  V <- solve(Aug)[seq_len(pk), seq_len(pk), drop = FALSE] / (N - 1)
  se <- sqrt(diag(V))
  # h2_i = rowSums(Lambda^2)_i has gradient 2 Lambda[i, ]; psi_i = 1 - h2_i shares its SE.
  G_h <- matrix(0, p, pk)
  for (i in seq_len(p)) G_h[i, (seq_len(k) - 1L) * p + i] <- 2 * L[i, ]
  list(loadings_se = matrix(se, p, k),
       uniquenesses_se = sqrt(rowSums((G_h %*% V) * G_h)))
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


test_that("information SEs match an independent correlation-structure reference", {
  L <- matrix(c(0.70, 0.60, 0.50, 0.10, 0.00, 0.20,
                0.10, 0.05, 0.20, 0.70, 0.60, 0.50), 6, 2)
  psi <- 1 - rowSums(L^2)
  fit_out <- list(unrot_loadings = L)

  ours <- EFAtools:::.se_information(fit_out, rot_info = NULL, N = 250, ci = 0.95,
                                     method = "ML")
  ref <- .ref_se_information(L, psi, N = 250)

  expect_equal(unname(ours$SE$unrot_loadings), ref$loadings_se, tolerance = 1e-5)
  expect_equal(unname(ours$SE$uniquenesses), ref$uniquenesses_se, tolerance = 1e-5)
})


test_that("the normal-theory Gamma reproduces the known variance of a correlation", {
  # The diagonal of the Olkin-Siotani covariance is the classical asymptotic variance of a
  # Pearson correlation, N Var(r_ij) = (1 - rho_ij^2)^2, which pins the formula's scale and its
  # pair ordering independently of the off-diagonal terms.
  R <- test_models$baseline$cormat
  prs <- utils::combn(ncol(R), 2L)
  G <- EFAtools:::.normal_theory_gamma(R, prs)

  expect_equal(unname(diag(G)), unname((1 - R[t(prs)]^2)^2), tolerance = 1e-12)
  expect_true(isSymmetric(G, tol = 1e-12))
  # It is the covariance of the sample correlations, so it must be positive definite at a
  # positive-definite Sigma -- the information Delta' Gamma^-1 Delta is built by inverting it.
  expect_gt(min(eigen(G, symmetric = TRUE, only.values = TRUE)$values), 0)
})


test_that("information uniqueness and loading SEs agree with the robust sandwich", {
  skip_on_cran()

  # `information` and `sandwich` build the same correlation-structure covariance from different
  # meats -- the normal-theory Gamma at the model-implied Sigma versus the fourth-moment (ADF)
  # Gamma of the observed correlations -- so on multivariate normal data they estimate the same
  # quantity and agree up to finite-sample meat noise. That makes the sandwich an independent
  # check on the uniqueness SEs within the correlation metric. lavaan's default fit is NOT the
  # comparator here: it analyses a covariance structure, so its uniqueness SEs carry the sampling
  # variability of a diagonal that is fixed at 1 in the model EFAtools fits (and under
  # `correlation = TRUE` lavaan derives the uniquenesses and reports their SE as 0).
  set.seed(11)
  R0 <- test_models$baseline$cormat
  X <- matrix(stats::rnorm(3000 * ncol(R0)), 3000) %*% chol(R0)
  colnames(X) <- colnames(R0)

  fi <- EFA(X, n_factors = 3, method = "ML", rotation = "none", se = "information")
  fs <- EFA(X, n_factors = 3, method = "ML", rotation = "none", se = "sandwich")

  expect_lt(abs(stats::median(fi$SE$uniquenesses / fs$SE$uniquenesses) - 1), 0.1)
  expect_lt(abs(stats::median(fi$SE$unrot_loadings / fs$SE$unrot_loadings) - 1), 0.1)
  expect_gt(stats::cor(c(fi$SE$unrot_loadings), c(fs$SE$unrot_loadings)), 0.9)
})


test_that("single-factor information loading SEs match lavaan's correlation structure", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  # With one factor there is no rotational indeterminacy, so the loadings are identified (up to
  # sign) the same way in both packages and their SEs are directly comparable. lavaan must be
  # asked for a CORRELATION structure (`correlation = TRUE`): its default fits a covariance
  # structure, which treats the analysed diagonal as carrying sampling variability and returns
  # loading SEs roughly twice as large. The residual gap is the N vs N - 1 scaling plus lavaan
  # evaluating its weight at the sample rather than the model-implied matrix.
  R <- stats::cor(GRiPS_raw)
  N <- nrow(GRiPS_raw)
  fit <- EFA(GRiPS_raw, n_factors = 1, method = "ML", rotation = "none", se = "information")

  mod <- paste0("f =~ ", paste(colnames(R), collapse = " + "))
  lav_fit <- lavaan::cfa(mod, data = as.data.frame(GRiPS_raw), std.lv = TRUE,
                         correlation = TRUE)
  pe <- lavaan::parameterEstimates(lav_fit)
  ld <- pe[pe$op == "=~", ]; ld <- ld[match(colnames(R), ld$rhs), ]

  # The residual few percent is the evaluation point: EFAtools uses the EXPECTED information, so
  # Gamma is taken at the model-implied Sigma, whereas lavaan evaluates it at the sample matrix.
  expect_equal(as.vector(fit$SE$unrot_loadings), ld$se, tolerance = 0.08)

  # ... and NOT lavaan's covariance-structure SEs, which are materially larger. Pinned as a
  # direction so a regression back to the covariance-structure information cannot pass silently.
  lav_cov <- lavaan::cfa(mod, data = as.data.frame(GRiPS_raw), std.lv = TRUE)
  se_cov <- lavaan::parameterEstimates(lav_cov)
  se_cov <- se_cov[se_cov$op == "=~", ]
  se_cov <- se_cov[match(colnames(R), se_cov$rhs), "se"]
  expect_gt(stats::median(se_cov / as.vector(fit$SE$unrot_loadings)), 1.2)
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

  expect_gt(cor(c(fa$SE$unrot_loadings), c(fb$SE$unrot_loadings)), 0.9)
  expect_lt(abs(median(fa$SE$unrot_loadings / fb$SE$unrot_loadings) - 1), 0.35)
})


test_that("information SEs populate the bootstrap SE/CI schema", {
  R <- test_models$baseline$cormat
  fit <- EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "none", se = "information")

  expect_setequal(names(fit$SE), c("unrot_loadings", "uniquenesses"))
  expect_setequal(names(fit$CI), c("unrot_loadings", "uniquenesses"))
  expect_null(fit$replicates)

  expect_equal(dim(fit$SE$unrot_loadings), dim(fit$unrot_loadings))
  expect_length(fit$SE$uniquenesses, ncol(R))
  expect_false(anyNA(fit$SE$unrot_loadings))

  # Wald intervals bracket the point estimates.
  expect_true(all(fit$CI$unrot_loadings$lower <= fit$unrot_loadings))
  expect_true(all(fit$unrot_loadings <= fit$CI$unrot_loadings$upper))
  expect_named(fit$CI$uniquenesses, c("lower", "upper"))
})


test_that("information SEs populate the rotated SE/CI schema under an oblique rotation", {
  R <- test_models$baseline$cormat
  fit <- EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "oblimin", se = "information")

  expect_setequal(names(fit$SE),
                  c("unrot_loadings", "uniquenesses", "rot_loadings",
                    "communalities", "Phi", "Structure"))
  expect_setequal(names(fit$CI),
                  c("unrot_loadings", "uniquenesses", "rot_loadings",
                    "communalities", "Phi", "Structure"))
  expect_null(fit$replicates)

  expect_equal(dim(fit$SE$rot_loadings), dim(fit$rot_loadings))
  expect_equal(dim(fit$SE$Structure), dim(fit$rot_loadings))
  expect_equal(dim(fit$SE$Phi), dim(fit$Phi))
  expect_length(fit$SE$communalities, ncol(R))
  expect_false(anyNA(fit$SE$rot_loadings))

  # The unit diagonal of Phi is fixed, so it carries no sampling variance.
  expect_true(all(diag(fit$SE$Phi) == 0))

  # Wald intervals bracket every rotated point estimate.
  expect_true(all(fit$CI$rot_loadings$lower <= fit$rot_loadings &
                    fit$rot_loadings <= fit$CI$rot_loadings$upper))
  expect_true(all(fit$CI$Structure$lower <= fit$Structure &
                    fit$Structure <= fit$CI$Structure$upper))
  expect_named(fit$CI$Phi, c("lower", "upper"))
})


test_that("information communality SEs equal the uniqueness SEs (one estimand)", {
  # h2_i = 1 - psi_i exactly, so the communality and uniqueness are a single estimand up to sign
  # and must carry the same standard error (and a mirrored Wald interval). Both are the ordinary
  # delta method on the loading covariance, so they now agree by construction -- under the former
  # covariance-structure information the two routes disagreed by up to 18%.
  R <- test_models$baseline$cormat
  fit <- EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "oblimin", se = "information")

  expect_equal(unname(fit$SE$communalities), unname(fit$SE$uniquenesses))
  # The Wald intervals mirror: communality upper = 1 - uniqueness lower (and lower = 1 - upper).
  expect_equal(unname(fit$CI$communalities$upper), unname(1 - fit$CI$uniquenesses$lower))
  expect_equal(unname(fit$CI$communalities$lower), unname(1 - fit$CI$uniquenesses$upper))
})


test_that("an unusable covariance NAs the communality SEs along with the uniqueness SEs", {
  # The post-solve PSD gate leaves a FINITE covariance next to NA marginal SEs (only the marginals
  # are NA-filled). The communalities are a delta method on that covariance, so they must fail
  # closed with it -- otherwise the object ships a finite `SE$communalities` beside an all-NA
  # `SE$uniquenesses` and an all-NA `vcov_unrot_loadings`, contradicting both the one-estimand
  # contract and the "no finite SE next to an NA covariance" invariant.
  R <- test_models$baseline$cormat
  fit <- EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "oblimin", se = "none")
  fit_out <- list(unrot_loadings = fit$unrot_loadings)
  rot_info <- list(rotation = "oblimin", rotmat = fit$rotmat, rot_loadings = fit$rot_loadings,
                   Phi = fit$Phi, normalize = fit$settings$normalize,
                   crit_args = list(gam = 0, delta = 0.01))

  out <- testthat::with_mocked_bindings(
    suppressWarnings(
      EFAtools:::.se_information(fit_out, rot_info, N = 500, ci = 0.95, method = "ML")),
    .is_psd = function(M) FALSE,
    .package = "EFAtools"
  )

  expect_true(all(is.na(out$SE$uniquenesses)))
  expect_true(all(is.na(out$SE$communalities)))
  expect_true(all(is.na(out$SE$rot_loadings)))
  expect_true(all(is.na(out$vcov_unrot_loadings)))
})


test_that("information SEs under an orthogonal rotation omit Phi and the structure matrix", {
  R <- test_models$baseline$cormat
  fit <- EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "varimax", se = "information")

  expect_setequal(names(fit$SE),
                  c("unrot_loadings", "uniquenesses", "rot_loadings", "communalities"))
  expect_null(fit$SE$Phi)
  expect_null(fit$SE$Structure)
  expect_equal(dim(fit$SE$rot_loadings), dim(fit$rot_loadings))
  expect_false(anyNA(fit$SE$rot_loadings))
  expect_false(anyNA(fit$SE$communalities))
})


test_that("information SEs are produced for every supported native rotation", {
  # The lavaan oracle below validates the SE magnitude for a representative rotation; this guards
  # the remaining native criteria against a wrong criterion mapping in `.rotation_se_method` or a
  # warm-start reproduction failure that would silently degrade their rotated SEs to NA.
  # One representative per family: a smooth CF / oblimin path (quartimax / quartimin) and the
  # bifactor general-factor-exempt path (bifactorT / bifactorQ). The bifactor REORDER path is
  # exercised by the next block; this block covers the non-reorder case.
  R <- test_models$baseline$cormat
  orth <- c("quartimax", "bifactorT")
  oblq <- c("quartimin", "bifactorQ")
  for (rot in c(orth, oblq)) {
    fit <- EFA(R, n_factors = 3, N = 500, method = "ML", rotation = rot, se = "information")
    expect_false(anyNA(fit$SE$rot_loadings), info = rot)
    expect_equal(dim(fit$SE$rot_loadings), dim(fit$rot_loadings), info = rot)
    expect_false(anyNA(fit$SE$communalities), info = rot)
    if (rot %in% oblq) {
      expect_false(anyNA(fit$SE$Phi), info = rot)
      expect_false(anyNA(fit$SE$Structure), info = rot)
      expect_true(all(diag(fit$SE$Phi) == 0), info = rot)
    } else {
      expect_null(fit$SE$Phi, info = rot)
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
    expect_false(anyNA(fit$SE$rot_loadings), info = rot)
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
  SEt_a <- fit$SE$rot_loadings[, al$perm, drop = FALSE]

  # Both packages reach the same (identification-invariant) rotated solution.
  expect_equal(max(abs(Lt_a - Ll)), 0, tolerance = 0.01)
  # Their delta-method loading SEs agree closely overall -- the same rotation Jacobian applied to
  # the same solution -- but not exactly, because `lavaan::efa()` propagates a COVARIANCE-structure
  # unrotated covariance (it has no correlation-structure EFA path) while EFAtools propagates the
  # correlation-structure one.
  expect_gt(stats::cor(c(SEt_a), c(SEl)), 0.93)
  expect_lt(mean(abs(SEt_a - SEl)), 0.004)

  # That difference has a signature, and it is the one that would disappear if the unrotated
  # covariance regressed to a covariance structure: it is confined to the SALIENT loadings (where
  # the spurious standardisation variance accumulates) and leaves the near-zero loadings alone.
  salient <- abs(Ll) > 0.4
  expect_lt(max((SEt_a / SEl)[salient]), 0.98)
  expect_lt(abs(stats::median((SEt_a / SEl)[!salient]) - 1), 0.05)

  # Factor-correlation SEs, matched by correlation value (the two packages order factors differently
  # so the off-diagonal positions do not line up).
  pe <- lavaan::parameterEstimates(lf)
  fc <- pe[pe$op == "~~" & pe$lhs != pe$rhs, ]
  ours <- data.frame(corr = fit$Phi[lower.tri(fit$Phi)],
                     se = fit$SE$Phi[lower.tri(fit$SE$Phi)])
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
  SEb <- fb$SE$rot_loadings[, al$perm, drop = FALSE]

  # Bootstrap rotated-loading SEs carry target-rotation (Procrustes) alignment noise that the
  # analytic delta method does not, so the two are different finite-sample estimators (Yuan &
  # Hayashi, 2006) that agree in overall magnitude (median ratio near 1) rather than cell by cell.
  expect_lt(abs(stats::median(fa$SE$rot_loadings / SEb) - 1), 0.3)

  SEs_b <- fb$SE$Structure[, al$perm, drop = FALSE]
  expect_lt(abs(stats::median(fa$SE$Structure / SEs_b) - 1), 0.3)

  # Communalities are rotation-invariant (no alignment), so they are compared cell by cell. The
  # analytic communality SE is the delta-method value it shares with the uniqueness
  # (h2 = 1 - psi); the bootstrap has no communality slot, so the reference is recomputed from the
  # replicate loadings. The two agree in magnitude (median ratio near 1); the normal-theory SE and
  # the bootstrap are different finite-sample estimators, so the per-variable pattern correlation is
  # positive but looser than the magnitude agreement.
  comm_b <- apply(fb$replicates$unrot_loadings, 3, function(L) rowSums(L^2))
  comm_b_se <- apply(comm_b, 1, stats::sd, na.rm = TRUE)
  expect_gt(stats::cor(fa$SE$communalities, comm_b_se), 0.4)
  expect_lt(abs(stats::median(fa$SE$communalities / comm_b_se) - 1), 0.4)
})


test_that("a Heywood case under a rotation yields NA rotated SEs with a classed warning", {
  # A communality above one drives a uniqueness below its zero boundary, where the solution sits
  # on the parameter-space boundary and the Wald intervals this path reports are not valid, so no
  # analytic SE -- rotated or unrotated -- is offered.
  L <- matrix(c(sqrt(1.11), 0.6, 0.5, 0.4, 0.2, 0.7,
                0.1, 0.05, 0.2, 0.3, 0.15, 0.25), 6, 2)   # h2[1] > 1, so psi[1] < 0
  fit_out <- list(unrot_loadings = L)
  rot_info <- list(rotation = "oblimin", rotmat = diag(2), rot_loadings = L,
                   Phi = diag(2), normalize = FALSE, crit_args = list(gam = 0, delta = 0.01))

  expect_warning(
    out <- EFAtools:::.se_information(fit_out, rot_info, N = 200, ci = 0.95, method = "ML"),
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

  expect_equal(raw$SE$unrot_loadings, cmat$SE$unrot_loadings, tolerance = 1e-5)
  expect_equal(raw$SE$uniquenesses, cmat$SE$uniquenesses, tolerance = 1e-5)
})


test_that("a Heywood case yields NA SEs with a classed warning", {
  # A communality above one drives a uniqueness below its zero boundary. The correlation-structure
  # information is still assembleable there, but the solution sits on the parameter-space boundary,
  # where a Wald interval is not valid, so the path declines to report one.
  L <- matrix(c(sqrt(1.11), 0.6, 0.5, 0.4), 4, 1)   # h2[1] = 1.11 > 1, so psi[1] < 0
  fit_out <- list(unrot_loadings = L)

  expect_warning(
    out <- EFAtools:::.se_information(fit_out, rot_info = NULL, N = 200, ci = 0.95,
                                      method = "ML"),
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
  # Sandwich SEs require an ordinal correlation method; a continuous (Pearson) matrix is rejected.
  expect_error(
    EFA(R, n_factors = 3, N = 500, method = "ML", rotation = "none", se = "sandwich"),
    class = "efa_se_unsupported"
  )
  expect_error(
    EFA(R, n_factors = 3, N = 500, method = "ULS", rotation = "none", se = "sandwich"),
    class = "efa_se_unsupported"
  )
  # A correlation matrix carries no sample size, which the analytic methods require.
  expect_error(
    EFA(R, n_factors = 3, method = "ML", rotation = "none", se = "information"),
    class = "efa_se_no_n"
  )
})


test_that("information SEs are rejected for non-Pearson correlation methods", {
  # The expected information is built from the normal-theory covariance of Pearson correlations.
  # A polychoric correlation additionally carries first-stage threshold estimation error and a rank
  # correlation is not a Pearson moment at all, so neither may reach this path silently -- the
  # sandwich is the supported route for both.
  x <- DOSPERT_raw[1:300, 1:8]

  for (cm in c("poly", "tetra", "spearman", "kendall")) {
    expect_error(
      EFA(x, n_factors = 2, method = "ML", rotation = "none", se = "information",
          cor_method = cm),
      class = "efa_se_unsupported",
      info = cm
    )
  }

  # Pearson and fiml stay available.
  expect_s3_class(
    EFA(x, n_factors = 2, method = "ML", rotation = "none", se = "information",
        cor_method = "pearson"),
    "EFA"
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
