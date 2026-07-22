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
  # The two factors carry deliberately different strengths, so the ML canonical variances
  # diag(Lambda' Psi^-1 Lambda) are well separated (4.23 vs 1.16) and the rotational gauge is
  # firmly identified. Equal strengths leave the orientation unpinned, and the unrotated loading
  # SEs then diverge along the gauge -- a regime the path deliberately withholds (see the gauge
  # tests below), which would leave nothing here to compare against the reference.
  L <- matrix(c(0.80, 0.75, 0.70, 0.10, 0.00, 0.20,
                0.10, 0.05, 0.20, 0.55, 0.50, 0.45), 6, 2)
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

  expect_setequal(names(fit$SE), c("unrot_loadings", "uniquenesses", "communalities"))
  expect_setequal(names(fit$CI), c("unrot_loadings", "uniquenesses", "communalities"))
  expect_null(fit$replicates)

  expect_equal(dim(fit$SE$unrot_loadings), dim(fit$unrot_loadings))
  expect_length(fit$SE$uniquenesses, ncol(R))
  expect_false(anyNA(fit$SE$unrot_loadings))

  # h2 = 1 - psi exactly, so the unrotated path reports the communality SEs it already holds,
  # with the mirrored Wald interval.
  expect_equal(unname(fit$SE$communalities), unname(fit$SE$uniquenesses))
  expect_equal(unname(fit$CI$communalities$upper), unname(1 - fit$CI$uniquenesses$lower))
  expect_equal(unname(fit$CI$communalities$lower), unname(1 - fit$CI$uniquenesses$upper))

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
  # every native criterion against a wrong criterion mapping in `.rotation_se_method` or a
  # warm-start reproduction failure that would silently degrade its rotated SEs to NA. All eleven
  # mapped rotations are covered: equamax (param = k / (2 p)) and geominT/geominQ (param = delta)
  # are the only entries whose Jacobian parameter is COMPUTED rather than fixed, so a transposed
  # k / p or a drifted delta default would hide there. The bifactor REORDER path is exercised by
  # the next block; this block covers the non-reorder case.
  R <- test_models$baseline$cormat
  orth <- c("varimax", "quartimax", "equamax", "bentlerT", "geominT", "bifactorT")
  oblq <- c("oblimin", "quartimin", "bentlerQ", "geominQ", "bifactorQ")
  for (rot in c(orth, oblq)) {
    # The native criterion rotations run `random_starts` random starts, so without a seed inside
    # the loop which local optimum is tested depends on the RNG state left by the preceding
    # iteration (and by preceding test blocks).
    set.seed(51)
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
    # The two group factors are equally strong by construction, so their canonical variances
    # nearly coincide and the unrotated orientation is weakly determined. The fit therefore
    # warns (`efa_se_unreliable`) and withholds the unrotated loading SEs -- that is the
    # documented gauge behaviour, not a failure, and the point of this test is that it must
    # not take the rotated SEs down with it. The warning is suppressed rather than asserted
    # by class because it sits on a numeric threshold and need not fire identically on every
    # BLAS; the substantive assertion below is what this test pins.
    fit <- suppressWarnings(
      EFA(X, n_factors = 3, method = "ML", rotation = rot, se = "information"))
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
  if (utils::packageVersion("lavaan") >= "0.7.0") {
    ef <- lavaan::efa(sample.cov = stats::cor(X), sample.nobs = N, nfactors = 3,
                      rotation = list("oblimin", row_weights = "none"),
                      rotation.se = "delta")
  } else {
    # lavaan < 0.7 used a separate rotation.args argument.
    ef <- lavaan::efa(sample.cov = stats::cor(X), sample.nobs = N, nfactors = 3,
                      rotation = "oblimin", rotation.se = "delta",
                      rotation.args = list(row.weights = "none"))
  }
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

  # TODO: `SE$Structure` has no external oracle here. `lavaan::efa()` reports no standard error
  # for the structure matrix, so validating it against lavaan would mean re-deriving
  # J_S = (Phi' (x) I_p) J_L + (I_k (x) L) J_Phi from lavaan's own parameter covariance inside the
  # test -- i.e. re-implementing the code under test. It is currently checked for magnitude
  # against a 300-replicate bootstrap in the block below (slow-gated) only.
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


test_that("a uniqueness pinned at the fitter's floor is treated as a Heywood case", {
  # The ML/ULS/DWLS fitters constrain the uniquenesses to [.uniqueness_floor, 1], so an improper
  # solution is pinned AT that floor and never reaches zero. A gate keyed on `psi <= 0` therefore
  # fires only for a hand-built covariance and never for a fitted one, letting a boundary solution
  # report a hugely inflated Wald standard error. Pin the boundary the fitters can actually reach.
  L <- matrix(c(sqrt(1 - EFAtools:::.uniqueness_floor), 0.6, 0.5, 0.4,
                0.00, 0.05, 0.20, 0.30), 4, 2)
  psi <- 1 - rowSums(L^2)
  expect_lte(min(psi), EFAtools:::.uniqueness_floor + sqrt(.Machine$double.eps))
  expect_gt(min(psi), 0)                       # strictly interior to the old psi <= 0 gate

  expect_warning(
    out <- EFAtools:::.se_information(list(unrot_loadings = L), rot_info = NULL,
                                      N = 250, ci = 0.95, method = "ML"),
    class = "efa_se_unreliable"
  )
  expect_true(all(is.na(out$SE$unrot_loadings)))
  expect_true(all(is.na(out$SE$uniquenesses)))
})


test_that("a Heywood case from an actual fit withholds the analytic SEs", {
  # The companion test above drives `.se_information()` with a hand-built psi, which bypasses the
  # optimiser -- and the optimiser is precisely what makes the boundary reachable only AT the floor
  # rather than below zero. Drive the public path on data that fits an improper solution, so the
  # gate is exercised against a uniqueness the fitter actually produced.
  set.seed(2)
  Lh <- matrix(0, 6, 2)
  Lh[1:3, 1] <- c(.95, .55, .50)
  Lh[4:6, 2] <- c(.90, .50, .45)
  Rh <- tcrossprod(Lh)
  diag(Rh) <- 1
  x <- matrix(stats::rnorm(150 * 6), 150, 6) %*% chol(Rh)
  colnames(x) <- paste0("V", seq_len(6))

  fit <- suppressMessages(suppressWarnings(
    efa_fit(x, n_factors = 2, estimator = "ML", rotation = "none", se = "information")
  ))
  skip_if(length(fit$heywood) == 0L, "this platform's optimiser did not land on the boundary")

  # The uniqueness is pinned AT the floor, strictly above the zero the old gate tested for.
  psi <- 1 - rowSums(unclass(fit$unrot_loadings)^2)
  expect_lte(min(psi), EFAtools:::.uniqueness_floor + sqrt(.Machine$double.eps))
  expect_gt(min(psi), 0)

  expect_true(all(is.na(fit$SE$unrot_loadings)))
  expect_true(all(is.na(fit$SE$uniquenesses)))
  expect_true(anyNA(fit$vcov_unrot_loadings))
})


test_that("a weakly determined rotational gauge withholds only the unrotated loading SEs", {
  # When two ML canonical variances diag(Lambda' Psi^-1 Lambda) coincide, the orientation within
  # that two-plane is arbitrary: the gauge-fixing constraint stops being a transversal to the gauge
  # orbit, the reflexive inverse is amplified by ||(C Z)^-1||, and the unrotated loadings have no
  # well-defined limiting distribution. The covariance stays finite and PSD, so it passes the
  # PSD gate -- without the transversal check the path would report a standard error of order 10 or
  # 100 for a parameter bounded in [-1, 1].
  L <- matrix(c(0.70, 0.60, 0.50, 0.10, 0.00, 0.20,
                0.10, 0.05, 0.20, 0.70, 0.60, 0.50), 6, 2)
  psi <- 1 - rowSums(L^2)
  canon <- diag(crossprod(L, L / psi))
  expect_lt(abs(diff(canon)) / max(canon), 1e-2)       # the two-plane really is degenerate

  Cmat <- EFAtools:::.se_sandwich_constraint(L, psi = psi)
  expect_lt(EFAtools:::.gauge_transversal(L, Cmat), EFAtools:::.gauge_transversal_floor)

  expect_warning(
    out <- EFAtools:::.se_information(list(unrot_loadings = L), rot_info = NULL,
                                      N = 250, ci = 0.95, method = "ML"),
    class = "efa_se_unreliable"
  )
  expect_true(all(is.na(out$SE$unrot_loadings)))
  # The communalities rowSums(Lambda^2) are invariant under Lambda -> Lambda T, so they survive.
  expect_false(anyNA(out$SE$communalities))
  expect_equal(unname(out$SE$uniquenesses), unname(out$SE$communalities))
  # The covariance itself is correct; only its gauge-dependent marginals are unusable.
  expect_false(anyNA(out$vcov_unrot_loadings))
})


test_that("a weakly determined gauge leaves the rotated SEs intact", {
  # The rotation Jacobian annihilates the gauge directions, so the rotated loadings, Phi and the
  # structure coefficients are gauge-invariant and must NOT be withheld along with the unrotated
  # ones -- otherwise the gate would destroy the very quantities that remain trustworthy.
  # Drive the whole public path so the rotation, its Jacobian and the gauge diagnostic all come
  # from the same fitted solution rather than a hand-built rot_info. Two disjoint blocks make
  # Lambda' Psi^-1 Lambda exactly diagonal, so the population solution IS the canonical one; the
  # 0.98 scaling puts the two canonical variances close together (2.371 / 2.203) without making
  # them equal, which is the near-degenerate regime the gate targets. At exactly equal strengths
  # the bordered matrix is exactly singular and the pre-existing singular-information path takes
  # over instead.
  Lp <- matrix(0, 6, 2)
  Lp[1:3, 1] <- c(0.72, 0.65, 0.60)
  Lp[4:6, 2] <- c(0.72, 0.65, 0.60) * 0.98
  R <- tcrossprod(Lp)
  diag(R) <- 1
  colnames(R) <- rownames(R) <- paste0("V", seq_len(6))

  expect_warning(
    fit <- efa_fit(R, n_factors = 2, N = 250, estimator = "ML", rotation = "oblimin",
                   se = "information"),
    class = "efa_se_unreliable"
  )
  expect_true(all(is.na(fit$SE$unrot_loadings)))
  expect_false(anyNA(fit$SE$rot_loadings))
  expect_false(anyNA(fit$SE$communalities))
  expect_false(anyNA(fit$SE$Phi))
  # ... and the surviving rotated SEs are of a plausible magnitude, not gauge-inflated.
  expect_lt(max(fit$SE$rot_loadings), 0.5)
  # The covariance the rotated SEs were propagated from stays intact, so consumers that pool in a
  # common gauge across fits (efa_mi) keep the gauge-invariant quantities.
  expect_false(anyNA(fit$vcov_unrot_loadings))
})


test_that("a well-separated gauge is not flagged and reports finite SEs", {
  # Guard the other direction: the transversal check must not fire on an ordinary solution.
  L <- matrix(c(0.80, 0.75, 0.70, 0.10, 0.00, 0.20,
                0.10, 0.05, 0.20, 0.55, 0.50, 0.45), 6, 2)
  psi <- 1 - rowSums(L^2)
  Cmat <- EFAtools:::.se_sandwich_constraint(L, psi = psi)
  expect_gt(EFAtools:::.gauge_transversal(L, Cmat), EFAtools:::.gauge_transversal_floor)

  expect_no_warning(
    out <- EFAtools:::.se_information(list(unrot_loadings = L), rot_info = NULL,
                                      N = 250, ci = 0.95, method = "ML"),
    class = "efa_se_unreliable"
  )
  expect_false(anyNA(out$SE$unrot_loadings))
  expect_lt(max(out$SE$unrot_loadings), 0.5)
})


test_that("the gauge transversal is the exact conditioning of the bordered inverse", {
  # `.gauge_transversal()` restricts the constraint to null(A) = {vec(Lambda S) : S antisymmetric},
  # using that basis in closed form rather than decomposing A. The claim rests on the gauge
  # direction leaving the modelled off-diagonal of Sigma = Lambda Lambda' unchanged, which is
  # checked here directly: the change along Lambda S must be second order in the step, so the
  # information A = Delta' Gamma^-1 Delta annihilates it whatever Gamma is.
  L <- matrix(c(0.70, 0.60, 0.50, 0.10, 0.00, 0.20,
                0.10, 0.05, 0.20, 0.66, 0.58, 0.47), 6, 2)
  psi <- 1 - rowSums(L^2)
  k <- ncol(L)
  nc <- k * (k - 1L) / 2L

  S <- matrix(c(0, -1, 1, 0), 2, 2)
  off <- function(M) { diag(M) <- 0; max(abs(M)) }
  step <- vapply(c(1e-2, 1e-3), function(e) off(tcrossprod(L + e * (L %*% S)) - tcrossprod(L)),
                 numeric(1))
  # Halving the step ten-fold cuts the deviation a hundred-fold: second order, no first-order part.
  expect_equal(step[1] / step[2], 100, tolerance = 1e-3)

  # The value agrees with a basis obtained by decomposing the information matrix instead. Build A
  # through the package's own assembly via the sandwich core at the efficient weight, then read its
  # null space off the smallest `nc` eigenvectors.
  Sigma <- tcrossprod(L); diag(Sigma) <- 1
  prs <- utils::combn(nrow(L), 2L)
  Gamma <- EFAtools:::.normal_theory_gamma(Sigma, prs)
  Delta <- t(vapply(seq_len(ncol(prs)), function(s) {
    G <- matrix(0, nrow(L), k)
    G[prs[1, s], ] <- L[prs[2, s], ]
    G[prs[2, s], ] <- G[prs[2, s], ] + L[prs[1, s], ]
    as.vector(G)
  }, numeric(length(L))))
  A <- crossprod(Delta, solve(Gamma) %*% Delta)
  ev <- eigen((A + t(A)) / 2, symmetric = TRUE)
  Z_eig <- ev$vectors[, seq(ncol(A) - nc + 1L, ncol(A)), drop = FALSE]

  Cmat <- EFAtools:::.se_sandwich_constraint(L, psi = psi)
  rn <- sqrt(rowSums(Cmat^2))
  expect_equal(EFAtools:::.gauge_transversal(L, Cmat),
               min(svd((Cmat / rn) %*% Z_eig)$d), tolerance = 1e-8)

  # A single factor has no rotational freedom, so there is no transversal to condition.
  expect_identical(
    EFAtools:::.gauge_transversal(L[, 1, drop = FALSE],
                                  EFAtools:::.se_sandwich_constraint(L[, 1, drop = FALSE])),
    Inf
  )
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
