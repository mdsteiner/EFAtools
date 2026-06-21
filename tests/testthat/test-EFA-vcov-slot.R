# Tests for the persisted unrotated loading vcov and correlation ACOV slots on EFA() results.
# Contracts pinned here:
#   * `vcov_unrot_loadings`: pk x pk numeric matrix; the FULL unrotated loading covariance the
#     marginal `SE$unrot_loadings` were derived from. Populated for se = "information" and
#     se = "sandwich" (in both unrotated and rotated paths), NA-filled at a Heywood/unreliable
#     fit, NULL for se = "none" and se = "np-boot".
#   * `Gamma`: n x n asymptotic covariance of the off-diagonal correlations on the
#     variance scale (n = p(p-1)/2, ordered by utils::combn(p, 2)). Populated for
#     se = "sandwich" only, NULL otherwise.
# Downstream MI pooling depends on these slots; the assertions below are the storage contract.

# Self-contained continuous fixture (skewed simple-structure factor model, ~unit variance per
# column). Independent of test-se-sandwich.R so the file can run in isolation.
make_adf_fixture <- function(N, loadings, seed) {
  set.seed(seed)
  p <- nrow(loadings)
  k <- ncol(loadings)
  skew <- function(n) (stats::rchisq(n, df = 1) - 1) / sqrt(2)
  Fmat <- matrix(skew(N * k), N, k)
  E <- matrix(skew(N * p), N, p)
  uniq <- sqrt(pmax(1 - rowSums(loadings^2), 0))
  X <- Fmat %*% t(loadings) + sweep(E, 2L, uniq, "*")
  colnames(X) <- paste0("V", seq_len(p))
  as.data.frame(X)
}


test_that("se = 'information' persists the full unrotated loading vcov (unrotated path)", {
  cormat <- test_models$baseline$cormat
  p <- nrow(cormat); k <- 3L; pk <- p * k

  fit <- EFA(cormat, n_factors = k, N = 500, method = "ML", rotation = "none",
             se = "information")

  expect_true("vcov_unrot_loadings" %in% names(fit))
  V <- fit$vcov_unrot_loadings
  expect_true(is.matrix(V) && is.numeric(V))
  expect_equal(dim(V), c(pk, pk))
  expect_true(isSymmetric(V, tol = 1e-10))

  # The persisted vcov must be the EXACT matrix that produced the marginal SEs:
  # sqrt(diag(V)) reshaped column-major == SE$unrot_loadings element-wise.
  diag_se <- sqrt(pmax(diag(V), 0))
  expect_equal(matrix(diag_se, p, k), unname(fit$SE$unrot_loadings), tolerance = 1e-12)
})


test_that("se = 'information' persists the unrotated loading vcov under rotation", {
  cormat <- test_models$baseline$cormat
  p <- nrow(cormat); k <- 3L; pk <- p * k

  fit <- EFA(cormat, n_factors = k, N = 500, method = "ML", rotation = "oblimin",
             se = "information")

  expect_true("vcov_unrot_loadings" %in% names(fit))
  V <- fit$vcov_unrot_loadings
  expect_equal(dim(V), c(pk, pk))
  expect_true(isSymmetric(V, tol = 1e-10))

  # Still the UNROTATED block: its diagonal must match the unrotated loading SEs, NOT the
  # rotated ones. This is what U2/M2 will pool across imputations.
  diag_se <- sqrt(pmax(diag(V), 0))
  expect_equal(matrix(diag_se, p, k), unname(fit$SE$unrot_loadings), tolerance = 1e-12)
})


test_that("se = 'sandwich' persists the unrotated robust loading vcov on the polychoric path", {
  skip_on_cran()

  dat <- DOSPERT_raw[, 1:8]
  p <- ncol(dat); k <- 1L; pk <- p * k

  fit <- EFA(dat, n_factors = k, cor_method = "poly", method = "ULS", rotation = "none",
             se = "sandwich")

  expect_true("vcov_unrot_loadings" %in% names(fit))
  V <- fit$vcov_unrot_loadings
  expect_equal(dim(V), c(pk, pk))
  expect_true(isSymmetric(V, tol = 1e-10))

  diag_se <- sqrt(pmax(diag(V), 0))
  expect_equal(matrix(diag_se, p, k), unname(fit$SE$unrot_loadings), tolerance = 1e-12)
})


test_that("se = 'sandwich' persists the unrotated robust loading vcov on the continuous path", {
  dat <- make_adf_fixture(N = 500,
                          loadings = matrix(c(.7, .65, .6, .75, .55, .5), ncol = 1L),
                          seed = 11)
  p <- ncol(dat); k <- 1L; pk <- p * k

  fit <- EFA(dat, n_factors = k, cor_method = "pearson", method = "ML", rotation = "none",
             se = "sandwich")

  V <- fit$vcov_unrot_loadings
  expect_equal(dim(V), c(pk, pk))
  expect_true(isSymmetric(V, tol = 1e-10))

  diag_se <- sqrt(pmax(diag(V), 0))
  expect_equal(matrix(diag_se, p, k), unname(fit$SE$unrot_loadings), tolerance = 1e-12)
})


test_that("se = 'sandwich' under rotation also persists the unrotated block", {
  L <- matrix(0, 8, 2)
  L[1:4, 1] <- c(.7, .65, .6, .55)
  L[5:8, 2] <- c(.7, .6, .65, .5)
  dat <- make_adf_fixture(N = 700, loadings = L, seed = 12)
  p <- ncol(dat); k <- 2L; pk <- p * k

  fit <- EFA(dat, n_factors = k, cor_method = "pearson", method = "ML",
             rotation = "oblimin", se = "sandwich")

  V <- fit$vcov_unrot_loadings
  expect_equal(dim(V), c(pk, pk))
  expect_true(isSymmetric(V, tol = 1e-10))
  diag_se <- sqrt(pmax(diag(V), 0))
  expect_equal(matrix(diag_se, p, k), unname(fit$SE$unrot_loadings), tolerance = 1e-12)
})


test_that("an unreliable information fit yields an NA-filled vcov of the correct shape, not NULL", {
  # Drive `.se_information_ml()` straight into its `na_out` branch with a synthetic loading
  # whose first communality exceeds one (psi[1] < 0 -> Psi^{-1} and the identification
  # constraint Lambda' Psi^{-1} Lambda are no longer defined). Going through EFA() is not a
  # reliable route here because the optimiser pulls the boundary off zero (e.g. psi[1] ~ 0.005
  # at the .999-loading fixture), so the SE path completes with finite numbers; the storage
  # contract under genuine unreliability is the one we need to pin.
  L <- matrix(c(sqrt(1.11), 0.6, 0.5, 0.4), 4, 1)   # h2[1] = 1.11 > 1, so psi[1] < 0
  fit_out <- list(unrot_loadings = L)
  p <- nrow(L); k <- ncol(L); pk <- p * k

  out <- suppressWarnings(EFAtools:::.se_information(fit_out, N = 200, ci = 0.95))

  expect_true("vcov_unrot_loadings" %in% names(out))
  V <- out$vcov_unrot_loadings
  # The unreliability gate must surface as a present, correctly shaped NA matrix so
  # pooling code can fail closed on anyNA() rather than on a missing slot.
  expect_false(is.null(V))
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(pk, pk))
  expect_true(all(is.na(V)))
})


test_that("an unreliable sandwich fit yields an NA-filled V_AA of the correct shape, not NULL", {
  # `.se_sandwich_core()` returns `na_core` (NA-filled V_AA, reliable = FALSE) whenever
  # Gamma is missing, the wrong shape, or carries an NA. Force that path with an
  # invalid Gamma; the wrapper must still attach a correctly shaped NA matrix on the
  # persisted slot.
  L <- matrix(c(.7, .65, .6, .55), 4, 1)
  fit_out <- list(unrot_loadings = L)
  p <- nrow(L); k <- ncol(L); pk <- p * k

  core <- EFAtools:::.se_sandwich_core(fit_out, N = 200, Gamma = NULL, method = "ULS")
  out <- suppressWarnings(EFAtools:::.se_sandwich_unrotated(fit_out, core, ci = 0.95))

  expect_true("vcov_unrot_loadings" %in% names(out))
  V <- out$vcov_unrot_loadings
  expect_false(is.null(V))
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(pk, pk))
  expect_true(all(is.na(V)))
})


test_that("a finite-but-not-PSD information vcov is NA-filled to match the marginal SEs", {
  # `.se_information_ml()` has TWO unreliability paths: (i) early Heywood / singular return
  # already covered above (matrix-of-NAs vcov); (ii) post-solve sqrt-gate where the bordered
  # inverse succeeds but yields a V whose diagonal carries a small negative variance, so the
  # marginal SEs are NA but `vcov = V` stays finite. The persisted slot must NA-fill on the
  # second path too, so a non-PSD V never ships next to NA marginal SEs and silently
  # propagates `sqrt(NaN)` downstream.
  L <- matrix(c(.7, .65, .6, .55, .5, .6), 6, 1)
  fit_out <- list(unrot_loadings = L)
  p <- nrow(L); k <- ncol(L); pk <- p * k
  q <- pk + p

  # Simulate the post-solve gate output: a finite (non-NA) bordered vcov alongside NA marginal
  # SEs. The persistence wrapper should still ship an all-NA pk x pk slot.
  fake_ml <- function(L, psi, N, R = NULL) {
    list(
      vcov = matrix(1.0, q, q),         # finite, deliberately non-PSD-looking
      loadings_se = matrix(NA_real_, p, k),
      uniquenesses_se = rep(NA_real_, p)
    )
  }

  out <- testthat::with_mocked_bindings(
    suppressWarnings(EFAtools:::.se_information(fit_out, N = 200, ci = 0.95)),
    .se_information_ml = fake_ml,
    .package = "EFAtools"
  )

  V <- out$vcov_unrot_loadings
  expect_false(is.null(V))
  expect_equal(dim(V), c(pk, pk))
  expect_true(all(is.na(V)))
})


test_that("a finite-but-not-PSD sandwich V_AA is NA-filled to match the marginal SEs", {
  # `.se_sandwich_core()` returns `V_AA = V_AA` (unchanged) even when `reliable = FALSE`,
  # only the marginal `loadings_se` is NA-filled in that branch. The persistence wrapper
  # must mirror the marginal NA into the persisted vcov.
  L <- matrix(c(.7, .65, .6, .55), 4, 1)
  fit_out <- list(unrot_loadings = L)
  p <- nrow(L); k <- ncol(L); pk <- p * k

  fake_core <- list(
    V_AA = matrix(1.0, pk, pk),         # finite numeric
    loadings_se = matrix(NA_real_, p, k),
    uniquenesses_se = rep(NA_real_, p),
    scaled_test = NULL,
    reliable = FALSE
  )
  out <- suppressWarnings(EFAtools:::.se_sandwich_unrotated(fit_out, fake_core, ci = 0.95))

  V <- out$vcov_unrot_loadings
  expect_false(is.null(V))
  expect_equal(dim(V), c(pk, pk))
  expect_true(all(is.na(V)))
})


test_that("the top-level name vector under analytic SE is pinned", {
  # The 9 `expect_named()` blocks in test-EFA.R only fire on se = "none" fits, where the
  # analytic-SE branch is skipped, so the order `..., settings, SE, CI, replicates,
  # vcov_unrot_loadings` is never asserted there. Pin it here so a future refactor that
  # reorders the slot assignments around `.compute_se_ci()` cannot silently drift.
  fit <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500, method = "ML",
             rotation = "none", se = "information")
  expect_named(fit, c("orig_R", "h2", "orig_eigen", "final_eigen", "iter",
                      "convergence", "heywood", "unrot_loadings", "vars_accounted",
                      "fit_indices", "model_implied_R", "residuals", "settings",
                      "vcov_unrot_loadings", "Gamma", "SE", "CI", "replicates"))
})


test_that("se = 'none' leaves vcov_unrot_loadings present but NULL", {
  fit <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500, method = "ML",
             rotation = "none", se = "none")

  expect_true("vcov_unrot_loadings" %in% names(fit))
  expect_null(fit$vcov_unrot_loadings)
})


test_that("se = 'np-boot' leaves vcov_unrot_loadings present but NULL", {
  # Tiny bootstrap; this test only asserts the storage convention.
  dat <- make_adf_fixture(N = 200,
                          loadings = matrix(c(.7, .6, .55, .65, .5, .5), ncol = 1L),
                          seed = 13)
  fit <- suppressWarnings(EFA(dat, n_factors = 1, cor_method = "pearson", method = "ML",
                              rotation = "none", se = "np-boot", b_boot = 20, seed = 1))

  expect_true("vcov_unrot_loadings" %in% names(fit))
  expect_null(fit$vcov_unrot_loadings)
})


test_that("se = 'sandwich' persists Gamma on the polychoric path with the right schema", {
  skip_on_cran()

  dat <- DOSPERT_raw[, 1:8]
  p <- ncol(dat); n_pairs <- p * (p - 1L) / 2L

  fit <- EFA(dat, n_factors = 1, cor_method = "poly", method = "ULS", rotation = "none",
             se = "sandwich")

  expect_true("Gamma" %in% names(fit))
  G <- fit$Gamma
  expect_true(is.matrix(G) && is.numeric(G))
  expect_equal(dim(G), c(n_pairs, n_pairs))
  expect_true(isSymmetric(G, tol = 1e-10))
  expect_true(all(diag(G) > 0))

  # Pair ordering matches utils::combn(p, 2) -- the dimnames are set by .polychoric() in
  # that order. Reconstruct the expected labels and confirm.
  vn <- colnames(dat)
  labels <- apply(utils::combn(vn, 2L), 2L, paste, collapse = "-")
  expect_equal(rownames(G), labels)
  expect_equal(colnames(G), labels)
})


test_that("se = 'sandwich' persists Gamma on the continuous Pearson path with the right schema", {
  dat <- make_adf_fixture(N = 500,
                          loadings = matrix(c(.7, .65, .6, .75, .55, .5), ncol = 1L),
                          seed = 14)
  p <- ncol(dat); n_pairs <- p * (p - 1L) / 2L

  fit <- EFA(dat, n_factors = 1, cor_method = "pearson", method = "ML", rotation = "none",
             se = "sandwich")

  G <- fit$Gamma
  expect_true(is.matrix(G) && is.numeric(G))
  expect_equal(dim(G), c(n_pairs, n_pairs))
  expect_true(isSymmetric(G, tol = 1e-10))
  expect_true(all(diag(G) > 0))

  vn <- colnames(dat)
  labels <- apply(utils::combn(vn, 2L), 2L, paste, collapse = "-")
  expect_equal(rownames(G), labels)
  expect_equal(colnames(G), labels)
})


test_that("Gamma is present but NULL outside the sandwich path", {
  cormat <- test_models$baseline$cormat

  fit_none <- EFA(cormat, n_factors = 3, N = 500, method = "ML", rotation = "none",
                  se = "none")
  fit_info <- EFA(cormat, n_factors = 3, N = 500, method = "ML", rotation = "none",
                  se = "information")
  fit_boot <- suppressWarnings(
    EFA(cormat, n_factors = 3, N = 500, method = "ML", rotation = "none",
        se = "np-boot", b_boot = 10, seed = 1)
  )

  for (fit in list(fit_none, fit_info, fit_boot)) {
    expect_true("Gamma" %in% names(fit))
    expect_null(fit$Gamma)
  }
})


test_that("Gamma matches lavaan's correlation NACOV up to the N scale (continuous Pearson)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if(utils::packageVersion("lavaan") < "0.6.9")

  L <- matrix(c(.7, .65, .6, .75, .55, .5), ncol = 1L)
  dat <- make_adf_fixture(N = 600, loadings = L, seed = 15)
  vn <- colnames(dat)
  N <- nrow(dat)

  fit <- EFA(dat, n_factors = 1, cor_method = "pearson", method = "ML", rotation = "none",
             se = "sandwich")

  # The persisted Gamma is the variance-scale ADF covariance; lavaan's NACOV is N * Gamma.
  mod <- paste0("f =~ ", paste(vn, collapse = " + "))
  lfit <- lavaan::cfa(mod, data = dat, estimator = "MLM", std.lv = TRUE,
                      correlation = TRUE, meanstructure = TRUE)
  Glav <- lavaan::lavInspect(lfit, "gamma")
  if (is.list(Glav)) Glav <- Glav[[1]]
  # With meanstructure = TRUE lavaan returns the means PLUS the off-diagonal correlation rows
  # (in utils::combn() order); subset to the ~~ rows for a clean comparison.
  cor_rows <- grep("~~", rownames(Glav))
  Glav_corr <- unname(Glav[cor_rows, cor_rows])
  expect_equal(nrow(Glav_corr), choose(length(vn), 2L))

  rel_F <- norm(N * unname(fit$Gamma) - Glav_corr, "F") / norm(Glav_corr, "F")
  expect_lt(rel_F, 1e-4)
})


test_that("Gamma matches lavaan's polychoric NACOV up to the N scale (ordinal)", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  dat <- as.data.frame(DOSPERT_raw[, 1:8])
  vn <- colnames(dat)
  N <- sum(stats::complete.cases(dat))

  fit <- EFA(dat, n_factors = 1, cor_method = "poly", method = "ULS", rotation = "none",
             se = "sandwich")

  # lavaan::efa() returns an `efaList`; lavInspect expects a single lavaan fit, so take [[1]].
  lf <- lavaan::efa(dat, nfactors = 1, ordered = vn, estimator = "DWLS", se = "robust.sem")
  lfit <- lf[[1]]
  Glav <- lavaan::lavInspect(lfit, "gamma")
  if (is.list(Glav)) Glav <- Glav[[1]]
  cor_rows <- grep("~~", rownames(Glav))
  Glav_corr <- unname(Glav[cor_rows, cor_rows])
  expect_equal(nrow(Glav_corr), choose(length(vn), 2L))

  rel_F <- norm(N * unname(fit$Gamma) - Glav_corr, "F") / norm(Glav_corr, "F")
  expect_lt(rel_F, 1e-3)
})


test_that("the round-trip sqrt(diag(vcov_unrot_loadings)) reproduces SE$unrot_loadings", {
  # Belt-and-braces: even after assembly into the EFA object, the persisted vcov stays
  # self-consistent with the SE slot -- so callers reading either get the same number.
  fit <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500, method = "ML",
             rotation = "oblimin", se = "information")
  p <- nrow(fit$unrot_loadings); k <- ncol(fit$unrot_loadings)

  reconstructed <- matrix(sqrt(pmax(diag(fit$vcov_unrot_loadings), 0)), p, k)
  expect_equal(reconstructed, unname(fit$SE$unrot_loadings), tolerance = 1e-12)
})
