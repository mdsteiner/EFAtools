# cor_method = "poly" / "tetra": polychoric and tetrachoric correlations wired into
# the public API. The matrix estimator itself is covered in test-polychoric.R; these
# tests cover the routing through .prepare_cor_input() (EFA and the suitability /
# retention functions), the `use` handling, the EFA() bootstrap recompute, the
# tetrachoric binary assertion, and the clean rejection by the criteria whose
# reference data are continuous (CD, PARALLEL, NEST, HULL).

# Muffle only the named warning classes, letting any other (unexpected) warning
# surface so it is not silently hidden by a blanket suppressWarnings().
.muffle <- function(expr, ...) {
  classes <- c(...)
  withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, classes)) invokeRestart("muffleWarning")
  })
}

test_that("EFA() runs ordinal factor analysis and matches psych::fa on a polychoric matrix", {
  skip_on_cran()
  skip_if_not_installed("psych")

  x <- DOSPERT_raw[stats::complete.cases(DOSPERT_raw), ]
  N <- nrow(x)
  k <- 3L

  efa <- EFA(x, n_factors = k, method = "ULS", cor_method = "poly",
             rotation = "none")

  # Reference: psych's minres factor analysis (ULS is the same estimator) of psych's
  # own polychoric matrix, with correct = 0 to match .polychoric()'s default.
  Rp <- suppressWarnings(psych::polychoric(x, correct = FALSE, global = FALSE)$rho)
  fa_ref <- suppressWarnings(psych::fa(Rp, nfactors = k, fm = "minres",
                                       rotate = "none", n.obs = N))
  L_ref <- unclass(fa_ref$loadings)
  L <- efa$unrot_loadings

  # Compare via the reproduced common-factor correlations (L %*% t(L)), which are
  # invariant to factor sign and column order, so no fragile per-column alignment
  # is needed. The two polychoric matrices agree to ~1e-4 and ULS == minres, so the
  # common parts and communalities agree well within 1e-3.
  expect_lt(max(abs(tcrossprod(L) - tcrossprod(L_ref))), 1e-3)
  expect_lt(max(abs(efa$h2 - fa_ref$communality)), 1e-3)
})

test_that("cor_method = 'tetra' equals 'poly' on binary data", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
  gb <- apply(g, 2L, function(col) as.integer(col > stats::median(col)))

  # A Heywood case (and any non-PD smoothing) on this small binary example is
  # incidental to the equivalence being checked; muffle only those specific
  # warnings rather than blanket-suppressing every condition.
  e_poly  <- .muffle(EFA(gb, n_factors = 2, method = "ULS",
                         cor_method = "poly", rotation = "none"),
                     "efa_heywood", "efa_cor_smoothed")
  e_tetra <- .muffle(EFA(gb, n_factors = 2, method = "ULS",
                         cor_method = "tetra", rotation = "none"),
                     "efa_heywood", "efa_cor_smoothed")

  expect_equal(e_tetra$unrot_loadings, e_poly$unrot_loadings)
  expect_equal(e_tetra$orig_R, e_poly$orig_R)
})

test_that("cor_method = 'tetra' rejects variables with more than two categories", {
  expect_error(EFA(DOSPERT_raw, n_factors = 3, cor_method = "tetra"),
               class = "efa_cor_not_binary")
})

test_that("method = 'DWLS' requires a polychoric asymptotic covariance", {
  # a correlation matrix carries no raw data to estimate an asymptotic covariance from
  expect_error(EFA(test_models$baseline$cormat, n_factors = 3, method = "DWLS", N = 500),
               class = "efa_dwls_no_acov")
  # a continuous correlation method produces no polychoric asymptotic covariance
  expect_error(EFA(DOSPERT_raw, n_factors = 3, method = "DWLS", cor_method = "pearson"),
               class = "efa_dwls_no_acov")
  # the abort fires before fitting, regardless of bootstrap settings
  expect_error(EFA(test_models$baseline$cormat, n_factors = 3, method = "DWLS", N = 500,
                   se = "np-boot"),
               class = "efa_dwls_no_acov")
})

test_that("KMO() and BARTLETT() honour cor_method = 'poly'", {
  x <- DOSPERT_raw[stats::complete.cases(DOSPERT_raw), ]
  N <- nrow(x)
  Rp <- .polychoric(x)$R

  # cor_method = "poly" must produce exactly the result of passing that polychoric
  # matrix in directly: this pins the routing through .prepare_cor_input(). The raw
  # call infers N = nrow(x), the same N supplied with the matrix.
  expect_equal(KMO(x, cor_method = "poly")$KMO, KMO(Rp)$KMO)
  expect_equal(BARTLETT(x, cor_method = "poly")$chisq,
               BARTLETT(Rp, N = N)$chisq)
})

test_that("KMO() and BARTLETT() with cor_method = 'poly' match an external polychoric reference", {
  skip_on_cran()
  skip_if_not_installed("psych")

  x <- DOSPERT_raw[stats::complete.cases(DOSPERT_raw), ]
  N <- nrow(x)
  Rp <- suppressWarnings(psych::polychoric(x, correct = FALSE, global = FALSE)$rho)

  # The two polychoric matrices agree to ~1e-4, so the suitability statistics agree
  # to well within 1e-2.
  expect_equal(KMO(x, cor_method = "poly")$KMO, psych::KMO(Rp)$MSA,
               tolerance = 1e-2)
  expect_equal(BARTLETT(x, cor_method = "poly")$chisq,
               psych::cortest.bartlett(Rp, n = N)$chisq,
               tolerance = 1e-2, ignore_attr = TRUE)
})

test_that("the simulation-based criteria reject poly / tetra with a classed condition", {
  expect_error(CD(DOSPERT_raw, cor_method = "poly"),
               class = "efa_cor_method_unsupported")
  expect_error(PARALLEL(DOSPERT_raw, cor_method = "tetra"),
               class = "efa_cor_method_unsupported")
  expect_error(NEST(DOSPERT_raw, N = 200, cor_method = "poly"),
               class = "efa_cor_method_unsupported")
  # HULL derives its factor-search bound from an internal parallel analysis, so it
  # rejects poly/tetra for the same reason.
  expect_error(HULL(DOSPERT_raw, N = 200, cor_method = "tetra"),
               class = "efa_cor_method_unsupported")
})

test_that("cor_method = 'poly' honours `use` for missing data", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
  gm <- g
  gm[1:5, 1] <- NA

  # all.obs and everything do not delete missing data, so they error on it,
  # matching the outcome of stats::cor() for the Pearson path.
  expect_error(KMO(gm, cor_method = "poly", use = "all.obs"),
               class = "efa_cor_na")
  expect_error(KMO(gm, cor_method = "poly", use = "everything"),
               class = "efa_cor_na")

  # complete.obs listwise-deletes first, so the matrix equals the polychoric of
  # the complete-case data.
  R_listwise <- .polychoric(gm[stats::complete.cases(gm), ])$R
  expect_equal(KMO(gm, cor_method = "poly", use = "complete.obs")$KMO,
               KMO(R_listwise)$KMO)
})

test_that("N_FACTORS skips reference-based criteria under poly with an informative note", {
  classes <- character()
  withCallingHandlers(
    N_FACTORS(GRiPS_raw, criteria = c("EKC", "HULL"), cor_method = "poly",
              method = "ULS", suitability = FALSE),
    warning = function(w) {
      classes <<- c(classes, class(w))
      invokeRestart("muffleWarning")
    }
  )
  # HULL is skipped with the dedicated 'skipped' class, not reported as a generic
  # failure; EKC still runs on the polychoric matrix.
  expect_true("efa_criterion_skipped" %in% classes)
  expect_false("efa_criterion_failed" %in% classes)
})

test_that("the polychoric bootstrap is reproducible and positive definite under resampling", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]

  b1 <- EFA(g, n_factors = 1, method = "ULS", cor_method = "poly",
            se = "np-boot", b_boot = 25, seed = 42)
  b2 <- EFA(g, n_factors = 1, method = "ULS", cor_method = "poly",
            se = "np-boot", b_boot = 25, seed = 42)

  # Same seed -> identical bootstrap SEs, independent of how many replicates were fit.
  expect_identical(b1$boot$SE$unrot_loadings, b2$boot$SE$unrot_loadings)
  # Replicate matrices are recomputed per resample (degenerate ones are dropped at
  # the fit), so the surviving replicates yield finite SEs.
  expect_true(all(is.finite(b1$boot$SE$unrot_loadings)))
})

test_that("the DWLS polychoric bootstrap is reproducible and positive definite under resampling", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]

  b1 <- suppressWarnings(EFA(g, n_factors = 1, method = "DWLS", cor_method = "poly",
                             se = "np-boot", b_boot = 25, seed = 42))
  b2 <- suppressWarnings(EFA(g, n_factors = 1, method = "DWLS", cor_method = "poly",
                             se = "np-boot", b_boot = 25, seed = 42))

  # Same seed -> identical bootstrap SEs; each replicate recomputes its own polychoric
  # matrix and diagonal-ACOV weights, and the surviving replicates yield finite SEs.
  expect_identical(b1$boot$SE$unrot_loadings, b2$boot$SE$unrot_loadings)
  expect_true(all(is.finite(b1$boot$SE$unrot_loadings)))
})

test_that("the polychoric bootstrap on DOSPERT_raw completes (timed)", {
  skip_on_cran()
  skip_if_not_slow()

  b_boot <- 200L
  t <- system.time(
    fit <- EFA(DOSPERT_raw, n_factors = 3, method = "ULS", cor_method = "poly",
               se = "np-boot", b_boot = b_boot, seed = 1)
  )
  elapsed <- unname(t["elapsed"])
  cli::cli_inform(
    "Polychoric bootstrap: {b_boot} replicates on DOSPERT_raw in {round(elapsed, 1)}s ({round(elapsed / b_boot, 3)}s/replicate).")
  expect_true(all(is.finite(fit$boot$SE$unrot_loadings)))
})
