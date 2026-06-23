# EFA_POOLED() reports the incremental fit indices CFI (Bentler, 1990) and TLI
# (Tucker & Lewis, 1973) as the average of the per-imputation indices. Each
# per-imputation index is computed in range by the per-model path, so the average
# stays in range and consistent with the component fits, and reduces to EFA()'s
# value in the single-imputation limit. Pooling the model and baseline chi-squares
# separately (the lavaan.mi/semTools convention) shrinks the two noncentralities by
# their own between-imputation ARIV and can drive these non-linear indices outside
# [0, 1] and against every per-imputation value; those separately pooled
# noncentralities remain available, for reference, in mi_diagnostics. These tests
# check the averaging convention.
#
# Assertion groups:
#   1. Single-imputation (identical-imputation) limit reduces to EFA().
#   2. Pooled CFI/TLI equal the mean of the per-imputation indices (not the
#      separately pooled (N - 1)-scale formula).
#   3. Pooled CFI/TLI stay in range and within the per-imputation range under
#      genuine multi-imputation fits where separate pooling did not.

# ---- Group 1 ---------------------------------------------------------------

test_that("identical imputations reproduce EFA()'s CFI/TLI exactly", {
  cormat <- test_models$baseline$cormat
  N <- 500L
  k <- 3L
  single <- EFA(cormat, n_factors = k, N = N, method = "ML", rotation = "none")
  pooled <- suppressMessages(EFA_POOLED(
    replicate(4, cormat, simplify = FALSE),
    n_factors = k, N = N, method = "ML", rotation = "none"
  ))

  # Identical imputations have zero between-imputation spread, so every per-
  # imputation CFI/TLI equals the single-fit value and their average must match
  # EFA() to machine precision.
  expect_equal(pooled$fit_indices$CFI, single$fit_indices$CFI, tolerance = 1e-8)
  expect_equal(pooled$fit_indices$TLI, single$fit_indices$TLI, tolerance = 1e-8)
})

# ---- Group 2 ---------------------------------------------------------------

test_that("pooled CFI/TLI equal the mean of the per-imputation indices", {
  set.seed(42)
  grips_list <- lapply(1:3, function(i) {
    GRiPS_raw[sample(seq_len(nrow(GRiPS_raw)), replace = TRUE), ]
  })
  pooled <- suppressWarnings(suppressMessages(
    EFA_POOLED(grips_list, n_factors = 1, method = "ML")
  ))
  fi <- pooled$fit_indices

  cfis <- vapply(pooled$fits, function(f) f$fit_indices$CFI, numeric(1))
  tlis <- vapply(pooled$fits, function(f) f$fit_indices$TLI, numeric(1))

  expect_equal(fi$CFI, mean(cfis[is.finite(cfis)]), tolerance = 1e-12)
  expect_equal(fi$TLI, mean(tlis[is.finite(tlis)]), tolerance = 1e-12)

  # The reported indices are the average, not the discarded separately pooled
  # (N - 1)-scale formula. That formula is still recoverable from the diagnostics
  # mi_diagnostics exposes, and here it gives a different number, so a revert to it
  # would be caught.
  md <- pooled$mi_diagnostics
  d_m <- max(md$chi_cfi - fi$df, 0)
  d_null <- max(md$chi_null_cfi - fi$df_null, 0)
  CFI_sep <- if (max(d_m, d_null) == 0) 1 else 1 - d_m / max(d_m, d_null)
  ratio <- md$chi_null_cfi / fi$df_null
  TLI_sep <- (ratio - md$chi_cfi / fi$df) / (ratio - 1)
  expect_false(isTRUE(all.equal(fi$CFI, CFI_sep)))
  expect_false(isTRUE(all.equal(fi$TLI, TLI_sep)))
})

# ---- Group 3 ---------------------------------------------------------------

test_that("pooled CFI/TLI stay in range and consistent under multi-imputation fits", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  # A 9-variable population with no clean two-factor structure: a k = 2 EFA
  # underfits, so each imputation lands at a moderate CFI/TLI well inside (0, 1).
  # Pooling the model and baseline chi-squares separately drove the pooled TLI
  # negative (and the pooled CFI far below every per-imputation value); averaging
  # the per-imputation indices keeps both in range and within the component range.
  set.seed(1)
  p <- 9L
  Sigma <- stats::cov2cor(tcrossprod(matrix(stats::rnorm(p * 4), p)) + diag(p))

  check_in_range <- function(data_list, k) {
    pooled <- suppressWarnings(suppressMessages(
      EFA_POOLED(data_list, n_factors = k, method = "ML", rotation = "none")
    ))
    fi <- pooled$fit_indices
    cfis <- vapply(pooled$fits, function(f) f$fit_indices$CFI, numeric(1))
    tlis <- vapply(pooled$fits, function(f) f$fit_indices$TLI, numeric(1))

    # CFI is always in [0, 1] (mean of in-range per-imputation values).
    expect_gte(fi$CFI, 0)
    expect_lte(fi$CFI, 1)
    # Both indices lie within the per-imputation range -- the consistency the
    # separately pooled statistics violated (its TLI fell below every component).
    expect_gte(fi$CFI, min(cfis) - 1e-8)
    expect_lte(fi$CFI, max(cfis) + 1e-8)
    expect_gte(fi$TLI, min(tlis) - 1e-8)
    expect_lte(fi$TLI, max(tlis) + 1e-8)
    invisible(pooled)
  }

  # Homogeneous N.
  imps_hom <- lapply(1:5, function(i) {
    as.data.frame(MASS::mvrnorm(300, rep(0, p), Sigma))
  })
  check_in_range(imps_hom, k = 2)

  # Heterogeneous N: the degenerate case where each separately pooled chi-square
  # floored to zero, returning the self-contradictory CFI = 1 / TLI = 0.
  Ns <- c(200, 400, 600, 250, 550)
  imps_het <- lapply(Ns, function(n) {
    as.data.frame(MASS::mvrnorm(n, rep(0, p), Sigma))
  })
  check_in_range(imps_het, k = 2)
})
