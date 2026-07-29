smt_cor <- efa_smt(test_models$baseline$cormat, N = 500)
smt_zero <- efa_smt(diag(nrow = 5, ncol = 5), N = 500)
smt_raw <- efa_smt(GRiPS_raw)

test_that("output class and dimensions are correct", {
  expect_s3_class(smt_cor, "efa_retention")
  expect_length(smt_cor, 6)
  expect_s3_class(smt_raw, "efa_retention")
  expect_length(smt_raw, 6)
  expect_s3_class(smt_zero, "efa_retention")
  expect_length(smt_zero, 6)

  expect_named(smt_cor$n_factors, c("chi", "RMSEA", "AIC"))
  expect_equal(.retention_record(smt_cor, "chi")$plot_type, "none")
})

test_that("number of factors are correct", {
  expect_equal(smt_cor$n_factors[["chi"]], 3)
  expect_equal(smt_cor$n_factors[["RMSEA"]], 2)
  expect_equal(smt_cor$n_factors[["AIC"]], 3)

  expect_equal(smt_raw$n_factors[["chi"]], 3)
  expect_equal(smt_raw$n_factors[["RMSEA"]], 1)
  expect_equal(smt_raw$n_factors[["AIC"]], 3)

  expect_equal(smt_zero$n_factors[["chi"]], 0)
  expect_equal(smt_zero$n_factors[["RMSEA"]], 0)
  expect_equal(smt_zero$n_factors[["AIC"]], 0)
})

test_that("p-values are correct", {
  # the chi record's y is c(p_null, ps_chi)
  chi_cor <- .retention_record(smt_cor, "chi")$y
  expect_lt(chi_cor[1], 0.05)
  expect_lt(chi_cor[2], 0.05)
  expect_lt(chi_cor[3], 0.05)
  expect_gte(chi_cor[4], 0.05)
  expect_gte(chi_cor[5], 0.05)

  chi_raw <- .retention_record(smt_raw, "chi")$y
  expect_lt(chi_raw[1], 0.05)
  expect_lt(chi_raw[2], 0.05)
  expect_lt(chi_raw[3], 0.05)
  expect_gte(chi_raw[4], 0.05)
  expect_gte(chi_raw[5], 0.05)

  chi_zero <- .retention_record(smt_zero, "chi")$y
  expect_gte(chi_zero[1], 0.05)
  expect_gte(chi_zero[2], 0.05)
})

test_that("RMSEA_LB and AIC values are correct", {
  rmsea_cor <- .retention_record(smt_cor, "RMSEA")$y
  rmsea_raw <- .retention_record(smt_raw, "RMSEA")$y
  rmsea_zero <- .retention_record(smt_zero, "RMSEA")$y

  expect_equal(rmsea_cor[1], 0.157856, tolerance = 1e-2)
  expect_equal(rmsea_raw[1], 0.461174, tolerance = 1e-2)
  expect_equal(rmsea_zero[1], 0, tolerance = 1e-2)

  expect_equal(rmsea_cor[-1], c(0.0567972, 0.0397977, rep(0, 10)),
               tolerance = 1e-2)
  expect_equal(rmsea_raw[-1], c(0.0354958, 0.0263924, rep(0, 2)),
               tolerance = 1e-2)
  expect_equal(rmsea_zero[-1], rep(0, 2), tolerance = 1e-2)

  aic_cor <- .retention_record(smt_cor, "AIC")$y
  aic_raw <- .retention_record(smt_raw, "AIC")$y
  aic_zero <- .retention_record(smt_zero, "AIC")$y

  expect_equal(aic_cor[1], 1867.28, tolerance = 0.1)
  expect_equal(aic_raw[1], 4998.06, tolerance = 0.1)
  expect_equal(aic_zero[1], -20, tolerance = 0.1)

  expect_equal(aic_cor[-1], c(133.20917, 13.14179, -80.25186, -79.07978,
                              -75.59415, -68.10822, -57.26154, -46.34868,
                              -36.42090, -26.76733, -17.03907, -5.13003),
               tolerance = 0.1)
  expect_equal(aic_raw[-1], c(19.610600, 7.631330, -4.512150, -1.868430),
               tolerance = 0.1)
  expect_equal(aic_zero[-1], c(-10, -2), tolerance = 0.1)
})

test_that("null-model statistics are computed from R and N, not the fitted model", {
  # The zero-factor (null) statistics depend only on the correlation matrix and
  # N. SMT must derive them directly so they stay finite even when the
  # max-factor model is degenerate (a Heywood / non-positive-definite case
  # leaves that model's fit_indices NA, which previously crashed the p_null and
  # RMSEA_LB_null comparisons).
  R <- test_models$baseline$cormat
  N <- 500

  # SMT smooths a non-positive-definite input before use (the trigger is a
  # smallest eigenvalue below .Machine$double.eps; see .prepare_cor_input). The
  # baseline matrix clears that threshold, so smoothing is a no-op and the R used
  # inside SMT equals the input here; enforce that so the recomputation is valid.
  expect_true(all(eigen(R, symmetric = TRUE, only.values = TRUE)$values >=
                    .Machine$double.eps))

  m <- ncol(R)
  chi_null <- .null_chisq(R, N)
  df_null <- (m^2 - m) / 2
  p_null <- stats::pchisq(chi_null, df_null, lower.tail = FALSE)
  # RMSEA is built on the uncorrected (N - 1) discrepancy scale, so the null-model
  # bound uses the uncorrected baseline chi-square (p_null and AIC_null keep the
  # Bartlett-corrected statistic).
  chi_null_rmsea <- .null_chisq(R, N, corrected = FALSE)
  RMSEA_LB_null <- sqrt(.rmsea_lambda(chi_null_rmsea, df_null, .95) / (df_null * (N - 1)))
  AIC_null <- chi_null - 2 * df_null

  # the null model is the first element of each criterion record's y vector
  null_out <- c(chi   = .retention_record(smt_cor, "chi")$y[1],
                RMSEA = .retention_record(smt_cor, "RMSEA")$y[1],
                AIC   = .retention_record(smt_cor, "AIC")$y[1])
  expect_equal(null_out, c(chi = p_null, RMSEA = RMSEA_LB_null, AIC = AIC_null))

  # finite null statistics in the SMT output are what keep the comparison guards
  # from crashing
  expect_false(anyNA(null_out))
})

test_that("SMT propagates NA null statistics (no crash) for tiny N", {
  # When N is small relative to the number of variables the Bartlett null
  # multiplier N - 1 - (2p + 5)/6 is non-positive, so .null_chisq() is NA. SMT
  # must carry that through to NA null statistics instead of crashing in the
  # RMSEA noncentrality search (.rmsea_lambda) with an opaque base error.
  R <- stats::cor(GRiPS_raw)  # 8 variables -> null multiplier <= 0 for N <= 4
  smt_tiny <- suppressWarnings(efa_smt(R, N = 4))
  expect_s3_class(smt_tiny, "efa_retention")
  expect_true(is.na(.retention_record(smt_tiny, "chi")$y[1]))
  expect_true(is.na(.retention_record(smt_tiny, "RMSEA")$y[1]))

  # Both sequential rules share one convention: an undefined value breaks the
  # strictly sequential test, so the search stops there and no number is suggested
  # rather than skipping ahead to a later model that meets the rule. Here that
  # happens at the very first (null) model.
  expect_true(is.na(smt_tiny$n_factors[["chi"]]))
  expect_true(is.na(smt_tiny$n_factors[["RMSEA"]]))
})

test_that("the chi-square sequence matches stats::factanal", {
  # Two independent ML optimisers on the same matrix: the agreement is tight but
  # not exact, and how tight depends on the BLAS/LAPACK the two share. Pin it
  # where the numerical environment is known rather than across every check
  # flavour (as the comparison-data pin in test-efa_cd.R does).
  skip_on_cran()
  R <- test_models$baseline$cormat
  N <- 500
  # the chi record's y is c(p_null, p_1, ..., p_max_fac)
  ps <- .retention_record(smt_cor, "chi")$y

  compared <- 0L
  for (k in seq_len(.det_max_factors(ncol(R)))) {
    fa <- stats::factanal(covmat = R, factors = k, n.obs = N)
    # Solutions with a uniqueness at the lower bound are boundary solutions: the
    # two optimisers can then attain different local optima (EFAtools reaches the
    # lower discrepancy on this matrix), so only interior solutions are comparable.
    if (min(fa$uniquenesses) <= .uniqueness_floor + 1e-8) next
    expect_equal(ps[k + 1], fa$PVAL[[1]], tolerance = 1e-6)
    compared <- compared + 1L
  }
  # guard against a vacuously passing loop
  expect_gt(compared, 0)
})

test_that("an inadmissible selected solution raises a classed warning", {
  # On this matrix the chi-square and AIC rules both select a 6-factor model that
  # has Heywood cases, while the RMSEA rule selects an admissible 4-factor model.
  expect_warning(smt_bad <- efa_smt(test_models$case_11b$cormat, N = 500),
                 class = "efa_smt_inadmissible")
  expect_equal(smt_bad$n_factors[["chi"]], 6)
  expect_equal(smt_bad$n_factors[["AIC"]], 6)

  # a clean sequence stays silent
  expect_no_warning(efa_smt(test_models$baseline$cormat, N = 500),
                    class = "efa_smt_inadmissible")
})

test_that("settings are returned correctly", {
  expect_named(smt_cor$settings, c("N", "use", "cor_method"))
  expect_named(smt_raw$settings, c("N", "use", "cor_method"))
  expect_named(smt_zero$settings, c("N", "use", "cor_method"))

  expect_equal(smt_cor$settings$N, 500)
  expect_equal(smt_raw$settings$N, 810)
  expect_equal(smt_zero$settings$N, 500)

  expect_equal(smt_cor$settings$use, "pairwise.complete.obs")
  expect_equal(smt_raw$settings$use, "pairwise.complete.obs")
  expect_equal(smt_zero$settings$use, "pairwise.complete.obs")

  expect_equal(smt_cor$settings$cor_method, "pearson")
  expect_equal(smt_raw$settings$cor_method, "pearson")
  expect_equal(smt_zero$settings$cor_method, "pearson")

})

# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

test_that("errors are thrown correctly", {
  expect_error(efa_smt(1:5), class = "efa_input_not_matrix")
  expect_error(efa_smt(test_models$baseline$cormat), class = "efa_n_required")
  expect_message(efa_smt(GRiPS_raw), class = "efa_cor_from_data")
  expect_warning(efa_smt(GRiPS_raw, N = 20), class = "efa_n_from_data")
  expect_error(efa_smt(dat_sing), class = "efa_cor_singular")
  expect_error(efa_smt(cor_sing, N = 10), class = "efa_cor_singular")
  expect_error(efa_smt(matrix(rnorm(50), ncol = 2)), class = "efa_smt_underidentified") # underidentified case
  expect_error(efa_smt(matrix(rnorm(60), ncol = 3)), class = "efa_smt_underidentified") # just identified case

  # SMT's sequential tests rest on a normal-theory chi-square that is not valid for
  # polychoric / tetrachoric correlations, so they are rejected rather than run.
  expect_error(efa_smt(GRiPS_raw, cor_method = "poly"),
               class = "efa_cor_method_unsupported")
  expect_error(efa_smt(GRiPS_raw, cor_method = "tetra"),
               class = "efa_cor_method_unsupported")
})

rm(smt_cor, smt_raw, smt_zero, x, y, z, dat_sing, cor_sing)
