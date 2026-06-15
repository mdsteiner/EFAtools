smt_cor <- SMT(test_models$baseline$cormat, N = 500)
smt_zero <- SMT(diag(nrow = 5, ncol = 5), N = 500)
smt_raw <- SMT(GRiPS_raw)

test_that("output class and dimensions are correct", {
  expect_s3_class(smt_cor, "efa_retention")
  expect_length(smt_cor, 7)
  expect_s3_class(smt_raw, "efa_retention")
  expect_length(smt_raw, 7)
  expect_s3_class(smt_zero, "efa_retention")
  expect_length(smt_zero, 7)

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

  expect_equal(rmsea_cor[1], 0.156644, tolerance = 1e-2)
  expect_equal(rmsea_raw[1], 0.460146, tolerance = 1e-2)
  expect_equal(rmsea_zero[1], 0, tolerance = 1e-2)

  expect_equal(rmsea_cor[-1], c(0.0560604, 0.0390112, rep(0, 10)),
               tolerance = 1e-2)
  expect_equal(rmsea_raw[-1], c(0.0352895, 0.0261428, rep(0, 2)),
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
  RMSEA_LB_null <- sqrt(.rmsea_lambda(chi_null, df_null, .95) / (df_null * (N - 1)))
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

burt <- matrix(c(1.00,  0.83,  0.81,  0.80,   0.71, 0.70, 0.54, 0.53,  0.59,  0.24, 0.13,
                 0.83,  1.00,  0.87,  0.62,   0.59, 0.44, 0.58, 0.44,  0.23,  0.45,  0.21,
                 0.81,  0.87,  1.00,  0.63,   0.37, 0.31, 0.30, 0.12,  0.33,  0.33,  0.36,
                 0.80,  0.62,  0.63,  1.00,   0.49, 0.54, 0.30, 0.28,  0.42,  0.29, -0.06,
                 0.71,  0.59,  0.37,  0.49,   1.00, 0.54, 0.34, 0.55,  0.40,  0.19, -0.10,
                 0.70,  0.44,  0.31,  0.54,   0.54, 1.00, 0.50, 0.51,  0.31,  0.11,  0.10,
                 0.54,  0.58,  0.30,  0.30,   0.34, 0.50, 1.00, 0.38,  0.29,  0.21,  0.08,
                 0.53,  0.44,  0.12,  0.28,   0.55, 0.51, 0.38, 1.00,  0.53,  0.10, -0.16,
                 0.59,  0.23,  0.33,  0.42,   0.40, 0.31, 0.29, 0.53,  1.00, -0.09, -0.10,
                 0.24,  0.45,  0.33,  0.29,   0.19, 0.11, 0.21, 0.10, -0.09,  1.00,  0.41,
                 0.13,  0.21,  0.36, -0.06,  -0.10, 0.10, 0.08, -0.16, -0.10, 0.41,  1.00),
               nrow = 11, ncol = 11)

test_that("errors are thrown correctly", {
  expect_error(SMT(1:5), class = "efa_input_not_matrix")
  expect_error(SMT(test_models$baseline$cormat), class = "efa_n_required")
  expect_message(SMT(GRiPS_raw), class = "efa_cor_from_data")
  expect_warning(SMT(GRiPS_raw, N = 20), class = "efa_n_from_data")
  expect_error(SMT(dat_sing), class = "efa_cor_singular")
  expect_error(SMT(cor_sing, N = 10), class = "efa_cor_singular")
  expect_error(SMT(matrix(rnorm(50), ncol = 2)), class = "efa_smt_underidentified") # underidentified case
  expect_error(SMT(matrix(rnorm(60), ncol = 3)), class = "efa_smt_underidentified") # just identified case
  # expect_warning(SMT(burt, N = 170), "Matrix was not positive definite, smoothing was done")
})

rm(smt_cor, smt_raw, smt_zero, x, y, z, dat_sing, cor_sing, burt)
