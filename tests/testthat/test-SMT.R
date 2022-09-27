smt_cor <- SMT(test_models$baseline$cormat, N = 500)
smt_zero <- SMT(diag(nrow = 5, ncol = 5), N = 500)
smt_raw <- SMT(GRiPS_raw)

test_that("output class and dimensions are correct", {
  expect_is(smt_cor, "SMT")
  expect_output(str(smt_cor), "List of 10")
  expect_is(smt_raw, "SMT")
  expect_output(str(smt_raw), "List of 10")
  expect_is(smt_zero, "SMT")
  expect_output(str(smt_zero), "List of 10")
})

test_that("number of factors are correct", {
  expect_equal(smt_cor$nfac_chi, 3)
  expect_equal(smt_cor$nfac_RMSEA, 2)
  expect_equal(smt_cor$nfac_AIC, 3)

  expect_equal(smt_raw$nfac_chi, 3)
  expect_equal(smt_raw$nfac_RMSEA, 1)
  expect_equal(smt_raw$nfac_AIC, 3)

  expect_equal(smt_zero$nfac_chi, 0)
  expect_equal(smt_zero$nfac_RMSEA, 0)
  expect_equal(smt_zero$nfac_AIC, 0)
})

test_that("p-values are correct", {
  expect_lt(smt_cor$p_null, 0.05)
  expect_lt(smt_raw$p_null, 0.05)
  expect_gte(smt_zero$p_null, 0.05)

  expect_lt(smt_cor$ps_chi[1],  0.05)
  expect_lt(smt_cor$ps_chi[2],  0.05)
  expect_gte(smt_cor$ps_chi[3],  0.05)
  expect_gte(smt_cor$ps_chi[4],  0.05)

  expect_lt(smt_raw$ps_chi[1],  0.05)
  expect_lt(smt_raw$ps_chi[2],  0.05)
  expect_gte(smt_raw$ps_chi[3],  0.05)
  expect_gte(smt_raw$ps_chi[4],  0.05)

  expect_gte(smt_zero$ps_chi[1], 0.05)
})

test_that("RMSEA_LB and AIC values are correct", {
  expect_equal(smt_cor$RMSEA_LB_null, 0.264456, tolerance = 1e-4)
  expect_equal(smt_raw$RMSEA_LB_null, 0.662764, tolerance = 1e-4)
  expect_equal(smt_zero$RMSEA_LB_null, 0, tolerance = 1e-4)

  expect_equal(smt_cor$RMSEA_LBs, c(0.05674033, 0.03975791, rep(0, 10)),
               tolerance = 1e-4)
  expect_equal(smt_raw$RMSEA_LBs, c(0.03547387, 0.02637614, rep(0, 2)),
               tolerance = 1e-4)
  expect_equal(smt_zero$RMSEA_LBs, rep(0, 2), tolerance = 1e-4)

  expect_equal(smt_cor$AIC_null, 5441.203, tolerance = 0.1)
  expect_equal(smt_raw$AIC_null, 10264.07, tolerance = 0.1)
  expect_equal(smt_zero$AIC_null, -20, tolerance = 0.1)

  expect_equal(smt_cor$AICs, c(139.36142, 17.28713, -78.02180, -77.23762,
                               -74.12983, -66.95662, -56.34764, -45.65766,
                               -35.95682, -25.891317, -16.952594, -5.103373),
               tolerance = 0.1)
  expect_equal(smt_raw$AICs, c(19.919177, 7.833463, -4.447209, -1.852061),
               tolerance = 0.1)
  expect_equal(smt_zero$AICs, c(-10, -2), tolerance = 0.1)
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
  expect_error(SMT(1:5), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n")
  expect_error(SMT(test_models$baseline$cormat), " Argument 'N' was NA. Either provide N or raw data.\n")
  expect_message(SMT(GRiPS_raw), " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n")
  expect_warning(SMT(GRiPS_raw, N = 20), " 'N' was set and data entered. Taking N from data.\n")
  expect_error(SMT(dat_sing), " Correlation matrix is singular, no further analyses are performed\n")
  expect_error(SMT(cor_sing, N = 10), " Correlation matrix is singular, no further analyses are performed\n")
  expect_error(SMT(matrix(rnorm(50), ncol = 2)), " The model is either underidentified or just identified with 1 factor already. SMTs cannot be performed. Please provide more indicators.\n") # underidentified case
  expect_error(SMT(matrix(rnorm(60), ncol = 3)), " The model is either underidentified or just identified with 1 factor already. SMTs cannot be performed. Please provide more indicators.\n") # just identified case
  # expect_warning(SMT(burt, N = 170), "Matrix was not positive definite, smoothing was done")
})

rm(smt_cor, smt_raw, smt_zero, x, y, z, dat_sing, cor_sing, burt)
