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

# Example from nearPD function from Matrix package
cor_nposdef <- matrix(c(1,     0.477, 0.644, 0.578, 0.651, 0.826,
                        0.477, 1,     0.516, 0.233, 0.682, 0.75,
                        0.644, 0.516, 1,     0.599, 0.581, 0.742,
                        0.478, 0.233, 0.599, 1,     0.741, 0.8,
                        0.651, 0.682, 0.581, 0.741, 1,     0.798,
                        0.826, 0.75,  0.742, 0.8,   0.798, 1),
                      nrow = 6, ncol = 6)

test_that("errors are thrown correctly", {
  expect_error(SMT(1:5), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n")
  expect_error(SMT(test_models$baseline$cormat), " Argument 'N' was NA. Either provide N or raw data.\n")
  expect_message(SMT(GRiPS_raw), " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n")
  expect_warning(SMT(GRiPS_raw, N = 20), " 'N' was set and data entered. Taking N from data.\n")
  expect_error(SMT(dat_sing), " Correlation matrix is singular, no further analyses are performed\n")
  expect_error(SMT(cor_sing, N = 10), " Correlation matrix is singular, no further analyses are performed\n")
  expect_error(SMT(matrix(rnorm(50), ncol = 2)), " The model is either underidentified or just identified with 1 factor already. SMTs cannot be performed. Please provide more indicators.\n") # underidentified case
  expect_error(SMT(matrix(rnorm(60), ncol = 3)), " The model is either underidentified or just identified with 1 factor already. SMTs cannot be performed. Please provide more indicators.\n") # just identified case
  expect_warning(SMT(cor_nposdef, N = 10), "Matrix was not positive definite, smoothing was done")
})

rm(smt_cor, smt_raw, smt_zero, x, y, z, dat_sing, cor_sing, cor_nposdef)
