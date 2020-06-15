smt_cor <- SMT(test_models$baseline$cormat, N = 500)
smt_raw <- SMT(GRiPS_raw)

set.seed(500)
smt_rand <- SMT(matrix(rnorm(200), ncol = 8))

test_that("output class and dimensions are correct", {
  expect_is(smt_cor, "SMT")
  expect_output(str(smt_cor), "List of 8")
  expect_is(smt_raw, "SMT")
  expect_output(str(smt_raw), "List of 8")
})

test_that("number of factors are correct", {
  expect_equal(smt_cor$nfac_chi, 3)
  expect_equal(smt_cor$nfac_RMSEA, 2)
  expect_equal(smt_cor$nfac_AIC, 3)

  expect_equal(smt_raw$nfac_chi, 3)
  expect_equal(smt_raw$nfac_RMSEA, 1)
  expect_equal(smt_raw$nfac_AIC, 3)

  expect_equal(smt_rand$nfac_chi, 0)
  expect_equal(smt_rand$nfac_RMSEA, NA)
  expect_equal(smt_rand$nfac_AIC, NA)
})

test_that("p-values are correct", {
  expect_lt(smt_cor$p_null, 0.05)
  expect_lt(smt_raw$p_null, 0.05)
  expect_gte(smt_rand$p_null, 0.05)

  expect_lt(smt_cor$ps_chi[1],  0.05)
  expect_lt(smt_cor$ps_chi[2],  0.05)
  expect_gte(smt_cor$ps_chi[3],  0.05)
  expect_gte(smt_cor$ps_chi[4],  0.05)

  expect_lt(smt_raw$ps_chi[1],  0.05)
  expect_lt(smt_raw$ps_chi[2],  0.05)
  expect_gte(smt_raw$ps_chi[3],  0.05)
  expect_gte(smt_raw$ps_chi[4],  0.05)

  expect_gte(smt_rand$ps_chi[1], 0.05)
})


test_that("settings are returned correctly", {
  expect_named(smt_cor$settings, c("N", "use", "cor_method"))
  expect_named(smt_raw$settings, c("N", "use", "cor_method"))
  expect_named(smt_rand$settings, c("N", "use", "cor_method"))

  expect_equal(smt_cor$settings$N, 500)
  expect_equal(smt_raw$settings$N, 810)
  expect_equal(smt_rand$settings$N, 25)

  expect_equal(smt_cor$settings$use, "pairwise.complete.obs")
  expect_equal(smt_raw$settings$use, "pairwise.complete.obs")
  expect_equal(smt_rand$settings$use, "pairwise.complete.obs")

  expect_equal(smt_cor$settings$cor_method, "pearson")
  expect_equal(smt_raw$settings$cor_method, "pearson")
  expect_equal(smt_rand$settings$cor_method, "pearson")

})

# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

test_that("errors are thrown correctly", {
  expect_error(SMT(1:5), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.")
  expect_error(SMT(test_models$baseline$cormat), " Argument 'N' was NA. Either provide N or raw data.")
  expect_message(SMT(GRiPS_raw), " 'x' was not a correlation matrix. Correlations are found from entered raw data.")
  expect_warning(SMT(GRiPS_raw, N = 20), " 'N' was set and data entered. Taking N from data.")
  expect_error(SMT(dat_sing), " Correlation matrix is singular, no further analyses are performed")
  expect_error(SMT(cor_sing, N = 10), " Correlation matrix is singular, no further analyses are performed")
  expect_warning(SMT(matrix(rnorm(100), ncol = 4)), " The maximum number of factors tested was 1. RMSEA lower bound and AIC criteria cannot be applied. nfac_RMSEA and nfac_AIC were set to NA.")
})

rm(smt_cor, smt_raw, smt_rand, x, y, z, dat_sing, cor_sing)
