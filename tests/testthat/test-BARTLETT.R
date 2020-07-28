bart_cor <- BARTLETT(test_models$baseline$cormat, N = 500)
bart_raw <- BARTLETT(GRiPS_raw)

set.seed(500)
bart_rand <- BARTLETT(matrix(rnorm(100), ncol = 4))

test_that("output class and dimensions are correct", {
  expect_is(bart_cor, "BARTLETT")
  expect_output(str(bart_cor), "List of 4")
  expect_is(bart_raw, "BARTLETT")
  expect_output(str(bart_raw), "List of 4")
})


test_that("p-values and df are correct", {
  expect_lt(bart_cor$p_value, 0.0001)
  expect_lt(bart_raw$p_value, 0.0001)
  expect_gt(bart_rand$p_value, 0.05)

  expect_equal(bart_cor$df, 153)
  expect_equal(bart_raw$df, 28)
  expect_equal(bart_rand$df, 6)
})

test_that("settings are returned correctly", {
  expect_named(bart_cor$settings, c("N", "use", "cor_method"))
  expect_named(bart_raw$settings, c("N", "use", "cor_method"))
  expect_named(bart_rand$settings, c("N", "use", "cor_method"))

  expect_equal(bart_cor$settings$N, 500)
  expect_equal(bart_raw$settings$N, 810)
  expect_equal(bart_rand$settings$N, 25)

  expect_equal(bart_cor$settings$use, "pairwise.complete.obs")
  expect_equal(bart_raw$settings$use, "pairwise.complete.obs")
  expect_equal(bart_rand$settings$use, "pairwise.complete.obs")

  expect_equal(bart_cor$settings$cor_method, "pearson")
  expect_equal(bart_raw$settings$cor_method, "pearson")
  expect_equal(bart_rand$settings$cor_method, "pearson")
})

# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y

dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors are thrown correctly", {
  expect_error(BARTLETT(1:5), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n")
  expect_error(BARTLETT(test_models$baseline$cormat), " Argument 'N' was NA, Bartlett's test could not be executed. Please provide either N or raw data.\n")
  expect_message(BARTLETT(GRiPS_raw), " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n")
  expect_warning(BARTLETT(GRiPS_raw, N = 20), " 'N' was set and data entered. Taking N from data.\n")
  expect_error(BARTLETT(dat_sing), " Correlation matrix is singular, Bartlett's test cannot be executed.\n")
  expect_error(BARTLETT(cor_sing, N = 10), " Correlation matrix is singular, Bartlett's test cannot be executed.\n")
  expect_warning(BARTLETT(cor_nposdef, N = 10), "Matrix was not positive definite, smoothing was done")
})

rm(bart_cor, bart_raw, bart_rand, x, y, z, dat_sing, cor_sing, cor_nposdef)


