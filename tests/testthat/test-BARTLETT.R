bart_cor <- BARTLETT(test_models$baseline$cormat, N = 500)
bart_raw <- BARTLETT(GRiPS_raw)
bart_rand <- BARTLETT(matrix(rnorm(100), ncol = 4))

test_that("output class and dimensions are correct", {
  expect_is(bart_cor, "BARTLETT")
  expect_output(str(bart_cor), "List of 4")
  expect_is(bart_raw, "BARTLETT")
  expect_output(str(bart_raw), "List of 4")
})


test_that("ch-square and p-values are correct", {
  expect_equal(bart_cor$p_value, 0, tolerance = 1e-10)
  expect_equal(bart_cor$p_value, 0, tolerance = 1e-10)
  expect_equal(bart_rand$p_value, 0, tolerance = 1e-10)

})

# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

test_that("errors are thrown correctly", {
  expect_error(BARTLETT(1:5), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.")
  expect_error(BARTLETT(test_models$baseline$cormat), " Argument 'N' was NA, Bartlett's test could not be executed. Please provide either N or raw data.")
  expect_warning(BARTLETT(GRiPS_raw, N = 20), " N was set and data entered. Taking N from data.")
  expect_error(BARTLETT(dat_sing), " Correlation matrix is singular, Bartlett's test cannot be executed.")
  expect_error(BARTLETT(cor_sing, N = 10), " Correlation matrix is singular, Bartlett's test cannot be executed.")
})

# REMOVE OBJECTS
rm(bart_cor, bart_raw, x, y, z, dat_sing, cor_sing)


