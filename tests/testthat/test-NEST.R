
nest_cor <- NEST(test_models$baseline$cormat, N = 500)
nest_raw <- NEST(GRiPS_raw)

test_that("output class and dimensions are correct", {
  expect_is(nest_cor, "NEST")
  expect_output(str(nest_cor), "List of 5")
  expect_is(nest_raw, "NEST")
  expect_output(str(nest_raw), "List of 5")
})


test_that("found eigenvalues are correct", {
  expect_equal(sum(nest_cor$eigenvalues),
               ncol(test_models$baseline$cormat))
  expect_equal(sum(nest_raw$eigenvalues), ncol(GRiPS_raw))
  expect_length(nest_cor$eigenvalues,
                ncol(test_models$baseline$cormat))
  expect_length(nest_raw$eigenvalues, ncol(GRiPS_raw))
})

test_that("reference eigenvalues are correct", {
  expect_length(nest_cor$references, nest_cor$n_factors + 1)
  expect_length(nest_raw$references, nest_raw$n_factors + 1)
  expect_lte(nest_cor$eigenvalues[nest_cor$n_factors + 1],
             nest_cor$references[nest_cor$n_factors + 1])
  expect_true(all(nest_cor$eigenvalues[1:nest_cor$n_factors] >
            nest_cor$references[1:nest_cor$n_factors]))
  expect_lte(nest_raw$eigenvalues[nest_raw$n_factors + 1],
             nest_raw$references[nest_raw$n_factors + 1])
  expect_true(all(nest_raw$eigenvalues[1:nest_raw$n_factors] >
            nest_raw$references[1:nest_raw$n_factors]))
})

test_that("identified number of factors is correct", {
  expect_equal(nest_cor$n_factors, 3)
  expect_equal(nest_raw$n_factors, 1)
})

# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z, rnorm(10), rnorm(10), rnorm(10)), ncol = 6)
cor_sing <- stats::cor(dat_sing)

cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors are thrown correctly", {
  expect_error(NEST(1:5), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n")
  expect_error(NEST(test_models$baseline$cormat), " Argument 'N' was NA but correlation matrix was entered. Please either provide N or raw data.\n")
  expect_message(NEST(GRiPS_raw), " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n")
  expect_warning(NEST(GRiPS_raw, N = 20), " 'N' was set and data entered. Taking N from data.\n")
  expect_error(NEST(dat_sing), " Correlation matrix is singular, no further analyses are performed\n")
  expect_error(NEST(cor_sing, N = 20), " Correlation matrix is singular, no further analyses are performed\n")
})

test_that("settings are returned correctly", {

  expect_equal(nest_cor$settings$N, 500)
  expect_equal(nest_raw$settings$N, 810)

  expect_equal(nest_cor$settings$use, "pairwise.complete.obs")
  expect_equal(nest_raw$settings$use, "pairwise.complete.obs")

  expect_equal(nest_cor$settings$cor_method, "pearson")
  expect_equal(nest_raw$settings$cor_method, "pearson")
})

rm(nest_cor, nest_raw, x, y, z, dat_sing, cor_sing, cor_nposdef)
