cd_grips <- CD(GRiPS_raw)

test_that("output class and dimensions are correct", {
  expect_is(cd_grips, "CD")
  expect_named(cd_grips, c("n_factors", "eigenvalues", "RMSE_eigenvalues",
                           "settings"))
  expect_is(cd_grips$RMSE_eigenvalues, "matrix")
})

test_that("CD returns the correct values", {
  expect_equal(cd_grips$n_factors, 1)
  expect_equal(sum(cd_grips$eigenvalues), 8)
})

test_that("errors etc. are thrown correctly", {
  expect_error(CD(1:10), " 'x' is neither a matrix nor a dataframe. Provide a dataframe or matrix with raw data.\n")
  expect_error(CD(test_models$baseline$cormat), " 'x' is a correlation matrix, but CD only works with raw data.\n")

  expect_warning(CD(GRiPS_raw, n_factors_max = 5), " n_factors_max was set to 5 but maximum possible factors to extract is 4 . Setting n_factors_max to 4 .\n")
})

rm(cd_grips)
