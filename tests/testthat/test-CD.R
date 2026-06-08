cd_grips <- CD(GRiPS_raw)

test_that("output class and dimensions are correct", {
  expect_s3_class(cd_grips, "CD")
  expect_named(cd_grips, c("n_factors", "eigenvalues", "RMSE_eigenvalues",
                           "settings"))
  checkmate::expect_matrix(cd_grips$RMSE_eigenvalues)
})

test_that("CD returns the correct values", {
  expect_equal(cd_grips$n_factors, 1)
  expect_equal(sum(cd_grips$eigenvalues), 8)
})

grips_na <- GRiPS_raw
grips_na[c(1,3,5), c(2, 4)] <- NA
test_that("errors etc. are thrown correctly", {
  expect_error(CD(1:10), class = "efa_input_not_matrix")
  expect_error(CD(test_models$baseline$cormat), class = "efa_cd_needs_raw")

  expect_warning(CD(GRiPS_raw, n_factors_max = 5), class = "efa_cd_max_factors")
  expect_warning(CD(grips_na, n_factors_max = 3), class = "efa_cd_missing_removed")
})

rm(cd_grips, grips_na)
