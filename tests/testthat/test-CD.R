# CD is stochastic; seed and use a small simulation for a stable, fast test
set.seed(123)
cd_grips <- CD(GRiPS_raw, N_pop = 1000, N_samples = 100)

test_that("output class and dimensions are correct", {
  expect_s3_class(cd_grips, "efa_retention")
  expect_length(cd_grips, 7)
  expect_named(cd_grips$n_factors, "CD")

  rec <- .retention_record(cd_grips, "CD")
  expect_equal(rec$plot_type, "eigen")
  expect_equal(rec$y_label, "RMSE eigenvalues")
  checkmate::expect_matrix(rec$rmse_eigenvalues)
})

test_that("CD returns the correct values", {
  expect_equal(cd_grips$n_factors[["CD"]], 1)

  rec <- .retention_record(cd_grips, "CD")
  expect_length(rec$y, cd_grips$n_factors[["CD"]])
  expect_length(rec$x, cd_grips$n_factors[["CD"]])
})

grips_na <- GRiPS_raw
grips_na[c(1,3,5), c(2, 4)] <- NA
test_that("errors etc. are thrown correctly", {
  expect_error(CD(1:10), class = "efa_input_not_matrix")
  expect_error(CD(test_models$baseline$cormat), class = "efa_cd_needs_raw")

  expect_warning(CD(GRiPS_raw, n_factors_max = 5, N_pop = 500, N_samples = 20), class = "efa_cd_max_factors")
  expect_warning(CD(grips_na, n_factors_max = 3, N_pop = 500, N_samples = 20), class = "efa_cd_missing_removed")
})

rm(cd_grips, grips_na)
