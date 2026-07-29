# CD is stochastic; seed and use a small simulation for a stable, fast test
set.seed(123)
cd_grips <- efa_cd(GRiPS_raw, N_pop = 1000, N_samples = 100)

test_that("output class and dimensions are correct", {
  expect_s3_class(cd_grips, "efa_retention")
  expect_length(cd_grips, 6)
  expect_named(cd_grips$n_factors, "CD")

  rec <- .retention_record(cd_grips, "CD")
  expect_equal(rec$plot_type, "eigen")
  expect_equal(rec$y_label, "RMSE eigenvalues")
  checkmate::expect_matrix(rec$rmse_eigenvalues)
})

test_that("CD returns the correct values", {
  expect_equal(cd_grips$n_factors[["CD"]], 1)

  rec <- .retention_record(cd_grips, "CD")
  # the RMSE curve runs one factor beyond the retained number (the first count
  # whose lack of significant improvement stopped the search)
  expect_length(rec$y, cd_grips$n_factors[["CD"]] + 1)
  expect_length(rec$x, cd_grips$n_factors[["CD"]] + 1)
})

test_that("the comparison-data RMSE curve is unchanged (regression)", {
  # Regression pin on the comparison-data draw: CD generates its populations
  # through the shared Ruscio-Kaczetow kernel, so pin the mean RMSE-eigenvalue
  # curve of the seeded run to catch any drift.
  #
  # The kernel is a stochastic search whose accept-or-shrink step compares two
  # floating-point RMSEs. Once the search plateaus that comparison is decided by
  # rounding, so the trajectory -- and with it this curve -- is reproducible only
  # within one BLAS; the retained factor count the method exists to report is
  # stable regardless and is asserted above. The pin therefore runs where the
  # BLAS is known rather than across every check flavour.
  skip_on_cran()
  rec <- .retention_record(cd_grips, "CD")
  # Only the tested factor counts (1 and 2) are populated; later columns stay
  # zero because the search stops after the first factor.
  expect_equal(
    colMeans(rec$rmse_eigenvalues),
    c(0.033821350175717364, 0.038807011125946637, 0, 0),
    tolerance = 1e-6
  )
})

grips_na <- GRiPS_raw
grips_na[c(1,3,5), c(2, 4)] <- NA
test_that("errors etc. are thrown correctly", {
  expect_error(efa_cd(1:10), class = "efa_input_not_matrix")
  expect_error(efa_cd(test_models$baseline$cormat), class = "efa_cd_needs_raw")

  expect_warning(efa_cd(GRiPS_raw, n_factors_max = 5, N_pop = 500, N_samples = 20), class = "efa_cd_max_factors")
  expect_warning(efa_cd(grips_na, n_factors_max = 3, N_pop = 500, N_samples = 20), class = "efa_cd_missing_removed")

  # With three indicators the one-factor model is just identified (df = 0), so no
  # over-identified model exists to compare against and CD must abort rather than
  # silently return zero factors. The check runs before any data are generated.
  expect_error(efa_cd(GRiPS_raw[, 1:3]), class = "efa_cd_min_indicators")
})

rm(cd_grips, grips_na)
