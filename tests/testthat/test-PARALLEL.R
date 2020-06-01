
pa_cor <- PARALLEL(test_models$baseline$cormat, N = 500)
pa_cor_pca <- PARALLEL(test_models$baseline$cormat, N = 500, eigen_type = "PCA")
pa_raw <- PARALLEL(GRiPS_raw)

test_that("output class and dimensions are correct", {
  expect_is(pa_cor, "PARALLEL")
  expect_output(str(pa_cor), "List of 7")
  expect_is(pa_raw, "PARALLEL")
  expect_output(str(pa_raw), "List of 7")
})


test_that("found eigenvalues are correct", {
  expect_equal(sum(pa_cor$eigenvalues_PCA[, 1]),
               ncol(test_models$baseline$cormat))
  expect_equal(sum(pa_cor$eigenvalues_PCA[, 2]),
               ncol(test_models$baseline$cormat))
  expect_gt(sum(pa_cor$eigenvalues_PCA[, 3]),
               ncol(test_models$baseline$cormat))
  expect_lt(sum(pa_cor$eigenvalues_SMC[, 1]),
            ncol(test_models$baseline$cormat))
  expect_lt(sum(pa_cor$eigenvalues_SMC[, 2]),
            ncol(test_models$baseline$cormat))
  expect_lt(sum(pa_cor$eigenvalues_EFA[, 1]),
               ncol(test_models$baseline$cormat))
  expect_lt(sum(pa_cor$eigenvalues_EFA[, 2]),
               ncol(test_models$baseline$cormat))

  expect_equal(sum(pa_raw$eigenvalues_PCA[, 1]),
               ncol(GRiPS_raw))
  expect_equal(sum(pa_raw$eigenvalues_PCA[, 2]),
               ncol(GRiPS_raw))
  expect_gt(sum(pa_raw$eigenvalues_PCA[, 3]),
            ncol(GRiPS_raw))
  expect_lt(sum(pa_raw$eigenvalues_SMC[, 1]),
            ncol(GRiPS_raw))
  expect_lt(sum(pa_raw$eigenvalues_SMC[, 2]),
            ncol(GRiPS_raw))
  expect_lt(sum(pa_raw$eigenvalues_EFA[, 1]),
            ncol(GRiPS_raw))
  expect_lt(sum(pa_raw$eigenvalues_EFA[, 2]),
            ncol(GRiPS_raw))

  expect_is(pa_cor$eigenvalues_PCA, "matrix")
  expect_is(pa_cor$eigenvalues_SMC, "matrix")
  expect_is(pa_cor$eigenvalues_EFA, "matrix")

  expect_is(pa_raw$eigenvalues_PCA, "matrix")
  expect_is(pa_raw$eigenvalues_SMC, "matrix")
  expect_is(pa_raw$eigenvalues_EFA, "matrix")

  expect_is(pa_cor_pca$eigenvalues_PCA, "matrix")
  expect_equal(c(pa_cor_pca$eigenvalues_SMC, pa_cor_pca$eigenvalues_EFA), c(NA, NA))

})


test_that("identified number of factors is correct", {
  expect_equal(pa_cor$n_fac_PCA, 3)
  expect_equal(pa_cor$n_fac_SMC, 3)
  expect_equal(pa_cor$n_fac_EFA, 7, tolerance = 1)

  expect_equal(pa_raw$n_fac_PCA, 1)
  expect_equal(pa_raw$n_fac_SMC, 1)
  expect_equal(pa_raw$n_fac_EFA, 3)

  expect_equal(pa_cor_pca$n_fac_PCA, 3)
  expect_equal(c(pa_cor_pca$n_fac_SMC, pa_cor_pca$n_fac_EFA), c(NA, NA))
})

test_that("errors are thrown correctly", {
  expect_message(PARALLEL(GRiPS_raw, eigen_type = "PCA"), "x was not a correlation matrix. Correlations are found from entered raw data.")
  expect_warning(PARALLEL(GRiPS_raw, N = 20, eigen_type = "PCA"), "N was set and data entered. Taking N from data.")
  expect_error(PARALLEL(test_models$baseline$cormat, eigen_type = "PCA"), '"N" was not set and could not be taken from data. Please specify N and try again.')
  expect_warning(PARALLEL(test_models$baseline$cormat, N = 500,
                          eigen_type = "PCA", decision_rule = "Crawford",
                          percent = 80), "decision_rule == 'Crawford' is specified, but 95 percentile was not used. Using Means instead. To use 'Crawford', make sure to specify percent = 95.")
})

rm(pa_cor, pa_cor_pca, pa_raw, na_cor)
