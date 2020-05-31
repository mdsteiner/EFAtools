
pa_cor <- PARALLEL(test_models$baseline$cormat, N = 500)
pa_cor_pca <- PARALLEL(test_models$baseline$cormat, N = 500, eigen_type = "PCA")
pa_raw <- PARALLEL(GRiPS_raw, N = 500)

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
  expect_equal(sum(ekc_raw$eigenvalues), ncol(GRiPS_raw))
  expect_length(ekc_cor$eigenvalues,
                ncol(test_models$baseline$cormat))
  expect_length(ekc_raw$eigenvalues, ncol(GRiPS_raw))
})

