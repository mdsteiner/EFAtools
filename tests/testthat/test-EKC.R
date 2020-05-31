
ekc_cor <- EKC(test_models$baseline$cormat, N = 500)
ekc_raw <- EKC(GRiPS_raw)

test_that("output class and dimensions are correct", {
  expect_is(ekc_cor, "EKC")
  expect_output(str(ekc_cor), "List of 4")
  expect_is(ekc_raw, "EKC")
  expect_output(str(ekc_raw), "List of 4")
})


test_that("found eigenvalues are correct", {
  expect_equal(sum(ekc_cor$eigenvalues),
               ncol(test_models$baseline$cormat))
  expect_equal(sum(ekc_raw$eigenvalues), ncol(GRiPS_raw))
  expect_length(ekc_cor$eigenvalues,
                ncol(test_models$baseline$cormat))
  expect_length(ekc_raw$eigenvalues, ncol(GRiPS_raw))
})

test_that("reference eigenvalues are correct", {
  expect_gt(ekc_cor$references[floor(ncol(test_models$baseline$cormat) / 2)], 1)
  expect_gt(ekc_raw$references[floor(ncol(GRiPS_raw) / 2)], 1)
  expect_length(ekc_cor$references,
                length(ekc_cor$eigenvalues))
  expect_length(ekc_raw$references, length(ekc_raw$eigenvalues))
})

test_that("identified number of factors is correct", {
  expect_equal(ekc_cor$n_factors, 2)
  expect_equal(ekc_raw$n_factors, 1)
})

test_that("message on input of raw data", {
  expect_message(EKC(GRiPS_raw), "x was not a correlation matrix. Correlations and N are found from entered raw data.")
})
