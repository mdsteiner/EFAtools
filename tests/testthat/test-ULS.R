ULS_test <- .ULS(test_models$baseline$cormat, n_factors = 3, N = 500)
ULS_test_1 <- .ULS(test_models$baseline$cormat, n_factors = 1, N = 500)

test_that("output class and dimensions are correct", {
  expect_s3_class(ULS_test$unrot_loadings, "LOADINGS")
  expect_output(str(ULS_test), "List of 11")
  expect_s3_class(ULS_test_1$unrot_loadings, "LOADINGS")
  expect_output(str(ULS_test_1), "List of 11")

})

test_that("outputs are correct", {
  expect_equal(ULS_test$orig_R, test_models$baseline$cormat)
  expect_equal(sum(ULS_test$orig_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(ULS_test$final_eigen), ncol(test_models$baseline$cormat))
  expect_equal(ULS_test$convergence, 0)

  expect_equal(ULS_test_1$orig_R, test_models$baseline$cormat)
  expect_equal(sum(ULS_test_1$orig_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(ULS_test_1$final_eigen), ncol(test_models$baseline$cormat))
  expect_equal(ULS_test_1$convergence, 0)
})

test_that("fit indices are returned correctly", {
  expect_output(str(ULS_test$fit_indices), "List of 15")
  expect_output(str(ULS_test_1$fit_indices), "List of 15")

  expect_type(ULS_test$fit_indices$chi, "double")
  expect_type(ULS_test$fit_indices$df, "double")
  expect_type(ULS_test$fit_indices$p_chi, "double")
  expect_type(ULS_test$fit_indices$CAF, "double")
  expect_type(ULS_test$fit_indices$CFI, "double")
  expect_type(ULS_test$fit_indices$RMSEA, "double")
  expect_type(ULS_test$fit_indices$RMSEA_LB, "double")
  expect_type(ULS_test$fit_indices$RMSEA_UB, "double")
  expect_type(ULS_test$fit_indices$AIC, "double")
  expect_type(ULS_test$fit_indices$BIC, "double")
  expect_type(ULS_test$fit_indices$Fm, "double")
  expect_type(ULS_test$fit_indices$chi_null, "double")
  expect_type(ULS_test$fit_indices$df_null, "double")
  expect_type(ULS_test$fit_indices$p_null, "double")

  expect_type(ULS_test_1$fit_indices$chi, "double")
  expect_type(ULS_test_1$fit_indices$df, "double")
  expect_type(ULS_test_1$fit_indices$p_chi, "double")
  expect_type(ULS_test_1$fit_indices$CAF, "double")
  expect_type(ULS_test_1$fit_indices$CFI, "double")
  expect_type(ULS_test_1$fit_indices$RMSEA, "double")
  expect_type(ULS_test_1$fit_indices$RMSEA_LB, "double")
  expect_type(ULS_test_1$fit_indices$RMSEA_UB, "double")
  expect_type(ULS_test_1$fit_indices$AIC, "double")
  expect_type(ULS_test_1$fit_indices$BIC, "double")
  expect_type(ULS_test_1$fit_indices$Fm, "double")
  expect_type(ULS_test_1$fit_indices$chi_null, "double")
  expect_type(ULS_test_1$fit_indices$df_null, "double")
  expect_type(ULS_test_1$fit_indices$p_null, "double")
})

rm(ULS_test, ULS_test_1)
