ULS_test <- .estimate_model(method = "ULS",test_models$baseline$cormat, n_factors = 3, N = 500)
ULS_test_1 <- .estimate_model(method = "ULS",test_models$baseline$cormat, n_factors = 1, N = 500)

test_that("output class and dimensions are correct", {
  expect_s3_class(ULS_test$unrot_loadings, "LOADINGS")
  expect_output(str(ULS_test), "List of 12")
  expect_s3_class(ULS_test_1$unrot_loadings, "LOADINGS")
  expect_output(str(ULS_test_1), "List of 12")

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
  expect_output(str(ULS_test$fit_indices), "List of 18")
  expect_output(str(ULS_test_1$fit_indices), "List of 18")

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

test_that("MINRES is accepted as a synonym for ULS", {
  uls <- suppressWarnings(EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                              method = "ULS"))
  minres <- suppressWarnings(EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                                 method = "MINRES"))

  # minimum residual and unweighted least squares are the same estimator
  expect_equal(minres$unrot_loadings, uls$unrot_loadings)
  expect_equal(minres$h2, uls$h2)
  expect_equal(minres$fit_indices, uls$fit_indices)
  # the alias resolves to the canonical method name
  expect_identical(minres$settings$method, "ULS")
})

test_that("SL() and EFA_AVERAGE() accept method = 'MINRES'", {
  EFA_mod <- suppressWarnings(EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
                                  method = "PAF", rotation = "promax"))
  expect_no_error(suppressWarnings(SL(EFA_mod, method = "MINRES")))
  expect_no_error(suppressWarnings(
    EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                method = "MINRES", show_progress = FALSE)))
})

rm(ULS_test, ULS_test_1)
