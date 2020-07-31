ML_test <- .ML(test_models$baseline$cormat, n_factors = 3, N = 500,
               start_method = "factanal")
ML_test_2 <- .ML(test_models$baseline$cormat, n_factors = 3, N = 500,
               start_method = "psych")

test_that("output class and dimensions are correct", {
  expect_is(ML_test$unrot_loadings, "LOADINGS")
  expect_output(str(ML_test), "List of 10")
  expect_is(ML_test_2$unrot_loadings, "LOADINGS")
  expect_output(str(ML_test_2), "List of 10")

})

test_that("outputs are correct", {
  expect_equal(ML_test$orig_R, test_models$baseline$cormat)
  expect_equal(sum(ML_test$orig_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(ML_test$final_eigen), ncol(test_models$baseline$cormat))
  expect_equal(ML_test$convergence, 0)
  expect_equal(ML_test_2$orig_R, test_models$baseline$cormat)
  expect_equal(sum(ML_test_2$orig_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(ML_test_2$final_eigen), ncol(test_models$baseline$cormat))
  expect_equal(ML_test_2$convergence, 0)
})

test_that("fit indices are returned correctly", {
  expect_output(str(ML_test$fit_indices), "List of 14")

  expect_is(ML_test$fit_indices$chi, "numeric")
  expect_is(ML_test$fit_indices$df, "numeric")
  expect_is(ML_test$fit_indices$p_chi, "numeric")
  expect_is(ML_test$fit_indices$CAF, "numeric")
  expect_is(ML_test$fit_indices$CFI, "numeric")
  expect_is(ML_test$fit_indices$RMSEA, "numeric")
  expect_is(ML_test$fit_indices$RMSEA_LB, "numeric")
  expect_is(ML_test$fit_indices$RMSEA_UB, "numeric")
  expect_is(ML_test$fit_indices$AIC, "numeric")
  expect_is(ML_test$fit_indices$BIC, "numeric")
  expect_is(ML_test$fit_indices$Fm, "numeric")
  expect_is(ML_test$fit_indices$chi_null, "numeric")
  expect_is(ML_test$fit_indices$df_null, "numeric")
  expect_is(ML_test$fit_indices$p_null, "numeric")

  expect_is(ML_test_2$fit_indices$chi, "numeric")
  expect_is(ML_test_2$fit_indices$df, "numeric")
  expect_is(ML_test_2$fit_indices$p_chi, "numeric")
  expect_is(ML_test_2$fit_indices$CAF, "numeric")
  expect_is(ML_test_2$fit_indices$CFI, "numeric")
  expect_is(ML_test_2$fit_indices$RMSEA, "numeric")
  expect_is(ML_test_2$fit_indices$RMSEA_LB, "numeric")
  expect_is(ML_test_2$fit_indices$RMSEA_UB, "numeric")
  expect_is(ML_test_2$fit_indices$AIC, "numeric")
  expect_is(ML_test_2$fit_indices$BIC, "numeric")
  expect_is(ML_test_2$fit_indices$Fm, "numeric")
  expect_is(ML_test_2$fit_indices$chi_null, "numeric")
  expect_is(ML_test_2$fit_indices$df_null, "numeric")
  expect_is(ML_test_2$fit_indices$p_null, "numeric")
})

test_that("settings are returned correctly", {
  expect_named(ML_test$settings, c("start_method"))
  expect_equal(ML_test$settings$start_method, "factanal")
  expect_named(ML_test_2$settings, c("start_method"))
  expect_equal(ML_test_2$settings$start_method, "psych")
})

rm(ML_test, ML_test_2)
