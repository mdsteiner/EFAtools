ML_test <- .ML(test_models$baseline$cormat, n_factors = 3, N = 500,
               start_method = "factanal")

test_that("output class and dimensions are correct", {
  expect_is(ML_test$unrot_loadings, "LOADINGS")
  expect_output(str(ML_test), "List of 10")

})

test_that("outputs are correct", {
  expect_equal(ML_test$orig_R, test_models$baseline$cormat)
  expect_equal(sum(ML_test$orig_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(ML_test$final_eigen), ncol(test_models$baseline$cormat))
  expect_equal(ML_test$convergence, 0)
})

test_that("fit indices are returned correctly", {
  expect_output(str(ML_test$fit_indices), "List of 12")

  expect_is(ML_test$fit_indices$chi, "numeric")
  expect_is(ML_test$fit_indices$df, "numeric")
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
})

# HIER STEHENGEBLIEBEN
test_that("settings are returned correctly", {
  expect_named(paf_efatools$settings, c("max_iter", "init_comm", "criterion",
                                        "criterion_type", "abs_eigen"))
  expect_named(paf_psych$settings, c("max_iter", "init_comm", "criterion",
                                     "criterion_type", "abs_eigen"))
  expect_named(paf_spss$settings, c("max_iter", "init_comm", "criterion",
                                    "criterion_type", "abs_eigen"))
  expect_named(paf_none$settings, c("max_iter", "init_comm", "criterion",
                                    "criterion_type", "abs_eigen"))

  expect_equal(paf_efatools$settings$max_iter, 300)
  expect_equal(paf_psych$settings$max_iter, 50)
  expect_equal(paf_spss$settings$max_iter, 25)
  expect_equal(paf_none$settings$max_iter, 500)

  expect_equal(paf_efatools$settings$init_comm, "mac")
  expect_equal(paf_psych$settings$init_comm, "smc")
  expect_equal(paf_spss$settings$init_comm, "smc")
  expect_equal(paf_none$settings$init_comm, "unity")

  expect_equal(paf_efatools$settings$criterion, 0.001)
  expect_equal(paf_psych$settings$criterion, 0.001)
  expect_equal(paf_spss$settings$criterion, 0.001)
  expect_equal(paf_none$settings$criterion, 1e-4)

  expect_equal(paf_efatools$settings$criterion_type, "sums")
  expect_equal(paf_psych$settings$criterion_type, "sums")
  expect_equal(paf_spss$settings$criterion_type, "max_individual")
  expect_equal(paf_none$settings$criterion_type, "sums")

  expect_equal(paf_efatools$settings$abs_eigen, TRUE)
  expect_equal(paf_psych$settings$abs_eigen, FALSE)
  expect_equal(paf_spss$settings$abs_eigen, TRUE)
  expect_equal(paf_none$settings$abs_eigen, TRUE)
})


test_that("warnings and errors are thrown correctly", {
  expect_error(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                    max_iter = 500, type = "none", init_comm = "unity",
                    criterion = 1e-4, criterion_type = "sums"), ' One of "init_comm", "criterion", "criterion_type", "abs_eigen", "max_iter" was NULL and no valid "type" was specified. Either use one of "EFAtools", "psych", or "SPSS" for type, or specify all other arguments')

  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "EFAtools", init_comm = "smc"), " Type and init_comm is specified. init_comm is used with value ' smc '. Results may differ from the specified type")
  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "EFAtools", criterion = 0.001), " Type and criterion is specified. criterion is used with value ' 0.001 '. Results may differ from the specified type")
  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "EFAtools", criterion_type = "sums"), " Type and criterion_type is specified. criterion_type is used with value ' sums '. Results may differ from the specified type")
  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "EFAtools", max_iter = 400), " Type and max_iter is specified. max_iter is used with value ' 400 '. Results may differ from the specified type")
  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "EFAtools", abs_eigen = TRUE), " Type and abs_eigen is specified. abs_eigen is used with value ' TRUE '. Results may differ from the specified type")

  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "psych", init_comm = "smc"), " Type and init_comm is specified. init_comm is used with value ' smc '. Results may differ from the specified type")
  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "psych", criterion = 0.001), " Type and criterion is specified. criterion is used with value ' 0.001 '. Results may differ from the specified type")
  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "psych", criterion_type = "sums"), " Type and criterion_type is specified. criterion_type is used with value ' sums '. Results may differ from the specified type")
  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "psych", max_iter = 400), " Type and max_iter is specified. max_iter is used with value ' 400 '. Results may differ from the specified type")
  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "psych", abs_eigen = TRUE), " Type and abs_eigen is specified. abs_eigen is used with value ' TRUE '. Results may differ from the specified type")

  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "SPSS", init_comm = "smc"), " Type and init_comm is specified. init_comm is used with value ' smc '. Results may differ from the specified type")
  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "SPSS", criterion = 0.001), " Type and criterion is specified. criterion is used with value ' 0.001 '. Results may differ from the specified type")
  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "SPSS", criterion_type = "sums"), " Type and criterion_type is specified. criterion_type is used with value ' sums '. Results may differ from the specified type")
  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "SPSS", max_iter = 400), " Type and max_iter is specified. max_iter is used with value ' 400 '. Results may differ from the specified type")
  expect_warning(.PAF(test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "SPSS", abs_eigen = TRUE), " Type and abs_eigen is specified. abs_eigen is used with value ' TRUE '. Results may differ from the specified type")
})

rm(paf_efatools, paf_psych, paf_spss, paf_none)
