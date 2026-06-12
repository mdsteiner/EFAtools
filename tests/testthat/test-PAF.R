paf_efatools <- .estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                     max_iter = NA, type = "EFAtools", init_comm = NA,
                     criterion = NA, criterion_type = NA, abs_eigen = NA)
paf_psych <- .estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                     max_iter = NA, type = "psych", init_comm = NA,
                     criterion = NA, criterion_type = NA, abs_eigen = NA)
paf_spss <- .estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                  max_iter = NA, type = "SPSS", init_comm = NA,
                  criterion = NA, criterion_type = NA, abs_eigen = NA)
paf_none <- .estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                 max_iter = 500, type = "none", init_comm = "unity",
                 criterion = 1e-4, criterion_type = "sum", abs_eigen = TRUE)
paf_mac_t <- .estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                 max_iter = 500, type = "none", init_comm = "mac",
                 criterion = 1e-4, criterion_type = "max_individual", abs_eigen = TRUE)
paf_mac_f <- .estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                  max_iter = 500, type = "none", init_comm = "mac",
                  criterion = 1e-4, criterion_type = "max_individual",
                  abs_eigen = FALSE)
paf_F1_t <- .estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 1, N = 500,
                max_iter = 500, type = "none", init_comm = "mac",
                criterion = 1e-4, criterion_type = "max_individual", abs_eigen = TRUE)
paf_F1_f <- .estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 1, N = 500,
               max_iter = 500, type = "none", init_comm = "mac",
               criterion = 1e-4, criterion_type = "max_individual", abs_eigen = FALSE)


test_that("output class and dimensions are correct", {
  expect_s3_class(paf_efatools$unrot_loadings, "LOADINGS")
  expect_s3_class(paf_psych$unrot_loadings, "LOADINGS")
  expect_s3_class(paf_spss$unrot_loadings, "LOADINGS")
  expect_s3_class(paf_none$unrot_loadings, "LOADINGS")
  expect_s3_class(paf_mac_t$unrot_loadings, "LOADINGS")
  expect_s3_class(paf_mac_f$unrot_loadings, "LOADINGS")
  expect_s3_class(paf_F1_t$unrot_loadings, "LOADINGS")
  expect_s3_class(paf_F1_f$unrot_loadings, "LOADINGS")


  expect_output(str(paf_efatools), "List of 15")
  expect_output(str(paf_psych), "List of 15")
  expect_output(str(paf_spss), "List of 15")
  expect_output(str(paf_none), "List of 15")
  expect_output(str(paf_mac_t), "List of 15")
  expect_output(str(paf_mac_f), "List of 15")
  expect_output(str(paf_F1_t), "List of 15")
  expect_output(str(paf_F1_f), "List of 15")
})

test_that("original correlation matrix and eigenvalues are correct", {
  expect_equal(paf_efatools$orig_R, test_models$baseline$cormat)
  expect_equal(paf_psych$orig_R, test_models$baseline$cormat)
  expect_equal(paf_spss$orig_R, test_models$baseline$cormat)
  expect_equal(paf_none$orig_R, test_models$baseline$cormat)
  expect_equal(paf_mac_t$orig_R, test_models$baseline$cormat)
  expect_equal(paf_mac_f$orig_R, test_models$baseline$cormat)
  expect_equal(paf_F1_t$orig_R, test_models$baseline$cormat)
  expect_equal(paf_F1_f$orig_R, test_models$baseline$cormat)

  expect_equal(sum(paf_efatools$orig_eigen), ncol(test_models$baseline$cormat))
  expect_equal(sum(paf_psych$orig_eigen), ncol(test_models$baseline$cormat))
  expect_equal(sum(paf_spss$orig_eigen), ncol(test_models$baseline$cormat))
  expect_equal(sum(paf_none$orig_eigen), ncol(test_models$baseline$cormat))
  expect_equal(sum(paf_mac_t$orig_eigen), ncol(test_models$baseline$cormat))
  expect_equal(sum(paf_mac_f$orig_eigen), ncol(test_models$baseline$cormat))
  expect_equal(sum(paf_F1_t$orig_eigen), ncol(test_models$baseline$cormat))
  expect_equal(sum(paf_F1_f$orig_eigen), ncol(test_models$baseline$cormat))

  expect_lt(sum(paf_efatools$init_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(paf_psych$init_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(paf_spss$init_eigen), ncol(test_models$baseline$cormat))
  expect_equal(sum(paf_none$orig_eigen), ncol(test_models$baseline$cormat))
  expect_equal(sum(paf_mac_t$orig_eigen), ncol(test_models$baseline$cormat))
  expect_equal(sum(paf_mac_f$orig_eigen), ncol(test_models$baseline$cormat))
  expect_equal(sum(paf_F1_t$orig_eigen), ncol(test_models$baseline$cormat))
  expect_equal(sum(paf_F1_f$orig_eigen), ncol(test_models$baseline$cormat))

  expect_lt(sum(paf_efatools$final_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(paf_psych$final_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(paf_spss$final_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(paf_none$final_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(paf_mac_t$final_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(paf_mac_f$final_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(paf_F1_t$final_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(paf_F1_f$final_eigen), ncol(test_models$baseline$cormat))
})

test_that("fit indices are returned correctly", {
  expect_output(str(paf_efatools$fit_indices), "List of 18")
  expect_output(str(paf_psych$fit_indices), "List of 18")
  expect_output(str(paf_spss$fit_indices), "List of 18")
  expect_output(str(paf_none$fit_indices), "List of 18")
  expect_output(str(paf_mac_t$fit_indices), "List of 18")
  expect_output(str(paf_mac_f$fit_indices), "List of 18")
  expect_output(str(paf_F1_t$fit_indices), "List of 18")
  expect_output(str(paf_F1_f$fit_indices), "List of 18")

  expect_equal(paf_efatools$fit_indices[c("chi", "p_chi", "CFI", "RMSEA",
                                          "RMSEA_LB", "RMSEA_UB", "AIC", "BIC",
                                          "Fm", "chi_null", "df_null", "p_null")],
               list(chi = NA, p_chi = NA, CFI = NA, RMSEA = NA, RMSEA_LB = NA,
                    RMSEA_UB = NA, AIC = NA, BIC = NA, Fm = NA, chi_null = NA,
                    df_null = NA, p_null = NA))
  expect_equal(paf_psych$fit_indices[c("chi", "p_chi", "CFI", "RMSEA",
                                       "RMSEA_LB", "RMSEA_UB", "AIC", "BIC",
                                       "Fm", "chi_null", "df_null", "p_null")],
               list(chi = NA, p_chi = NA, CFI = NA, RMSEA = NA, RMSEA_LB = NA,
                    RMSEA_UB = NA, AIC = NA, BIC = NA, Fm = NA, chi_null = NA,
                    df_null = NA, p_null = NA))
  expect_equal(paf_spss$fit_indices[c("chi", "p_chi", "CFI", "RMSEA",
                                      "RMSEA_LB", "RMSEA_UB", "AIC", "BIC",
                                      "Fm", "chi_null", "df_null", "p_null")],
               list(chi = NA, p_chi = NA, CFI = NA, RMSEA = NA, RMSEA_LB = NA,
                    RMSEA_UB = NA, AIC = NA, BIC = NA, Fm = NA, chi_null = NA,
                    df_null = NA, p_null = NA))
  expect_equal(paf_none$fit_indices[c("chi", "p_chi", "CFI", "RMSEA",
                                      "RMSEA_LB", "RMSEA_UB", "AIC", "BIC",
                                      "Fm", "chi_null", "df_null", "p_null")],
               list(chi = NA, p_chi = NA, CFI = NA, RMSEA = NA, RMSEA_LB = NA,
                    RMSEA_UB = NA, AIC = NA, BIC = NA, Fm = NA, chi_null = NA,
                    df_null = NA, p_null = NA))
  expect_equal(paf_mac_t$fit_indices[c("chi", "p_chi", "CFI", "RMSEA",
                                       "RMSEA_LB", "RMSEA_UB", "AIC", "BIC",
                                       "Fm", "chi_null", "df_null", "p_null")],
               list(chi = NA, p_chi = NA, CFI = NA, RMSEA = NA, RMSEA_LB = NA,
                    RMSEA_UB = NA, AIC = NA, BIC = NA, Fm = NA, chi_null = NA,
                    df_null = NA, p_null = NA))
  expect_equal(paf_mac_f$fit_indices[c("chi", "p_chi", "CFI", "RMSEA",
                                       "RMSEA_LB", "RMSEA_UB", "AIC", "BIC",
                                       "Fm", "chi_null", "df_null", "p_null")],
               list(chi = NA, p_chi = NA, CFI = NA, RMSEA = NA, RMSEA_LB = NA,
                    RMSEA_UB = NA, AIC = NA, BIC = NA, Fm = NA, chi_null = NA,
                    df_null = NA, p_null = NA))
  expect_equal(paf_F1_t$fit_indices[c("chi", "p_chi", "CFI", "RMSEA",
                                      "RMSEA_LB", "RMSEA_UB", "AIC", "BIC",
                                      "Fm", "chi_null", "df_null", "p_null")],
               list(chi = NA, p_chi = NA, CFI = NA, RMSEA = NA, RMSEA_LB = NA,
                    RMSEA_UB = NA, AIC = NA, BIC = NA, Fm = NA, chi_null = NA,
                    df_null = NA, p_null = NA))
  expect_equal(paf_F1_f$fit_indices[c("chi", "p_chi", "CFI", "RMSEA",
                                      "RMSEA_LB", "RMSEA_UB", "AIC", "BIC",
                                      "Fm", "chi_null", "df_null", "p_null")],
               list(chi = NA, p_chi = NA, CFI = NA, RMSEA = NA, RMSEA_LB = NA,
                    RMSEA_UB = NA, AIC = NA, BIC = NA, Fm = NA, chi_null = NA,
                    df_null = NA, p_null = NA))

  expect_type(paf_efatools$fit_indices$df, "double")
  expect_type(paf_psych$fit_indices$df, "double")
  expect_type(paf_spss$fit_indices$df, "double")
  expect_type(paf_none$fit_indices$df, "double")
  expect_type(paf_mac_t$fit_indices$df, "double")
  expect_type(paf_mac_f$fit_indices$df, "double")
  expect_type(paf_F1_t$fit_indices$df, "double")
  expect_type(paf_F1_f$fit_indices$df, "double")

  expect_type(paf_efatools$fit_indices$CAF, "double")
  expect_type(paf_psych$fit_indices$CAF, "double")
  expect_type(paf_spss$fit_indices$CAF, "double")
  expect_type(paf_none$fit_indices$CAF, "double")
  expect_type(paf_mac_t$fit_indices$CAF, "double")
  expect_type(paf_mac_f$fit_indices$CAF, "double")
  expect_type(paf_F1_t$fit_indices$CAF, "double")
  expect_type(paf_F1_f$fit_indices$CAF, "double")

})

test_that("settings are returned correctly", {
  expect_named(paf_efatools$settings, c("max_iter", "init_comm", "criterion",
                                        "criterion_type", "abs_eigen"))
  expect_named(paf_psych$settings, c("max_iter", "init_comm", "criterion",
                                        "criterion_type", "abs_eigen"))
  expect_named(paf_spss$settings, c("max_iter", "init_comm", "criterion",
                                        "criterion_type", "abs_eigen"))
  expect_named(paf_none$settings, c("max_iter", "init_comm", "criterion",
                                    "criterion_type", "abs_eigen"))
  expect_named(paf_mac_t$settings, c("max_iter", "init_comm", "criterion",
                                    "criterion_type", "abs_eigen"))
  expect_named(paf_mac_f$settings, c("max_iter", "init_comm", "criterion",
                                   "criterion_type", "abs_eigen"))
  expect_named(paf_F1_t$settings, c("max_iter", "init_comm", "criterion",
                                   "criterion_type", "abs_eigen"))
  expect_named(paf_F1_f$settings, c("max_iter", "init_comm", "criterion",
                                   "criterion_type", "abs_eigen"))

  expect_equal(paf_efatools$settings$max_iter, 300)
  expect_equal(paf_psych$settings$max_iter, 50)
  expect_equal(paf_spss$settings$max_iter, 25)
  expect_equal(paf_none$settings$max_iter, 500)
  expect_equal(paf_mac_t$settings$max_iter, 500)
  expect_equal(paf_mac_f$settings$max_iter, 500)
  expect_equal(paf_F1_t$settings$max_iter, 500)
  expect_equal(paf_F1_f$settings$max_iter, 500)

  expect_equal(paf_efatools$settings$init_comm, "smc")
  expect_equal(paf_psych$settings$init_comm, "smc")
  expect_equal(paf_spss$settings$init_comm, "smc")
  expect_equal(paf_none$settings$init_comm, "unity")
  expect_equal(paf_mac_t$settings$init_comm, "mac")
  expect_equal(paf_mac_f$settings$init_comm, "mac")
  expect_equal(paf_F1_t$settings$init_comm, "mac")
  expect_equal(paf_F1_f$settings$init_comm, "mac")

  expect_equal(paf_efatools$settings$criterion, 0.001)
  expect_equal(paf_psych$settings$criterion, 0.001)
  expect_equal(paf_spss$settings$criterion, 0.001)
  expect_equal(paf_none$settings$criterion, 1e-4)
  expect_equal(paf_mac_t$settings$criterion, 1e-4)
  expect_equal(paf_mac_f$settings$criterion, 1e-4)
  expect_equal(paf_F1_t$settings$criterion, 1e-4)
  expect_equal(paf_F1_f$settings$criterion, 1e-4)

  expect_equal(paf_efatools$settings$criterion_type, "sum")
  expect_equal(paf_psych$settings$criterion_type, "sum")
  expect_equal(paf_spss$settings$criterion_type, "max_individual")
  expect_equal(paf_none$settings$criterion_type, "sum")
  expect_equal(paf_mac_t$settings$criterion_type, "max_individual")
  expect_equal(paf_mac_f$settings$criterion_type, "max_individual")
  expect_equal(paf_F1_t$settings$criterion_type, "max_individual")
  expect_equal(paf_F1_f$settings$criterion_type, "max_individual")

  expect_equal(paf_efatools$settings$abs_eigen, TRUE)
  expect_equal(paf_psych$settings$abs_eigen, FALSE)
  expect_equal(paf_spss$settings$abs_eigen, TRUE)
  expect_equal(paf_none$settings$abs_eigen, TRUE)
  expect_equal(paf_mac_t$settings$abs_eigen, TRUE)
  expect_equal(paf_mac_f$settings$abs_eigen, FALSE)
  expect_equal(paf_F1_t$settings$abs_eigen, TRUE)
  expect_equal(paf_F1_f$settings$abs_eigen, FALSE)
})


test_that("warnings and errors are thrown correctly", {
  expect_error(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                    max_iter = 500, type = "none", init_comm = "unity",
                    criterion = 1e-4, criterion_type = "sum"), class = "efa_type_none")

  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                    type = "EFAtools", init_comm = "smc"), class = "efa_type_override")
  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                    type = "EFAtools", criterion = 0.001), class = "efa_type_override")
  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "EFAtools", criterion_type = "sum"), class = "efa_type_override")
  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "EFAtools", max_iter = 400), class = "efa_type_override")
  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "EFAtools", abs_eigen = TRUE), class = "efa_type_override")

  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "psych", init_comm = "smc"), class = "efa_type_override")
  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "psych", criterion = 0.001), class = "efa_type_override")
  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "psych", criterion_type = "sum"), class = "efa_type_override")
  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "psych", max_iter = 400), class = "efa_type_override")
  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "psych", abs_eigen = TRUE), class = "efa_type_override")

  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "SPSS", init_comm = "smc"), class = "efa_type_override")
  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "SPSS", criterion = 0.001), class = "efa_type_override")
  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "SPSS", criterion_type = "sum"), class = "efa_type_override")
  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "SPSS", max_iter = 400), class = "efa_type_override")
  expect_warning(.estimate_model(method = "PAF",test_models$baseline$cormat, n_factors = 3, N = 500,
                      type = "SPSS", abs_eigen = TRUE), class = "efa_type_override")
  expect_error(.estimate_model(method = "PAF",IDS2_R, n_factors = 7, N = 1991, max_iter = 500, type = "none",
                    init_comm = "smc", criterion = 1e-4, criterion_type = "sum",
                    abs_eigen = FALSE), "Negative Eigenvalues detected; cannot compute communality estimates. Try again with init_comm = 'unity' or 'mac'")
  expect_error(.estimate_model(method = "PAF",IDS2_R, n_factors = 7, N = 1991, max_iter = 500, type = "none",
                    init_comm = "smc", criterion = 1e-4,
                    criterion_type = "max_individual",
                    abs_eigen = FALSE), "Negative Eigenvalues detected; cannot compute communality estimates. Try again with init_comm = 'unity' or 'mac'")
})

rm(paf_efatools, paf_psych, paf_spss, paf_none, paf_mac_t, paf_mac_f, paf_F1_t,
   paf_F1_f)

