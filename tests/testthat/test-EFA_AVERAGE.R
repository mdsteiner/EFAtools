efa_def <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500)
efa_ml <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                      method = "ML")
efa_uls <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                      method = "ULS")

efa_all <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                       method = c("PAF", "ML", "ULS"),
                       type = c("none", "EFAtools", "psych", "SPSS"),
                       salience_threshold = .2)
efa_all_np <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                       method = c("PAF", "ML", "ULS"),
                       type = c("none", "EFAtools", "psych", "SPSS"),
                       salience_threshold = .2, show_progress = FALSE)

efa_all_oblq <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                            method = c("PAF", "ML", "ULS"),
                            type = c("none", "EFAtools", "psych", "SPSS"),
                            rotation = "oblique")
efa_all_orth <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                            method = c("PAF", "ML", "ULS"),
                            type = c("none", "EFAtools", "psych", "SPSS"),
                            rotation = "orthogonal")
efa_all_none <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                            method = c("PAF", "ML", "ULS"),
                            type = c("none", "EFAtools", "psych", "SPSS"),
                            rotation = "none")

efa_all_md <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                            method = c("PAF", "ML", "ULS"),
                            type = c("none", "EFAtools", "psych", "SPSS"),
                            rotation = "oblique", averaging = "median")
efa_all_tm <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                          method = c("PAF", "ML", "ULS"),
                          type = c("none", "EFAtools", "psych", "SPSS"),
                          rotation = "oblique", averaging = "mean",
                          trim = .2)
efa_raw <- EFA_AVERAGE(GRiPS_raw, n_factors = 1, rotation = "none")
efa_raw_p <- EFA_AVERAGE(GRiPS_raw, n_factors = 2, rotation = "promax")

test_that("output class and dimensions are correct", {
  expect_is(efa_def, "EFA_AVERAGE")
  expect_is(efa_ml, "EFA_AVERAGE")
  expect_is(efa_uls, "EFA_AVERAGE")
  expect_is(efa_all, "EFA_AVERAGE")
  expect_is(efa_all_oblq, "EFA_AVERAGE")
  expect_is(efa_all_orth, "EFA_AVERAGE")
  expect_is(efa_all_none, "EFA_AVERAGE")
  expect_is(efa_all_md, "EFA_AVERAGE")
  expect_is(efa_all_tm, "EFA_AVERAGE")
  expect_is(efa_raw, "EFA_AVERAGE")
  expect_is(efa_raw_p, "EFA_AVERAGE")

  expect_is(efa_def$loadings$average, "LOADINGS")
  expect_is(efa_ml$loadings$average, "LOADINGS")
  expect_is(efa_uls$loadings$average, "LOADINGS")
  expect_is(efa_all$loadings$average, "LOADINGS")
  expect_is(efa_all_oblq$loadings$average, "LOADINGS")
  expect_is(efa_all_orth$loadings$average, "LOADINGS")
  expect_is(efa_all_none$loadings$average, "LOADINGS")
  expect_is(efa_all_md$loadings$average, "LOADINGS")
  expect_is(efa_all_tm$loadings$average, "LOADINGS")
  expect_is(efa_raw$loadings$average, "LOADINGS")
  expect_is(efa_raw_p$loadings$average, "LOADINGS")

  expect_named(efa_def, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_ml, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_uls, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_all, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_all_oblq, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_all_orth, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_all_none, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_all_md, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_all_tm, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_raw, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
  expect_named(efa_raw_p, c("orig_R", "h2", "loadings", "Phi", "ind_fac_corres",
                          "vars_accounted", "fit_indices", "implementations_grid",
                          "efa_list", "settings"))
})

test_that("settings are returned correctly", {
  expect_named(efa_def$settings, c("method", "rotation", "type", "n_factors", "N",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "varimax_type", "normalize",
                                   "k_promax", "k_simplimax", "P_type",
                                   "precision", "start_method", "use",
                                   "cor_method", "max_iter", "averaging",
                                   "trim", "salience_threshold"))
  expect_named(efa_ml$settings, c("method", "rotation", "type", "n_factors", "N",
                                  "init_comm", "criterion", "criterion_type",
                                  "abs_eigen", "varimax_type", "normalize",
                                  "k_promax", "k_simplimax", "P_type",
                                  "precision", "start_method", "use",
                                  "cor_method", "max_iter", "averaging",
                                  "trim", "salience_threshold"))
  expect_named(efa_uls$settings, c("method", "rotation", "type", "n_factors", "N",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "varimax_type", "normalize",
                                   "k_promax", "k_simplimax", "P_type",
                                   "precision", "start_method", "use",
                                   "cor_method", "max_iter", "averaging",
                                   "trim", "salience_threshold"))
  expect_named(efa_all$settings, c("method", "rotation", "type", "n_factors", "N",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "varimax_type", "normalize",
                                   "k_promax", "k_simplimax", "P_type",
                                   "precision", "start_method", "use",
                                   "cor_method", "max_iter", "averaging",
                                   "trim", "salience_threshold"))
  expect_named(efa_all_oblq$settings, c("method", "rotation", "type", "n_factors", "N",
                                        "init_comm", "criterion", "criterion_type",
                                        "abs_eigen", "varimax_type", "normalize",
                                        "k_promax", "k_simplimax", "P_type",
                                        "precision", "start_method", "use",
                                        "cor_method", "max_iter", "averaging",
                                        "trim", "salience_threshold"))
  expect_named(efa_all_orth$settings, c("method", "rotation", "type", "n_factors", "N",
                                        "init_comm", "criterion", "criterion_type",
                                        "abs_eigen", "varimax_type", "normalize",
                                        "k_promax", "k_simplimax", "P_type",
                                        "precision", "start_method", "use",
                                        "cor_method", "max_iter", "averaging",
                                        "trim", "salience_threshold"))
  expect_named(efa_all_none$settings, c("method", "rotation", "type", "n_factors", "N",
                                        "init_comm", "criterion", "criterion_type",
                                        "abs_eigen", "varimax_type", "normalize",
                                        "k_promax", "k_simplimax", "P_type",
                                        "precision", "start_method", "use",
                                        "cor_method", "max_iter", "averaging",
                                        "trim", "salience_threshold"))
  expect_named(efa_all_md$settings, c("method", "rotation", "type", "n_factors", "N",
                                      "init_comm", "criterion", "criterion_type",
                                      "abs_eigen", "varimax_type", "normalize",
                                      "k_promax", "k_simplimax", "P_type",
                                      "precision", "start_method", "use",
                                      "cor_method", "max_iter", "averaging",
                                      "trim", "salience_threshold"))
  expect_named(efa_all_tm$settings, c("method", "rotation", "type", "n_factors", "N",
                                      "init_comm", "criterion", "criterion_type",
                                      "abs_eigen", "varimax_type", "normalize",
                                      "k_promax", "k_simplimax", "P_type",
                                      "precision", "start_method", "use",
                                      "cor_method", "max_iter", "averaging",
                                      "trim", "salience_threshold"))
  expect_named(efa_raw$settings, c("method", "rotation", "type", "n_factors", "N",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "varimax_type", "normalize",
                                   "k_promax", "k_simplimax", "P_type",
                                   "precision", "start_method", "use",
                                   "cor_method", "max_iter", "averaging",
                                   "trim", "salience_threshold"))
  expect_named(efa_raw_p$settings, c("method", "rotation", "type", "n_factors", "N",
                                     "init_comm", "criterion", "criterion_type",
                                     "abs_eigen", "varimax_type", "normalize",
                                     "k_promax", "k_simplimax", "P_type",
                                     "precision", "start_method", "use",
                                     "cor_method", "max_iter", "averaging",
                                     "trim", "salience_threshold"))


  expect_equal(efa_def$settings$method, "PAF")
  expect_equal(efa_ml$settings$method, "ML")
  expect_equal(efa_uls$settings$method, "ULS")
  expect_equal(efa_all$settings$method, c("PAF", "ML", "ULS"))
  expect_equal(efa_all_oblq$settings$method, c("PAF", "ML", "ULS"))
  expect_equal(efa_all_orth$settings$method, c("PAF", "ML", "ULS"))
  expect_equal(efa_all_none$settings$method, c("PAF", "ML", "ULS"))
  expect_equal(efa_all_md$settings$method, c("PAF", "ML", "ULS"))
  expect_equal(efa_all_tm$settings$method, c("PAF", "ML", "ULS"))
  expect_equal(efa_raw$settings$method, "PAF")
  expect_equal(efa_raw_p$settings$method, "PAF")

  expect_equal(efa_def$settings$rotation, "promax")
  expect_equal(efa_ml$settings$rotation, "promax")
  expect_equal(efa_uls$settings$rotation, "promax")
  expect_equal(efa_all$settings$rotation, "promax")
  expect_equal(efa_all_oblq$settings$rotation, "oblique")
  expect_equal(efa_all_orth$settings$rotation, "orthogonal")
  expect_equal(efa_all_none$settings$rotation, "none")
  expect_equal(efa_all_md$settings$rotation, "oblique")
  expect_equal(efa_all_tm$settings$rotation, "oblique")
  expect_equal(efa_raw$settings$rotation, "none")
  expect_equal(efa_raw_p$settings$rotation, "promax")

  expect_equal(efa_def$settings$type, "none")
  expect_equal(efa_ml$settings$type, "none")
  expect_equal(efa_uls$settings$type, "none")
  expect_equal(efa_all$settings$type, c("none", "EFAtools", "psych", "SPSS"))
  expect_equal(efa_all_oblq$settings$type, c("none", "EFAtools", "psych", "SPSS"))
  expect_equal(efa_all_orth$settings$type, c("none", "EFAtools", "psych", "SPSS"))
  expect_equal(efa_all_none$settings$type, c("none", "EFAtools", "psych", "SPSS"))
  expect_equal(efa_all_md$settings$type, c("none", "EFAtools", "psych", "SPSS"))
  expect_equal(efa_all_tm$settings$type, c("none", "EFAtools", "psych", "SPSS"))
  expect_equal(efa_raw$settings$type, "none")
  expect_equal(efa_raw_p$settings$type, "none")

  expect_equal(efa_def$settings$n_factors, 3)
  expect_equal(efa_ml$settings$n_factors, 3)
  expect_equal(efa_uls$settings$n_factors, 3)
  expect_equal(efa_all$settings$n_factors, 3)
  expect_equal(efa_all_oblq$settings$n_factors, 3)
  expect_equal(efa_all_orth$settings$n_factors, 3)
  expect_equal(efa_all_none$settings$n_factors, 3)
  expect_equal(efa_all_md$settings$n_factors, 3)
  expect_equal(efa_all_tm$settings$n_factors, 3)
  expect_equal(efa_raw$settings$n_factors, 1)
  expect_equal(efa_raw_p$settings$n_factors, 2)

  expect_equal(efa_def$settings$N, 500)
  expect_equal(efa_ml$settings$N, 500)
  expect_equal(efa_uls$settings$N, 500)
  expect_equal(efa_all$settings$N, 500)
  expect_equal(efa_all_oblq$settings$N, 500)
  expect_equal(efa_all_orth$settings$N, 500)
  expect_equal(efa_all_none$settings$N, 500)
  expect_equal(efa_all_md$settings$N, 500)
  expect_equal(efa_all_tm$settings$N, 500)
  expect_equal(efa_raw$settings$N, 810)
  expect_equal(efa_raw_p$settings$N, 810)

  expect_equal(efa_def$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_ml$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_uls$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_all$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_all_oblq$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_all_orth$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_all_none$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_all_md$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_all_tm$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_raw$settings$init_comm, c("smc", "mac", "unity"))
  expect_equal(efa_raw_p$settings$init_comm, c("smc", "mac", "unity"))

  expect_equal(efa_def$settings$criterion, 0.001)
  expect_equal(efa_ml$settings$criterion, 0.001)
  expect_equal(efa_uls$settings$criterion, 0.001)
  expect_equal(efa_all$settings$criterion, 0.001)
  expect_equal(efa_all_oblq$settings$criterion, 0.001)
  expect_equal(efa_all_orth$settings$criterion, 0.001)
  expect_equal(efa_all_none$settings$criterion, 0.001)
  expect_equal(efa_all_md$settings$criterion, 0.001)
  expect_equal(efa_all_tm$settings$criterion, 0.001)
  expect_equal(efa_raw$settings$criterion, 0.001)
  expect_equal(efa_raw_p$settings$criterion, 0.001)

  expect_equal(efa_def$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_ml$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_uls$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_all$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_all_oblq$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_all_orth$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_all_none$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_all_md$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_all_tm$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_raw$settings$criterion_type, c("sum", "max_individual"))
  expect_equal(efa_raw_p$settings$criterion_type, c("sum", "max_individual"))

  expect_equal(efa_def$settings$abs_eigen, TRUE)
  expect_equal(efa_ml$settings$abs_eigen, TRUE)
  expect_equal(efa_uls$settings$abs_eigen, TRUE)
  expect_equal(efa_all$settings$abs_eigen, TRUE)
  expect_equal(efa_all_oblq$settings$abs_eigen, TRUE)
  expect_equal(efa_all_orth$settings$abs_eigen, TRUE)
  expect_equal(efa_all_none$settings$abs_eigen, TRUE)
  expect_equal(efa_all_md$settings$abs_eigen, TRUE)
  expect_equal(efa_all_tm$settings$abs_eigen, TRUE)
  expect_equal(efa_raw$settings$abs_eigen, TRUE)
  expect_equal(efa_raw_p$settings$abs_eigen, TRUE)

  expect_equal(efa_def$settings$abs_eigen, TRUE)
  expect_equal(efa_ml$settings$abs_eigen, TRUE)
  expect_equal(efa_uls$settings$abs_eigen, TRUE)
  expect_equal(efa_all$settings$abs_eigen, TRUE)
  expect_equal(efa_all_oblq$settings$abs_eigen, TRUE)
  expect_equal(efa_all_orth$settings$abs_eigen, TRUE)
  expect_equal(efa_all_none$settings$abs_eigen, TRUE)
  expect_equal(efa_all_md$settings$abs_eigen, TRUE)
  expect_equal(efa_all_tm$settings$abs_eigen, TRUE)
  expect_equal(efa_raw$settings$abs_eigen, TRUE)
  expect_equal(efa_raw_p$settings$abs_eigen, TRUE)

  expect_equal(efa_def$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_ml$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_uls$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_all$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_all_oblq$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_all_orth$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_all_none$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_all_md$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_all_tm$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_raw$settings$varimax_type, c("svd", "kaiser"))
  expect_equal(efa_raw_p$settings$varimax_type, c("svd", "kaiser"))

  expect_equal(efa_def$settings$normalize, TRUE)
  expect_equal(efa_ml$settings$normalize, TRUE)
  expect_equal(efa_uls$settings$normalize, TRUE)
  expect_equal(efa_all$settings$normalize, TRUE)
  expect_equal(efa_all_oblq$settings$normalize, TRUE)
  expect_equal(efa_all_orth$settings$normalize, TRUE)
  expect_equal(efa_all_none$settings$normalize, TRUE)
  expect_equal(efa_all_md$settings$normalize, TRUE)
  expect_equal(efa_all_tm$settings$normalize, TRUE)
  expect_equal(efa_raw$settings$normalize, TRUE)
  expect_equal(efa_raw_p$settings$normalize, TRUE)

  expect_equal(efa_def$settings$k_promax, 2:4)
  expect_equal(efa_ml$settings$k_promax, 2:4)
  expect_equal(efa_uls$settings$k_promax, 2:4)
  expect_equal(efa_all$settings$k_promax, 2:4)
  expect_equal(efa_all_oblq$settings$k_promax, 2:4)
  expect_equal(efa_all_orth$settings$k_promax, 2:4)
  expect_equal(efa_all_none$settings$k_promax, 2:4)
  expect_equal(efa_all_md$settings$k_promax, 2:4)
  expect_equal(efa_all_tm$settings$k_promax, 2:4)
  expect_equal(efa_raw$settings$k_promax, 2:4)
  expect_equal(efa_raw_p$settings$k_promax, 2:4)

  expect_equal(efa_def$settings$k_simplimax, 18)
  expect_equal(efa_ml$settings$k_simplimax, 18)
  expect_equal(efa_uls$settings$k_simplimax, 18)
  expect_equal(efa_all$settings$k_simplimax, 18)
  expect_equal(efa_all_oblq$settings$k_simplimax, 18)
  expect_equal(efa_all_orth$settings$k_simplimax, 18)
  expect_equal(efa_all_none$settings$k_simplimax, 18)
  expect_equal(efa_all_md$settings$k_simplimax, 18)
  expect_equal(efa_all_tm$settings$k_simplimax, 18)
  expect_equal(efa_raw$settings$k_simplimax, 8)
  expect_equal(efa_raw_p$settings$k_simplimax, 8)

  expect_equal(efa_def$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_ml$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_uls$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_all$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_all_oblq$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_all_orth$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_all_none$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_all_md$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_all_tm$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_raw$settings$P_type, c("norm", "unnorm"))
  expect_equal(efa_raw_p$settings$P_type, c("norm", "unnorm"))

  expect_equal(efa_def$settings$precision, 1e-5)
  expect_equal(efa_ml$settings$precision, 1e-5)
  expect_equal(efa_uls$settings$precision, 1e-5)
  expect_equal(efa_all$settings$precision, 1e-5)
  expect_equal(efa_all_oblq$settings$precision, 1e-5)
  expect_equal(efa_all_orth$settings$precision, 1e-5)
  expect_equal(efa_all_none$settings$precision, 1e-5)
  expect_equal(efa_all_md$settings$precision, 1e-5)
  expect_equal(efa_all_tm$settings$precision, 1e-5)
  expect_equal(efa_raw$settings$precision, 1e-5)
  expect_equal(efa_raw_p$settings$precision, 1e-5)

  expect_equal(efa_def$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_ml$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_uls$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_all$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_all_oblq$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_all_orth$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_all_none$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_all_md$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_all_tm$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_raw$settings$start_method, c("psych", "factanal"))
  expect_equal(efa_raw_p$settings$start_method, c("psych", "factanal"))

  expect_equal(efa_def$settings$use, "pairwise.complete.obs")
  expect_equal(efa_ml$settings$use, "pairwise.complete.obs")
  expect_equal(efa_uls$settings$use, "pairwise.complete.obs")
  expect_equal(efa_all$settings$use, "pairwise.complete.obs")
  expect_equal(efa_all_oblq$settings$use, "pairwise.complete.obs")
  expect_equal(efa_all_orth$settings$use, "pairwise.complete.obs")
  expect_equal(efa_all_none$settings$use, "pairwise.complete.obs")
  expect_equal(efa_all_md$settings$use, "pairwise.complete.obs")
  expect_equal(efa_all_tm$settings$use, "pairwise.complete.obs")
  expect_equal(efa_raw$settings$use, "pairwise.complete.obs")
  expect_equal(efa_raw_p$settings$use, "pairwise.complete.obs")

  expect_equal(efa_def$settings$cor_method, "pearson")
  expect_equal(efa_ml$settings$cor_method, "pearson")
  expect_equal(efa_uls$settings$cor_method, "pearson")
  expect_equal(efa_all$settings$cor_method, "pearson")
  expect_equal(efa_all_oblq$settings$cor_method, "pearson")
  expect_equal(efa_all_orth$settings$cor_method, "pearson")
  expect_equal(efa_all_none$settings$cor_method, "pearson")
  expect_equal(efa_all_md$settings$cor_method, "pearson")
  expect_equal(efa_all_tm$settings$cor_method, "pearson")
  expect_equal(efa_raw$settings$cor_method, "pearson")
  expect_equal(efa_raw_p$settings$cor_method, "pearson")

  expect_equal(efa_def$settings$max_iter, 10000)
  expect_equal(efa_ml$settings$max_iter, 10000)
  expect_equal(efa_uls$settings$max_iter, 10000)
  expect_equal(efa_all$settings$max_iter, 10000)
  expect_equal(efa_all_oblq$settings$max_iter, 10000)
  expect_equal(efa_all_orth$settings$max_iter, 10000)
  expect_equal(efa_all_none$settings$max_iter, 10000)
  expect_equal(efa_all_md$settings$max_iter, 10000)
  expect_equal(efa_all_tm$settings$max_iter, 10000)
  expect_equal(efa_raw$settings$max_iter, 10000)
  expect_equal(efa_raw_p$settings$max_iter, 10000)

  expect_equal(efa_def$settings$averaging, "mean")
  expect_equal(efa_ml$settings$averaging, "mean")
  expect_equal(efa_uls$settings$averaging, "mean")
  expect_equal(efa_all$settings$averaging, "mean")
  expect_equal(efa_all_oblq$settings$averaging, "mean")
  expect_equal(efa_all_orth$settings$averaging, "mean")
  expect_equal(efa_all_none$settings$averaging, "mean")
  expect_equal(efa_all_md$settings$averaging, "median")
  expect_equal(efa_all_tm$settings$averaging, "mean")
  expect_equal(efa_raw$settings$averaging, "mean")
  expect_equal(efa_raw_p$settings$averaging, "mean")

  expect_equal(efa_def$settings$trim, 0)
  expect_equal(efa_ml$settings$trim, 0)
  expect_equal(efa_uls$settings$trim, 0)
  expect_equal(efa_all$settings$trim, 0)
  expect_equal(efa_all_oblq$settings$trim, 0)
  expect_equal(efa_all_orth$settings$trim, 0)
  expect_equal(efa_all_none$settings$trim, 0)
  expect_equal(efa_all_md$settings$trim, 0)
  expect_equal(efa_all_tm$settings$trim, 0.2)
  expect_equal(efa_raw$settings$trim, 0)
  expect_equal(efa_raw_p$settings$trim, 0)

  expect_equal(efa_def$settings$salience_threshold, 0.3)
  expect_equal(efa_ml$settings$salience_threshold, 0.3)
  expect_equal(efa_uls$settings$salience_threshold, 0.3)
  expect_equal(efa_all$settings$salience_threshold, 0.2)
  expect_equal(efa_all_oblq$settings$salience_threshold, 0.3)
  expect_equal(efa_all_orth$settings$salience_threshold, 0.3)
  expect_equal(efa_all_none$settings$salience_threshold, 0.3)
  expect_equal(efa_all_md$settings$salience_threshold, 0.3)
  expect_equal(efa_all_tm$settings$salience_threshold, 0.3)
  expect_equal(efa_raw$settings$salience_threshold, 0.3)
  expect_equal(efa_raw_p$settings$salience_threshold, 0.3)


})


# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors are thrown correctly", {
  expect_error(EFA_AVERAGE(1:5), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n")
  expect_message(EFA_AVERAGE(GRiPS_raw, n_factors = 2, method = "PAF", type = c("EFAtools", "psych")),
                 " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n")
  expect_warning(EFA_AVERAGE(GRiPS_raw, n_factors = 2, method = "PAF", type = c("EFAtools", "psych"),
                             N = 20),
                 " 'N' was set and data entered. Taking N from data.\n")
  expect_error(EFA_AVERAGE(dat_sing, n_factors = 1),
               " Correlation matrix is singular, no further analyses are performed\n")
  expect_error(EFA_AVERAGE(cor_sing, N = 10, n_factors = 1),
               " Correlation matrix is singular, no further analyses are performed\n")
  expect_error(EFA_AVERAGE(matrix(rnorm(30), ncol = 3), n_factors = 2),
               " The model is underidentified. Please enter a lower number of factors or use a larger number of indicators and try again.\n")
  expect_warning(EFA_AVERAGE(matrix(rnorm(30), ncol = 3), n_factors = 1,
                             method = "PAF", type = c("EFAtools", "psych")),
                 " The model is just identified (df = 0). We suggest to try again with a lower number of factors or a larger number of indicators.\n", fixed = TRUE)
  expect_warning(EFA_AVERAGE(cor_nposdef, n_factors = 1, N = 10, method = "PAF",
                     type = c("EFAtools", "psych")), "Matrix was not positive definite, smoothing was done")
  expect_message(EFA_AVERAGE(GRiPS_raw, n_factors = 1, method = "PAF", type = c("EFAtools", "psych")),
                 " 'n_factors' is 1, but rotation != 'none'. Setting rotation to 'none' to avoid many warnings, as 1-factor solutions cannot be rotated.\n")
  expect_warning(EFA_AVERAGE(GRiPS_raw, n_factors = 1, method = "PAF", type = c("EFAtools"),
                             rotation = "none"),
                 " There was only one combination of arguments, returning normal EFA output.\n")
})

rm(efa_def, efa_ml, efa_uls, efa_all, efa_all_oblq, efa_all_orth, efa_all_none,
   efa_all_md, efa_all_tm, efa_raw, efa_raw_p)



