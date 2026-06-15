# EFA_AVERAGE runs its grid of EFAs through future.apply with future.seed = TRUE,
# so the random starts of the oblique rotation engines draw from the session RNG.
# Seed once here so the fixtures below - and the print snapshot built from them -
# are reproducible regardless of what ran earlier in the test run.
set.seed(42)

efa_def <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                       show_progress = FALSE)
efa_ml <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                      method = "ML", show_progress = FALSE)
efa_uls <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                      method = "ULS", show_progress = FALSE)

efa_all <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                       method = c("PAF", "ML", "ULS"),
                       type = c("none", "EFAtools", "psych", "SPSS"),
                       salience_threshold = .2, show_progress = FALSE)
efa_all_oblq <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                            method = c("PAF", "ML", "ULS"),
                            type = c("none", "EFAtools", "psych", "SPSS"),
                            rotation = "oblique", show_progress = FALSE)
efa_all_orth <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                            method = c("PAF", "ML", "ULS"),
                            type = c("none", "EFAtools", "psych", "SPSS"),
                            rotation = "orthogonal", show_progress = FALSE)
efa_all_none <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                            method = c("PAF", "ML", "ULS"),
                            type = c("none", "EFAtools", "psych", "SPSS"),
                            rotation = "none", show_progress = FALSE)

efa_all_md <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                            method = c("PAF", "ML", "ULS"),
                            type = c("none", "EFAtools", "psych", "SPSS"),
                            rotation = "oblique", averaging = "median", show_progress = FALSE)
efa_all_tm <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                          method = c("PAF", "ML", "ULS"),
                          type = c("none", "EFAtools", "psych", "SPSS"),
                          rotation = "oblique", averaging = "mean",
                          trim = .2, show_progress = FALSE)
efa_raw <- EFA_AVERAGE(GRiPS_raw, n_factors = 1, rotation = "none", show_progress = FALSE)
efa_raw_p <- EFA_AVERAGE(GRiPS_raw, n_factors = 2, rotation = "promax", show_progress = FALSE)

test_that("output class and dimensions are correct", {
  expect_s3_class(efa_def, "EFA_AVERAGE")
  expect_s3_class(efa_ml, "EFA_AVERAGE")
  expect_s3_class(efa_uls, "EFA_AVERAGE")
  expect_s3_class(efa_all, "EFA_AVERAGE")
  expect_s3_class(efa_all_oblq, "EFA_AVERAGE")
  expect_s3_class(efa_all_orth, "EFA_AVERAGE")
  expect_s3_class(efa_all_none, "EFA_AVERAGE")
  expect_s3_class(efa_all_md, "EFA_AVERAGE")
  expect_s3_class(efa_all_tm, "EFA_AVERAGE")
  expect_s3_class(efa_raw, "EFA_AVERAGE")
  expect_s3_class(efa_raw_p, "EFA_AVERAGE")

  expect_s3_class(efa_def$loadings$average, "LOADINGS")
  expect_s3_class(efa_ml$loadings$average, "LOADINGS")
  expect_s3_class(efa_uls$loadings$average, "LOADINGS")
  expect_s3_class(efa_all$loadings$average, "LOADINGS")
  expect_s3_class(efa_all_oblq$loadings$average, "LOADINGS")
  expect_s3_class(efa_all_orth$loadings$average, "LOADINGS")
  expect_s3_class(efa_all_none$loadings$average, "LOADINGS")
  expect_s3_class(efa_all_md$loadings$average, "LOADINGS")
  expect_s3_class(efa_all_tm$loadings$average, "LOADINGS")
  expect_s3_class(efa_raw$loadings$average, "LOADINGS")
  expect_s3_class(efa_raw_p$loadings$average, "LOADINGS")

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
  expect_error(EFA_AVERAGE(1:5, show_progress = FALSE), class = "efa_input_not_matrix")
  expect_message(EFA_AVERAGE(GRiPS_raw, n_factors = 2, method = "PAF", type = c("EFAtools", "psych"), show_progress = FALSE),
                 class = "efa_cor_from_data")
  expect_warning(EFA_AVERAGE(GRiPS_raw, n_factors = 2, method = "PAF", type = c("EFAtools", "psych"),
                             N = 20, show_progress = FALSE),
                 class = "efa_n_from_data")
  expect_error(EFA_AVERAGE(dat_sing, n_factors = 1, show_progress = FALSE),
               class = "efa_cor_singular")
  expect_error(EFA_AVERAGE(cor_sing, N = 10, n_factors = 1, show_progress = FALSE),
               class = "efa_cor_singular")
  expect_error(EFA_AVERAGE(matrix(rnorm(30), ncol = 3), n_factors = 2, show_progress = FALSE),
               class = "efa_underidentified")
  expect_warning(EFA_AVERAGE(matrix(rnorm(30), ncol = 3), n_factors = 1,
                             method = "PAF", type = c("EFAtools", "psych"), show_progress = FALSE),
                 class = "efa_just_identified")
  expect_warning(EFA_AVERAGE(cor_nposdef, n_factors = 1, N = 10, method = "PAF",
                     type = c("EFAtools", "psych"), show_progress = FALSE), class = "efa_cor_smoothed")
  expect_message(EFA_AVERAGE(GRiPS_raw, n_factors = 1, method = "PAF", type = c("EFAtools", "psych"), show_progress = FALSE),
                 class = "efa_avg_single_factor_rotation")
  expect_warning(EFA_AVERAGE(GRiPS_raw, n_factors = 1, method = "PAF", type = c("EFAtools"),
                             rotation = "none", show_progress = FALSE),
                 class = "efa_avg_single_combination")
})

test_that("an all-failed averaging grid returns an empty (NA) result", {
  # When every solution fails (here all runs hit max_iter and do not converge),
  # the averaged result is NA rather than an error or an average over an empty set.
  res <- suppressWarnings(EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3,
                                      N = 500, max_iter = 1, method = "PAF",
                                      type = "none", rotation = "none",
                                      show_progress = FALSE))
  expect_s3_class(res, "EFA_AVERAGE")
  expect_true(all((res$implementations_grid$converged != 0) %in% TRUE))
  expect_true(all(is.na(res$h2)))
  expect_true(all(is.na(res$loadings)))
})

test_that("an all-Heywood averaging grid returns an empty (NA) result", {
  # A 3-variable matrix whose single-factor solution implies a communality > 1, so
  # every fitted solution is a Heywood case and is excluded; with nothing left to
  # average, the communalities and loadings come back as NA rather than an error.
  m <- matrix(c(1, .8, .8,  .8, 1, .5,  .8, .5, 1), 3, 3)
  res <- suppressWarnings(EFA_AVERAGE(m, n_factors = 1, N = 200, method = "PAF",
                                      type = "none", rotation = "none",
                                      show_progress = FALSE))
  expect_s3_class(res, "EFA_AVERAGE")
  expect_true(all(res$implementations_grid$heywood))
  expect_true(all(is.na(res$h2)))
  expect_true(all(is.na(res$loadings)))
})

test_that("print output is stable", {
  skip_on_cran()
  local_reproducible_output()

  # default view (oblique -> Phi section, average + range)
  expect_snapshot(print(efa_def, plot = FALSE), transform = scrub_num_pct)

  # all statistics, no Phi (orthogonal-free rotation = "none")
  expect_snapshot(print(efa_all_none, stat = c("average", "sd", "min", "max"),
                        plot = FALSE), transform = scrub_num_pct)

  # median averaging (Md / Median labels)
  expect_snapshot(print(efa_all_md, plot = FALSE), transform = scrub_num_pct)
})

test_that("plot returns a ggplot", {
  skip_on_cran()
  expect_s3_class(plot(efa_def), "ggplot")
})

test_that("a vector-valued precision is recycled across the grid", {
  # A vector precision is expanded into the grid; each EFA must receive its own
  # scalar value rather than the whole vector (which would fail every fit and
  # return an all-NA result).
  res <- suppressWarnings(
    EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                method = "ML", rotation = "promax", type = "none",
                precision = c(1e-5, 1e-3), show_progress = FALSE))
  expect_s3_class(res, "EFA_AVERAGE")
  expect_setequal(unique(res$implementations_grid$precision), c(1e-5, 1e-3))
  expect_false(any(res$implementations_grid$errors))
  expect_false(all(is.na(res$loadings$average)))
  expect_true(all(is.finite(res$loadings$average)))
})

test_that("problematic solutions are summarised in a single warning", {
  # All runs hit max_iter and are excluded, so one summary warning is raised
  # rather than one per model.
  expect_warning(
    EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                max_iter = 1, method = "PAF", type = "none", rotation = "none",
                show_progress = FALSE),
    class = "efa_avg_excluded_solutions"
  )
})

test_that("printing does not plot by default", {
  expect_false(formals(getS3method("print", "EFA_AVERAGE"))$plot)
})

test_that("admissibility is reported as an outcome, not a varied setting", {
  skip_on_cran()
  local_reproducible_output()
  obj <- efa_def
  # Force mixed admissibility so the column would surface in the varied-settings
  # list if it were (incorrectly) treated as a setting.
  obj$implementations_grid$admissible[1] <- !obj$implementations_grid$admissible[1]
  out <- paste(cli::ansi_strip(utils::capture.output(print(obj, plot = FALSE))),
               collapse = " ")
  # "admissible" appears only in the convergence summary, never in the settings list.
  expect_equal(unname(lengths(gregexpr("admissible", out, fixed = TRUE))), 1L)
})

test_that("averaged fit indices match the per-model grid means", {
  # efa_def uses the untrimmed mean, so each averaged index equals the column
  # mean over the included models (and the new residual indices are averaged too).
  fi <- efa_def$fit_indices
  grid <- efa_def$implementations_grid
  expect_equal(fi$average[fi$index == "caf"], mean(grid$caf, na.rm = TRUE))
  expect_equal(fi$average[fi$index == "srmr"], mean(grid$srmr, na.rm = TRUE))
})

rm(efa_def, efa_ml, efa_uls, efa_all, efa_all_oblq, efa_all_orth, efa_all_none,
   efa_all_md, efa_all_tm, efa_raw, efa_raw_p)



