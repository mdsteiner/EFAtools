efa_cor <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500)
efa_raw <- EFA(GRiPS_raw, n_factors = 1)

# different types
efa_psych <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                 type = "psych", rotation = "promax")
efa_spss <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                type = "SPSS", rotation = "promax")

# different methods
efa_ml <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
              method = "ML")
efa_uls <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
               method = "ULS")

# different rotation methods from GPA rotation package (orthogonal and oblique)
efa_equa <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                rotation = "equamax")
efa_quart <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                 rotation = "quartimin")

# PAF with promax rotation without a specified type
efa_none <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                type = "none", method = "PAF", rotation = "promax",
                max_iter = 500, init_comm = "unity", criterion = 1e-4,
                criterion_type = "sum", abs_eigen = FALSE, k = 3,
                P_type = "unnorm", precision= 1e-5, order_type = "eigen",
                varimax_type = "svd")

# create correlation matrices from population models
cormat_zero <- population_models$loadings$baseline %*% population_models$phis_3$zero %*% t(population_models$loadings$baseline)
diag(cormat_zero) <- 1

cormat_moderate <- population_models$loadings$baseline %*% population_models$phis_3$moderate %*% t(population_models$loadings$baseline)
diag(cormat_moderate) <- 1

efa_paf_zero <- EFA(cormat_zero, 3, 500, rotation = "varimax")
efa_ml_zero <- EFA(cormat_zero, 3, 500, method = "ML", rotation = "varimax")
efa_uls_zero <- EFA(cormat_zero, 3, 500, method = "ULS", rotation = "varimax")

efa_paf_moderate <- EFA(cormat_moderate, 3, 500, rotation = "promax")
efa_ml_moderate <- EFA(cormat_moderate, 3, 500, method = "ML",
                       rotation = "promax")
efa_uls_moderate <- EFA(cormat_moderate, 3, 500, method = "ULS",
                        rotation = "promax")


test_that("output class and dimensions are correct", {
  # Objects carry the new `efa` class first and keep the legacy `EFA` string for
  # inherits()/back-compat; the loading matrices are dual-classed the same way.
  expect_identical(class(efa_cor), c("efa", "EFA"))
  expect_s3_class(efa_raw, "EFA")
  expect_s3_class(efa_psych, "EFA")
  expect_s3_class(efa_spss, "EFA")
  expect_s3_class(efa_ml, "EFA")
  expect_s3_class(efa_uls, "EFA")
  expect_s3_class(efa_equa, "EFA")
  expect_s3_class(efa_quart, "EFA")
  expect_s3_class(efa_none, "EFA")

  expect_identical(class(efa_cor$unrot_loadings), c("efa_loadings", "LOADINGS"))
  expect_s3_class(efa_raw$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_psych$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_spss$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_ml$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_uls$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_equa$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_quart$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_none$unrot_loadings, "LOADINGS")

  expect_identical(class(efa_psych$rot_loadings), c("efa_loadings", "LOADINGS"))
  expect_s3_class(efa_spss$rot_loadings, "LOADINGS")
  expect_s3_class(efa_equa$rot_loadings, "LOADINGS")
  expect_s3_class(efa_quart$rot_loadings, "LOADINGS")
  expect_s3_class(efa_none$rot_loadings, "LOADINGS")

  expect_named(efa_cor, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                          "final_eigen", "iter", "convergence", "heywood", "unrot_loadings",
                          "vars_accounted", "fit_indices", "model_implied_R",
                          "residuals", "settings", "vcov_unrot_loadings", "Gamma"))
  expect_named(efa_raw, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                          "final_eigen", "iter", "convergence", "heywood", "unrot_loadings",
                          "vars_accounted", "fit_indices", "model_implied_R",
                          "residuals", "settings", "vcov_unrot_loadings", "Gamma"))
  expect_named(efa_psych, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                            "final_eigen", "iter", "convergence", "heywood", "unrot_loadings",
                            "vars_accounted", "fit_indices", "model_implied_R",
                            "residuals", "rot_loadings",
                            "Phi", "Structure", "rotmat", "vars_accounted_rot",
                            "settings", "vcov_unrot_loadings", "Gamma"))
  expect_named(efa_spss, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                           "final_eigen", "iter", "convergence", "heywood", "unrot_loadings",
                           "vars_accounted", "fit_indices", "model_implied_R",
                           "residuals", "rot_loadings",
                           "Phi", "Structure", "rotmat", "vars_accounted_rot",
                           "settings", "vcov_unrot_loadings", "Gamma"))
  expect_named(efa_ml, c("orig_R", "h2", "orig_eigen", "final_eigen", "iter",
                         "convergence", "heywood", "unrot_loadings", "vars_accounted",
                         "fit_indices", "model_implied_R",
                         "residuals", "settings", "vcov_unrot_loadings", "Gamma"))
  expect_named(efa_uls, c("orig_R", "h2", "orig_eigen", "final_eigen", "iter",
                         "convergence", "heywood", "unrot_loadings", "vars_accounted",
                         "fit_indices", "model_implied_R",
                         "residuals", "settings", "vcov_unrot_loadings", "Gamma"))
  expect_named(efa_equa, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                           "final_eigen", "iter", "convergence", "heywood", "unrot_loadings",
                           "vars_accounted", "fit_indices", "model_implied_R",
                           "residuals", "rot_loadings",
                           "rotmat", "vars_accounted_rot",
                           "settings", "vcov_unrot_loadings", "Gamma"))
  expect_named(efa_quart, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                            "final_eigen", "iter", "convergence", "heywood", "unrot_loadings",
                            "vars_accounted", "fit_indices", "model_implied_R",
                            "residuals", "rot_loadings",
                            "Phi", "Structure", "rotmat", "vars_accounted_rot",
                            "settings", "vcov_unrot_loadings", "Gamma"))
  expect_named(efa_none, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                           "final_eigen", "iter", "convergence", "heywood", "unrot_loadings",
                           "vars_accounted", "fit_indices", "model_implied_R",
                           "residuals", "rot_loadings",
                           "Phi", "Structure", "rotmat", "vars_accounted_rot",
                           "settings", "vcov_unrot_loadings", "Gamma"))
})

test_that("settings are returned correctly", {
  expect_named(efa_cor$settings, c("estimator", "method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "input_type",
                                   "cor_method_used", "se", "b_boot", "ci", "seed", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen"))
  expect_named(efa_raw$settings, c("estimator", "method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "input_type",
                                   "cor_method_used", "se", "b_boot", "ci", "seed", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen"))
  expect_named(efa_psych$settings, c("estimator", "method", "rotation", "type", "n_factors",
                                     "N", "use", "cor_method", "input_type",
                                     "cor_method_used", "se", "b_boot", "ci", "seed", "max_iter",
                                     "init_comm", "criterion", "criterion_type",
                                     "abs_eigen", "normalize", "P_type", "precision",
                                     "order_type", "varimax_type", "k"))
  expect_named(efa_spss$settings, c("estimator", "method", "rotation", "type", "n_factors",
                                    "N", "use", "cor_method", "input_type",
                                    "cor_method_used", "se", "b_boot", "ci", "seed", "max_iter",
                                    "init_comm", "criterion", "criterion_type",
                                    "abs_eigen", "normalize", "P_type", "precision",
                                    "order_type", "varimax_type", "k"))
  expect_named(efa_ml$settings, c("estimator", "method", "rotation", "type", "n_factors",
                                    "N", "use", "cor_method", "input_type",
                                    "cor_method_used", "se", "b_boot", "ci", "seed", "start_method"))
  expect_named(efa_uls$settings, c("estimator", "method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "input_type",
                                   "cor_method_used", "se", "b_boot", "ci", "seed"))
  expect_named(efa_equa$settings, c("estimator", "method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "input_type",
                                   "cor_method_used", "se", "b_boot", "ci", "seed", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "normalize", "precision",
                                   "order_type", "randomStarts", "rotation_diagnostics"))
  expect_named(efa_quart$settings, c("estimator", "method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "input_type",
                                   "cor_method_used", "se", "b_boot", "ci", "seed", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "normalize", "precision",
                                   "order_type", "k", "randomStarts", "rotation_diagnostics"))
  expect_named(efa_none$settings, c("estimator", "method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "input_type",
                                   "cor_method_used", "se", "b_boot", "ci", "seed", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "normalize", "P_type", "precision",
                                   "order_type", "varimax_type", "k"))

  expect_equal(efa_cor$settings$estimator, "PAF")
  expect_equal(efa_raw$settings$estimator, "PAF")
  expect_equal(efa_psych$settings$estimator, "PAF")
  expect_equal(efa_spss$settings$estimator, "PAF")
  expect_equal(efa_ml$settings$estimator, "ML")
  expect_equal(efa_uls$settings$estimator, "ULS")
  expect_equal(efa_equa$settings$estimator, "PAF")
  expect_equal(efa_quart$settings$estimator, "PAF")
  expect_equal(efa_none$settings$estimator, "PAF")

  expect_equal(efa_cor$settings$rotation, "none")
  expect_equal(efa_raw$settings$rotation, "none")
  expect_equal(efa_psych$settings$rotation, "promax")
  expect_equal(efa_spss$settings$rotation, "promax")
  expect_equal(efa_ml$settings$rotation, "none")
  expect_equal(efa_uls$settings$rotation, "none")
  expect_equal(efa_equa$settings$rotation, "equamax")
  expect_equal(efa_quart$settings$rotation, "quartimin")
  expect_equal(efa_none$settings$rotation, "promax")

  expect_equal(efa_cor$settings$type, "EFAtools")
  expect_equal(efa_raw$settings$type, "EFAtools")
  expect_equal(efa_psych$settings$type, "psych")
  expect_equal(efa_spss$settings$type, "SPSS")
  expect_equal(efa_ml$settings$type, "EFAtools")
  expect_equal(efa_uls$settings$type, "EFAtools")
  expect_equal(efa_equa$settings$type, "EFAtools")
  expect_equal(efa_quart$settings$type, "EFAtools")
  expect_equal(efa_none$settings$type, "none")

  expect_equal(efa_cor$settings$n_factors, 3)
  expect_equal(efa_raw$settings$n_factors, 1)
  expect_equal(efa_psych$settings$n_factors, 3)
  expect_equal(efa_spss$settings$n_factors, 3)
  expect_equal(efa_ml$settings$n_factors, 3)
  expect_equal(efa_uls$settings$n_factors, 3)
  expect_equal(efa_equa$settings$n_factors, 3)
  expect_equal(efa_quart$settings$n_factors, 3)
  expect_equal(efa_none$settings$n_factors, 3)

  expect_equal(efa_cor$settings$N, 500)
  expect_equal(efa_raw$settings$N, 810)
  expect_equal(efa_psych$settings$N, 500)
  expect_equal(efa_spss$settings$N, 500)
  expect_equal(efa_ml$settings$N, 500)
  expect_equal(efa_uls$settings$N, 500)
  expect_equal(efa_equa$settings$N, 500)
  expect_equal(efa_quart$settings$N, 500)
  expect_equal(efa_none$settings$N, 500)

  expect_equal(efa_cor$settings$use, "pairwise.complete.obs")
  expect_equal(efa_raw$settings$use, "pairwise.complete.obs")
  expect_equal(efa_psych$settings$use, "pairwise.complete.obs")
  expect_equal(efa_spss$settings$use, "pairwise.complete.obs")
  expect_equal(efa_ml$settings$use, "pairwise.complete.obs")
  expect_equal(efa_uls$settings$use, "pairwise.complete.obs")
  expect_equal(efa_equa$settings$use, "pairwise.complete.obs")
  expect_equal(efa_quart$settings$use, "pairwise.complete.obs")
  expect_equal(efa_none$settings$use, "pairwise.complete.obs")

  expect_equal(efa_cor$settings$cor_method, "pearson")
  expect_equal(efa_raw$settings$cor_method, "pearson")
  expect_equal(efa_psych$settings$cor_method, "pearson")
  expect_equal(efa_spss$settings$cor_method, "pearson")
  expect_equal(efa_ml$settings$cor_method, "pearson")
  expect_equal(efa_uls$settings$cor_method, "pearson")
  expect_equal(efa_equa$settings$cor_method, "pearson")
  expect_equal(efa_quart$settings$cor_method, "pearson")
  expect_equal(efa_none$settings$cor_method, "pearson")

  expect_equal(efa_cor$settings$max_iter, 300)
  expect_equal(efa_raw$settings$max_iter, 300)
  expect_equal(efa_psych$settings$max_iter, 50)
  expect_equal(efa_spss$settings$max_iter, 25)
  expect_equal(efa_equa$settings$max_iter, 300)
  expect_equal(efa_quart$settings$max_iter, 300)
  expect_equal(efa_none$settings$max_iter, 500)

  expect_equal(efa_cor$settings$init_comm, "smc")
  expect_equal(efa_raw$settings$init_comm, "smc")
  expect_equal(efa_psych$settings$init_comm, "smc")
  expect_equal(efa_spss$settings$init_comm, "smc")
  expect_equal(efa_equa$settings$init_comm, "smc")
  expect_equal(efa_quart$settings$init_comm, "smc")
  expect_equal(efa_none$settings$init_comm, "unity")

  expect_equal(efa_cor$settings$criterion, 0.001)
  expect_equal(efa_raw$settings$criterion,  0.001)
  expect_equal(efa_psych$settings$criterion,  0.001)
  expect_equal(efa_spss$settings$criterion,  0.001)
  expect_equal(efa_equa$settings$criterion,  0.001)
  expect_equal(efa_quart$settings$criterion,  0.001)
  expect_equal(efa_none$settings$criterion,  1e-4)

  expect_equal(efa_cor$settings$criterion_type, "sum")
  expect_equal(efa_raw$settings$criterion_type, "sum")
  expect_equal(efa_psych$settings$criterion_type, "sum")
  expect_equal(efa_spss$settings$criterion_type, "max_individual")
  expect_equal(efa_equa$settings$criterion_type, "sum")
  expect_equal(efa_quart$settings$criterion_type, "sum")
  expect_equal(efa_none$settings$criterion_type, "sum")

  expect_equal(efa_cor$settings$abs_eigen, TRUE)
  expect_equal(efa_raw$settings$abs_eigen,  TRUE)
  expect_equal(efa_psych$settings$abs_eigen, FALSE)
  expect_equal(efa_spss$settings$abs_eigen,  TRUE)
  expect_equal(efa_equa$settings$abs_eigen, TRUE)
  expect_equal(efa_quart$settings$abs_eigen,  TRUE)
  expect_equal(efa_none$settings$abs_eigen, FALSE)

  expect_equal(efa_psych$settings$normalize, TRUE)
  expect_equal(efa_spss$settings$normalize, TRUE)
  expect_equal(efa_equa$settings$normalize, TRUE)
  expect_equal(efa_quart$settings$normalize, TRUE)
  expect_equal(efa_none$settings$normalize, TRUE)

  expect_equal(efa_psych$settings$P_type, "unnorm")
  expect_equal(efa_spss$settings$P_type, "norm")
  expect_equal(efa_none$settings$P_type, "unnorm")

  expect_equal(efa_psych$settings$precision, 1e-05)
  expect_equal(efa_spss$settings$precision, 1e-05)
  expect_equal(efa_equa$settings$precision, 1e-05)
  expect_equal(efa_quart$settings$precision, 1e-05)
  expect_equal(efa_none$settings$precision, 1e-05)

  expect_equal(efa_psych$settings$order_type, "eigen")
  expect_equal(efa_spss$settings$order_type, "ss_factors")
  expect_equal(efa_equa$settings$order_type, "eigen")
  expect_equal(efa_quart$settings$order_type, "eigen")
  expect_equal(efa_none$settings$order_type, "eigen")

  expect_equal(efa_psych$settings$varimax_type, "svd")
  expect_equal(efa_spss$settings$varimax_type, "kaiser")
  expect_equal(efa_none$settings$varimax_type, "svd")

  expect_equal(efa_psych$settings$k, 4)
  expect_equal(efa_spss$settings$k, 4)
  expect_equal(efa_none$settings$k, 3)

  expect_equal(efa_ml$settings$start_method, "psych")
})

test_that("settings record the seed and how the correlations were obtained", {
  cormat <- test_models$baseline$cormat

  # a correlation matrix consumes no correlation method, so the resolved record is NA
  # while the requested value is kept unchanged
  s_cor <- efa_fit(cormat, n_factors = 3, N = 500)$settings
  expect_identical(s_cor$input_type, "correlation")
  expect_identical(s_cor$cor_method_used, NA_character_)
  expect_identical(s_cor$cor_method, "pearson")
  expect_null(s_cor$seed)

  # a method requested but never run is still echoed as requested, and still not
  # reported as the one that ran
  s_tetra <- efa_fit(cormat, n_factors = 3, N = 500, cor_method = "tetra")$settings
  expect_identical(s_tetra$cor_method, "tetra")
  expect_identical(s_tetra$cor_method_used, NA_character_)

  # raw data records the method that actually produced the correlations, and the seed
  s_raw <- suppressMessages(
    efa_fit(GRiPS_raw, n_factors = 1, cor_method = "spearman", seed = 42))$settings
  expect_identical(s_raw$input_type, "raw")
  expect_identical(s_raw$cor_method_used, "spearman")
  expect_equal(s_raw$seed, 42)
})

test_that("factor analyses are performed correctly", {
  expect_equal(mean(abs(efa_compare(efa_paf_zero$rot_loadings,
                                   population_models$loadings$baseline,
                                   plot = FALSE)$diff)), 0, tolerance = .01)
  expect_equal(mean(abs(efa_compare(efa_ml_zero$rot_loadings,
                                   population_models$loadings$baseline,
                                   plot = FALSE)$diff)), 0, tolerance = .01)
  expect_equal(mean(abs(efa_compare(efa_uls_zero$rot_loadings,
                                   population_models$loadings$baseline,
                                   plot = FALSE)$diff)), 0, tolerance = .01)
  # expect_equal(mean(abs(efa_compare(efa_paf_moderate$rot_loadings,
  #                                  population_models$loadings$baseline,
  #                                  plot = FALSE)$diff)), 0, tolerance = .01)
  expect_equal(mean(abs(efa_compare(efa_ml_moderate$rot_loadings,
                                   population_models$loadings$baseline,
                                   plot = FALSE)$diff)), 0, tolerance = .01)
  expect_equal(mean(abs(efa_compare(efa_uls_moderate$rot_loadings,
                                   population_models$loadings$baseline,
                                   plot = FALSE)$diff)), 0, tolerance = .01)
})

test_that("errors are thrown correctly", {
  expect_error(EFA(1:5), class = "efa_input_not_matrix")
  expect_error(EFA(cor_nposdef, n_factors = 1, N = 10, type = "SPSS"), class = "efa_cor_not_posdef")
  expect_message(EFA(GRiPS_raw, n_factors = 1), class = "efa_cor_from_data")
  expect_warning(EFA(GRiPS_raw, N = 20, n_factors = 1), class = "efa_n_from_data")
  expect_error(EFA(sing_raw, n_factors = 1), class = "efa_cor_singular")
  expect_error(EFA(sing_cor, N = sing_N, n_factors = 1), class = "efa_cor_singular")
  # a criterion of 1 is not a convergence tolerance. EFA() repacks it into the estimation
  # control, which rejects it up front -- as it already did for any larger value -- so the
  # invalid call is caught at the argument rather than deep inside the fit.
  expect_error(EFA(test_models$baseline$cormat, n_factors = 3, N = 500, method = "PAF", criterion = 1),
               class = "efa_control_input")
  # Identification is a property of the dimensions, not of the values, but fitting a
  # 10 x 3 sample does depend on the draw: whether the extraction Heywoods or runs out of
  # iterations varies from one realization to the next. Seed the fixture and let those two
  # incidental conditions through, so the assertion stays about identification alone.
  set.seed(1)
  dat_ident <- matrix(rnorm(30), ncol = 3)
  expect_warning(EFA(dat_ident, n_factors = 2), class = "efa_underidentified")
  expect_warning(
    suppressWarnings(
      EFA(dat_ident, n_factors = 1),
      classes = c("efa_heywood", "efa_paf_nonconvergence")
    ),
    class = "efa_just_identified"
  )
  expect_warning(EFA(test_models$baseline$cormat, n_factors = 3, method = "ML"), class = "efa_fit_na_n")
  expect_warning(EFA(test_models$baseline$cormat, n_factors = 3, method = "ULS"), class = "efa_fit_na_n")
  expect_warning(
    suppressWarnings(
      EFA(cor_nposdef, n_factors = 1, N = 10),
      classes = c("efa_just_identified", "efa_heywood")
    ),
    class = "efa_cor_smoothed"
  )
  expect_warning(EFA(test_models$baseline$cormat, n_factors = 3, method = "ML", N = 500, type = "SPSS"), class = "efa_spss_method_untested")
  expect_warning(EFA(test_models$baseline$cormat, n_factors = 3, method = "ULS", N = 500, type = "SPSS"), class = "efa_spss_method_untested")
  expect_warning(EFA(test_models$baseline$cormat, n_factors = 3, rotation = "oblimin", N = 500, type = "SPSS"), class = "efa_spss_rotation_untested")
  expect_warning(EFA(test_models$baseline$cormat, n_factors = 3, rotation = "quartimax", N = 500, type = "SPSS"), class = "efa_spss_rotation_untested")
})

test_that("an unusable N is rejected or reported, and never left silent", {
  R <- test_models$baseline$cormat

  # N = 0 is an impossible sample size, not an unknown one: it goes through the classed
  # argument-validation path. NA remains the supported way of saying "sample size unknown",
  # and still warns.
  expect_error(efa_fit(R, n_factors = 3, N = 0, estimator = "ML"),
               class = "efa_invalid_argument")
  expect_warning(efa_fit(R, n_factors = 3, N = NA, estimator = "ML"),
                 class = "efa_fit_na_n")

  # A valid but very small N relative to p turns Bartlett's small-sample multiplier
  # non-positive, so the chi-square block is undefined. That was previously silent and
  # indistinguishable from a PAF fit; it now warns exactly once, and the residual summaries
  # and df it leaves behind are still populated.
  expect_warning(fit_small <- efa_fit(R, n_factors = 3, N = 8, estimator = "ML"),
                 class = "efa_fit_indices_undefined")
  expect_true(is.na(fit_small$fit_indices$chi))
  expect_true(is.na(fit_small$fit_indices$CFI))
  expect_true(is.na(fit_small$fit_indices$RMSEA))
  expect_true(is.finite(fit_small$fit_indices$SRMR))
  expect_equal(fit_small$fit_indices$df, .efa_df(ncol(R), 3))

  # PAF reports no chi-square at any N, so the small-sample warning would be noise there.
  expect_no_warning(efa_fit(R, n_factors = 3, N = 8, estimator = "PAF"),
                    class = "efa_fit_indices_undefined")

  # An underidentified model leaves the same block undefined, and its warning must not claim
  # that no fit indices at all were computed: the residual summaries are returned, and are
  # exactly the numbers a reader would otherwise misread as a flawless solution.
  expect_warning(
    fit_ui <- suppressWarnings(efa_fit(R, n_factors = 17, N = 500, estimator = "ML"),
                               classes = c("efa_heywood", "efa_nonconvergence")),
    class = "efa_underidentified"
  )
  expect_lt(fit_ui$fit_indices$df, 0)
  expect_true(is.na(fit_ui$fit_indices$chi))
  expect_true(is.finite(fit_ui$fit_indices$SRMR))
  expect_true(is.finite(fit_ui$fit_indices$CAF))
})

test_that("a singular correlation matrix is rejected for every type except psych + PAF", {
  # The psych preset relaxes the singularity check because psych::smc() falls back
  # to a pseudo-inverse for the PAF starting communalities. ML and ULS have no such
  # fallback -- they solve R for their starting values, and the ML discrepancy needs
  # log|R| -- so the check must stay on for them under every preset, and the abort
  # must be the package's classed one rather than a raw LAPACK error.
  for (ty in c("EFAtools", "psych", "SPSS")) {
    for (est in c("PAF", "ML", "ULS")) {
      if (ty == "psych" && est == "PAF") next
      expect_error(
        efa_fit(sing_cor, n_factors = 1, N = sing_N, estimator = est,
                estimate_control = estimate_control(type = ty)),
        class = "efa_cor_singular"
      )
    }
  }

  # psych + PAF keeps the documented pseudo-inverse behaviour and still fits.
  expect_no_error(
    suppressWarnings(
      efa_fit(sing_cor, n_factors = 1, N = sing_N, estimator = "PAF",
              estimate_control = estimate_control(type = "psych"))
    )
  )
})

test_that("EFA rejects n_factors >= number of variables", {
  dat <- matrix(rnorm(30), ncol = 3)  # 3 variables

  # n_factors == n_variables: out of bounds for ML, degenerate for ULS/PAF
  expect_error(EFA(dat, n_factors = 3, method = "ML"),  class = "efa_too_many_factors")
  expect_error(EFA(dat, n_factors = 3, method = "ULS"), class = "efa_too_many_factors")
  expect_error(EFA(dat, n_factors = 3, method = "PAF"), class = "efa_too_many_factors")

  # n_factors > n_variables
  expect_error(EFA(dat, n_factors = 4, method = "ML"),  class = "efa_too_many_factors")

  # the bootstrap path reuses n_factors and is guarded before any resampling runs
  expect_error(EFA(dat, n_factors = 3, method = "ML", se = "np-boot"),
               class = "efa_too_many_factors")
})

test_that("EFA rejects n_factors below 1 and a missing n_factors", {
  cm <- test_models$baseline$cormat

  # Several factor retention criteria legitimately return 0 ("no factor worth
  # extracting"), so the retain-then-fit idiom can hand efa_fit() a zero. It must
  # be rejected in R with a classed condition, not by the C++ bound check.
  for (est in c("PAF", "ML", "ULS")) {
    expect_error(efa_fit(cm, n_factors = 0, N = 500, estimator = est),
                 class = "efa_too_few_factors")
  }
  expect_error(EFA(cm, n_factors = 0, N = 500), class = "efa_too_few_factors")

  # Omitting n_factors points at the function whose job is to answer it.
  expect_error(efa_fit(cm, N = 500), class = "efa_missing_n_factors")
  expect_error(EFA(cm, N = 500), class = "efa_missing_n_factors")
})

test_that("print.efa output is stable (PAF, promax)", {
  local_reproducible_output()

  expect_snapshot(print(efa_psych), transform = scrub_num)
})

test_that("print.efa output is stable (ML, promax)", {
  local_reproducible_output()

  expect_snapshot(print(efa_ml_moderate), transform = scrub_num)
})

test_that("format.efa is the source of truth and honours the colour state", {
  # print() is exactly cat(format(x), sep = "\n"), so the two agree line for line.
  expect_identical(utils::capture.output(print(efa_psych)), format(efa_psych))

  old <- options(cli.num_colors = 256)
  on.exit(options(old), add = TRUE)

  # With colours on the report embeds ANSI ...
  expect_true(any(grepl("\033", format(efa_psych), fixed = TRUE)))

  # ... and with colours off it is plain.
  options(cli.num_colors = 1)
  expect_false(any(grepl("\033", format(efa_psych), fixed = TRUE)))
})

test_that("the model header follows the console width", {
  # The header names the analysis and is emitted verbatim, so a "setting = 'value'" token is
  # never split; it is packed to the console around those tokens instead.
  header <- function(w, colors = 1) {
    out <- withr::with_options(list(cli.width = w, cli.num_colors = colors),
                               cli::ansi_strip(format(efa_psych)))
    wrapped_item(out, "^EFA performed")
  }

  for (w in c(120L, 80L, 60L)) {
    lines <- header(w)
    expect_true(all(cli::ansi_nchar(lines, type = "width") <= w))
    # wrapping moves the line break only: the same sentence, with both tokens intact
    expect_identical(paste(trimws(lines), collapse = " "),
                     "EFA performed with estimator = 'PAF' and rotation = 'promax'.")
    expect_false(any(grepl("= '[^']*$", lines)))
  }

  # a narrow console does split the line, and does so identically with colours on
  expect_length(header(60L), 2L)
  expect_identical(header(60L), header(60L, colors = 256))
})

test_that("summary.efa output is stable (PAF, promax)", {
  local_reproducible_output()

  expect_snapshot(print(summary(efa_psych)), transform = scrub_num)
})

test_that("summary.efa output is stable (ML, promax)", {
  local_reproducible_output()

  expect_snapshot(print(summary(efa_ml_moderate)), transform = scrub_num)
})

test_that("format.summary.efa is the source of truth and honours the colour state", {
  s <- summary(efa_psych)

  # print() is exactly cat(format(x), sep = "\n"), so the two agree line for line.
  expect_identical(utils::capture.output(print(s)), format(s))

  old <- options(cli.num_colors = 256)
  on.exit(options(old), add = TRUE)

  # With colours on the report embeds ANSI ...
  expect_true(any(grepl("\033", format(s), fixed = TRUE)))

  # ... and with colours off it is plain.
  options(cli.num_colors = 1)
  expect_false(any(grepl("\033", format(s), fixed = TRUE)))
})
test_that("SRMR is the printed residual fit index while RMSR remains returned", {
  residuals <- matrix(c(0, .10, -.20,
                        .10, 0, .30,
                        -.20, .30, 0), nrow = 3, byrow = TRUE)
  p <- nrow(residuals)
  rmsr <- .rmsr(residuals)
  srmr <- sqrt(sum(residuals[upper.tri(residuals)]^2) / (p * (p + 1) / 2))

  expect_equal(srmr, rmsr * sqrt((p - 1) / (p + 1)))
  expect_lt(srmr, rmsr)

  expect_true(all(c("RMSR", "SRMR") %in% names(efa_ml$fit_indices)))
  out <- cli::ansi_strip(format(efa_ml))
  summary_out <- cli::ansi_strip(format(summary(efa_ml)))
  expect_false(any(grepl("^RMSR\\b", out)))
  expect_false(any(grepl("^RMSR\\b", summary_out)))
  expect_true(any(grepl("^SRMR:", out)))
  expect_true(any(grepl("^SRMR:", summary_out)))
})

test_that("print/summary.efa omit the inapplicable tables for a rotated single factor", {
  local_reproducible_output()

  efa_1fac <- suppressWarnings(
    EFA(test_models$baseline$cormat, n_factors = 1, N = 500, method = "PAF",
        rotation = "promax")
  )

  # A single factor cannot be rotated, so the rotation-only outputs are absent and
  # their print sections are skipped (rather than rendering a stray NA).
  expect_null(efa_1fac$Phi)
  expect_null(efa_1fac$vars_accounted_rot)

  # The variance table also drops the three cumulative/common-variance rows here: with one
  # factor they would only repeat the two above them. This shape is part of the documented
  # return, because user code indexing a row by name has to allow for it.
  expect_identical(rownames(efa_1fac$vars_accounted), c("SS loadings", "Prop Tot Var"))

  expect_snapshot(print(efa_1fac), transform = scrub_num)
  expect_snapshot(print(summary(efa_1fac)), transform = scrub_num)
})

test_that("EFA records rotation diagnostics for multistart rotations", {
  # A native gradient-projection rotation records the per-start diagnostics.
  efa_geo <- suppressWarnings(
    EFA(test_models$baseline$cormat, n_factors = 3, N = 500, rotation = "geominQ")
  )
  rd <- efa_geo$settings$rotation_diagnostics
  expect_type(rd, "list")
  expect_named(rd, c("converged", "n_starts_total", "n_optimized", "n_converged",
                     "n_distinct_minima", "criterion_spread", "criterion_best"))
  # the convergence status of the solution that is actually returned, so a saved fit can be
  # audited for it after the fact rather than only through the transient warning
  expect_true(rd$converged)
  # 100 random starts plus the rational (identity) start; only a few of them are ever
  # optimized under the screen-and-triage strategy.
  expect_equal(rd$n_starts_total, 101L)
  expect_lte(rd$n_optimized, rd$n_starts_total)
  # The distinct-optima count is taken over the converged starts only, so it is at least
  # one (the rational start converges here) and never exceeds the converged count, which
  # in turn can never exceed the number of optimized starts.
  expect_gte(rd$n_converged, 1L)
  expect_lte(rd$n_converged, rd$n_optimized)
  expect_gte(rd$n_distinct_minima, 1L)
  expect_lte(rd$n_distinct_minima, rd$n_converged)
  expect_true(is.finite(rd$criterion_best))
  expect_gte(rd$criterion_spread, 0)

  # simplimax also records diagnostics. It runs full multistart, so every requested random
  # start is optimized and the two counts coincide.
  efa_simp <- suppressWarnings(
    EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
        rotation = "simplimax", randomStarts = 20)
  )
  rd_simp <- efa_simp$settings$rotation_diagnostics
  expect_equal(rd_simp$n_starts_total, 21L)
  expect_equal(rd_simp$n_optimized, 21L)
  expect_true(is.finite(rd_simp$criterion_best))

  # Rotations that do not use random starts carry no diagnostics.
  expect_null(efa_psych$settings$rotation_diagnostics)   # promax
  expect_null(efa_cor$settings$rotation_diagnostics)     # rotation = "none"

  # summary() surfaces the diagnostic only for the multistart rotations.
  local_reproducible_output()
  geo_summary <- utils::capture.output(print(summary(efa_geo)))
  expect_true(any(grepl("Rotation local optima", geo_summary, fixed = TRUE)))
  # The line reports the starts that were actually optimized alongside the total, so it
  # cannot overstate the evidence by counting screened-and-discarded starts as optima.
  expect_true(any(grepl(
    paste0(rd$n_distinct_minima, " distinct from ", rd$n_optimized,
           " of ", rd$n_starts_total, " starts"),
    geo_summary, fixed = TRUE)))
  promax_summary <- utils::capture.output(print(summary(efa_psych)))
  expect_false(any(grepl("Rotation local optima", promax_summary, fixed = TRUE)))
})

test_that("the rotation diagnostics record the returned solution's own convergence", {
  # The returned solution is the start with the lowest criterion value, which need not be a
  # converged one: a start that stopped short of the tolerance at a lower value beats a
  # converged start at a higher one. `converged` therefore has to be carried through from
  # the winning start rather than read off the aggregate count, or such a solution would be
  # presented as converged because its rivals were. Exercised on the helper so the
  # disagreement is deterministic instead of depending on which start a seed happens to win.
  rd <- .rotation_diagnostics(all_values = c(0.5, 0.1), all_converged = c(1L, 0L),
                              randomStarts = 1, converged = FALSE)
  expect_false(rd$converged)
  expect_equal(rd$n_converged, 1L)
  expect_equal(rd$criterion_best, 0.1)

  # end to end: with a one-iteration budget no start can converge, and the object records
  # that instead of leaving it in the transient warning only
  capped <- suppressWarnings(
    EFA(test_models$baseline$cormat, n_factors = 3, N = 500, rotation = "geominQ",
        maxit = 1)
  )
  expect_false(capped$settings$rotation_diagnostics$converged)
})

test_that("residuals.efa is a pure extractor", {
  # returns the residual matrix with no printing side effect
  expect_output(residuals(efa_psych), NA)
  expect_identical(residuals(efa_psych), efa_psych$residuals)
  expect_identical(residuals(efa_psych, "raw"), efa_psych$residuals)

  # standardized residuals need bootstrap SEs, which efa_psych does not have
  expect_error(residuals(efa_psych, "standardized"),
               class = "efa_no_standardized_residuals")

  # an unknown type is rejected
  expect_error(residuals(efa_psych, "bogus"))
})

# A near-unit first loading pins item1's uniqueness at the ML/ULS estimation
# boundary, giving a genuine improper solution with every communality still below 1.
boundary_heywood_cormat <- function() {
  L <- c(.999, .7, .6, .5)
  R <- L %*% t(L)
  diag(R) <- 1
  dimnames(R) <- list(paste0("item", seq_along(L)), paste0("item", seq_along(L)))
  R
}

test_that("Heywood cases are detected, warned, and recorded", {
  # Over-extracting the baseline model with PAF yields a Heywood case (V14).
  expect_warning(
    suppressWarnings(
      EFA(test_models$baseline$cormat, 6, N = 500, method = "PAF"),
      classes = "efa_paf_nonconvergence"
    ),
    class = "efa_heywood"
  )

  efa_hey <- suppressWarnings(EFA(test_models$baseline$cormat, 6, N = 500,
                                  method = "PAF"))
  expect_gt(length(efa_hey$heywood), 0)
  expect_true(all(efa_hey$h2[efa_hey$heywood] >= 1))

  # A proper solution records no Heywood cases and emits no Heywood warning.
  efa_clean <- suppressWarnings(EFA(test_models$baseline$cormat, 3, N = 500,
                                    method = "ML"))
  expect_length(efa_clean$heywood, 0)
  # ... and it is still a named integer vector when empty, as documented
  expect_named(efa_clean$heywood)

  # ML/ULS constrain the uniquenesses to a lower bound, so an improper solution
  # keeps the communality just below 1 and instead pins a uniqueness at that
  # bound. A near-unit first loading forces this boundary case; it must be flagged
  # for both estimators even though no communality reaches 1.
  R_bnd <- boundary_heywood_cormat()

  for (m in c("ML", "ULS")) {
    expect_warning(EFA(R_bnd, 1, N = 200, method = m), class = "efa_heywood")
    efa_bnd <- suppressWarnings(EFA(R_bnd, 1, N = 200, method = m))
    expect_gt(length(efa_bnd$heywood), 0)
    expect_true(all(efa_bnd$h2[efa_bnd$heywood] < 1))
    # heywood is documented as a *named* integer vector: the optimiser-based
    # fitters return the uniquenesses as a one-column matrix, whose dim attribute
    # would otherwise displace the variable names when the two criteria are
    # combined, leaving the warning below with nothing to name.
    expect_named(efa_bnd$heywood, "item1")
  }

  # The same matrix with a clearly proper structure is not flagged.
  L_ok <- c(.9, .7, .6, .5)
  R_ok <- L_ok %*% t(L_ok)
  diag(R_ok) <- 1
  expect_length(suppressWarnings(EFA(R_ok, 1, N = 200, method = "ML"))$heywood, 0)
})

test_that("the Heywood warning names the affected variables", {
  local_reproducible_output()

  R_bnd <- boundary_heywood_cormat()

  # ML and ULS reach the Heywood case through a uniqueness pinned at the
  # estimation boundary; the message must name the variable rather than report a
  # bare count.
  expect_snapshot({
    for (est in c("ML", "ULS")) {
      cat("--", est, "--\n")
      withCallingHandlers(
        efa_fit(R_bnd, 1, N = 200, estimator = est),
        efa_heywood = function(w) cat(conditionMessage(w), "\n"),
        warning = function(w) invokeRestart("muffleWarning")
      )
    }
  })
})

test_that("point-estimate non-convergence is surfaced for the iterative estimators", {
  R <- test_models$baseline$cormat

  # Force the convergence code via a mock so the assertions do not depend on the
  # well-behaved baseline fits actually converging. EFA() must warn on a non-zero
  # code and stay silent on a zero code, for ML, ULS, and DWLS. heywood is cleared
  # so the only candidate warning is the convergence one. Each mock is scoped to its
  # own local() frame.
  for (m in c("ML", "ULS")) {
    local({
      fit <- suppressWarnings(.estimate_model(R, method = m, n_factors = 3, N = 500,
                                              start_method = "psych"))
      fit$heywood <- integer(0)
      testthat::local_mocked_bindings(.estimate_model = function(...) fit)

      fit$convergence <- 1L
      expect_warning(EFA(R, 3, N = 500, method = m), class = "efa_nonconvergence")

      # A zero code must not warn (the mock makes this deterministic).
      fit$convergence <- 0L
      conv_classes <- character(0)
      withCallingHandlers(
        EFA(R, 3, N = 500, method = m),
        warning = function(w) {
          conv_classes <<- c(conv_classes, class(w))
          invokeRestart("muffleWarning")
        }
      )
      expect_false("efa_nonconvergence" %in% conv_classes)
    })
  }

  # DWLS needs raw ordinal data with cor_method = "poly". The mocked fitter makes
  # the real estimation irrelevant, so a minimal fit carrying only the fields EFA()
  # reads after the fit suffices (the polychoric matrix is still built inside EFA()).
  local({
    xo <- DOSPERT_raw[, 1:6]
    testthat::local_mocked_bindings(
      .estimate_model = function(...) list(convergence = 1L, heywood = integer(0))
    )
    expect_warning(EFA(xo, 2, method = "DWLS", cor_method = "poly"),
                   class = "efa_nonconvergence")
  })

  # PAF surfaces non-convergence through its own warning class, never the
  # optimiser-estimator one.
  paf_classes <- character(0)
  withCallingHandlers(
    EFA(R, 3, N = 500, method = "PAF", type = "none", init_comm = "smc",
        criterion = 1e-3, criterion_type = "sum", abs_eigen = TRUE, max_iter = 1),
    warning = function(w) {
      paf_classes <<- c(paf_classes, class(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true("efa_paf_nonconvergence" %in% paf_classes)
  expect_false("efa_nonconvergence" %in% paf_classes)
})

test_that("the non-convergence explanation is branched on the convergence code", {
  # Code 1 is the optimiser's iteration limit; 51 and 52 come from L-BFGS-B's line
  # search and mean something different, so they must not be explained as an iteration
  # limit. Anything else falls back to the generic explanation.
  cap <- .nonconvergence_cause(1L)
  line_search <- .nonconvergence_cause(52L)
  generic <- .nonconvergence_cause(10L)

  for (msg in list(cap, line_search, generic)) {
    expect_type(msg, "character")
    expect_length(msg, 1L)
    expect_true(nzchar(msg))
  }

  # the two L-BFGS-B line-search codes share one explanation, distinct from the others
  expect_identical(.nonconvergence_cause(51L), line_search)
  expect_false(identical(cap, line_search))
  expect_false(identical(cap, generic))
  expect_false(identical(line_search, generic))

  # the branch is what the warning carries: a mocked code-52 fit must not be explained
  # with the iteration-limit text. The message is compared with its line breaks
  # collapsed, so console wrapping cannot decide the assertion.
  testthat::local_reproducible_output()
  squash <- function(x) gsub("\\s+", " ", paste(x, collapse = " "))

  R <- test_models$baseline$cormat
  fit <- suppressWarnings(.estimate_model(R, method = "ULS", n_factors = 3, N = 500))
  fit$heywood <- integer(0)
  fit$convergence <- 52L
  testthat::local_mocked_bindings(.estimate_model = function(...) fit)

  w <- expect_warning(EFA(R, 3, N = 500, method = "ULS"),
                      class = "efa_nonconvergence")
  msg <- squash(conditionMessage(w))
  expect_true(grepl(squash(line_search), msg, fixed = TRUE))
  expect_false(grepl(squash(cap), msg, fixed = TRUE))
})

test_that("the print banner reports optimiser non-convergence", {
  testthat::local_reproducible_output()
  R <- test_models$baseline$cormat
  efa <- suppressWarnings(EFA(R, 3, N = 500, method = "ML"))

  # Set the code explicitly so the assertions do not depend on the fixture's actual
  # convergence. A zero code shows no banner; a non-zero code shows the generic
  # optimiser banner for ML/ULS/DWLS.
  efa$convergence <- 0L
  expect_false(any(grepl("The optimiser did not converge", format(efa), fixed = TRUE)))
  efa$convergence <- 1L
  expect_true(any(grepl("The optimiser did not converge", format(efa), fixed = TRUE)))

  # PAF keeps its iteration-specific banner.
  paf <- suppressWarnings(
    EFA(R, 3, N = 500, method = "PAF", type = "none", init_comm = "smc",
        criterion = 1e-3, criterion_type = "sum", abs_eigen = TRUE, max_iter = 1)
  )
  expect_true(any(grepl("Maximum number of iterations reached without convergence",
                        format(paf), fixed = TRUE)))

  # An object carrying no convergence code (e.g. one saved by an older version) still
  # shows the banner via the iteration-count fallback (iter >= max_iter).
  paf$convergence <- NULL
  expect_true(any(grepl("Maximum number of iterations reached without convergence",
                        format(paf), fixed = TRUE)))
})

rm(efa_cor, efa_raw, efa_psych, efa_spss, efa_ml, efa_uls, efa_equa, efa_quart,
   efa_none, cormat_zero, cormat_moderate, efa_paf_zero, efa_ml_zero, efa_uls_zero,
   efa_paf_moderate, efa_ml_moderate, efa_uls_moderate)

test_that("print.efa argument validators raise classed conditions", {
  efa <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500)

  expect_error(print(efa, cutoff = -1), class = "efa_print_invalid_cutoff")
  expect_error(print(efa, digits = 1.5), class = "efa_print_invalid_digits")
  expect_error(print(efa, max_name_length = 0),
               class = "efa_print_invalid_max_name_length")
  expect_error(print(efa, diagnostics_top_n = -1),
               class = "efa_print_invalid_diagnostics_top_n")
  expect_error(print(efa, residual_cutoff = -1),
               class = "efa_print_invalid_residual_cutoff")
  expect_error(print(efa, residual_top_n = 0),
               class = "efa_print_invalid_residual_top_n")
  expect_error(print(efa, show_structure = NA),
               class = "efa_print_invalid_show_structure")
  expect_error(print(efa, show_loading_legend = NA),
               class = "efa_print_invalid_show_loading_legend")
  expect_error(print(efa, cross_loading_cutoff = -1),
               class = "efa_print_invalid_cross_loading_cutoff")
  expect_error(print(efa, min_primary_gap = -1),
               class = "efa_print_invalid_min_primary_gap")
  expect_error(print(efa, min_salient_per_factor = 0),
               class = "efa_print_invalid_min_salient_per_factor")
  expect_error(print(efa, max_factors_per_block = 0),
               class = "efa_print_invalid_max_factors_per_block")
  expect_error(print(efa, show_mi_diagnostics = NA),
               class = "efa_print_invalid_show_mi_diagnostics")
})


