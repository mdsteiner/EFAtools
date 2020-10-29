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
  expect_is(efa_cor, "EFA")
  expect_is(efa_raw, "EFA")
  expect_is(efa_psych, "EFA")
  expect_is(efa_spss, "EFA")
  expect_is(efa_ml, "EFA")
  expect_is(efa_uls, "EFA")
  expect_is(efa_equa, "EFA")
  expect_is(efa_quart, "EFA")
  expect_is(efa_none, "EFA")

  expect_is(efa_cor$unrot_loadings, "LOADINGS")
  expect_is(efa_raw$unrot_loadings, "LOADINGS")
  expect_is(efa_psych$unrot_loadings, "LOADINGS")
  expect_is(efa_spss$unrot_loadings, "LOADINGS")
  expect_is(efa_ml$unrot_loadings, "LOADINGS")
  expect_is(efa_uls$unrot_loadings, "LOADINGS")
  expect_is(efa_equa$unrot_loadings, "LOADINGS")
  expect_is(efa_quart$unrot_loadings, "LOADINGS")
  expect_is(efa_none$unrot_loadings, "LOADINGS")

  expect_is(efa_psych$rot_loadings, "LOADINGS")
  expect_is(efa_spss$rot_loadings, "LOADINGS")
  expect_is(efa_equa$rot_loadings, "LOADINGS")
  expect_is(efa_quart$rot_loadings, "LOADINGS")
  expect_is(efa_none$rot_loadings, "LOADINGS")

  expect_named(efa_cor, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                          "final_eigen", "iter", "convergence", "unrot_loadings",
                          "vars_accounted", "fit_indices", "settings"))
  expect_named(efa_raw, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                          "final_eigen", "iter", "convergence", "unrot_loadings",
                          "vars_accounted", "fit_indices", "settings"))
  expect_named(efa_psych, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                            "final_eigen", "iter", "convergence", "unrot_loadings",
                            "vars_accounted", "fit_indices", "rot_loadings",
                            "Phi", "Structure", "rotmat", "vars_accounted_rot",
                            "settings"))
  expect_named(efa_spss, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                           "final_eigen", "iter", "convergence", "unrot_loadings",
                           "vars_accounted", "fit_indices", "rot_loadings",
                           "Phi", "Structure", "rotmat", "vars_accounted_rot",
                           "settings"))
  expect_named(efa_ml, c("orig_R", "h2", "orig_eigen", "final_eigen", "iter",
                         "convergence", "unrot_loadings", "vars_accounted",
                         "fit_indices", "settings"))
  expect_named(efa_uls, c("orig_R", "h2", "orig_eigen", "final_eigen", "iter",
                         "convergence", "unrot_loadings", "vars_accounted",
                         "fit_indices", "settings"))
  expect_named(efa_equa, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                           "final_eigen", "iter", "convergence", "unrot_loadings",
                           "vars_accounted", "fit_indices", "rot_loadings",
                           "rotmat", "vars_accounted_rot",
                           "settings"))
  expect_named(efa_quart, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                            "final_eigen", "iter", "convergence", "unrot_loadings",
                            "vars_accounted", "fit_indices", "rot_loadings",
                            "Phi", "Structure", "rotmat", "vars_accounted_rot",
                            "settings"))
  expect_named(efa_none, c("orig_R", "h2_init", "h2", "orig_eigen", "init_eigen",
                           "final_eigen", "iter", "convergence", "unrot_loadings",
                           "vars_accounted", "fit_indices", "rot_loadings",
                           "Phi", "Structure", "rotmat", "vars_accounted_rot",
                           "settings"))
})

test_that("settings are returned correctly", {
  expect_named(efa_cor$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen"))
  expect_named(efa_raw$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen"))
  expect_named(efa_psych$settings, c("method", "rotation", "type", "n_factors",
                                     "N", "use", "cor_method", "max_iter",
                                     "init_comm", "criterion", "criterion_type",
                                     "abs_eigen", "normalize", "P_type", "precision",
                                     "order_type", "varimax_type", "k"))
  expect_named(efa_spss$settings, c("method", "rotation", "type", "n_factors",
                                    "N", "use", "cor_method", "max_iter",
                                    "init_comm", "criterion", "criterion_type",
                                    "abs_eigen", "normalize", "P_type", "precision",
                                    "order_type", "varimax_type", "k"))
  expect_named(efa_ml$settings, c("method", "rotation", "type", "n_factors",
                                    "N", "use", "cor_method", "start_method"))
  expect_named(efa_uls$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method"))
  expect_named(efa_equa$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "normalize", "precision",
                                   "order_type"))
  expect_named(efa_quart$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "normalize", "precision",
                                   "order_type", "k"))
  expect_named(efa_none$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "normalize", "P_type", "precision",
                                   "order_type", "varimax_type", "k"))

  expect_equal(efa_cor$settings$method, "PAF")
  expect_equal(efa_raw$settings$method, "PAF")
  expect_equal(efa_psych$settings$method, "PAF")
  expect_equal(efa_spss$settings$method, "PAF")
  expect_equal(efa_ml$settings$method, "ML")
  expect_equal(efa_uls$settings$method, "ULS")
  expect_equal(efa_equa$settings$method, "PAF")
  expect_equal(efa_quart$settings$method, "PAF")
  expect_equal(efa_none$settings$method, "PAF")

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

test_that("factor analyses are performed correctly", {
  expect_equal(mean(abs(COMPARE(efa_paf_zero$rot_loadings,
                               population_models$loadings$baseline,
                               plot = FALSE)$diff)), 0, tolerance = .01)
  expect_equal(mean(abs(COMPARE(efa_ml_zero$rot_loadings,
                               population_models$loadings$baseline,
                               plot = FALSE)$diff)), 0, tolerance = .01)
  expect_equal(mean(abs(COMPARE(efa_uls_zero$rot_loadings,
                               population_models$loadings$baseline,
                               plot = FALSE)$diff)), 0, tolerance = .01)
  # expect_equal(mean(abs(COMPARE(efa_paf_moderate$rot_loadings,
  #                              population_models$loadings$baseline,
  #                              plot = FALSE)$diff)), 0, tolerance = .01)
  expect_equal(mean(abs(COMPARE(efa_ml_moderate$rot_loadings,
                               population_models$loadings$baseline,
                               plot = FALSE)$diff)), 0, tolerance = .01)
  expect_equal(mean(abs(COMPARE(efa_uls_moderate$rot_loadings,
                               population_models$loadings$baseline,
                               plot = FALSE)$diff)), 0, tolerance = .01)
})

# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors are thrown correctly", {
  expect_error(EFA(1:5), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n")
  expect_message(EFA(GRiPS_raw, n_factors = 1), " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n")
  expect_warning(EFA(GRiPS_raw, N = 20, n_factors = 1), " 'N' was set and data entered. Taking N from data.\n")
  expect_error(EFA(dat_sing, n_factors = 1), " Correlation matrix is singular, no further analyses are performed\n")
  expect_error(EFA(cor_sing, N = 10, n_factors = 1), " Correlation matrix is singular, no further analyses are performed\n")
  expect_error(EFA(matrix(rnorm(30), ncol = 3), n_factors = 2), " The model is underidentified. Please enter a lower number of factors or use a larger number of indicators and try again.\n")
  expect_warning(EFA(matrix(rnorm(30), ncol = 3), n_factors = 1), " The model is just identified (df = 0). We suggest to try again with a lower number of factors or a larger number of indicators.\n", fixed = TRUE)
  expect_warning(EFA(test_models$baseline$cormat, n_factors = 3, method = "ML"), " Argument 'N' was NA, not all fit indices could be computed. To get all fit indices, either provide N or raw data.\n")
  expect_warning(EFA(test_models$baseline$cormat, n_factors = 3, method = "ULS"), " Argument 'N' was NA, not all fit indices could be computed. To get all fit indices, either provide N or raw data.\n")
  expect_warning(EFA(cor_nposdef, n_factors = 1, N = 10), "Matrix was not positive definite, smoothing was done")
})

rm(efa_cor, efa_raw, efa_psych, efa_spss, efa_ml, efa_uls, efa_equa, efa_quart,
   efa_none, cormat_zero, cormat_moderate, efa_paf_zero, efa_ml_zero, efa_uls_zero,
   efa_paf_moderate, efa_ml_moderate, efa_uls_moderate, x, y, z, dat_sing, cor_sing,
   cor_nposdef)


