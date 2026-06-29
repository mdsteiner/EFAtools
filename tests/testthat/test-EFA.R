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
  expect_s3_class(efa_cor, "EFA")
  expect_s3_class(efa_raw, "EFA")
  expect_s3_class(efa_psych, "EFA")
  expect_s3_class(efa_spss, "EFA")
  expect_s3_class(efa_ml, "EFA")
  expect_s3_class(efa_uls, "EFA")
  expect_s3_class(efa_equa, "EFA")
  expect_s3_class(efa_quart, "EFA")
  expect_s3_class(efa_none, "EFA")

  expect_s3_class(efa_cor$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_raw$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_psych$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_spss$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_ml$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_uls$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_equa$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_quart$unrot_loadings, "LOADINGS")
  expect_s3_class(efa_none$unrot_loadings, "LOADINGS")

  expect_s3_class(efa_psych$rot_loadings, "LOADINGS")
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
  expect_named(efa_cor$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "se", "b_boot", "ci", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen"))
  expect_named(efa_raw$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "se", "b_boot", "ci", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen"))
  expect_named(efa_psych$settings, c("method", "rotation", "type", "n_factors",
                                     "N", "use", "cor_method", "se", "b_boot", "ci", "max_iter",
                                     "init_comm", "criterion", "criterion_type",
                                     "abs_eigen", "normalize", "P_type", "precision",
                                     "order_type", "varimax_type", "k"))
  expect_named(efa_spss$settings, c("method", "rotation", "type", "n_factors",
                                    "N", "use", "cor_method", "se", "b_boot", "ci", "max_iter",
                                    "init_comm", "criterion", "criterion_type",
                                    "abs_eigen", "normalize", "P_type", "precision",
                                    "order_type", "varimax_type", "k"))
  expect_named(efa_ml$settings, c("method", "rotation", "type", "n_factors",
                                    "N", "use", "cor_method", "se", "b_boot", "ci", "start_method"))
  expect_named(efa_uls$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "se", "b_boot", "ci"))
  expect_named(efa_equa$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "se", "b_boot", "ci", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "normalize", "precision",
                                   "order_type", "randomStarts", "rotation_diagnostics"))
  expect_named(efa_quart$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "se", "b_boot", "ci", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen", "normalize", "precision",
                                   "order_type", "k", "randomStarts", "rotation_diagnostics"))
  expect_named(efa_none$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "se", "b_boot", "ci", "max_iter",
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
  expect_error(EFA(1:5), class = "efa_input_not_matrix")
  expect_error(EFA(cor_nposdef, n_factors = 1, N = 10, type = "SPSS"), class = "efa_cor_not_posdef")
  expect_message(EFA(GRiPS_raw, n_factors = 1), class = "efa_cor_from_data")
  expect_warning(EFA(GRiPS_raw, N = 20, n_factors = 1), class = "efa_n_from_data")
  expect_error(EFA(dat_sing, n_factors = 1), class = "efa_cor_singular")
  expect_error(EFA(cor_sing, N = 10, n_factors = 1), class = "efa_cor_singular")
  expect_error(EFA(test_models$baseline$cormat, n_factors = 3, N = 500, method = "PAF", criterion = 1),
               class = "efa_criterion_too_large")
  expect_warning(EFA(matrix(rnorm(30), ncol = 3), n_factors = 2), class = "efa_underidentified")
  expect_warning(EFA(matrix(rnorm(30), ncol = 3), n_factors = 1), class = "efa_just_identified")
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

test_that("print.EFA output is stable (PAF, promax)", {
  local_reproducible_output()

  expect_snapshot(print(efa_psych), transform = scrub_num)
})

test_that("print.EFA output is stable (ML, promax)", {
  local_reproducible_output()

  expect_snapshot(print(efa_ml_moderate), transform = scrub_num)
})

test_that("format.EFA is the source of truth and honours the colour state", {
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

test_that("summary.EFA output is stable (PAF, promax)", {
  local_reproducible_output()

  expect_snapshot(print(summary(efa_psych)), transform = scrub_num)
})

test_that("summary.EFA output is stable (ML, promax)", {
  local_reproducible_output()

  expect_snapshot(print(summary(efa_ml_moderate)), transform = scrub_num)
})

test_that("format.summary.EFA is the source of truth and honours the colour state", {
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

test_that("print/summary.EFA omit the inapplicable tables for a rotated single factor", {
  local_reproducible_output()

  efa_1fac <- suppressWarnings(
    EFA(test_models$baseline$cormat, n_factors = 1, N = 500, method = "PAF",
        rotation = "promax")
  )

  # A single factor cannot be rotated, so the rotation-only outputs are absent and
  # their print sections are skipped (rather than rendering a stray NA).
  expect_null(efa_1fac$Phi)
  expect_null(efa_1fac$vars_accounted_rot)

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
  expect_named(rd, c("n_starts", "n_converged", "n_distinct_minima",
                     "criterion_spread", "criterion_best"))
  expect_equal(rd$n_starts, 100)
  # The distinct-optima count is taken over the converged starts only, so it is at least
  # one (the rational start converges here) and never exceeds the converged count.
  expect_gte(rd$n_converged, 1L)
  expect_gte(rd$n_distinct_minima, 1L)
  expect_lte(rd$n_distinct_minima, rd$n_converged)
  expect_true(is.finite(rd$criterion_best))
  expect_gte(rd$criterion_spread, 0)

  # simplimax also records diagnostics; the requested random starts are stored verbatim.
  efa_simp <- suppressWarnings(
    EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
        rotation = "simplimax", randomStarts = 20)
  )
  rd_simp <- efa_simp$settings$rotation_diagnostics
  expect_equal(rd_simp$n_starts, 20)
  expect_true(is.finite(rd_simp$criterion_best))

  # Rotations that do not use random starts carry no diagnostics.
  expect_null(efa_psych$settings$rotation_diagnostics)   # promax
  expect_null(efa_cor$settings$rotation_diagnostics)     # rotation = "none"

  # summary() surfaces the diagnostic only for the multistart rotations.
  local_reproducible_output()
  geo_summary <- utils::capture.output(print(summary(efa_geo)))
  expect_true(any(grepl("Rotation local optima", geo_summary, fixed = TRUE)))
  promax_summary <- utils::capture.output(print(summary(efa_psych)))
  expect_false(any(grepl("Rotation local optima", promax_summary, fixed = TRUE)))
})

test_that("residuals.EFA is a pure extractor", {
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

  # ML/ULS constrain the uniquenesses to a lower bound, so an improper solution
  # keeps the communality just below 1 and instead pins a uniqueness at that
  # bound. A near-unit first loading forces this boundary case; it must be flagged
  # for both estimators even though no communality reaches 1.
  L_bnd <- c(.999, .7, .6, .5)
  R_bnd <- L_bnd %*% t(L_bnd)
  diag(R_bnd) <- 1

  for (m in c("ML", "ULS")) {
    expect_warning(EFA(R_bnd, 1, N = 200, method = m), class = "efa_heywood")
    efa_bnd <- suppressWarnings(EFA(R_bnd, 1, N = 200, method = m))
    expect_gt(length(efa_bnd$heywood), 0)
    expect_true(all(efa_bnd$h2[efa_bnd$heywood] < 1))
  }

  # The same matrix with a clearly proper structure is not flagged.
  L_ok <- c(.9, .7, .6, .5)
  R_ok <- L_ok %*% t(L_ok)
  diag(R_ok) <- 1
  expect_length(suppressWarnings(EFA(R_ok, 1, N = 200, method = "ML"))$heywood, 0)
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
   efa_paf_moderate, efa_ml_moderate, efa_uls_moderate, x, y, z, dat_sing, cor_sing,
   cor_nposdef)

test_that("print.EFA argument validators raise classed conditions", {
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


