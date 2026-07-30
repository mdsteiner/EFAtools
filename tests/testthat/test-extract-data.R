# Unit tests for .extract_data() (R/averaging.R), which stacks the per-run efa_fit()
# solutions of efa_average() into the loading/Phi/fit-index arrays the averaging step
# consumes. Covered across three fixtures: a clean run, one carrying a non-converged
# fit, and a rotated run.

efa_list <- list(EFA(test_models$baseline$cormat, n_factors = 3, N = 500),
                 EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                     method = "ML"),
                 EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                     method = "ULS"),
                 suppressWarnings(EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                     max_iter = 2)))


ext_a <- .extract_data(efa_list, test_models$baseline$cormat, 3, 4, "none", .3)

efa_list_er <- list(EFA(test_models$baseline$cormat, n_factors = 3, N = 500),
                 EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                     method = "ML"),
                 EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                     method = "ULS"),
                 try(suppressWarnings(EFA(test_models$baseline$cormat, n_factors = 15, N = 500,
                         type = "psych")),
                     silent = TRUE))

ext_er <- .extract_data(efa_list_er, test_models$baseline$cormat, 3, 4, "none", .3)

efa_list_rot <- list(EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                         rotation = "promax"),
                 EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                     method = "ML", rotation = "promax"),
                 EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                     method = "ULS", rotation = "promax"))

ext_rot <- .extract_data(efa_list_rot, test_models$baseline$cormat, 3, 3, "promax",
                         .3)

test_that(".extract_data works", {
  ### tests for ext_a with one non-convergence; no rotation; no error
  expect_type(ext_a, "list")
  expect_named(ext_a, c("L", "L_corres", "phi", "extract_phi", "h2",
                        "vars_accounted", "for_grid"))
  expect_named(ext_a$for_grid, c("errors", "error_m", "converged", "heywood",
                                 "admissible", "chisq", "p_chi", "caf", "cfi",
                                 "rmsea", "aic", "bic", "srmr", "tli", "ecvi",
                                 "rmsr"))
  checkmate::expect_array(ext_a$L)
  expect_equal(dim(ext_a$L), c(ncol(test_models$baseline$cormat), 3, 3))
  checkmate::expect_array(ext_a$L_corres)
  expect_equal(dim(ext_a$L_corres), c(ncol(test_models$baseline$cormat), 3, 3))
  expect_equal(ext_a$phi, NA)
  expect_equal(ext_a$extract_phi, FALSE)
  checkmate::expect_matrix(ext_a$h2)
  expect_equal(dim(ext_a$h2), c(4, ncol(test_models$baseline$cormat)))
  expect_equal(ext_a$for_grid$errors, rep(FALSE, 4))
  expect_true(all(is.na(ext_a$for_grid$error_m)))
  expect_equal(ext_a$for_grid$converged, c(0, 0, 0, 1))
  expect_equal(ext_a$for_grid$heywood, c(FALSE, FALSE, FALSE, NA))
  expect_equal(sum(is.na(ext_a$for_grid$chisq)), 2)
  expect_type(ext_a$for_grid$chisq, "double")
  expect_equal(sum(is.na(ext_a$for_grid$p_chi)), 2)
  expect_type(ext_a$for_grid$p_chi, "double")
  expect_equal(round(ext_a$for_grid$caf, 2), c(0.5, 0.5, 0.5, NA))
  expect_equal(ext_a$for_grid$cfi > .95, c(NA, TRUE, TRUE, NA))
  expect_equal(ext_a$for_grid$rmsea < .05, c(NA, TRUE, TRUE, NA))
  expect_equal(sign(ext_a$for_grid$aic), c(NA, -1, -1, NA))
  expect_equal(sign(ext_a$for_grid$bic), c(NA, -1, -1, NA))
  # SRMR/RMSR are residual-based (also computed for PAF); TLI/ECVI need a
  # chi-square and are NA for PAF. The non-converged run contributes NA to all.
  expect_type(ext_a$for_grid$srmr, "double")
  expect_equal(sum(is.na(ext_a$for_grid$srmr)), 1)
  expect_type(ext_a$for_grid$rmsr, "double")
  expect_equal(sum(is.na(ext_a$for_grid$rmsr)), 1)
  expect_type(ext_a$for_grid$tli, "double")
  expect_equal(sum(is.na(ext_a$for_grid$tli)), 2)
  expect_type(ext_a$for_grid$ecvi, "double")
  expect_equal(sum(is.na(ext_a$for_grid$ecvi)), 2)
  checkmate::expect_array(ext_a$vars_accounted)
  expect_equal(dim(ext_a$vars_accounted), c(3, 3, 3))

  ### tests for ext_er with one error; no rotation
  expect_type(ext_er, "list")
  expect_named(ext_er, c("L", "L_corres", "phi", "extract_phi", "h2",
                         "vars_accounted", "for_grid"))
  expect_named(ext_er$for_grid, c("errors", "error_m", "converged", "heywood",
                                  "admissible", "chisq", "p_chi", "caf", "cfi",
                                  "rmsea", "aic", "bic", "srmr", "tli", "ecvi",
                                  "rmsr"))
  checkmate::expect_array(ext_er$L)
  expect_equal(dim(ext_er$L), c(ncol(test_models$baseline$cormat), 3, 3))
  checkmate::expect_array(ext_er$L_corres)
  expect_equal(dim(ext_er$L_corres), c(ncol(test_models$baseline$cormat), 3, 3))
  expect_equal(ext_er$phi, NA)
  expect_equal(ext_er$extract_phi, FALSE)
  checkmate::expect_matrix(ext_er$h2)
  expect_equal(dim(ext_er$h2), c(4, ncol(test_models$baseline$cormat)))
  expect_equal(ext_er$for_grid$errors, c(rep(FALSE, 3), TRUE))
  expect_equal(sum(is.na(ext_er$for_grid$error_m)), 3)
  expect_equal(ext_er$for_grid$converged, c(0, 0, 0, NA))
  expect_equal(ext_er$for_grid$heywood, c(FALSE, FALSE, FALSE, NA))
  expect_equal(sum(is.na(ext_er$for_grid$chisq)), 2)
  expect_type(ext_er$for_grid$chisq, "double")
  expect_equal(sum(is.na(ext_er$for_grid$p_chi)), 2)
  expect_type(ext_er$for_grid$p_chi, "double")
  expect_equal(round(ext_er$for_grid$caf, 2), c(0.5, 0.5, 0.5, NA))
  expect_equal(ext_er$for_grid$cfi > .95, c(NA, TRUE, TRUE, NA))
  expect_equal(ext_er$for_grid$rmsea < .05, c(NA, TRUE, TRUE, NA))
  expect_equal(sign(ext_er$for_grid$aic), c(NA, -1, -1, NA))
  expect_equal(sign(ext_er$for_grid$bic), c(NA, -1, -1, NA))
  expect_type(ext_er$for_grid$srmr, "double")
  expect_equal(sum(is.na(ext_er$for_grid$srmr)), 1)
  expect_type(ext_er$for_grid$rmsr, "double")
  expect_equal(sum(is.na(ext_er$for_grid$rmsr)), 1)
  expect_type(ext_er$for_grid$tli, "double")
  expect_equal(sum(is.na(ext_er$for_grid$tli)), 2)
  expect_type(ext_er$for_grid$ecvi, "double")
  expect_equal(sum(is.na(ext_er$for_grid$ecvi)), 2)
  checkmate::expect_array(ext_er$vars_accounted)
  expect_equal(dim(ext_er$vars_accounted), c(3, 3, 3))


  ### tests for ext_rot with no errors; promax rotation
  expect_type(ext_rot, "list")
  expect_named(ext_rot, c("L", "L_corres", "phi", "extract_phi", "h2",
                          "vars_accounted", "for_grid"))
  expect_named(ext_rot$for_grid, c("errors", "error_m", "converged", "heywood",
                                   "admissible", "chisq", "p_chi", "caf", "cfi",
                                   "rmsea", "aic", "bic", "srmr", "tli", "ecvi",
                                   "rmsr"))
  checkmate::expect_array(ext_rot$L)
  expect_equal(dim(ext_rot$L), c(ncol(test_models$baseline$cormat), 3, 3))
  checkmate::expect_array(ext_rot$L_corres)
  expect_equal(dim(ext_rot$L_corres), c(ncol(test_models$baseline$cormat), 3, 3))
  checkmate::expect_array(ext_rot$phi)
  expect_equal(dim(ext_rot$phi), c(3, 3, 3))
  checkmate::expect_matrix(ext_rot$h2)
  expect_equal(dim(ext_rot$h2), c(3, ncol(test_models$baseline$cormat)))
  expect_equal(ext_rot$for_grid$errors, rep(FALSE, 3))
  expect_equal(ext_rot$for_grid$converged, c(0, 0, 0))
  expect_equal(ext_rot$for_grid$heywood, c(FALSE, FALSE, FALSE))
  expect_equal(sum(is.na(ext_rot$for_grid$chisq)), 1)
  expect_type(ext_rot$for_grid$chisq, "double")
  expect_equal(sum(is.na(ext_rot$for_grid$p_chi)), 1)
  expect_type(ext_rot$for_grid$p_chi, "double")
  expect_equal(round(ext_rot$for_grid$caf, 2), c(0.5, 0.5, 0.5))
  expect_equal(ext_rot$for_grid$cfi > .95, c(NA, TRUE, TRUE))
  expect_equal(ext_rot$for_grid$rmsea < .05, c(NA, TRUE, TRUE))
  expect_equal(sign(ext_rot$for_grid$aic), c(NA, -1, -1))
  expect_equal(sign(ext_rot$for_grid$bic), c(NA, -1, -1))
  expect_type(ext_rot$for_grid$srmr, "double")
  expect_equal(sum(is.na(ext_rot$for_grid$srmr)), 0)
  expect_type(ext_rot$for_grid$rmsr, "double")
  expect_equal(sum(is.na(ext_rot$for_grid$rmsr)), 0)
  expect_type(ext_rot$for_grid$tli, "double")
  expect_equal(sum(is.na(ext_rot$for_grid$tli)), 1)
  expect_type(ext_rot$for_grid$ecvi, "double")
  expect_equal(sum(is.na(ext_rot$for_grid$ecvi)), 1)
  checkmate::expect_array(ext_rot$vars_accounted)
  expect_equal(dim(ext_rot$vars_accounted), c(3, 3, 3))
})


rm(efa_list, ext_a, efa_list_er, ext_er, efa_list_rot, ext_rot)
