test_that(".numformat works", {
  expect_equal(.numformat(0.23), " .23")
  expect_equal(.numformat(0.2345, digits = 3), " .234")
  expect_equal(.numformat(0.2345, digits = 3, print_zero = TRUE), " 0.234")
})


efa_temp <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500)
efa_pro <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                rotation = "promax")
test_that(".compute_vars works", {
  expect_is(.compute_vars(efa_temp$unrot_loadings,
                             efa_temp$unrot_loadings), "matrix")
  expect_is(.compute_vars(efa_pro$rot_loadings,
                             efa_pro$unrot_loadings,
                             efa_pro$Phi), "matrix")
})

x_base <- population_models$loadings$baseline
x_NA <- population_models$loadings$baseline
x_NA[1, 3] <- NA
y_base <- x_base[, c(3,2,1)]
y_NA <- y_base
y_NA[2, 2] <- NA

test_that(".factor_congruence works", {
  expect_is(.factor_congruence(x_base, y_base), "matrix")
  expect_equal(sum(.factor_congruence(x_base, y_base)), 3)
  expect_warning(.factor_congruence(x_NA, y_NA, na.rm = FALSE), " Input contained missing values. Check your data or rerun with na.rm = TRUE.\n")
  expect_warning(.factor_congruence(x_NA, y_NA), " Input contained missing values. Analysis is performed on complete cases.\n")
})

efa_ml <- suppressWarnings(EFA(cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100),
                                     rnorm(100), rnorm(100)), 3, N = 500,
                               method = "ML"))
efa_uls <- suppressWarnings(EFA(cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100),
                    rnorm(100), rnorm(100)), 3, method = "ULS"))
efa_paf <- suppressWarnings(EFA(cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100),
                    rnorm(100), rnorm(100)), 3, method = "PAF"))

gof_ml <- .gof(efa_ml$unrot_loadings, efa_ml$orig_R, efa_ml$settings$N,
               "ML", efa_ml$fit_indices$Fm)
gof_uls <- .gof(efa_uls$unrot_loadings, efa_uls$orig_R, efa_uls$settings$N,
               "ULS", efa_uls$fit_indices$Fm)
gof_paf <- .gof(efa_paf$unrot_loadings, efa_paf$orig_R, efa_paf$settings$N,
               "PAF", NA)
m <- 6 # n variables
q <- 3 # n factors
test_that(".gof works", {
  expect_is(gof_ml, "list")
  expect_named(gof_ml,
               c("chi", "df", "p_chi", "CAF", "CFI", "RMSEA", "RMSEA_LB",
                 "RMSEA_UB", "AIC", "BIC", "Fm", "chi_null", "df_null",
                 "p_null"))
  expect_lt(gof_ml$p_chi, .05)
  expect_equal(gof_ml$CFI, 1)
  expect_equal(gof_ml$RMSEA, 0)
  expect_equal(gof_ml$CAF, .5, tolerance = .1)
  expect_equal(gof_ml$df, ((m - q)**2 - (m + q)) / 2)

  expect_is(gof_uls, "list")
  expect_named(gof_uls,
               c("chi", "df", "p_chi", "CAF", "CFI", "RMSEA", "RMSEA_LB",
                 "RMSEA_UB", "AIC", "BIC", "Fm", "chi_null", "df_null",
                 "p_null"))
  expect_lt(gof_uls$p_chi, .05)
  expect_equal(gof_uls$CFI, 1)
  expect_equal(gof_uls$RMSEA, 0)
  expect_equal(gof_uls$CAF, .5, tolerance = .1)
  expect_equal(gof_uls$df, ((m - q)**2 - (m + q)) / 2)

  expect_is(gof_paf, "list")
  expect_named(gof_paf,
               c("chi", "df", "p_chi", "CAF", "CFI", "RMSEA", "RMSEA_LB",
                 "RMSEA_UB", "AIC", "BIC", "Fm", "chi_null", "df_null",
                 "p_null"))
  expect_equal(gof_paf$chi, NA)
  expect_equal(gof_paf$p_chi, NA)
  expect_equal(gof_paf$CFI, NA)
  expect_equal(gof_paf$RMSEA, NA)
  expect_equal(gof_paf$CAF, .5, tolerance = .1)
  expect_equal(gof_paf$df, ((m - q)**2 - (m + q)) / 2)
  expect_equal(gof_paf$chi_null, NA)
  expect_equal(gof_paf$df_null, NA)
  expect_equal(gof_paf$p_null, NA)
})




test_that(".is_cormat works", {
  expect_equal(.is_cormat(cor(cbind(rnorm(100), rnorm(100)))), TRUE)
  expect_equal(.is_cormat(cbind(rnorm(100), rnorm(100))), FALSE)
  expect_equal(.is_cormat(cbind(rnorm(2), rnorm(2))), FALSE)
  expect_equal(.is_cormat(cbind(c(1, NA), rnorm(2))), FALSE)
  expect_equal(.is_cormat(matrix(c(1, .1, .3, 1), ncol = 2)), TRUE)
  expect_error(.is_cormat(matrix(c(1, NA, .3, 1), ncol = 2)), ' "x" is likely a correlation matrix but contains missing values. Please check the entered data.\n')
})

q_p <- .det_max_factors(8) + 1
test_that(".det_max_factors works", {
  expect_is(.det_max_factors(8), "numeric")
  expect_lte(((8 - q_p)**2 - (8 + q_p)) / 2, 0)
  expect_equal(.det_max_factors(0), 0)
  expect_equal(.det_max_factors(1), 0)
  expect_equal(.det_max_factors(2), 0)
  expect_equal(.det_max_factors(3), 0)
  expect_gt(.det_max_factors(4), 0)
})

dat_unname <- population_models$loadings$case_1a
dimnames(dat_unname) <- NULL

dat_unname_2 <- population_models$loadings$case_1a
colnames(dat_unname_2) <- NULL

test_that(".get_compare_matrix works", {
  expect_equal(capture_output(cat(.get_compare_matrix(population_models$loadings$case_1a))), "  \t F1  \t F2  \t F3  \nV1\t 0.600\t 0.000\t 0.000\nV2\t 0.600\t 0.000\t 0.000\nV3\t 0.000\t 0.600\t 0.000\nV4\t 0.000\t 0.600\t 0.000\nV5\t 0.000\t 0.000\t 0.600\nV6\t 0.000\t 0.000\t 0.600")
  expect_equal(capture_output(cat(.get_compare_matrix(dat_unname))), "  \t F1  \t F2  \t F3  \nV1\t 0.600\t 0.000\t 0.000\nV2\t 0.600\t 0.000\t 0.000\nV3\t 0.000\t 0.600\t 0.000\nV4\t 0.000\t 0.600\t 0.000\nV5\t 0.000\t 0.000\t 0.600\nV6\t 0.000\t 0.000\t 0.600")
  expect_equal(capture_output(cat(.get_compare_matrix(dat_unname_2))), "  \t F1  \t F2  \t F3  \nV1\t 0.600\t 0.000\t 0.000\nV2\t 0.600\t 0.000\t 0.000\nV3\t 0.000\t 0.600\t 0.000\nV4\t 0.000\t 0.600\t 0.000\nV5\t 0.000\t 0.000\t 0.600\nV6\t 0.000\t 0.000\t 0.600")
  expect_equal(capture_output(cat(.get_compare_matrix(population_models$loadings$case_1a, gof = TRUE))), " F1  \t F2  \t F3  \n 0.600\t 0.000\t 0.000\n 0.600\t 0.000\t 0.000\n 0.000\t 0.600\t 0.000\n 0.000\t 0.600\t 0.000\n 0.000\t 0.000\t 0.600\n 0.000\t 0.000\t 0.600")
  expect_equal(capture_output(cat(.get_compare_matrix(population_models$loadings$case_1a,
                                                      n_char = 1, gof = FALSE))), " \t F1  \t F2  \t F3  \nV\t 0.600\t 0.000\t 0.000\nV\t 0.600\t 0.000\t 0.000\nV\t 0.000\t 0.600\t 0.000\nV\t 0.000\t 0.600\t 0.000\nV\t 0.000\t 0.000\t 0.600\nV\t 0.000\t 0.000\t 0.600")
})

test_that(".get_compare_vector works", {
  expect_equal(capture_output(cat(.get_compare_vector(population_models$loadings$case_1a))), " 0.600   0.600   0.000   0.000   0.000   0.000   0.000\n 0.000   0.600   0.600   0.000   0.000   0.000   0.000\n 0.000   0.000   0.600   0.600")
})

test_that(".decimals works", {
  expect_is(.decimals(8), "numeric")
  expect_equal(.decimals(8), 0)
  expect_is(.decimals(8), "numeric")
  expect_error(.decimals("a"), " 'x' is of class 'character' but must be a numeric vector or matrix\n")
})

test_that(".settings_string works", {
  expect_equal(capture_output(cat(.settings_string(efa_ml$settings), sep = "")), "ML, none, EFAtools, 3, 100, pairwise.complete.obs, pearson, and psych")
  expect_equal(capture_output(cat(.settings_string(c("a", "b")), sep = "")), "a and b")
  expect_equal(capture_output(cat(.settings_string(c("a")), sep = "")), "a")
})

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
                 try(EFA(test_models$baseline$cormat, n_factors = 15, N = 500),
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
  expect_is(ext_a, "list")
  expect_named(ext_a, c("L", "L_corres", "phi", "extract_phi", "h2",
                        "vars_accounted", "for_grid"))
  expect_named(ext_a$for_grid, c("errors", "error_m", "converged", "heywood",
                                 "admissible", "chisq", "p_chi", "caf", "cfi",
                                 "rmsea", "aic", "bic"))
  expect_is(ext_a$L, "array")
  expect_equal(dim(ext_a$L), c(ncol(test_models$baseline$cormat), 3, 3))
  expect_is(ext_a$L_corres, "array")
  expect_equal(dim(ext_a$L_corres), c(ncol(test_models$baseline$cormat), 3, 3))
  expect_equal(ext_a$phi, NA)
  expect_equal(ext_a$extract_phi, FALSE)
  expect_is(ext_a$h2, "matrix")
  expect_equal(dim(ext_a$h2), c(4, ncol(test_models$baseline$cormat)))
  expect_equal(ext_a$for_grid$errors, rep(FALSE, 4))
  expect_true(all(is.na(ext_a$for_grid$error_m)))
  expect_equal(ext_a$for_grid$converged, c(0, 0, 0, 1))
  expect_equal(ext_a$for_grid$heywood, c(FALSE, FALSE, FALSE, NA))
  expect_equal(sum(is.na(ext_a$for_grid$chisq)), 2)
  expect_is(ext_a$for_grid$chisq, "numeric")
  expect_equal(sum(is.na(ext_a$for_grid$p_chi)), 2)
  expect_is(ext_a$for_grid$p_chi, "numeric")
  expect_equal(round(ext_a$for_grid$caf, 2), c(0.5, 0.5, 0.5, NA))
  expect_equal(ext_a$for_grid$cfi > .95, c(NA, TRUE, TRUE, NA))
  expect_equal(ext_a$for_grid$rmsea < .05, c(NA, TRUE, TRUE, NA))
  expect_equal(sign(ext_a$for_grid$aic), c(NA, -1, -1, NA))
  expect_equal(sign(ext_a$for_grid$bic), c(NA, -1, -1, NA))
  expect_is(ext_a$vars_accounted, "array")
  expect_equal(dim(ext_a$vars_accounted), c(3, 3, 3))

  ### tests for ext_er with one error; no rotation
  expect_is(ext_er, "list")
  expect_named(ext_er, c("L", "L_corres", "phi", "extract_phi", "h2",
                         "vars_accounted", "for_grid"))
  expect_named(ext_er$for_grid, c("errors", "error_m", "converged", "heywood",
                                  "admissible", "chisq", "p_chi", "caf", "cfi",
                                  "rmsea", "aic", "bic"))
  expect_is(ext_er$L, "array")
  expect_equal(dim(ext_er$L), c(ncol(test_models$baseline$cormat), 3, 3))
  expect_is(ext_er$L_corres, "array")
  expect_equal(dim(ext_er$L_corres), c(ncol(test_models$baseline$cormat), 3, 3))
  expect_equal(ext_er$phi, NA)
  expect_equal(ext_er$extract_phi, FALSE)
  expect_is(ext_er$h2, "matrix")
  expect_equal(dim(ext_er$h2), c(4, ncol(test_models$baseline$cormat)))
  expect_equal(ext_er$for_grid$errors, c(rep(FALSE, 3), TRUE))
  expect_equal(sum(is.na(ext_er$for_grid$error_m)), 3)
  expect_equal(ext_er$for_grid$converged, c(0, 0, 0, NA))
  expect_equal(ext_er$for_grid$heywood, c(FALSE, FALSE, FALSE, NA))
  expect_equal(sum(is.na(ext_er$for_grid$chisq)), 2)
  expect_is(ext_er$for_grid$chisq, "numeric")
  expect_equal(sum(is.na(ext_er$for_grid$p_chi)), 2)
  expect_is(ext_er$for_grid$p_chi, "numeric")
  expect_equal(round(ext_er$for_grid$caf, 2), c(0.5, 0.5, 0.5, NA))
  expect_equal(ext_er$for_grid$cfi > .95, c(NA, TRUE, TRUE, NA))
  expect_equal(ext_er$for_grid$rmsea < .05, c(NA, TRUE, TRUE, NA))
  expect_equal(sign(ext_er$for_grid$aic), c(NA, -1, -1, NA))
  expect_equal(sign(ext_er$for_grid$bic), c(NA, -1, -1, NA))
  expect_is(ext_er$vars_accounted, "array")
  expect_equal(dim(ext_er$vars_accounted), c(3, 3, 3))


  ### tests for ext_rot with no errors; promax rotation
  expect_is(ext_rot, "list")
  expect_named(ext_rot, c("L", "L_corres", "phi", "extract_phi", "h2",
                          "vars_accounted", "for_grid"))
  expect_named(ext_rot$for_grid, c("errors", "error_m", "converged", "heywood",
                                   "admissible", "chisq", "p_chi", "caf", "cfi",
                                   "rmsea", "aic", "bic"))
  expect_is(ext_rot$L, "array")
  expect_equal(dim(ext_rot$L), c(ncol(test_models$baseline$cormat), 3, 3))
  expect_is(ext_rot$L_corres, "array")
  expect_equal(dim(ext_rot$L_corres), c(ncol(test_models$baseline$cormat), 3, 3))
  expect_is(ext_rot$phi, "array")
  expect_equal(dim(ext_rot$phi), c(3, 3, 3))
  expect_is(ext_rot$h2, "matrix")
  expect_equal(dim(ext_rot$h2), c(3, ncol(test_models$baseline$cormat)))
  expect_equal(ext_rot$for_grid$errors, rep(FALSE, 3))
  expect_equal(ext_rot$for_grid$converged, c(0, 0, 0))
  expect_equal(ext_rot$for_grid$heywood, c(FALSE, FALSE, FALSE))
  expect_equal(sum(is.na(ext_rot$for_grid$chisq)), 1)
  expect_is(ext_rot$for_grid$chisq, "numeric")
  expect_equal(sum(is.na(ext_rot$for_grid$p_chi)), 1)
  expect_is(ext_rot$for_grid$p_chi, "numeric")
  expect_equal(round(ext_rot$for_grid$caf, 2), c(0.5, 0.5, 0.5))
  expect_equal(ext_rot$for_grid$cfi > .95, c(NA, TRUE, TRUE))
  expect_equal(ext_rot$for_grid$rmsea < .05, c(NA, TRUE, TRUE))
  expect_equal(sign(ext_rot$for_grid$aic), c(NA, -1, -1))
  expect_equal(sign(ext_rot$for_grid$bic), c(NA, -1, -1))
  expect_is(ext_rot$vars_accounted, "array")
  expect_equal(dim(ext_rot$vars_accounted), c(3, 3, 3))
})


av_mean_NA <- .average_values(L = array(c(rep(1, 9), rep(3, 9), rep(4, 9)),
                                          c(3, 3, 3)),
                                L_corres = array(c(1, 0, 0, 0, 1, 0, 0, 0, 1,
                                                   1, 0, 0, 0, 0, 1, 0, 1, 0,
                                                   1, 0, 0, 0, 1, 0, 0, 0, 1),
                                                 c(3, 3, 3)),
                                vars_accounted = array(c(rep(1, 9), rep(3, 9), rep(4, 9)),
                                                       c(3, 3, 3)),
                                h2 = matrix(rep(c(1, 3, 4), each = 3), ncol = 3, byrow = TRUE),
                                phi = NA,
                                extract_phi = FALSE,
                                averaging = "mean",
                                trim = 0,
                                for_grid = data.frame(chisq = c(1, 3, 4),
                                                      p_chi = c(1, 3, 4),
                                                      caf = c(1, 3, 4),
                                                      cfi = c(1, 3, 4),
                                                      rmsea = c(1, 3, 4),
                                                      aic = c(1, 3, 4),
                                                      bic= c(1, 3, 4)),
                                df = 5, ind_names = paste0("Ind", 1:3))

av_mean_NA_t01 <- .average_values(L = array(c(rep(1, 9), rep(3, 9), rep(4, 9)),
                                           c(3, 3, 3)),
                                 L_corres = array(c(1, 0, 0, 0, 1, 0, 0, 0, 1,
                                                    1, 0, 0, 0, 0, 1, 0, 1, 0,
                                                    1, 0, 0, 0, 1, 0, 0, 0, 1),
                                                  c(3, 3, 3)),
                                 vars_accounted = array(c(rep(1, 9), rep(3, 9), rep(4, 9)),
                                                        c(3, 3, 3)),
                                 h2 = matrix(rep(c(1, 3, 4), each = 3), ncol = 3, byrow = TRUE),
                                 phi = NA,
                                 extract_phi = FALSE,
                                 averaging = "mean",
                                 trim = .5,
                                 for_grid = data.frame(chisq = c(1, 3, 4),
                                                       p_chi = c(1, 3, 4),
                                                       caf = c(1, 3, 4),
                                                       cfi = c(1, 3, 4),
                                                       rmsea = c(1, 3, 4),
                                                       aic = c(1, 3, 4),
                                                       bic= c(1, 3, 4)),
                                 df = 5, ind_names = paste0("Ind", 1:3))

av_mean <- .average_values(L = array(c(rep(1, 9), rep(3, 9), rep(4, 9)),
                                           c(3, 3, 3)),
                                 L_corres = array(c(1, 0, 0, 0, 1, 0, 0, 0, 1,
                                                    1, 0, 0, 0, 0, 1, 0, 1, 0,
                                                    1, 0, 0, 0, 1, 0, 0, 0, 1),
                                                  c(3, 3, 3)),
                              vars_accounted = array(c(rep(1, 9), rep(3, 9), rep(4, 9)),
                                                     c(3, 3, 3)),
                                 h2 = matrix(rep(c(1, 3, 4), each = 3), ncol = 3, byrow = TRUE),
                                 phi = array(c(rep(1, 9), rep(3, 9), rep(4, 9)),
                                             c(3, 3, 3)),
                                 extract_phi = TRUE,
                                 averaging = "mean",
                                 trim = 0,
                                 for_grid = data.frame(chisq = c(1, 3, 4),
                                                       p_chi = c(1, 3, 4),
                                                       caf = c(1, 3, 4),
                                                       cfi = c(1, 3, 4),
                                                       rmsea = c(1, 3, 4),
                                                       aic = c(1, 3, 4),
                                                       bic= c(1, 3, 4)),
                                 df = 5, ind_names = paste0("Ind", 1:3))
av_median <- .average_values(L = array(c(rep(1, 9), rep(3, 9), rep(4, 9)),
                                        c(3, 3, 3)),
                              L_corres = array(c(1, 0, 0, 0, 1, 0, 0, 0, 1,
                                                 1, 0, 0, 0, 0, 1, 0, 1, 0,
                                                 1, 0, 0, 0, 1, 0, 0, 0, 1),
                                               c(3, 3, 3)),
                              vars_accounted = array(c(rep(1, 9), rep(3, 9), rep(4, 9)),
                                                     c(3, 3, 3)),
                              h2 = matrix(rep(c(1, 3, 4), each = 3), ncol = 3, byrow = TRUE),
                              phi = array(c(rep(1, 9), rep(3, 9), rep(4, 9)),
                                          c(3, 3, 3)),
                              extract_phi = TRUE,
                              averaging = "median",
                              trim = 0,
                              for_grid = data.frame(chisq = c(1, 3, 4),
                                                    p_chi = c(1, 3, 4),
                                                    caf = c(1, 3, 4),
                                                    cfi = c(1, 3, 4),
                                                    rmsea = c(1, 3, 4),
                                                    aic = c(1, 3, 4),
                                                    bic= c(1, 3, 4)),
                              df = 5, ind_names = paste0("Ind", 1:3))

av_median_NA <- .average_values(L = array(c(rep(1, 9), rep(3, 9), rep(4, 9)),
                                           c(3, 3, 3)),
                                 L_corres = array(c(1, 0, 0, 0, 1, 0, 0, 0, 1,
                                                    1, 0, 0, 0, 0, 1, 0, 1, 0,
                                                    1, 0, 0, 0, 1, 0, 0, 0, 1),
                                                  c(3, 3, 3)),
                                 vars_accounted = array(c(rep(1, 9), rep(3, 9), rep(4, 9)),
                                                        c(3, 3, 3)),
                                 h2 = matrix(rep(c(1, 3, 4), each = 3), ncol = 3, byrow = TRUE),
                                 phi = NA,
                                 extract_phi = FALSE,
                                 averaging = "median",
                                 trim = 0.1,
                                 for_grid = data.frame(chisq = c(1, 3, 4),
                                                       p_chi = c(1, 3, 4),
                                                       caf = c(1, 3, 4),
                                                       cfi = c(1, 3, 4),
                                                       rmsea = c(1, 3, 4),
                                                       aic = c(1, 3, 4),
                                                       bic= c(1, 3, 4)),
                                 df = 5, ind_names = paste0("Ind", 1:3))

test_that(".average_values works", {
  ### tests for av_mean_NA with extract_phi = FALSE and trim = 0
  expect_is(av_mean_NA, "list")
  expect_named(av_mean_NA, c("h2", "loadings", "phi", "vars_accounted",
                              "ind_fac_corres", "fit_indices"))
  expect_is(av_mean_NA$h2, "list")
  expect_named(av_mean_NA$h2, c("average", "sd", "min", "max", "range"))
  expect_is(av_mean_NA$h2$average, "numeric")
  expect_equal(unname(round(av_mean_NA$h2$average, 2)), rep(2.67, 3))
  expect_named(av_mean_NA$h2$average, paste0("Ind", 1:3))
  expect_equal(unname(av_mean_NA$h2$sd), rep(1.527525, 3), tolerance = .01)
  expect_equal(unname(av_mean_NA$h2$min), rep(1, 3))
  expect_equal(unname(av_mean_NA$h2$max), rep(4, 3))
  expect_is(av_mean_NA$loadings, "list")
  expect_named(av_mean_NA$loadings, c("average", "sd", "min", "max", "range"))
  expect_is(av_mean_NA$loadings$average, "LOADINGS")
  expect_equal(unclass(round(av_mean_NA$loadings$average, 2)), matrix(rep(2.67, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(av_mean_NA$loadings$sd, matrix(rep(1.527525, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))), tolerance = .01)
  expect_equal(unclass(av_mean_NA$loadings$min), matrix(rep(1, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(unclass(av_mean_NA$loadings$max), matrix(rep(4, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(av_mean_NA$loadings$range, matrix(rep(3, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))), tolerance = .01)
  expect_is(av_mean_NA$vars_accounted, "list")
  expect_named(av_mean_NA$vars_accounted, c("average", "sd", "min", "max", "range"))
  expect_is(av_mean_NA$vars_accounted$average, "matrix")
  expect_equal(round(av_mean_NA$vars_accounted$average, 2), matrix(rep(2.67, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_mean_NA$vars_accounted$sd, matrix(rep(1.527525, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))), tolerance = .01)
  expect_equal(av_mean_NA$vars_accounted$min, matrix(rep(1, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_mean_NA$vars_accounted$max, matrix(rep(4, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_mean_NA$vars_accounted$range, matrix(rep(3, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))), tolerance = .01)
  expect_equal(av_mean_NA$phi, NA)
  expect_is(av_mean_NA$ind_fac_corres, "matrix")
  expect_equal(round(av_mean_NA$ind_fac_corres, 2),
               matrix(c(1, 0, 0, 0, .67, .33, 0, .33, .67), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_is(av_mean_NA$fit_indices, "data.frame")
  expect_named(av_mean_NA$fit_indices, c("index", "average", "sd", "range",
                                          "min", "max"))
  expect_is(av_mean_NA$fit_indices$index, "character")
  expect_equal(av_mean_NA$fit_indices$index, c("chisq", "p_chi", "caf", "cfi",
                                             "rmsea", "aic", "bic", "df"))
  expect_is(av_mean_NA$fit_indices$average, "numeric")
  expect_equal(round(av_mean_NA$fit_indices$average, 2), c(rep(2.67, 7), 5))
  expect_is(av_mean_NA$fit_indices$sd, "numeric")
  expect_equal(round(av_mean_NA$fit_indices$sd, 2), c(rep(1.53, 7), 5))
  expect_is(av_mean_NA$fit_indices$range, "numeric")
  expect_equal(av_mean_NA$fit_indices$range, c(rep(3, 7), 5))
  expect_is(av_mean_NA$fit_indices$min, "numeric")
  expect_equal(av_mean_NA$fit_indices$min, c(rep(1, 7), 5))
  expect_is(av_mean_NA$fit_indices$max, "numeric")
  expect_equal(av_mean_NA$fit_indices$max, c(rep(4, 7), 5))


  ### tests for av_mean_NA_t01 with extract_phi = FALSE and trim = .10
  expect_is(av_mean_NA_t01, "list")
  expect_named(av_mean_NA_t01, c("h2", "loadings", "phi", "vars_accounted",
                                  "ind_fac_corres", "fit_indices"))
  expect_is(av_mean_NA_t01$h2, "list")
  expect_named(av_mean_NA_t01$h2, c("average", "sd", "min", "max", "range"))
  expect_is(av_mean_NA_t01$h2$average, "numeric")
  expect_equal(unname(round(av_mean_NA_t01$h2$average, 2)), rep(3, 3))
  expect_named(av_mean_NA_t01$h2$average, paste0("Ind", 1:3))
  expect_equal(unname(av_mean_NA_t01$h2$sd), rep(1.527525, 3), tolerance = .01)
  expect_equal(unname(av_mean_NA_t01$h2$min), rep(1, 3))
  expect_equal(unname(av_mean_NA_t01$h2$max), rep(4, 3))
  expect_is(av_mean_NA_t01$loadings, "list")
  expect_named(av_mean_NA_t01$loadings, c("average", "sd", "min", "max", "range"))
  expect_is(av_mean_NA_t01$loadings$average, "LOADINGS")
  expect_equal(unclass(round(av_mean_NA_t01$loadings$average, 2)), matrix(rep(3, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(av_mean_NA_t01$loadings$sd, matrix(rep(1.527525, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))), tolerance = .01)
  expect_equal(unclass(av_mean_NA_t01$loadings$min), matrix(rep(1, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(unclass(av_mean_NA_t01$loadings$max), matrix(rep(4, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_is(av_mean_NA_t01$vars_accounted, "list")
  expect_named(av_mean_NA_t01$vars_accounted, c("average", "sd", "min", "max", "range"))
  expect_is(av_mean_NA_t01$vars_accounted$average, "matrix")
  expect_equal(round(av_mean_NA_t01$vars_accounted$average, 2), matrix(rep(3, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_mean_NA_t01$vars_accounted$sd, matrix(rep(1.527525, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))), tolerance = .01)
  expect_equal(av_mean_NA_t01$vars_accounted$min, matrix(rep(1, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_mean_NA_t01$vars_accounted$max, matrix(rep(4, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_mean_NA_t01$vars_accounted$range, matrix(rep(3, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))), tolerance = .01)
  expect_equal(av_mean_NA_t01$phi, NA)
  expect_is(av_mean_NA_t01$ind_fac_corres, "matrix")
  expect_equal(round(av_mean_NA_t01$ind_fac_corres, 2),
               matrix(c(1, 0, 0, 0, .67, .33, 0, .33, .67), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_is(av_mean_NA_t01$fit_indices, "data.frame")
  expect_named(av_mean_NA_t01$fit_indices, c("index", "average", "sd", "range",
                                          "min", "max"))
  expect_is(av_mean_NA_t01$fit_indices$index, "character")
  expect_equal(av_mean_NA_t01$fit_indices$index, c("chisq", "p_chi", "caf", "cfi",
                                                "rmsea", "aic", "bic", "df"))
  expect_is(av_mean_NA_t01$fit_indices$average, "numeric")
  expect_equal(round(av_mean_NA_t01$fit_indices$average, 2), c(rep(3, 7), 5))
  expect_is(av_mean_NA_t01$fit_indices$sd, "numeric")
  expect_equal(round(av_mean_NA_t01$fit_indices$sd, 2), c(rep(1.53, 7), 5))
  expect_is(av_mean_NA_t01$fit_indices$range, "numeric")
  expect_equal(av_mean_NA_t01$fit_indices$range, c(rep(3, 7), 5))
  expect_is(av_mean_NA_t01$fit_indices$min, "numeric")
  expect_equal(av_mean_NA_t01$fit_indices$min, c(rep(1, 7), 5))
  expect_is(av_mean_NA_t01$fit_indices$max, "numeric")
  expect_equal(av_mean_NA_t01$fit_indices$max, c(rep(4, 7), 5))


  ### tests for av_mean with extract_phi = TRUE (only affected output tested)
  expect_is(av_mean$phi, "list")
  expect_is(av_mean$phi$average, "matrix")
  expect_equal(round(av_mean$phi$average, 2), matrix(rep(2.67, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))
  expect_is(av_mean$phi$sd, "matrix")
  expect_equal(round(av_mean$phi$sd, 2), matrix(rep(1.53, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))
  expect_is(av_mean$phi$min, "matrix")
  expect_equal(av_mean$phi$min, matrix(rep(1, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))
  expect_is(av_mean$phi$max, "matrix")
  expect_equal(av_mean$phi$max, matrix(rep(4, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))


  ### tests for av_median_NA with extract_phi = FALSE
  expect_is(av_median_NA, "list")
  expect_named(av_median_NA, c("h2", "loadings", "phi", "vars_accounted",
                                "ind_fac_corres", "fit_indices"))
  expect_is(av_median_NA$h2, "list")
  expect_named(av_median_NA$h2, c("average", "sd", "min", "max", "range"))
  expect_is(av_median_NA$h2$average, "numeric")
  expect_equal(unname(round(av_median_NA$h2$average, 2)), rep(3, 3))
  expect_named(av_median_NA$h2$average, paste0("Ind", 1:3))
  expect_equal(unname(av_median_NA$h2$sd), rep(1.527525, 3), tolerance = .01)
  expect_equal(unname(av_median_NA$h2$min), rep(1, 3))
  expect_equal(unname(av_median_NA$h2$max), rep(4, 3))
  expect_is(av_median_NA$loadings, "list")
  expect_named(av_median_NA$loadings, c("average", "sd", "min", "max", "range"))
  expect_is(av_median_NA$loadings$average, "LOADINGS")
  expect_equal(unclass(round(av_median_NA$loadings$average, 2)), matrix(rep(3, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(av_median_NA$loadings$sd, matrix(rep(1.527525, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))), tolerance = .01)
  expect_equal(unclass(av_median_NA$loadings$min), matrix(rep(1, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(unclass(av_median_NA$loadings$max), matrix(rep(4, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_is(av_median_NA$vars_accounted, "list")
  expect_named(av_median_NA$vars_accounted, c("average", "sd", "min", "max", "range"))
  expect_is(av_median_NA$vars_accounted$average, "matrix")
  expect_equal(round(av_median_NA$vars_accounted$average, 2), matrix(rep(3, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_median_NA$vars_accounted$sd, matrix(rep(1.527525, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))), tolerance = .01)
  expect_equal(av_median_NA$vars_accounted$min, matrix(rep(1, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_median_NA$vars_accounted$max, matrix(rep(4, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_median_NA$vars_accounted$range, matrix(rep(3, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))), tolerance = .01)
  expect_equal(av_median_NA$phi, NA)
  expect_is(av_median_NA$ind_fac_corres, "matrix")
  expect_equal(round(av_median_NA$ind_fac_corres, 2),
               matrix(c(1, 0, 0, 0, .67, .33, 0, .33, .67), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_is(av_median_NA$fit_indices, "data.frame")
  expect_named(av_median_NA$fit_indices, c("index", "average", "sd", "range",
                                              "min", "max"))
  expect_is(av_median_NA$fit_indices$index, "character")
  expect_equal(av_median_NA$fit_indices$index, c("chisq", "p_chi", "caf", "cfi",
                                                    "rmsea", "aic", "bic", "df"))
  expect_is(av_median_NA$fit_indices$average, "numeric")
  expect_equal(round(av_median_NA$fit_indices$average, 2), c(rep(3, 7), 5))
  expect_is(av_median_NA$fit_indices$sd, "numeric")
  expect_equal(round(av_median_NA$fit_indices$sd, 2), c(rep(1.53, 7), 5))
  expect_is(av_median_NA$fit_indices$range, "numeric")
  expect_equal(av_median_NA$fit_indices$range, c(rep(3, 7), 5))
  expect_is(av_median_NA$fit_indices$min, "numeric")
  expect_equal(av_median_NA$fit_indices$min, c(rep(1, 7), 5))
  expect_is(av_median_NA$fit_indices$max, "numeric")
  expect_equal(av_median_NA$fit_indices$max, c(rep(4, 7), 5))


  ### tests for av_median with extract_phi = TRUE (only affected output tested)
  expect_is(av_median$phi, "list")
  expect_is(av_median$phi$average, "matrix")
  expect_equal(round(av_median$phi$average, 2), matrix(rep(3, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))
  expect_is(av_median$phi$sd, "matrix")
  expect_equal(round(av_median$phi$sd, 2), matrix(rep(1.53, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))
  expect_is(av_median$phi$min, "matrix")
  expect_equal(av_median$phi$min, matrix(rep(1, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))
  expect_is(av_median$phi$max, "matrix")
  expect_equal(av_median$phi$max, matrix(rep(4, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))


})

arr_re_NA <- .array_reorder(L = array(c(rep(.6, 6), rep(0, 12),
                                     rep(0, 6), rep(.6, 6), rep(0, 6),
                                     rep(0, 6), rep(0, 6), rep(.6, 6),
                                     rep(.6, 6), rep(0, 12),
                                     rep(0, 6), rep(.6, 6), rep(0, 6),
                                     rep(0, 6), rep(0, 6), rep(.6, 6),
                                     rep(0, 6), rep(.6, 6), rep(0, 6),
                                     rep(-.6, 6), rep(0, 12),
                                     rep(0, 6), rep(0, 6), rep(.6, 6)),
                                   c(18, 3, 3)),
                            vars_accounted = array(c(rep(.2, 3),
                                                     rep(.3, 3),
                                                     rep(.4, 3),
                                                     rep(.2, 3),
                                                     rep(.3, 3),
                                                     rep(.4, 3),
                                                     rep(.3, 3),
                                                     rep(.2, 3),
                                                     rep(.4, 3)),
                                                   c(3, 3, 3)),
                         L_corres = array(as.numeric(abs(c(rep(.6, 6), rep(0, 12),
                                            rep(0, 6), rep(.6, 6), rep(0, 6),
                                            rep(0, 6), rep(0, 6), rep(.6, 6),
                                            rep(.6, 6), rep(0, 12),
                                            rep(0, 6), rep(.6, 6), rep(0, 6),
                                            rep(0, 6), rep(0, 6), rep(.6, 6),
                                            rep(0, 6), rep(.6, 6), rep(0, 6),
                                            rep(.6, 6), rep(0, 12),
                                            rep(0, 6), rep(0, 6), rep(.6, 6)))> .3),
                                          c(18, 3, 3)),
                         phi = NA, extract_phi = FALSE, n_factors = 3)

arr_re <- .array_reorder(L = array(c(rep(.6, 6), rep(0, 12),
                                        rep(0, 6), rep(.6, 6), rep(0, 6),
                                        rep(0, 6), rep(0, 6), rep(.6, 6),
                                        rep(.6, 6), rep(0, 12),
                                        rep(0, 6), rep(.6, 6), rep(0, 6),
                                        rep(0, 6), rep(0, 6), rep(.6, 6),
                                        rep(0, 6), rep(.6, 6), rep(0, 6),
                                        rep(-.6, 6), rep(0, 12),
                                        rep(0, 6), rep(0, 6), rep(.6, 6)),
                                      c(18, 3, 3)),
                         vars_accounted = array(c(rep(.2, 3),
                                                  rep(.3, 3),
                                                  rep(.4, 3),
                                                  rep(.2, 3),
                                                  rep(.3, 3),
                                                  rep(.4, 3),
                                                  rep(.3, 3),
                                                  rep(.2, 3),
                                                  rep(.4, 3)),
                                                c(3, 3, 3)),
                            L_corres = array(as.numeric(abs(c(rep(.6, 6), rep(0, 12),
                                                              rep(0, 6), rep(.6, 6), rep(0, 6),
                                                              rep(0, 6), rep(0, 6), rep(.6, 6),
                                                              rep(.6, 6), rep(0, 12),
                                                              rep(0, 6), rep(.6, 6), rep(0, 6),
                                                              rep(0, 6), rep(0, 6), rep(.6, 6),
                                                              rep(0, 6), rep(.6, 6), rep(0, 6),
                                                              rep(.6, 6), rep(0, 12),
                                                              rep(0, 6), rep(0, 6), rep(.6, 6)))> .3),
                                             c(18, 3, 3)),
                            phi = array(rep(c(1, .3, .4, .3, 1, .2, .4, .2, 1), 3), c(3, 3, 3)),
                         extract_phi = TRUE, n_factors = 3)
test_that(".array_reorder works", {
  ### tests for arr_re_NA with phi = NA and extract_phi = FALSE
  expect_is(arr_re_NA, "list")
  expect_named(arr_re_NA, c("L", "L_corres", "phi", "vars_accounted"))
  expect_is(arr_re_NA$L, "array")
  expect_equal(dim(arr_re_NA$L), c(18, 3, 3))
  expect_equal(arr_re_NA$L,
               array(c(rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6),
                       rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6),
                       rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6)),
                      c(18, 3, 3)))
  expect_is(arr_re_NA$L_corres, "array")
  expect_equal(dim(arr_re_NA$L_corres), c(18, 3, 3))
  expect_equal(arr_re_NA$L_corres,
               array(as.numeric(c(rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6),
                       rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6),
                       rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6)) > .3),
                     c(18, 3, 3)))
  expect_equal(arr_re_NA$phi, NA)
  expect_is(arr_re_NA$vars_accounted, "array")
  expect_equal(arr_re_NA$vars_accounted,
               array(c(rep(.2, 3),
                       rep(.3, 3),
                       rep(.4, 3),
                       rep(.2, 3),
                       rep(.3, 3),
                       rep(.4, 3),
                       rep(.2, 3),
                       rep(.3, 3),
                       rep(.4, 3)),
                     c(3, 3, 3)))

  ### tests for arr_re with phi = array() and extract_phi = TRUE
  expect_is(arr_re, "list")
  expect_named(arr_re, c("L", "L_corres", "phi", "vars_accounted"))
  expect_is(arr_re$L, "array")
  expect_equal(dim(arr_re$L), c(18, 3, 3))
  expect_equal(arr_re$L,
               array(c(rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6),
                       rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6),
                       rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6)),
                     c(18, 3, 3)))
  expect_is(arr_re$L_corres, "array")
  expect_equal(dim(arr_re$L_corres), c(18, 3, 3))
  expect_equal(arr_re$L_corres,
               array(as.numeric(c(rep(.6, 6), rep(0, 12),
                                  rep(0, 6), rep(.6, 6), rep(0, 6),
                                  rep(0, 6), rep(0, 6), rep(.6, 6),
                                  rep(.6, 6), rep(0, 12),
                                  rep(0, 6), rep(.6, 6), rep(0, 6),
                                  rep(0, 6), rep(0, 6), rep(.6, 6),
                                  rep(.6, 6), rep(0, 12),
                                  rep(0, 6), rep(.6, 6), rep(0, 6),
                                  rep(0, 6), rep(0, 6), rep(.6, 6)) > .3),
                     c(18, 3, 3)))
  expect_is(arr_re$phi, "array")
  expect_equal(dim(arr_re$phi), c(3, 3, 3))
  expect_equal(arr_re$phi,
               array(c(1, .3, .4, .3, 1, .2, .4, .2, 1,
                       1, .3, .4, .3, 1, .2, .4, .2, 1,
                       1, .3, .2, .3, 1, .4, .2, .4, 1), c(3, 3, 3)))
  expect_is(arr_re$vars_accounted, "array")
  expect_equal(arr_re$vars_accounted,
               array(c(rep(.2, 3),
                       rep(.3, 3),
                       rep(.4, 3),
                       rep(.2, 3),
                       rep(.3, 3),
                       rep(.4, 3),
                       rep(.2, 3),
                       rep(.3, 3),
                       rep(.4, 3)),
                     c(3, 3, 3)))

})

### test .oblq_grid

obl_grid_1 <- .oblq_grid(c("PAF"), c("smc", "mac"), .001,
                         c("sum", "max_individual"), c(FALSE, TRUE),
                         NA, c("promax", "simplimax", "oblimin"),
                         c(3, 4), TRUE, c("norm", "unnorm"), 1e-5, c("kaiser", "svd"),
                         30)
obl_grid_2 <- .oblq_grid(c("PAF"), c("smc", "mac"), .001,
                         c("sum", "max_individual"), c(FALSE, TRUE),
                         NA, c("simplimax", "oblimin"),
                         c(3, 4), TRUE, c("norm", "unnorm"), 1e-5, c("kaiser", "svd"),
                         30)
obl_grid_3 <- .oblq_grid(c("PAF"), c("smc", "mac"), .001,
                         c("sum", "max_individual"), c(FALSE, TRUE),
                         NA, c("simplimax", "oblimin"),
                         NA, TRUE, NA, 1e-5, NA, 30)
obl_grid_4 <- .oblq_grid("ML", NA, NA, NA, NA, c("psych", "factanal"),
                         "oblimin", NA, TRUE, NA, 1e-5, NA, NA)

test_that(".oblq_grid works", {
  ### tests for arr_re_NA with phi = NA and extract_phi = FALSE
  expect_is(obl_grid_1, "data.frame")
  expect_named(obl_grid_1, c("method", "init_comm", "criterion", "criterion_type",
                            "abs_eigen", "start_method", "rotation", "k_promax",
                            "normalize", "P_type", "precision", "varimax_type",
                            "k_simplimax"))
  expect_equal(nrow(obl_grid_1), 80)
  expect_equal(sum(is.na(obl_grid_1$k_simplimax)), 72)
  expect_equal(sum(is.na(obl_grid_1$k_promax)), 16)
  expect_equal(sum(is.na(obl_grid_1$varimax_type)), 16)
  expect_equal(unique(obl_grid_1$rotation), c("promax", "simplimax", "oblimin"))

  expect_is(obl_grid_2, "data.frame")
  expect_named(obl_grid_2, c("method", "init_comm", "criterion", "criterion_type",
                             "abs_eigen", "start_method", "rotation", "k_promax",
                             "normalize", "P_type", "precision", "varimax_type",
                             "k_simplimax"))
  expect_equal(nrow(obl_grid_2), 16)
  expect_equal(sum(is.na(obl_grid_2$k_promax)), 16)
  expect_equal(sum(is.na(obl_grid_2$varimax_type)), 16)
  expect_equal(unique(obl_grid_2$rotation), c("simplimax", "oblimin"))

  expect_equal(obl_grid_2, obl_grid_3)

  expect_is(obl_grid_4, "data.frame")
  expect_named(obl_grid_4, c("method", "init_comm", "criterion", "criterion_type",
                             "abs_eigen", "start_method", "rotation", "k_promax",
                             "normalize", "P_type", "precision", "varimax_type",
                             "k_simplimax"))
  expect_equal(nrow(obl_grid_4), 2)
  expect_equal(sum(is.na(obl_grid_4$k_simplimax)), 2)
  expect_equal(sum(is.na(obl_grid_4$k_promax)), 2)
  expect_equal(sum(is.na(obl_grid_4$varimax_type)), 2)
  expect_equal(unique(obl_grid_4$rotation), c("oblimin"))
  expect_equal(sum(is.na(obl_grid_4$init_comm)), 2)


})
### test .orth_grid

orth_grid_1 <- .orth_grid(c("PAF"), c("smc", "mac"), .001,
                         c("sum", "max_individual"), c(FALSE, TRUE),
                         NA, c("varimax", "quartimax"),
                         TRUE, 1e-5, c("kaiser", "svd"))
orth_grid_2 <- .orth_grid(c("PAF"), c("smc", "mac"), .001,
                          c("sum", "max_individual"), c(FALSE, TRUE),
                          NA, c("quartimax"), TRUE, 1e-5,
                          c("kaiser", "svd"))
orth_grid_3 <- .orth_grid(c("PAF"), c("smc", "mac"), .001,
                          c("sum", "max_individual"), c(FALSE, TRUE),
                          NA, c("quartimax"), TRUE, 1e-5, NA)
orth_grid_4 <- .orth_grid("ML", NA, NA, NA, NA, c("psych", "factanal"),
                         "quartimax", TRUE, 1e-5, NA)

test_that(".orth_grid works", {
  ### tests for arr_re_NA with phi = NA and extract_phi = FALSE
  expect_is(orth_grid_1, "data.frame")
  expect_named(orth_grid_1, c("method", "init_comm", "criterion", "criterion_type",
                             "abs_eigen", "start_method", "rotation", "k_promax",
                             "normalize", "P_type", "precision", "varimax_type",
                             "k_simplimax"))
  expect_equal(nrow(orth_grid_1), 24)
  expect_equal(sum(is.na(orth_grid_1$varimax_type)), 8)
  expect_equal(sum(is.na(orth_grid_1$k_promax)), 24)
  expect_equal(unique(orth_grid_1$rotation), c("varimax", "quartimax"))

  expect_is(orth_grid_2, "data.frame")
  expect_named(orth_grid_2, c("method", "init_comm", "criterion", "criterion_type",
                             "abs_eigen", "start_method", "rotation", "k_promax",
                             "normalize", "P_type", "precision", "varimax_type",
                             "k_simplimax"))
  expect_equal(nrow(orth_grid_2), 8)
  expect_equal(sum(is.na(orth_grid_2$varimax_type)), 8)
  expect_equal(sum(is.na(orth_grid_2$k_promax)), 8)
  expect_equal(unique(orth_grid_2$rotation), c("quartimax"))

  expect_equal(orth_grid_2, orth_grid_3)

  expect_is(orth_grid_4, "data.frame")
  expect_named(orth_grid_4, c("method", "init_comm", "criterion", "criterion_type",
                             "abs_eigen", "start_method", "rotation", "k_promax",
                             "normalize", "P_type", "precision", "varimax_type",
                             "k_simplimax"))
  expect_equal(nrow(orth_grid_4), 2)
  expect_equal(sum(is.na(orth_grid_4$k_simplimax)), 2)
  expect_equal(sum(is.na(orth_grid_4$k_promax)), 2)
  expect_equal(sum(is.na(orth_grid_4$varimax_type)), 2)
  expect_equal(unique(orth_grid_4$rotation), c("quartimax"))
  expect_equal(sum(is.na(orth_grid_4$init_comm)), 2)

})


### test .type_grid

tg_ob <- .type_grid("PAF", c("smc", "mac"), .001,
                   c("sum", "max_individual"), c(FALSE, TRUE),
                   NA, "oblique", c(3, 4), TRUE, c("norm", "unnorm"),
                   1e-5, c("kaiser", "svd"), 30)
tg_ob2 <- .type_grid("PAF", c("smc", "mac"), .001,
                    c("sum", "max_individual"), c(FALSE, TRUE),
                    NA, c("oblimin", "promax"), c(3, 4), TRUE, c("norm", "unnorm"),
                    1e-5, c("kaiser", "svd"), 30)
tg_orth <- .type_grid("PAF", c("smc", "mac"), .001,
                    c("sum", "max_individual"), c(FALSE, TRUE),
                    NA, "orthogonal", c(3, 4), TRUE, c("norm", "unnorm"),
                    1e-5, c("kaiser", "svd"), 30)
tg_orth2 <- .type_grid("PAF", c("smc", "mac"), .001,
                      c("sum", "max_individual"), c(FALSE, TRUE),
                      NA, c("varimax", "quartimax"), c(3, 4), TRUE,
                      c("norm", "unnorm"), 1e-5, c("kaiser", "svd"), 30)
tg_nn <- .type_grid("PAF", c("smc", "mac"), .001,
                      c("sum", "max_individual"), c(FALSE, TRUE),
                      NA, "none", c(3, 4), TRUE, c("norm", "unnorm"),
                      1e-5, c("kaiser", "svd"), 30)

test_that(".type_grid works", {
  ### test errors
  expect_error(.type_grid("PAF", NA, NA, NA, NA, NA, c("oblique", "none"), NA,
                          NA, NA, NA, NA, NA))
  expect_error(.type_grid("PAF", NA, NA, NA, NA, NA, c("oblique", "varimax"), NA,
                          NA, NA, NA, NA, NA))
  expect_error(.type_grid("PAF", NA, NA, NA, NA, NA, c("orthogonal", "varimax"),
                          NA, NA, NA, NA, NA, NA))
  expect_error(.type_grid("PAF", NA, NA, NA, NA, NA, c("promax", "varimax"),
                          NA, NA, NA, NA, NA, NA),
               " 'rotation' contains both oblique rotations and orthogonal rotations, but can only average rotations of the same kind. Oblique rotations are 'promax', 'oblimin', 'quartimin', 'simplimax', 'bentlerQ', 'geominQ', and 'bifactorQ'. Orthogonal rotations are 'varimax', 'quartimax', 'equamax', 'bentlerT', 'geominT', and 'bifactorT'.\n")

  expect_is(tg_ob, "data.frame")
  expect_named(tg_ob, c("method", "init_comm", "criterion", "criterion_type",
                              "abs_eigen", "start_method", "rotation", "k_promax",
                              "normalize", "P_type", "precision", "varimax_type",
                              "k_simplimax"))
  expect_equal(nrow(tg_ob), 112)
  expect_equal(sum(is.na(tg_ob$varimax_type)),
               nrow(tg_ob) - sum(tg_ob$rotation == "promax"))
  expect_equal(sum(is.na(tg_ob$k_promax)),
               nrow(tg_ob) - sum(tg_ob$rotation == "promax"))
  expect_equal(sum(is.na(tg_ob$P_type)),
               nrow(tg_ob) - sum(tg_ob$rotation == "promax"))
  expect_equal(sum(is.na(tg_ob$k_simplimax)),
               nrow(tg_ob) - sum(tg_ob$rotation == "simplimax"))
  expect_equal(sort(unique(tg_ob$rotation)),
               sort(c("promax", "oblimin", "quartimin", "simplimax",
                 "bentlerQ", "geominQ", "bifactorQ")))

  expect_is(tg_ob2, "data.frame")
  expect_named(tg_ob2, c("method", "init_comm", "criterion", "criterion_type",
                        "abs_eigen", "start_method", "rotation", "k_promax",
                        "normalize", "P_type", "precision", "varimax_type",
                        "k_simplimax"))
  expect_equal(nrow(tg_ob2), 72)
  expect_equal(sum(is.na(tg_ob2$varimax_type)),
               nrow(tg_ob2) - sum(tg_ob2$rotation == "promax"))
  expect_equal(sum(is.na(tg_ob2$k_promax)),
               nrow(tg_ob2) - sum(tg_ob2$rotation == "promax"))
  expect_equal(sum(is.na(tg_ob2$P_type)),
               nrow(tg_ob2) - sum(tg_ob2$rotation == "promax"))
  expect_equal(sum(is.na(tg_ob2$k_simplimax)),
               nrow(tg_ob2) - sum(tg_ob2$rotation == "simplimax"))
  expect_equal(unique(tg_ob2$rotation),
               c("promax", "oblimin"))

  expect_is(tg_orth, "data.frame")
  expect_named(tg_orth, c("method", "init_comm", "criterion", "criterion_type",
                        "abs_eigen", "start_method", "rotation", "k_promax",
                        "normalize", "P_type", "precision", "varimax_type",
                        "k_simplimax"))
  expect_equal(nrow(tg_orth), 56)
  expect_equal(sum(is.na(tg_orth$varimax_type)),
               nrow(tg_orth) - sum(tg_orth$rotation == "varimax"))
  expect_equal(sort(unique(tg_orth$rotation)),
               sort(c("varimax", "quartimax", "equamax",
                      "bentlerT", "geominT", "bifactorT")))

  expect_is(tg_orth2, "data.frame")
  expect_named(tg_orth2, c("method", "init_comm", "criterion", "criterion_type",
                         "abs_eigen", "start_method", "rotation", "k_promax",
                         "normalize", "P_type", "precision", "varimax_type",
                         "k_simplimax"))
  expect_equal(nrow(tg_orth2), 24)
  expect_equal(sum(is.na(tg_orth2$varimax_type)),
               nrow(tg_orth2) - sum(tg_orth2$rotation == "varimax"))
  expect_equal(unique(tg_orth2$rotation),
               c("varimax", "quartimax"))


  expect_is(tg_nn, "data.frame")
  expect_named(tg_nn, c("method", "init_comm", "criterion", "criterion_type",
                           "abs_eigen", "start_method", "rotation", "k_promax",
                           "normalize", "P_type", "precision", "varimax_type",
                           "k_simplimax"))
  expect_equal(nrow(tg_nn), 8)
  expect_true(all(is.na(tg_nn$varimax_type)))
  expect_true(all(tg_nn$rotation == "none"))
  expect_true(all(is.na(tg_nn$k_promax)))
  expect_true(all(is.na(tg_nn$normalize)))
  expect_true(all(is.na(tg_nn$P_type)))
  expect_true(all(is.na(tg_nn$precision)))
  expect_true(all(is.na(tg_nn$k_simplimax)))

})



rm(efa_pro, efa_temp, x_base, y_base, efa_ml, efa_uls, efa_paf, gof_ml, gof_uls,
   gof_paf, m, q, q_p, dat_unname, dat_unname_2, efa_list, ext_a, efa_list_er,
   ext_er, efa_list_rot, ext_rot, av_mean_NA, av_mean_NA_t01, av_mean,
   av_median, av_median_NA, arr_re_NA, arr_re, obl_grid_1, obl_grid_2,
   obl_grid_3, obl_grid_4, orth_grid_1, orth_grid_2, orth_grid_3, orth_grid_4,
   tg_ob, tg_ob2, tg_orth, tg_orth2, tg_nn)
