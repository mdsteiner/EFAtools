test_that(".numformat works", {
  expect_equal(.numformat(0.23), " .23")
  expect_equal(.numformat(0.2345, digits = 3), " .234")
  expect_equal(.numformat(0.2345, digits = 3, print_zero = TRUE), " 0.234")
})


efa_temp <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500)
efa_pro <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                rotation = "promax")
test_that(".compute_vars works", {
  checkmate::expect_matrix(
    .compute_vars(efa_temp$unrot_loadings, efa_temp$unrot_loadings)
  )
  checkmate::expect_matrix(
    .compute_vars(efa_pro$rot_loadings, efa_pro$unrot_loadings, efa_pro$Phi)
  )
})

x_base <- population_models$loadings$baseline
x_NA <- population_models$loadings$baseline
x_NA[1, 3] <- NA
y_base <- x_base[, c(3,2,1)]
y_NA <- y_base
y_NA[2, 2] <- NA

test_that(".factor_congruence works", {
  checkmate::expect_matrix(.factor_congruence(x_base, y_base))
  expect_equal(sum(.factor_congruence(x_base, y_base)), 3)
  expect_warning(.factor_congruence(x_NA, y_NA, na.rm = FALSE), class = "efa_missing_check")
  expect_warning(.factor_congruence(x_NA, y_NA), class = "efa_missing_complete")
})

test_that(".factor_congruence has (i, j) = cos(x_i, y_j) orientation", {
  # Asymmetric inputs so the congruence matrix differs from its transpose; this
  # pins the orientation the COMPARE / averaging reordering depends on.
  A <- matrix(c(2, 0, 1,
                0, 1, 0,
                0, 0, 3,
                1, 2, 0), nrow = 4, byrow = TRUE)
  B <- matrix(c(0, 1, 0,
                3, 0, 1,
                0, 0, 2,
                1, 1, 0), nrow = 4, byrow = TRUE)
  # Reference computed entrywise from column dot products, independent of the
  # matrix-algebra implementation.
  ref <- sapply(seq_len(ncol(B)), function(j)
    sapply(seq_len(ncol(A)), function(i)
      sum(A[, i] * B[, j]) / sqrt(sum(A[, i]^2) * sum(B[, j]^2))))
  expect_false(isTRUE(all.equal(ref, t(ref))))
  expect_equal(unname(.factor_congruence(A, B)), ref)
})

test_that(".factor_congruence aborts when too few complete cases remain", {
  x_few <- population_models$loadings$baseline
  x_few[-1, 1] <- NA   # only the first row has complete data
  expect_error(
    suppressWarnings(.factor_congruence(x_few, x_base)),
    class = "efa_too_few_complete"
  )
})

set.seed(42)
efa_ml <- suppressWarnings(EFA(cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100),
                                     rnorm(100), rnorm(100)), 3, N = 500,
                               method = "ML"))
efa_uls <- suppressWarnings(EFA(cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100),
                    rnorm(100), rnorm(100)), 3, method = "ULS"))
efa_paf <- suppressWarnings(EFA(cbind(rnorm(100), rnorm(100), rnorm(100), rnorm(100),
                    rnorm(100), rnorm(100)), 3, method = "PAF"))

gof_ml <- .gof(efa_ml$unrot_loadings, efa_ml$orig_R, efa_ml$settings$N,
               "ML", efa_ml$fit_indices$Fm)
w_ml <- tryCatch(.gof(efa_ml$unrot_loadings, efa_ml$orig_R, efa_ml$settings$N,
                      "ML", efa_ml$fit_indices$Fm),
                  error=function(e) e,
                  warning=function(w) w)
gof_uls <- .gof(efa_uls$unrot_loadings, efa_uls$orig_R, efa_uls$settings$N,
               "ULS", efa_uls$fit_indices$Fm)
w_uls <- tryCatch(.gof(efa_uls$unrot_loadings, efa_uls$orig_R, efa_uls$settings$N,
                    "ULS", efa_uls$fit_indices$Fm),
               error=function(e) e,
               warning=function(w) w)

gof_paf <- .gof(efa_paf$unrot_loadings, efa_paf$orig_R, efa_paf$settings$N,
               "PAF", NA)
w_paf <- tryCatch(.gof(efa_paf$unrot_loadings, efa_paf$orig_R, efa_paf$settings$N,
                       "PAF", NA),
                  error=function(e) e,
                  warning=function(w) w)
m <- 6 # n variables
q <- 3 # n factors
test_that(".gof works", {
  expect_type(gof_ml, "list")
  expect_named(gof_ml,
               c("chi", "df", "p_chi", "CAF", "RMSR", "SRMR", "CFI", "TLI",
                 "RMSEA", "RMSEA_LB", "RMSEA_UB", "AIC", "BIC", "ECVI", "Fm",
                 "chi_null", "df_null", "p_null"))
  expect_lt(gof_ml$p_chi, .05)
  expect_equal(gof_ml$CFI, 1)
  expect_equal(gof_ml$RMSEA, 0)
  if (inherits(w_ml, "warning")) {
    expect_equal(gof_ml$CAF, 0, tolerance = .01)
  } else {
    expect_equal(gof_ml$CAF, .5, tolerance = .1)
  }
  expect_equal(gof_ml$df, ((m - q)**2 - (m + q)) / 2)

  expect_type(gof_uls, "list")
  expect_named(gof_uls,
               c("chi", "df", "p_chi", "CAF", "RMSR", "SRMR", "CFI", "TLI",
                 "RMSEA", "RMSEA_LB", "RMSEA_UB", "AIC", "BIC", "ECVI", "Fm",
                 "chi_null", "df_null", "p_null"))
  expect_lt(gof_uls$p_chi, .05)
  expect_equal(gof_uls$CFI, 1)
  expect_equal(gof_uls$RMSEA, 0)
  if (inherits(w_uls, "warning")) {
    expect_equal(gof_uls$CAF, 0, tolerance = .01)
  } else {
    expect_equal(gof_uls$CAF, .5, tolerance = .1)
  }

  expect_equal(gof_uls$df, ((m - q)**2 - (m + q)) / 2)

  expect_type(gof_paf, "list")
  expect_named(gof_paf,
               c("chi", "df", "p_chi", "CAF", "RMSR", "SRMR", "CFI", "TLI",
                 "RMSEA", "RMSEA_LB", "RMSEA_UB", "AIC", "BIC", "ECVI", "Fm",
                 "chi_null", "df_null", "p_null"))
  expect_equal(gof_paf$chi, NA)
  expect_equal(gof_paf$p_chi, NA)
  expect_equal(gof_paf$CFI, NA)
  expect_equal(gof_paf$RMSEA, NA)
  if (inherits(w_paf, "warning")) {
    expect_equal(gof_paf$CAF, 0, tolerance = .01)
  } else {
    expect_equal(gof_paf$CAF, .5, tolerance = .1)
  }
  expect_equal(gof_paf$df, ((m - q)**2 - (m + q)) / 2)
  expect_equal(gof_paf$chi_null, NA)
  expect_equal(gof_paf$df_null, NA)
  expect_equal(gof_paf$p_null, NA)
})


test_that(".compute_caf returns 0 (with warning) when KMO is not computable", {
  # Hollow residual matrix: solve() succeeds (not a try-error) but the inverse
  # has a negative diagonal, so KMO is NaN. CAF must fall back to 0, not NaN.
  delta_hat <- matrix(c(0, .5, .5, .5, 0, .5, .5, .5, 0), 3)
  expect_warning(caf <- .compute_caf(delta_hat), class = "efa_caf_failed")
  expect_equal(caf, 0)
})




test_that(".is_cormat works", {
  expect_equal(.is_cormat(cor(cbind(rnorm(100), rnorm(100)))), TRUE)
  expect_equal(.is_cormat(cbind(rnorm(100), rnorm(100))), FALSE)
  expect_equal(.is_cormat(cbind(rnorm(2), rnorm(2))), FALSE)
  expect_equal(.is_cormat(cbind(c(1, NA, .57, .85))), FALSE)
  expect_equal(.is_cormat(matrix(c(1, .1, .3, 1), ncol = 2)), FALSE)
  expect_error(.is_cormat(matrix(c(1, NA, NA, 1), ncol = 2)),
               class = "efa_cormat_has_na")
})

q_p <- .det_max_factors(8) + 1
test_that(".det_max_factors works", {
  expect_type(.det_max_factors(8), "double")
  expect_lte(((8 - q_p)**2 - (8 + q_p)) / 2, 0)
  expect_equal(.det_max_factors(0), 0)
  expect_equal(.det_max_factors(1), 0)
  expect_equal(.det_max_factors(2), 0)
  expect_equal(.det_max_factors(3), 0)
  expect_gt(.det_max_factors(4), 0)
})


test_that(".decimals works", {
  expect_type(.decimals(8), "double")
  expect_equal(.decimals(8), 0)
  expect_type(.decimals(8), "double")
  expect_error(.decimals("a"), class = "efa_not_numeric")
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
  expect_type(av_mean_NA, "list")
  expect_named(av_mean_NA, c("h2", "loadings", "phi", "vars_accounted",
                              "ind_fac_corres", "fit_indices"))
  expect_type(av_mean_NA$h2, "list")
  expect_named(av_mean_NA$h2, c("average", "sd", "min", "max", "range"))
  expect_type(av_mean_NA$h2$average, "double")
  expect_equal(unname(round(av_mean_NA$h2$average, 2)), rep(2.67, 3))
  expect_named(av_mean_NA$h2$average, paste0("Ind", 1:3))
  expect_equal(unname(av_mean_NA$h2$sd), rep(1.527525, 3), tolerance = .01)
  expect_equal(unname(av_mean_NA$h2$min), rep(1, 3))
  expect_equal(unname(av_mean_NA$h2$max), rep(4, 3))
  expect_type(av_mean_NA$loadings, "list")
  expect_named(av_mean_NA$loadings, c("average", "sd", "min", "max", "range"))
  expect_s3_class(av_mean_NA$loadings$average, "LOADINGS")
  expect_equal(unclass(round(av_mean_NA$loadings$average, 2)), matrix(rep(2.67, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(av_mean_NA$loadings$sd, matrix(rep(1.527525, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))), tolerance = .01)
  expect_equal(unclass(av_mean_NA$loadings$min), matrix(rep(1, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(unclass(av_mean_NA$loadings$max), matrix(rep(4, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(av_mean_NA$loadings$range, matrix(rep(3, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))), tolerance = .01)
  expect_type(av_mean_NA$vars_accounted, "list")
  expect_named(av_mean_NA$vars_accounted, c("average", "sd", "min", "max", "range"))
  checkmate::expect_matrix(av_mean_NA$vars_accounted$average)
  expect_equal(round(av_mean_NA$vars_accounted$average, 2), matrix(rep(2.67, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_mean_NA$vars_accounted$sd, matrix(rep(1.527525, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))), tolerance = .01)
  expect_equal(av_mean_NA$vars_accounted$min, matrix(rep(1, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_mean_NA$vars_accounted$max, matrix(rep(4, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_mean_NA$vars_accounted$range, matrix(rep(3, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))), tolerance = .01)
  expect_equal(av_mean_NA$phi, NA)
  checkmate::expect_matrix(av_mean_NA$ind_fac_corres)
  expect_equal(round(av_mean_NA$ind_fac_corres, 2),
               matrix(c(1, 0, 0, 0, .67, .33, 0, .33, .67), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_s3_class(av_mean_NA$fit_indices, "data.frame")
  expect_named(av_mean_NA$fit_indices, c("index", "average", "sd", "range",
                                          "min", "max"))
  expect_type(av_mean_NA$fit_indices$index, "character")
  expect_equal(av_mean_NA$fit_indices$index, c("chisq", "p_chi", "caf", "cfi",
                                             "rmsea", "aic", "bic", "df"))
  expect_type(av_mean_NA$fit_indices$average, "double")
  expect_equal(round(av_mean_NA$fit_indices$average, 2), c(rep(2.67, 7), 5))
  expect_type(av_mean_NA$fit_indices$sd, "double")
  expect_equal(round(av_mean_NA$fit_indices$sd, 2), c(rep(1.53, 7), 5))
  expect_type(av_mean_NA$fit_indices$range, "double")
  expect_equal(av_mean_NA$fit_indices$range, c(rep(3, 7), 5))
  expect_type(av_mean_NA$fit_indices$min, "double")
  expect_equal(av_mean_NA$fit_indices$min, c(rep(1, 7), 5))
  expect_type(av_mean_NA$fit_indices$max, "double")
  expect_equal(av_mean_NA$fit_indices$max, c(rep(4, 7), 5))


  ### tests for av_mean_NA_t01 with extract_phi = FALSE and trim = .10
  expect_type(av_mean_NA_t01, "list")
  expect_named(av_mean_NA_t01, c("h2", "loadings", "phi", "vars_accounted",
                                  "ind_fac_corres", "fit_indices"))
  expect_type(av_mean_NA_t01$h2, "list")
  expect_named(av_mean_NA_t01$h2, c("average", "sd", "min", "max", "range"))
  expect_type(av_mean_NA_t01$h2$average, "double")
  expect_equal(unname(round(av_mean_NA_t01$h2$average, 2)), rep(3, 3))
  expect_named(av_mean_NA_t01$h2$average, paste0("Ind", 1:3))
  expect_equal(unname(av_mean_NA_t01$h2$sd), rep(1.527525, 3), tolerance = .01)
  expect_equal(unname(av_mean_NA_t01$h2$min), rep(1, 3))
  expect_equal(unname(av_mean_NA_t01$h2$max), rep(4, 3))
  expect_type(av_mean_NA_t01$loadings, "list")
  expect_named(av_mean_NA_t01$loadings, c("average", "sd", "min", "max", "range"))
  expect_s3_class(av_mean_NA_t01$loadings$average, "LOADINGS")
  expect_equal(unclass(round(av_mean_NA_t01$loadings$average, 2)), matrix(rep(3, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(av_mean_NA_t01$loadings$sd, matrix(rep(1.527525, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))), tolerance = .01)
  expect_equal(unclass(av_mean_NA_t01$loadings$min), matrix(rep(1, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(unclass(av_mean_NA_t01$loadings$max), matrix(rep(4, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_type(av_mean_NA_t01$vars_accounted, "list")
  expect_named(av_mean_NA_t01$vars_accounted, c("average", "sd", "min", "max", "range"))
  checkmate::expect_matrix(av_mean_NA_t01$vars_accounted$average)
  expect_equal(round(av_mean_NA_t01$vars_accounted$average, 2), matrix(rep(3, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_mean_NA_t01$vars_accounted$sd, matrix(rep(1.527525, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))), tolerance = .01)
  expect_equal(av_mean_NA_t01$vars_accounted$min, matrix(rep(1, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_mean_NA_t01$vars_accounted$max, matrix(rep(4, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_mean_NA_t01$vars_accounted$range, matrix(rep(3, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))), tolerance = .01)
  expect_equal(av_mean_NA_t01$phi, NA)
  checkmate::expect_matrix(av_mean_NA_t01$ind_fac_corres)
  expect_equal(round(av_mean_NA_t01$ind_fac_corres, 2),
               matrix(c(1, 0, 0, 0, .67, .33, 0, .33, .67), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_s3_class(av_mean_NA_t01$fit_indices, "data.frame")
  expect_named(av_mean_NA_t01$fit_indices, c("index", "average", "sd", "range",
                                          "min", "max"))
  expect_type(av_mean_NA_t01$fit_indices$index, "character")
  expect_equal(av_mean_NA_t01$fit_indices$index, c("chisq", "p_chi", "caf", "cfi",
                                                "rmsea", "aic", "bic", "df"))
  expect_type(av_mean_NA_t01$fit_indices$average, "double")
  expect_equal(round(av_mean_NA_t01$fit_indices$average, 2), c(rep(3, 7), 5))
  expect_type(av_mean_NA_t01$fit_indices$sd, "double")
  expect_equal(round(av_mean_NA_t01$fit_indices$sd, 2), c(rep(1.53, 7), 5))
  expect_type(av_mean_NA_t01$fit_indices$range, "double")
  expect_equal(av_mean_NA_t01$fit_indices$range, c(rep(3, 7), 5))
  expect_type(av_mean_NA_t01$fit_indices$min, "double")
  expect_equal(av_mean_NA_t01$fit_indices$min, c(rep(1, 7), 5))
  expect_type(av_mean_NA_t01$fit_indices$max, "double")
  expect_equal(av_mean_NA_t01$fit_indices$max, c(rep(4, 7), 5))


  ### tests for av_mean with extract_phi = TRUE (only affected output tested)
  expect_type(av_mean$phi, "list")
  checkmate::expect_matrix(av_mean$phi$average)
  expect_equal(round(av_mean$phi$average, 2), matrix(rep(2.67, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))
  checkmate::expect_matrix(av_mean$phi$sd)
  expect_equal(round(av_mean$phi$sd, 2), matrix(rep(1.53, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))
  checkmate::expect_matrix(av_mean$phi$min)
  expect_equal(av_mean$phi$min, matrix(rep(1, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))
  checkmate::expect_matrix(av_mean$phi$max)
  expect_equal(av_mean$phi$max, matrix(rep(4, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))


  ### tests for av_median_NA with extract_phi = FALSE
  expect_type(av_median_NA, "list")
  expect_named(av_median_NA, c("h2", "loadings", "phi", "vars_accounted",
                                "ind_fac_corres", "fit_indices"))
  expect_type(av_median_NA$h2, "list")
  expect_named(av_median_NA$h2, c("average", "sd", "min", "max", "range"))
  expect_type(av_median_NA$h2$average, "double")
  expect_equal(unname(round(av_median_NA$h2$average, 2)), rep(3, 3))
  expect_named(av_median_NA$h2$average, paste0("Ind", 1:3))
  expect_equal(unname(av_median_NA$h2$sd), rep(1.527525, 3), tolerance = .01)
  expect_equal(unname(av_median_NA$h2$min), rep(1, 3))
  expect_equal(unname(av_median_NA$h2$max), rep(4, 3))
  expect_type(av_median_NA$loadings, "list")
  expect_named(av_median_NA$loadings, c("average", "sd", "min", "max", "range"))
  expect_s3_class(av_median_NA$loadings$average, "LOADINGS")
  expect_equal(unclass(round(av_median_NA$loadings$average, 2)), matrix(rep(3, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(av_median_NA$loadings$sd, matrix(rep(1.527525, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))), tolerance = .01)
  expect_equal(unclass(av_median_NA$loadings$min), matrix(rep(1, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_equal(unclass(av_median_NA$loadings$max), matrix(rep(4, 9), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_type(av_median_NA$vars_accounted, "list")
  expect_named(av_median_NA$vars_accounted, c("average", "sd", "min", "max", "range"))
  checkmate::expect_matrix(av_median_NA$vars_accounted$average)
  expect_equal(round(av_median_NA$vars_accounted$average, 2), matrix(rep(3, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_median_NA$vars_accounted$sd, matrix(rep(1.527525, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))), tolerance = .01)
  expect_equal(av_median_NA$vars_accounted$min, matrix(rep(1, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_median_NA$vars_accounted$max, matrix(rep(4, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))))
  expect_equal(av_median_NA$vars_accounted$range, matrix(rep(3, 9), ncol = 3, dimnames = list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))), tolerance = .01)
  expect_equal(av_median_NA$phi, NA)
  checkmate::expect_matrix(av_median_NA$ind_fac_corres)
  expect_equal(round(av_median_NA$ind_fac_corres, 2),
               matrix(c(1, 0, 0, 0, .67, .33, 0, .33, .67), ncol = 3, dimnames = list(paste0("Ind", 1:3), paste0("F", 1:3))))
  expect_s3_class(av_median_NA$fit_indices, "data.frame")
  expect_named(av_median_NA$fit_indices, c("index", "average", "sd", "range",
                                              "min", "max"))
  expect_type(av_median_NA$fit_indices$index, "character")
  expect_equal(av_median_NA$fit_indices$index, c("chisq", "p_chi", "caf", "cfi",
                                                    "rmsea", "aic", "bic", "df"))
  expect_type(av_median_NA$fit_indices$average, "double")
  expect_equal(round(av_median_NA$fit_indices$average, 2), c(rep(3, 7), 5))
  expect_type(av_median_NA$fit_indices$sd, "double")
  expect_equal(round(av_median_NA$fit_indices$sd, 2), c(rep(1.53, 7), 5))
  expect_type(av_median_NA$fit_indices$range, "double")
  expect_equal(av_median_NA$fit_indices$range, c(rep(3, 7), 5))
  expect_type(av_median_NA$fit_indices$min, "double")
  expect_equal(av_median_NA$fit_indices$min, c(rep(1, 7), 5))
  expect_type(av_median_NA$fit_indices$max, "double")
  expect_equal(av_median_NA$fit_indices$max, c(rep(4, 7), 5))


  ### tests for av_median with extract_phi = TRUE (only affected output tested)
  expect_type(av_median$phi, "list")
  checkmate::expect_matrix(av_median$phi$average)
  expect_equal(round(av_median$phi$average, 2), matrix(rep(3, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))
  checkmate::expect_matrix(av_median$phi$sd)
  expect_equal(round(av_median$phi$sd, 2), matrix(rep(1.53, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))
  checkmate::expect_matrix(av_median$phi$min)
  expect_equal(av_median$phi$min, matrix(rep(1, 9), ncol = 3, dimnames = list(paste0("F", 1:3), paste0("F", 1:3))))
  checkmate::expect_matrix(av_median$phi$max)
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
  expect_type(arr_re_NA, "list")
  expect_named(arr_re_NA, c("L", "L_corres", "phi", "vars_accounted"))
  checkmate::expect_array(arr_re_NA$L)
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
  checkmate::expect_array(arr_re_NA$L_corres)
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
  checkmate::expect_array(arr_re_NA$vars_accounted)
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
  expect_type(arr_re, "list")
  expect_named(arr_re, c("L", "L_corres", "phi", "vars_accounted"))
  checkmate::expect_array(arr_re$L)
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
  checkmate::expect_array(arr_re$L_corres)
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
  checkmate::expect_array(arr_re$phi)
  expect_equal(dim(arr_re$phi), c(3, 3, 3))
  expect_equal(arr_re$phi,
               array(c(1, .3, .4, .3, 1, .2, .4, .2, 1,
                       1, .3, .4, .3, 1, .2, .4, .2, 1,
                       1, -.3, -.2, -.3, 1, .4, -.2, .4, 1), c(3, 3, 3)))
  checkmate::expect_array(arr_re$vars_accounted)
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

test_that(".array_reorder uses an optimal permutation alignment", {
  # Two solutions engineered so that a greedy per-row congruence match maps both
  # target factors 1 and 3 onto source column 1 (order 1, 2, 1): column 1 would
  # be duplicated and column 3 dropped, corrupting the average. The optimal
  # (linear-sum-assignment) alignment instead recovers the true permutation.
  L1 <- cbind(c(.9, .8, 0, 0, 0, 0),
              c(0, 0, .9, .8, 0, 0),
              c(0, 0, 0, 0, .9, .8))
  Ln <- cbind(c(.2, .1, 0, 0, .9, .8),
              c(0, 0, .9, .8, 0, 0),
              c(0, 0, 0, 0, .2, .1))

  # The greedy which.max alignment used previously collapses to a non-permutation.
  greedy_order <- apply(abs(.factor_congruence(L1, Ln, skip_checks = TRUE)), 1,
                        which.max)
  expect_equal(greedy_order, c(1, 2, 1))

  L_arr <- array(c(L1, Ln), dim = c(6, 3, 2))
  corres_arr <- array(as.numeric(abs(c(L1, Ln)) > .3), dim = c(6, 3, 2))
  va_arr <- array(1, dim = c(3, 3, 2))

  res <- .array_reorder(vars_accounted = va_arr, L = L_arr, L_corres = corres_arr,
                        phi = NA, extract_phi = FALSE, n_factors = 3)

  aligned <- res$L[, , 2]
  # All three source factors are retained (no duplicated/dropped column) ...
  expect_equal(nrow(unique(round(t(aligned), 8))), 3)
  # ... and the recovered order is the true permutation, leaving Ln in place.
  expect_equal(aligned, Ln)
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
  expect_s3_class(obl_grid_1, "data.frame")
  expect_named(obl_grid_1, c("method", "init_comm", "criterion", "criterion_type",
                            "abs_eigen", "start_method", "rotation", "k_promax",
                            "normalize", "P_type", "precision", "varimax_type",
                            "k_simplimax"))
  expect_equal(nrow(obl_grid_1), 80)
  expect_equal(sum(is.na(obl_grid_1$k_simplimax)), 72)
  expect_equal(sum(is.na(obl_grid_1$k_promax)), 16)
  expect_equal(sum(is.na(obl_grid_1$varimax_type)), 16)
  expect_equal(unique(obl_grid_1$rotation), c("promax", "simplimax", "oblimin"))

  expect_s3_class(obl_grid_2, "data.frame")
  expect_named(obl_grid_2, c("method", "init_comm", "criterion", "criterion_type",
                             "abs_eigen", "start_method", "rotation", "k_promax",
                             "normalize", "P_type", "precision", "varimax_type",
                             "k_simplimax"))
  expect_equal(nrow(obl_grid_2), 16)
  expect_equal(sum(is.na(obl_grid_2$k_promax)), 16)
  expect_equal(sum(is.na(obl_grid_2$varimax_type)), 16)
  expect_equal(unique(obl_grid_2$rotation), c("simplimax", "oblimin"))

  expect_equal(obl_grid_2, obl_grid_3)

  expect_s3_class(obl_grid_4, "data.frame")
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
  expect_s3_class(orth_grid_1, "data.frame")
  expect_named(orth_grid_1, c("method", "init_comm", "criterion", "criterion_type",
                             "abs_eigen", "start_method", "rotation", "k_promax",
                             "normalize", "P_type", "precision", "varimax_type",
                             "k_simplimax"))
  expect_equal(nrow(orth_grid_1), 24)
  expect_equal(sum(is.na(orth_grid_1$varimax_type)), 8)
  expect_equal(sum(is.na(orth_grid_1$k_promax)), 24)
  expect_equal(unique(orth_grid_1$rotation), c("varimax", "quartimax"))

  expect_s3_class(orth_grid_2, "data.frame")
  expect_named(orth_grid_2, c("method", "init_comm", "criterion", "criterion_type",
                             "abs_eigen", "start_method", "rotation", "k_promax",
                             "normalize", "P_type", "precision", "varimax_type",
                             "k_simplimax"))
  expect_equal(nrow(orth_grid_2), 8)
  expect_equal(sum(is.na(orth_grid_2$varimax_type)), 8)
  expect_equal(sum(is.na(orth_grid_2$k_promax)), 8)
  expect_equal(unique(orth_grid_2$rotation), c("quartimax"))

  expect_equal(orth_grid_2, orth_grid_3)

  expect_s3_class(orth_grid_4, "data.frame")
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
               class = "efa_rotation_mismatch")

  expect_s3_class(tg_ob, "data.frame")
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

  expect_s3_class(tg_ob2, "data.frame")
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

  expect_s3_class(tg_orth, "data.frame")
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

  expect_s3_class(tg_orth2, "data.frame")
  expect_named(tg_orth2, c("method", "init_comm", "criterion", "criterion_type",
                         "abs_eigen", "start_method", "rotation", "k_promax",
                         "normalize", "P_type", "precision", "varimax_type",
                         "k_simplimax"))
  expect_equal(nrow(tg_orth2), 24)
  expect_equal(sum(is.na(tg_orth2$varimax_type)),
               nrow(tg_orth2) - sum(tg_orth2$rotation == "varimax"))
  expect_equal(unique(tg_orth2$rotation),
               c("varimax", "quartimax"))


  expect_s3_class(tg_nn, "data.frame")
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


### test .rotation_family

test_that(".rotation_family works", {
  expect_identical(.rotation_family("none"), "none")

  # representative names hardcoded so the test catches a rotation placed in the
  # wrong family vector (including the special-cased varimax/promax and the
  # geomin/bentler/bifactor T-vs-Q pairs)
  expect_identical(.rotation_family("varimax"), "orthogonal")
  expect_identical(.rotation_family("geominT"), "orthogonal")
  expect_identical(.rotation_family("bifactorT"), "orthogonal")
  expect_identical(.rotation_family("promax"), "oblique")
  expect_identical(.rotation_family("geominQ"), "oblique")
  expect_identical(.rotation_family("bifactorQ"), "oblique")

  # and every canonical name dispatches to its family
  for (rot in .orth_rotations) {
    expect_identical(.rotation_family(rot), "orthogonal")
  }

  for (rot in .oblq_rotations) {
    expect_identical(.rotation_family(rot), "oblique")
  }

  expect_error(.rotation_family("bogus"), class = "efa_unknown_rotation")
})



rm(efa_pro, efa_temp, x_base, y_base, efa_ml, efa_uls, efa_paf, gof_ml, gof_uls,
   gof_paf, m, q, q_p, efa_list, ext_a, efa_list_er,
   ext_er, efa_list_rot, ext_rot, av_mean_NA, av_mean_NA_t01, av_mean,
   av_median, av_median_NA, arr_re_NA, arr_re, obl_grid_1, obl_grid_2,
   obl_grid_3, obl_grid_4, orth_grid_1, orth_grid_2, orth_grid_3, orth_grid_4,
   tg_ob, tg_ob2, tg_orth, tg_orth2, tg_nn)
