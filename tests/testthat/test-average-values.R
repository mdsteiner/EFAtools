# Unit tests for .average_values() (R/averaging.R), which reduces the stacked
# efa_average() arrays to their average, range, and dispersion summaries. Inputs are
# literal arrays with hand-chosen values, so the expected means, medians, trimmed
# means, and NA propagation can be written out exactly.

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
  expect_equal(round(av_mean_NA$fit_indices$sd, 2), c(rep(1.53, 7), 0))
  expect_type(av_mean_NA$fit_indices$range, "double")
  expect_equal(av_mean_NA$fit_indices$range, c(rep(3, 7), 0))
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
  expect_equal(round(av_mean_NA_t01$fit_indices$sd, 2), c(rep(1.53, 7), 0))
  expect_type(av_mean_NA_t01$fit_indices$range, "double")
  expect_equal(av_mean_NA_t01$fit_indices$range, c(rep(3, 7), 0))
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
  expect_equal(round(av_median_NA$fit_indices$sd, 2), c(rep(1.53, 7), 0))
  expect_type(av_median_NA$fit_indices$range, "double")
  expect_equal(av_median_NA$fit_indices$range, c(rep(3, 7), 0))
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


rm(av_mean_NA, av_mean_NA_t01, av_mean, av_median, av_median_NA)
