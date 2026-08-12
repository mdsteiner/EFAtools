# Unit tests for .average_values() (R/averaging.R), which reduces the stacked
# efa_average() arrays to their average, range, and dispersion summaries. Inputs are
# literal arrays with hand-chosen values, so the expected means, medians, trimmed
# means, and NA propagation can be written out exactly.

# Five solutions per cell, all carrying the same five values. They are chosen so the
# three averaging modes are mutually distinguishable: the untrimmed mean is 6, trim = .2
# drops the most extreme value at each end and averages the remaining three to 4, and the
# median is 3. Three values per cell could not tell them apart -- with n = 3 nothing is
# trimmed until trim reaches .5, where the trimmed mean is exactly the median.
av_values <- c(1, 2, 3, 7, 17)

# Indicator-to-factor correspondence: four solutions assign item 2 to F2 and item 3 to
# F3, the fifth swaps them, so the averaged correspondence is .8 / .2. This summary is a
# plain proportion and is never trimmed.
av_corres_plain <- c(1, 0, 0, 0, 1, 0, 0, 0, 1)
av_corres_swap <- c(1, 0, 0, 0, 0, 1, 0, 1, 0)

av_fit_names <- c("chisq", "p_chi", "caf", "cfi", "rmsea", "aic", "bic")
av_for_grid <- as.data.frame(
  matrix(rep(av_values, times = length(av_fit_names)), nrow = length(av_values),
         dimnames = list(NULL, av_fit_names))
)
av_h2 <- matrix(rep(av_values, times = 3), nrow = length(av_values), ncol = 3)

# One .average_values() call. Every argument is fixed except the three the fixtures
# actually vary, so what each fixture tests is visible on its own line below.
av_fixture <- function(extract_phi = FALSE, averaging = "mean", trim = 0,
                       for_grid = av_for_grid, h2 = av_h2) {
  stacked <- array(rep(av_values, each = 9), c(3, 3, 5))
  .average_values(
    L = stacked,
    L_corres = array(c(av_corres_plain, av_corres_swap,
                       rep(av_corres_plain, times = 3)), c(3, 3, 5)),
    vars_accounted = stacked,
    h2 = h2,
    phi = if (isTRUE(extract_phi)) stacked else NA,
    extract_phi = extract_phi,
    averaging = averaging,
    trim = trim,
    for_grid = for_grid,
    df = 5,
    ind_names = paste0("Ind", 1:3)
  )
}

av_mean_NA <- av_fixture(extract_phi = FALSE, averaging = "mean", trim = 0)
av_mean_NA_t20 <- av_fixture(extract_phi = FALSE, averaging = "mean", trim = .2)
av_mean <- av_fixture(extract_phi = TRUE, averaging = "mean", trim = 0)
av_median <- av_fixture(extract_phi = TRUE, averaging = "median", trim = 0)
av_median_NA <- av_fixture(extract_phi = FALSE, averaging = "median", trim = .2)

# Dimnames and the closed-form dispersion of av_values, written once.
av_ind_dn <- list(paste0("Ind", 1:3), paste0("F", 1:3))
av_var_dn <- list(c("SS loadings", "Prop Tot Var", "Prop Comm Var"), paste0("F", 1:3))
av_fac_dn <- list(paste0("F", 1:3), paste0("F", 1:3))
# sum((av_values - 6)^2) = 172, over n - 1 = 4
av_sd <- sqrt(43)


test_that(".average_values works", {
  ### tests for av_mean_NA with extract_phi = FALSE and trim = 0
  expect_type(av_mean_NA, "list")
  expect_named(av_mean_NA, c("h2", "loadings", "phi", "vars_accounted",
                              "ind_fac_corres", "fit_indices"))
  expect_type(av_mean_NA$h2, "list")
  expect_named(av_mean_NA$h2, c("average", "sd", "min", "max", "range"))
  expect_type(av_mean_NA$h2$average, "double")
  expect_equal(unname(av_mean_NA$h2$average), rep(6, 3))
  expect_named(av_mean_NA$h2$average, paste0("Ind", 1:3))
  expect_equal(unname(av_mean_NA$h2$sd), rep(av_sd, 3))
  expect_equal(unname(av_mean_NA$h2$min), rep(1, 3))
  expect_equal(unname(av_mean_NA$h2$max), rep(17, 3))
  expect_type(av_mean_NA$loadings, "list")
  expect_named(av_mean_NA$loadings, c("average", "sd", "min", "max", "range"))
  expect_s3_class(av_mean_NA$loadings$average, "LOADINGS")
  expect_equal(unclass(av_mean_NA$loadings$average), matrix(rep(6, 9), ncol = 3, dimnames = av_ind_dn))
  expect_equal(av_mean_NA$loadings$sd, matrix(rep(av_sd, 9), ncol = 3, dimnames = av_ind_dn))
  expect_equal(unclass(av_mean_NA$loadings$min), matrix(rep(1, 9), ncol = 3, dimnames = av_ind_dn))
  expect_equal(unclass(av_mean_NA$loadings$max), matrix(rep(17, 9), ncol = 3, dimnames = av_ind_dn))
  expect_equal(av_mean_NA$loadings$range, matrix(rep(16, 9), ncol = 3, dimnames = av_ind_dn))
  expect_type(av_mean_NA$vars_accounted, "list")
  expect_named(av_mean_NA$vars_accounted, c("average", "sd", "min", "max", "range"))
  checkmate::expect_matrix(av_mean_NA$vars_accounted$average)
  expect_equal(av_mean_NA$vars_accounted$average, matrix(rep(6, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_mean_NA$vars_accounted$sd, matrix(rep(av_sd, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_mean_NA$vars_accounted$min, matrix(rep(1, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_mean_NA$vars_accounted$max, matrix(rep(17, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_mean_NA$vars_accounted$range, matrix(rep(16, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_mean_NA$phi, NA)
  checkmate::expect_matrix(av_mean_NA$ind_fac_corres)
  expect_equal(av_mean_NA$ind_fac_corres,
               matrix(c(1, 0, 0, 0, .8, .2, 0, .2, .8), ncol = 3, dimnames = av_ind_dn))
  expect_s3_class(av_mean_NA$fit_indices, "data.frame")
  expect_named(av_mean_NA$fit_indices, c("index", "average", "sd", "range",
                                          "min", "max"))
  expect_type(av_mean_NA$fit_indices$index, "character")
  expect_equal(av_mean_NA$fit_indices$index, c(av_fit_names, "df"))
  expect_type(av_mean_NA$fit_indices$average, "double")
  expect_equal(av_mean_NA$fit_indices$average, c(rep(6, 7), 5))
  expect_type(av_mean_NA$fit_indices$sd, "double")
  expect_equal(av_mean_NA$fit_indices$sd, c(rep(av_sd, 7), 0))
  expect_type(av_mean_NA$fit_indices$range, "double")
  expect_equal(av_mean_NA$fit_indices$range, c(rep(16, 7), 0))
  expect_type(av_mean_NA$fit_indices$min, "double")
  expect_equal(av_mean_NA$fit_indices$min, c(rep(1, 7), 5))
  expect_type(av_mean_NA$fit_indices$max, "double")
  expect_equal(av_mean_NA$fit_indices$max, c(rep(17, 7), 5))


  ### tests for av_mean_NA_t20 with extract_phi = FALSE and trim = .20
  expect_type(av_mean_NA_t20, "list")
  expect_named(av_mean_NA_t20, c("h2", "loadings", "phi", "vars_accounted",
                                  "ind_fac_corres", "fit_indices"))
  expect_type(av_mean_NA_t20$h2, "list")
  expect_named(av_mean_NA_t20$h2, c("average", "sd", "min", "max", "range"))
  expect_type(av_mean_NA_t20$h2$average, "double")
  expect_equal(unname(av_mean_NA_t20$h2$average), rep(4, 3))
  expect_named(av_mean_NA_t20$h2$average, paste0("Ind", 1:3))
  # only the average is trimmed; the dispersion summaries still describe all five values
  expect_equal(unname(av_mean_NA_t20$h2$sd), rep(av_sd, 3))
  expect_equal(unname(av_mean_NA_t20$h2$min), rep(1, 3))
  expect_equal(unname(av_mean_NA_t20$h2$max), rep(17, 3))
  expect_type(av_mean_NA_t20$loadings, "list")
  expect_named(av_mean_NA_t20$loadings, c("average", "sd", "min", "max", "range"))
  expect_s3_class(av_mean_NA_t20$loadings$average, "LOADINGS")
  expect_equal(unclass(av_mean_NA_t20$loadings$average), matrix(rep(4, 9), ncol = 3, dimnames = av_ind_dn))
  expect_equal(av_mean_NA_t20$loadings$sd, matrix(rep(av_sd, 9), ncol = 3, dimnames = av_ind_dn))
  expect_equal(unclass(av_mean_NA_t20$loadings$min), matrix(rep(1, 9), ncol = 3, dimnames = av_ind_dn))
  expect_equal(unclass(av_mean_NA_t20$loadings$max), matrix(rep(17, 9), ncol = 3, dimnames = av_ind_dn))
  expect_type(av_mean_NA_t20$vars_accounted, "list")
  expect_named(av_mean_NA_t20$vars_accounted, c("average", "sd", "min", "max", "range"))
  checkmate::expect_matrix(av_mean_NA_t20$vars_accounted$average)
  expect_equal(av_mean_NA_t20$vars_accounted$average, matrix(rep(4, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_mean_NA_t20$vars_accounted$sd, matrix(rep(av_sd, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_mean_NA_t20$vars_accounted$min, matrix(rep(1, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_mean_NA_t20$vars_accounted$max, matrix(rep(17, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_mean_NA_t20$vars_accounted$range, matrix(rep(16, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_mean_NA_t20$phi, NA)
  checkmate::expect_matrix(av_mean_NA_t20$ind_fac_corres)
  expect_equal(av_mean_NA_t20$ind_fac_corres,
               matrix(c(1, 0, 0, 0, .8, .2, 0, .2, .8), ncol = 3, dimnames = av_ind_dn))
  expect_s3_class(av_mean_NA_t20$fit_indices, "data.frame")
  expect_named(av_mean_NA_t20$fit_indices, c("index", "average", "sd", "range",
                                          "min", "max"))
  expect_type(av_mean_NA_t20$fit_indices$index, "character")
  expect_equal(av_mean_NA_t20$fit_indices$index, c(av_fit_names, "df"))
  expect_type(av_mean_NA_t20$fit_indices$average, "double")
  expect_equal(av_mean_NA_t20$fit_indices$average, c(rep(4, 7), 5))
  expect_type(av_mean_NA_t20$fit_indices$sd, "double")
  expect_equal(av_mean_NA_t20$fit_indices$sd, c(rep(av_sd, 7), 0))
  expect_type(av_mean_NA_t20$fit_indices$range, "double")
  expect_equal(av_mean_NA_t20$fit_indices$range, c(rep(16, 7), 0))
  expect_type(av_mean_NA_t20$fit_indices$min, "double")
  expect_equal(av_mean_NA_t20$fit_indices$min, c(rep(1, 7), 5))
  expect_type(av_mean_NA_t20$fit_indices$max, "double")
  expect_equal(av_mean_NA_t20$fit_indices$max, c(rep(17, 7), 5))


  ### tests for av_mean with extract_phi = TRUE (only affected output tested)
  expect_type(av_mean$phi, "list")
  checkmate::expect_matrix(av_mean$phi$average)
  expect_equal(av_mean$phi$average, matrix(rep(6, 9), ncol = 3, dimnames = av_fac_dn))
  checkmate::expect_matrix(av_mean$phi$sd)
  expect_equal(av_mean$phi$sd, matrix(rep(av_sd, 9), ncol = 3, dimnames = av_fac_dn))
  checkmate::expect_matrix(av_mean$phi$min)
  expect_equal(av_mean$phi$min, matrix(rep(1, 9), ncol = 3, dimnames = av_fac_dn))
  checkmate::expect_matrix(av_mean$phi$max)
  expect_equal(av_mean$phi$max, matrix(rep(17, 9), ncol = 3, dimnames = av_fac_dn))


  ### tests for av_median_NA with extract_phi = FALSE
  expect_type(av_median_NA, "list")
  expect_named(av_median_NA, c("h2", "loadings", "phi", "vars_accounted",
                                "ind_fac_corres", "fit_indices"))
  expect_type(av_median_NA$h2, "list")
  expect_named(av_median_NA$h2, c("average", "sd", "min", "max", "range"))
  expect_type(av_median_NA$h2$average, "double")
  expect_equal(unname(av_median_NA$h2$average), rep(3, 3))
  expect_named(av_median_NA$h2$average, paste0("Ind", 1:3))
  expect_equal(unname(av_median_NA$h2$sd), rep(av_sd, 3))
  expect_equal(unname(av_median_NA$h2$min), rep(1, 3))
  expect_equal(unname(av_median_NA$h2$max), rep(17, 3))
  expect_type(av_median_NA$loadings, "list")
  expect_named(av_median_NA$loadings, c("average", "sd", "min", "max", "range"))
  expect_s3_class(av_median_NA$loadings$average, "LOADINGS")
  expect_equal(unclass(av_median_NA$loadings$average), matrix(rep(3, 9), ncol = 3, dimnames = av_ind_dn))
  expect_equal(av_median_NA$loadings$sd, matrix(rep(av_sd, 9), ncol = 3, dimnames = av_ind_dn))
  expect_equal(unclass(av_median_NA$loadings$min), matrix(rep(1, 9), ncol = 3, dimnames = av_ind_dn))
  expect_equal(unclass(av_median_NA$loadings$max), matrix(rep(17, 9), ncol = 3, dimnames = av_ind_dn))
  expect_type(av_median_NA$vars_accounted, "list")
  expect_named(av_median_NA$vars_accounted, c("average", "sd", "min", "max", "range"))
  checkmate::expect_matrix(av_median_NA$vars_accounted$average)
  expect_equal(av_median_NA$vars_accounted$average, matrix(rep(3, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_median_NA$vars_accounted$sd, matrix(rep(av_sd, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_median_NA$vars_accounted$min, matrix(rep(1, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_median_NA$vars_accounted$max, matrix(rep(17, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_median_NA$vars_accounted$range, matrix(rep(16, 9), ncol = 3, dimnames = av_var_dn))
  expect_equal(av_median_NA$phi, NA)
  checkmate::expect_matrix(av_median_NA$ind_fac_corres)
  expect_equal(av_median_NA$ind_fac_corres,
               matrix(c(1, 0, 0, 0, .8, .2, 0, .2, .8), ncol = 3, dimnames = av_ind_dn))
  expect_s3_class(av_median_NA$fit_indices, "data.frame")
  expect_named(av_median_NA$fit_indices, c("index", "average", "sd", "range",
                                              "min", "max"))
  expect_type(av_median_NA$fit_indices$index, "character")
  expect_equal(av_median_NA$fit_indices$index, c(av_fit_names, "df"))
  expect_type(av_median_NA$fit_indices$average, "double")
  expect_equal(av_median_NA$fit_indices$average, c(rep(3, 7), 5))
  expect_type(av_median_NA$fit_indices$sd, "double")
  expect_equal(av_median_NA$fit_indices$sd, c(rep(av_sd, 7), 0))
  expect_type(av_median_NA$fit_indices$range, "double")
  expect_equal(av_median_NA$fit_indices$range, c(rep(16, 7), 0))
  expect_type(av_median_NA$fit_indices$min, "double")
  expect_equal(av_median_NA$fit_indices$min, c(rep(1, 7), 5))
  expect_type(av_median_NA$fit_indices$max, "double")
  expect_equal(av_median_NA$fit_indices$max, c(rep(17, 7), 5))


  ### tests for av_median with extract_phi = TRUE (only affected output tested)
  expect_type(av_median$phi, "list")
  checkmate::expect_matrix(av_median$phi$average)
  expect_equal(av_median$phi$average, matrix(rep(3, 9), ncol = 3, dimnames = av_fac_dn))
  checkmate::expect_matrix(av_median$phi$sd)
  expect_equal(av_median$phi$sd, matrix(rep(av_sd, 9), ncol = 3, dimnames = av_fac_dn))
  checkmate::expect_matrix(av_median$phi$min)
  expect_equal(av_median$phi$min, matrix(rep(1, 9), ncol = 3, dimnames = av_fac_dn))
  checkmate::expect_matrix(av_median$phi$max)
  expect_equal(av_median$phi$max, matrix(rep(17, 9), ncol = 3, dimnames = av_fac_dn))


})


test_that("trim moves the mean and is inert under the median", {
  # The point of the five-value fixture: a trimmed mean that is neither the untrimmed
  # mean nor the median. Dropping trim from the mean() call, or mis-scaling it, would
  # collapse one of these differences.
  expect_false(isTRUE(all.equal(av_mean_NA_t20$h2$average, av_mean_NA$h2$average)))
  expect_false(isTRUE(all.equal(av_mean_NA_t20$h2$average, av_median_NA$h2$average)))
  expect_false(isTRUE(all.equal(av_mean_NA$h2$average, av_median_NA$h2$average)))

  # trim is a parameter of mean() alone, so the median route ignores it entirely
  expect_equal(av_median_NA$h2$average, av_median$h2$average)
  expect_equal(unclass(av_median_NA$loadings$average),
               unclass(av_median$loadings$average))
  expect_equal(av_median_NA$fit_indices$average, av_median$fit_indices$average)

  # trimming changes only the average, never the dispersion summaries
  expect_equal(av_mean_NA_t20$fit_indices$sd, av_mean_NA$fit_indices$sd)
  expect_equal(av_mean_NA_t20$fit_indices$min, av_mean_NA$fit_indices$min)
  expect_equal(av_mean_NA_t20$fit_indices$max, av_mean_NA$fit_indices$max)
})


test_that("a summary with nothing to summarise is NA and raises no warning", {
  # A PAF-only grid carries no chi-square-based fit index, so those columns are missing
  # for every solution -- the routine case. base min()/max() would warn there and return
  # -Inf/Inf, which must never reach the reported range.
  no_chisq <- av_for_grid
  no_chisq$chisq <- NA_real_
  no_h2 <- av_h2
  no_h2[, 2] <- NA_real_

  expect_no_warning(empty <- av_fixture(for_grid = no_chisq, h2 = no_h2))

  chisq_row <- empty$fit_indices$index == "chisq"
  expect_true(is.na(empty$fit_indices$average[chisq_row]))
  expect_true(is.na(empty$fit_indices$min[chisq_row]))
  expect_true(is.na(empty$fit_indices$max[chisq_row]))
  expect_true(is.na(empty$fit_indices$range[chisq_row]))

  expect_true(is.na(empty$h2$min[[2]]))
  expect_true(is.na(empty$h2$max[[2]]))
  expect_true(is.na(empty$h2$range[[2]]))

  # the indices that do have values are unaffected
  expect_equal(empty$fit_indices$min[empty$fit_indices$index == "caf"], 1)
  expect_equal(empty$fit_indices$max[empty$fit_indices$index == "caf"], 17)
  expect_equal(unname(empty$h2$min[c(1, 3)]), c(1, 1))
})


rm(av_mean_NA, av_mean_NA_t20, av_mean, av_median, av_median_NA, av_fixture,
   av_values, av_corres_plain, av_corres_swap, av_fit_names, av_for_grid, av_h2,
   av_ind_dn, av_var_dn, av_fac_dn, av_sd)
