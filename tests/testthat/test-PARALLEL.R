
pa_cor <- PARALLEL(test_models$baseline$cormat, N = 500)
pa_cor_pca <- PARALLEL(test_models$baseline$cormat, N = 500, eigen_type = "PCA")
pa_raw <- PARALLEL(GRiPS_raw)
pa_nodat <- PARALLEL(N = 20, n_vars = 5)
pa_craw <- PARALLEL(test_models$baseline$cormat, N = 500, eigen_type = "PCA",
                    decision_rule = "crawford")
pa_perc <- PARALLEL(test_models$baseline$cormat, N = 500, eigen_type = "PCA",
                    decision_rule = "percentile")

test_that("output class and dimensions are correct", {
  expect_is(pa_cor, "PARALLEL")
  expect_output(str(pa_cor), "List of 7")
  expect_is(pa_raw, "PARALLEL")
  expect_output(str(pa_raw), "List of 7")
  expect_is(pa_cor_pca, "PARALLEL")
  expect_output(str(pa_cor_pca), "List of 7")
  expect_is(pa_nodat, "PARALLEL")
  expect_output(str(pa_nodat), "List of 7")
  expect_is(pa_craw, "PARALLEL")
  expect_output(str(pa_craw), "List of 7")
  expect_is(pa_perc, "PARALLEL")
  expect_output(str(pa_perc), "List of 7")
})


test_that("found eigenvalues are correct", {
  expect_equal(sum(pa_cor$eigenvalues_PCA[, 1]),
               ncol(test_models$baseline$cormat))
  expect_equal(sum(pa_cor$eigenvalues_PCA[, 2]),
               ncol(test_models$baseline$cormat))
  expect_gt(sum(pa_cor$eigenvalues_PCA[, 3]),
               ncol(test_models$baseline$cormat))
  expect_lt(sum(pa_cor$eigenvalues_SMC[, 1]),
            ncol(test_models$baseline$cormat))
  expect_lt(sum(pa_cor$eigenvalues_SMC[, 2]),
            ncol(test_models$baseline$cormat))
  expect_lt(sum(pa_cor$eigenvalues_EFA[, 1]),
               ncol(test_models$baseline$cormat))
  expect_lt(sum(pa_cor$eigenvalues_EFA[, 2]),
               ncol(test_models$baseline$cormat))

  expect_equal(sum(pa_raw$eigenvalues_PCA[, 1]),
               ncol(GRiPS_raw))
  expect_equal(sum(pa_raw$eigenvalues_PCA[, 2]),
               ncol(GRiPS_raw))
  expect_gt(sum(pa_raw$eigenvalues_PCA[, 3]),
            ncol(GRiPS_raw))
  expect_lt(sum(pa_raw$eigenvalues_SMC[, 1]),
            ncol(GRiPS_raw))
  expect_lt(sum(pa_raw$eigenvalues_SMC[, 2]),
            ncol(GRiPS_raw))
  expect_lt(sum(pa_raw$eigenvalues_EFA[, 1]),
            ncol(GRiPS_raw))
  expect_lt(sum(pa_raw$eigenvalues_EFA[, 2]),
            ncol(GRiPS_raw))

  expect_is(pa_cor$eigenvalues_PCA, "matrix")
  expect_is(pa_cor$eigenvalues_SMC, "matrix")
  expect_is(pa_cor$eigenvalues_EFA, "matrix")

  expect_is(pa_raw$eigenvalues_PCA, "matrix")
  expect_is(pa_raw$eigenvalues_SMC, "matrix")
  expect_is(pa_raw$eigenvalues_EFA, "matrix")

  expect_is(pa_nodat$eigenvalues_PCA, "matrix")
  expect_is(pa_nodat$eigenvalues_SMC, "matrix")
  expect_is(pa_nodat$eigenvalues_EFA, "matrix")

  expect_is(pa_cor_pca$eigenvalues_PCA, "matrix")
  expect_equal(c(pa_cor_pca$eigenvalues_SMC, pa_cor_pca$eigenvalues_EFA), c(NA, NA))

  expect_is(pa_craw$eigenvalues_PCA, "matrix")
  expect_equal(c(pa_craw$eigenvalues_SMC, pa_craw$eigenvalues_EFA), c(NA, NA))

  expect_is(pa_perc$eigenvalues_PCA, "matrix")
  expect_equal(c(pa_perc$eigenvalues_SMC, pa_perc$eigenvalues_EFA), c(NA, NA))

})


test_that("identified number of factors is correct", {
  expect_equal(pa_cor$n_fac_PCA, 3)
  expect_equal(pa_cor$n_fac_SMC, 3)
  expect_equal(pa_cor$n_fac_EFA, 7, tolerance = 2)

  expect_equal(pa_raw$n_fac_PCA, 1)
  expect_equal(pa_raw$n_fac_SMC, 1)
  expect_equal(pa_raw$n_fac_EFA, 3, tolerance = 2)

  expect_equal(c(pa_nodat$n_fac_PCA, pa_nodat$n_fac_SMC, pa_nodat$n_fac_EFA),
               rep(NA, 3))

  expect_equal(pa_cor_pca$n_fac_PCA, 3)
  expect_equal(c(pa_cor_pca$n_fac_SMC, pa_cor_pca$n_fac_EFA), c(NA, NA))

  expect_equal(pa_craw$n_fac_PCA, 3)
  expect_equal(c(pa_craw$n_fac_SMC, pa_craw$n_fac_EFA), c(NA, NA))

  expect_equal(pa_perc$n_fac_PCA, 3)
  expect_equal(c(pa_perc$n_fac_SMC, pa_perc$n_fac_EFA), c(NA, NA))
})


test_that("settings are returned correctly", {
  expect_named(pa_cor$settings, c("x_dat", "N", "n_vars", "n_datasets", "percent",
                                  "eigen_type", "use", "cor_method", "decision_rule",
                                  "n_factors"))
  expect_named(pa_raw$settings, c("x_dat", "N", "n_vars", "n_datasets", "percent",
                                  "eigen_type", "use", "cor_method", "decision_rule",
                                  "n_factors"))
  expect_named(pa_cor_pca$settings, c("x_dat", "N", "n_vars", "n_datasets", "percent",
                                      "eigen_type", "use", "cor_method",
                                      "decision_rule", "n_factors"))
  expect_named(pa_nodat$settings, c("x_dat", "N", "n_vars", "n_datasets", "percent",
                                      "eigen_type", "use", "cor_method",
                                      "decision_rule", "n_factors"))
  expect_named(pa_craw$settings, c("x_dat", "N", "n_vars", "n_datasets", "percent",
                                      "eigen_type", "use", "cor_method",
                                      "decision_rule", "n_factors"))
  expect_named(pa_perc$settings, c("x_dat", "N", "n_vars", "n_datasets", "percent",
                                      "eigen_type", "use", "cor_method",
                                      "decision_rule", "n_factors"))

  expect_true(pa_cor$settings$x_dat)
  expect_true(pa_raw$settings$x_dat)
  expect_true(pa_cor_pca$settings$x_dat)
  expect_false(pa_nodat$settings$x_dat)
  expect_true(pa_craw$settings$x_dat)
  expect_true(pa_perc$settings$x_dat)

  expect_equal(pa_cor$settings$N, 500)
  expect_equal(pa_raw$settings$N, 810)
  expect_equal(pa_cor_pca$settings$N, 500)
  expect_equal(pa_nodat$settings$N, 20)
  expect_equal(pa_craw$settings$N, 500)
  expect_equal(pa_perc$settings$N, 500)

  expect_equal(pa_cor$settings$n_vars, 18)
  expect_equal(pa_raw$settings$n_vars, 8)
  expect_equal(pa_cor_pca$settings$n_vars, 18)
  expect_equal(pa_nodat$settings$n_vars, 5)
  expect_equal(pa_craw$settings$n_vars, 18)
  expect_equal(pa_perc$settings$n_vars, 18)

  expect_equal(pa_cor$settings$n_datasets, 1000)
  expect_equal(pa_raw$settings$n_datasets, 1000)
  expect_equal(pa_cor_pca$settings$n_datasets, 1000)
  expect_equal(pa_nodat$settings$n_datasets, 1000)
  expect_equal(pa_craw$settings$n_datasets, 1000)
  expect_equal(pa_perc$settings$n_datasets, 1000)

  expect_equal(pa_cor$settings$percent, 95)
  expect_equal(pa_raw$settings$percent, 95)
  expect_equal(pa_cor_pca$settings$percent, 95)
  expect_equal(pa_nodat$settings$percent, 95)
  expect_equal(pa_craw$settings$percent, 95)
  expect_equal(pa_perc$settings$percent, 95)

  expect_equal(pa_cor$settings$eigen_type, c("PCA", "SMC", "EFA"))
  expect_equal(pa_raw$settings$eigen_type, c("PCA", "SMC", "EFA"))
  expect_equal(pa_cor_pca$settings$eigen_type, "PCA")
  expect_equal(pa_nodat$settings$eigen_type, c("PCA", "SMC", "EFA"))
  expect_equal(pa_craw$settings$eigen_type, "PCA")
  expect_equal(pa_perc$settings$eigen_type, "PCA")

  expect_equal(pa_cor$settings$use, "pairwise.complete.obs")
  expect_equal(pa_raw$settings$use, "pairwise.complete.obs")
  expect_equal(pa_cor_pca$settings$use, "pairwise.complete.obs")
  expect_equal(pa_nodat$settings$use, "pairwise.complete.obs")
  expect_equal(pa_craw$settings$use, "pairwise.complete.obs")
  expect_equal(pa_perc$settings$use, "pairwise.complete.obs")

  expect_equal(pa_cor$settings$cor_method, "pearson")
  expect_equal(pa_raw$settings$cor_method, "pearson")
  expect_equal(pa_cor_pca$settings$cor_method, "pearson")
  expect_equal(pa_nodat$settings$cor_method, "pearson")
  expect_equal(pa_craw$settings$cor_method, "pearson")
  expect_equal(pa_perc$settings$cor_method, "pearson")

  expect_equal(pa_cor$settings$decision_rule, "means")
  expect_equal(pa_raw$settings$decision_rule, "means")
  expect_equal(pa_cor_pca$settings$decision_rule, "means")
  expect_equal(pa_nodat$settings$decision_rule, "means")
  expect_equal(pa_craw$settings$decision_rule, "crawford")
  expect_equal(pa_perc$settings$decision_rule, "percentile")

  expect_equal(pa_cor$settings$n_factors, 1)
  expect_equal(pa_raw$settings$n_factors, 1)
  expect_equal(pa_cor_pca$settings$n_factors, 1)
  expect_equal(pa_nodat$settings$n_factors, 1)
  expect_equal(pa_craw$settings$n_factors, 1)
  expect_equal(pa_perc$settings$n_factors, 1)

})

# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

burt <- matrix(c(1.00,  0.83,  0.81,  0.80,   0.71, 0.70, 0.54, 0.53,  0.59,  0.24, 0.13,
                 0.83,  1.00,  0.87,  0.62,   0.59, 0.44, 0.58, 0.44,  0.23,  0.45,  0.21,
                 0.81,  0.87,  1.00,  0.63,   0.37, 0.31, 0.30, 0.12,  0.33,  0.33,  0.36,
                 0.80,  0.62,  0.63,  1.00,   0.49, 0.54, 0.30, 0.28,  0.42,  0.29, -0.06,
                 0.71,  0.59,  0.37,  0.49,   1.00, 0.54, 0.34, 0.55,  0.40,  0.19, -0.10,
                 0.70,  0.44,  0.31,  0.54,   0.54, 1.00, 0.50, 0.51,  0.31,  0.11,  0.10,
                 0.54,  0.58,  0.30,  0.30,   0.34, 0.50, 1.00, 0.38,  0.29,  0.21,  0.08,
                 0.53,  0.44,  0.12,  0.28,   0.55, 0.51, 0.38, 1.00,  0.53,  0.10, -0.16,
                 0.59,  0.23,  0.33,  0.42,   0.40, 0.31, 0.29, 0.53,  1.00, -0.09, -0.10,
                 0.24,  0.45,  0.33,  0.29,   0.19, 0.11, 0.21, 0.10, -0.09,  1.00,  0.41,
                 0.13,  0.21,  0.36, -0.06,  -0.10, 0.10, 0.08, -0.16, -0.10, 0.41,  1.00),
               nrow = 11, ncol = 11)


#sim_NA <- data.frame(rnorm(30), rnorm(30), rnorm(30), rep("a", 30))

test_that("errors are thrown correctly", {
  expect_error(PARALLEL(1:5), " 'x' is neither NULL, nor a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data or leave x at NULL.\n")
  expect_warning(suppressMessages(PARALLEL(GRiPS_raw, n_vars = 5)), " n_vars was set and data entered. Taking n_vars from data\n")
  expect_warning(suppressMessages(PARALLEL(GRiPS_raw, N = 20, eigen_type = "PCA")), " 'N' was set and data entered. Taking N from data.\n")
  expect_error(suppressMessages(PARALLEL(N = 500)), ' "n_vars" was not set and could not be taken from data. Please specify n_vars and try again.\n')
  expect_message(PARALLEL(GRiPS_raw, eigen_type = "PCA"), " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n")
  expect_error(PARALLEL(test_models$baseline$cormat, eigen_type = "PCA"), '"N" was not set and could not be taken from data. Please specify N and try again.\n')
  expect_warning(PARALLEL(test_models$baseline$cormat, N = 500,
                          eigen_type = "PCA", decision_rule = "crawford",
                          percent = 80), "decision_rule == 'crawford' is specified, but 95 percentile was not used. Using means instead. To use 'crawford', make sure to specify percent = 95.\n")
  expect_error(PARALLEL(dat_sing), " Correlation matrix is singular, parallel analysis is not possible.\n")
  expect_error(PARALLEL(cor_sing, N = 10), " Correlation matrix is singular, parallel analysis is not possible.\n")
  expect_warning(PARALLEL(burt, N = 100), "Matrix was not positive definite, smoothing was done")
  expect_error(PARALLEL(test_models$baseline$cormat, N = 15), ' "N" was smaller than or equal to the number of variables but must be larger.\n')
  expect_error(PARALLEL(test_models$baseline$cormat, N = 18), ' "N" was smaller than or equal to the number of variables but must be larger.\n')
})

rm(pa_cor, pa_cor_pca, pa_raw, pa_nodat, pa_craw, pa_perc, x, y, z, dat_sing,
   cor_sing, burt)
