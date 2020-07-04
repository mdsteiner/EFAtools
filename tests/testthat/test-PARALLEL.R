
pa_cor <- PARALLEL(test_models$baseline$cormat, N = 500)
pa_cor_pca <- PARALLEL(test_models$baseline$cormat, N = 500, eigen_type = "PCA")
pa_raw <- PARALLEL(GRiPS_raw)

test_that("output class and dimensions are correct", {
  expect_is(pa_cor, "PARALLEL")
  expect_output(str(pa_cor), "List of 7")
  expect_is(pa_raw, "PARALLEL")
  expect_output(str(pa_raw), "List of 7")
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

  expect_is(pa_cor_pca$eigenvalues_PCA, "matrix")
  expect_equal(c(pa_cor_pca$eigenvalues_SMC, pa_cor_pca$eigenvalues_EFA), c(NA, NA))

})


test_that("identified number of factors is correct", {
  expect_equal(pa_cor$n_fac_PCA, 3)
  expect_equal(pa_cor$n_fac_SMC, 3)
  expect_equal(pa_cor$n_fac_EFA, 7, tolerance = 2)

  expect_equal(pa_raw$n_fac_PCA, 1)
  expect_equal(pa_raw$n_fac_SMC, 1)
  expect_equal(pa_raw$n_fac_EFA, 3, tolerance = 2)

  expect_equal(pa_cor_pca$n_fac_PCA, 3)
  expect_equal(c(pa_cor_pca$n_fac_SMC, pa_cor_pca$n_fac_EFA), c(NA, NA))
})

# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)


test_that("errors are thrown correctly", {
  expect_warning(PARALLEL(GRiPS_raw, N = 20, eigen_type = "PCA"), " 'N' was set and data entered. Taking N from data.\n")
  expect_message(PARALLEL(GRiPS_raw, eigen_type = "PCA"), " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n")
  expect_error(PARALLEL(test_models$baseline$cormat, eigen_type = "PCA"), '"N" was not set and could not be taken from data. Please specify N and try again.\n')
  expect_warning(PARALLEL(test_models$baseline$cormat, N = 500,
                          eigen_type = "PCA", decision_rule = "crawford",
                          percent = 80), "decision_rule == 'crawford' is specified, but 95 percentile was not used. Using means instead. To use 'crawford', make sure to specify percent = 95.\n")
  expect_error(PARALLEL(dat_sing), " Correlation matrix is singular, parallel analysis is not possible.\n")
  expect_error(PARALLEL(cor_sing, N = 10), " Correlation matrix is singular, parallel analysis is not possible.\n")
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

  expect_true(pa_cor$settings$x_dat)
  expect_true(pa_raw$settings$x_dat)
  expect_true(pa_cor_pca$settings$x_dat)

  expect_equal(pa_cor$settings$N, 500)
  expect_equal(pa_raw$settings$N, 810)
  expect_equal(pa_cor_pca$settings$N, 500)

  expect_equal(pa_cor$settings$n_vars, 18)
  expect_equal(pa_raw$settings$n_vars, 8)
  expect_equal(pa_cor_pca$settings$n_vars, 18)

  expect_equal(pa_cor$settings$n_datasets, 1000)
  expect_equal(pa_raw$settings$n_datasets, 1000)
  expect_equal(pa_cor_pca$settings$n_datasets, 1000)

  expect_equal(pa_cor$settings$percent, 95)
  expect_equal(pa_raw$settings$percent, 95)
  expect_equal(pa_cor_pca$settings$percent, 95)

  expect_equal(pa_cor$settings$eigen_type, c("PCA", "SMC", "EFA"))
  expect_equal(pa_raw$settings$eigen_type, c("PCA", "SMC", "EFA"))
  expect_equal(pa_cor_pca$settings$eigen_type, "PCA")

  expect_equal(pa_cor$settings$use, "pairwise.complete.obs")
  expect_equal(pa_raw$settings$use, "pairwise.complete.obs")
  expect_equal(pa_cor_pca$settings$use, "pairwise.complete.obs")

  expect_equal(pa_cor$settings$cor_method, "pearson")
  expect_equal(pa_raw$settings$cor_method, "pearson")
  expect_equal(pa_cor_pca$settings$cor_method, "pearson")

  expect_equal(pa_cor$settings$decision_rule, "means")
  expect_equal(pa_raw$settings$decision_rule, "means")
  expect_equal(pa_cor_pca$settings$decision_rule, "means")

  expect_equal(pa_cor$settings$n_factors, 1)
  expect_equal(pa_raw$settings$n_factors, 1)
  expect_equal(pa_cor_pca$settings$n_factors, 1)

})

rm(pa_cor, pa_cor_pca, pa_raw, x, y, z, dat_sing, cor_sing)
