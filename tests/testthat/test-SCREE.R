scree_cor <- SCREE(test_models$baseline$cormat)
scree_cor_smc <- SCREE(test_models$baseline$cormat, eigen_type = "SMC")
scree_raw <- SCREE(GRiPS_raw)
# Check with an argument passed to "EFA"
scree_efa_ml<- SCREE(test_models$baseline$cormat, eigen_type = "EFA", method = "ML")

test_that("output class and dimensions are correct", {
  expect_is(scree_cor, "SCREE")
  expect_output(str(scree_cor), "List of 4")
  expect_is(scree_cor_smc, "SCREE")
  expect_output(str(scree_cor_smc), "List of 4")
  expect_is(scree_raw, "SCREE")
  expect_output(str(scree_raw), "List of 4")
  expect_is(scree_efa_ml, "SCREE")
  expect_output(str(scree_efa_ml), "List of 4")
})

test_that("found eigenvalues are correct", {
  expect_equal(sum(scree_cor$eigen_PCA), ncol(test_models$baseline$cormat))
  expect_lt(sum(scree_cor$eigen_SMC), ncol(test_models$baseline$cormat))
  expect_lt(sum(scree_cor$eigen_EFA), ncol(test_models$baseline$cormat))

  expect_equal(sum(scree_raw$eigen_PCA), ncol(GRiPS_raw))
  expect_lt(sum(scree_raw$eigen_SMC), ncol(GRiPS_raw))
  expect_lt(sum(scree_raw$eigen_EFA), ncol(GRiPS_raw))

  expect_lt(sum(scree_raw$eigen_SMC), ncol(test_models$baseline$cormat))
  expect_equal(c(scree_cor_smc$eigen_PCA, scree_cor_smc$eigen_EFA), c(NA, NA))

  expect_lt(sum(scree_raw$eigen_EFA), ncol(test_models$baseline$cormat))
  expect_equal(c(scree_efa_ml$eigen_PCA, scree_efa_ml$eigen_SMC), c(NA, NA))
})


# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors are thrown correctly", {
  expect_error(SCREE(1:5), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n")
  expect_message(SCREE(GRiPS_raw, eigen_type = "PCA"), " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n")
  expect_error(SCREE(dat_sing), " Correlation matrix is singular, no further analyses are performed.\n")
  expect_error(SCREE(cor_sing, N = 10), " Correlation matrix is singular, no further analyses are performed.\n")
  expect_warning(SCREE(cor_nposdef, N = 10), "Matrix was not positive definite, smoothing was done")
})

test_that("settings are returned correctly", {
  expect_named(scree_cor$settings, c("eigen_type", "use", "cor_method",
                                     "n_factors"))
  expect_named(scree_raw$settings, c("eigen_type", "use", "cor_method",
                                     "n_factors"))
  expect_named(scree_cor_smc$settings, c("eigen_type", "use", "cor_method",
                                         "n_factors"))
  expect_named(scree_efa_ml$settings, c("eigen_type", "use", "cor_method",
                                        "n_factors"))

  expect_equal(scree_cor$settings$eigen_type, c("PCA", "SMC", "EFA"))
  expect_equal(scree_raw$settings$eigen_type, c("PCA", "SMC", "EFA"))
  expect_equal(scree_cor_smc$settings$eigen_type, "SMC")
  expect_equal(scree_efa_ml$settings$eigen_type, "EFA")

  expect_equal(scree_cor$settings$use, "pairwise.complete.obs")
  expect_equal(scree_raw$settings$use, "pairwise.complete.obs")
  expect_equal(scree_cor_smc$settings$use, "pairwise.complete.obs")
  expect_equal(scree_efa_ml$settings$use, "pairwise.complete.obs")

  expect_equal(scree_cor$settings$cor_method, "pearson")
  expect_equal(scree_raw$settings$cor_method, "pearson")
  expect_equal(scree_cor_smc$settings$cor_method, "pearson")
  expect_equal(scree_efa_ml$settings$cor_method, "pearson")

  expect_equal(scree_cor$settings$n_factors, 1)
  expect_equal(scree_raw$settings$n_factors, 1)
  expect_equal(scree_cor_smc$settings$n_factors, 1)
  expect_equal(scree_efa_ml$settings$n_factors, 1)

})

rm(scree_cor, scree_cor_smc, scree_raw, scree_efa_ml, x, y, z, dat_sing, cor_sing,
   cor_nposdef)
