kgc_cor <- KGC(test_models$baseline$cormat)
kgc_cor_smc <- KGC(test_models$baseline$cormat, eigen_type = "SMC")
kgc_raw <- KGC(GRiPS_raw)
# Check with an argument passed to "EFA"
kgc_efa_ml<- KGC(test_models$baseline$cormat, eigen_type = "EFA", method = "ML")

test_that("output class and dimensions are correct", {
  expect_is(kgc_cor, "KGC")
  expect_output(str(kgc_cor), "List of 7")
  expect_is(kgc_cor_smc, "KGC")
  expect_output(str(kgc_cor_smc), "List of 7")
  expect_is(kgc_raw, "KGC")
  expect_output(str(kgc_raw), "List of 7")
  expect_is(kgc_efa_ml, "KGC")
  expect_output(str(kgc_efa_ml), "List of 7")
})

test_that("found eigenvalues are correct", {
  expect_equal(sum(kgc_cor$eigen_PCA), ncol(test_models$baseline$cormat))
  expect_lt(sum(kgc_cor$eigen_SMC), ncol(test_models$baseline$cormat))
  expect_lt(sum(kgc_cor$eigen_EFA), ncol(test_models$baseline$cormat))

  expect_equal(sum(kgc_raw$eigen_PCA), ncol(GRiPS_raw))
  expect_lt(sum(kgc_raw$eigen_SMC), ncol(GRiPS_raw))
  expect_lt(sum(kgc_raw$eigen_EFA), ncol(GRiPS_raw))

  expect_lt(sum(kgc_raw$eigen_SMC), ncol(test_models$baseline$cormat))
  expect_equal(c(kgc_cor_smc$eigen_PCA, kgc_cor_smc$eigen_EFA), c(NA, NA))

  expect_lt(sum(kgc_raw$eigen_EFA), ncol(test_models$baseline$cormat))
  expect_equal(c(kgc_efa_ml$eigen_PCA, kgc_efa_ml$eigen_SMC), c(NA, NA))
})

test_that("identified number of factors is correct", {
  expect_equal(kgc_cor$n_fac_PCA, 3)
  expect_equal(kgc_cor$n_fac_SMC, 1)
  expect_equal(kgc_cor$n_fac_EFA, 1)

  expect_equal(kgc_raw$n_fac_PCA, 1)
  expect_equal(kgc_raw$n_fac_SMC, 1)
  expect_equal(kgc_raw$n_fac_EFA, 1)

  expect_equal(kgc_cor_smc$n_fac_SMC, 1)
  expect_equal(c(kgc_cor_smc$n_fac_PCA, kgc_cor_smc$n_fac_EFA), c(NA, NA))

  expect_equal(kgc_efa_ml$n_fac_EFA, 1)
  expect_equal(c(kgc_efa_ml$n_fac_PCA, kgc_efa_ml$n_fac_SMC), c(NA, NA))
})

# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors are thrown correctly", {
  expect_error(KGC(1:5), " 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n")
  expect_message(KGC(GRiPS_raw, eigen_type = "PCA"), " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n")
  expect_error(KGC(dat_sing), " Correlation matrix is singular, no further analyses are performed.\n")
  expect_error(KGC(cor_sing, N = 10), " Correlation matrix is singular, no further analyses are performed.\n")
  expect_warning(KGC(cor_nposdef, N = 10), "Matrix was not positive definite, smoothing was done")
})

test_that("settings are returned correctly", {
  expect_named(kgc_cor$settings, c("eigen_type", "use", "cor_method", "n_factors"))
  expect_named(kgc_raw$settings, c("eigen_type", "use", "cor_method", "n_factors"))
  expect_named(kgc_cor_smc$settings, c("eigen_type", "use", "cor_method", "n_factors"))
  expect_named(kgc_efa_ml$settings, c("eigen_type", "use", "cor_method", "n_factors"))

  expect_equal(kgc_cor$settings$eigen_type, c("PCA", "SMC", "EFA"))
  expect_equal(kgc_raw$settings$eigen_type, c("PCA", "SMC", "EFA"))
  expect_equal(kgc_cor_smc$settings$eigen_type, "SMC")
  expect_equal(kgc_efa_ml$settings$eigen_type, "EFA")

  expect_equal(kgc_cor$settings$use, "pairwise.complete.obs")
  expect_equal(kgc_raw$settings$use, "pairwise.complete.obs")
  expect_equal(kgc_cor_smc$settings$use, "pairwise.complete.obs")
  expect_equal(kgc_efa_ml$settings$use, "pairwise.complete.obs")

  expect_equal(kgc_cor$settings$cor_method, "pearson")
  expect_equal(kgc_raw$settings$cor_method, "pearson")
  expect_equal(kgc_cor_smc$settings$cor_method, "pearson")
  expect_equal(kgc_efa_ml$settings$cor_method, "pearson")

  expect_equal(kgc_cor$settings$n_factors, 1)
  expect_equal(kgc_raw$settings$n_factors, 1)
  expect_equal(kgc_cor_smc$settings$n_factors, 1)
  expect_equal(kgc_efa_ml$settings$n_factors, 1)

})

rm(kgc_cor, kgc_cor_smc, kgc_raw, kgc_efa_ml, x, y, z, dat_sing, cor_sing,
   cor_nposdef)
