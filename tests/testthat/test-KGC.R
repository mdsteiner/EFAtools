kgc_cor <- KGC(test_models$baseline$cormat)
kgc_cor_smc <- KGC(test_models$baseline$cormat, eigen_type = "SMC")
kgc_raw <- KGC(GRiPS_raw)
# Check with an argument passed to "EFA"
kgc_efa_ml<- KGC(test_models$baseline$cormat, eigen_type = "EFA", method = "ML")

test_that("output class and dimensions are correct", {
  expect_s3_class(kgc_cor, "efa_retention")
  expect_length(kgc_cor, 7)
  expect_s3_class(kgc_cor_smc, "efa_retention")
  expect_length(kgc_cor_smc, 7)
  expect_s3_class(kgc_raw, "efa_retention")
  expect_s3_class(kgc_efa_ml, "efa_retention")

  expect_named(kgc_cor$n_factors, c("PCA", "SMC", "EFA"))
  expect_named(kgc_cor_smc$n_factors, "SMC")
  expect_named(kgc_efa_ml$n_factors, "EFA")
  expect_equal(.retention_record(kgc_cor, "PCA")$plot_type, "eigen")
  expect_equal(.retention_record(kgc_cor, "PCA")$threshold, 1)
})

test_that("found eigenvalues are correct", {
  expect_equal(sum(.retention_record(kgc_cor, "PCA")$y),
               ncol(test_models$baseline$cormat))
  expect_lt(sum(.retention_record(kgc_cor, "SMC")$y),
            ncol(test_models$baseline$cormat))
  expect_lt(sum(.retention_record(kgc_cor, "EFA")$y),
            ncol(test_models$baseline$cormat))

  expect_equal(sum(.retention_record(kgc_raw, "PCA")$y), ncol(GRiPS_raw))
  expect_lt(sum(.retention_record(kgc_raw, "SMC")$y), ncol(GRiPS_raw))
  expect_lt(sum(.retention_record(kgc_raw, "EFA")$y), ncol(GRiPS_raw))

  # Only the requested eigenvalue type produces a record
  expect_null(.retention_record(kgc_cor_smc, "PCA"))
  expect_null(.retention_record(kgc_cor_smc, "EFA"))
  expect_null(.retention_record(kgc_efa_ml, "PCA"))
  expect_null(.retention_record(kgc_efa_ml, "SMC"))
})

test_that("identified number of factors is correct", {
  expect_equal(kgc_cor$n_factors[["PCA"]], 3)
  expect_equal(kgc_cor$n_factors[["SMC"]], 1)
  expect_equal(kgc_cor$n_factors[["EFA"]], 1)

  expect_equal(kgc_raw$n_factors[["PCA"]], 1)
  expect_equal(kgc_raw$n_factors[["SMC"]], 1)
  expect_equal(kgc_raw$n_factors[["EFA"]], 1)

  expect_equal(kgc_cor_smc$n_factors[["SMC"]], 1)
  expect_false(any(c("PCA", "EFA") %in% names(kgc_cor_smc$n_factors)))

  expect_equal(kgc_efa_ml$n_factors[["EFA"]], 1)
  expect_false(any(c("PCA", "SMC") %in% names(kgc_efa_ml$n_factors)))
})

# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors are thrown correctly", {
  expect_error(KGC(1:5), class = "efa_input_not_matrix")
  expect_message(KGC(GRiPS_raw, eigen_type = "PCA"), class = "efa_cor_from_data")
  expect_error(KGC(dat_sing), class = "efa_cor_singular")
  expect_error(KGC(cor_sing, N = 10), class = "efa_cor_singular")
  expect_warning(KGC(cor_nposdef, N = 10), class = "efa_cor_smoothed")
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
