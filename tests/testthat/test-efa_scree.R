scree_cor <- efa_scree(test_models$baseline$cormat)
scree_cor_smc <- efa_scree(test_models$baseline$cormat, eigen_type = "SMC")
scree_raw <- efa_scree(GRiPS_raw)
# Check with an argument passed to "EFA"
scree_efa_ml<- efa_scree(test_models$baseline$cormat, eigen_type = "EFA", estimator = "ML")

test_that("output class and dimensions are correct", {
  expect_s3_class(scree_cor, "efa_retention")
  expect_length(scree_cor, 6)
  expect_s3_class(scree_cor_smc, "efa_retention")
  expect_length(scree_cor_smc, 6)
  expect_s3_class(scree_raw, "efa_retention")
  expect_s3_class(scree_efa_ml, "efa_retention")

  expect_named(scree_cor$n_factors, c("PCA", "SMC", "EFA"))
  expect_named(scree_cor_smc$n_factors, "SMC")
  expect_named(scree_efa_ml$n_factors, "EFA")
  # the scree plot is purely visual: no numeric suggestion
  expect_true(all(is.na(scree_cor$n_factors)))
  expect_equal(.retention_record(scree_cor, "PCA")$plot_type, "eigen")
})

test_that("found eigenvalues are correct", {
  expect_equal(sum(.retention_record(scree_cor, "PCA")$y),
               ncol(test_models$baseline$cormat))
  expect_lt(sum(.retention_record(scree_cor, "SMC")$y),
            ncol(test_models$baseline$cormat))
  expect_lt(sum(.retention_record(scree_cor, "EFA")$y),
            ncol(test_models$baseline$cormat))

  expect_equal(sum(.retention_record(scree_raw, "PCA")$y), ncol(GRiPS_raw))
  expect_lt(sum(.retention_record(scree_raw, "SMC")$y), ncol(GRiPS_raw))
  expect_lt(sum(.retention_record(scree_raw, "EFA")$y), ncol(GRiPS_raw))

  # Only the requested eigenvalue type produces a record
  expect_null(.retention_record(scree_cor_smc, "PCA"))
  expect_null(.retention_record(scree_cor_smc, "EFA"))
  expect_null(.retention_record(scree_efa_ml, "PCA"))
  expect_null(.retention_record(scree_efa_ml, "SMC"))
})


test_that("errors are thrown correctly", {
  expect_error(efa_scree(1:5), class = "efa_input_not_matrix")
  expect_message(efa_scree(GRiPS_raw, eigen_type = "PCA"), class = "efa_cor_from_data")
  expect_error(efa_scree(sing_raw), class = "efa_cor_singular")
  expect_error(efa_scree(sing_cor, N = sing_N), class = "efa_cor_singular")
  expect_warning(efa_scree(cor_nposdef, N = 10), class = "efa_cor_smoothed")
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

rm(scree_cor, scree_cor_smc, scree_raw, scree_efa_ml)
