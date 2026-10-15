kgc_cor <- efa_kgc(test_models$baseline$cormat)
kgc_cor_smc <- efa_kgc(test_models$baseline$cormat, eigen_type = "SMC")
kgc_raw <- efa_kgc(GRiPS_raw)
# Check with an argument passed to "EFA"
kgc_efa_ml<- efa_kgc(test_models$baseline$cormat, eigen_type = "EFA", estimator = "ML")

test_that("output class and dimensions are correct", {
  expect_s3_class(kgc_cor, "efa_retention")
  expect_length(kgc_cor, 6)
  expect_s3_class(kgc_cor_smc, "efa_retention")
  expect_length(kgc_cor_smc, 6)
  expect_s3_class(kgc_raw, "efa_retention")
  expect_s3_class(kgc_efa_ml, "efa_retention")

  expect_named(kgc_cor$n_factors, c("PCA", "SMC", "EFA"))
  expect_named(kgc_cor_smc$n_factors, "SMC")
  expect_named(kgc_efa_ml$n_factors, "EFA")
  expect_equal(.retention_record(kgc_cor, "PCA")$plot_type, "eigen")
  # the plotted cutoff of each eigenvalue type
  expect_equal(vapply(kgc_cor$results, function(r) r$threshold, numeric(1)), c(1, 0, 0))
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
  # SMC and EFA count the positive eigenvalues of the reduced matrix. The eigenvalue
  # closest to 0 is at least .004 away from it in every fixture, far above any BLAS noise.
  expect_equal(kgc_cor$n_factors[["PCA"]], 3)
  expect_equal(kgc_cor$n_factors[["SMC"]], 8)
  expect_equal(kgc_cor$n_factors[["EFA"]], 7)

  expect_equal(kgc_raw$n_factors[["PCA"]], 1)
  expect_equal(kgc_raw$n_factors[["SMC"]], 3)
  expect_equal(kgc_raw$n_factors[["EFA"]], 4)

  expect_equal(kgc_cor_smc$n_factors[["SMC"]], 8)
  expect_false(any(c("PCA", "EFA") %in% names(kgc_cor_smc$n_factors)))

  expect_equal(kgc_efa_ml$n_factors[["EFA"]], 7)
  expect_false(any(c("PCA", "SMC") %in% names(kgc_efa_ml$n_factors)))
})

test_that("PCA counts the eigenvalues above 1, SMC and EFA those above 0", {
  # Six variables that all correlate .1 have the eigenvalues 1.5 and .9 (five times). The
  # SMC-reduced matrix subtracts 1 - SMC = .9 * 1.5 / 1.4 from each, which leaves 15 / 28
  # and -.9 / 14 (five times). So PCA retains 1 (the cutoff 0 would give 6), and SMC
  # retains 1 (the cutoff 1 would give 0).
  R <- matrix(0.1, 6, 6)
  diag(R) <- 1
  kgc <- efa_kgc(R, eigen_type = c("PCA", "SMC"))
  expect_equal(.retention_record(kgc, "SMC")$y, c(15 / 28, rep(-0.9 / 14, 5)),
               tolerance = 1e-10)
  expect_equal(kgc$n_factors, c(PCA = 1, SMC = 1))

  # the EFA eigenvalues follow the cutoff 0 too
  rec <- .retention_record(kgc_cor, "EFA")
  expect_equal(kgc_cor$n_factors[["EFA"]], sum(rec$y > 0))
  expect_gt(kgc_cor$n_factors[["EFA"]], sum(rec$y > 1))

  # An eigenvalue at the cutoff is not counted, whatever round-off the eigensolver adds.
  # One correlated pair among six variables gives the eigenvalues 1.5, .5, and 1 (four
  # times), and the SMC-reduced eigenvalues .75, -.25, and 0 (four times).
  R2 <- diag(6)
  R2[1, 2] <- R2[2, 1] <- 0.5
  expect_equal(efa_kgc(R2, eigen_type = c("PCA", "SMC"))$n_factors, c(PCA = 1, SMC = 1))
})

test_that("the PCA suggestion matches nFactors::nScree", {
  skip_if_not_installed("nFactors")
  # nScree's nkaiser is the eigenvalues-greater-than-one count on the unaltered
  # correlation matrix, i.e. the PCA variant of the Kaiser-Guttman criterion.
  for (cmat in list(test_models$baseline$cormat, stats::cor(GRiPS_raw))) {
    expect_equal(
      efa_kgc(cmat, eigen_type = "PCA")$n_factors[["PCA"]],
      nFactors::nScree(cmat, model = "components")$Components$nkaiser
    )
  }
})

test_that("errors are thrown correctly", {
  expect_error(efa_kgc(1:5), class = "efa_input_not_matrix")
  expect_message(efa_kgc(GRiPS_raw, eigen_type = "PCA"), class = "efa_cor_from_data")
  expect_error(efa_kgc(sing_raw), class = "efa_cor_singular")
  expect_error(efa_kgc(sing_cor, N = sing_N), class = "efa_cor_singular")
  expect_warning(efa_kgc(cor_nposdef, N = 10), class = "efa_cor_smoothed")
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

rm(kgc_cor, kgc_cor_smc, kgc_raw, kgc_efa_ml)
