
pa_cor <- PARALLEL(test_models$baseline$cormat, N = 500)
pa_cor_pca <- PARALLEL(test_models$baseline$cormat, N = 500, eigen_type = "PCA")
pa_raw <- PARALLEL(GRiPS_raw)
pa_nodat <- PARALLEL(N = 20, n_vars = 5)
pa_craw <- PARALLEL(test_models$baseline$cormat, N = 500, eigen_type = "PCA",
                    decision_rule = "crawford")
pa_perc <- PARALLEL(test_models$baseline$cormat, N = 500, eigen_type = "PCA",
                    decision_rule = "percentile")

test_that("output class and dimensions are correct", {
  expect_s3_class(pa_cor, "efa_retention")
  expect_length(pa_cor, 7)
  expect_s3_class(pa_raw, "efa_retention")
  expect_length(pa_raw, 7)
  expect_s3_class(pa_cor_pca, "efa_retention")
  expect_s3_class(pa_nodat, "efa_retention")
  expect_s3_class(pa_craw, "efa_retention")
  expect_s3_class(pa_perc, "efa_retention")

  expect_named(pa_cor$n_factors, c("PCA", "SMC", "EFA"))
  expect_named(pa_cor_pca$n_factors, "PCA")
  expect_named(pa_nodat$n_factors, c("PCA", "SMC", "EFA"))
  expect_equal(.retention_record(pa_cor, "PCA")$plot_type, "eigen")
})


test_that("found eigenvalues are correct", {
  # real eigenvalues form the solid line (record $y)
  expect_equal(sum(.retention_record(pa_cor, "PCA")$y),
               ncol(test_models$baseline$cormat))
  expect_lt(sum(.retention_record(pa_cor, "SMC")$y),
            ncol(test_models$baseline$cormat))
  expect_lt(sum(.retention_record(pa_cor, "EFA")$y),
            ncol(test_models$baseline$cormat))

  # simulated eigenvalues form the dashed reference series (record $references)
  expect_equal(sum(.retention_record(pa_cor, "PCA")$references$Means),
               ncol(test_models$baseline$cormat))
  expect_gt(sum(.retention_record(pa_cor, "PCA")$references[["95 Percentile"]]),
            ncol(test_models$baseline$cormat))
  expect_lt(sum(.retention_record(pa_cor, "SMC")$references$Means),
            ncol(test_models$baseline$cormat))
  expect_lt(sum(.retention_record(pa_cor, "EFA")$references$Means),
            ncol(test_models$baseline$cormat))

  expect_equal(sum(.retention_record(pa_raw, "PCA")$y), ncol(GRiPS_raw))
  expect_named(.retention_record(pa_raw, "PCA")$references,
               c("Means", "95 Percentile"))
  expect_equal(sum(.retention_record(pa_raw, "PCA")$references$Means),
               ncol(GRiPS_raw))
  expect_gt(sum(.retention_record(pa_raw, "PCA")$references[["95 Percentile"]]),
            ncol(GRiPS_raw))
  expect_lt(sum(.retention_record(pa_raw, "SMC")$y), ncol(GRiPS_raw))
  expect_lt(sum(.retention_record(pa_raw, "EFA")$y), ncol(GRiPS_raw))

  # references is a named two-series list (means + percentile)
  refs <- .retention_record(pa_cor, "PCA")$references
  checkmate::expect_list(refs, len = 2)
  expect_named(refs, c("Means", "95 Percentile"))

  # only the requested eigenvalue type produces a record
  expect_null(.retention_record(pa_cor_pca, "SMC"))
  expect_null(.retention_record(pa_cor_pca, "EFA"))

  # no real data: no real-eigenvalue series, but the simulated references remain
  expect_null(.retention_record(pa_nodat, "PCA")$y)
  expect_named(.retention_record(pa_nodat, "PCA")$references,
               c("Means", "95 Percentile"))
})


test_that("identified number of factors is correct", {
  expect_equal(pa_cor$n_factors[["PCA"]], 3)
  expect_equal(pa_cor$n_factors[["SMC"]], 3)
  expect_equal(pa_cor$n_factors[["EFA"]], 7, tolerance = 2)

  expect_equal(pa_raw$n_factors[["PCA"]], 1)
  expect_equal(pa_raw$n_factors[["SMC"]], 1)
  expect_equal(pa_raw$n_factors[["EFA"]], 3, tolerance = 2)

  expect_true(all(is.na(pa_nodat$n_factors)))

  expect_equal(pa_cor_pca$n_factors[["PCA"]], 3)
  expect_equal(pa_craw$n_factors[["PCA"]], 3)
  expect_equal(pa_perc$n_factors[["PCA"]], 3)
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
  expect_error(PARALLEL(1:5), class = "efa_input_not_matrix")
  expect_warning(suppressMessages(PARALLEL(GRiPS_raw, n_vars = 5)), class = "efa_nvars_from_data")
  expect_warning(suppressMessages(PARALLEL(GRiPS_raw, N = 20, eigen_type = "PCA")), class = "efa_n_from_data")
  expect_error(suppressMessages(PARALLEL(N = 500)), class = "efa_nvars_required")
  expect_message(PARALLEL(GRiPS_raw, eigen_type = "PCA"), class = "efa_cor_from_data")
  expect_error(PARALLEL(test_models$baseline$cormat, eigen_type = "PCA"), class = "efa_n_required")
  expect_warning(PARALLEL(test_models$baseline$cormat, N = 500,
                          eigen_type = "PCA", decision_rule = "crawford",
                          percent = 80), class = "efa_parallel_crawford")
  expect_error(PARALLEL(dat_sing), class = "efa_cor_singular")
  expect_error(PARALLEL(cor_sing, N = 10), class = "efa_cor_singular")
  expect_warning(PARALLEL(burt, N = 100, eigen_type = "PCA"), class = "efa_cor_smoothed")
  expect_error(PARALLEL(test_models$baseline$cormat, N = 15), class = "efa_n_too_small")
  expect_error(PARALLEL(test_models$baseline$cormat, N = 18), class = "efa_n_too_small")
})

test_that(".parallel_chunks splits exactly into non-negative integer chunks", {
  # the documented degenerate case: 11 datasets across 7 workers must not yield a
  # negative chunk (the old round-and-backfill produced c(2, 2, 2, 2, 2, 2, -1))
  chunks_11_7 <- .parallel_chunks(11, 7)
  expect_length(chunks_11_7, 7)
  expect_equal(sum(chunks_11_7), 11)
  expect_true(all(chunks_11_7 >= 0))

  # a range of awkward (n_datasets, n_cores) pairs: always exact, never negative
  for (n_cores in 1:8) {
    for (n_datasets in c(1, 3, 11, 100, 1000)) {
      chunks <- .parallel_chunks(n_datasets, n_cores)
      expect_length(chunks, n_cores)
      expect_equal(sum(chunks), n_datasets)
      expect_true(all(chunks >= 0))
    }
  }

  # a single worker (the sequential default) takes all datasets
  expect_equal(.parallel_chunks(1000, 1), 1000)
})


test_that(".parallel_summarise uses stats::quantile for the percentile series", {
  set.seed(42)
  eig_vals <- matrix(stats::rnorm(1000 * 4), nrow = 1000, ncol = 4)

  res <- .parallel_summarise(eig_vals, percent = 95, n_vars = 4)
  expect_equal(dim(res), c(4, 2))
  expect_equal(res[, 1], colMeans(eig_vals))
  for (root in 1:4) {
    expect_equal(res[root, 2],
                 stats::quantile(eig_vals[, root], probs = 0.95, names = FALSE))
  }

  # multiple percentiles produce one column each
  res2 <- .parallel_summarise(eig_vals, percent = c(50, 95), n_vars = 4)
  expect_equal(dim(res2), c(4, 3))
  for (root in 1:4) {
    expect_equal(res2[root, 2:3],
                 stats::quantile(eig_vals[, root], probs = c(0.5, 0.95),
                                 names = FALSE))
  }

  # a NA/NaN simulated eigenvalue is tolerated, not turned into a hard error
  eig_na <- eig_vals
  eig_na[1, 1] <- NA
  expect_no_error(res_na <- .parallel_summarise(eig_na, percent = 95, n_vars = 4))
  expect_equal(res_na[1, 2],
               stats::quantile(eig_na[, 1], probs = 0.95, names = FALSE,
                               na.rm = TRUE))
})

rm(pa_cor, pa_cor_pca, pa_raw, pa_nodat, pa_craw, pa_perc, x, y, z, dat_sing,
   cor_sing, burt)
