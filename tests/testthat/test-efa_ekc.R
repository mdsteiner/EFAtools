
ekc_cor <- efa_ekc(test_models$baseline$cormat, N = 500)
ekc_raw <- efa_ekc(GRiPS_raw)
ekc_both <- efa_ekc(test_models$baseline$cormat, N = 500,
                type = c("BvA2017", "AM2019"))

test_that("output class and dimensions are correct", {
  expect_s3_class(ekc_cor, "efa_retention")
  expect_length(ekc_cor, 6)
  expect_s3_class(ekc_raw, "efa_retention")
  expect_length(ekc_raw, 6)

  expect_equal(unname(ekc_cor$criterion["id"]), "EKC")
  expect_length(ekc_cor$results, 1)
  expect_length(ekc_both$results, 2)
  expect_named(ekc_both$n_factors, c("BvA2017", "AM2019"))
  expect_equal(ekc_cor$results[[1]]$plot_type, "eigen")
})


test_that("found eigenvalues are correct", {
  expect_equal(sum(ekc_cor$results[[1]]$y),
               ncol(test_models$baseline$cormat))
  expect_equal(sum(ekc_raw$results[[1]]$y), ncol(GRiPS_raw))
  expect_length(ekc_cor$results[[1]]$y,
                ncol(test_models$baseline$cormat))
  expect_length(ekc_raw$results[[1]]$y, ncol(GRiPS_raw))
})

test_that("reference eigenvalues are correct", {
  expect_equal(ekc_cor$results[[1]]$reference[floor(ncol(test_models$baseline$cormat) / 2)], 1)
  expect_equal(ekc_raw$results[[1]]$reference[floor(ncol(GRiPS_raw) / 2)], 1)
  expect_length(ekc_cor$results[[1]]$reference,
                length(ekc_cor$results[[1]]$y))
  expect_length(ekc_raw$results[[1]]$reference, length(ekc_raw$results[[1]]$y))
})

test_that("identified number of factors is correct", {
  expect_equal(ekc_cor$n_factors[["BvA2017"]], 3)
  expect_equal(ekc_raw$n_factors[["BvA2017"]], 1)
})

test_that("the BvA2017 reference series reproduces the published equation", {
  # Braeken and van Assen (2017, p. 454) define the j-th reference eigenvalue as
  # max(((J - sum_{i < j} lambda_i) / (J - j + 1)) * (1 + sqrt(J / N))^2, 1),
  # transcribed here as an explicit loop -- an independent expression of the
  # equation the vectorised cumsum form in efa_ekc() implements. The equation holds
  # for any N, so one nominal value is used for both matrices.
  N <- 500
  for (cmat in list(test_models$baseline$cormat, stats::cor(GRiPS_raw))) {
    lambda <- eigen(cmat, symmetric = TRUE, only.values = TRUE)$values
    J <- ncol(cmat)

    ref <- numeric(J)
    for (j in seq_len(J)) {
      used <- if (j == 1) 0 else sum(lambda[seq_len(j - 1)])
      ref[j] <- max((J - used) / (J - j + 1) * (1 + sqrt(J / N))^2, 1)
    }

    out <- efa_ekc(cmat, N = N)
    expect_equal(.retention_record(out, "BvA2017")$reference, ref)
    # the paper's stopping rule: retain the leading run of eigenvalues above their
    # reference value
    expect_equal(out$n_factors[["BvA2017"]], sum(cumprod(lambda > ref)))
  }
})

test_that("AM2019 yields a well-defined (non-NA) number of factors", {
  # The no-crossing case is mapped to p (retain all). For any valid correlation
  # matrix the smallest eigenvalue is <= 1 <= the floored last reference, so a
  # crossing always exists and the AM2019 count is never NA.
  for (cmat in list(test_models$baseline$cormat, test_models$case_11b$cormat)) {
    nf <- efa_ekc(cmat, N = 500, type = "AM2019")$n_factors[["AM2019"]]
    expect_false(is.na(nf))
    expect_gte(nf, 0)
    expect_lte(nf, ncol(cmat))
  }

  nf_raw <- suppressMessages(efa_ekc(GRiPS_raw, type = "AM2019"))$n_factors[["AM2019"]]
  expect_false(is.na(nf_raw))
  expect_gte(nf_raw, 0)
  expect_lte(nf_raw, ncol(GRiPS_raw))
})

# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z, rnorm(10), rnorm(10), rnorm(10)), ncol = 6)
cor_sing <- stats::cor(dat_sing)

cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors are thrown correctly", {
  expect_error(efa_ekc(1:5), class = "efa_input_not_matrix")
  expect_error(efa_ekc(test_models$baseline$cormat), class = "efa_n_required")
  expect_message(efa_ekc(GRiPS_raw), class = "efa_cor_from_data")
  expect_warning(efa_ekc(GRiPS_raw, N = 20), class = "efa_n_from_data")
  expect_error(efa_ekc(dat_sing), class = "efa_cor_singular")
  expect_error(efa_ekc(cor_sing, N = 20), class = "efa_cor_singular")
  expect_warning(efa_ekc(cor_nposdef, N = 20), class = "efa_cor_smoothed")
})

test_that("settings are returned correctly", {
  expect_named(ekc_cor$settings, c("use", "cor_method", "N", "type"))
  expect_named(ekc_raw$settings, c("use", "cor_method", "N", "type"))

  expect_equal(ekc_cor$settings$N, 500)
  expect_equal(ekc_raw$settings$N, 810)

  expect_equal(ekc_cor$settings$use, "pairwise.complete.obs")
  expect_equal(ekc_raw$settings$use, "pairwise.complete.obs")

  expect_equal(ekc_cor$settings$cor_method, "pearson")
  expect_equal(ekc_raw$settings$cor_method, "pearson")
})

rm(ekc_cor, ekc_raw, ekc_both, x, y, z, dat_sing, cor_sing, cor_nposdef)
