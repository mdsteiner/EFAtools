# NEST simulates reference eigenvalues at every candidate factor count, which makes
# the file-top fixtures ~10s. Skipped by default; opt in with
# `Sys.setenv(EFATOOLS_TEST_SLOW = "true")`. See helper-slow.R.
if (is_slow_test()) {
# seed the reference-data simulation so the retained-factor counts are reproducible
set.seed(42)
nest_cor <- NEST(test_models$baseline$cormat, N = 500)
nest_raw <- NEST(GRiPS_raw)
}  # is_slow_test()

test_that("output class and dimensions are correct", {
  skip_if_not_slow()
  expect_s3_class(nest_cor, "efa_retention")
  expect_length(nest_cor, 7)
  expect_s3_class(nest_raw, "efa_retention")
  expect_length(nest_raw, 7)
  expect_equal(.retention_record(nest_cor, "NEST")$plot_type, "none")
})


test_that("found eigenvalues are correct", {
  skip_if_not_slow()
  expect_equal(sum(.retention_record(nest_cor, "NEST")$y),
               ncol(test_models$baseline$cormat))
  expect_equal(sum(.retention_record(nest_raw, "NEST")$y), ncol(GRiPS_raw))
  expect_length(.retention_record(nest_cor, "NEST")$y,
                ncol(test_models$baseline$cormat))
  expect_length(.retention_record(nest_raw, "NEST")$y, ncol(GRiPS_raw))
})

test_that("reference eigenvalues are correct", {
  skip_if_not_slow()
  nf_cor <- nest_cor$n_factors[["NEST"]]
  eig_cor <- .retention_record(nest_cor, "NEST")$y
  ref_cor <- .retention_record(nest_cor, "NEST")$reference
  nf_raw <- nest_raw$n_factors[["NEST"]]
  eig_raw <- .retention_record(nest_raw, "NEST")$y
  ref_raw <- .retention_record(nest_raw, "NEST")$reference

  expect_length(ref_cor, nf_cor + 1)
  expect_length(ref_raw, nf_raw + 1)
  expect_lte(eig_cor[nf_cor + 1], ref_cor[nf_cor + 1])
  expect_true(all(eig_cor[1:nf_cor] > ref_cor[1:nf_cor]))
  expect_lte(eig_raw[nf_raw + 1], ref_raw[nf_raw + 1])
  expect_true(all(eig_raw[1:nf_raw] > ref_raw[1:nf_raw]))
})

test_that("identified number of factors is correct", {
  skip_if_not_slow()
  expect_equal(nest_cor$n_factors[["NEST"]], 3)
  expect_equal(nest_raw$n_factors[["NEST"]], 1)
})

# Create singular correlation matrix for tests
set.seed(7)
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z, rnorm(10), rnorm(10), rnorm(10)), ncol = 6)
cor_sing <- stats::cor(dat_sing)

cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors are thrown correctly", {
  skip_if_not_slow()
  expect_error(NEST(1:5), class = "efa_input_not_matrix")
  expect_error(NEST(test_models$baseline$cormat), class = "efa_n_required")
  expect_message(NEST(GRiPS_raw), class = "efa_cor_from_data")
  expect_warning(NEST(GRiPS_raw, N = 20), class = "efa_n_from_data")
  expect_error(NEST(dat_sing), class = "efa_cor_singular")
  expect_error(NEST(cor_sing, N = 20), class = "efa_cor_singular")
})

test_that("a Heywood case in the reference model raises a classed error", {
  skip_if_not_slow()
  # a positive-definite matrix whose one-factor model has a communality above 1:
  # the reference data cannot be simulated (negative uniqueness), so NEST must
  # abort with a classed condition instead of an unclassed C++ error
  R_hey <- diag(4)
  R_hey[1, 2:4] <- R_hey[2:4, 1] <- .9
  R_hey[2, 3:4] <- R_hey[3:4, 2] <- .78
  R_hey[3, 4] <- R_hey[4, 3] <- .78
  expect_gt(min(eigen(R_hey, symmetric = TRUE, only.values = TRUE)$values), 0)
  # the reference EFA legitimately warns about the Heywood case; assert the abort
  expect_error(suppressWarnings(NEST(R_hey, N = 200)), class = "efa_nest_heywood")
})

test_that("the last accepted model is retained at the no-stop boundary", {
  skip_if_not_slow()
  # Four variables with two clear factors: NEST accepts the 1st empirical
  # eigenvalue (vs the identity reference) and the 2nd (vs the one-factor
  # reference). The search then runs out of testable factors -- floor(.8 * nvar)
  # is 3, but the (nf - 1)-factor reference model must stay over-identified
  # (df > 0), which caps it at 2. With no rejection, the retained count must be
  # the last *accepted* model (2), not the loop's just-past index (1).
  R_bound <- matrix(c(1, .7, .1, .1,
                      .7, 1, .1, .1,
                      .1, .1, 1, .7,
                      .1, .1, .7, 1), 4, 4)

  set.seed(123)
  nest_bound <- NEST(R_bound, N = 200, n_datasets = 300)
  rec_bound <- .retention_record(nest_bound, "NEST")

  expect_equal(nest_bound$n_factors[["NEST"]], 2)
  # every tested empirical eigenvalue exceeded its reference (no rejection), so
  # the reference series has one entry per retained factor (not n_factors + 1)
  expect_length(rec_bound$reference, nest_bound$n_factors[["NEST"]])
  expect_true(all(rec_bound$y[seq_len(nest_bound$n_factors[["NEST"]])] >
                    rec_bound$reference[seq_len(nest_bound$n_factors[["NEST"]])]))

  # The reference fits stay over-identified: an uncapped floor(.8 * nvar) = 3
  # would fit a 2-factor (df < 0) reference for 4 variables; the df > 0 bound
  # prevents that, so no under-identified reference fit is triggered.
  expect_no_warning(
    {
      set.seed(123)
      NEST(R_bound, N = 200, n_datasets = 300)
    },
    class = "efa_underidentified"
  )
})

test_that("settings are returned correctly", {
  skip_if_not_slow()

  expect_equal(nest_cor$settings$N, 500)
  expect_equal(nest_raw$settings$N, 810)

  expect_equal(nest_cor$settings$use, "pairwise.complete.obs")
  expect_equal(nest_raw$settings$use, "pairwise.complete.obs")

  expect_equal(nest_cor$settings$cor_method, "pearson")
  expect_equal(nest_raw$settings$cor_method, "pearson")
})

if (is_slow_test()) rm(nest_cor, nest_raw)
rm(x, y, z, dat_sing, cor_sing, cor_nposdef)
