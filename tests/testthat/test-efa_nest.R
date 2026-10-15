# NEST simulates reference eigenvalues at every candidate factor count, which makes
# the file-top fixtures ~10s. Skipped by default; opt in with
# `Sys.setenv(EFATOOLS_TEST_SLOW = "true")`. See helper-slow.R.
if (is_slow_test()) {
# seed the reference-data simulation so the retained-factor counts are reproducible
set.seed(42)
nest_cor <- efa_nest(test_models$baseline$cormat, N = 500)
nest_raw <- efa_nest(GRiPS_raw)
}  # is_slow_test()

test_that("output class and dimensions are correct", {
  skip_if_not_slow()
  expect_s3_class(nest_cor, "efa_retention")
  expect_length(nest_cor, 6)
  expect_s3_class(nest_raw, "efa_retention")
  expect_length(nest_raw, 6)
  expect_equal(.retention_record(nest_cor, "NEST")$plot_type, "eigen")
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

  # the reference series is padded with NA to the length of the empirical
  # eigenvalues; only the tested positions carry a value
  expect_length(ref_cor, length(eig_cor))
  expect_length(ref_raw, length(eig_raw))
  expect_equal(sum(!is.na(ref_cor)), nf_cor + 1)
  expect_equal(sum(!is.na(ref_raw)), nf_raw + 1)
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

test_that("errors are thrown correctly", {
  skip_if_not_slow()
  expect_error(efa_nest(1:5), class = "efa_input_not_matrix")
  expect_error(efa_nest(test_models$baseline$cormat), class = "efa_n_required")
  expect_message(efa_nest(GRiPS_raw), class = "efa_cor_from_data")
  expect_warning(efa_nest(GRiPS_raw, N = 20), class = "efa_n_from_data")
  expect_error(efa_nest(sing_raw), class = "efa_cor_singular")
  expect_error(efa_nest(sing_cor, N = sing_N), class = "efa_cor_singular")
})

test_that("NEST requires more observations than variables", {
  # the accepted call below simulates reference data, so pin the stream: the RNG
  # state reaching this test otherwise depends on which earlier blocks were skipped
  set.seed(321)
  R <- diag(5)
  expect_error(efa_nest(R, N = 0, n_datasets = 1), "Must be >= 1")
  expect_error(efa_nest(R, N = 1, n_datasets = 1), class = "efa_n_too_small")
  expect_error(efa_nest(R, N = 5, n_datasets = 1), class = "efa_n_too_small")
  expect_s3_class(efa_nest(R, N = 6, n_datasets = 1), "efa_retention")
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
  expect_error(suppressWarnings(efa_nest(R_hey, N = 200)), class = "efa_nest_heywood")
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
  nest_bound <- efa_nest(R_bound, N = 200, n_datasets = 300)
  rec_bound <- .retention_record(nest_bound, "NEST")

  expect_equal(nest_bound$n_factors[["NEST"]], 2)
  # every tested empirical eigenvalue exceeded its reference (no rejection), so
  # the reference series has one value per retained factor (not n_factors + 1),
  # padded with NA up to the number of empirical eigenvalues
  expect_length(rec_bound$reference, length(rec_bound$y))
  expect_equal(sum(!is.na(rec_bound$reference)),
               nest_bound$n_factors[["NEST"]])
  expect_true(all(rec_bound$y[seq_len(nest_bound$n_factors[["NEST"]])] >
                    rec_bound$reference[seq_len(nest_bound$n_factors[["NEST"]])]))

  # The reference fits stay over-identified: an uncapped floor(.8 * nvar) = 3
  # would fit a 2-factor (df < 0) reference for 4 variables; the df > 0 bound
  # prevents that, so no under-identified reference fit is triggered.
  expect_no_warning(
    {
      set.seed(123)
      efa_nest(R_bound, N = 200, n_datasets = 300)
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

# The reference correlation matrices are drawn by .simulate_cfm_eigen from their Wishart
# sampling distribution. These fast tests run unconditionally: the first two check the
# draw against stats::rWishart() and pin it, the third drives efa_nest() itself so the
# wiring from each reference model into the kernel is covered without the slow fixtures.
test_that("the reference draw is a Wishart draw on the stats::rWishart() stream", {
  # Under one seed, the nf-th reference eigenvalue must equal a recompute from
  # stats::rWishart() with N - 1 degrees of freedom and the reference-model covariance,
  # rescaled with stats::cov2cor(): the kernel uses the same Bartlett draw order. This
  # checks the null (identity) model (nf = 1) and a one-factor model (nf = 2). The
  # tolerance covers the different linear-algebra routines of the two computations.
  N <- 60; p <- 4; nd <- 4
  L <- matrix(c(.7, .6, .5, .4), p, 1)
  Psi <- as.vector(1 - L^2)

  wishart_eigen <- function(Sigma, nf) {
    vapply(seq_len(nd), function(i) {
      W <- stats::rWishart(1L, N - 1, Sigma)[, , 1]
      eigen(stats::cov2cor(W), symmetric = TRUE, only.values = TRUE)$values[nf]
    }, numeric(1))
  }

  set.seed(11)
  got_null <- .simulate_cfm_eigen(1L, N, matrix(numeric(0), p, 0), rep(1, p), nd)
  set.seed(11)
  expect_equal(as.vector(got_null), wishart_eigen(diag(p), 1L), tolerance = 1e-10)

  set.seed(12)
  got_one <- .simulate_cfm_eigen(2L, N, L, Psi, nd)
  set.seed(12)
  expect_equal(as.vector(got_one), wishart_eigen(tcrossprod(L) + diag(Psi), 2L),
               tolerance = 1e-10)
})

test_that("the factor-model reference draw matches its pinned reference eigenvalues", {
  # A one-factor reference model (the nf = 2 draw path). The reference eigenvalues are
  # pinned; the draw must reproduce them under this seed. expect_equal (not identical)
  # because the linear algebra downstream of the platform-independent RNG draw varies by
  # ~1e-13 across BLAS implementations, whereas any change to the draw stream is O(1).
  N <- 50; nd <- 5
  L <- matrix(c(.7, .6, .5, .4), 4, 1)
  Psi <- 1 - L^2

  set.seed(123)
  got <- .simulate_cfm_eigen(2L, N, L, as.vector(Psi), nd)

  expect_equal(
    as.vector(got),
    c(1.01450926035357, 0.972389824811252, 0.894813770054858,
      0.860823287614562, 0.938033057140992)
  )
})

test_that("NEST recovers a clear two-factor structure end to end", {
  # Drives efa_nest() itself so the wiring from each reference model into .simulate_cfm_eigen
  # -- the loadings and uniquenesses efa_nest() passes it -- is exercised by a fast test. Two
  # well-separated factors keep the retained count stable at a small number of simulated
  # datasets; the nf = 3 step builds a two-factor reference, covering a multi-column
  # loading matrix in the draw.
  R2 <- matrix(.1, 6, 6)
  R2[1:3, 1:3] <- .6
  R2[4:6, 4:6] <- .6
  diag(R2) <- 1

  set.seed(2024)
  nest_2f <- efa_nest(R2, N = 300, n_datasets = 100)

  expect_s3_class(nest_2f, "efa_retention")
  expect_equal(nest_2f$n_factors[["NEST"]], 2)
})

if (is_slow_test()) rm(nest_cor, nest_raw)
