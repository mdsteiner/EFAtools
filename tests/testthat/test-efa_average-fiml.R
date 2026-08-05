# cor_method = "fiml" propagated to the sibling entry points: efa_average() averages a
# two-stage FIML correlation (estimated once and reused across the grid), and the
# simulation-/eigenvalue-based retention and suitability functions leave "fiml" out of
# their choices. The EM engine and the EFA() routing themselves are covered elsewhere
# (test-fiml-moments.R, test-EFA-fiml.R).

# A missing-at-random fixture: a clean two-factor population, with column 1 fully observed
# and driving the missingness in the others (so the mechanism depends only on observed data).
# The deviates are coloured by the Cholesky factor rather than by an eigendecomposition: this
# population covariance has repeated eigenvalues, whose eigenvector basis is undetermined and
# hence settled by rounding, while chol() is unique for a positive definite matrix, so the same
# seed reproduces the same sample on every LAPACK build.
avg_fiml_mar_data <- function(n = 800, seed = 4242) {
  set.seed(seed)
  L <- matrix(0, 6, 2)
  L[1:3, 1] <- 0.7
  L[4:6, 2] <- 0.7
  S <- tcrossprod(L); diag(S) <- 1
  X <- matrix(stats::rnorm(n * 6), n) %*% chol(S)
  colnames(X) <- paste0("V", seq_len(6))
  X[X[, 1] >  0.8, 2] <- NA
  X[X[, 1] < -0.8, 3] <- NA
  X[X[, 1] >  1.2, 4] <- NA
  X[X[, 1] < -1.2, 5] <- NA
  X
}

test_that("efa_average with cor_method = 'fiml' analyses the two-stage FIML correlation", {
  skip_on_cran()

  X <- avg_fiml_mar_data()

  avg <- suppressMessages(suppressWarnings(
    efa_average(X, n_factors = 2, estimator = c("ML", "ULS"), rotation = "none",
                cor_method = "fiml", show_progress = FALSE)
  ))
  expect_s3_class(avg, "efa_average")

  # The analysed matrix is the two-stage FIML correlation: the EM saturated covariance,
  # standardized, estimated once up front from the raw data with missing values.
  em <- .fiml_em_moments(X)
  R_fiml <- stats::cov2cor(em$sigma)
  expect_equal(unclass(avg$orig_R), R_fiml, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(avg$settings$N, em$n)

  # Routing raw data through cor_method = "fiml" is consistent with passing the precomputed
  # FIML correlation matrix through the same grid: both analyse the same matrix with the same
  # (deterministic, unrotated ML/ULS) solutions, so the averaged loadings coincide. (That the
  # EM is run only once for the whole grid is pinned separately below.)
  avg_ref <- suppressMessages(suppressWarnings(
    efa_average(R_fiml, N = em$n, n_factors = 2, estimator = c("ML", "ULS"),
                rotation = "none", show_progress = FALSE)
  ))
  expect_equal(unclass(avg$loadings$average), unclass(avg_ref$loadings$average),
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("efa_average runs the FIML EM once for the whole grid, not per solution", {
  skip_on_cran()

  X <- avg_fiml_mar_data()

  # Count the EM calls during one efa_average run. The saturated moments are estimated once up
  # front in .prepare_cor_input() and the resulting correlation is reused by every solution in
  # the grid, so the engine is hit exactly once; a regression to per-cell EM would call it once
  # per grid cell. The per-cell EFA() calls receive that correlation matrix, so they never
  # re-enter the EM regardless of the future plan.
  real_em <- .fiml_em_moments
  n_em <- 0L
  testthat::local_mocked_bindings(
    .fiml_em_moments = function(...) {
      n_em <<- n_em + 1L
      real_em(...)
    }
  )

  invisible(suppressMessages(suppressWarnings(
    efa_average(X, n_factors = 2, estimator = c("ML", "ULS"), rotation = "none",
                cor_method = "fiml", show_progress = FALSE)
  )))

  expect_identical(n_em, 1L)
})

test_that("retention and suitability functions reject cor_method = 'fiml'", {
  # FIML is a single-fit missing-data correlation; the simulation-/eigenvalue-based
  # retention and suitability criteria do not support it, so it is left out of their
  # cor_method choices and rejected up front (before any computation) by the shared
  # choice-matching layer. Asserting that layer's condition class, rather than merely that
  # something failed, is what keeps the test tied to the rejection it is named for: any other
  # abort would otherwise satisfy it silently.
  expect_error(N_FACTORS(GRiPS_raw, cor_method = "fiml"), class = "efa_bad_choice")
  expect_error(KMO(GRiPS_raw, cor_method = "fiml"), class = "efa_bad_choice")
})
