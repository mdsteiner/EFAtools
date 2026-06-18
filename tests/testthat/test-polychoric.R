# Polychoric / tetrachoric correlation matrix (.polychoric + the C++ backend).
#
# The two-step estimator is validated against the conditional-ML two-step references
# polycor::polychor(ML = FALSE) and psych::polychoric(correct = FALSE, global = FALSE);
# on complete data the pooled thresholds used here coincide with the per-pair thresholds
# those functions derive, so the estimates agree. The bivariate-normal rectangle that
# drives the likelihood is cross-checked against mnormt::sadmvn. Heavy reference loops
# over the large item sets run only under EFATOOLS_TEST_SLOW.

# Reference two-step polychoric matrix via polycor, pairwise.
.ref_polychor <- function(x) {
  p <- ncol(x)
  m <- diag(p)
  for (i in seq_len(p - 1L)) for (j in (i + 1L):p) {
    m[i, j] <- m[j, i] <- suppressWarnings(polycor::polychor(x[, i], x[, j], ML = FALSE))
  }
  m
}

# Deterministic raw ordinal data whose two-step polychoric matrix is not positive
# definite: four items are noisy monotone functions of the common score z1 + z2, so the
# columns are near-collinear and disattenuation pushes the smallest eigenvalue clearly
# negative (~ -0.07). Built without quantile cut points so it cannot tie / collapse a
# category across RNG streams.
.nonpd_ordinal <- function() {
  set.seed(1)
  n <- 120L
  z1 <- sample(1:5, n, TRUE)
  z2 <- sample(1:5, n, TRUE)
  core <- z1 + z2
  derive <- function() {
    v <- core + sample(c(-1L, 0L, 1L), n, TRUE)
    as.integer(pmin(6L, pmax(1L, v - 1L)))
  }
  cbind(z1, z2, derive(), derive(), derive(), derive())
}

test_that(".bvn_rect_cpp matches mnormt::sadmvn rectangle probabilities", {
  skip_on_cran()
  skip_if_not_installed("mnormt")

  cuts <- c(-Inf, -2.5, -1, 0, 1, 2.5, Inf)
  worst <- 0
  for (rho in c(0, 0.3, 0.5, 0.75, 0.9, 0.95, -0.8)) {
    V <- matrix(c(1, rho, rho, 1), 2L)
    for (ia in seq_len(length(cuts) - 1L)) for (ib in seq_len(length(cuts) - 1L)) {
      ours <- .bvn_rect_cpp(cuts[ia], cuts[ia + 1L], cuts[ib], cuts[ib + 1L], rho)
      ref <- as.numeric(mnormt::sadmvn(c(cuts[ia], cuts[ib]),
                                       c(cuts[ia + 1L], cuts[ib + 1L]),
                                       c(0, 0), V))
      worst <- max(worst, abs(ours - ref))
    }
  }
  expect_lt(worst, 1e-7)
})

test_that("polychoric matrix matches polycor and psych on GRiPS", {
  skip_on_cran()
  skip_if_not_installed("polycor")

  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
  ours <- .polychoric(g)$R

  expect_lt(max(abs(ours - .ref_polychor(g))), 1e-4)

  psych_rho <- suppressWarnings(
    psych::polychoric(g, correct = FALSE, global = FALSE)$rho)
  expect_lt(max(abs(ours - unname(psych_rho))), 1e-4)
})

test_that("tetrachoric (2 categories) matches polycor and psych", {
  skip_on_cran()
  skip_if_not_installed("polycor")

  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
  gb <- apply(g, 2L, function(col) as.integer(col > stats::median(col)))
  ours <- .polychoric(gb)$R

  expect_lt(max(abs(ours - .ref_polychor(gb))), 1e-4)

  tet_rho <- suppressWarnings(psych::tetrachoric(gb, correct = FALSE)$rho)
  expect_lt(max(abs(ours - unname(tet_rho))), 1e-4)
})

test_that("the matrix is a valid correlation matrix with named thresholds", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
  res <- .polychoric(g)

  expect_equal(dim(res$R), c(ncol(g), ncol(g)))
  expect_equal(dimnames(res$R), list(colnames(g), colnames(g)))
  expect_equal(diag(res$R), rep(1, ncol(g)), ignore_attr = TRUE)
  expect_true(isSymmetric(res$R))
  expect_true(all(res$R >= -1 & res$R <= 1))

  # one threshold vector per variable, length = (#categories - 1)
  expect_named(res$thresholds, colnames(g))
  n_cat <- apply(g, 2L, function(col) length(unique(col)))
  expect_equal(lengths(res$thresholds), n_cat - 1L, ignore_attr = TRUE)
})

test_that("the result is deterministic and independent of thread count", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
  expect_identical(.polychoric(g, n_threads = 1L)$R, .polychoric(g, n_threads = 1L)$R)
  expect_identical(.polychoric(g, n_threads = 1L)$R, .polychoric(g, n_threads = 4L)$R)
})

test_that("a constant column is rejected with a classed condition", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
  g[, 1L] <- 3L
  expect_error(.polychoric(g), class = "efa_cor_constant_col")
})

test_that("a non-positive-definite matrix is left alone unless nearest_pd is requested", {
  x <- .nonpd_ordinal()

  raw <- .polychoric(x)$R
  expect_lt(min(eigen(raw, symmetric = TRUE, only.values = TRUE)$values), 0)

  expect_warning(adj <- .polychoric(x, nearest_pd = TRUE)$R, class = "efa_cor_smoothed")
  expect_gt(min(eigen(adj, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_equal(diag(adj), rep(1, ncol(x)), ignore_attr = TRUE)
})

test_that("missing data is handled pairwise-complete", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
  gm <- g
  gm[1:20, 1] <- NA          # scattered missingness; every pair still overlaps
  gm[21:40, 2] <- NA
  res <- .polychoric(gm)
  expect_false(anyNA(res$R))
  expect_true(isSymmetric(res$R))
  expect_equal(diag(res$R), rep(1, ncol(gm)), ignore_attr = TRUE)
})

test_that("a pair with no overlapping complete cases is rejected with a classed condition", {
  v6 <- c(1L, 2L, 3L, 1L, 2L, 3L)
  na6 <- rep(NA_integer_, 6L)
  # columns 1 and 2 are never observed on the same row -> that pair is uncomputable;
  # column 3 overlaps both so only the 1-2 entry is NA.
  m <- cbind(c(v6, na6), c(na6, v6), rep(1:3, 4L))
  expect_error(.polychoric(m), class = "efa_cor_na")
})

test_that("the empty-cell continuity correction runs and changes the estimate", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
  res <- .polychoric(g, correct = 0.5)
  expect_false(anyNA(res$R))
  expect_true(isSymmetric(res$R))
  expect_equal(diag(res$R), rep(1, ncol(g)), ignore_attr = TRUE)
  # GRiPS has empty corner cells, so the pseudo-counts must move at least one estimate.
  expect_false(isTRUE(all.equal(res$R, .polychoric(g, correct = 0)$R)))
})

test_that("a likely-continuous (many-category) variable warns", {
  x <- cbind(rep(1:12, length.out = 100L), rep(1:3, length.out = 100L))
  expect_warning(.polychoric(x), class = "efa_cor_many_categories")
})

test_that("invalid correct or n_threads are rejected with a classed condition", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
  expect_error(.polychoric(g, correct = -1), class = "efa_cor_bad_arg")
  expect_error(.polychoric(g, n_threads = 0L), class = "efa_cor_bad_arg")
})

test_that("polychoric matches polycor on DOSPERT and UPPS", {
  skip_on_cran()
  skip_if_not_slow()
  skip_if_not_installed("polycor")

  for (data in list(DOSPERT_raw, UPPS_raw)) {
    x <- data[stats::complete.cases(data), ]
    ours <- .polychoric(x)$R
    expect_lt(max(abs(ours - .ref_polychor(x))), 1e-4)
    psych_rho <- suppressWarnings(
      psych::polychoric(x, correct = FALSE, global = FALSE)$rho)
    expect_lt(max(abs(ours - unname(psych_rho))), 1e-4)
  }
})
