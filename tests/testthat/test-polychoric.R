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

# Listwise-complete GRiPS and its default polychoric matrix, computed once and
# reused by the structural blocks that would otherwise rebuild them. Tests that
# need non-default arguments (`correct`, `nearest_pd`, `n_threads`, `acov`) or
# that mutate the data still build their own.
g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
poly_g <- .polychoric(g)

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
  # The 12-node Gauss-Legendre rule reproduces the reference bivariate-normal integrator to
  # better than 1e-6 over this grid (worst case near rho = 0.95); this is far tighter than the
  # ~1e-4 agreement the polychoric matrix itself needs.
  expect_lt(worst, 1e-6)
})

test_that("polychoric matrix matches polycor and psych on GRiPS", {
  skip_on_cran()
  skip_if_not_installed("polycor")

  ours <- poly_g$R

  expect_lt(max(abs(ours - .ref_polychor(g))), 1e-4)

  psych_rho <- suppressWarnings(
    psych::polychoric(g, correct = FALSE, global = FALSE)$rho)
  expect_lt(max(abs(ours - unname(psych_rho))), 1e-4)
})

test_that("tetrachoric (2 categories) matches polycor and psych", {
  skip_on_cran()
  skip_if_not_installed("polycor")

  gb <- apply(g, 2L, function(col) as.integer(col > stats::median(col)))
  ours <- .polychoric(gb)$R

  expect_lt(max(abs(ours - .ref_polychor(gb))), 1e-4)

  tet_rho <- suppressWarnings(psych::tetrachoric(gb, correct = FALSE)$rho)
  expect_lt(max(abs(ours - unname(tet_rho))), 1e-4)
})

test_that("the matrix is a valid correlation matrix with named thresholds", {
  res <- poly_g

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
  expect_identical(.polychoric(g, n_threads = 1L)$R, .polychoric(g, n_threads = 1L)$R)
  expect_identical(.polychoric(g, n_threads = 1L)$R, .polychoric(g, n_threads = 4L)$R)
})

test_that("a constant column is rejected with a classed condition", {
  # local copy: this block mutates the column to force the constant-column path
  gx <- g
  gx[, 1L] <- 3L
  expect_error(.polychoric(gx), class = "efa_cor_constant_col")
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

# --- Asymptotic covariance of the polychoric correlations (acov = "diag" / "full") -------
#
# The two-step ACOV accounts for the estimated thresholds (Muthen, 1984; Joreskog, 1994);
# it is validated against lavaan's NACOV and an independent Monte-Carlo simulation. It is
# NOT checked against polycor::polychor(std.err = TRUE): that returns the threshold-FIXED
# conditional variance 1/(N * I_rho), which omits threshold-estimation uncertainty and runs
# ~7-23% below the true variance, so it is the wrong reference for the threshold-corrected
# quantity DWLS and robust SEs require (and that diag = diag(full) mandates).

# lavaan's variance-scale NACOV for the off-diagonal correlations, labelled "Vi-Vj" to match
# .polychoric()'s pair order. lavaan reports the per-observation Gamma, so divide by N.
.lav_acov <- function(x) {
  df <- as.data.frame(lapply(as.data.frame(x), ordered))
  fit <- lavaan::lavCor(df, ordered = names(df), output = "fit", se = "robust.sem")
  G <- lavaan::lavInspect(fit, "gamma")
  ci <- grep("~~", rownames(G))
  out <- G[ci, ci, drop = FALSE] / nrow(x)
  lab <- gsub("~~", "-", rownames(out))
  dimnames(out) <- list(lab, lab)
  out
}

# Independent reference for one rho's asymptotic variance: the same two-step (Muthen/Joreskog)
# outer-product sandwich, but the bread/meat assembly (A11/A21/A22/IF) is rebuilt from scratch
# in plain R, and the rho-derivative uses the bivariate-normal density at the four corners
# (Plackett, 1954) -- a different formula than the backend's conditioning integral. This
# catches an error in the C++ influence-function assembly that diag == diag(full) cannot. The
# cell probabilities (via .bvn_rect_cpp) and the thresholds (qnorm of the cumulative
# proportions) are deliberately the SAME as the backend, so this guards the assembly while the
# lavaan cross-checks independently guard the cell-probability quadrature. Inputs are 0-based
# category codes for one complete-data pair.
.if_diag_pair <- function(a, b) {
  N <- length(a)
  Ki <- max(a) + 1L; Kj <- max(b) + 1L
  ti <- stats::qnorm(cumsum(tabulate(a + 1L, Ki)) / N)[-Ki]
  tj <- stats::qnorm(cumsum(tabulate(b + 1L, Kj)) / N)[-Kj]
  rcut <- c(-Inf, ti, Inf); ccut <- c(-Inf, tj, Inf)
  rho <- .polychoric(cbind(a, b))$R[1, 2]; s <- sqrt(1 - rho^2)
  dbn <- function(xv, yv) if (is.finite(xv) && is.finite(yv)) {
    exp(-(xv^2 - 2 * rho * xv * yv + yv^2) / (2 * (1 - rho^2))) / (2 * pi * s)
  } else 0
  P <- dP <- matrix(0, Ki, Kj)
  for (aa in seq_len(Ki)) for (bb in seq_len(Kj)) {
    P[aa, bb]  <- .bvn_rect_cpp(rcut[aa], rcut[aa + 1L], ccut[bb], ccut[bb + 1L], rho)
    dP[aa, bb] <- dbn(rcut[aa + 1L], ccut[bb + 1L]) - dbn(rcut[aa], ccut[bb + 1L]) -
                  dbn(rcut[aa + 1L], ccut[bb]) + dbn(rcut[aa], ccut[bb])
  }
  Pf <- pmax(P, .Machine$double.xmin); dxr <- dP / Pf
  n_ab <- matrix(0, Ki, Kj)
  for (r in seq_len(N)) n_ab[a[r] + 1L, b[r] + 1L] <- n_ab[a[r] + 1L, b[r] + 1L] + 1
  A22 <- sum(n_ab * dxr^2)
  A21i <- vapply(seq_len(Ki - 1L), function(k) {
    Bc <- stats::pnorm((ccut[-1L] - rho * ti[k]) / s) -
          stats::pnorm((ccut[-(Kj + 1L)] - rho * ti[k]) / s)
    stats::dnorm(ti[k]) * sum(Bc * (n_ab[k, ] * dxr[k, ] / Pf[k, ] -
                                    n_ab[k + 1L, ] * dxr[k + 1L, ] / Pf[k + 1L, ]))
  }, numeric(1))
  A21j <- vapply(seq_len(Kj - 1L), function(k) {
    Bc <- stats::pnorm((rcut[-1L] - rho * tj[k]) / s) -
          stats::pnorm((rcut[-(Ki + 1L)] - rho * tj[k]) / s)
    stats::dnorm(tj[k]) * sum(Bc * (n_ab[, k] * dxr[, k] / Pf[, k] -
                                    n_ab[, k + 1L] * dxr[, k + 1L] / Pf[, k + 1L]))
  }, numeric(1))
  bread <- function(K, tau, m) {                     # IFth = S (A11)^{-1}, A11 = S' diag(m) S
    S <- matrix(0, K, K - 1L)
    for (aa in 0:(K - 1L)) {
      ipm <- if (m[aa + 1L] > 0) N / m[aa + 1L] else 0
      if (aa <= K - 2L) S[aa + 1L, aa + 1L] <- S[aa + 1L, aa + 1L] + stats::dnorm(tau[aa + 1L]) * ipm
      if (aa >= 1L)     S[aa + 1L, aa]       <- S[aa + 1L, aa]       - stats::dnorm(tau[aa]) * ipm
    }
    S %*% solve(t(S) %*% diag(m) %*% S)
  }
  Ti <- as.numeric(bread(Ki, ti, tabulate(a + 1L, Ki)) %*% A21i)
  Tj <- as.numeric(bread(Kj, tj, tabulate(b + 1L, Kj)) %*% A21j)
  IF <- (dxr - outer(Ti, rep(1, Kj)) - outer(rep(1, Ki), Tj)) / A22
  sum(n_ab * IF^2)
}

test_that("acov diag/full are well-formed and mutually consistent", {
  x <- DOSPERT_raw[stats::complete.cases(DOSPERT_raw), 1:6]
  p <- ncol(x)
  pstar <- p * (p - 1L) / 2L
  labels <- apply(utils::combn(colnames(x), 2L), 2L, paste, collapse = "-")

  d <- .polychoric(x, acov = "diag")
  f <- .polychoric(x, acov = "full")

  # shapes + labelling
  expect_length(d$acov, pstar)
  expect_named(d$acov, labels)
  expect_equal(dim(f$acov), c(pstar, pstar))
  expect_equal(dimnames(f$acov), list(labels, labels))

  # the diagonal is a valid set of variances; the cheap diag path equals diag(full) (a
  # structural check that the per-cell influence is scattered/cross-multiplied correctly)
  expect_true(all(is.finite(d$acov)))
  expect_true(all(d$acov > 0))
  expect_equal(unname(d$acov), unname(diag(f$acov)), tolerance = 1e-10)

  # full Gamma is symmetric and positive-semidefinite (it is a Gram matrix by construction)
  expect_true(isSymmetric(f$acov))
  expect_gt(min(eigen(f$acov, symmetric = TRUE, only.values = TRUE)$values), -1e-8)

  # default is matrix-only; a bad level is rejected (base match.arg, as elsewhere in the pkg)
  expect_null(.polychoric(x)$acov)
  expect_error(.polychoric(x, acov = "nope"))
})

test_that("acov diag matches an independent influence-function computation", {
  # Simulated equicorrelated ordinal items with well-populated categories (so the threshold
  # bread is non-singular); compares the backend's diagonal to the from-scratch R reference.
  set.seed(3)
  N <- 1500L; pp <- 4L
  Lc <- chol(matrix(0.4, pp, pp) + diag(0.6, pp))
  Z <- matrix(stats::rnorm(N * pp), N, pp) %*% Lc
  thr <- list(stats::qnorm(c(.3, .65)), stats::qnorm(c(.25, .5, .75)),
              stats::qnorm(c(.2, .5, .8)), stats::qnorm(c(.4, .7)))
  x <- vapply(seq_len(pp), function(j) findInterval(Z[, j], thr[[j]]), integer(N))

  d <- .polychoric(x, acov = "diag")$acov
  codes <- apply(x, 2L, function(col) match(col, sort(unique(col))) - 1L)
  ref <- numeric(length(d)); t <- 1L
  for (i in seq_len(pp - 1L)) for (j in (i + 1L):pp) {
    ref[t] <- .if_diag_pair(codes[, i], codes[, j]); t <- t + 1L
  }
  expect_equal(unname(d), ref, tolerance = 1e-5)
})

test_that("acov uses the listwise-complete rows consistently under missingness", {
  g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), 1:6]
  gm <- g
  gm[1:60, 1] <- NA
  gm[61:110, 2] <- NA

  res <- .polychoric(gm, acov = "full")
  lw <- .polychoric(gm[stats::complete.cases(gm), ])

  # the returned matrix and the ACOV come from the SAME (listwise) case set
  expect_equal(res$R, lw$R)
  expect_true(all(is.finite(res$acov)))
  expect_equal(diag(res$acov), .polychoric(gm, acov = "diag")$acov, ignore_attr = TRUE)
})

test_that("acov aborts when there are too few listwise-complete cases", {
  m <- 60L
  # the two columns are never observed on the same row: no listwise-complete case exists
  X <- cbind(c(rep(1:3, length.out = m), rep(NA, m)),
             c(rep(NA, m), rep(1:3, length.out = m)))
  expect_error(.polychoric(X, acov = "diag"), class = "efa_cor_no_complete_cases")
  expect_error(.polychoric(X, acov = "full"), class = "efa_cor_no_complete_cases")
  # without an ACOV the same data is handled pairwise, failing as an uncomputable pair instead
  expect_error(.polychoric(X), class = "efa_cor_na")
})

test_that("acov full and diag match lavaan's NACOV on DOSPERT and UPPS", {
  skip_on_cran()
  skip_if_not_installed("lavaan")

  for (data in list(DOSPERT_raw, UPPS_raw)) {
    x <- data[stats::complete.cases(data), 1:6]
    f <- .polychoric(x, acov = "full")$acov
    lav <- .lav_acov(x)
    ord <- match(rownames(lav), rownames(f))
    ours <- f[ord, ord]

    # The residual is our two-step rho vs lavaan's optimiser plus the cancellation-free vs
    # corner cell-probability derivatives. Check the variances relatively and the (tiny,
    # ~1e-4) off-diagonal covariances absolutely -- mixing the two scales in one absolute
    # cap would let the variance-scale diagonal dominate. Measured agreement is ~3e-5 (diag,
    # relative) and ~1e-8 (off-diagonal, absolute), so these still catch a real sign/scale
    # error while leaving headroom against lavaan/BLAS version drift.
    expect_lt(max(abs(diag(ours) - diag(lav)) / diag(lav)), 1e-3)
    expect_lt(max(abs(ours[upper.tri(ours)] - lav[upper.tri(lav)])), 1e-7)
  }
})

test_that("acov full matches lavaan's NACOV on the full item sets", {
  skip_on_cran()
  skip_if_not_slow()
  skip_if_not_installed("lavaan")

  for (data in list(DOSPERT_raw, UPPS_raw)) {
    x <- data[stats::complete.cases(data), ]
    f <- .polychoric(x, acov = "full")$acov
    lav <- .lav_acov(x)
    ord <- match(rownames(lav), rownames(f))
    ours <- f[ord, ord]
    expect_lt(max(abs(diag(ours) - diag(lav)) / diag(lav)), 2e-3)
    expect_lt(max(abs(ours[upper.tri(ours)] - lav[upper.tri(lav)])), 1e-6)
  }
})

test_that("acov diag recovers the Monte-Carlo sampling variance of rho", {
  skip_on_cran()
  skip_if_not_slow()

  # Discretise a known bivariate normal into well-populated 4-category items; the analytic
  # diagonal variance must track the empirical Var(rho-hat) over many fresh samples. The
  # tolerance covers the Monte-Carlo error (~2% at R = 4000) plus the small outer-product
  # bias, with margin against RNG/platform drift.
  rho_true <- 0.5
  thr <- stats::qnorm(c(.25, .5, .75))
  disc <- function(z) findInterval(z, thr)
  N <- 2000L
  Lc <- chol(matrix(c(1, rho_true, rho_true, 1), 2L))

  set.seed(7)
  Z <- matrix(stats::rnorm(2 * N), N, 2L) %*% Lc
  analytic <- .polychoric(cbind(disc(Z[, 1]), disc(Z[, 2])), acov = "diag")$acov[[1]]

  set.seed(99)
  R <- 4000L
  rho_hat <- numeric(R)
  for (r in seq_len(R)) {
    Zb <- matrix(stats::rnorm(2 * N), N, 2L) %*% Lc
    rho_hat[r] <- .polychoric(cbind(disc(Zb[, 1]), disc(Zb[, 2])))$R[1, 2]
  }
  expect_lt(abs(analytic / stats::var(rho_hat) - 1), 0.15)
})
