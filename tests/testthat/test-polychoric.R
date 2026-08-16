# Polychoric / tetrachoric correlation matrix (.polychoric + the C++ backend).
#
# The two-step estimator is validated against the conditional-ML two-step references
# polycor::polychor(ML = FALSE) and psych::polychoric(correct = FALSE, global = FALSE);
# the thresholds match those functions' per-pair thresholds (the full-column marginals on
# complete data, the pairwise-complete marginals under missing data), so the estimates agree.
# Near-collinear pairs, where those fixed-quadrature references themselves lose accuracy, are
# instead checked against an mnormt two-step gold. The bivariate-normal rectangle that drives
# the likelihood is cross-checked against mnormt::sadmvn. Heavy reference loops over the large
# item sets run only under EFATOOLS_TEST_SLOW.

# Reference two-step polychoric matrix via polycor, pairwise.
.ref_polychor <- function(x) {
  p <- ncol(x)
  m <- diag(p)
  for (i in seq_len(p - 1L)) for (j in (i + 1L):p) {
    m[i, j] <- m[j, i] <- suppressWarnings(polycor::polychor(x[, i], x[, j], ML = FALSE))
  }
  m
}

# Exact two-step polychoric negative log-likelihood for one pair, integrating each cell with
# mnormt's adaptive bivariate-normal rule. Returned as a closure over the fixed thresholds and
# the contingency table, so a candidate estimate can be scored on the same objective
# .gold_polychor maximises.
.gold_nll <- function(a, b) {
  ai <- match(a, sort(unique(a))) - 1L; bi <- match(b, sort(unique(b))) - 1L
  Ki <- max(ai) + 1L; Kj <- max(bi) + 1L; N <- length(ai)
  ti <- stats::qnorm(cumsum(tabulate(ai + 1L, Ki)) / N)[-Ki]
  tj <- stats::qnorm(cumsum(tabulate(bi + 1L, Kj)) / N)[-Kj]
  rc <- c(-Inf, ti, Inf); cc <- c(-Inf, tj, Inf)
  nab <- matrix(0, Ki, Kj)
  for (r in seq_len(N)) nab[ai[r] + 1L, bi[r] + 1L] <- nab[ai[r] + 1L, bi[r] + 1L] + 1
  function(rho) {
    V <- matrix(c(1, rho, rho, 1), 2L); s <- 0
    for (i in seq_len(Ki)) for (j in seq_len(Kj)) if (nab[i, j] > 0)
      s <- s + nab[i, j] * log(max(as.numeric(mnormt::sadmvn(
        c(rc[i], cc[j]), c(rc[i + 1L], cc[j + 1L]), c(0, 0), V)), 1e-300))
    -s
  }
}

# Exact two-step polychoric for one pair: the maximiser of .gold_nll. This is the oracle in the
# near-collinear regime, where the fixed-quadrature reference estimators (polycor / psych)
# themselves lose accuracy.
.gold_polychor <- function(a, b) {
  stats::optimize(.gold_nll(a, b), c(-0.9999, 0.9999), tol = 1e-9)$minimum
}

# Expand a contingency table of counts into raw two-column ordinal data (row, column codes).
.expand_table <- function(tab) {
  do.call(rbind, lapply(which(tab > 0), function(k) {
    a <- ((k - 1L) %% nrow(tab)) + 1L
    b <- ((k - 1L) %/% nrow(tab)) + 1L
    cbind(rep(a, tab[a, b]), rep(b, tab[a, b]))
  }))
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
# need non-default arguments (`nearest_pd`, `acov`) or that mutate the data still
# build their own.
g <- GRiPS_raw[stats::complete.cases(GRiPS_raw), ]
poly_g <- .polychoric(g)

test_that(".bvn_rect_cpp matches mnormt::sadmvn rectangle probabilities", {
  skip_on_cran()
  skip_if_not_installed("mnormt")

  cuts <- c(-Inf, -2.5, -1, 0, 1, 2.5, Inf)
  worst <- 0
  for (rho in c(0, 0.3, 0.5, 0.75, 0.9, 0.95, 0.99, -0.8)) {
    V <- matrix(c(1, rho, rho, 1), 2L)
    for (ia in seq_len(length(cuts) - 1L)) for (ib in seq_len(length(cuts) - 1L)) {
      ours <- .bvn_rect_cpp(cuts[ia], cuts[ia + 1L], cuts[ib], cuts[ib + 1L], rho)
      ref <- as.numeric(mnormt::sadmvn(c(cuts[ia], cuts[ib]),
                                       c(cuts[ia + 1L], cuts[ib + 1L]),
                                       c(0, 0), V))
      worst <- max(worst, abs(ours - ref))
    }
  }
  # The quadrature reproduces the reference bivariate-normal integrator to better than 1e-6 over
  # this grid: the 12-node rule for |rho| up to 0.95, and the finer rule (used once |rho| exceeds
  # the refinement threshold) at rho = 0.99, where the conditional transition is narrow. This is
  # far tighter than the ~1e-4 agreement the polychoric matrix itself needs.
  expect_lt(worst, 1e-6)
})

test_that(".bvn_rect_cpp keeps its relative accuracy in the far tails", {
  # The absolute check above cannot see a 1e-45 rectangle collapse to zero, and that is exactly
  # what the opposite-corner cells of a strongly correlated pair are. Such a cell can be
  # OBSERVED -- a few careless or extreme responses put it there -- and then its probability
  # enters the log-likelihood and its score the asymptotic covariance, so what matters is
  # relative accuracy at magnitudes far below any absolute tolerance. The upper-tail rectangle
  # is the demanding one: its conditional band has to be taken from the complementary tail,
  # since as a plain difference of normal CDFs it is 1 - 1 = 0 and the rectangle vanishes
  # however large it really is.
  #
  # The oracle is the (X, Y) -> (-X, -Y) symmetry of the bivariate normal, which maps an
  # upper-tail rectangle onto its mirror image in the LOWER tail -- where a CDF difference
  # cannot cancel, because both terms are small. The identity is exact, needs no external
  # integrator (mnormt itself returns 0 for some of these rectangles), and holds to machine
  # precision here because the two sides are the same quadrature rule on mirrored bands.
  #
  # Compared as a RATIO, not as two values: expect_equal() on quantities this small would
  # compare them absolutely -- any two numbers below 1e-45 are within any usable tolerance of
  # each other -- and the assertion would hold whatever the upper tail returned. On the ratio the
  # tolerance means what it says, and the pre-fix values fail it by 6e-11 at rho = 0.9 and by
  # 7.7e-06 at rho = 0.95, where the cancelled band is inaccurate rather than absent.
  cut7 <- stats::qnorm(seq_len(6) / 7)
  # rho and the corner offset k; deeper than (0.999, 2) the rectangle genuinely underflows,
  # and at |rho| > 0.95 the finer quadrature rule is used, so both rules are covered.
  grid <- list(c(0.9, 1), c(0.9, 2), c(0.95, 1), c(0.95, 2),
               c(0.99, 1), c(0.99, 2), c(0.999, 2))
  for (gr in grid) {
    rho <- gr[1]; k <- gr[2]
    a <- cut7[k]; b <- cut7[7L - k]
    up <- .bvn_rect_cpp(-Inf, a, b, Inf, rho)          # X <= a, Y > b   (upper Y tail)
    lo <- .bvn_rect_cpp(-a, Inf, -Inf, -b, rho)        # X > -a, Y <= -b (its mirror image)
    expect_gt(up, 0)                                   # the cancelled band collapses to exactly 0
    expect_equal(up / lo, 1, tolerance = 1e-12)
  }
})

test_that("near-collinear pairs match the exact two-step estimate", {
  skip_on_cran()
  skip_if_not_installed("mnormt")

  # Near |rho| = 1 the contingency table has near-impossible off-diagonal cells whose
  # probability underflows; a naive expected-information step false-converges at the warm start,
  # and a fixed-order quadrature under-resolves the narrow conditional band. The estimator must
  # still recover the exact two-step value (and the reference estimators polycor / psych are
  # themselves biased here, so the oracle is the mnormt two-step).
  cut7 <- stats::qnorm(seq_len(6) / 7)
  set.seed(101)
  N <- 4000
  for (rt in c(0.97, 0.99, 0.995, 0.999)) {
    L <- chol(matrix(c(1, rt, rt, 1), 2L))
    Z <- matrix(stats::rnorm(2 * N), N, 2L) %*% L
    a <- findInterval(Z[, 1], cut7); b <- findInterval(Z[, 2], cut7)
    expect_equal(.polychoric(cbind(a, b))$R[1, 2], .gold_polychor(a, b), tolerance = 1e-4)
  }
})

test_that("binary pairs below the refinement gate hold the documented 1e-4 accuracy", {
  skip_on_cran()
  skip_if_not_installed("mnormt")

  # The band just below the gate at which the finer quadrature rule takes over is where the
  # fixed 12-node rule is least accurate, and a two-by-two table is its worst case: with only
  # one threshold per variable the whole likelihood rests on four wide bands whose conditional
  # transition the base rule under-resolves. The help page promises about 1e-4 there (against
  # about 1e-5 elsewhere), so that is what is pinned. The bound is set a few times the
  # measured discrepancy -- about 3e-5 at rho = .85 and 1.5e-4 at rho = .92 across seeds --
  # because it must hold on any BLAS, while still catching a loss of the refinement rule or a
  # coarser base rule, which cost an order of magnitude.
  gate_tol <- 5e-4
  for (rt in c(0.85, 0.92, -0.85, -0.92)) {
    for (seed in c(11, 12)) {
      set.seed(seed)
      L <- chol(matrix(c(1, rt, rt, 1), 2L))
      Z <- matrix(stats::rnorm(2 * 1500), 1500, 2L) %*% L
      a <- as.integer(Z[, 1] > 0); b <- as.integer(Z[, 2] > 0)
      ours <- .polychoric(cbind(a, b))$R[1, 2]
      expect_lt(abs(ours - .gold_polychor(a, b)), gate_tol)
      # Reversing one variable's coding reflects the latent variable about its threshold,
      # so the estimate must come back negated -- the accuracy above is a quadrature
      # residual, not an asymmetry between the two signs. Measured at exactly zero
      # difference here, but pinned loosely enough that a last-ulp difference in the
      # complementary threshold on another platform's qnorm cannot fail it, and still
      # six orders tighter than the accuracy band it has to distinguish.
      expect_equal(.polychoric(cbind(a, 1L - b))$R[1, 2], -ours, tolerance = 1e-10)
    }
  }

  # Unequal marginals put the single threshold off centre, which moves the bands the rule has
  # to resolve; the same accuracy has to hold there.
  set.seed(31)
  L <- chol(matrix(c(1, 0.92, 0.92, 1), 2L))
  Z <- matrix(stats::rnorm(2 * 1500), 1500, 2L) %*% L
  a <- as.integer(Z[, 1] > 0.6); b <- as.integer(Z[, 2] > 0.6)
  expect_lt(abs(.polychoric(cbind(a, b))$R[1, 2] - .gold_polychor(a, b)), gate_tol)
})

test_that("an empty pairwise-complete marginal category does not corrupt the estimate", {
  skip_on_cran()
  skip_if_not_installed("polycor")

  set.seed(303)
  N <- 1500
  L <- chol(matrix(c(1, .5, .5, 1), 2L))
  Z <- matrix(stats::rnorm(2 * N), N, 2L) %*% L
  a <- findInterval(Z[, 1], stats::qnorm(c(.2, .5, .8)))   # categories 0..3
  b <- findInterval(Z[, 2], stats::qnorm(c(.3, .7)))       # categories 0..2

  # Drop b wherever a is in its top (then bottom) category, so that category is absent from the
  # pair's pairwise-complete table while still present in a's full column. The per-pair
  # cumulative proportion then lands exactly on 1 (or 0), which must map to an infinite cut and
  # a zero-width empty band -- not a NaN that silently corrupts the estimate.
  for (empty in c(3L, 0L)) {
    bm <- b; bm[a == empty] <- NA
    ours <- .polychoric(cbind(a, bm))$R[1, 2]
    expect_false(is.na(ours))
    expect_equal(ours, suppressWarnings(polycor::polychor(a, bm, ML = FALSE)),
                 tolerance = 1e-3)
  }
})

test_that("an observed near-impossible cell does not bias a strongly-correlated pair", {
  skip_on_cran()
  skip_if_not_installed("mnormt")

  cut7 <- stats::qnorm(seq_len(6) / 7)
  set.seed(404)
  N <- 5000
  L <- chol(matrix(c(1, 0.995, 0.995, 1), 2L))
  Z <- matrix(stats::rnorm(2 * N), N, 2L) %*% L
  a <- findInterval(Z[, 1], cut7); b <- findInterval(Z[, 2], cut7)
  # A few careless responses fall in the two opposite corners, whose model probabilities at the
  # optimum are about 3e-45 and 1e-44. Both are far above the smallest representable double, but
  # the upper-right one is reached only from the complementary tail: taken as a plain difference
  # of normal CDFs its band is 1 - 1 = 0 exactly, which annihilates the cell and biases the pair
  # by ~4e-3 (the two mirror-image corners then contribute quite differently despite carrying the
  # same three cases each). Computed from the nearer tail, both corners contribute and the pair
  # tracks the exact two-step value to ~1e-5.
  a[1:3] <- 0L; b[1:3] <- 6L
  a[4:6] <- 6L; b[4:6] <- 0L

  # The asymptotic variance of this pair is available too: it is estimated in the interior, so it
  # is neither a boundary nor a continuity-corrected pair, and nothing in its influence function
  # overflows once the corner probabilities are accurate. Before that it came out NaN, which made
  # the pair unusable for DWLS and withheld every robust standard error involving it.
  # (The empty interior cells of this table trip the separate sparse-cell diagnostic, which is
  # computed from the observed counts and is not what is under test here.)
  d <- suppressWarnings(.polychoric(cbind(a, b), acov = "diag"),
                        classes = "efa_cor_sparse_cells")
  f <- suppressWarnings(.polychoric(cbind(a, b), acov = "full"),
                        classes = "efa_cor_sparse_cells")
  expect_equal(d$R[1, 2], .gold_polychor(a, b), tolerance = 1e-4)
  expect_false(d$at_bound[[1]])
  expect_false(d$zero_corrected[[1]])
  expect_true(is.finite(d$acov[[1]]) && d$acov[[1]] > 0)
  expect_equal(unname(d$acov), unname(diag(f$acov)), tolerance = 1e-10)
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

  # Pinned against the UNcorrected references even though binary pairs are continuity-corrected:
  # no pair of this median split has an empty cell, so the correction never fires and the two
  # conventions coincide here (psych gives the identical matrix for correct = 0.5 and FALSE).
  gb <- apply(g, 2L, function(col) as.integer(col > stats::median(col)))
  poly_gb <- .polychoric(gb)
  expect_false(any(poly_gb$zero_corrected))
  ours <- poly_gb$R

  expect_lt(max(abs(ours - .ref_polychor(gb))), 1e-4)

  tet_rho <- suppressWarnings(psych::tetrachoric(gb, correct = FALSE)$rho)
  expect_lt(max(abs(ours - unname(tet_rho))), 1e-4)
})

test_that("the matrix is a valid, named correlation matrix", {
  res <- poly_g

  expect_equal(dim(res$R), c(ncol(g), ncol(g)))
  expect_equal(dimnames(res$R), list(colnames(g), colnames(g)))
  expect_equal(diag(res$R), rep(1, ncol(g)), ignore_attr = TRUE)
  expect_true(isSymmetric(res$R))
  expect_true(all(res$R >= -1 & res$R <= 1))
})

test_that("the result is deterministic", {
  # Compared against the file-top `poly_g`, an independent earlier call with the same
  # arguments on the same object: two invocations separated by other work are a stronger
  # determinism probe than two adjacent ones, and one matrix fewer is built.
  expect_identical(poly_g$R, .polychoric(g)$R)
})

test_that("a constant column is rejected with a classed condition", {
  # local copy: this block mutates the column to force the constant-column path
  gx <- g
  gx[, 1L] <- 3L
  expect_error(.polychoric(gx), class = "efa_cor_constant_col")
})

test_that("a single variable is rejected at the R level with a classed condition", {
  # The C++ backend keeps the same check, but it raises an unclassed Rcpp error, so the
  # abort must happen in the wrapper.
  expect_error(.polychoric(g[, 1L, drop = FALSE]), class = "efa_cor_too_few_vars")
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

test_that("the boundary estimate of a Frechet-bound pair attains the maximum likelihood", {
  skip_if_not_installed("mnormt")
  # This table is the comonotone coupling of its margins -- the Frechet upper bound (Frechet,
  # 1951). On the cumulative scale that bound is the copula M(u, v) = min(u, v), which is exactly
  # the bivariate normal copula at rho = 1, and the two-step thresholds fix the model's margins to
  # this table's margins. So at rho = 1 the model reproduces the table exactly, the log-likelihood
  # attains its saturated value, and the maximiser is the boundary of the parameter space. The
  # estimator reports the boundary of the domain it maximises over (+/- 0.9999, as in polycor);
  # the likelihood is numerically flat over the whole approach to it, so an optimiser left to
  # search there stops at an arbitrary, platform-dependent point instead.
  #
  # What is checked here is the statistical claim behind that: the reported value fits at least as
  # well as the mnormt two-step gold, i.e. it really is a maximiser and not merely a bound. The
  # 1e-6 slack is measured, not inherited from sadmvn's abseps: at a perfect fit every cell's
  # n_ab / P_ab equals N, so an abseps-sized (1e-6) probability error would move this nll by as
  # much as sum(n_ab / P_ab) * abseps = 1.5e-3. sadmvn is far more accurate than it promises on
  # these smooth rectangles -- the measured spread over the whole plateau is ~6e-14 -- leaving
  # 1e-6 seven orders above the noise and three below the ~1.4e-3 a mislocated estimate pays.
  x <- .expand_table(rbind(c(41, 0), c(13, 0), c(28, 0), c(171, 47)))
  ours <- suppressWarnings(.polychoric(x))$R[1, 2]
  expect_identical(ours, 0.9999)
  nll <- .gold_nll(x[, 1], x[, 2])
  gold <- .gold_polychor(x[, 1], x[, 2])
  expect_lte(nll(ours), nll(gold) + 1e-6)
})

test_that("a Frechet-bound pair is detected exactly and reported at the boundary", {
  # The detector compares the table cell-by-cell against both bound couplings of its own
  # marginals, built from the cumulative margins by inclusion-exclusion. Every quantity is an
  # integer count, so the verdict involves no floating point and is identical on every platform --
  # which is the whole point: the estimate and the asymptotic variance of such a pair are boundary
  # quantities that an optimiser cannot pin down (the variance of this fixture has been observed
  # anywhere between 1e-13 and 1e+26 depending on the platform).
  upper <- rbind(c(41, 0), c(13, 0), c(28, 0), c(171, 47))

  # Upper bound (perfectly ordered) and its mirror, the lower bound (perfectly reversed): the
  # sign of the estimate follows the bound, and neither has an asymptotic variance.
  for (case in list(list(tab = upper,          rho = 0.9999),
                    list(tab = upper[, 2:1],   rho = -0.9999))) {
    res <- suppressWarnings(.polychoric(.expand_table(case$tab), acov = "diag"))
    expect_identical(res$R[1, 2], case$rho)
    expect_true(res$at_bound[[1L]])
    expect_identical(unname(res$acov), NA_real_)
  }

  # A 2x2 table with a zero cell is at a bound too, but it is repaired by the continuity
  # correction rather than reported (see the tetrachoric tests below). Two empty cells leave
  # nothing to repair -- the correction is defined for a single empty cell -- so those tables
  # stay on the boundary path, in both directions.
  for (case in list(list(tab = rbind(c(40, 0), c(0, 200)), rho = 0.9999),
                    list(tab = rbind(c(0, 40), c(200, 0)), rho = -0.9999))) {
    res <- suppressWarnings(.polychoric(.expand_table(case$tab), acov = "diag"))
    expect_identical(res$R[1, 2], case$rho)
    expect_true(res$at_bound[[1L]])
    expect_false(res$zero_corrected[[1L]])
    expect_identical(unname(res$acov), NA_real_)
  }

  # One count off the bound and the pair is ordinary again: not flagged, estimated in the
  # interior, and with a finite positive variance. There is no grey zone between the two.
  off <- suppressWarnings(.polychoric(.expand_table(rbind(c(40, 2), c(60, 198))), acov = "diag"))
  expect_false(off$at_bound[[1L]])
  expect_false(off$zero_corrected[[1L]])
  expect_lt(abs(off$R[1, 2]), 0.95)
  expect_true(is.finite(off$acov[[1L]]) && off$acov[[1L]] > 0)
})

test_that("a Frechet-bound pair is reproducible under a one-count perturbation", {
  # Moving a single observation between populated cells leaves the table at the bound (the
  # marginals change, but the monotone rearrangement of the new marginals is again this table), so
  # every member of this family is the same boundary problem. It is precisely where the old
  # arbitrary-stopping-point estimate made the asymptotic variance vary by 39 orders of magnitude
  # across the family, non-monotonically, and differ between platforms. All of them must now give
  # the identical answer. The probe is only meaningful among perturbations that stay on the same
  # side of the bound -- one that moved a count INTO a structurally empty cell would leave the
  # family, and is expected to be estimated in the interior instead (covered above).
  tab <- rbind(c(41, 0), c(13, 0), c(28, 0), c(171, 47))
  for (k in which(tab > 0)) {
    for (step in c(-1L, 1L)) {
      t2 <- tab
      t2[k] <- t2[k] + step
      res <- suppressWarnings(.polychoric(.expand_table(t2), acov = "diag"))
      expect_identical(res$R[1, 2], 0.9999)
      expect_identical(unname(res$acov), NA_real_)
    }
  }
})

test_that("a binary pair with one empty cell is continuity-corrected, not reported at the bound", {
  # Every 2x2 table with an empty cell is the extreme coupling of its own margins, so without a
  # correction its estimate would be +/-1 whatever the underlying correlation -- which for binary
  # data is the common case, not a corner case. A continuity correction of 0.5 is added to the
  # empty cell instead, so the pair is estimated in the interior. The sign still follows the
  # corner the empty cell sits in.
  for (case in list(list(tab = rbind(c(40, 0), c(60, 200)),  sign = 1),
                    list(tab = rbind(c(0, 40), c(200, 60)),  sign = -1),
                    list(tab = rbind(c(150, 50), c(100, 0)), sign = -1))) {
    res <- suppressWarnings(.polychoric(.expand_table(case$tab), acov = "diag"))
    expect_false(res$at_bound[[1L]])
    expect_true(res$zero_corrected[[1L]])
    expect_identical(sign(res$R[1, 2]), case$sign)
    expect_lt(abs(res$R[1, 2]), 0.9999)
    # The correction repairs the point estimate only: the sandwich would treat the nudged counts
    # as data, and its coverage in this regime is not trustworthy, so no variance is reported.
    expect_identical(unname(res$acov), NA_real_)
  }

  # The mirror image must give the mirrored estimate, so the correction cannot depend on which
  # corner is empty. Compared numerically rather than bit-wise: the two tables put the correction
  # in different cells and so take slightly different routes through the iteration, which leaves
  # them a single ulp apart.
  up <- suppressWarnings(.polychoric(.expand_table(rbind(c(40, 0), c(60, 200)))))$R[1, 2]
  dn <- suppressWarnings(.polychoric(.expand_table(rbind(c(0, 40), c(200, 60)))))$R[1, 2]
  expect_equal(up, -dn, tolerance = 1e-12)
})

test_that("the continuity correction preserves the table margins", {
  # The correction adds 0.5 to the empty cell AND to its diagonal opposite, and takes 0.5 out of
  # the other two, so every row and column sum is unchanged. That is what keeps the shared
  # per-variable thresholds valid, and it is checked here independently: with the margins fixed a
  # 2x2 has one degree of freedom, so the two-step estimate is the rho solving
  # Phi2(t1, t2, rho) = (a - 0.5) / N with t1, t2 taken from the ORIGINAL margins. Agreement to
  # quadrature accuracy shows the margins -- hence the thresholds -- really were left alone.
  tab <- rbind(c(40, 0), c(60, 200))
  N <- sum(tab)
  t1 <- stats::qnorm(sum(tab[1, ]) / N)
  t2 <- stats::qnorm(sum(tab[, 1]) / N)
  target <- (tab[1, 1] - 0.5) / N
  gold <- stats::uniroot(function(rho) .bvn_rect_cpp(-Inf, t1, -Inf, t2, rho) - target,
                         c(-0.9999, 0.9999), tol = 1e-10)$root
  expect_equal(suppressWarnings(.polychoric(.expand_table(tab)))$R[1, 2], gold, tolerance = 1e-4)
})

test_that("the continuity correction preserves the margins under missing data too", {
  # With missing data each pair uses the thresholds of its OWN pairwise-complete cases, derived
  # from the contingency table -- the corrected one. That is the case the margin-preserving design
  # exists for: a correction that shifted the margins would move the thresholds with them. Build a
  # pair whose complete cases form a 2x2 with one empty cell, add missing values in a third
  # variable so the local-threshold branch is taken, and check the pair still matches the estimate
  # implied by its ORIGINAL margins.
  tab <- rbind(c(40, 0), c(60, 200))
  xy <- .expand_table(tab)
  x <- cbind(xy, rep(c(1L, 2L, NA), length.out = nrow(xy)))
  colnames(x) <- c("a", "b", "z")

  res <- suppressWarnings(.polychoric(x))
  expect_true(res$zero_corrected[["a-b"]])

  N <- sum(tab)
  t1 <- stats::qnorm(sum(tab[1, ]) / N)
  t2 <- stats::qnorm(sum(tab[, 1]) / N)
  gold <- stats::uniroot(
    function(rho) .bvn_rect_cpp(-Inf, t1, -Inf, t2, rho) - (tab[1, 1] - 0.5) / N,
    c(-0.9999, 0.9999), tol = 1e-10)$root
  expect_equal(res$R[1, 2], gold, tolerance = 1e-4)
})

test_that("the corrected tetrachoric matches lavaan, which applies the same correction", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  # lavaan corrects a 2x2 zero cell by default (`zero.add = c(0.5, 0)`) and keeps the margins
  # (`zero.keep.margins`), which is exactly the convention implemented here -- so the estimates
  # must agree. With `zero.add = 0` lavaan instead returns its own `maxcor` of 0.999, confirming
  # that the two estimators differ only by the correction and not in the estimator itself.
  tab <- rbind(c(40, 0), c(60, 200))
  d <- as.data.frame(.expand_table(tab))
  names(d) <- c("a", "b")
  d[] <- lapply(d, ordered)

  ours <- suppressWarnings(.polychoric(.expand_table(tab)))$R[1, 2]
  ref <- lavaan::lavCor(d, ordered = names(d), output = "cor")[1, 2]
  expect_equal(ours, ref, tolerance = 1e-4)
  # Regression pin, so a change in the correction is caught even without lavaan installed.
  expect_equal(ours, 0.9276913, tolerance = 1e-6)
})

# Collect the class strings of every warning raised while evaluating `expr`, so one expected class
# can be asserted (and another denied) when several conditions fire on the same data.
.poly_warn_classes <- function(expr) {
  cls <- character()
  withCallingHandlers(
    expr,
    warning = function(w) { cls <<- c(cls, class(w)); invokeRestart("muffleWarning") }
  )
  cls
}

# Variables with structurally empty cells whose expected counts are large (32), estimated well
# inside the interior (rho = 0.82) rather than at a Frechet bound: a common ordinal score nudged up
# by each variable's own 0/1 pattern and capped at the top category. Reaching the lowest category
# needs an unnudged zero score and the highest a nudged or already-top score, so cells (0, 2) and
# (2, 0) are both empty in every pair's table -- a pattern no extreme coupling has, which is what
# keeps the pairs off the bound.
.sparse_interior <- function(n_var = 2L) {
  r <- 0:383L
  u <- (r %/% 64L) %% 3L
  out <- vapply(seq_len(n_var),
                function(k) pmin(2L, u + bitwAnd(bitwShiftR(r, k - 1L), 1L)),
                integer(length(r)))
  colnames(out) <- paste0("v", seq_len(n_var))
  out
}

test_that("a structurally empty cell warns when an asymptotic covariance is requested", {
  # An empty cell with a non-negligible expected count is flagged; an otherwise dense ordinal set
  # whose only empty cells are rare corners is not. The covariance is still returned.
  sparse <- .sparse_interior()
  expect_warning(.polychoric(sparse, acov = "diag"), class = "efa_cor_sparse_cells")
  expect_no_warning(
    .polychoric(DOSPERT_raw[, 1:8], acov = "diag"),
    class = "efa_cor_sparse_cells"
  )
  # The per-replicate bootstrap recompute (label_acov = FALSE) skips the diagnostic entirely.
  expect_no_warning(
    .polychoric(sparse, acov = "diag", label_acov = FALSE),
    class = "efa_cor_sparse_cells"
  )
  # A pair at a Frechet bound has empty cells too, but it is reported as a boundary solution
  # instead: that warning says strictly more (there is no asymptotic covariance at all, rather
  # than an unreliable one), so repeating the weaker sparse-cell diagnosis for it would only
  # obscure the diagnosis.
  bound <- .expand_table(rbind(c(41, 0), c(13, 0), c(28, 0), c(171, 47)))
  cls <- .poly_warn_classes(.polychoric(bound, acov = "diag"))
  expect_true("efa_cor_boundary" %in% cls)
  expect_false("efa_cor_sparse_cells" %in% cls)
  # The bootstrap recompute is silent about the boundary too.
  expect_no_warning(.polychoric(bound, acov = "diag", label_acov = FALSE))
})

test_that("the sparse-cell warning names the offending pairs and caps the list", {
  testthat::local_reproducible_output()

  # The warning itself is the subject here, so catch it rather than snapshotting whatever
  # else the fixture emits.
  sparse_msg <- function(x) {
    w <- tryCatch(.polychoric(x, acov = "diag"),
                  efa_cor_sparse_cells = function(w) w)
    expect_s3_class(w, "efa_cor_sparse_cells")
    conditionMessage(w)
  }

  one <- .sparse_interior()
  colnames(one) <- c("risk", "gamble")
  expect_snapshot(cat(sparse_msg(one)))

  # The same construction over six items: cells (0, 2) and (2, 0) are structurally empty in all
  # 15 tables, each with an expected count of 32, and the nudges keep every pair well short of
  # perfect dependence (rho = 0.82) so none of them lands on a Frechet bound -- which the naive
  # monotone construction would. The reported list is capped at five plus a count.
  nudged <- .sparse_interior(6L)
  expect_snapshot(cat(sparse_msg(nudged)))
})

test_that("the boundary warning names the offending pairs and caps the list", {
  testthat::local_reproducible_output()

  boundary_msg <- function(x) {
    w <- tryCatch(.polychoric(x, acov = "diag"), efa_cor_boundary = function(w) w)
    expect_s3_class(w, "efa_cor_boundary")
    conditionMessage(w)
  }

  upper <- rbind(c(41, 0), c(13, 0), c(28, 0), c(171, 47))
  one <- .expand_table(upper)
  colnames(one) <- c("risk", "gamble")
  expect_snapshot(cat(boundary_msg(one)))

  # A perfectly REVERSED table is estimated at the negative endpoint, so the snapshot must show
  # -0.9999. Reporting the positive value here, or describing the pair as perfectly ordered, would
  # contradict the correlation the same call returns.
  rev <- .expand_table(upper[, 2:1])
  colnames(rev) <- c("risk", "gamble")
  expect_snapshot(cat(boundary_msg(rev)))

  # Both directions in one matrix: the report names both values rather than picking one. `gamble`
  # is monotone increasing in `risk` and `mirror` is its reverse, so the two pairs involving
  # `mirror` sit on the opposite bound from the `risk`-`gamble` pair (the `gamble`-`mirror` table
  # has two empty cells, so it is not continuity-corrected and stays a boundary pair).
  both <- cbind(one, mirror = max(one[, 2L]) - one[, 2L])
  expect_snapshot(cat(boundary_msg(both)))
})

test_that("the boundary warning quotes the estimate, not the nearest-PD projection of it", {
  # The reported value comes from the backend's own boundary estimate rather than from the returned
  # matrix, because `nearest_pd` projects the matrix afterwards and moves the entry off the
  # endpoint. Without this the message would quote the projected value (measured here: 0.9954538
  # instead of 0.9999), so the check is what keeps the estimate and the message from drifting apart.
  x <- cbind(.expand_table(rbind(c(41, 0), c(13, 0), c(28, 0), c(171, 47))),
             .sparse_interior(2L)[seq_len(300L), ])
  colnames(x) <- paste0("v", 1:4)

  projected <- suppressWarnings(.polychoric(x, nearest_pd = TRUE))$R[1, 2]
  expect_lt(projected, 0.9999)                      # the projection really did move it

  # The projection also warns that it happened; that is asserted elsewhere, so muffle it here and
  # let only the boundary report through.
  w <- tryCatch(
    suppressWarnings(.polychoric(x, nearest_pd = TRUE), classes = "efa_cor_smoothed"),
    efa_cor_boundary = function(w) w)
  expect_s3_class(w, "efa_cor_boundary")
  expect_snapshot(cat(conditionMessage(w)))
})

test_that("the continuity-correction warning names the offending pairs and caps the list", {
  testthat::local_reproducible_output()

  zero_cell_msg <- function(x) {
    w <- tryCatch(.polychoric(x, acov = "diag"), efa_cor_zero_cell = function(w) w)
    expect_s3_class(w, "efa_cor_zero_cell")
    conditionMessage(w)
  }

  # Seven binary items that are all monotone in one another (an item is answered 1 only if every
  # easier item is), so every one of the 21 tables has a single empty cell. Being 2x2 they take
  # the continuity correction rather than the boundary report, which is what separates this case
  # from the larger table above. The reported list is capped at five plus a count.
  guttman <- vapply(1:7, function(k) as.integer((0:279L) %% 8L >= k), integer(280L))
  colnames(guttman) <- paste0("i", 1:7)
  expect_snapshot(cat(zero_cell_msg(guttman)))

  # The pairs are corrected, not bounded, and none of them carries a variance.
  res <- suppressWarnings(.polychoric(guttman, acov = "diag"))
  expect_true(all(res$zero_corrected))
  expect_false(any(res$at_bound))
  expect_true(all(is.na(res$acov)))
})

test_that("a likely-continuous (many-category) variable warns", {
  x <- cbind(rep(1:12, length.out = 100L), rep(1:3, length.out = 100L))
  expect_warning(.polychoric(x), class = "efa_cor_many_categories")
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

test_that("a Frechet-bound pair has no asymptotic covariance and contaminates no other pair", {
  # A boundary pair has no asymptotic variance to report: the rho-information vanishes identically
  # at the bound, so the influence function (score - threshold correction) / A22 is 0/0 -- both
  # parts underflow to zero well before the boundary is reached, and the limit of the ratio is
  # infinite. Recomputed from scratch at high accuracy, A22 falls from 1.6e-03 at rho = 0.90 to
  # 3.8e-16 at 0.97 and to exactly 0 beyond, so any finite value here is an artefact of the
  # probability floor rather than a variance. It is fail-closed to NA at both levels, which is what
  # makes the branch every consumer takes independent of the platform.
  #
  # The other five pairs are ordinary and must come through untouched: the boundary pair's
  # influence column is zeroed before the crossprod precisely so its NAs cannot spread into them.
  x <- cbind(.expand_table(rbind(c(41, 0), c(13, 0), c(28, 0), c(171, 47))),
             .sparse_interior(2L)[seq_len(300L), ])
  colnames(x) <- paste0("v", 1:4)
  d <- suppressWarnings(.polychoric(x, acov = "diag"))
  f <- suppressWarnings(.polychoric(x, acov = "full"))$acov

  bnd <- d$at_bound
  expect_identical(unname(which(bnd)), 1L)                 # the v1-v2 pair, in combn order
  expect_true(all(is.na(d$acov[bnd])))
  expect_true(all(is.finite(d$acov[!bnd]) & d$acov[!bnd] > 0))

  # The full level withholds the boundary pair's whole row and column -- its covariance with every
  # other pair rests on the same vanishing information -- and leaves the rest estimable.
  expect_true(all(is.na(f[bnd, ])) && all(is.na(f[, bnd])))
  expect_true(all(is.finite(f[!bnd, !bnd])))

  # The cheap diagonal still equals diag(full), NA entries included: the two are computed by
  # different routes (a direct cell sum versus scatter-then-crossprod).
  expect_equal(unname(d$acov), unname(diag(f)), tolerance = 1e-6)
})

test_that("a continuity-corrected pair has no asymptotic covariance at either level", {
  # A corrected pair is withheld for a different reason than a boundary pair -- its estimate is
  # interior, but it comes from a nudged table the sandwich would treat as data -- and it is
  # withheld at BOTH levels, not just the cheap diagonal. The other pairs must still come through
  # estimable, which is what the zeroed influence column before the crossprod is for.
  x <- cbind(.expand_table(rbind(c(40, 0), c(60, 200))),
             .sparse_interior(2L)[seq_len(300L), ])
  colnames(x) <- paste0("v", 1:4)
  d <- suppressWarnings(.polychoric(x, acov = "diag"))
  f <- suppressWarnings(.polychoric(x, acov = "full"))$acov

  zc <- d$zero_corrected
  expect_identical(unname(which(zc)), 1L)                  # the v1-v2 pair, in combn order
  expect_false(any(d$at_bound))
  expect_true(all(is.na(d$acov[zc])))
  expect_true(all(is.finite(d$acov[!zc]) & d$acov[!zc] > 0))

  expect_true(all(is.na(f[zc, ])) && all(is.na(f[, zc])))
  expect_true(all(is.finite(f[!zc, !zc])))
  expect_equal(unname(d$acov), unname(diag(f)), tolerance = 1e-6)
})

test_that("a near-collinear pair with structurally empty tail cells still has a variance", {
  skip_on_cran()

  # A pair correlated this strongly has empty cells whose model probability is astronomically
  # small -- at rho = 0.999 some are below the smallest representable double outright. Those
  # cells carry no cases, so they must contribute exactly nothing to the influence function: the
  # asymptotic variance has to come out finite and positive, and the cheap diagonal has to agree
  # with the crossprod route. The pair is estimated in the interior (its table is off both
  # Frechet bounds and it is not binary), so nothing upstream withholds the variance and the
  # whole assembly is exercised. The estimates themselves are pinned against the exact two-step
  # value elsewhere; this covers the ACOV those same tables never requested.
  #
  # Same seed and same four rt values as that estimate test, so these are the same four tables:
  # each draw depends on the ones before it, and a shorter loop would silently substitute
  # different data.
  cut7 <- stats::qnorm(seq_len(6) / 7)
  set.seed(101)
  N <- 4000
  for (rt in c(0.97, 0.99, 0.995, 0.999)) {
    L <- chol(matrix(c(1, rt, rt, 1), 2L))
    Z <- matrix(stats::rnorm(2 * N), N, 2L) %*% L
    a <- findInterval(Z[, 1], cut7); b <- findInterval(Z[, 2], cut7)
    d <- suppressWarnings(.polychoric(cbind(a, b), acov = "diag"),
                          classes = "efa_cor_sparse_cells")
    f <- suppressWarnings(.polychoric(cbind(a, b), acov = "full"),
                          classes = "efa_cor_sparse_cells")
    expect_false(d$at_bound[[1]])
    expect_false(d$zero_corrected[[1]])
    expect_true(is.finite(d$acov[[1]]) && d$acov[[1]] > 0)
    expect_equal(unname(d$acov), unname(diag(f$acov)), tolerance = 1e-10)
  }
})

# Deterministic ordinal data whose alphabetical label order differs from the response order:
# a bivariate normal (rho = 0.7) cut into four categories 1..4, then labelled in the TRUE
# order never < rarely < often < always. Alphabetical sorting scrambles that to
# always, never, often, rarely, so an unordered factor / character column would attenuate rho.
.labelled_ordinal <- function() {
  set.seed(123)
  N <- 600L
  z1 <- stats::rnorm(N)
  z2 <- 0.7 * z1 + sqrt(1 - 0.7^2) * stats::rnorm(N)
  br <- stats::qnorm(c(.25, .5, .75))
  a <- findInterval(z1, br) + 1L
  b <- findInterval(z2, br) + 1L
  lab <- c("never", "rarely", "often", "always")        # the true response order
  list(
    num = data.frame(a = a, b = b),
    ord = data.frame(a = ordered(lab[a], levels = lab), b = ordered(lab[b], levels = lab)),
    uno = data.frame(a = factor(lab[a]), b = factor(lab[b])),
    chr = data.frame(a = lab[a], b = lab[b], stringsAsFactors = FALSE)
  )
}

test_that("ordered factors match numeric codes but unordered / character columns are rejected", {
  d <- .labelled_ordinal()

  # An ordered factor carries the response order, so it gives exactly the numeric-code answer.
  expect_equal(.polychoric(d$ord)$R, .polychoric(d$num)$R)

  # An unordered factor or character column has no response order; data.matrix() would rank it
  # alphabetically and silently attenuate the correlation, so it is refused up front.
  expect_error(.polychoric(d$uno), class = "efa_cor_unordered_factor")
  expect_error(.polychoric(d$chr), class = "efa_cor_unordered_factor")
  expect_error(.polychoric(as.matrix(d$chr)), class = "efa_cor_unordered_factor")
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

  res <- suppressWarnings(
    .polychoric(gm, acov = "full"),
    classes = "efa_cor_sparse_cells"
  )
  lw <- .polychoric(gm[stats::complete.cases(gm), ])

  # the returned matrix and the ACOV come from the SAME (listwise) case set
  expect_equal(res$R, lw$R)
  expect_true(all(is.finite(res$acov)))
  expect_equal(
    diag(res$acov),
    suppressWarnings(
      .polychoric(gm, acov = "diag"),
      classes = "efa_cor_sparse_cells"
    )$acov,
    ignore_attr = TRUE
  )
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
    f <- suppressWarnings(
      .polychoric(x, acov = "full"),
      classes = "efa_cor_sparse_cells"
    )$acov
    lav <- suppressWarnings(.lav_acov(x))
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
    f <- suppressWarnings(
      .polychoric(x, acov = "full"),
      classes = "efa_cor_sparse_cells"
    )$acov
    lav <- suppressWarnings(.lav_acov(x))
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
