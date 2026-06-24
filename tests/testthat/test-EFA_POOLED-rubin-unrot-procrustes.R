# Procrustes-aligned analytic-SE pooling in EFA_POOLED().
#
# Companion to test-EFA_POOLED-rubin-unrot.R, focused on
# align_unrotated = "procrustes". The Procrustes orthogonal transform Q_d
# mixes loading columns, so the per-imputation full unrotated vcov is
# propagated through the Kronecker identity
#
#     Var(vec(L_d Q_d)) = (Q_d' (x) I_p) V_d (Q_d (x) I_p)
#
# (column-major vec) before per-element Rubin pooling. Q_d is treated as fixed
# (Schoenemann 1966; lavaan / mice convention).
#
# Assertion groups:
#   1. Identity-imputation invariance.
#   2. Signed-permutation orbit: procrustes pool recovers the anchor SE.
#   3. Helper kronecker check: closed-form aligned diagonals match the full
#      Kronecker propagation and a Monte-Carlo simulation.
#   4. Continuous-rotation orbit: procrustes pool recovers the anchor SE.
#   5. Cross-check that the procrustes pool and the signed_tucker_congruence
#      pool agree on a signed-permutation orbit (both reduce to the anchor SE).
#   6. NA-vcov gating: NA-filled vcov on any imputation aborts with
#      `efa_pooled_unreliable_vcov`; per-row NA in the helper still
#      NA-propagates when the helper is called directly.
#   7. Classed `efa_pooled_no_vcov` abort when any component fit is missing
#      `vcov_unrot_loadings` under `align_unrotated = "procrustes"`.

cormat <- test_models$baseline$cormat
p_vars <- ncol(cormat)
N_id   <- 500L
k      <- 3L
m_id   <- 5L
identical_list <- replicate(m_id, cormat, simplify = FALSE)

# -- Synthetic-fit helpers ----------------------------------------------------

# Build the (p k x p k) Kronecker that maps vec(L) -> vec(L S) under
# column-major vec: vec(L S) = (S' (x) I_p) vec(L). Used to transport the
# anchor's vcov_unrot_loadings into a known gauge so the test inputs satisfy
# the same Kronecker identity the production helper relies on.
.k_lhs <- function(S, p) kronecker(t(S), diag(p))

# Marginal SEs of (L S) implied by V_d in the gauged frame.
.row_marginal_se <- function(V, p, k) {
  out <- matrix(NA_real_, p, k)
  for (i in seq_len(p)) {
    idx <- (seq_len(k) - 1L) * p + i
    out[i, ] <- sqrt(pmax(diag(V[idx, idx, drop = FALSE]), 0))
  }
  out
}

# Apply a k x k gauge S to a base fit: replace the unrotated loadings,
# unrot_loading SEs, and full vcov to reflect the (L S, (S' (x) I) V (S (x) I))
# transform. Uniquenesses (rotation-invariant) are kept identical to the base.
.gauge_fit <- function(base_fit, S) {
  L_base <- unclass(base_fit$unrot_loadings)
  p <- nrow(L_base)
  kk <- ncol(L_base)

  L_g <- L_base %*% S

  V_base <- base_fit$vcov_unrot_loadings
  J <- .k_lhs(S, p)
  V_g <- J %*% V_base %*% t(J)
  V_g <- (V_g + t(V_g)) / 2

  SE_g <- .row_marginal_se(V_g, p, kk)
  dimnames(SE_g) <- dimnames(L_base)

  fit <- base_fit
  fit$unrot_loadings <- structure(L_g,
    class    = class(base_fit$unrot_loadings),
    dimnames = dimnames(base_fit$unrot_loadings)
  )
  fit$SE$unrot_loadings <- SE_g
  fit$vcov_unrot_loadings <- V_g
  fit
}

# Rotation matrix in factor space (k x k orthogonal). For k > 2, a single
# Givens rotation in the (1, 2) plane is enough to mix columns; the identity
# elsewhere keeps the test focused on the Kronecker propagation rather than
# k-specific structure.
.givens <- function(theta, k, i = 1L, j = 2L) {
  R <- diag(k)
  c_ <- cos(theta); s_ <- sin(theta)
  R[i, i] <-  c_; R[j, j] <- c_
  R[i, j] <- -s_; R[j, i] <- s_
  R
}

# Signed permutation matrix S with S[perm[j], j] = signs[j] (so L S permutes
# and sign-flips columns of L).
.signed_perm <- function(perm, signs) {
  k <- length(perm)
  S <- matrix(0, k, k)
  for (j in seq_len(k)) S[perm[j], j] <- signs[j]
  S
}

# ---- Group 1 ---------------------------------------------------------------

test_that("identity imputations: procrustes pool equals the per-fit SE and FMI is zero", {
  pooled <- suppressMessages(EFA_POOLED(
    identical_list, n_factors = k, N = N_id, method = "ML",
    rotation = "none", se = "information",
    align_unrotated = "procrustes"
  ))

  oracle <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                rotation = "none", se = "information")

  expect_identical(pooled$settings$align_unrotated, "procrustes")
  expect_identical(pooled$settings$se, "information")

  expect_equal(pooled$SE$unrot_loadings, oracle$SE$unrot_loadings, tolerance = 1e-10)
  expect_equal(pooled$SE$uniquenesses,   oracle$SE$uniquenesses,   tolerance = 1e-10)
  expect_equal(unclass(pooled$unrot_loadings), unclass(oracle$unrot_loadings),
               tolerance = 1e-10)

  # Q_d = I_k under identical imputations only up to LAPACK SVD numerical
  # noise. Tolerate that noise; the signed_tucker_congruence pool gives an
  # exact 0 because it permutes integer column indices, not floats.
  fmi <- pooled$MI$unrot_loadings$FMI
  riv <- pooled$MI$unrot_loadings$RIV
  expect_true(all(abs(fmi[is.finite(fmi)]) < 1e-10))
  expect_true(all(abs(riv[is.finite(riv)]) < 1e-10))
})

# ---- Group 2 ---------------------------------------------------------------

test_that("procrustes pool over a signed-permutation orbit returns the anchor SE", {
  base_fit <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                  rotation = "none", se = "information")
  SE_base <- base_fit$SE$unrot_loadings
  psi_se  <- base_fit$SE$uniquenesses

  gauges <- list(
    diag(k),
    .signed_perm(c(2L, 1L, 3L), c(-1,  1,  1)),
    .signed_perm(c(3L, 2L, 1L), c( 1, -1,  1)),
    .signed_perm(c(1L, 3L, 2L), c( 1,  1, -1)),
    .signed_perm(c(2L, 3L, 1L), c(-1, -1,  1))
  )
  gauge_fits <- lapply(gauges, .gauge_fit, base_fit = base_fit)

  unrot_list <- lapply(gauge_fits, function(f) unclass(f$unrot_loadings))
  aligned <- EFAtools:::.efa_pooled_align_unrotated_list(
    unrot_list,
    align_unrotated = "procrustes",
    return_meta = TRUE
  )

  pool <- EFAtools:::.efa_pooled_analytic_pool(
    fits = gauge_fits,
    unrot_loadings_aligned = aligned$loadings,
    align_meta = aligned$meta,
    ci = 0.95,
    align_unrotated = "procrustes"
  )

  expect_equal(unname(pool$SE$unrot_loadings), unname(SE_base), tolerance = 1e-10)
  expect_equal(unname(pool$SE$uniquenesses),   unname(psi_se),  tolerance = 1e-12)
})

# ---- Group 3 ---------------------------------------------------------------

test_that(".efa_pooled_propagate_procrustes_vcov matches the full kronecker formula and a Monte-Carlo estimate", {
  skip_if_not_installed("MASS")
  set.seed(20260621)
  p <- 4L
  kk <- 2L

  # PSD V of dimension (p k) x (p k); rescaling keeps numbers in O(1e-2) so
  # the MC variances do not get drowned by row-block magnitude mismatches.
  A <- matrix(stats::rnorm((p * kk)^2), p * kk, p * kk)
  V <- crossprod(A) / (p * kk)

  theta <- pi / 7
  Q <- .givens(theta, kk)

  helper <- EFAtools:::.efa_pooled_propagate_procrustes_vcov(V, Q, p, kk)

  # Brute-force Kronecker reference: V_aligned = (Q' (x) I_p) V (Q (x) I_p);
  # row-marginal variances are the (i + (j - 1) p)-th diagonal element.
  J <- kronecker(t(Q), diag(p))
  V_aligned <- J %*% V %*% t(J)
  ref <- matrix(diag(V_aligned), p, kk)
  expect_equal(helper, sqrt(pmax(ref, 0)), tolerance = 1e-12)

  # Monte-Carlo sanity: draw vec(L) ~ N(0, V), align by Q, take the empirical
  # row-marginal variance. Use a sample of 5e4; the MC standard error of a
  # variance estimate is sqrt(2 / n) ~ 6e-3, so 2e-2 is a safe tolerance.
  n_draws <- 50000L
  draws <- MASS::mvrnorm(n_draws, mu = rep(0, p * kk), Sigma = V)
  mc_var <- matrix(0, p, kk)
  for (b in seq_len(n_draws)) {
    L_b <- matrix(draws[b, ], p, kk)
    mc_var <- mc_var + (L_b %*% Q)^2
  }
  mc_var <- mc_var / n_draws
  expect_equal(mc_var, ref, tolerance = 2e-2)
})

# ---- Group 4 ---------------------------------------------------------------

test_that("procrustes pool recovers the anchor SE under a continuous rotation gauge", {
  base_fit <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                  rotation = "none", se = "information")
  SE_base <- base_fit$SE$unrot_loadings

  rotations <- list(
    diag(k),
    .givens(pi / 7,  k, 1L, 2L),
    .givens(pi / 5,  k, 2L, 3L),
    .givens(pi / 11, k, 1L, 3L)
  )
  gauge_fits <- lapply(rotations, .gauge_fit, base_fit = base_fit)

  unrot_list <- lapply(gauge_fits, function(f) unclass(f$unrot_loadings))
  aligned <- EFAtools:::.efa_pooled_align_unrotated_list(
    unrot_list,
    align_unrotated = "procrustes",
    return_meta = TRUE
  )

  pool <- EFAtools:::.efa_pooled_analytic_pool(
    fits = gauge_fits,
    unrot_loadings_aligned = aligned$loadings,
    align_meta = aligned$meta,
    ci = 0.95,
    align_unrotated = "procrustes"
  )

  # All aligned vcovs collapse to V_base, so B = 0 and the pool returns the
  # anchor's marginal SE matrix up to floating-point noise from the SVD-based
  # Q recovery and Kronecker round-trip.
  expect_equal(unname(pool$SE$unrot_loadings), unname(SE_base), tolerance = 1e-8)
})

# ---- Group 5 ---------------------------------------------------------------

test_that("procrustes pool and signed_tucker_congruence pool agree on the signed-permutation orbit", {
  base_fit <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                  rotation = "none", se = "information")

  gauges <- list(
    diag(k),
    .signed_perm(c(2L, 1L, 3L), c(-1,  1,  1)),
    .signed_perm(c(3L, 2L, 1L), c( 1, -1,  1)),
    .signed_perm(c(1L, 3L, 2L), c( 1,  1, -1)),
    .signed_perm(c(2L, 3L, 1L), c(-1, -1,  1))
  )
  gauge_fits <- lapply(gauges, .gauge_fit, base_fit = base_fit)

  unrot_list <- lapply(gauge_fits, function(f) unclass(f$unrot_loadings))

  aligned_stc <- EFAtools:::.efa_pooled_align_unrotated_list(
    unrot_list, align_unrotated = "signed_tucker_congruence", return_meta = TRUE
  )
  pool_stc <- EFAtools:::.efa_pooled_analytic_pool(
    fits = gauge_fits,
    unrot_loadings_aligned = aligned_stc$loadings,
    align_meta = aligned_stc$meta,
    ci = 0.95,
    align_unrotated = "signed_tucker_congruence"
  )

  aligned_pr <- EFAtools:::.efa_pooled_align_unrotated_list(
    unrot_list, align_unrotated = "procrustes", return_meta = TRUE
  )
  pool_pr <- EFAtools:::.efa_pooled_analytic_pool(
    fits = gauge_fits,
    unrot_loadings_aligned = aligned_pr$loadings,
    align_meta = aligned_pr$meta,
    ci = 0.95,
    align_unrotated = "procrustes"
  )

  expect_equal(unname(pool_pr$SE$unrot_loadings),
               unname(pool_stc$SE$unrot_loadings), tolerance = 1e-9)
  expect_equal(unname(pool_pr$SE$uniquenesses),
               unname(pool_stc$SE$uniquenesses),   tolerance = 1e-12)
})

# ---- Group 6 ---------------------------------------------------------------

test_that("NA-filled vcov_unrot_loadings on any imputation aborts with classed condition", {
  base_fit <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                  rotation = "none", se = "information")

  gauges <- list(diag(k),
                 .signed_perm(c(2L, 1L, 3L), c(-1, 1, 1)),
                 .signed_perm(c(1L, 3L, 2L), c( 1, 1,-1)))
  gauge_fits <- lapply(gauges, .gauge_fit, base_fit = base_fit)

  # `EFA()` NA-fills the whole `vcov_unrot_loadings` matrix when the analytic
  # covariance is unreliable (Heywood / singular bordered information); mimic
  # that here. The dispatcher must abort rather than silently produce all-NA
  # pooled SEs (the propagation helper would otherwise NA-skip every row and
  # the fail-closed mask would blank every pooled element with no diagnostic).
  gauge_fits[[3]]$vcov_unrot_loadings[] <- NA_real_

  unrot_list <- lapply(gauge_fits, function(f) unclass(f$unrot_loadings))
  aligned <- EFAtools:::.efa_pooled_align_unrotated_list(
    unrot_list, align_unrotated = "procrustes", return_meta = TRUE
  )

  expect_error(
    EFAtools:::.efa_pooled_analytic_pool(
      fits = gauge_fits,
      unrot_loadings_aligned = aligned$loadings,
      align_meta = aligned$meta,
      ci = 0.95,
      align_unrotated = "procrustes"
    ),
    class = "efa_pooled_unreliable_vcov"
  )
})

test_that(".efa_pooled_propagate_procrustes_vcov NA-fills rows whose vcov block has NAs", {
  set.seed(20260622)
  p <- 4L
  kk <- 2L
  A <- matrix(stats::rnorm((p * kk)^2), p * kk, p * kk)
  V <- crossprod(A) / (p * kk)
  Q <- .givens(pi / 7, kk)

  i_bad <- 2L
  idx_bad <- (seq_len(kk) - 1L) * p + i_bad
  V[idx_bad, idx_bad[1]] <- NA_real_

  out <- EFAtools:::.efa_pooled_propagate_procrustes_vcov(V, Q, p, kk)
  expect_true(all(is.na(out[i_bad, ])))
  other_rows <- setdiff(seq_len(p), i_bad)
  expect_true(all(!is.na(out[other_rows, , drop = FALSE])))
})

# ---- Group 7 ---------------------------------------------------------------

test_that("missing vcov_unrot_loadings under align_unrotated = 'procrustes' aborts with classed condition", {
  base_fit <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                  rotation = "none", se = "information")

  gauges <- list(diag(k), .signed_perm(c(2L, 1L, 3L), c(-1, 1, 1)))
  gauge_fits <- lapply(gauges, .gauge_fit, base_fit = base_fit)
  gauge_fits[[2]]$vcov_unrot_loadings <- NULL

  unrot_list <- lapply(gauge_fits, function(f) unclass(f$unrot_loadings))
  aligned <- EFAtools:::.efa_pooled_align_unrotated_list(
    unrot_list, align_unrotated = "procrustes", return_meta = TRUE
  )

  expect_error(
    EFAtools:::.efa_pooled_analytic_pool(
      fits = gauge_fits,
      unrot_loadings_aligned = aligned$loadings,
      align_meta = aligned$meta,
      ci = 0.95,
      align_unrotated = "procrustes"
    ),
    class = "efa_pooled_no_vcov"
  )
})

test_that("procrustes pool returns correctly shaped slots", {
  pooled <- suppressMessages(EFA_POOLED(
    identical_list, n_factors = k, N = N_id, method = "ML",
    rotation = "none", se = "information",
    align_unrotated = "procrustes"
  ))

  loading_dn <- dimnames(unclass(pooled$unrot_loadings))

  expect_identical(dim(pooled$SE$unrot_loadings),         c(p_vars, k))
  expect_identical(dimnames(pooled$SE$unrot_loadings),    loading_dn)
  expect_identical(dim(pooled$CI$unrot_loadings$lower),   c(p_vars, k))
  expect_identical(dim(pooled$CI$unrot_loadings$upper),   c(p_vars, k))
  expect_identical(dimnames(pooled$CI$unrot_loadings$lower), loading_dn)

  expect_setequal(names(pooled$MI$unrot_loadings), c("RIV", "FMI", "df"))
  expect_identical(dim(pooled$MI$unrot_loadings$RIV), c(p_vars, k))
  expect_identical(dim(pooled$MI$unrot_loadings$FMI), c(p_vars, k))
  expect_identical(dim(pooled$MI$unrot_loadings$df),  c(p_vars, k))

  expect_length(pooled$SE$uniquenesses,             p_vars)
  expect_length(pooled$CI$uniquenesses$lower,       p_vars)
  expect_length(pooled$CI$uniquenesses$upper,       p_vars)
  expect_setequal(names(pooled$MI$uniquenesses),    c("RIV", "FMI", "df"))

  expect_null(pooled$replicates)
})
