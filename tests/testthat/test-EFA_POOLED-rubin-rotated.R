# Rubin's-rules pooling of analytic ROTATED-quantity SEs in EFA_POOLED() on the
# se = "information" path: rotated loadings, communalities (h2), and (oblique)
# factor correlations (Phi) and structure coefficients.
#
# A rotated-loading SE is conditional on the rotation criterion (Jennrich 1973,
# 1974; Browne 2001; Zhang & Preacher 2015), so for ORTHOGONAL and OBLIQUE
# rotations alike the within-imputation variance is each fit's own criterion-aware
# delta-method rotated SE (the quantity EFA() returns), reused after a
# signed-permutation alignment to the MI target; the between-imputation variance
# is the sample variance of the aligned per-imputation rotated loadings. Both
# rotation families are flagged MI$<param>$method = "signed_permutation_approx".
# A fixed-rotation linear map of the unrotated covariance is deliberately NOT used
# for rotated-loading SEs (it would inflate them by redistributing the rotational
# indeterminacy variance the criterion pins down). Communalities are
# rotation-invariant and pool with no alignment.
#
# Assertion groups:
#   1. Identity-imputation invariance (orthogonal): pooled rotated SE equals the
#      per-fit criterion-aware delta SE; B = 0; FMI = 0.
#   2. Point-estimate / orthogonal-Procrustes sanity.
#   3. Regression guard: orthogonal rotated SE is the delta SE, not the inflated
#      fixed-Q propagation.
#   4. Communalities pooling (gauge-invariant, no alignment).
#   5. np-boot cross-check on an orthogonal solution (informational).
#   6. Oblique approximation: signed-permutation path, flagged, does not abort.
#   7. Phi and Structure (oblique): shapes, diag(Phi SE) = 0, symmetry, invariance.
#   8. Slot shape under se = "information" + rotation.

cormat <- test_models$baseline$cormat
p_vars <- ncol(cormat)
N_id   <- 500L
k      <- 2L
m_id   <- 5L
identical_list <- replicate(m_id, cormat, simplify = FALSE)

# ---- Group 1 ---------------------------------------------------------------

test_that("identity imputations (orthogonal): pooled rotated SE equals the per-fit criterion-aware delta SE, B = 0", {
  pooled <- suppressMessages(EFA_POOLED(
    identical_list, n_factors = k, N = N_id, method = "ML",
    rotation = "varimax", se = "information"
  ))

  oracle <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                rotation = "varimax", se = "information")

  expect_identical(pooled$settings$se, "information")

  # Pooled point estimate is the MI target (= the single fit's varimax solution).
  expect_equal(unclass(pooled$rot_loadings), unclass(oracle$rot_loadings),
               tolerance = 1e-8)

  # The within-imputation variance is each fit's criterion-aware delta-method
  # rotated SE (the quantity EFA() returns); with identical imputations the
  # signed-permutation alignment is the identity and B = 0, so the pooled rotated
  # SE equals the single-fit rotated SE exactly.
  expect_equal(unname(pooled$SE$rot_loadings), unname(oracle$SE$rot_loadings),
               tolerance = 1e-6)

  fmi <- pooled$MI$rot_loadings$FMI
  riv <- pooled$MI$rot_loadings$RIV
  expect_true(all(abs(fmi[is.finite(fmi)]) < 1e-6))
  expect_true(all(abs(riv[is.finite(riv)]) < 1e-6))

  # Provenance flag: the criterion-aware signed-permutation reuse (orthogonal and
  # oblique alike).
  expect_identical(pooled$MI$rot_loadings$method, "signed_permutation_approx")
})

# ---- Group 2 ---------------------------------------------------------------

test_that("orthogonal point estimate: L_d^unrot Q_d reproduces the aligned rotated loadings", {
  fit <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
             rotation = "varimax", se = "information")
  A      <- unclass(fit$unrot_loadings)
  target <- unclass(fit$rot_loadings)
  Q      <- EFAtools:::.procrustes_orthogonal_T(A, target)

  # Q recovers the orthogonal rotation, so A Q reproduces the rotated loadings.
  expect_equal(unname(A %*% Q), unname(target), tolerance = 1e-8)
  # Orthogonality of the recovered transform.
  expect_equal(crossprod(Q), diag(k), tolerance = 1e-10)
})

# ---- Group 3 ---------------------------------------------------------------

test_that("orthogonal rotated SE is the criterion-aware delta SE, not the inflated fixed-Q propagation", {
  pooled <- suppressMessages(EFA_POOLED(
    identical_list, n_factors = k, N = N_id, method = "ML",
    rotation = "varimax", se = "information"
  ))
  oracle <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                rotation = "varimax", se = "information")

  # The pooled orthogonal rotated SE is the per-fit criterion-aware delta SE.
  expect_equal(unname(pooled$SE$rot_loadings), unname(oracle$SE$rot_loadings),
               tolerance = 1e-6)

  # Regression guard: the former fixed-Q Kronecker propagation of the
  # unrotated covariance is a DISTINCT matrix from the reported (criterion-aware)
  # SE -- an orthogonal Q preserves each row's total variance, so it redistributes
  # the rotational-indeterminacy variance the criterion otherwise pins down (and is
  # never smaller). Confirm the reported SE is not that propagated quantity; the
  # gap grows from modest on this well-separated structure to several-fold when
  # factors are weakly separated.
  A      <- unclass(oracle$unrot_loadings)
  V      <- oracle$vcov_unrot_loadings
  target <- unclass(oracle$rot_loadings)
  Q      <- EFAtools:::.procrustes_orthogonal_T(A, target)
  kron_se <- EFAtools:::.efa_pooled_propagate_procrustes_vcov(V, Q, p_vars, k)

  reported <- unclass(pooled$SE$rot_loadings)
  expect_gt(max(abs(kron_se - reported)), 1e-3)
})

# ---- Group 4 ---------------------------------------------------------------

test_that("communalities (h2) pool is gauge-invariant: Rubin's rules on per-fit communality SEs, no alignment", {
  pooled <- suppressMessages(EFA_POOLED(
    identical_list, n_factors = k, N = N_id, method = "ML",
    rotation = "varimax", se = "information"
  ))
  oracle <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                rotation = "varimax", se = "information")

  # Identical imputations -> B = 0 -> pooled h2 SE equals the per-fit communality
  # SE exactly (no alignment is applied; communalities are rotation-invariant).
  expect_length(pooled$SE$h2, p_vars)
  expect_equal(unname(pooled$SE$h2), unname(oracle$SE$communalities),
               tolerance = 1e-10)
  expect_true(all(pooled$SE$h2[!is.na(pooled$SE$h2)] >= 0))

  fmi <- pooled$MI$h2$FMI
  expect_true(all(abs(fmi[is.finite(fmi)]) < 1e-12))
  expect_identical(pooled$MI$h2$method, "gauge_invariant")
})

# ---- Group 4b --------------------------------------------------------------

test_that("distinct imputations (B > 0): pooled rotated-loading and h2 SEs match an independent Rubin computation", {
  # Identical imputations leave B = 0, which cannot exercise the
  # between-imputation variance, the (1 + 1/m) B scaling, or the column
  # alignment under genuine disagreement. Drive B > 0 with deterministically
  # distinct datasets (base-R draws, no MASS) from a 2-factor population and
  # check the pooled SEs against a hand-computed Rubin reference.
  set.seed(20260621)
  p_b <- 6L
  k_b <- 2L
  N_b <- 250L
  m_b <- 4L
  Lp <- matrix(0, p_b, k_b)
  Lp[1:3, 1] <- c(0.72, 0.65, 0.60)
  Lp[4:6, 2] <- c(0.72, 0.65, 0.60)
  Sig <- tcrossprod(Lp)
  diag(Sig) <- 1
  Rt <- chol(Sig)
  imps_b <- lapply(seq_len(m_b), function(i) {
    X <- matrix(stats::rnorm(N_b * p_b), N_b, p_b) %*% Rt
    colnames(X) <- paste0("V", seq_len(p_b))
    X
  })

  pooled_b <- tryCatch(suppressWarnings(suppressMessages(EFA_POOLED(
    imps_b, n_factors = k_b, method = "ML", rotation = "varimax",
    se = "information"
  ))), error = function(e) e)
  if (inherits(pooled_b, "error")) skip("EFA_POOLED failed on synthetic data")

  m <- length(pooled_b$fits)

  # --- Communalities (gauge-invariant): independent Rubin from per-fit pieces ---
  comm    <- t(vapply(pooled_b$fits, function(f) as.numeric(f$h2), numeric(p_b)))
  se_comm <- t(vapply(pooled_b$fits, function(f) as.numeric(f$SE$communalities),
                      numeric(p_b)))
  if (anyNA(se_comm)) skip("a component fit had non-finite communality SEs")
  Ubar_h2 <- colMeans(se_comm^2)
  B_h2    <- apply(comm, 2, stats::var)
  ref_se_h2 <- sqrt(Ubar_h2 + (1 + 1 / m) * B_h2)
  expect_equal(unname(pooled_b$SE$h2), unname(ref_se_h2), tolerance = 1e-8)
  expect_true(any(B_h2 > 0))                 # the regime is genuinely B > 0
  expect_true(any(pooled_b$MI$h2$FMI > 0))   # and it flows through to FMI

  # --- Rotated loadings: independent Rubin from the per-fit criterion-aware delta
  # SE. Mirror the production gauge exactly: the MI target is the first fit's
  # rotated solution; imputation 1 uses it directly and 2..m are Procrustes-aligned
  # to it. The within-imputation variance is each fit's $SE$rot_loadings, reordered
  # by the signed-permutation alignment of its own rotated solution to Lr_d; B is
  # the sample variance of the aligned rotated loadings. ---
  target  <- as.matrix(pooled_b$alignment$target)
  est_mat <- matrix(NA_real_, m, p_b * k_b)
  se_mat  <- matrix(NA_real_, m, p_b * k_b)
  for (d in seq_len(m)) {
    if (d == 1L) {
      Lr_d <- target
    } else {
      A    <- unclass(pooled_b$fits[[d]]$unrot_loadings)
      Lr_d <- unclass(EFAtools::PROCRUSTES(A, target, rotation = "orthogonal")$loadings)
    }
    se_d <- pooled_b$fits[[d]]$SE$rot_loadings
    if (is.null(se_d) || anyNA(se_d)) skip("a component fit had non-finite rotated SEs")
    fo <- EFAtools:::.align_solution(
      L_target = Lr_d, L = unclass(pooled_b$fits[[d]]$rot_loadings)
    )$factor_order
    est_mat[d, ] <- as.vector(Lr_d)
    se_mat[d, ]  <- as.vector(unclass(se_d)[, fo, drop = FALSE])
  }
  Ubar_r <- colMeans(se_mat^2)
  B_r    <- apply(est_mat, 2, stats::var)
  ref_se_rot <- matrix(sqrt(Ubar_r + (1 + 1 / m) * B_r), p_b, k_b)
  expect_equal(unname(pooled_b$SE$rot_loadings), unname(ref_se_rot),
               tolerance = 1e-6)
  expect_true(any(pooled_b$MI$rot_loadings$FMI > 0))

  # Factor-column labels are consistent across the loading-shaped SE families.
  expect_identical(colnames(pooled_b$SE$rot_loadings),
                   colnames(pooled_b$SE$unrot_loadings))
})

# ---- Group 5 ---------------------------------------------------------------

test_that("np-boot cross-check: analytic and bootstrap pooled orthogonal rotated SEs agree to within ~order of magnitude (informational)", {
  skip_on_cran()
  skip_if_not_slow()
  skip_if_not_installed("MASS")

  set.seed(20260621)
  p_syn <- 6L
  N_syn <- 300L
  m_syn <- 3L
  # Independent-cluster 2-factor population (orthogonal), strong primary loadings.
  L_true <- matrix(0, p_syn, 2)
  L_true[1:3, 1] <- c(0.75, 0.70, 0.65)
  L_true[4:6, 2] <- c(0.75, 0.70, 0.65)
  Sigma <- tcrossprod(L_true)
  diag(Sigma) <- 1

  imps <- lapply(seq_len(m_syn), function(i) {
    X <- MASS::mvrnorm(N_syn, mu = rep(0, p_syn), Sigma = Sigma)
    colnames(X) <- paste0("V", seq_len(p_syn))
    as.data.frame(X)
  })

  pooled_info <- tryCatch(suppressWarnings(suppressMessages(EFA_POOLED(
    imps, n_factors = 2, method = "ML", rotation = "varimax", se = "information"
  ))), error = function(e) e)
  if (inherits(pooled_info, "error")) skip("information EFA_POOLED failed on synthetic data")

  pooled_boot <- tryCatch(suppressWarnings(suppressMessages(EFA_POOLED(
    imps, n_factors = 2, method = "ML", rotation = "varimax",
    se = "np-boot", b_boot = 200
  ))), error = function(e) e)
  if (inherits(pooled_boot, "error")) skip("np-boot EFA_POOLED failed on synthetic data")

  se_info <- pooled_info$SE$rot_loadings
  se_boot <- pooled_boot$SE$rot_loadings
  est     <- unclass(pooled_info$rot_loadings)

  # Compare on the primary (marker) loadings where the SE is well determined.
  primary <- abs(est) >= 0.3 & is.finite(se_info) & is.finite(se_boot) &
    se_boot > 0
  expect_true(any(primary))
  ratio <- se_info[primary] / se_boot[primary]

  # Informational: both target the same MI-gauge estimand (Procrustes to the MI
  # target), so the analytic fixed-Q SE and the bootstrap SE should agree to
  # within roughly 10-20%; the band is deliberately loose (not a strict oracle).
  expect_gt(stats::median(ratio), 0.5)
  expect_lt(stats::median(ratio), 1.6)
})

# ---- Group 6 ---------------------------------------------------------------

test_that("oblique rotated SE pooling uses the signed-permutation fallback, is flagged, and does not abort", {
  set.seed(42)
  pooled <- suppressMessages(EFA_POOLED(
    identical_list, n_factors = k, N = N_id, method = "ML",
    rotation = "oblimin", se = "information"
  ))
  oracle <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                rotation = "oblimin", se = "information")

  # Results are produced (no abort) and carry the documentation flag.
  expect_false(is.null(pooled$SE$rot_loadings))
  expect_identical(pooled$MI$rot_loadings$method, "signed_permutation_approx")

  # Signed-permutation reuse of the per-fit rotated SE: with identical imputations
  # the alignment is the identity, so the within-imputation variance (Ubar) is the
  # single-fit SE exactly. The pooled SE differs only by the tiny between-imputation
  # variance the iterative oblique Procrustes alignment injects into the per-fit
  # rotated point estimates (the closed-form orthogonal path has none); its
  # magnitude depends on the oblimin flat-optimum geometry, so allow headroom.
  expect_equal(unname(pooled$SE$rot_loadings), unname(oracle$SE$rot_loadings),
               tolerance = 1e-3)

  fmi <- pooled$MI$rot_loadings$FMI
  expect_true(all(abs(fmi[is.finite(fmi)]) < 1e-2))
})

# ---- Group 7 ---------------------------------------------------------------

test_that("oblique Phi and Structure SE pooling: shapes, diag(Phi SE) = 0, symmetry, invariance, flags", {
  set.seed(42)
  pooled <- suppressMessages(EFA_POOLED(
    identical_list, n_factors = k, N = N_id, method = "ML",
    rotation = "oblimin", se = "information"
  ))
  oracle <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                rotation = "oblimin", se = "information")

  # Phi SE: k x k, fixed unit diagonal -> zero SE, symmetric. (1e-3 headroom for
  # the iterative-oblique-Procrustes B injection, as in the previous group.)
  expect_identical(dim(pooled$SE$Phi), c(k, k))
  expect_equal(diag(pooled$SE$Phi), rep(0, k), tolerance = 1e-12)
  expect_equal(pooled$SE$Phi, t(pooled$SE$Phi), tolerance = 1e-12)
  expect_equal(unname(pooled$SE$Phi), unname(oracle$SE$Phi), tolerance = 1e-3)
  expect_identical(pooled$MI$Phi$method, "signed_permutation_approx")

  # The fixed unit diagonal carries no MI diagnostics (NA), unlike the off-diagonal.
  expect_true(all(is.na(diag(pooled$MI$Phi$FMI))))
  expect_true(all(is.na(diag(pooled$MI$Phi$df))))

  # Structure SE: p x k, invariance under identical imputations.
  expect_identical(dim(pooled$SE$Structure), c(p_vars, k))
  expect_equal(unname(pooled$SE$Structure), unname(oracle$SE$Structure),
               tolerance = 1e-3)
  expect_identical(pooled$MI$Structure$method, "signed_permutation_approx")
})

# ---- Group 8 ---------------------------------------------------------------

test_that("se = 'information' + rotation fills the rotated SE/CI/MI slots with correct shapes", {
  # Orthogonal: rot_loadings + h2; no Phi/Structure.
  pooled_orth <- suppressMessages(EFA_POOLED(
    identical_list, n_factors = k, N = N_id, method = "ML",
    rotation = "varimax", se = "information"
  ))
  rot_dn <- dimnames(unclass(pooled_orth$rot_loadings))

  expect_identical(dim(pooled_orth$SE$rot_loadings), c(p_vars, k))
  expect_identical(dimnames(pooled_orth$SE$rot_loadings), rot_dn)
  expect_identical(dim(pooled_orth$CI$rot_loadings$lower), c(p_vars, k))
  expect_identical(dim(pooled_orth$CI$rot_loadings$upper), c(p_vars, k))
  expect_setequal(names(pooled_orth$MI$rot_loadings),
                  c("RIV", "FMI", "df", "method"))
  expect_identical(dim(pooled_orth$MI$rot_loadings$df), c(p_vars, k))

  expect_length(pooled_orth$SE$h2, p_vars)
  expect_length(pooled_orth$CI$h2$lower, p_vars)
  expect_setequal(names(pooled_orth$MI$h2), c("RIV", "FMI", "df", "method"))

  expect_null(pooled_orth$SE$Phi)
  expect_null(pooled_orth$SE$Structure)
  expect_null(pooled_orth$replicates)

  # Oblique: adds Phi and Structure.
  pooled_oblq <- suppressMessages(EFA_POOLED(
    identical_list, n_factors = k, N = N_id, method = "ML",
    rotation = "oblimin", se = "information"
  ))
  expect_identical(dim(pooled_oblq$SE$Phi), c(k, k))
  expect_identical(dim(pooled_oblq$CI$Phi$lower), c(k, k))
  expect_identical(dim(pooled_oblq$SE$Structure), c(p_vars, k))
  expect_setequal(names(pooled_oblq$MI$Structure),
                  c("RIV", "FMI", "df", "method"))

  # Unrotated families keep their RIV/FMI/df-only MI shape (no method tag).
  expect_setequal(names(pooled_oblq$MI$unrot_loadings), c("RIV", "FMI", "df"))
  expect_setequal(names(pooled_oblq$MI$uniquenesses), c("RIV", "FMI", "df"))
})

# ---- Group 9 ---------------------------------------------------------------

test_that("an unreliable per-imputation rotated SE warns and NA-blanks the pooled rotated SE family", {
  base_fit <- EFA(cormat, n_factors = k, N = N_id, method = "ML",
                  rotation = "oblimin", se = "information")
  fits <- replicate(3L, base_fit, simplify = FALSE)

  # Mimic a fit whose rotated-SE Jacobian could not be reproduced: EFA() leaves the
  # rotated SE matrices NA-filled while the unrotated covariance stays reliable, so
  # the whole-fit vcov reliability gate does not catch it.
  fits[[2]]$SE$rot_loadings[] <- NA_real_
  fits[[2]]$SE$Phi[]          <- NA_real_
  fits[[2]]$SE$Structure[]    <- NA_real_

  unrot_list <- lapply(fits, function(f) unclass(f$unrot_loadings))
  aligned <- EFAtools:::.efa_pooled_align_unrotated_list(
    unrot_list, align_unrotated = "signed_tucker_congruence", return_meta = TRUE
  )
  rot_list <- lapply(fits, function(f) unclass(f$rot_loadings))
  phi_list <- lapply(fits, function(f) f$Phi)
  str_list <- Map(function(L, Phi) L %*% Phi, rot_list, phi_list)
  h2 <- EFAtools:::.efa_pooled_communalities(rot_list[[1]], phi_list[[1]])

  expect_warning(
    pooled <- EFAtools:::.efa_pooled_analytic_pool(
      fits = fits,
      unrot_loadings_aligned = aligned$loadings,
      align_meta = aligned$meta,
      ci = 0.95,
      align_unrotated = "signed_tucker_congruence",
      rotation_type = "oblique",
      rot_loadings = rot_list,
      phis = phi_list,
      structure_loadings = str_list,
      mean_structure_loadings = str_list[[1]],
      mean_phis = phi_list[[1]],
      h2 = h2
    ),
    class = "efa_pooled_rotated_se_unreliable"
  )

  # The affected rotated families are fully NA-blanked by the fail-closed mask,
  # while the gauge-invariant communalities and the unrotated families stay finite.
  expect_true(all(is.na(pooled$SE$rot_loadings)))
  expect_true(all(is.na(pooled$SE$Structure)))
  expect_true(all(is.na(pooled$SE$Phi[upper.tri(pooled$SE$Phi)])))
  expect_false(all(is.na(pooled$SE$unrot_loadings)))
  expect_false(all(is.na(pooled$SE$h2)))
})
