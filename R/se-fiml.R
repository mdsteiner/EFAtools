# Corrected two-stage standard errors and the corrected (Satorra-Bentler) two-stage chi-square for
# cor_method = "fiml". The factor model is fitted to the EM-estimated saturated correlation R
# (Stage 1), so the naive Stage-2 standard errors -- which treat R as if it were complete raw data
# -- are inconsistent under missingness. The corrected covariance is the two-stage sandwich
#   Omega_tilde = (D' H D)^- D' H Omega_delta H D (D' H D)^-,
# with D = d sigma_offdiag / d vec(Lambda) the model derivative, H the Stage-2 weight at the
# solution (the normal-theory ML weight for ML, the identity for ULS), and Omega_delta the
# correlation-scale asymptotic covariance of the saturated FIML estimates
# (.fiml_saturated_acov()$cor). Likewise the plain two-stage
# likelihood-ratio chi-square is not asymptotically chi-square(df) under the two-stage estimator;
# it is rescaled by the Satorra-Bentler correction built from the same H and Omega_delta.
# Two-stage standard errors / rescaled statistic: Yuan & Bentler (2000); Savalei & Bentler (2009,
# SEM); Yuan, Marshall & Bentler (2002, Psychometrika 67:95-121). MAR/ignorability: Little & Rubin
# (2002). The construction mirrors the ordinal/continuous robust sandwich (.se_sandwich_core), with
# H fixed to the normal-theory ML weight and the meat supplied by the FIML saturated covariance.

# Shared building blocks of the two-stage sandwich and the rescaled chi-square: the off-diagonal
# model Jacobian D, the Stage-2 weight V matched to the extraction method (the normal-theory ML
# weight for ML, the identity for ULS; as in .se_sandwich_core), and the bordered (constrained
# generalised) inverse of the bread A = D' V D. `method` also selects the gauge the unrotated
# loadings are reported in (Lambda' Psi^-1 Lambda for ML, Lambda' Lambda otherwise; Lawley &
# Maxwell); the projector D' V D and the rescaled chi-square are gauge-invariant. Returns NULL when
# the model-implied correlation or the bordered system is singular.
.fiml_sandwich_pieces <- function(L, method) {
  L <- unclass(L)
  p <- nrow(L)
  k <- ncol(L)
  pk <- p * k

  pairs <- utils::combn(p, 2L)
  pi <- pairs[1, ]
  pj <- pairs[2, ]
  n <- ncol(pairs)

  # Model Jacobian Delta = d sigma_ij / d vec(Lambda) (n x pk). For pair (i, j),
  # d sigma_ij / d Lambda[a, f] = (a == i) Lambda[j, f] + (a == j) Lambda[i, f]; rows in
  # utils::combn(p, 2) order, columns in column-major vec(Lambda) order (as in .se_sandwich_core).
  Delta <- matrix(0, n, pk)
  for (f in seq_len(k)) {
    Delta[cbind(seq_len(n), (f - 1L) * p + pi)] <- L[pj, f]
    Delta[cbind(seq_len(n), (f - 1L) * p + pj)] <- L[pi, f]
  }

  # Stage-2 weight V, matched to the extraction method exactly as .se_sandwich_core() does: ULS
  # minimises the unweighted off-diagonal residual sum of squares (identity weight); ML uses the
  # normal-theory GLS weight V = 1/2 (Sigma^-1 (x) Sigma^-1) restricted to the off-diagonal pairs,
  # at the model-implied correlation Sigma = Lambda Lambda' (unit diagonal). Held as a diagonal
  # (vdiag) or a full matrix (Vmat), the other NULL, so the meat and the rescaled chi-square apply
  # the same weight the loadings were fitted with.
  if (method == "ULS") {
    vdiag <- rep(1, n)
    Vmat <- NULL
  } else {
    Sigma <- tcrossprod(L)
    diag(Sigma) <- 1
    P <- tryCatch(solve(Sigma), error = function(e) NULL)
    if (is.null(P)) return(NULL)
    vdiag <- NULL
    Vmat <- matrix(0, n, n)
    for (s in seq_len(n)) {
      a <- pi[s]
      b <- pj[s]
      Vmat[, s] <- 0.5 * (P[pi, a] * P[pj, b] + P[pi, b] * P[pj, a])
    }
  }
  is_diag <- is.null(Vmat)

  # Bread A = D' V D (singular by the k(k-1)/2 rotational gauge freedoms); border with the
  # gauge-fixing constraint matched to the solution's orientation and take the leading pk-block of
  # the augmented inverse (the reflexive generalised inverse of A).
  VD <- if (is_diag) vdiag * Delta else Vmat %*% Delta
  A <- crossprod(Delta, VD)
  psi <- 1 - rowSums(L^2)
  use_ltpil <- k >= 2L && !anyNA(psi) && all(psi > 0) && method == "ML"
  Cmat <- if (use_ltpil) .se_sandwich_constraint(L, psi = psi) else .se_sandwich_constraint(L)
  Abread <- .bordered_inverse_block(A, Cmat, pk)
  if (is.null(Abread)) return(NULL)

  list(Delta = Delta, vdiag = vdiag, Vmat = Vmat, VD = VD, Abread = Abread, pairs = pairs)
}

# Robust two-stage loading covariance V_AA (p*k x p*k), the unrotated loading/uniqueness SEs, and a
# reliability flag, from the fitted loadings and the FIML moments. Returns NA SEs (reliable = FALSE)
# at a singular bordered system, a degenerate saturated covariance, or a non-PSD V_AA.
.se_fiml_core <- function(fit_out, fiml, method) {

  L <- unclass(fit_out$unrot_loadings)
  p <- nrow(L)
  k <- ncol(L)
  pk <- p * k

  na_core <- list(
    V_AA = matrix(NA_real_, pk, pk),
    loadings_se = matrix(NA_real_, p, k),
    uniquenesses_se = rep(NA_real_, p),
    reliable = FALSE
  )

  pieces <- .fiml_sandwich_pieces(L, method)
  if (is.null(pieces)) return(na_core)

  # Stage-1 saturated correlation asymptotic covariance Omega_delta (the meat), on the variance
  # scale Var(r-hat) and in utils::combn() order. Reuse the cached `acov_cor` when EFA() already
  # built it for the scaled chi-square (the point-estimate analytic-SE path); else build it here.
  Omega <- fiml$acov_cor
  if (is.null(Omega)) {
    Omega <- tryCatch(.fiml_saturated_acov(fiml$data, fiml$mu, fiml$sigma)$cor,
                      error = function(e) NULL)
  }
  if (is.null(Omega) || anyNA(Omega) || nrow(Omega) != p * (p - 1L) / 2L) return(na_core)

  # Two-stage sandwich V = A^- (V D)' Omega_delta (V D) A^-. Omega_delta is already the covariance of
  # the saturated estimates (inverse observed information), so it enters directly -- no N scaling
  # (contrast the polychoric/ADF sandwich, whose Gamma enters on the unit asymptotic-variance scale).
  VD <- pieces$VD
  Abread <- pieces$Abread
  meat <- crossprod(VD, Omega %*% VD)
  V_AA <- Abread %*% meat %*% Abread
  V_AA <- (V_AA + t(V_AA)) / 2

  # The covariance must be positive semidefinite; a non-finite entry or negative eigenvalue would
  # otherwise corrupt the marginal, uniqueness, and rotated SEs that read the full V_AA.
  reliable <- .is_psd(V_AA)
  loadings_se <- if (reliable) sqrt(diag(V_AA)) else rep(NA_real_, pk)
  # Uniqueness SE = communality SE (psi_i = 1 - rowSums(Lambda^2)_i), via the shared gradient.
  uniq_se <- if (!reliable) rep(NA_real_, p) else .communality_se(L, V_AA)

  list(
    V_AA = V_AA,
    loadings_se = matrix(loadings_se, p, k),
    uniquenesses_se = uniq_se,
    reliable = reliable
  )
}

# Corrected two-stage SE/CI for cor_method = "fiml" (se = "information" / "sandwich" both route
# here -- the naive Stage-2 information SE is never shipped). Fills the standard analytic SE/CI
# schema. The unrotated branch reuses the robust unrotated wrapper; the rotated branch reuses the
# information-SE rotation propagation by supplying the robust covariance through its `se0` contract
# (so promax/simplimax stay bootstrap-only, exactly as elsewhere). NA-fills with the classed
# efa_se_unreliable warning at a Heywood case / non-PSD covariance.
.se_fiml <- function(fit_out, rot_info, N, ci, fiml, method) {

  core <- .se_fiml_core(fit_out, fiml, method)

  if (is.null(rot_info)) {
    .se_sandwich_unrotated(fit_out, core, ci)
  } else {
    se0 <- list(vcov = core$V_AA,
                loadings_se = core$loadings_se,
                uniquenesses_se = core$uniquenesses_se)
    .se_information_rotated(fit_out, rot_info, N, ci, se0 = se0)
  }
}

# Satorra-Bentler-corrected two-stage chi-square block for cor_method = "fiml", in the schema
# .apply_scaled_test() consumes. The factor model is fitted to the saturated EM correlation R, so
# the corrected statistic is the residual-based two-stage statistic on the off-diagonal correlation
# structure: the normal-theory ML quadratic form rescaled by the saturated FIML covariance
# Omega_delta (the meat) -- exactly the .scaled_chisq() machinery the robust (sandwich) path uses,
# with the FIML covariance in place of the polychoric/ADF Gamma. Returns NULL for PAF (no
# discrepancy), a just/under-identified model, a missing N, or a degenerate/non-PD saturated
# covariance, where the plain likelihood-ratio / NA fallback is reported instead.
.fiml_scaled_test <- function(L, R, N, method, df, m, fiml) {

  if (!(method %in% c("ML", "ULS")) || is.na(N) ||
      is.null(df) || is.na(df) || df <= 0) {
    return(NULL)
  }

  L <- unclass(L)
  p <- nrow(L)

  pieces <- .fiml_sandwich_pieces(L, method)
  if (is.null(pieces)) return(NULL)

  # Reuse the cached saturated correlation covariance when EFA() pre-built it for the SE sandwich
  # (the point-estimate analytic-SE path); else build it here (every bootstrap replicate, and the
  # se = "none"/"np-boot" point estimate, omit the cache and land in this branch).
  Omega <- fiml$acov_cor
  if (is.null(Omega)) {
    Omega <- tryCatch(.fiml_saturated_acov(fiml$data, fiml$mu, fiml$sigma)$cor,
                      error = function(e) NULL)
  }
  if (is.null(Omega) || anyNA(Omega) || nrow(Omega) != p * (p - 1L) / 2L) return(NULL)

  # .scaled_chisq() takes Gamma on the unit asymptotic-variance scale (N * Var(rho-hat); = lavaan's
  # NACOV), so scale the variance-scale FIML covariance by N, exactly as the polychoric/ADF path
  # does. The minimal fit_out shape supplies the loadings, the analysis correlation, and df.
  fo <- list(unrot_loadings = L, orig_R = R, fit_indices = list(df = df))
  .scaled_chisq(fo, N * Omega, pieces$pairs, pieces$VD, pieces$vdiag, pieces$Vmat,
                pieces$Abread, N)
}
