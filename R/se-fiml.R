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
# (2002). The construction shares its geometry with the ordinal/continuous robust sandwich -- both
# take D, H, the bordered bread and the gauge constraint from `.offdiag_sandwich_pieces()` -- and
# differs only in the meat, which here is the FIML saturated covariance.

# Robust two-stage loading covariance V_AA (p*k x p*k), the unrotated loading/uniqueness SEs, and
# the reliability flags, from the fitted loadings and the FIML moments. The model Jacobian, the
# Stage-2 weight (the normal-theory ML weight for ML, the identity for ULS), the bordered bread and
# the gauge constraint come from `.offdiag_sandwich_pieces()`, shared with the polychoric/ADF and
# expected-information sandwiches; only the meat differs. Returns NA SEs (reliable = FALSE) at a
# singular bordered system, a degenerate saturated covariance, or a non-PSD V_AA.
.se_fiml_core <- function(fit_out, fiml, N, method) {

  L <- unclass(fit_out$unrot_loadings)
  p <- nrow(L)
  k <- ncol(L)
  pk <- p * k

  # `gauge_reliable = TRUE` on the NA path, as on the other analytic paths: everything is already
  # withheld and the gauge is not what withheld it, so the caller's gauge-specific message must not
  # be attached to this branch.
  na_core <- list(
    V_AA = matrix(NA_real_, pk, pk),
    loadings_se = matrix(NA_real_, p, k),
    uniquenesses_se = rep(NA_real_, p),
    reliable = FALSE,
    gauge_reliable = TRUE
  )

  pieces <- .offdiag_sandwich_pieces(L, method)
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

  # Two-stage sandwich V = A^- (V D)' Omega_delta (V D) A^- / (N - 1). Omega_delta is the covariance
  # of the saturated estimates (the inverse observed information, so on the variance scale
  # Var(r-hat)); it is put on the unit asymptotic-variance scale by N and the covariance divided by
  # N - 1, exactly as the polychoric/ADF sandwich and the expected information do. That N - 1 is the
  # package-wide small-sample convention (it is also the scale `.gof()` and `.fiml_scaled_test()`
  # work on), and it makes the complete-data limit of this sandwich agree with the Pearson path
  # instead of differing from it by a systematic factor of N / (N - 1).
  VD <- pieces$VD
  Abread <- pieces$Abread
  meat <- crossprod(VD, (N * Omega) %*% VD)
  V_AA <- Abread %*% meat %*% Abread / (N - 1)
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
    reliable = reliable,
    # Conditioning of the rotational gauge the unrotated loadings are reported in, from the shared
    # geometry. Reported here so a degenerate orientation is caught on this path exactly as on the
    # Pearson/polychoric one, rather than shipping an inflated unrotated loading SE as informative.
    gauge_reliable = pieces$gauge_reliable
  )
}

# Corrected two-stage SE/CI for cor_method = "fiml" (se = "information" / "sandwich" both route
# here -- the naive Stage-2 information SE is never shipped). Fills the standard analytic SE/CI
# schema. The unrotated branch reuses the robust unrotated wrapper; the rotated branch reuses the
# information-SE rotation propagation by supplying the robust covariance through its `se0` contract
# (so promax/simplimax stay bootstrap-only, exactly as elsewhere). NA-fills with the classed
# efa_se_unreliable warning at a Heywood case / non-PSD covariance.
.se_fiml <- function(fit_out, rot_info, N, ci, fiml, method) {

  core <- .se_fiml_core(fit_out, fiml, N, method)

  # The two-stage sandwich withholds at a boundary solution on the same test, and through the same
  # shared policy, as every other analytic path. The corrected chi-square is untouched: it is built
  # in `.gof()` rather than here, and is a discrepancy-function quantity the boundary does not
  # invalidate.
  core <- .withhold_at_boundary(core, fit_out$unrot_loadings)

  if (is.null(rot_info)) {
    .se_sandwich_unrotated(fit_out, core, ci)
  } else {
    se0 <- list(vcov = core$V_AA,
                loadings_se = core$loadings_se,
                uniquenesses_se = core$uniquenesses_se,
                gauge_reliable = core$gauge_reliable)
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

  # The rescaled statistic is gauge-invariant, so the gauge diagnostic is not read here; skip it.
  # This runs once per bootstrap replicate, where the discarded decomposition would be pure cost.
  pieces <- .offdiag_sandwich_pieces(L, method, gauge = FALSE)
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
