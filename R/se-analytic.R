# Dispatch standard-error/confidence-interval computation on the requested method. The
# nonparametric bootstrap aggregates resampled fits; the analytic methods derive SEs from
# the fitted model itself. Every branch returns the same list(SE, CI, replicates) schema, so
# the print and summary methods are agnostic to how the SEs were produced. The sandwich branch
# additionally returns a `scaled_test` element (the robust scaled chi-square), which efa_fit()
# strips off and folds into the fit indices. `gamma` (the polychoric ACOV meat) and `method`
# are used only by the sandwich. `fiml` (the two-stage EM moments + raw data) is non-NULL only
# for cor_method = "fiml", where both analytic settings route to the corrected two-stage sandwich.
.compute_se_ci <- function(fit_out, L_rot, se_method, boot_fits = NULL,
                           boot_rot = "none", ci = .95, b = NULL, N = NULL,
                           rot_info = NULL, gamma = NULL, method = NULL,
                           fiml = NULL) {

  # cor_method = "fiml": the EM correlation is a two-stage estimate, so both analytic settings
  # ("information"/"sandwich") return the corrected two-stage sandwich SE built on the saturated
  # FIML covariance. The naive Stage-2 information SE (treating the EM correlation as complete data)
  # is inconsistent under missingness and is never shipped. The bootstrap path is unchanged.
  if (!is.null(fiml) && se_method %in% c("information", "sandwich")) {
    return(.se_fiml(fit_out, rot_info, N, ci, fiml, method))
  }

  switch(se_method,
    "np-boot" = .boot_se_ci(fit_out, L_rot, boot_fits, boot_rot, ci, b),
    # rot_info is non-NULL exactly for an analytic-SE fit under a real rotation; the unrotated
    # path (rotation = "none") leaves it NULL. Both are handled inside .se_information().
    "information" = .se_information(fit_out, rot_info, N, ci, method),
    "sandwich" = .se_sandwich_dispatch(fit_out, rot_info, N, ci, gamma, method),
    NULL
  )
}


# Leading `keep`-block of the inverse of the bordered matrix [A C'; C 0] -- the constrained
# (reflexive generalised) inverse of a gauge-singular information matrix A under the identification
# constraints whose Jacobian is `Cmat` (k(k-1)/2 rows, zero rows for a single factor). Returns NULL
# if the augmented system is singular. Shared by the expected-information and sandwich SEs.
.bordered_inverse_block <- function(A, Cmat, keep) {
  nc <- nrow(Cmat)
  Aug <- if (nc > 0L) {
    rbind(cbind(A, t(Cmat)), cbind(Cmat, matrix(0, nc, nc)))
  } else {
    A
  }
  inv <- tryCatch(solve(Aug), error = function(e) NULL)
  if (is.null(inv)) return(NULL)
  inv[seq_len(keep), seq_len(keep), drop = FALSE]
}

# Delta-method SEs of the communalities h2_i = rowSums(Lambda^2)_i from a loading covariance V
# (p*k x p*k over column-major vec(Lambda)). The gradient of h2_i is 2 Lambda[i, ], nonzero only in
# variable i's loading columns. Shared by every analytic path: psi_i = 1 - h2_i is a function of
# Lambda, not a free parameter, so this single route supplies the coinciding uniqueness and
# communality SEs from whichever loading covariance the path built.
.communality_se <- function(L, V) {
  p <- nrow(L)
  k <- ncol(L)
  G_h <- matrix(0, p, p * k)
  for (i in seq_len(p)) G_h[i, (seq_len(k) - 1L) * p + i] <- 2 * L[i, ]
  sqrt(pmax(rowSums((G_h %*% V) * G_h), 0))
}

# Skeleton gradient of the rotational identification constraint for the factor pair (u, v): a
# p x k matrix with column u set to X[, v] and column v to X[, u]. The off-diagonal (u, v) entry
# of X' Lambda has exactly this gradient in vec(Lambda) -- with X = Lambda it is the gradient of
# off-diag(Lambda' Lambda), and with X = Psi^-1 Lambda the direct part of off-diag(Lambda' Psi^-1
# Lambda); the caller adds the psi(Lambda) chain-rule term. Used by `.se_sandwich_constraint()`,
# which builds the gauge constraint for every analytic path.
.gauge_grad <- function(X, u, v) {
  g <- matrix(0, nrow(X), ncol(X))
  g[, u] <- X[, v]
  g[, v] <- X[, u]
  g
}

# TRUE if M is a (numerically) positive-semidefinite covariance: finite, with a non-negative
# diagonal and no eigenvalue below a small tolerance scaled by its largest diagonal entry. Used to
# gate the analytic parameter covariances before their square roots are reported: checking only the
# diagonal would pass a covariance that is non-PSD off the diagonal (where `pmax(., 0)` would then
# silently floor a negative rotated variance to an understated SE rather than flag it).
.is_psd <- function(M) {
  if (!all(is.finite(M))) return(FALSE)
  d <- diag(M)
  if (any(d < 0)) return(FALSE)
  sym <- (M + t(M)) / 2
  ev <- tryCatch(min(eigen(sym, symmetric = TRUE, only.values = TRUE)$values),
                 error = function(e) NA_real_)
  is.finite(ev) && ev >= -1e-8 * max(d, 1)
}

# TRUE when the solution sits on the parameter-space boundary: a uniqueness at (or below) its lower
# boundary, i.e. the Heywood case `.finalize_fit()` flags. The Wald approximation every analytic
# standard error rests on is not valid there, so this is the single gate on which the analytic paths
# withhold. The ML and ULS fitters constrain the uniquenesses to [.uniqueness_floor, 1] (DWLS does
# not -- see R/DWLS.R -- and can land below zero), so an improper solution from those fitters is
# pinned AT that floor and never reaches zero: a gate keyed on `psi <= 0` would fire only for a
# hand-built covariance and never for a fitted one. Test the boundary the fitters can actually
# reach.
.at_uniqueness_boundary <- function(L) {
  psi <- 1 - rowSums(unclass(L)^2)
  anyNA(psi) || any(psi <= .uniqueness_floor + sqrt(.Machine$double.eps))
}

# The hint the analytic paths attach to `efa_se_unreliable` when the boundary above is what withheld
# the standard errors. Held next to the predicate, and shared by the rotated and unrotated message
# sites, so the wording cannot drift between them. It states a total withholding, so a caller must
# only reach for it once it has established that nothing was reported.
.se_boundary_hint <- paste0(
  "The solution is a Heywood case -- a uniqueness sits at its lower boundary -- where the Wald ",
  "approximation these standard errors rest on is not valid, so none is reported. Use ",
  "{.code se = \"np-boot\"}, or extract fewer factors."
)

# Row/column labels for a vec(Lambda)-ordered loading covariance: "<variable>_<factor>" in the same
# column-major order the block is assembled in, so the documented ordering can be read off the
# object. NULL when the loadings carry no dimnames, in which case the block ships unlabelled rather
# than mislabelled.
.vec_loading_labels <- function(L) {
  rn <- rownames(L)
  cn <- colnames(L)
  if (is.null(rn) || is.null(cn)) return(NULL)
  as.vector(outer(rn, cn, paste, sep = "_"))
}


# Analytic expected-information SEs and Wald CIs for the ML solution, rotated or unrotated.
#
# EFAtools fits a CORRELATION structure: the analysed matrix has a diagonal fixed at exactly 1, and
# the uniquenesses psi = 1 - rowSums(Lambda^2) are a function of Lambda rather than free parameters.
# The Fisher information is therefore the one for theta = vec(Lambda) over the off-diagonal
# correlations alone,
#
#   I(theta) = Delta' Gamma^{-1} Delta,
#
# with Delta = d sigma_offdiag / d vec(Lambda) and Gamma the normal-theory asymptotic covariance of
# the sample correlations evaluated at the model-implied Sigma (Cudeck, 1989; Ogasawara, 1998).
# Assembling the covariance-structure information over (vec(Lambda), psi) instead would attribute
# Wishart sampling variability to a diagonal that carries none, and leak it into the loading block.
#
# I(theta) is exactly the Godambe sandwich of `.se_sandwich_core()` at the optimal (efficient)
# weight V = Gamma^{-1}, where the meat Delta' V Gamma V Delta collapses back onto the bread
# Delta' V Delta. Routing through the sandwich therefore reuses the model Jacobian, the gauge
# constraint, the bordered inverse, the PSD gate and the rotation propagation rather than
# duplicating them. No scaled chi-square is requested: the normal-theory weight is the one the ML
# discrepancy already minimises, so the ordinary chi-square stands.
#
# Cost: Gamma is n x n over the n = p(p-1)/2 pairs and is inverted once, so this is O(n^3) rather
# than the O((p k + p)^3) of a covariance-structure assembly -- a one-shot 0.1 s at p = 18 and under
# a second at p = 40, far below a bootstrap.
.se_information <- function(fit_out, rot_info, N, ci, method) {

  L <- unclass(fit_out$unrot_loadings)

  # At a boundary solution (`.at_uniqueness_boundary()`) the Wald intervals this path reports are
  # not valid, so withhold the covariance: the core then NA-fills through its usual unusable-Gamma
  # branch and the shared efa_se_unreliable warning fires, rather than reporting a boundary standard
  # error. Short-circuiting here also skips assembling and inverting an n x n Gamma the fit would
  # never use; `.se_sandwich_dispatch()` applies the same gate for the paths that still need the
  # core (it builds the scaled chi-square, which stays reportable at a boundary fit).
  Gamma <- if (.at_uniqueness_boundary(L)) {
    NULL
  } else {
    # Model-implied correlation matrix (unit diagonal by construction).
    Sigma <- tcrossprod(L)
    diag(Sigma) <- 1
    # `.se_sandwich_core()` takes Gamma on the variance scale Var(rho-hat) and restores the unit
    # scale itself, so undo the unit scaling `.normal_theory_gamma()` returns.
    .normal_theory_gamma(Sigma, utils::combn(nrow(L), 2L)) / N
  }

  .se_sandwich_dispatch(fit_out, rot_info, N, ci, Gamma, method,
                        optimal_weight = TRUE, scaled = FALSE)
}


# Map a rotation name to the warm-start criterion family and tuning argument used by the compiled
# `.rotation_se_jacobian()`, mirroring the engine selection in `.orth_engines`/`.oblq_engines`
# (R/rotate_model.R). varimax shares the Crawford-Ferguson criterion at kappa = 1 / p, quartimax at
# kappa = 0, equamax at kappa = k / (2 p); quartimin is oblimin at gam = 0. The geomin offset
# `delta` and oblimin `gam` are taken from the resolved criterion arguments (the same `.gpf_crit`
# defaults the rotation itself used). promax and simplimax have no usable analytic Jacobian and are
# rejected before reaching here, so an unrecognised rotation returns NULL.
.rotation_se_method <- function(rotation, p, k, crit_args = list()) {
  gam <- if (is.null(crit_args$gam)) 0 else crit_args$gam
  delta <- if (is.null(crit_args$delta)) 0.01 else crit_args$delta
  switch(rotation,
    varimax   = list(method = "cf",       param = 1 / p,       oblique = FALSE),
    quartimax = list(method = "cf",       param = 0,           oblique = FALSE),
    equamax   = list(method = "cf",       param = k / (2 * p), oblique = FALSE),
    bentlerT  = list(method = "bentler",  param = 0,           oblique = FALSE),
    geominT   = list(method = "geomin",   param = delta,       oblique = FALSE),
    bifactorT = list(method = "bifactor", param = 0,           oblique = FALSE),
    oblimin   = list(method = "oblimin",  param = gam,         oblique = TRUE),
    quartimin = list(method = "oblimin",  param = 0,           oblique = TRUE),
    bentlerQ  = list(method = "bentler",  param = 0,           oblique = TRUE),
    geominQ   = list(method = "geomin",   param = delta,       oblique = TRUE),
    bifactorQ = list(method = "bifactor", param = 0,           oblique = TRUE),
    NULL)
}

# Wald lower/upper interval around an estimate from its standard error. `est` is unclassed so the
# returned bounds are plain numeric matrices/vectors (the rotated estimates carry the "LOADINGS"
# class).
.wald_ci <- function(est, se, z) {
  est <- unclass(est)
  list(lower = est - z * se, upper = est + z * se)
}


# Analytic standard errors and Wald CIs for an obliquely or orthogonally ROTATED solution. The
# caller's covariance of the unrotated loadings (`se0`) is propagated through the rotation by the
# delta method; every analytic path -- expected information, robust sandwich, and the corrected
# two-stage FIML sandwich -- supplies its own, so the propagation is written once. The rotation
# maps the unrotated loadings A
# to the rotated pattern L = A T(A)^{-1 T} (oblique) / A T(A) (orthogonal) and, for oblique
# rotations, the factor correlations Phi = T(A)' T(A); both depend on A alone (not the
# uniquenesses), so only the p*k loading block V of the parameter covariance is needed. The
# rotation Jacobians d vec(L) / d vec(A) and d vec(Phi) / d vec(A) are obtained by
# finite-differencing a warm-started re-rotation of A (`.rotation_se_jacobian()`), which re-solves
# the rotation from the converged transformation so the optimum tracks each perturbation smoothly.
# This is the delta-method standard error of Jennrich (1973), matching lavaan's rotation.se =
# "delta".
# Because rotated quantities are identification-invariant they are comparable across packages
# (unlike the unrotated loadings, whose orientation is identification-dependent).
#
#   Var(vec L)   = J_L V J_L'                                   -> rotated-loading SEs
#   Var(vec Phi) = J_Phi V J_Phi'                              -> factor-correlation SEs (oblique)
#   S = L Phi:  J_S = (Phi' (x) I_p) J_L + (I_k (x) L) J_Phi;  Var(vec S) = J_S V J_S'  (oblique)
#   h2_i = rowSums(A^2)_i = 1 - psi_i (rotation-invariant); SE(h2_i) by the delta method on V
#
# Fills the same SE/CI schema as the bootstrap (rot_loadings/Phi/Structure plus the unrotated
# loadings and uniquenesses); there are no replicate arrays. Falls back to NA rotated SEs (with a
# classed warning) at a Heywood case, a singular information matrix, or a rotation the warm start
# cannot reproduce (e.g. a non-converged transformation). References: Jennrich (1973); Zhang &
# Preacher (2015).
.se_information_rotated <- function(fit_out, rot_info, N, ci, se0) {

  A <- unclass(fit_out$unrot_loadings)
  p <- nrow(A)
  k <- ncol(A)
  pk <- p * k
  psi <- 1 - rowSums(A^2)
  z <- stats::qnorm(1 - (1 - ci) / 2)

  # Unrotated pieces (always reported, mirroring the bootstrap's unrot_loadings) and the loading
  # covariance block propagated through the rotation. Each analytic path supplies its own covariance
  # in `se0` (uniform $vcov/$loadings_se/$uniquenesses_se schema), so the propagation below is
  # identical for all of them and is purely a function of the loading covariance V.
  V <- se0$vcov[seq_len(pk), seq_len(pk), drop = FALSE]
  SE_unrot <- se0$loadings_se
  SE_psi <- se0$uniquenesses_se
  dimnames(SE_unrot) <- dimnames(A)
  names(SE_psi) <- rownames(A)
  # The caller NA's the unrotated SEs whenever the parameter covariance is non-finite OR carries a
  # negative variance on its diagonal (a near-degenerate orientation `solve()` still inverts). Reuse
  # that as the single reliability signal for the rotated quantities too: they are propagated from
  # the same covariance, so if it was too ill-conditioned for the unrotated SEs it must not yield
  # finite rotated SEs either. NA the covariance here, before anything reads it: on the non-PSD path
  # the caller ships a FINITE covariance next to NA marginal SEs, and every quantity derived from it
  # below (the communality delta method, the persisted slot) must fail closed with them rather than
  # report a standard error the caller already judged untrustworthy.
  info_reliable <- !anyNA(SE_unrot)
  if (!info_reliable) V[] <- NA_real_

  # A degenerate rotational gauge is a strictly weaker signal: it diverges only along the gauge
  # directions, which the rotation Jacobian annihilates, so the rotated loadings, Phi, the structure
  # coefficients and the communalities all stay valid and V must NOT be NA'd here. Withhold the
  # unrotated loadings alone -- and only after `info_reliable` has been read, so a degenerate gauge
  # cannot masquerade as an unusable covariance and take the rotated quantities down with it.
  # Absent (rather than FALSE) means the path supplied no gauge diagnostic, e.g. the two-stage FIML
  # sandwich in R/se-fiml.R, which builds its own `se0`; those paths keep their previous behaviour.
  gauge_degenerate <- isFALSE(se0$gauge_reliable)
  if (gauge_degenerate) SE_unrot[] <- NA_real_

  spec <- .rotation_se_method(rot_info$rotation, p, k, rot_info$crit_args)
  rotmat <- rot_info$rotmat
  L_rot <- unclass(rot_info$rot_loadings)
  Phi_pt <- rot_info$Phi
  # Orthogonal rotations and the single-factor fall-back carry no factor correlations, so they
  # report rotated loadings and communalities only.
  oblique <- !is.null(Phi_pt)

  # The communalities are the exact complement of the uniquenesses (h2_i = rowSums(A^2)_i =
  # 1 - psi_i = (L Phi L')_ii), so h2_i and psi_i are one estimand up to sign and must share a
  # standard error. Both are the ordinary delta method on the loading covariance V through the
  # gradient 2 A[i, ] of h2_i (`.communality_se()`); because every path's `se0$uniquenesses_se`
  # comes from that same call on that same V, the two agree exactly.
  h2 <- rowSums(A^2)
  names(h2) <- rownames(A)
  SE_h2 <- .communality_se(A, V)
  names(SE_h2) <- rownames(A)

  na_mat <- function() {
    m <- matrix(NA_real_, p, k)
    dimnames(m) <- dimnames(L_rot)
    m
  }
  SE_rot <- na_mat()
  SE_Phi <- if (oblique) {
    matrix(NA_real_, k, k, dimnames = dimnames(Phi_pt))
  } else {
    NULL
  }
  SE_S <- if (oblique) na_mat() else NULL
  S_pt <- if (oblique) `dimnames<-`(L_rot %*% Phi_pt, dimnames(L_rot)) else NULL

  can_rotate <- !is.null(spec) && is.matrix(rotmat) && !anyNA(rotmat) &&
    k >= 2L && !is.null(L_rot) && info_reliable

  if (can_rotate) {
    normalize <- isTRUE(rot_info$normalize)
    # The bifactor criterion exempts a fixed column as the general factor, but `.reflect_and_order`
    # may have moved the general factor out of the first column. It is the one column that loads on
    # every variable -- so the column whose smallest absolute loading is largest, the others being
    # group factors with near-zero out-group loadings. Identify it so the re-rotation optimises the
    # same criterion the point estimate did; the other criteria ignore general_col.
    general_col <- if (identical(spec$method, "bifactor")) {
      which.max(apply(abs(L_rot), 2, min)) - 1L
    } else {
      0L
    }
    # Forward difference of the warm-started re-rotation (the stencil lavaan's delta method uses);
    # the whole p*k finite-difference loop runs in compiled code (.rotation_se_jacobian()) to avoid
    # p*k round trips to R. This is a one-shot cost per fit: a few hundredths of a second for a
    # typical problem (p = 18, k = 3) and about 1 s for a large one (p = 40, k = 8) -- far below a
    # bootstrap, and an order of magnitude faster than the equivalent pure-R numeric differentiation.
    eps <- 1e-4
    jac <- .rotation_se_jacobian(A, rotmat, spec$method, spec$param, normalize,
                                 spec$oblique, eps, general_col)
    # The re-rotation at the unperturbed A must reproduce the reported rotated loadings; a gross
    # mismatch flags a criterion/parameter mismatch or a transformation outside the criterion's
    # basin, for which the Jacobian is not trustworthy. The tolerance absorbs the small drift
    # between the reported optimum and a re-solve to a tighter tolerance -- and, for varimax, between
    # the SVD/SPSS varimax algorithm and the Crawford-Ferguson re-rotation at a flat optimum -- while
    # still catching a genuine basin/criterion mismatch (which moves loadings an order of magnitude
    # more).
    if (isTRUE(jac$valid) && max(abs(jac$base_loadings - L_rot)) < 5e-2) {
      J_L <- jac$J_L
      SE_rot <- matrix(sqrt(pmax(rowSums((J_L %*% V) * J_L), 0)), p, k)
      dimnames(SE_rot) <- dimnames(L_rot)
      if (spec$oblique) {
        J_Phi <- jac$J_Phi
        SE_Phi <- matrix(sqrt(pmax(rowSums((J_Phi %*% V) * J_Phi), 0)), k, k)
        diag(SE_Phi) <- 0       # the unit diagonal of Phi is fixed, so it has no variance
        SE_Phi <- (SE_Phi + t(SE_Phi)) / 2
        dimnames(SE_Phi) <- dimnames(Phi_pt)
        # Structure S = L Phi: dvec(S) = (Phi' (x) I_p) dvec(L) + (I_k (x) L) dvec(Phi); the
        # rotation matrices are the re-rotated point estimates, consistent with the Jacobians.
        L0 <- jac$base_loadings
        Phi0 <- jac$base_Phi
        J_S <- kronecker(t(Phi0), diag(p)) %*% J_L + kronecker(diag(k), L0) %*% J_Phi
        SE_S <- matrix(sqrt(pmax(rowSums((J_S %*% V) * J_S), 0)), p, k)
        dimnames(SE_S) <- dimnames(L_rot)
      }
    }
  } else if (info_reliable && k < 2L) {
    # No rotation was actually performed (single factor): the rotated loadings equal the unrotated
    # ones, so their SEs are the unrotated loading SEs; there is no Phi or structure matrix.
    SE_rot <- SE_unrot
    dimnames(SE_rot) <- dimnames(L_rot)
  }

  SE <- list(unrot_loadings = SE_unrot, uniquenesses = SE_psi,
             rot_loadings = SE_rot, communalities = SE_h2)
  CI <- list(unrot_loadings = .wald_ci(A, SE_unrot, z),
             uniquenesses = .wald_ci(psi, SE_psi, z),
             rot_loadings = .wald_ci(L_rot, SE_rot, z),
             communalities = .wald_ci(h2, SE_h2, z))
  if (oblique) {
    SE$Phi <- SE_Phi
    SE$Structure <- SE_S
    CI$Phi <- .wald_ci(Phi_pt, SE_Phi, z)
    CI$Structure <- .wald_ci(S_pt, SE_S, z)
  }

  if (anyNA(SE_unrot) || anyNA(SE_psi) || anyNA(SE_rot) || anyNA(SE_h2) ||
      (oblique && (anyNA(SE_Phi) || anyNA(SE_S)))) {
    cli::cli_warn(
      c("Analytic standard errors could not be computed for all parameters.",
        "i" = if (gauge_degenerate && info_reliable) {
          "The factor solution's rotational orientation is only weakly determined (two canonical variances nearly coincide), so the unrotated loadings have no well-defined standard error. The rotated loadings and communalities are gauge-invariant and unaffected."
        } else if (!info_reliable && .at_uniqueness_boundary(A)) {
          # Named rather than listed among the possible causes: `.se_sandwich_dispatch()` withholds
          # on exactly this test, so when it holds it can be reported as a fact about the solution.
          # `!info_reliable` is what makes that attribution safe: the two-stage FIML sandwich builds
          # its own `se0` and carries no boundary gate, so a boundary FIML fit can arrive here with
          # finite unrotated SEs and only the rotation Jacobian failed -- the hint's claim that
          # nothing was reported would be false, and the last clause is the one that applies.
          .se_boundary_hint
        } else {
          "This occurs when the parameter covariance is singular, or when the rotation could not be reproduced for the standard-error Jacobian."
        }),
      class = "efa_se_unreliable"
    )
  }

  # Surface the pk x pk unrotated loading vcov (whichever analytic covariance `se0` carried) so it
  # can be persisted on the EFA object alongside the rotated SEs. It was already NA-filled above
  # when the unrotated SEs were unreliable, so the slot is consistent with the marginal SEs and
  # downstream consumers can fail closed on `anyNA()`. A degenerate gauge deliberately does NOT
  # NA-fill it: the covariance is finite, PSD and correct, and only its gauge-dependent marginals
  # are unusable (see the matching note in `.se_sandwich_unrotated()`).
  lab <- .vec_loading_labels(A)
  if (!is.null(lab)) dimnames(V) <- list(lab, lab)

  list(SE = SE, CI = CI, replicates = NULL, vcov_unrot_loadings = V)
}


# Robust (sandwich) standard errors and a scaled chi-square for the ordinal/polychoric path.
# The factor model is fitted by minimising a weighted off-diagonal discrepancy
# F(theta) = (s - sigma(theta))' V (s - sigma(theta)) over the free loadings theta = vec(Lambda),
# where s = vech of the off-diagonal sample (polychoric/tetrachoric) correlations,
# sigma_ij(theta) = (Lambda Lambda')_ij, and V is the estimator's weight matrix (DWLS: diagonal
# inverse asymptotic variances; ULS: identity; ML: the normal-theory GLS weight). When V is not
# the inverse of the true asymptotic covariance of s (which it never is for ordinal data), the
# naive inverse-information SEs are biased, and the Godambe sandwich
#   Var(theta) = (1/(N-1)) A^- (Delta' V Gamma V Delta) A^-,   A = Delta' V Delta,
# with Gamma the (threshold-adjusted) asymptotic covariance of s, gives consistent SEs (Browne,
# 1984; Muthen, 1984; Satorra & Bentler, 1994). This mirrors lavaan's se = "robust.sem". The
# rotational over-parameterisation of A is handled by bordering it with the identification-
# constraint Jacobian, exactly as the expected-information path does. The bordered loading
# covariance is propagated through a rotation by the same Jacobian machinery as the information
# SEs (.se_information_rotated), so promax/simplimax stay bootstrap-only.
.se_sandwich_dispatch <- function(fit_out, rot_info, N, ci, Gamma, method,
                                  optimal_weight = FALSE, scaled = TRUE) {

  core <- .se_sandwich_core(fit_out, N, Gamma, method, optimal_weight, scaled)

  # A boundary solution invalidates the Wald interval however the covariance was estimated, so the
  # robust path withholds on the same test the expected-information path uses. The gate sits here
  # rather than inside the core for two reasons. The core's own `psi <= 0` test asks a different
  # question -- whether Psi^-1 exists for the gauge constraint -- and must not be moved to the floor
  # (see the note there). And the scaled chi-square the core builds is a discrepancy-function
  # quantity rather than a Wald one: the boundary does not invalidate it, and for DWLS it is the
  # only chi-square block the fit has, so the standard errors are withheld and it is kept. The
  # persisted covariance follows from these flags -- both consumers NA-fill it when the marginal SEs
  # are unusable -- and `gauge_reliable` is reset because the gauge is not what withheld them, so
  # the caller must not attach its gauge-specific message.
  if (.at_uniqueness_boundary(fit_out$unrot_loadings)) {
    core$loadings_se[] <- NA_real_
    core$uniquenesses_se[] <- NA_real_
    core$reliable <- FALSE
    core$gauge_reliable <- TRUE
  }

  res <- if (is.null(rot_info)) {
    .se_sandwich_unrotated(fit_out, core, ci)
  } else {
    # Reuse the information-SE rotation propagation verbatim by supplying the robust loading
    # covariance in place of the expected-information one (identical $vcov schema). Because the
    # rotation Jacobian annihilates the rotational gauge directions, the rotated SEs are
    # identification-invariant and do not depend on which bordering constraint produced V_AA.
    se0 <- list(vcov = core$V_AA,
                loadings_se = core$loadings_se,
                uniquenesses_se = core$uniquenesses_se,
                gauge_reliable = core$gauge_reliable)
    .se_information_rotated(fit_out, rot_info, N, ci, se0 = se0)
  }

  res$scaled_test <- core$scaled_test
  res
}

# Identification-constraint Jacobian for the rotational gauge: the gradient, with respect to the
# column-major vec(Lambda), of the off-diagonal entries the estimator's canonical orientation sets
# to zero. nc = k(k-1)/2 rows, p*k columns; for k = 1 there is no gauge freedom and the matrix has
# zero rows. The constraint must match the orientation the loadings are reported in so the unrotated
# loading SEs are scaled in that gauge (the rotated SEs are gauge-invariant and the same for any
# transversal choice). With `psi = NULL` (eigen-based ULS/DWLS) it fixes off-diag(Lambda'Lambda);
# with `psi` supplied (ML, where psi = 1 - rowSums(Lambda^2)) it fixes off-diag(Lambda' Psi^-1
# Lambda) (Lawley & Maxwell), whose gradient carries a chain-rule term through psi(Lambda).
.se_sandwich_constraint <- function(L, psi = NULL) {
  p <- nrow(L)
  k <- ncol(L)
  nc <- k * (k - 1L) / 2L
  Cmat <- matrix(0, nc, p * k)
  if (k > 1L) {
    Astar <- if (is.null(psi)) NULL else L / psi   # Psi^-1 Lambda for the ML gauge
    uv <- 0L
    for (v in 2:k) {
      for (u in seq_len(v - 1L)) {
        uv <- uv + 1L
        grad <- if (is.null(psi)) {
          .gauge_grad(L, u, v)
        } else {
          # d off-diag(Lambda' Psi^-1 Lambda)_uv / d vec(Lambda): the direct part plus the term
          # from psi_a = 1 - rowSums(Lambda^2)_a depending on Lambda (d psi_a / d lambda_ac = -2
          # lambda_ac), which contributes 2 (Astar[, u] Astar[, v]) Lambda across all columns.
          .gauge_grad(Astar, u, v) + 2 * (Astar[, u] * Astar[, v]) * L
        }
        Cmat[uv, ] <- as.vector(grad)
      }
    }
  }
  Cmat
}

# Smallest gauge-transversal conditioning a solution can have before its rotational orientation is
# treated as degenerate. The quantity `.gauge_transversal()` returns is a cosine in [0, 1] -- the
# smallest principal angle between the gauge-fixing constraint's row space and the gauge directions
# it has to pin down -- so this is a dimensionless angle, not a parameter scale. Blow-ups in the
# reported standard errors set in below about 0.05 and well-determined solutions sit above 0.2, so
# the floor flags only a clearly degenerate orientation.
.gauge_transversal_floor <- 0.05

# Conditioning of the gauge-fixing constraint as a transversal to the rotational gauge orbit.
#
# The bordered inverse of `A` under the constraint `Cmat` is the constrained (reflexive generalised)
# inverse (Rao & Mitra, 1971, ch. 3; Silvey, 1959)
#   A^- = (I - Z (C Z)^-1 C) A^+ (I - C' (C Z)^-1' Z'),
# where Z spans null(A). The whole conditioning of the augmented system therefore sits in the
# k(k-1)/2 x k(k-1)/2 matrix C Z, and the reported unrotated loading covariance is amplified by
# ||(C Z)^-1||. The null space is available in closed form: the gauge directions are dLambda =
# Lambda S for antisymmetric S, since d(Lambda Lambda') = Lambda (S + S') Lambda' = 0, so no
# eigendecomposition of A is needed.
#
# C Z degenerates when the canonical orientation is not pinned down. Dropping the chain-rule term,
# C Z is diagonal with entries (d_u - d_v) for d = diag(Lambda' Psi^-1 Lambda): when two of the
# canonical variances that identify the ML solution collide (Lawley & Maxwell, 1971, sec. 2.3), the
# orientation within that two-plane is arbitrary and the unrotated loadings have no well-defined
# limiting distribution -- the classic rotational indeterminacy of the unrotated factor solution.
# The divergence is genuine, not a numerical artefact -- but it is confined to the gauge, so it must
# be reported rather than reported around: the rotated loadings, Phi and the structure coefficients
# are gauge-invariant (the rotation Jacobian annihilates these directions) and the communalities
# rowSums(Lambda^2) are likewise invariant, so all of them stay valid and are left untouched by the
# caller.
#
# Returns the smallest singular value of the row-normalised constraint restricted to the gauge
# directions -- dimensionless and in [0, 1] -- or Inf when there is no gauge freedom to fix (k = 1).
.gauge_transversal <- function(L, Cmat) {

  k <- ncol(L)
  if (nrow(Cmat) == 0L || k < 2L) return(Inf)
  if (!all(is.finite(L)) || !all(is.finite(Cmat))) return(0)

  # Z = [vec(L S_uv)] over the antisymmetric basis, in the same (u, v) order as the constraint rows.
  # L %*% S_uv is just two sign-flipped column copies (column u becomes -L[, v], column v becomes
  # L[, u]), so fill it directly rather than forming S and multiplying.
  Z <- matrix(0, length(L), nrow(Cmat))
  col <- 0L
  for (v in 2:k) {
    for (u in seq_len(v - 1L)) {
      col <- col + 1L
      Zc <- matrix(0, nrow(L), k)
      Zc[, u] <- -L[, v]
      Zc[, v] <- L[, u]
      Z[, col] <- as.vector(Zc)
    }
  }
  # A rank-deficient Z means the loadings themselves are collinear across factors, so the gauge is
  # not even a well-defined orbit; treat it as maximally degenerate.
  qrZ <- qr(Z)
  if (qrZ$rank < ncol(Z)) return(0)
  Z <- qr.Q(qrZ)

  # Normalise the constraint rows so the singular values measure orientation only: the constraint is
  # a set of equations fixed at zero, so its overall row scaling is arbitrary and must not enter.
  rn <- sqrt(rowSums(Cmat^2))
  if (any(rn <= 0)) return(0)

  min(svd((Cmat / rn) %*% Z)$d)
}

# Normal-theory asymptotic covariance of the off-diagonal Pearson correlations, on the unit scale
# (N * Cov(r-hat)) and in the `pairs` (utils::combn(p, 2)) order the sandwich machinery uses
# throughout. `Sigma` is the correlation matrix it is evaluated at -- the model-implied one for the
# expected-information SEs.
#
# For pairs (i, j) and (k, l) (Pearson & Filon, 1898; Olkin & Siotani, 1976; Steiger, 1980):
#   N Cov(r_ij, r_kl) = 1/2 rho_ij rho_kl (rho_ik^2 + rho_il^2 + rho_jk^2 + rho_jl^2)
#                       + rho_ik rho_jl + rho_il rho_jk
#                       - rho_ij (rho_ik rho_il + rho_jk rho_jl)
#                       - rho_kl (rho_ik rho_jk + rho_il rho_jl),
# which reduces to the familiar (1 - rho^2)^2 on the diagonal. Assembled by whole-matrix indexing
# (`Sigma[ii, ii]` and friends are already n x n in pair order) rather than looping over the n^2
# pair combinations; the result is symmetric analytically and is symmetrised against rounding.
.normal_theory_gamma <- function(Sigma, pairs) {

  ii <- pairs[1, ]
  jj <- pairs[2, ]
  n <- ncol(pairs)

  # Drop the variable dimnames while sub-setting: `Sigma[ii, ii]` would otherwise label the pair
  # axis with repeated variable names ("V1", "V1", ...), which mislabels a pair-indexed matrix and
  # would propagate into the covariance built from it.
  Sigma <- unname(Sigma)
  r_ij <- Sigma[cbind(ii, jj)]
  Rik <- Sigma[ii, ii, drop = FALSE]
  Ril <- Sigma[ii, jj, drop = FALSE]
  Rjk <- Sigma[jj, ii, drop = FALSE]
  Rjl <- Sigma[jj, jj, drop = FALSE]

  # rho_ij indexes the row pair, rho_kl the column pair: `r_ij * M` scales rows (column-major
  # recycling) and `rep(r_ij, each = n) * M` scales columns.
  G <- 0.5 * outer(r_ij, r_ij) * (Rik^2 + Ril^2 + Rjk^2 + Rjl^2) +
    Rik * Rjl + Ril * Rjk -
    r_ij * (Rik * Ril + Rjk * Rjl) -
    rep(r_ij, each = n) * (Rik * Rjk + Ril * Rjl)

  (G + t(G)) / 2
}

# Build the robust loading covariance V_AA (p*k x p*k), the unrotated loading/uniqueness SEs, and
# the scaled chi-square from the fitted loadings and the polychoric ACOV `Gamma`. `Gamma` enters on
# the variance scale Var(rho-hat) and is converted here to the unit asymptotic-variance scale
# (N * Var) used by the WLS/sandwich formulas (= lavaan's NACOV). Returns NA SEs (reliable = FALSE)
# at a singular bordered information matrix or an unusable covariance.
#
# The linear algebra is O(n^2 q) for the meat and O(n^3) for the chi-square trace (n = p(p-1)/2
# pairs, q = p*k), and is kept in R: it is a one-shot cost of a few hundredths of a second for a
# typical problem and about 0.3 s (p = 30) to 1.5 s (p = 40, n = 780 pairs) for a large one --
# negligible next to the polychoric estimate it builds on, and far below a bootstrap.
#
# `optimal_weight = TRUE` replaces the estimator weight with V = Gamma^-1, at which the sandwich
# collapses onto its bread and returns the Fisher information of the correlation structure (the
# se = "information" path; see `.se_information()`). `scaled = FALSE` skips the scaled chi-square,
# which only the robust paths report.
.se_sandwich_core <- function(fit_out, N, Gamma, method,
                              optimal_weight = FALSE, scaled = TRUE) {

  L <- unclass(fit_out$unrot_loadings)
  p <- nrow(L)
  k <- ncol(L)
  pk <- p * k

  # `gauge_reliable = TRUE` on the NA path: everything is already withheld, and the gauge is not
  # what withheld it, so the caller's gauge-specific message must not be attached to this branch.
  na_core <- list(
    V_AA = matrix(NA_real_, pk, pk),
    loadings_se = matrix(NA_real_, p, k),
    uniquenesses_se = rep(NA_real_, p),
    scaled_test = NULL,
    reliable = FALSE,
    gauge_reliable = TRUE
  )

  if (is.null(Gamma) || anyNA(Gamma) || nrow(Gamma) != p * (p - 1L) / 2L) {
    return(na_core)
  }

  # Variance scale (Var(rho-hat)) -> unit asymptotic-variance scale (N * Var); the 1/(N-1) on the
  # final covariance and the N on the chi-square then follow lavaan's robust.sem conventions.
  Gamma <- N * Gamma

  pairs <- utils::combn(p, 2L)
  pi <- pairs[1, ]
  pj <- pairs[2, ]
  n <- ncol(pairs)

  # Model Jacobian Delta = d sigma_offdiag / d vec(Lambda) (n x pk). For pair (i, j),
  # d sigma_ij / d Lambda[a, f] = (a == i) Lambda[j, f] + (a == j) Lambda[i, f]; rows in
  # utils::combn(p, 2) order, columns in column-major vec(Lambda) order.
  Delta <- matrix(0, n, pk)
  for (f in seq_len(k)) {
    Delta[cbind(seq_len(n), (f - 1L) * p + pi)] <- L[pj, f]
    Delta[cbind(seq_len(n), (f - 1L) * p + pj)] <- L[pi, f]
  }

  # Estimator weight V (unit scale). DWLS: diagonal inverse variances; ULS: identity; ML: the
  # normal-theory GLS weight 1/2 (Sigma^-1 (x) Sigma^-1) restricted to the off-diagonal pairs, at
  # the model-implied correlation matrix Sigma = Lambda Lambda' (unit diagonal).
  if (optimal_weight) {
    # Efficient weight V = Gamma^-1. The meat Delta' V Gamma V Delta then equals the bread
    # Delta' V Delta = Delta' Gamma^-1 Delta, so V_AA reduces to the (bordered) inverse of the
    # correlation-structure Fisher information.
    Vmat <- tryCatch(solve(Gamma), error = function(e) NULL)
    if (is.null(Vmat)) return(na_core)
    Vmat <- (Vmat + t(Vmat)) / 2
    vdiag <- NULL
  } else if (method == "DWLS") {
    vdiag <- 1 / diag(Gamma)
    Vmat <- NULL
  } else if (method == "ULS") {
    vdiag <- rep(1, n)
    Vmat <- NULL
  } else {
    Sigma <- tcrossprod(L)
    diag(Sigma) <- 1
    P <- tryCatch(solve(Sigma), error = function(e) NULL)
    if (is.null(P)) return(na_core)
    vdiag <- NULL
    Vmat <- matrix(0, n, n)
    for (s in seq_len(n)) {
      a <- pi[s]
      b <- pj[s]
      Vmat[, s] <- 0.5 * (P[pi, a] * P[pj, b] + P[pi, b] * P[pj, a])
    }
  }
  is_diag <- is.null(Vmat)

  # VD = V Delta; bread A = Delta' V Delta (singular by the k(k-1)/2 rotational gauge freedoms).
  VD <- if (is_diag) vdiag * Delta else Vmat %*% Delta
  A <- crossprod(Delta, VD)

  # Border A with the gauge-fixing constraint so the augmented system is invertible; the leading
  # pk-block of its inverse is the reflexive generalised inverse of A. The constraint fixes the
  # rotational orientation, so it must match the one the solution's loadings are reported in: the
  # unrotated loading SEs are scaled in that gauge (the rotated SEs and scaled chi-square are
  # gauge-invariant, so the choice does not affect them). Rather than assume the orientation from the
  # estimator label, detect it from the loadings: an eigen-based ULS/DWLS solution is
  # Lambda'Lambda-diagonal, an ML solution is Lambda'Psi^-1 Lambda-diagonal (Lawley & Maxwell). The
  # orientation the solution is in leaves a (scale-free) relative off-diagonal at the optimiser floor
  # while the other is O(0.01-1), so the smaller ratio identifies the gauge -- and detecting keeps it
  # tied to the actual solution, so a future estimator with either identification is handled without
  # special-casing. Two cases the loadings cannot resolve fall back to a fixed choice: the
  # Lambda'Psi^-1 Lambda gauge needs Psi^-1, so it is unavailable where that is undefined (psi <= 0)
  # or at a non-finite solution -- routed to Lambda'Lambda; and homogeneous uniquenesses (Psi
  # proportional to I) make BOTH orientations diagonal, so the gauge is undetermined by the loadings
  # (yet the two still give different SEs through the chain-rule term), broken by the estimator's
  # identification. A single factor has no rotational freedom (the constraint is empty either way).
  #
  # The test here is psi <= 0 and deliberately NOT the `.uniqueness_floor` boundary that
  # `.se_information()` treats as a Heywood case. The two ask different questions: that one asks
  # whether the solution sits on the parameter-space boundary (where a Wald interval is invalid),
  # this one only whether Psi^-1 exists. At a uniqueness pinned at the floor Psi^-1 is large but
  # perfectly well defined, and the solution is still oriented Lambda'Psi^-1 Lambda-diagonal -- so
  # forcing the Lambda'Lambda constraint there would mismatch the orientation the loadings are
  # reported in and hand `.gauge_transversal()` a near-tangent transversal (0.002 against 0.16 for
  # the matching constraint on a boundary ML fit), withholding standard errors that are fine.
  psi <- 1 - rowSums(L^2)
  use_ltpil_gauge <- if (k < 2L || anyNA(psi) || any(psi <= 0)) {
    FALSE
  } else {
    rel_off <- function(M) max(abs(M[upper.tri(M)])) / max(abs(diag(M)), .Machine$double.eps)
    r_ltl <- rel_off(crossprod(L))
    r_ltpil <- rel_off(crossprod(L, L / psi))
    # Both orientations near-diagonal (Psi proportional to I): the loadings cannot distinguish the
    # gauge, so use the estimator's known identification; otherwise the smaller ratio wins.
    if (max(r_ltl, r_ltpil) < 1e-4) method == "ML" else r_ltpil < r_ltl
  }
  Cmat <- if (use_ltpil_gauge) {
    .se_sandwich_constraint(L, psi = psi)
  } else {
    .se_sandwich_constraint(L)
  }
  Abread <- .bordered_inverse_block(A, Cmat, pk)
  if (is.null(Abread)) return(na_core)

  # Conditioning of that constraint as a transversal to the gauge orbit. A near-degenerate
  # transversal amplifies the unrotated loading covariance by ||(C Z)^-1|| without making it
  # non-finite or non-PSD, so it passes the `.is_psd()` gate below and would otherwise ship a
  # silently meaningless standard error. Flagged separately from `reliable` because the divergence
  # lives purely in the gauge: everything gauge-invariant the caller derives from V_AA stays valid.
  # Computed after the bordered inverse so an exactly singular augmented system -- which returns
  # through `na_core` and discards this flag -- does not pay for the decomposition.
  gauge_reliable <- .gauge_transversal(L, Cmat) >= .gauge_transversal_floor

  # Meat Delta' V Gamma V Delta = (V Delta)' Gamma (V Delta); robust covariance A^- meat A^- /(N-1).
  Gamma_theta <- crossprod(VD, Gamma %*% VD)
  V_AA <- (Abread %*% Gamma_theta %*% Abread) / (N - 1)

  # The covariance must be positive semidefinite: a non-finite entry or a negative eigenvalue (even
  # with a still-positive diagonal) would otherwise corrupt the marginal, uniqueness, and rotated
  # SEs, which read the full V_AA.
  reliable <- .is_psd(V_AA)
  loadings_se <- if (reliable) sqrt(diag(V_AA)) else rep(NA_real_, pk)

  # Uniqueness SE = communality SE (psi_i = 1 - rowSums(Lambda^2)_i), via the shared gradient.
  uniq_se <- if (!reliable) rep(NA_real_, p) else .communality_se(L, V_AA)

  scaled_test <- if (!reliable || !scaled) NULL else {
    .scaled_chisq(fit_out, Gamma, pairs, VD, vdiag, Vmat, Abread, N)
  }

  list(
    V_AA = V_AA,
    loadings_se = matrix(loadings_se, p, k),
    uniquenesses_se = uniq_se,
    scaled_test = scaled_test,
    reliable = reliable,
    gauge_reliable = gauge_reliable
  )
}

# Scaled chi-square test statistics (Satorra & Bentler, 1994; Asparouhov & Muthen, 2010) for the
# weighted off-diagonal fit. T = N * (s - sigma)' V (s - sigma) is the (unscaled) fit statistic;
# U = V - V Delta A^- Delta' V is the residual projector; c1 = tr(U Gamma), c2 = tr((U Gamma)^2)
# its trace coefficients (Gamma on the unit scale). The projector spans only the off-diagonal
# correlation residuals, so this is the two-stage correlation-structure correction (Browne, 1984)
# and is not identical to the full WLSMV statistic of lavaan/Mplus, which also projects the
# thresholds. Returns the scaled-shifted, mean-adjusted, and mean-and-variance-adjusted statistics
# plus the scaled baseline statistic for the robust CFI/TLI/RMSEA. NULL when the model is
# just-identified (df <= 0) or the traces degenerate.
.scaled_chisq <- function(fit_out, Gamma, pairs, VD, vdiag, Vmat, Abread, N) {

  L <- unclass(fit_out$unrot_loadings)
  df <- fit_out$fit_indices$df
  if (is.null(df) || is.na(df) || df <= 0) return(NULL)

  n <- ncol(pairs)
  is_diag <- is.null(Vmat)
  apply_V <- function(x) if (is_diag) vdiag * x else as.vector(Vmat %*% x)

  # The sample correlations s and residuals e are read from the (possibly smoothed) analysis matrix,
  # while Gamma is the asymptotic covariance of the observed un-projected polychoric correlations --
  # the same convention the DWLS point estimate and lavaan use (the weights/ACOV describe the
  # observed correlations even when the matrix is projected to positive definiteness).
  R <- unclass(fit_out$orig_R)
  s <- R[t(pairs)]                       # off-diagonal sample correlations
  e <- (R - tcrossprod(L))[t(pairs)]     # residuals r_ij - (Lambda Lambda')_ij

  Tstat <- N * sum(e * apply_V(e))

  # Residual projector U = V - V Delta A^- Delta' V. The dense V Delta A^- Delta' V term makes U
  # dense anyway, so add the (diagonal or full) weight onto it in place instead of materialising a
  # dense V.
  U <- -(VD %*% Abread %*% t(VD))
  if (is_diag) diag(U) <- diag(U) + vdiag else U <- U + Vmat
  c1 <- sum(U * Gamma)
  UG <- U %*% Gamma
  c2 <- sum(UG * t(UG))
  if (!is.finite(c1) || !is.finite(c2) || c1 <= 0 || c2 <= 0) return(NULL)

  # Independence baseline (all off-diagonal correlations fixed to 0): residual is s, projector the
  # baseline weight V0. V0 is diagonal for every estimator -- it equals the (model-independent)
  # weight for DWLS/ULS, and 1/2 I for ML (the normal-theory weight at Sigma = I) -- so the baseline
  # traces use only its diagonal.
  v0 <- if (is_diag) vdiag else rep(0.5, n)
  Tbase <- N * sum(s * (v0 * s))
  df_base <- n
  c1_base <- sum(v0 * diag(Gamma))       # tr(diag(v0) Gamma)
  V0G <- v0 * Gamma                      # diag(v0) %*% Gamma (row scaling)
  c2_base <- sum(V0G * t(V0G))           # tr((diag(v0) Gamma)^2)

  scaled <- .scaled_variants(Tstat, df, c1, c2)
  chi_null <- if (is.finite(c1_base) && is.finite(c2_base) && c2_base > 0) {
    .scaled_variants(Tbase, df_base, c1_base, c2_base)$T_ss
  } else {
    NA_real_
  }

  list(
    m = nrow(L),
    df = df,
    chi = scaled$T_ss,                   # scaled-shifted (WLSMV default) -> the reported chi-square
    chi_scaling = scaled$a,
    chi_shift = scaled$b,
    chi_unscaled = Tstat,
    chi_mean_adjusted = scaled$T_mean,
    chi_mean_var = scaled$T_mv,
    df_mean_var = scaled$df_mv,
    chi_null = chi_null,
    df_null = df_base
  )
}

# The three scalings of a fit statistic T (Satorra & Bentler, 1994; Asparouhov & Muthen, 2010)
# from its trace coefficients c1 = tr(U Gamma), c2 = tr((U Gamma)^2) and degrees of freedom df:
# mean-adjusted (T * df/c1, reference df); scaled-shifted a*T + b with a = sqrt(df/c2),
# b = df - a*c1 (reference df, the WLSMV default); mean-and-variance-adjusted T * c1/c2 with the
# Satterthwaite df* = c1^2/c2.
.scaled_variants <- function(Tstat, df, c1, c2) {
  a <- sqrt(df / c2)
  list(
    T_mean = Tstat * df / c1,
    a = a,
    b = df - a * c1,
    T_ss = a * Tstat + (df - a * c1),
    df_mv = c1^2 / c2,
    T_mv = Tstat * c1 / c2
  )
}

# Unrotated analytic SE/CI wrapper: fills the uniform list(SE, CI, replicates) schema from the
# core's loading/uniqueness SEs. Shared by every analytic path -- expected information, robust
# sandwich, and the corrected two-stage FIML sandwich -- which differ only in the covariance the
# core built (`.se_information_rotated()` plays the same role under a rotation).
.se_sandwich_unrotated <- function(fit_out, core, ci) {

  L <- unclass(fit_out$unrot_loadings)
  SE_L <- core$loadings_se
  SE_psi <- core$uniquenesses_se
  dimnames(SE_L) <- dimnames(L)
  names(SE_psi) <- rownames(L)

  # A degenerate rotational gauge inflates the unrotated loading SEs without bound while leaving
  # every gauge-invariant quantity intact, so withhold the unrotated loadings alone and keep the
  # communalities and uniquenesses, which are invariant under Lambda -> Lambda T and unaffected.
  gauge_degenerate <- isFALSE(core$gauge_reliable)
  if (gauge_degenerate) SE_L[] <- NA_real_

  if (anyNA(SE_L) || anyNA(SE_psi)) {
    cli::cli_warn(
      c("Analytic standard errors could not be computed for all parameters.",
        "i" = if (gauge_degenerate) {
          "The factor solution's rotational orientation is only weakly determined (two canonical variances nearly coincide), so the unrotated loadings have no well-defined standard error. The communalities and uniquenesses are unaffected; apply a rotation, or use {.code se = \"np-boot\"}, for loading-level uncertainty."
        } else if (anyNA(SE_L) && .at_uniqueness_boundary(L)) {
          # As in `.se_information_rotated()`: the hint claims a total withholding, so it is only
          # attached once the loading SEs are actually gone. They and the uniqueness SEs are
          # all-or-nothing together today, but a path that withheld only one of them would
          # otherwise be described by a hint that its own output contradicts.
          .se_boundary_hint
        } else {
          "This occurs when the bordered information matrix is singular, or when the asymptotic covariance is not usable."
        }),
      class = "efa_se_unreliable"
    )
  }

  # h2_i = 1 - psi_i exactly, so the communality is the uniqueness up to sign and shares its
  # standard error verbatim (the delta method the core already ran through the gradient
  # 2 L[i, ]). Report it here too, so an unrotated fit does not force the user to know the
  # identity to obtain communality intervals.
  h2 <- rowSums(L^2)
  names(h2) <- rownames(L)
  psi <- 1 - h2
  SE_h2 <- SE_psi
  z <- stats::qnorm(1 - (1 - ci) / 2)

  # `.se_sandwich_core()` always returns a numeric V_AA, even when `reliable = FALSE` (it is
  # only the marginal `loadings_se` that gets NA-filled in the unreliable branch). NA-fill the
  # persisted covariance to match, so a finite-but-not-PSD V_AA does not silently ship next to
  # NA SEs and propagate sqrt(NaN) into downstream pooling.
  # A degenerate gauge does NOT NA-fill the persisted block. The documented contract for that slot
  # is "NA-filled if the analytic covariance is unreliable (a Heywood case or a singular bordered
  # information matrix)", and a weakly determined orientation is neither: the covariance is finite,
  # PSD and -- verified against the closed-form reflexive inverse -- correct. What is unusable is
  # reading its diagonal as the standard error of a parameter the gauge does not pin down. Keeping
  # it finite also leaves consumers that pool in a COMMON gauge across fits (`efa_mi()`) able to
  # recover the gauge-invariant quantities, which a wholesale NA fill would take down with it.
  V_AA <- core$V_AA
  if (!isTRUE(core$reliable)) V_AA[] <- NA_real_
  lab <- .vec_loading_labels(L)
  if (!is.null(lab)) dimnames(V_AA) <- list(lab, lab)

  list(
    SE = list(unrot_loadings = SE_L, uniquenesses = SE_psi,
              communalities = SE_h2),
    CI = list(unrot_loadings = .wald_ci(L, SE_L, z),
              uniquenesses = .wald_ci(psi, SE_psi, z),
              communalities = .wald_ci(h2, SE_h2, z)),
    replicates = NULL,
    vcov_unrot_loadings = V_AA
  )
}

# Fold the robust scaled chi-square into the fit indices, overwriting the chi-square-derived block
# (.gof() leaves it undefined for DWLS and reports the non-robust discrepancy for ML/ULS). The
# reported chi-square is the scaled-shifted (WLSMV-default) statistic; the mean-adjusted and
# mean-and-variance-adjusted statistics and the scaling/shift are added as extra fields. CFI/TLI/
# RMSEA come from the scaled model and baseline statistics via the shared .chi_fit_indices().
.apply_scaled_test <- function(fit_indices, st, N) {

  if (is.null(st)) return(fit_indices)

  # The scaled (sandwich) model and baseline statistics are already on a comparable scale, so they
  # serve directly as the CFI/TLI/RMSEA noncentrality inputs (chi_cfi / chi_null_cfi).
  idx <- .chi_fit_indices(st$chi, st$df, st$chi_null, st$df_null, N, st$m, ci = TRUE,
                          chi_cfi = st$chi, chi_null_cfi = st$chi_null)

  fit_indices$chi <- st$chi
  fit_indices$df <- st$df
  fit_indices$p_chi <- idx$p_chi
  fit_indices$CFI <- idx$CFI
  fit_indices$TLI <- idx$TLI
  fit_indices$RMSEA <- idx$RMSEA
  fit_indices$RMSEA_LB <- idx$RMSEA_LB
  fit_indices$RMSEA_UB <- idx$RMSEA_UB
  # AIC/BIC/ECVI are likelihood-ratio chi-square information criteria; they have no standard
  # interpretation when built on the moment-scaled (Satorra-Bentler) statistic, so leave them NA
  # rather than report a misleading model-comparison number.
  fit_indices$AIC <- NA_real_
  fit_indices$BIC <- NA_real_
  fit_indices$ECVI <- NA_real_
  fit_indices$chi_null <- st$chi_null
  fit_indices$df_null <- st$df_null
  fit_indices$p_null <- idx$p_null

  fit_indices$chi_scaled_type <- "scaled.shifted"
  fit_indices$chi_scaling <- st$chi_scaling
  fit_indices$chi_shift <- st$chi_shift
  fit_indices$chi_unscaled <- st$chi_unscaled
  fit_indices$chi_mean_adjusted <- st$chi_mean_adjusted
  fit_indices$chi_mean_var <- st$chi_mean_var
  fit_indices$df_mean_var <- st$df_mean_var

  fit_indices
}


.boot_se_ci <- function(fit_target, L_rot, boot_fit, boot_rot, ci, b) {

  l_ci <- (1 - ci) / 2
  ps <- c(l_ci, ci + l_ci)

  ### calculate stats for unrot loadings and gof measures
  L_unrot <- fit_target$unrot_loadings

  ncol_L <- ncol(L_unrot)
  nrow_L <- nrow(L_unrot)
  colnam_L <- colnames(L_unrot)
  rownam_L <- rownames(L_unrot)

  L_unrot_boot <- array(NA_real_, c(nrow_L, ncol_L, b),
                        dimnames = list(rownam_L, colnam_L,
                                        NULL))
  # The bootstrap aggregates only the numeric fit indices: a scaled chi-square fit (the
  # cor_method = "fiml" two-stage statistic, or any se = "sandwich" fit) adds the character
  # `chi_scaled_type` tag and the extra scaled-statistic components, which the quantile/SD
  # aggregation below cannot consume. Aligning each replicate to these names also lets a
  # replicate whose scaled statistic degenerated (and so lacks a component) contribute NA
  # there rather than shifting every column.
  gof_names <- names(fit_target$fit_indices)[
    vapply(fit_target$fit_indices, is.numeric, logical(1))]
  # The RMSEA confidence bounds are dropped: the replicate fits are run with ci = FALSE, so
  # they are NA in every replicate, and a bootstrap standard error of an interval bound is
  # not a meaningful quantity in the first place.
  gof_names <- setdiff(gof_names, c("RMSEA_LB", "RMSEA_UB"))
  gof_boot <- matrix(NA_real_, ncol = length(gof_names), nrow = b,
                     dimnames = list(NULL, gof_names))
  residuals_boot <- array(NA_real_, c(nrow_L, nrow_L, b),
                          dimnames = list(rownam_L, rownam_L,
                                          NULL))

  # Track replicates that could not be fit (NULL) or aligned, so they are
  # excluded from the bootstrap statistics rather than aborting the whole call.
  failed <- vapply(boot_fit, is.null, logical(1))

  for (boot_i in seq_len(b)) {

    if (failed[boot_i]) next

    # save aligned loading matrix
    aligned <- tryCatch(
      .align_solution(L_unrot, boot_fit[[boot_i]]$unrot_loadings),
      error = function(e) NULL
    )
    if (is.null(aligned)) {
      failed[boot_i] <- TRUE
      next
    }

    L_unrot_boot[,, boot_i] <- aligned$loadings
    fi_i <- boot_fit[[boot_i]]$fit_indices
    gof_boot[boot_i, ] <- vapply(gof_names, function(nm) {
      v <- fi_i[[nm]]
      if (is.null(v) || !is.numeric(v)) NA_real_ else as.numeric(v[1L])
    }, numeric(1))
    residuals_boot[,, boot_i] <- boot_fit[[boot_i]]$residuals

  }

  n_failed <- sum(failed)
  # Reported on the object so the effective B behind the standard errors stays recoverable
  # from a saved fit, not only from the warning below.
  valid_reps <- b - n_failed
  if (n_failed == b) {
    cli::cli_abort(
      c("All {b} bootstrap replicates failed; no bootstrap standard errors could be computed.",
        "i" = "The resampled correlation matrices may be degenerate; try more observations or fewer factors."),
      class = "efa_boot_all_failed"
    )
  }
  if (n_failed > 0) {
    cli::cli_warn(
      c("{n_failed} bootstrap replicate{?s} failed and {?was/were} excluded.",
        "i" = "Bootstrap standard errors and confidence intervals are based on {b - n_failed} replicate{?s}."),
      class = "efa_boot_replicate_failed"
    )
  }
  # A bootstrap standard error is the dispersion across replicates, so fewer than two leaves it
  # undefined (`sd()` of a single value is NA) and collapses the percentile interval onto that one
  # replicate. `b_boot` is bounded below at 2, so this is reachable only when replicates fail at
  # run time -- including the extreme case where all but one of a large `b_boot` do.
  if (valid_reps < 2L) {
    cli::cli_warn(
      c("Bootstrap standard errors are not defined from {valid_reps} usable replicate{?s}.",
        "i" = "The standard errors are {.val {NA}} and the confidence bounds collapse onto the single replicate; raise {.arg b_boot}, or check why the replicate fits failed."),
      class = "efa_se_unreliable"
    )
  }

  # se = sd of bootstrap replications (Zhang, 2014, Estimating Standard Errors
  # in Exploratory Factor Analysis)
  L_unrot_se_ci <- .array_se_ci(L_unrot_boot, ps)
  gof_se_ci <- .array_se_ci(gof_boot, ps, M = 2)
  residuals_se_ci <- .array_se_ci(residuals_boot, ps)
  names(gof_se_ci$se) <- gof_names
  names(gof_se_ci$ci$lower) <- gof_names
  names(gof_se_ci$ci$upper) <- gof_names

  if(boot_rot == "oblique") {

    colnam_L <- colnames(L_rot)

    L_rot_boot <- array(NA_real_, c(nrow_L, ncol_L, b),
                        dimnames = list(rownam_L, colnam_L,
                                        NULL))
    Phi_rot_boot <- array(NA_real_, c(ncol_L, ncol_L, b),
                          dimnames = list(colnam_L, colnam_L,
                                          NULL))
    Structure_boot <- array(NA_real_, c(nrow_L, ncol_L, b),
                        dimnames = list(rownam_L, colnam_L,
                                        NULL))

    failed_rot <- 0

    # Align every successfully fit replicate to the target in a single compiled
    # call over the loading cube, rather than one efa_procrustes() round trip per
    # replicate. Replicates that failed to fit or could not be aligned to the
    # point estimate are excluded from the cube and stay NA in the output arrays.
    keep <- which(!failed)
    if (length(keep) > 0) {
      A_cube <- array(NA_real_, c(nrow_L, ncol_L, length(keep)))
      for (j in seq_along(keep)) {
        A_cube[, , j] <- boot_fit[[keep[j]]]$unrot_loadings
      }

      aligned <- .oblique_procrustes_batch(A_cube, L_rot, random_starts = 5)

      for (j in seq_along(keep)) {
        boot_i <- keep[j]
        # Exclude a replicate only when no valid alignment could be produced. The
        # best multi-start fit is kept even if it did not formally converge: it is
        # the lowest-objective alignment available and its loadings are well-defined.
        if (!isTRUE(aligned$valid[j])) {
          failed_rot <- failed_rot + 1
          next
        }
        L_j <- matrix(aligned$loadings[, , j], nrow_L, ncol_L)
        Phi_j <- matrix(aligned$Phi[, , j], ncol_L, ncol_L)
        L_rot_boot[, , boot_i] <- L_j
        Phi_rot_boot[, , boot_i] <- Phi_j
        Structure_boot[, , boot_i] <- L_j %*% Phi_j
      }
    }

    valid_rot <- b - n_failed - failed_rot
    if (failed_rot > 0) {
      cli::cli_warn(c("{failed_rot} target rotation{?s} in the bootstrap procedure could not be aligned.",
                    "i" = "Bootstrap SE and CI of rotated loadings, factor correlations and structure coefficients are based on {valid_rot} bootstrap sample{?s}."),
                    class = "efa_boot_rotation_failed")
    }

    L_rot_se_ci <- .array_se_ci(L_rot_boot, ps)
    Phi_rot_se_ci <- .array_se_ci(Phi_rot_boot, ps)
    Structure_se_ci <- .array_se_ci(Structure_boot, ps)

    out <- list(
      SE = list(
        unrot_loadings = L_unrot_se_ci$se,
        rot_loadings = L_rot_se_ci$se,
        Phi = Phi_rot_se_ci$se,
        Structure = Structure_se_ci$se,
        fit_indices = gof_se_ci$se,
        residuals = residuals_se_ci$se,
        valid_replicates = valid_reps,
        valid_target_rotations = valid_rot
      ),
      CI = list(
        unrot_loadings = L_unrot_se_ci$ci,
        rot_loadings = L_rot_se_ci$ci,
        Phi = Phi_rot_se_ci$ci,
        Structure = Structure_se_ci$ci,
        fit_indices = gof_se_ci$ci,
        residuals = residuals_se_ci$ci
      ),
      replicates = list(
        unrot_loadings = L_unrot_boot,
        rot_loadings = L_rot_boot,
        Phi = Phi_rot_boot,
        Structure = Structure_boot,
        fit_indices = gof_boot,
        residuals = residuals_boot
      )
    )

  } else if (boot_rot == "orthogonal") {
    colnam_L <- colnames(L_rot)

    L_rot_boot <- array(NA_real_, c(nrow_L, ncol_L, b),
                        dimnames = list(rownam_L, colnam_L,
                                        NULL))

    failed_rot <- 0

    for (boot_i in seq_len(b)) {

      if (failed[boot_i]) next

      # save target-rotated loading matrix
      aligned_i <- tryCatch(
        efa_procrustes(boot_fit[[boot_i]]$unrot_loadings,
                       Target = L_rot, rotation = "orthogonal"),
        error = function(e) NULL
      )
      if (is.null(aligned_i)) {
        failed_rot <- failed_rot + 1
        next
      }
      L_rot_boot[,, boot_i] <- aligned_i$loadings
    }

    valid_rot <- b - n_failed - failed_rot
    if (failed_rot > 0) {
      cli::cli_warn(c("{failed_rot} target rotation{?s} in the bootstrap procedure could not be aligned.",
                    "i" = "Bootstrap SE and CI of rotated loadings are based on {valid_rot} bootstrap sample{?s}."),
                    class = "efa_boot_rotation_failed")
    }

    L_rot_se_ci <- .array_se_ci(L_rot_boot, ps)

    out <- list(
      SE = list(
        unrot_loadings = L_unrot_se_ci$se,
        rot_loadings = L_rot_se_ci$se,
        fit_indices = gof_se_ci$se,
        residuals = residuals_se_ci$se,
        valid_replicates = valid_reps,
        valid_target_rotations = valid_rot
      ),
      CI = list(
        unrot_loadings = L_unrot_se_ci$ci,
        rot_loadings = L_rot_se_ci$ci,
        fit_indices = gof_se_ci$ci,
        residuals = residuals_se_ci$ci
      ),
      replicates = list(
        unrot_loadings = L_unrot_boot,
        rot_loadings = L_rot_boot,
        fit_indices = gof_boot,
        residuals = residuals_boot
      )
    )


  } else {

    out <- list(
      SE = list(
        unrot_loadings = L_unrot_se_ci$se,
        fit_indices = gof_se_ci$se,
        residuals = residuals_se_ci$se,
        valid_replicates = valid_reps
      ),
      CI = list(
        unrot_loadings = L_unrot_se_ci$ci,
        fit_indices = gof_se_ci$ci,
        residuals = residuals_se_ci$ci
      ),
      replicates = list(
        unrot_loadings = L_unrot_boot,
        fit_indices = gof_boot,
        residuals = residuals_boot
      )
    )

  }


  out

}

.array_se_ci <- function(x, probs, M = c(1, 2)) {
  se <- apply(x, M, stats::sd, na.rm = TRUE)
  ci <- lapply(probs, function(p) {
    apply(x, M, stats::quantile, probs = p, na.rm = TRUE)
  })
  names(ci) <- c("lower", "upper")

  list(
    se = se,
    ci = ci
  )
}
