.efa_pooled_analytic_pool <- function(fits,
                                      unrot_loadings_aligned,
                                      align_meta,
                                      ci = 0.95,
                                      align_unrotated = c("signed_tucker_congruence",
                                                          "none", "procrustes"),
                                      rotation_type = c("none", "orthogonal",
                                                        "oblique"),
                                      rot_loadings = NULL,
                                      phis = NULL,
                                      structure_loadings = NULL,
                                      mean_structure_loadings = NULL,
                                      mean_phis = NULL,
                                      h2 = NULL) {
  # Rubin's-rules pool of analytic per-imputation marginal SEs for the unrotated
  # loadings and uniquenesses. Serves component fits run with se = "information":
  # each fit carries $SE$unrot_loadings (p x k) and $SE$uniquenesses (length p),
  # and there are no bootstrap replicate cubes to estimate within-imputation
  # variance from. (The sandwich route is handled separately by
  # .efa_pooled_mi2s() and never reaches this helper.)
  #
  # Under `align_unrotated = "signed_tucker_congruence"` or `"none"`, the
  # per-imputation gauge transform is a signed permutation that preserves
  # marginal SEs up to a column permutation, so $SE$unrot_loadings is reused
  # directly after applying the alignment's `factor_order`.
  #
  # Under `align_unrotated = "procrustes"`, the per-imputation orthogonal Q_d
  # MIXES loading columns, so marginal SEs alone no longer determine the aligned
  # marginal SE matrix. The full unrotated covariance persisted on each fit
  # (`vcov_unrot_loadings`, populated by `se = "information"` and
  # `se = "sandwich"`) is propagated by
  # `.efa_pooled_propagate_procrustes_vcov()`, which implements
  # diag(Q' V_d^(i) Q) per row i. Q_d is treated as fixed.
  #
  # Per-element NA propagation is fail-closed: if any imputation has NA for an
  # element (in either Q_d or SE_d), every pooled output at that position (SE,
  # CI bounds, RIV, FMI, df) is NA.
  #
  # The returned shape mirrors .efa_pooled_bootstrap_pool() so the efa_mi()
  # assembly block can route either pool through the same downstream code:
  # list(SE, CI, replicates = NULL, MI, n_boot = NULL). Empty inputs or fits
  # missing the analytic SE slot return NULL so the caller falls through to the
  # existing silent-downgrade branch. Fits missing the full vcov slot under the
  # procrustes branch abort with a classed condition (fail-closed; the marginal
  # SE alone cannot represent the aligned gauge).
  align_unrotated <- match.arg(align_unrotated)
  rotation_type <- match.arg(rotation_type)
  if (length(fits) == 0L) return(NULL)
  has_se <- vapply(fits, function(x) {
    !is.null(x$SE) &&
      !is.null(x$SE$unrot_loadings) &&
      !is.null(x$SE$uniquenesses)
  }, logical(1L))
  if (!all(has_se)) {
    if (any(has_se)) {
      cli::cli_warn(
        c("Analytic SE slots were found for only some imputations; pooled analytic SEs/CIs were not computed.",
          "i" = "Inspect each fit's {.code $SE} slot to identify the imputations with missing analytic SEs."),
        class = "efa_pooled_analytic_partial_se"
      )
    }
    return(NULL)
  }

  # Whole-fit reliability gate (every alignment method). `EFA()` NA-fills the
  # entire `vcov_unrot_loadings` when the analytic covariance is unreliable (a
  # Heywood case or a singular bordered information matrix), NA-filling the
  # marginal SE matrices in lockstep and firing `efa_se_unreliable` at fit time.
  # Pooling those would yield an all-NA SE bundle with no diagnostic, so abort
  # here; `efa_mi()` converts this into a classed `efa_pooled_se_unavailable`
  # warning and falls back to point-estimate-only pooling. The gate keys on the
  # whole-matrix vcov NA fill, not on individual `$SE` elements, so a genuine
  # per-element NA (one imputation NA at a single cell) still fail-closes per
  # element via the mask below rather than dropping the entire pool.
  vcov_reliable <- vapply(fits, function(x) {
    is.null(x$vcov_unrot_loadings) || !anyNA(x$vcov_unrot_loadings)
  }, logical(1L))
  if (!all(vcov_reliable)) {
    bad <- which(!vcov_reliable)
    n_bad <- length(bad)
    cli::cli_abort(
      c("Pooled analytic standard errors require a reliable unrotated loading covariance on every imputation.",
        "i" = "{cli::qty(n_bad)}Imputation{?s} {.val {bad}} {cli::qty(n_bad)}carr{?ies/y} an NA-filled {.code $vcov_unrot_loadings} (Heywood case or singular information). Drop or re-fit {cli::qty(n_bad)}{?it/them}."),
      class = "efa_pooled_unreliable_vcov"
    )
  }

  if (align_unrotated == "procrustes") {
    # The Procrustes branch propagates the full covariance through Q_d, so the
    # slot must be present (not merely reliable) on every imputation.
    has_vcov <- vapply(fits, function(x) {
      !is.null(x$vcov_unrot_loadings)
    }, logical(1L))
    if (!all(has_vcov)) {
      bad <- which(!has_vcov)
      n_bad <- length(bad)
      cli::cli_abort(
        c("{.code align_unrotated = \"procrustes\"} with {.code se = \"information\"} requires the full unrotated loading covariance on every imputation.",
          "i" = "{cli::qty(n_bad)}Imputation{?s} {.val {bad}} {cli::qty(n_bad)}{?is/are} missing {.code $vcov_unrot_loadings}; re-fit those imputations or use {.code align_unrotated = \"signed_tucker_congruence\"}."),
        class = "efa_pooled_no_vcov"
      )
    }
  }

  m <- length(fits)
  L_target <- as.matrix(unrot_loadings_aligned[[1]])
  p_vars   <- nrow(L_target)
  k_fac    <- ncol(L_target)
  loading_dn <- dimnames(L_target)
  var_names  <- rownames(L_target)

  # Per-imputation aligned loadings and aligned SEs.
  # - signed_tucker_congruence / none: column-permute the marginal SE matrix by
  #   the alignment's factor_order.
  # - procrustes: propagate the full vcov through Q_d row-by-row to obtain the
  #   aligned marginal variances, then take sqrt to land back on the same
  #   variance-of-vec(L_d Q_d) scale as Ubar_L (which squares SE_mat).
  Q_mat   <- matrix(NA_real_, nrow = m, ncol = p_vars * k_fac)
  SE_mat  <- matrix(NA_real_, nrow = m, ncol = p_vars * k_fac)
  Q_psi   <- matrix(NA_real_, nrow = m, ncol = p_vars)
  SE_psi_mat <- matrix(NA_real_, nrow = m, ncol = p_vars)
  # Imputations whose covariance block the canonical-gauge branch needed but did
  # not find; collected so the silent all-NA block below can be reported.
  gauge_no_vcov <- integer(0)
  for (d in seq_len(m)) {
    Q_mat[d, ] <- as.vector(as.matrix(unrot_loadings_aligned[[d]]))

    if (identical(align_meta[[d]]$type, "procrustes")) {
      V_d  <- fits[[d]]$vcov_unrot_loadings
      Q_d  <- align_meta[[d]]$Q
      se_d <- .efa_pooled_propagate_procrustes_vcov(V_d, Q_d, p_vars, k_fac)
      # Reading the covariance block bypasses the marginal-SE NA mask, so apply
      # it here too: an SE the component fit could not certify must not come
      # back finite. Q_d mixes columns, so one NA contaminates its whole row.
      se_marg <- fits[[d]]$SE$unrot_loadings
      if (!is.null(se_marg)) {
        se_d[apply(is.na(as.matrix(se_marg)), 1, any), ] <- NA_real_
      }
      SE_mat[d, ] <- as.vector(se_d)
    } else if (!is.null(align_meta[[d]]$C)) {
      # The canonical-gauge rotation mixes loading columns, so the marginal SEs
      # cannot merely be permuted: the full unrotated covariance is propagated
      # through the composed orthogonal transform E_d C, where E_d is the signed
      # permutation that aligned imputation d. Fail-closed at the imputation
      # level -- a missing or non-finite covariance block leaves this row NA and
      # the per-element NA mask blanks the affected pooled outputs -- rather than
      # aborting, so a component fit without the block still pools everything else.
      fo <- align_meta[[d]]$factor_order
      sg <- align_meta[[d]]$factor_sign
      if (is.null(fo) || is.null(sg)) {
        cli::cli_abort(
          "Per-imputation alignment metadata is missing for imputation {d}; analytic-SE pooling requires a column-permutation/sign alignment.",
          class = "efa_pooled_analytic_align_meta_missing"
        )
      }
      E <- matrix(0, k_fac, k_fac)
      E[cbind(fo, seq_len(k_fac))] <- sg
      V_d <- fits[[d]]$vcov_unrot_loadings
      if (is.null(V_d)) {
        gauge_no_vcov <- c(gauge_no_vcov, d)
      } else {
        se_aligned <- .efa_pooled_propagate_procrustes_vcov(
          as.matrix(V_d), E %*% align_meta[[d]]$C, p_vars, k_fac)
        # Keep the marginal-SE NA mask meaningful. The propagation reads the
        # covariance block, so an unreliable marginal SE would otherwise be
        # silently replaced by a finite value. Because the gauge rotation mixes
        # columns, one NA contaminates its whole row.
        se_marg <- fits[[d]]$SE$unrot_loadings
        if (!is.null(se_marg)) {
          bad_row <- apply(is.na(as.matrix(se_marg)), 1, any)
          se_aligned[bad_row, ] <- NA_real_
        }
        SE_mat[d, ] <- as.vector(se_aligned)
      }
    } else {
      se_d <- as.matrix(fits[[d]]$SE$unrot_loadings)
      fo <- align_meta[[d]]$factor_order
      if (is.null(fo)) {
        # All currently-supported producers of `align_meta` populate
        # `factor_order` (the `none` and `signed_tucker_congruence` branches of
        # `.efa_pooled_align_unrotated_list`); procrustes meta is routed into
        # the V_d / Q_d branch above by the `type` test and never reaches
        # here. Catch any future bypass at the helper boundary rather than
        # silently falling back to an identity permutation that would pool
        # SEs against unaligned loadings.
        cli::cli_abort(
          "Per-imputation alignment metadata is missing for imputation {d}; analytic-SE pooling requires a column-permutation/sign alignment.",
          class = "efa_pooled_analytic_align_meta_missing"
        )
      }
      SE_mat[d, ] <- as.vector(se_d[, fo, drop = FALSE])
    }

    psi_d <- .efa_pooled_communalities(unrot_loadings_aligned[[d]], NULL)
    Q_psi[d, ] <- 1 - as.numeric(psi_d)
    SE_psi_mat[d, ] <- as.numeric(fits[[d]]$SE$uniquenesses)
  }

  if (length(gauge_no_vcov)) {
    # Without the covariance block the gauge rotation cannot be propagated, so
    # those rows stay NA and the per-element mask blanks the whole unrotated
    # loading family. Say so: an all-NA block returned in silence is not
    # diagnosable, and every sibling branch here signals its downgrade.
    cli::cli_warn(
      c("The unrotated loading covariance is missing for {cli::qty(length(gauge_no_vcov))}imputation{?s} {.val {gauge_no_vcov}}; pooled unrotated loading SEs/CIs were not computed.",
        "i" = "The pooled solution is returned in the extraction's canonical gauge, which mixes loading columns, so the full {.code $vcov_unrot_loadings} is required rather than the marginal SEs alone."),
      class = "efa_pooled_gauge_vcov_missing"
    )
  }

  # A component fit whose analytic covariance was unreliable at the loading level carries
  # an NA-filled $SE$unrot_loadings while keeping a finite $vcov_unrot_loadings, so the
  # whole-fit gate above passes and the fail-closed mask then blanks the entire pooled
  # unrotated SE block. Name the affected imputations, as the rotated family does: an
  # all-NA block returned in silence reads as a computation that succeeded. The
  # imputations already reported by the missing-covariance warning are excluded so a
  # single downgrade is not announced twice.
  unrot_se_unreliable <- setdiff(
    which(apply(SE_mat, 1L, function(r) all(is.na(r)))),
    gauge_no_vcov
  )
  if (length(unrot_se_unreliable) > 0L) {
    n_bad <- length(unrot_se_unreliable)
    cli::cli_warn(
      c("Pooled unrotated standard errors could not be produced for {n_bad} imputation{?s}.",
        "i" = "Affected: {.val {unrot_se_unreliable}} (each left an NA-filled unrotated standard error, typically the upstream {.code efa_se_unreliable} signal).",
        "i" = "The pooled unrotated loading standard errors are returned as {.code NA}."),
      class = "efa_pooled_unrotated_se_unreliable"
    )
  }

  # Fail-closed NA mask and Rubin pooling are shared with the rotated families
  # below via .efa_pooled_analytic_marginal().
  alpha <- 1 - ci
  pool_L   <- .efa_pooled_analytic_marginal(Q_mat, SE_mat,     m, alpha)
  pool_psi <- .efa_pooled_analytic_marginal(Q_psi, SE_psi_mat, m, alpha)

  reshape_loading <- function(v) {
    M <- matrix(v, p_vars, k_fac)
    dimnames(M) <- loading_dn
    M
  }
  reshape_psi <- function(v) {
    out <- as.numeric(v)
    if (!is.null(var_names)) names(out) <- var_names
    out
  }

  unrot_asm <- .efa_pooled_assemble_family(pool_L,   reshape_loading)
  psi_asm   <- .efa_pooled_assemble_family(pool_psi, reshape_psi)

  SE <- list(unrot_loadings = unrot_asm$SE, uniquenesses = psi_asm$SE)
  CI <- list(unrot_loadings = unrot_asm$CI, uniquenesses = psi_asm$CI)
  MI <- list(unrot_loadings = unrot_asm$MI, uniquenesses = psi_asm$MI)

  ## Rotated quantities -------------------------------------------------------
  ## Pool the rotated loadings, communalities, and (oblique) factor correlations
  ## and structure coefficients, each in a common MI rotational gauge. A
  ## rotated-loading standard error is conditional on the rotation criterion
  ## (Archer & Jennrich 1973; Jennrich 1973, 1974; Zhang, Preacher, & Jennrich
  ## 2012; Zhang & Preacher 2015), so for orthogonal
  ## and oblique rotations alike the within-imputation variance is each fit's own
  ## criterion-aware delta-method rotated SE (the quantity EFA() returns), reused
  ## after a signed-permutation alignment to the MI target. The between-imputation
  ## spread of the aligned point estimates supplies B, so the pooled SE shares the
  ## estimand of the np-boot pool, which likewise aligns each replicate to the
  ## target. This is a deliberate approximation -- the per-fit SE is conditional on
  ## that fit's rotation optimum rather than the common gauge -- and is flagged in
  ## MI$<param>$method = "signed_permutation_approx". Communalities are
  ## rotation-invariant and pool without any alignment.
  if (rotation_type != "none") {
    # Phi and Structure exist only for a genuine oblique rotation (a single-factor
    # solution carries no factor correlations); gate on k_fac so the oblique branch
    # never dereferences the NULL Phi/Structure slots of a one-factor fit.
    oblique <- rotation_type == "oblique" && k_fac >= 2L

    q_rot  <- matrix(NA_real_, m, p_vars * k_fac)
    se_rot <- matrix(NA_real_, m, p_vars * k_fac)
    q_h2   <- matrix(NA_real_, m, p_vars)
    se_h2  <- matrix(NA_real_, m, p_vars)
    if (oblique) {
      vlen   <- k_fac * (k_fac + 1L) / 2L
      q_phi  <- matrix(NA_real_, m, vlen)
      se_phi <- matrix(NA_real_, m, vlen)
      q_str  <- matrix(NA_real_, m, p_vars * k_fac)
      se_str <- matrix(NA_real_, m, p_vars * k_fac)
    }

    for (d in seq_len(m)) {
      Lr_d  <- as.matrix(rot_loadings[[d]])
      Phi_d <- if (oblique) as.matrix(phis[[d]]) else NULL
      q_rot[d, ] <- .efa_pooled_vec(Lr_d)
      q_h2[d, ]  <- .efa_pooled_communalities(Lr_d, Phi_d)
      se_h2[d, ] <- as.numeric(fits[[d]]$SE$communalities)

      # Within-imputation variance for the rotated loadings = this fit's
      # criterion-aware delta-method rotated SE, reused after a signed-permutation
      # alignment to THIS imputation's pooled rotated loading Lr_d. Aligning to
      # Lr_d (rather than the global target) keeps the SE column order matched to
      # the order q_rot/q_phi/q_str use; aligning to the global target can pick a
      # different order when factors are weakly separated, pairing a cell's pooled
      # estimate with another factor's within-imputation variance. Only the column
      # order is used (SE magnitudes are sign-invariant).
      fo <- .align_solution(
        L_target = Lr_d,
        L        = as.matrix(fits[[d]]$rot_loadings)
      )$factor_order
      se_rot[d, ] <- .efa_pooled_vec(
        as.matrix(fits[[d]]$SE$rot_loadings)[, fo, drop = FALSE]
      )

      if (oblique) {
        # Apply the same column alignment to the per-fit Phi and Structure SEs.
        se_phi[d, ] <- .vech(
          as.matrix(fits[[d]]$SE$Phi)[fo, fo, drop = FALSE]
        )
        q_phi[d, ]  <- .vech(Phi_d)
        se_str[d, ] <- .efa_pooled_vec(
          as.matrix(fits[[d]]$SE$Structure)[, fo, drop = FALSE]
        )
        q_str[d, ]  <- .efa_pooled_vec(as.matrix(structure_loadings[[d]]))
      }
    }

    # A component fit whose rotated-SE Jacobian could not be reproduced carries an
    # all-NA $SE$rot_loadings (the upstream efa_se_unreliable signal) even when its
    # unrotated covariance is fine, so the whole-fit reliability gate above does not
    # catch it. Such a fit leaves an all-NA row in se_rot, which the fail-closed mask
    # in .efa_pooled_analytic_marginal() then propagates to every cell of the pooled
    # rotated-loading (and, for oblique fits, Phi/Structure) SE family. Name the
    # affected imputations so the resulting all-NA rotated SE block is diagnosable
    # rather than silent.
    rot_se_unreliable <- which(apply(se_rot, 1L, function(r) all(is.na(r))))
    if (length(rot_se_unreliable) > 0L) {
      n_bad <- length(rot_se_unreliable)
      affected <- if (oblique) {
        "rotated-loading, factor-correlation, and structure-coefficient"
      } else {
        "rotated-loading"
      }
      cli::cli_warn(
        c("Pooled rotated standard errors could not be produced for {n_bad} imputation{?s}.",
          "i" = "Affected: {.val {rot_se_unreliable}} (each carried an NA-filled rotated standard error, the upstream {.code efa_se_unreliable} signal).",
          "i" = "The pooled {affected} standard errors are returned as {.code NA}."),
        class = "efa_pooled_rotated_se_unreliable"
      )
    }

    pool_rot <- .efa_pooled_analytic_marginal(q_rot, se_rot, m, alpha)
    pool_h2  <- .efa_pooled_analytic_marginal(q_h2, se_h2, m, alpha,
                                              est_override = as.numeric(h2))

    # Rotated loadings and structure are p x k like the unrotated loadings, so they
    # share the same factor-column labels (reshape_loading carries loading_dn).
    rot_asm <- .efa_pooled_assemble_family(
      pool_rot, reshape_loading,
      method = "signed_permutation_approx"
    )
    SE$rot_loadings <- rot_asm$SE; CI$rot_loadings <- rot_asm$CI
    MI$rot_loadings <- rot_asm$MI

    h2_asm <- .efa_pooled_assemble_family(pool_h2, reshape_psi,
                                          method = "gauge_invariant")
    SE$communalities <- h2_asm$SE; CI$communalities <- h2_asm$CI
    MI$communalities <- h2_asm$MI
    # `communalities` is the name EFA() uses, so it is the canonical one here too;
    # `h2` is what this route has always returned and is kept as an alias of it. The
    # MI diagnostics carry no alias: the printed FMI/RIV summary walks the whole MI
    # tree, so a duplicated family would enter those statistics twice.
    SE$h2 <- SE$communalities; CI$h2 <- CI$communalities

    if (oblique) {
      phi_dn <- dimnames(as.matrix(phis[[1]]))
      unvech_phi <- function(v) {
        M <- .unvech(v, k_fac)
        dimnames(M) <- phi_dn
        M
      }
      # Center the Phi CI on the same symmetrised pooled estimate reported as $Phi
      # (mean_phis, always present on this oblique branch), matching the est_override
      # the h2 and Structure families use; without it the off-diagonal CI would
      # center on the unsymmetrised column mean of vech(Phi), which can differ from
      # the reported (symmetrised) estimate.
      pool_phi <- .efa_pooled_analytic_marginal(
        q_phi, se_phi, m, alpha,
        est_override = .vech(mean_phis)
      )
      phi_asm <- .efa_pooled_assemble_family(pool_phi, unvech_phi,
                                             method = "signed_permutation_approx")
      # The fixed unit diagonal of Phi is not an estimated parameter: zero its SE
      # and NA its MI diagnostics so the fixed cells stay out of the printed
      # FMI/RIV summaries.
      diag(phi_asm$SE) <- 0
      diag(phi_asm$MI$RIV) <- NA_real_
      diag(phi_asm$MI$FMI) <- NA_real_
      diag(phi_asm$MI$df)  <- NA_real_
      SE$Phi <- phi_asm$SE; CI$Phi <- phi_asm$CI; MI$Phi <- phi_asm$MI

      pool_str <- .efa_pooled_analytic_marginal(
        q_str, se_str, m, alpha,
        est_override = .efa_pooled_vec(mean_structure_loadings)
      )
      str_asm <- .efa_pooled_assemble_family(pool_str, reshape_loading,
                                             method = "signed_permutation_approx")
      SE$Structure <- str_asm$SE; CI$Structure <- str_asm$CI
      MI$Structure <- str_asm$MI
    }
  }

  list(SE = SE, CI = CI, replicates = NULL, MI = MI, n_boot = NULL)
}


.efa_pooled_align_unrot_boot <- function(L_boot, target, align_unrotated) {
  # Align one bootstrap unrotated loading matrix to the pooled unrotated target.
  if (align_unrotated == "none") {
    return(as.matrix(L_boot))
  }

  if (align_unrotated == "procrustes") {
    return(.change_class(efa_procrustes(
      A = as.matrix(L_boot),
      Target = as.matrix(target),
      rotation = "orthogonal"
    )$loadings, "matrix"))
  }

  .align_solution(
    L_target = as.matrix(target),
    L = as.matrix(L_boot)
  )$loadings
}

.efa_pooled_align_rot_boot <- function(L_boot_unrot, final_target, rotation_type,
                                       procrustes_args = list()) {
  # Re-rotate/re-align the bootstrap unrotated loadings to the final MI target.
  # This is the key step that makes within-imputation bootstrap variances refer
  # to the same estimand as the pooled rotated solution.
  if (rotation_type == "none") {
    return(NULL)
  }

  pr <- do.call(
    efa_procrustes,
    c(list(A = as.matrix(L_boot_unrot),
           Target = as.matrix(final_target),
           rotation = rotation_type),
      procrustes_args)
  )

  # Drop the replicate only when no valid alignment was produced. A best
  # multi-start fit that did not formally converge is still the lowest-objective
  # alignment available, with well-defined loadings, and is retained.
  if (isFALSE(pr$valid)) {
    return(NULL)
  }

  list(
    Lambda = .change_class(pr$loadings, "matrix"),
    Phi = if (!is.null(pr$Phi)) as.matrix(pr$Phi) else NULL
  )
}

.efa_pooled_bootstrap_pool <- function(fits,
                                       orig_R_list,
                                       unrot_loadings_aligned,
                                       mean_unrot_loadings,
                                       rot_loadings = NULL,
                                       phis = NULL,
                                       structure_loadings = NULL,
                                       mean_structure_loadings = NULL,
                                       final_target = NULL,
                                       rotation_type = c("none", "orthogonal", "oblique"),
                                       align_unrotated = c("signed_tucker_congruence", "none", "procrustes"),
                                       procrustes_args = list(),
                                       h2,
                                       residuals,
                                       alpha = 0.05) {

  rotation_type <- match.arg(rotation_type)
  align_unrotated <- match.arg(align_unrotated)

  has_boot <- .efa_pooled_has_replicates(fits)
  if (!all(has_boot)) {
    if (any(has_boot)) {
      cli::cli_warn("Bootstrap replicates were found for only some imputations; pooled bootstrap SEs/CIs were not computed.",
                    class = "efa_pooled_partial_boot")
    }
    return(NULL)
  }

  m <- length(fits)
  # Coerce the loop-invariant alignment targets once; mean_unrot_loadings is a
  # LOADINGS object whose as.matrix() coercion is otherwise repeated on every
  # bootstrap replicate inside the d/b loops.
  mean_unrot_target <- as.matrix(mean_unrot_loadings)
  final_target_m <- if (is.null(final_target)) NULL else as.matrix(final_target)
  p_vars <- nrow(mean_unrot_target)
  k <- ncol(mean_unrot_target)
  loading_dimnames <- dimnames(mean_unrot_target)

  # A single factor cannot be rotated and has no factor correlations, so an
  # oblique request on one factor degenerates to the no-Phi case (matching a
  # single-fit EFA()); align such a solution by sign only.
  oblique <- rotation_type == "oblique" && k >= 2L
  boot_proc_rotation <- if (oblique) "oblique" else "orthogonal"

  B_vec <- vapply(fits, function(x) dim(x$replicates$unrot_loadings)[3], integer(1))
  if (length(unique(B_vec)) > 1L) {
    B_use <- min(B_vec)
    cli::cli_warn("The number of bootstrap replicates differs across imputations; using the minimum number available in all imputations.",
                  class = "efa_pooled_unequal_boot")
  } else {
    B_use <- B_vec[[1]]
  }

  if (!is.finite(B_use) || B_use < 2L) {
    cli::cli_warn("At least two bootstrap replicates per imputation are required for pooled bootstrap SEs/CIs.",
                  class = "efa_pooled_min_boot")
    return(NULL)
  }

  q_unrot <- lapply(unrot_loadings_aligned, .efa_pooled_vec)
  boot_unrot <- vector("list", m)

  q_rot <- boot_rot <- NULL
  q_phi <- boot_phi <- NULL
  q_structure <- boot_structure <- NULL
  q_h2 <- vector("list", m)
  boot_h2 <- vector("list", m)
  q_residuals <- vector("list", m)
  boot_residuals <- vector("list", m)

  if (rotation_type != "none") {
    q_rot <- lapply(rot_loadings, .efa_pooled_vec)
    boot_rot <- vector("list", m)
  }
  if (oblique) {
    q_phi <- lapply(phis, .vech)
    boot_phi <- vector("list", m)
    q_structure <- lapply(structure_loadings, .efa_pooled_vec)
    boot_structure <- vector("list", m)
  }

  nonconv_procrustes <- integer(m)
  boot_failures <- integer(m)

  for (d in seq_len(m)) {
    arr <- fits[[d]]$replicates
    unrot_arr <- arr$unrot_loadings

    boot_unrot[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars * k)

    if (rotation_type != "none") {
      boot_rot[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars * k)
    }
    if (oblique) {
      boot_phi[[d]] <- matrix(NA_real_, nrow = B_use, ncol = k * (k + 1) / 2)
      boot_structure[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars * k)
    }

    boot_h2[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars)
    boot_residuals[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars * p_vars)

    if (oblique) {
      q_h2[[d]] <- .efa_pooled_communalities(rot_loadings[[d]], phis[[d]])
      q_residuals[[d]] <- .efa_pooled_vec(.efa_pooled_residual_from_solution(
        orig_R_list[[d]], rot_loadings[[d]], phis[[d]]
      ))
    } else if (rotation_type != "none") {
      q_h2[[d]] <- .efa_pooled_communalities(rot_loadings[[d]], NULL)
      q_residuals[[d]] <- .efa_pooled_vec(.efa_pooled_residual_from_solution(
        orig_R_list[[d]], rot_loadings[[d]], NULL
      ))
    } else {
      q_h2[[d]] <- .efa_pooled_communalities(unrot_loadings_aligned[[d]], NULL)
      q_residuals[[d]] <- .efa_pooled_vec(.efa_pooled_residual_from_solution(
        orig_R_list[[d]], unrot_loadings_aligned[[d]], NULL
      ))
    }

    for (b in seq_len(B_use)) {
      Lb_unrot0 <- unrot_arr[, , b]

      # A component EFA NA-fills a bootstrap replicate it could not fit; skip it
      # (its pre-initialized rows stay NA and are dropped by the NA-tolerant
      # pooling) and tally it, rather than letting the NA reach the alignment /
      # Procrustes step and abort the whole pooled bootstrap.
      if (!all(is.finite(Lb_unrot0))) {
        boot_failures[d] <- boot_failures[d] + 1L
        next
      }

      Lb_unrot <- .efa_pooled_align_unrot_boot(
        L_boot = Lb_unrot0,
        target = mean_unrot_target,
        align_unrotated = align_unrotated
      )
      boot_unrot[[d]][b, ] <- .efa_pooled_vec(Lb_unrot)

      if (rotation_type != "none") {
        # Do not use arr$rot_loadings here. Those matrices were aligned to the
        # imputation-specific EFA target. MI pooling requires every bootstrap
        # replicate to be represented in the final MI target space.
        rot_b <- .efa_pooled_align_rot_boot(
          L_boot_unrot = Lb_unrot0,
          final_target = final_target_m,
          rotation_type = boot_proc_rotation,
          procrustes_args = procrustes_args
        )

        if (!is.null(rot_b)) {
          boot_rot[[d]][b, ] <- .efa_pooled_vec(rot_b$Lambda)

          if (oblique) {
            boot_phi[[d]][b, ] <- .vech(rot_b$Phi)
            boot_structure[[d]][b, ] <- .efa_pooled_vec(rot_b$Lambda %*% rot_b$Phi)
            boot_h2[[d]][b, ] <- .efa_pooled_communalities(rot_b$Lambda, rot_b$Phi)
          } else {
            boot_h2[[d]][b, ] <- .efa_pooled_communalities(rot_b$Lambda, NULL)
          }
        } else {
          nonconv_procrustes[d] <- nonconv_procrustes[d] + 1L
        }
      } else {
        boot_h2[[d]][b, ] <- .efa_pooled_communalities(Lb_unrot, NULL)
      }

      if (!is.null(arr$residuals)) {
        boot_residuals[[d]][b, ] <- .efa_pooled_vec(arr$residuals[, , b])
      }
    }
  }

  # Insufficient implies some failed, so test it first and emit a single warning.
  if (any(B_use - boot_failures < 2L)) {
    cli::cli_warn("Too few valid bootstrap replicates remain after dropping failures; pooled bootstrap SEs/CIs could not be computed.",
                  class = "efa_pooled_boot_insufficient")
    return(NULL)
  }

  if (any(boot_failures > 0L)) {
    cli::cli_warn("Some bootstrap replicates failed in the component EFAs; pooled bootstrap SEs/CIs are based on the valid replicates.",
                  class = "efa_pooled_boot_failed")
  }

  if (any(nonconv_procrustes > 0L)) {
    cli::cli_warn("Some bootstrap Procrustes rotations did not converge; pooled bootstrap SEs/CIs are based on the valid aligned replicates.",
                  class = "efa_pooled_boot_nonconv")
  }

  pool_unrot <- .efa_pooled_rubin_pool(q_unrot, boot_unrot, alpha = alpha)
  pool_h2 <- .efa_pooled_rubin_pool(
    q_h2, boot_h2, alpha = alpha,
    est_override = as.vector(h2)
  )
  pool_residuals <- .efa_pooled_rubin_pool(
    q_residuals, boot_residuals, alpha = alpha,
    est_override = .efa_pooled_vec(residuals)
  )

  # Reshapers for the two rectangular families pooled below (loading-shaped and
  # correlation-matrix-shaped); Phi is expanded from its vech instead.
  as_loading <- \(v) matrix(v, nrow = p_vars, ncol = k)

  unrot_res <- .efa_pooled_rubin_result(pool_unrot, as_loading, loading_dimnames)
  residual_res <- .efa_pooled_rubin_result(
    pool_residuals, \(v) matrix(v, nrow = p_vars, ncol = p_vars),
    dimnames(as.matrix(residuals))
  )

  SE <- list(
    communalities = pool_h2$se,
    unrot_loadings = unrot_res$se,
    residuals = residual_res$se
  )

  CI <- list(
    communalities = list(lower = pool_h2$ci$lower, upper = pool_h2$ci$upper),
    unrot_loadings = unrot_res$ci,
    residuals = residual_res$ci
  )

  arrays <- list(
    unrot_loadings = boot_unrot,
    h2 = boot_h2,
    residuals = boot_residuals
  )

  # Same canonical-name/alias split as the analytic route above.
  SE$h2 <- SE$communalities; CI$h2 <- CI$communalities

  MI <- list(
    unrot_loadings = list(RIV = pool_unrot$RIV, FMI = pool_unrot$FMI, df = pool_unrot$df),
    communalities = list(RIV = pool_h2$RIV, FMI = pool_h2$FMI, df = pool_h2$df),
    residuals = list(RIV = pool_residuals$RIV, FMI = pool_residuals$FMI, df = pool_residuals$df),
    # Per-imputation alignment outcomes. `boot_failures` are replicates the
    # component EFA could not fit (NA-filled, dropped before rotation);
    # `bootstrap_rotation_failures` are replicates whose Procrustes alignment to
    # the MI target did not produce a valid rotation. The valid count subtracts
    # both, so it never overstates the replicates that actually entered the pool.
    bootstrap_source_failures = boot_failures,
    bootstrap_rotation_failures = nonconv_procrustes,
    bootstrap_rotation_valid = B_use - boot_failures - nonconv_procrustes
  )

  if (rotation_type != "none") {
    pool_rot <- .efa_pooled_rubin_pool(q_rot, boot_rot, alpha = alpha)
    rot_res <- .efa_pooled_rubin_result(pool_rot, as_loading, loading_dimnames)
    SE$rot_loadings <- rot_res$se
    CI$rot_loadings <- rot_res$ci
    arrays$rot_loadings <- boot_rot
    MI$rot_loadings <- list(RIV = pool_rot$RIV, FMI = pool_rot$FMI, df = pool_rot$df)
  }

  if (oblique) {
    pool_phi <- .efa_pooled_rubin_pool(q_phi, boot_phi, alpha = alpha)
    phi_res <- .efa_pooled_rubin_result(
      pool_phi, \(v) .unvech(v, k), dimnames(as.matrix(phis[[1]]))
    )
    SE$Phi <- phi_res$se
    CI$Phi <- phi_res$ci
    arrays$Phi <- boot_phi
    MI$Phi <- list(RIV = pool_phi$RIV, FMI = pool_phi$FMI, df = pool_phi$df)

    pool_structure <- .efa_pooled_rubin_pool(
      q_structure, boot_structure, alpha = alpha,
      est_override = .efa_pooled_vec(mean_structure_loadings)
    )
    structure_res <- .efa_pooled_rubin_result(
      pool_structure, as_loading, loading_dimnames
    )
    SE$Structure <- structure_res$se
    CI$Structure <- structure_res$ci
    arrays$Structure <- boot_structure
    MI$Structure <- list(
      RIV = pool_structure$RIV,
      FMI = pool_structure$FMI,
      df = pool_structure$df
    )
  }

  ## Fit-index summaries ------------------------------------------------------
  ## fit_indices_descriptive: Rubin-Wald MI summaries of the per-imputation
  ## bootstrap fit indices. Chi-square-type fit is pooled with the D2 rule and its
  ## reference distribution supplies the point-estimate uncertainty, so no
  ## bootstrap-percentile interval is constructed for the pooled fit statistic.

  fit_q_names <- Reduce(intersect, lapply(fits, function(x) names(x$fit_indices)))
  fit_q_names <- fit_q_names[vapply(fits[[1]]$fit_indices[fit_q_names], function(x) {
    is.numeric(x) && length(x) == 1L
  }, logical(1))]

  # The bootstrap does not aggregate every fit index -- the RMSEA confidence bounds have no
  # per-replicate value, because the replicate fits run with ci = FALSE -- so the replicate
  # matrix is narrower than the point-estimate fit-index list. Take the pooled quantities from
  # the columns each fit's replicate matrix actually carries, rather than assuming the two line
  # up.
  boot_fit_arrs <- lapply(fits, function(x) {
    fit_arr <- x$replicates$fit_indices
    if (is.null(fit_arr)) return(NULL)
    fit_arr <- as.matrix(fit_arr)
    if (is.null(colnames(fit_arr))) {
      # An object saved before the replicate matrix was labelled carries no dimnames; there the
      # fit-index names are the only available labels, and only if they match by position.
      nms <- names(x$fit_indices)
      if (length(nms) != ncol(fit_arr)) return(NULL)
      colnames(fit_arr) <- nms
    }
    fit_arr
  })
  fit_q_names <- intersect(
    fit_q_names,
    Reduce(intersect, lapply(boot_fit_arrs, function(a) {
      if (is.null(a)) character(0) else colnames(a)
    }))
  )

  can_describe_fit_boot <- length(fit_q_names) > 0L &&
    all(vapply(boot_fit_arrs, function(a) {
      !is.null(a) && nrow(a) >= B_use
    }, logical(1)))

  if (can_describe_fit_boot) {
    q_fit <- lapply(fits, function(x) unlist(x$fit_indices[fit_q_names]))
    boot_fit <- lapply(boot_fit_arrs, function(a) {
      a[seq_len(B_use), fit_q_names, drop = FALSE]
    })

    pool_fit_desc <- .efa_pooled_rubin_pool(q_fit, boot_fit, alpha = alpha)
    SE$fit_indices_descriptive <- stats::setNames(pool_fit_desc$se, fit_q_names)
    CI$fit_indices_descriptive <- list(
      lower = stats::setNames(pool_fit_desc$ci$lower, fit_q_names),
      upper = stats::setNames(pool_fit_desc$ci$upper, fit_q_names)
    )
    arrays$fit_indices_descriptive <- boot_fit
    MI$fit_indices_descriptive <- list(
      RIV = stats::setNames(pool_fit_desc$RIV, fit_q_names),
      FMI = stats::setNames(pool_fit_desc$FMI, fit_q_names),
      df = stats::setNames(pool_fit_desc$df, fit_q_names)
    )
  }

  list(
    SE = SE,
    CI = CI,
    replicates = arrays,
    MI = MI,
    n_boot = B_use
  )
}
