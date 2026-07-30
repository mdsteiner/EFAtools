.efa_pooled_medoid_anchor <- function(unrot_loadings) {
  # Index of the imputation that is closest, in aligned squared Frobenius
  # distance, to all the others. Anchoring the alignment there instead of on
  # whichever imputation happens to be first makes the pooled solution invariant
  # to the order of `data_list`: the medoid is a property of the set. Cost is
  # O(D^2) cheap column matchings.
  #
  # Aligning to a single arbitrarily chosen imputation is the usual approach --
  # it is what the reference-target route below does, and what psych's
  # fa.pooled() does -- but it makes the result depend on which imputation that
  # is. The alternative in the literature is to drop the reference and iterate a
  # Generalized Procrustes consensus instead (Gower 1975; van Ginkel &
  # Kroonenberg 2014); choosing a central reference keeps the pooled estimate a
  # plain average of matched solutions, which the consensus centroid is not.
  D <- length(unrot_loadings)
  if (D < 2L) return(1L)
  Ls <- lapply(unrot_loadings, as.matrix)
  # The distance is symmetric -- the column matching of (b, a) is the inverse of
  # the matching of (a, b), and an orthogonal factor leaves the norm alone -- so
  # only the upper triangle needs matching, halving the O(D^2) work.
  dm <- matrix(0, D, D)
  for (a in seq_len(D - 1L)) {
    A <- Ls[[a]]
    for (b in seq.int(a + 1L, D)) {
      v <- sum((unclass(.align_solution(L_target = A, L = Ls[[b]])$loadings) - A)^2)
      dm[a, b] <- v
      dm[b, a] <- v
    }
  }
  tot <- rowSums(dm)

  # Ties are not a numerical edge case here: with two imputations every
  # candidate is tied by construction, and identical imputations tie for any D.
  # Breaking them on list position would put back exactly the order dependence
  # the medoid removes, so break them on the candidates' own content -- the
  # squared Frobenius norm, which no alignment can change and which does not
  # depend on where the imputation sits in the list. A candidate set that is
  # still tied after that falls back to the lowest position, which is order
  # dependent -- it needs two imputations with the same total distance *and* the
  # same norm, which outside of exactly duplicated solutions does not happen.
  scale <- max(1, max(abs(tot)))
  cand <- which(tot <= min(tot) + 1e-10 * scale)
  if (length(cand) == 1L) return(cand)
  ss <- vapply(cand, function(i) sum(Ls[[i]]^2), numeric(1))
  cand[which.max(round(ss, 10))]
}

.efa_pooled_canonical_gauge <- function(M, anchor) {
  # The extraction identifies the unrotated solution by making L' Psi^-1 L
  # diagonal with descending diagonal (Lawley & Maxwell 1971), so a single fit's
  # unrotated loadings are fixed up to column sign. Averaging aligned
  # per-imputation solutions does not preserve that constraint, which leaves the
  # pooled matrix in a gauge no component fit uses and makes it incomparable
  # element-by-element with any single-fit unrotated solution. Returns the k x k
  # orthogonal C that puts M C back in the canonical gauge. Psi is taken from the
  # pooled matrix itself, so the gauge is consistent with the loadings actually
  # reported.
  #
  # C is common to every imputation and orthogonal, so it leaves all
  # gauge-invariant quantities (communalities, total variance accounted for, the
  # model-implied correlation matrix and hence the residuals and RMSR) exactly
  # unchanged; only the split of variance across factors moves.
  k <- ncol(M)
  # Which constraint the extraction uses depends on the estimator: principal-axis
  # extractions diagonalise L'L, maximum likelihood diagonalises L' Psi^-1 L.
  # Rather than assume one, read it off `anchor`, which is an untouched
  # single-fit solution and therefore already satisfies whichever holds.
  gram  <- function(L) crossprod(L)
  # Psi is reconstructed from the loadings, and it has to be floored at the same
  # value the fitting routines use. A maximum likelihood solution with a
  # uniqueness pinned at the boundary diagonalises L' Psi^-1 L against the
  # *floored* psi, whereas 1 - rowSums(L^2) comes back a whisker below it; the
  # inverse multiplies that tiny gap by ~1/floor and turns it into a defect large
  # enough to fail the tolerance below, which would silently switch the
  # re-gauging off exactly for the boundary solutions. Flooring here reproduces
  # the psi the extraction actually used.
  wgram <- function(L) crossprod(L, L / pmax(1 - rowSums(L^2), .uniqueness_floor))
  defect <- function(A) {
    dg <- sum(abs(diag(A)))
    if (!is.finite(dg) || dg <= 0) return(Inf)
    sum(abs(A[upper.tri(A)])) / dg
  }
  d_g <- defect(gram(anchor))
  d_w <- defect(wgram(anchor))
  # Loose enough to clear the extraction's own convergence noise, tight enough
  # that the constraint the estimator does *not* use does not pass. The margin is
  # narrower than it looks: over a sweep of 492 fits (PAF/ML/ULS/MINRES/DWLS x
  # pearson/poly/tetra) the smallest defect for a *used* constraint was 1.7e-3
  # away on the ML side, while the *unused* constraint came as low as 2.5e-4
  # (DWLS with polychoric correlations). Do not widen this to 1e-3 on the
  # assumption that the two are an order of magnitude apart -- they are not.
  # Ties are settled by taking the smaller defect, not by the tolerance alone.
  tol <- 1e-4
  A <- if (is.finite(d_g) && d_g <= tol && d_g <= d_w) {
    gram(M)
  } else if (is.finite(d_w) && d_w <= tol) {
    # A row at or past the Heywood boundary carries no usable uniqueness, so leave an
    # improper average in the gauge the alignment produced rather than letting a
    # degenerate row steer the weighting. `wgram()` floors psi at `.uniqueness_floor`,
    # which already bounds the weight (at 1 / 0.005 = 200, not the unbounded blow-up an
    # unfloored Psi^-1 would give), so this guard is the second line of defence and not
    # the thing that keeps the weight finite -- do not remove the floor on the strength
    # of it.
    if (any(1 - rowSums(M^2) <= 1e-6)) return(diag(k))
    wgram(M)
  } else {
    # The component fits sit in no gauge this recognises (an improper solution,
    # say). Imposing a foreign one would be worse than leaving the average alone.
    return(diag(k))
  }
  # If the average never left the gauge -- identical or near-identical
  # imputations -- there is nothing to correct, and rotating anyway would move
  # the result by the extraction's convergence tolerance for no gain.
  if (defect(A) <= tol) return(diag(k))
  e <- eigen((A + t(A)) / 2, symmetric = TRUE)
  # Two (near-)equal eigenvalues leave their eigenvectors undetermined: every
  # rotation inside that eigenspace diagonalises A equally well, so the gauge
  # would be settled by rounding noise and could come out differently under a
  # different BLAS. Leave such a solution as the alignment produced it.
  # The sensitivity of an eigenvector goes as (machine epsilon)/gap, so the
  # threshold has to be well above epsilon to bound it: at a relative gap of
  # 1e-6 the worst-case movement is around 1e-10, comfortably below anything a
  # user would see, whereas a threshold at 1e-8 would still admit bases that
  # move in the sixth decimal between one BLAS and another.
  if (k > 1L) {
    ev <- e$values
    if (min(abs(diff(ev))) <= 1e-6 * max(abs(ev))) return(diag(k))
  }
  C <- e$vectors
  # Orient each column to agree with the anchor rather than by an internal rule,
  # so the pooled solution keeps the sign convention its component fits use. The
  # anchor is the medoid, so this stays independent of the imputation order.
  MC <- M %*% C
  s <- sign(colSums(MC * anchor))
  s[s == 0] <- 1
  C %*% diag(s, nrow = k)
}

.efa_pooled_align_unrotated_list <- function(unrot_loadings,
                                             align_unrotated = c("signed_tucker_congruence", "none", "procrustes")) {
  # Unrotated factor axes are arbitrary up to ordering and signs. This helper
  # puts them into a common orientation before simple arithmetic averaging.
  # Unlike the rotated solution below, this step should not seek simple
  # structure; it only removes indeterminacy in the unrotated axes.
  #
  # Returns the aligned loadings together with a parallel list of per-imputation
  # alignment metadata, so that downstream consumers can apply the same column
  # permutation (and, when signed gauges matter, sign vector) to other
  # per-imputation matrices that share the loading gauge, e.g. marginal-SE
  # matrices used by Rubin's-rules pooling of analytic SEs.
  align_unrotated <- match.arg(align_unrotated)

  D <- length(unrot_loadings)
  if (align_unrotated == "none") {
    k_first <- if (D > 0L) ncol(as.matrix(unrot_loadings[[1]])) else 0L
    identity_meta <- list(type         = "none",
                          factor_order = seq_len(k_first),
                          factor_sign  = rep(1, k_first))
    meta <- rep(list(identity_meta), D)
    return(list(loadings = unrot_loadings, meta = meta))
  }

  # The signed-permutation default anchors on the medoid, so the pooled solution
  # does not depend on the order of the imputation list. `procrustes` keeps its
  # documented first-imputation anchor.
  anchor <- if (align_unrotated == "signed_tucker_congruence") {
    .efa_pooled_medoid_anchor(unrot_loadings)
  } else {
    1L
  }
  target <- as.matrix(unrot_loadings[[anchor]])
  k <- ncol(target)

  out  <- vector("list", D)
  meta <- vector("list", D)
  out[[anchor]] <- target
  meta[[anchor]] <- if (align_unrotated == "signed_tucker_congruence") {
    list(type         = "signed_tucker_congruence",
         factor_order = seq_len(k),
         factor_sign  = rep(1, k))
  } else {
    # Procrustes anchor gauge: Q_1 = I_k, so the aligned analytic-SE propagation
    # reduces to the identity on the first imputation.
    list(type = "procrustes", Q = diag(k))
  }

  for (d in setdiff(seq_along(unrot_loadings), anchor)) {
    Ld <- as.matrix(unrot_loadings[[d]])

    if (align_unrotated == "signed_tucker_congruence") {
      # Preserve the axes up to a one-to-one permutation/sign change based on
      # Tucker congruence. .align_solution() implements the assignment and sign
      # correction without applying a continuous rotation.
      al <- .align_solution(L_target = target, L = Ld)
      out[[d]]  <- al$loadings
      meta[[d]] <- list(type         = "signed_tucker_congruence",
                        factor_order = as.integer(al$factor_order),
                        factor_sign  = as.numeric(al$factor_sign))
    } else if (align_unrotated == "procrustes") {
      # Continuous orthogonal alignment: Q_d minimises ||L_d Q_d - L_1||_F over
      # orthogonal Q. The transform is retained alongside the aligned loadings
      # so the analytic-SE pool can propagate the per-imputation full unrotated
      # covariance via the Kronecker identity Var(vec(L_d Q_d)) = (Q_d' (x) I_p)
      # V_d (Q_d (x) I_p).
      pr <- efa_procrustes(A = Ld, Target = target, rotation = "orthogonal")
      out[[d]]  <- .change_class(pr$loadings, "matrix")
      meta[[d]] <- list(type = "procrustes", Q = pr$T)
    }
  }

  if (align_unrotated == "signed_tucker_congruence") {
    # Put the pooled solution back into the gauge every component fit uses. The
    # same orthogonal C is applied to each imputation, so the arithmetic mean of
    # the returned list is itself canonical and the per-imputation estimates stay
    # in the gauge of the pooled estimate they are pooled against.
    C <- .efa_pooled_canonical_gauge(.average_matrices(out), target)
    # Record the rotation only when it does something. Leaving `C` absent keeps
    # the cheaper marginal-SE alignment in `.efa_pooled_analytic_pool()`, which
    # needs no covariance block, for every case that does not need re-gauging.
    if (max(abs(C - diag(k))) > 1e-8) {
      dn <- dimnames(target)
      for (d in seq_len(D)) {
        out[[d]] <- structure(out[[d]] %*% C, dimnames = dn)
        meta[[d]]$C <- C
      }
    }
  }

  list(loadings = out, meta = meta)
}

.efa_pooled_rmsea_ci <- function(chi, df, N, level = .90) {
  # RMSEA confidence bounds for a pooled test statistic: the same noncentral
  # chi-square inversion .chi_fit_indices() uses for a single fit (.rmsea_lambda();
  # Browne & Cudeck, 1992), at an arbitrary confidence `level` and on the pooled
  # (N - 1) scale. The chi-square supplied here should already be the pooled
  # statistic. An undefined bound is reported as NA rather than aborting.
  if (is.na(chi) || is.na(df) || is.na(N) || df <= 0 || N <= 1) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  denom <- df * (N - 1)
  if (!is.finite(denom) || denom <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  alpha <- 1 - level
  # Cap the bounds at 1 like the non-pooled .chi_fit_indices() RMSEA CI, so all
  # pooled routes (information, bootstrap, MI2S) report RMSEA in [0, 1].
  c(lower = min(sqrt(.rmsea_lambda(chi, df, 1 - alpha / 2) / denom), 1),
    upper = min(sqrt(.rmsea_lambda(chi, df, alpha / 2) / denom), 1))
}

.efa_pooled_D2 <- function(chis, df) {
  # D2 pools complete-data chi-square statistics across imputations. The
  # returned chi value is the asymptotic chi-square approximation df * F_D2,
  # which is used downstream for RMSEA/CFI.
  chis <- chis[is.finite(chis)]
  M <- length(chis)

  if (M < 2L || !is.finite(df) || df <= 0) {
    return(NULL)
  }

  Tbar <- mean(chis)

  if (Tbar <= 0) {
    return(list(
      F = 0,
      chi = 0,
      df1 = df,
      df2 = Inf,
      p = 1,
      ARIV = 0,
      FMI = 0,
      Tbar = Tbar,
      M = M
    ))
  }

  # Average relative increase in variance (Li, Meng, Raghunathan & Rubin, 1991):
  # the between-imputation variance of the sqrt-transformed statistics is centred
  # on mean(sqrt(chi^2)) = stats::var(sqrt(chis)), not on sqrt(mean(chi^2)).
  ARIV <- (1 + 1 / M) * stats::var(sqrt(chis))

  if (ARIV <= .Machine$double.eps) {
    F_D2 <- Tbar / df
    df2 <- Inf
    p_val <- stats::pchisq(Tbar, df = df, lower.tail = FALSE)
  } else {
    F_D2 <- (Tbar / df - ((M + 1) / (M - 1)) * ARIV) / (1 + ARIV)
    df2 <- df^(-3 / M) * (M - 1) * (1 + 1 / ARIV)^2
    p_val <- if (F_D2 < 0) 1 else stats::pf(F_D2, df1 = df, df2 = df2, lower.tail = FALSE)
  }

  chi_D2 <- max(0, df * F_D2)
  FMI <- ARIV / (1 + ARIV)

  list(
    F = F_D2,
    chi = chi_D2,
    df1 = df,
    df2 = df2,
    p = p_val,
    ARIV = ARIV,
    FMI = FMI,
    Tbar = Tbar,
    M = M
  )
}

.efa_pooled_fit_indices <- function(fits,
                                    pooled_R,
                                    residuals,
                                    RMSR,
                                    N,
                                    Ns,
                                    pool_method = "D2",
                                    rmsea_ci_level = .90) {

  fit_list <- .extract_list_object(fits, "fit_indices")
  keep <- !vapply(fit_list, is.null, logical(1))
  fit_list <- fit_list[keep]

  # Per-imputation N aligned to the kept fits. Each imputation's chi-square was
  # computed in .gof() with that imputation's own N, so the rescaling to the
  # common (N - 1) scale used by the CFI/TLI diagnostics below must use the same
  # per-imputation N.
  Ns_kept <- as.numeric(Ns)[keep]

  p_vars <- nrow(pooled_R)
  df <- NA_real_
  chis <- numeric(0)
  chis_null <- numeric(0)

  if (length(fit_list) > 0L) {
    dfs <- vapply(fit_list, function(x) {
      if (!is.null(x$df)) x$df else NA_real_
    }, numeric(1))
    finite_dfs <- unique(stats::na.omit(dfs))
    if (length(finite_dfs) > 1L) {
      cli::cli_abort("Cannot D2-pool chi-square fit because the imputation-specific dfs differ.",
                     class = "efa_pooled_chisq_df")
    }
    if (length(finite_dfs) == 1L) {
      df <- finite_dfs[[1L]]
    }

    chis <- vapply(fit_list, function(x) {
      if (!is.null(x$chi)) x$chi else NA_real_
    }, numeric(1))
    chis_null <- vapply(fit_list, function(x) {
      if (!is.null(x$chi_null)) x$chi_null else NA_real_
    }, numeric(1))
  }

  # D2-pool a vector of per-imputation chi-squares at `dfree` degrees of freedom,
  # or NULL when pooling does not apply (non-D2 request, no statistics, or a
  # degenerate df). Shared by the reported model/baseline statistics and their
  # (N - 1)-scaled diagnostic counterparts (chi_cfi / chi_null_cfi).
  .pool_chi <- function(chi_stats, dfree) {
    if (identical(pool_method, "D2") && length(chi_stats) > 0L &&
        is.finite(dfree) && dfree > 0) {
      .efa_pooled_D2(chi_stats, dfree)
    } else {
      NULL
    }
  }

  D2 <- .pool_chi(chis, df)
  if (!is.null(D2)) {
    chi <- D2$chi
    p_chi <- D2$p
  } else {
    chi <- NA_real_
    p_chi <- NA_real_
  }

  df_null <- p_vars * (p_vars - 1) / 2
  D2_null <- .pool_chi(chis_null, df_null)

  if (!is.null(D2_null)) {
    chi_null <- D2_null$chi
    p_null <- D2_null$p
  } else if (is.finite(N) && N > 1) {
    chi_null <- .null_chisq(pooled_R, N)
    p_null <- stats::pchisq(chi_null, df_null, lower.tail = FALSE)
  } else {
    chi_null <- NA_real_
    p_null <- NA_real_
  }

  # The incremental indices CFI and TLI are reported as the average of the
  # per-imputation indices rather than formed from separately pooled model/baseline
  # chi-squares; see the rationale in the function documentation (@details). Average
  # via .efa_pooled_col_means(), the same finite-mean primitive that centres the
  # bootstrap fit_indices_descriptive CI, so the reported point estimate and that CI
  # share one definition.
  cfis <- vapply(fit_list, function(x) {
    if (!is.null(x$CFI)) x$CFI else NA_real_
  }, numeric(1))
  tlis <- vapply(fit_list, function(x) {
    if (!is.null(x$TLI)) x$TLI else NA_real_
  }, numeric(1))
  # Columns are positional (col_means drops dimnames): column 1 is CFI, column 2 TLI.
  inc_means <- .efa_pooled_col_means(cbind(cfis, tlis))
  CFI <- inc_means[[1L]]
  TLI <- inc_means[[2L]]

  # The model and baseline chi-squares D2-pooled on the common (N - 1) noncentrality
  # scale, the basis on which lavaan.mi/semTools form pooled incremental indices.
  # chi_cfi is the (N - 1)-scaled pooled model statistic that drives the reported RMSEA
  # below; chi_null_cfi is exposed via mi_diagnostics so the pooled fit can be reconciled
  # against that reference (the reported CFI/TLI are averaged per imputation, not formed
  # from these). Each imputation's Bartlett-corrected chi is mapped back to the
  # (N - 1) scale -- chi = F * mult implies F * (N - 1) = chi * (N - 1) / mult, with
  # .bartlett_mult() shared with .gof()/.null_chisq() so the multiplier cannot drift
  # -- before pooling, since D2 is non-linear in a constant rescaling. The model and
  # baseline use different multipliers (the factor-count term (2q)/3 enters only the
  # model), so a small-N imputation can pass one finiteness mask (ok_m / ok_n) and
  # not the other; the model mask ok_m gates the RMSEA statistic, ok_n the chi_null_cfi
  # diagnostic. q is read from the first fit, valid because .efa_pooled_check_fits()
  # pins an equal factor count across imputations.
  chi_cfi <- NA_real_
  chi_null_cfi <- NA_real_
  if (identical(pool_method, "D2") && length(chis) > 0L && is.finite(df) && df > 0) {
    q <- ncol(as.matrix(fits[[1L]]$unrot_loadings))
    mult_model <- .bartlett_mult(Ns_kept, p_vars, q)
    mult_null <- .bartlett_mult(Ns_kept, p_vars)
    chis_cfi <- rep(NA_real_, length(chis))
    ok_m <- is.finite(Ns_kept) & is.finite(mult_model) & mult_model > 0
    chis_cfi[ok_m] <- chis[ok_m] * (Ns_kept[ok_m] - 1) / mult_model[ok_m]
    chis_null_cfi <- rep(NA_real_, length(chis_null))
    ok_n <- is.finite(Ns_kept) & is.finite(mult_null) & mult_null > 0
    chis_null_cfi[ok_n] <- chis_null[ok_n] * (Ns_kept[ok_n] - 1) / mult_null[ok_n]

    D2_cfi <- .pool_chi(chis_cfi, df)
    if (!is.null(D2_cfi)) chi_cfi <- D2_cfi$chi
    D2_null_cfi <- .pool_chi(chis_null_cfi, df_null)
    if (!is.null(D2_null_cfi)) chi_null_cfi <- D2_null_cfi$chi
  }

  # AIC, BIC and ECVI are chi-square-derived information criteria with no standard
  # interpretation once the component statistic is a scaled one (the Satorra-Bentler
  # corrected two-stage statistic of Yuan, Marshall & Bentler, 2002, and its
  # scaled-shifted variants), so .gof() withholds them from every fit that reports
  # such a statistic. A D2 pool cannot rehabilitate them, so take the components'
  # decision rather than recomputing from the pooled chi-square: a component that
  # produced a chi-square but no AIC withheld it deliberately. Reading it off the
  # returned indices covers the `chi_scaled_type` path (.apply_scaled_test() NAs all
  # three) and the tagless FIML fallbacks alike -- the plain-LRT and just-identified
  # cases at R/fit-indices.R:507-515 never reach .apply_scaled_test() and so carry no
  # tag, but they withhold the criteria just the same. A fit with no chi-square at all
  # (PAF, DWLS) is not a withholding: its pooled chi is NA, which zeroes these anyway.
  withheld_components <- any(vapply(fit_list, function(x) {
    !is.null(x$chi) && is.finite(x$chi) && (is.null(x$AIC) || !is.finite(x$AIC))
  }, logical(1)))

  if (!withheld_components && is.finite(chi) && is.finite(df) && is.finite(N) && N > 1) {
    n_params <- p_vars * (p_vars + 1) / 2 - df
    ECVI <- (chi + 2 * n_params) / (N - 1)
  } else {
    ECVI <- NA_real_
  }

  # RMSEA is built on the uncorrected (N - 1) discrepancy scale (as in
  # .chi_fit_indices()), so it uses the (N - 1)-scaled pooled model statistic chi_cfi,
  # not the Bartlett-corrected chi (which drives the reported model chi-square test).
  rmsea_denom <- df * (N - 1)
  if (is.finite(chi_cfi) && is.finite(df) && is.finite(N) &&
      df > 0 && N > 1 && is.finite(rmsea_denom) && rmsea_denom > 0) {
    RMSEA <- .rmsea_point(chi_cfi, df, N)
    rmsea_ci <- .efa_pooled_rmsea_ci(chi_cfi, df, N, level = rmsea_ci_level)
  } else {
    RMSEA <- NA_real_
    rmsea_ci <- c(lower = NA_real_, upper = NA_real_)
  }

  # These mirror the chi-square-derived quantities used elsewhere in EFAtools.
  # Under MI/D2 pooling they are descriptive only, not likelihood-based MI AIC/BIC,
  # and they are withheld entirely when the component fits withheld them (see
  # `withheld_components` above).
  AIC <- if (!withheld_components && is.finite(chi) && is.finite(df)) {
    chi - 2 * df
  } else {
    NA_real_
  }
  BIC <- if (!withheld_components && is.finite(chi) && is.finite(df) && is.finite(N)) {
    chi - log(N) * df
  } else {
    NA_real_
  }

  ## CAF in EFAtools is 1 - KMO(delta_hat), with diagonal set to 1.
  # Ensure symmetry of residuals by (residuals + t(residuals)) / 2. Asymmetry
  # can arise due to floating point imprecision from averaging the matrices.
  delta_hat <- (residuals + t(residuals)) / 2
  diag(delta_hat) <- 1
  CAF <- .compute_caf(delta_hat)

  ## SRMR (standardized RMR; Bentler, 1995) over the pooled residuals
  SRMR <- .srmr(residuals)

  # The fit indices and the pooling diagnostics are returned separately: the pooled
  # object carries the diagnostics in its own top-level `mi_diagnostics` slot, so the
  # `fit_indices` slot keeps the same all-scalar shape a single efa_fit() returns.
  list(
    fit_indices = list(
      chi = chi,
      df = df,
      p_chi = p_chi,
      CAF = CAF,
      CFI = CFI,
      TLI = TLI,
      RMSEA = RMSEA,
      RMSEA_LB = rmsea_ci[["lower"]],
      RMSEA_UB = rmsea_ci[["upper"]],
      AIC = AIC,
      BIC = BIC,
      ECVI = ECVI,
      RMSR = RMSR,
      SRMR = SRMR,
      chi_null = chi_null,
      df_null = df_null,
      p_null = p_null,
      pool_method = pool_method
    ),
    mi_diagnostics = list(
      D2_F = if (!is.null(D2)) D2$F else NA_real_,
      D2_df1 = if (!is.null(D2)) D2$df1 else NA_real_,
      D2_df2 = if (!is.null(D2)) D2$df2 else NA_real_,
      D2_chi_asymptotic = if (!is.null(D2)) D2$chi else NA_real_,
      chi_bar_naive = if (!is.null(D2)) D2$Tbar else if (length(chis) > 0L) mean(chis, na.rm = TRUE) else NA_real_,
      D2_null_F = if (!is.null(D2_null)) D2_null$F else NA_real_,
      D2_null_df1 = if (!is.null(D2_null)) D2_null$df1 else NA_real_,
      D2_null_df2 = if (!is.null(D2_null)) D2_null$df2 else NA_real_,
      D2_null_chi_asymptotic = if (!is.null(D2_null)) D2_null$chi else NA_real_,
      chi_null_bar_naive = if (!is.null(D2_null)) D2_null$Tbar else if (length(chis_null) > 0L) mean(chis_null, na.rm = TRUE) else NA_real_,
      ARIV = if (!is.null(D2)) D2$ARIV else NA_real_,
      FMI = if (!is.null(D2)) D2$FMI else NA_real_,
      ARIV_null = if (!is.null(D2_null)) D2_null$ARIV else NA_real_,
      FMI_null = if (!is.null(D2_null)) D2_null$FMI else NA_real_,
      # Pooled model and baseline chi-squares on the common (N - 1) noncentrality
      # scale -- the lavaan.mi/semTools basis for pooled incremental indices, kept
      # for reconciliation against that reference. The reported CFI/TLI average the
      # per-imputation indices instead; these are diagnostic and distinct from the
      # reported chi / chi_null, which keep their Bartlett multipliers.
      chi_cfi = chi_cfi,
      chi_null_cfi = chi_null_cfi,
      m = if (!is.null(D2)) D2$M else length(fits)
    )
  )
}

## -----------------------------------------------------------------------------
## Rubin pooling helpers, shared by the analytic and bootstrap routes
## -----------------------------------------------------------------------------

.efa_pooled_has_replicates <- function(fits) {
  # Bootstrap MI pooling requires the raw unrotated bootstrap loading replicates.
  # Scalar bootstrap SEs/CIs are not enough because the replicates must be
  # re-expressed in the final MI target space.
  vapply(fits, function(x) {
    !is.null(x$replicates) &&
      !is.null(x$replicates$unrot_loadings)
  }, logical(1))
}

.efa_pooled_vec <- function(M) {
  # Column-major vectorization so matrix estimates can be Rubin-pooled.
  as.vector(as.matrix(M))
}

.efa_pooled_communalities <- function(Lambda, Phi = NULL) {
  # Communalities are the diagonal of the common-factor reproduced matrix; read it
  # directly instead of forming the full p x p product:
  # diag(L L') = rowSums(L^2) and diag(L Phi L') = rowSums((L Phi) * L).
  Lambda <- as.matrix(Lambda)
  if (is.null(Phi)) {
    unname(rowSums(Lambda^2))
  } else {
    unname(rowSums((Lambda %*% as.matrix(Phi)) * Lambda))
  }
}

.efa_pooled_model_implied <- function(Lambda, Phi = NULL) {
  # Correlation-model implied matrix: common-factor part plus uniquenesses,
  # implemented by setting the diagonal of Lambda Phi Lambda' to one.
  Lambda <- as.matrix(Lambda)
  common_R <- if (is.null(Phi)) {
    Lambda %*% t(Lambda)
  } else {
    Lambda %*% as.matrix(Phi) %*% t(Lambda)
  }
  implied <- common_R
  diag(implied) <- 1
  implied
}

.efa_pooled_residual_from_solution <- function(R, Lambda, Phi = NULL) {
  # Residuals for one imputation in the same target space as its aligned solution.
  E <- as.matrix(R) - .efa_pooled_model_implied(Lambda, Phi)
  diag(E) <- 0
  E
}

.efa_pooled_make_ci <- function(est, se, df, alpha) {
  # Wald-type confidence intervals using Rubin degrees of freedom when finite.
  # `est` may be a plug-in estimate (communalities, Structure, symmetrised Phi)
  # supplied via the callers' `est_override`, in which case the interval is
  # centred on the plug-in while its half-width (se, df) comes from the pooled
  # column-mean's between/within variance; the two coincide when the plug-in
  # equals that column mean and differ only by the (typically small) plug-in gap.
  crit <- ifelse(
    is.finite(df),
    stats::qt(1 - alpha / 2, df = df),
    stats::qnorm(1 - alpha / 2)
  )
  list(
    lower = est - crit * se,
    upper = est + crit * se
  )
}


.efa_pooled_col_means <- function(x) {
  # Column means that return NA for columns without finite values.
  x <- as.matrix(x)
  if (ncol(x) < 1L) {
    return(numeric(0))
  }
  vapply(seq_len(ncol(x)), function(j) {
    vals <- x[, j]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 1L) {
      NA_real_
    } else {
      mean(vals)
    }
  }, numeric(1))
}

.efa_pooled_col_vars <- function(x) {
  # Column variances that return NA for columns with fewer than two finite values.
  # Rubin pooling only needs diagonal within- and between-imputation variances,
  # so computing full covariance matrices is unnecessary and can be prohibitively
  # memory-intensive for residual matrices.
  x <- as.matrix(x)
  if (ncol(x) < 1L) {
    return(numeric(0))
  }
  if (nrow(x) < 2L) {
    return(rep(NA_real_, ncol(x)))
  }
  vapply(seq_len(ncol(x)), function(j) {
    vals <- x[, j]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2L) {
      NA_real_
    } else {
      stats::var(vals)
    }
  }, numeric(1))
}

.efa_pooled_rubin_core <- function(Qbar, Ubar, B, m, N_pool = Inf,
                                   alpha = 0.05) {
  # Element-wise Rubin-style pooling core, shared by the bootstrap (analytic SEs
  # absent) and analytic (per-imputation marginal SEs available) paths.
  #
  # Qbar  pooled point estimate per parameter (vector or matrix-as-vector).
  # Ubar  mean within-imputation variance (vector).
  # B     between-imputation variance (sample, m-1 denominator) (vector).
  # m     number of imputations.
  # N_pool complete-data sample size for Barnard-Rubin small-sample df. When
  #        non-finite, the plain Rubin (1987) df ν_old = (m-1)(1 + 1/r)^2 is
  #        used (lavaan.mi convention). When finite, the Barnard & Rubin
  #        (1999, Biometrika 86) adjustment is applied:
  #           ν_obs = ((ν_com + 1)/(ν_com + 3)) ν_com (1 - λ),  ν_com = N_pool - 1
  #           ν_BR  = 1 / (1/ν_old + 1/ν_obs)
  #        Every pooling route in the package reports the asymptotically normal
  #        reference and so passes N_pool = Inf; the finite-N_pool branch is a
  #        tested extension point, kept because it is the only place the
  #        small-sample reference is implemented. It reproduces
  #        mice::pool.scalar(n = N_pool) up to the documented lambda-vs-gamma FMI
  #        convention below, and is exercised by the test suite.
  # alpha 1 - confidence level for the Wald interval.
  #
  # Returns a list(se, ci, RIV, FMI, df). Bootstrap callers receive the same
  # df as before (N_pool defaults to Inf for byte-identical behavior).
  Tot <- Ubar + (1 + 1 / m) * B
  Tot[is.finite(Tot) & Tot < 0 & abs(Tot) < sqrt(.Machine$double.eps)] <- 0
  se <- sqrt(Tot)

  r <- rep(0, length(se))
  has_U <- is.finite(Ubar) & Ubar > 0
  r[has_U] <- ((1 + 1 / m) * B[has_U]) / Ubar[has_U]
  r[!has_U & is.finite(B) & B > 0] <- Inf
  r[!is.finite(B) & !is.infinite(r)] <- NA_real_

  df_old <- rep(Inf, length(se))
  finite_pos_r <- is.finite(r) & r > 0
  df_old[finite_pos_r] <- (m - 1) * (1 + 1 / r[finite_pos_r])^2
  df_old[is.infinite(r)] <- m - 1
  df_old[is.na(r)] <- NA_real_

  # Reported FMI is Rubin's asymptotic lambda = r / (1 + r) (the limiting fraction
  # of missing information; equals lavaan.mi's `fmi`), not mice's finite-m gamma
  # = (r + 2/(df_old + 3)) / (r + 1), which adds a small-m correction. The two
  # agree as m grows; lambda is the conventional EFA/SEM-MI reporting choice.
  fmi <- r / (1 + r)
  fmi[is.infinite(r)] <- 1
  fmi[is.nan(fmi)] <- 0

  if (is.finite(N_pool) && N_pool > 1) {
    nu_com <- N_pool - 1
    nu_obs <- ((nu_com + 1) / (nu_com + 3)) * nu_com * (1 - fmi)
    # ν_BR = 1/(1/ν_old + 1/ν_obs). When ν_old is Inf (no between-variance),
    # 1/ν_old = 0 and ν_BR collapses to ν_obs (the complete-data t reference).
    df <- 1 / (1 / df_old + 1 / nu_obs)
    df[is.infinite(df_old)] <- nu_obs[is.infinite(df_old)]
    # When fmi -> 1 (Ubar = 0, B > 0), ν_obs -> 0 and the BR formula collapses
    # to 0, sending qt(., 0) to NaN. Fall back to the plain-Rubin ν_old in that
    # regime so the Wald interval is well-defined (wide, but finite).
    no_obs <- is.finite(nu_obs) & nu_obs == 0
    df[no_obs] <- df_old[no_obs]
    df[is.na(df_old)] <- NA_real_
  } else {
    df <- df_old
  }

  ci <- .efa_pooled_make_ci(Qbar, se, df, alpha)

  list(
    se  = se,
    Tot = Tot,
    ci  = ci,
    RIV = r,
    FMI = fmi,
    df  = df
  )
}

.efa_pooled_rubin_pool <- function(q_list, boot_mat_list, alpha = 0.05,
                                   est_override = NULL) {
  # Rubin-style pooling for one parameter family.
  #
  # q_list contains one aligned point-estimate vector per imputation. boot_mat_list
  # contains aligned bootstrap replicates for estimating the within-imputation
  # variances U_d. Only diagonal variances are required for the returned SEs, CIs,
  # RIV, FMI, and df, so the implementation intentionally avoids building full
  # p x p covariance matrices for each parameter vector.
  m <- length(q_list)
  q_mat <- do.call(rbind, q_list)
  q_bar <- .efa_pooled_col_means(q_mat)

  U_diag_mat <- do.call(rbind, lapply(boot_mat_list, .efa_pooled_col_vars))
  U_bar <- colMeans(U_diag_mat)
  B_mi <- .efa_pooled_col_vars(q_mat)

  est <- if (is.null(est_override)) q_bar else est_override
  core <- .efa_pooled_rubin_core(est, U_bar, B_mi, m, N_pool = Inf, alpha = alpha)

  list(
    est = est,
    se  = core$se,
    ci  = core$ci,
    RIV = core$RIV,
    FMI = core$FMI,
    df  = core$df
  )
}

.efa_pooled_propagate_procrustes_vcov <- function(V, Q, p, k) {
  # Marginal SEs of L_d Q_d implied by the per-imputation full unrotated
  # covariance V_d (p k x p k over column-major vec(L_d)) and a fixed orthogonal
  # Q_d (k x k). Under the column-major identity vec(L_d Q_d) = (Q_d' (x) I_p)
  # vec(L_d), V_d^aligned = (Q_d' (x) I_p) V_d (Q_d (x) I_p). For element (i, j)
  # only the variance of L_d Q_d at row i is required, so the implementation
  # extracts the k x k row sub-block V_d^(i) at indices (0:(k-1)) p + i and
  # returns sqrt(diag(Q_d' V_d^(i) Q_d)) per row. Cost is O(p k^3).
  #
  # Q_d is treated as fixed (Schoenemann 1966; mirrors the lavaan / mice
  # convention of not differentiating the per-imputation orthogonal Procrustes
  # solution back through the SVD on vec(L_d)).
  #
  # NA propagation is fail-closed at the row level: if any element of the row
  # block of V or of Q is non-finite, the entire row of the returned p x k
  # marginal-SE matrix is NA so the per-element NA mask in the analytic pool
  # can blank out the affected pooled outputs.
  q_finite <- all(is.finite(Q))
  out <- matrix(NA_real_, p, k)
  if (!q_finite) return(out)
  for (i in seq_len(p)) {
    idx <- (seq_len(k) - 1L) * p + i
    block <- V[idx, idx, drop = FALSE]
    if (!all(is.finite(block))) next
    # diag(Q' (block Q)) = colSums(Q * (block Q)); avoids forming the full k x k
    # product only to keep its diagonal.
    BQ <- block %*% Q
    out[i, ] <- sqrt(pmax(colSums(Q * BQ), 0))
  }
  out
}

.efa_pooled_analytic_marginal <- function(q_mat, se_mat, m, alpha,
                                          est_override = NULL) {
  # Rubin pool of one parameter family from aligned per-imputation point estimates
  # (rows of `q_mat`) and aligned per-imputation marginal SEs (rows of `se_mat`, so
  # the within-imputation variance U_d is `se_mat^2`). Only diagonal variances are
  # needed. Fail-closed: any per-imputation NA in q or se for an element blanks the
  # pooled outputs (SE, CI, RIV, FMI, df) at that element. `est_override` supplies
  # the pooled point estimate for plug-in quantities (communalities, structure)
  # whose returned estimate is computed from the pooled solution, not the mean of
  # the per-imputation values.
  Qbar <- .efa_pooled_col_means(q_mat)
  Ubar <- .efa_pooled_col_means(se_mat^2)
  B    <- .efa_pooled_col_vars(q_mat)

  na <- apply(is.na(q_mat) | is.na(se_mat), 2, any)
  Qbar[na] <- NA_real_; Ubar[na] <- NA_real_; B[na] <- NA_real_

  est <- if (is.null(est_override)) Qbar else est_override
  # Plain Rubin (1987) df (N_pool = Inf): analytic information-matrix loadings are
  # asymptotically normal (Wald-z), so the complete-data reference is normal
  # (nu_com = Inf) and the Barnard-Rubin small-sample adjustment does not apply -- it
  # reduces to plain Rubin. This matches lavaan.mi, which reports df -> Inf for
  # asymptotically-normal SEM parameters.
  .efa_pooled_rubin_core(est, Ubar, B, m, N_pool = Inf, alpha = alpha)
}

.efa_pooled_assemble_family <- function(pool, reshape, method = NULL) {
  # Reshape one Rubin-pooled family's SE, Wald CI bounds, and MI diagnostics with a
  # single reshape function (loading matrix, uniqueness/communality vector, or
  # vech-to-symmetric). `method` tags the gauge-alignment provenance when supplied.
  mi <- list(RIV = reshape(pool$RIV), FMI = reshape(pool$FMI), df = reshape(pool$df))
  if (!is.null(method)) mi$method <- method
  list(
    SE = reshape(pool$se),
    CI = list(lower = reshape(pool$ci$lower), upper = reshape(pool$ci$upper)),
    MI = mi
  )
}

.efa_pooled_rubin_result <- function(pool, reshape, dimnames = NULL) {
  # Convert vectorized Rubin output (estimate, SE, and Wald CI bounds) back into one
  # parameter family's own shape with a single `reshape` function -- a rectangular
  # matrix or a vech-to-symmetric expansion -- the way .efa_pooled_assemble_family()
  # does for the analytic route's SE/CI/MI triple.
  shape <- function(v) {
    out <- reshape(v)
    if (!is.null(dimnames)) dimnames(out) <- dimnames
    out
  }

  list(est = shape(pool$est), se = shape(pool$se),
       ci = list(lower = shape(pool$ci$lower), upper = shape(pool$ci$upper)))
}
