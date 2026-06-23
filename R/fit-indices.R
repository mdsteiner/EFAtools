# Goodness-of-fit indices for a fitted factor solution: explained-variance accounting,
# RMSR, the KMO-based CAF, and the chi-square / CFI / RMSEA / AIC / BIC block from `.gof()`.

#' Compute explained variances from loadings
#'
#' From unrotated loadings compute the communalities and uniquenesses for total
#' variance. Compute explained variances per factor from rotated loadings (and
#' factor intercorrelations Phi if oblique rotation was used).
#'
#' @param L_unrot matrix. Unrotated factor loadings.
#' @param L_rot matrix. Rotated factor loadings.
#' @param Phi matrix. Factor intercorrelations. Provide only if oblique rotation
#'  is used.
#'
#' @return A matrix with sum of squared loadings, proportion explained variance
#'  from total variance per factor, same as previous but cumulative, Proportion
#'  of explained variance from total explained variance, and same as previous but
#'  cumulative.
.compute_vars <- function(L_unrot, L_rot, Phi = NULL) {

  if (is.null(Phi)) {
    # sum of squared loadings per factor; colSums() already returns the scalar for
    # a single-factor solution, so no special case is needed.
    vars <- colSums(L_rot^2)
  } else {
    # compute variance proportions
    vars <- diag(Phi %*% t(L_rot) %*% L_rot)
  }

  # Compute the explained variances. The code is based on the psych::fac() function
  # total variance (sum of communalities and uniquenesses)
  h2 <- rowSums(L_unrot^2)            # diag(L L'), without forming the full p x p product
  var_total <- sum(h2 + (1 - h2))
  vars_explained <- rbind(`SS loadings` = vars)
  vars_explained <- rbind(vars_explained, `Prop Tot Var` = vars / var_total)

  if (ncol(L_rot) > 1) {
    vars_explained <- rbind(vars_explained,
                            `Cum Prop Tot Var` = cumsum(vars / var_total))
    vars_explained <- rbind(vars_explained,
                            `Prop Comm Var` = vars / sum(vars))
    vars_explained <- rbind(vars_explained,
                            `Cum Prop Comm Var` = cumsum(vars / sum(vars)))
  }

  vars_explained
}

.rmsr <- function(residuals, upper = TRUE) {
  E <- as.matrix(residuals)
  if (isTRUE(upper)) {
    vals <- E[upper.tri(E, diag = FALSE)]
  } else {
    vals <- E[row(E) != col(E)]
  }
  sqrt(mean(vals^2, na.rm = TRUE))
}

.compute_caf <- function(delta_hat) {

  delta_hat_KMO <- try(.compute_kmo(delta_hat)$KMO, silent = TRUE)

  if (inherits(delta_hat_KMO, "try-error") || is.na(delta_hat_KMO)) {
    CAF <- 0
    cli::cli_warn(
      c("CAF could not be computed; it was set to {.val 0} (the worst value).",
        "i" = "Inspect the results carefully."),
      class = "efa_caf_failed"
    )
  } else {
    CAF <- 1 - delta_hat_KMO
  }

  CAF

}

# Bartlett's (1951) small-sample multiplier for the ML/ULS discrepancy chi-square,
# N - 1 - (2p + 5)/6 - (2q)/3 with p variables and q factors. q = 0 gives the
# independence-model (sphericity) multiplier. Shared by the model chi-square in
# .gof(), the baseline chi-square in .null_chisq(), and EFA_POOLED()'s pooled
# (N - 1) rescaling so all three use one definition.
.bartlett_mult <- function(N, p, q = 0) N - 1 - (2 * p + 5) / 6 - (2 * q) / 3

# Independence-model (baseline) chi-square: Bartlett's test of sphericity, the Bartlett-
# corrected ML discrepancy of the null model (model-implied matrix = identity),
# -log|R| * (N - 1 - (2p + 5)/6). Reported as the baseline statistic and reused by HULL.
# (The incremental indices CFI/TLI in .gof() place the model and baseline on a common
# (N - 1) scale separately; see .chi_fit_indices().)
.null_chisq <- function(R, N, ld = determinant(R, logarithm = TRUE)) {
  p <- ncol(R)
  # Use the log-determinant directly: det(R) underflows to 0 for large,
  # near-singular (but still positive-definite) matrices, which would otherwise
  # turn the statistic into NA and propagate into .gof(), SMT(), and BARTLETT().
  # `ld` defaults to that decomposition; callers that already have it (e.g. .gof()
  # for the model chi-square) pass it in to avoid recomputing determinant(R).
  if (ld$sign <= 0 || !is.finite(ld$modulus)) return(NA_real_)
  # Guard the Bartlett multiplier: for small N relative to p it turns
  # non-positive, which would flip the baseline chi-square sign and propagate a
  # spurious statistic into .gof(), SMT(), and BARTLETT().
  mult <- .bartlett_mult(N, p)
  if (is.na(mult) || mult <= 0) return(NA_real_)
  -as.numeric(ld$modulus) * mult
}

# Noncentrality parameter for an RMSEA confidence bound: solve pchisq(chi, df, ncp) = goal
# for ncp via stats::uniroot (the 90% CI lower bound uses goal = .95, the upper bound goal =
# .05; Browne & Cudeck, 1992). Returns 0 when the model chi-square already lies below the
# target quantile, so the bound collapses to 0. Shared by .gof() and SMT().
.rmsea_lambda <- function(chi, df, goal) {
  # An undefined chi-square (e.g. the null model for a tiny N, where .null_chisq()
  # returns NA) has no noncentrality bound; propagate NA rather than failing the
  # `if (pchisq(NA) >= goal)` test with "missing value where TRUE/FALSE needed".
  if (is.na(chi)) return(NA_real_)
  if (stats::pchisq(chi, df = df, ncp = 0) >= goal) {
    p_chi_fun <- function(x, val, df, goal) goal - stats::pchisq(val, df, ncp = x)
    stats::uniroot(f = p_chi_fun, interval = c(1e-10, 10000), val = chi, df = df,
                   goal = goal, extendInt = "upX", maxiter = 100L)$root
  } else {
    0
  }
}

# CFI (Bentler, 1990) and TLI (Tucker & Lewis, 1973) from model and baseline
# noncentralities on a common scaling constant: chi_cfi and chi_null_cfi, each with
# its df subtracted and floored at 0 so an over-fitting model cannot deflate the
# index and it stays in [0, 1] (matches lavaan::fitMeasures()). Returns NA for both
# when either statistic is undefined. Called from .chi_fit_indices() (the per-model
# path); EFA_POOLED() reports the average of the per-imputation indices this produces
# rather than calling it on pooled statistics.
.incremental_fit <- function(chi_cfi, df, chi_null_cfi, df_null) {
  if (is.na(chi_cfi) || is.na(chi_null_cfi)) {
    return(list(CFI = NA_real_, TLI = NA_real_))
  }
  d_m <- max(chi_cfi - df, 0)
  d_null <- max(chi_null_cfi - df_null, 0)
  CFI <- if (df == 0 || max(d_m, d_null) == 0) 1 else 1 - d_m / max(d_m, d_null)
  if (df == 0) {
    TLI <- 1
  } else {
    ratio_null <- chi_null_cfi / df_null
    # ratio_null == 1 leaves the index undefined (0/0); report a perfect 1.
    TLI <- if (ratio_null == 1) 1 else (ratio_null - chi_cfi / df) / (ratio_null - 1)
  }
  list(CFI = CFI, TLI = TLI)
}

# Chi-square-derived fit indices: the p-values, CFI (Bentler, 1990), TLI/NNFI (Tucker &
# Lewis, 1973), RMSEA with its analytic 90% bounds (Browne & Cudeck, 1992), AIC, BIC, and
# ECVI (Browne & Cudeck, 1989). Computed from a model chi-square and an independence-
# baseline chi-square with their degrees of freedom; `m` is the number of variables, `ci`
# toggles the (uniroot-solved) RMSEA bounds. Returns NA throughout when the model chi-
# square is undefined. Shared by .gof() (the ML/ULS discrepancy chi-square) and the
# scaled-chi-square sandwich-SE path, which passes its scaled model and baseline statistics
# in, so both report the block by exactly the same formulas.
.chi_fit_indices <- function(chi, df, chi_null, df_null, N, m, ci,
                             chi_cfi = chi, chi_null_cfi = chi_null) {

  if (is.na(chi)) {
    return(list(p_chi = NA, CFI = NA, TLI = NA, RMSEA = NA, RMSEA_LB = NA,
                RMSEA_UB = NA, AIC = NA, BIC = NA, ECVI = NA, p_null = NA))
  }

  p_null <- stats::pchisq(chi_null, df_null, lower.tail = F)
  # A just-identified model (df = 0) has no chi-square test: pchisq(chi, 0) returns 0 for
  # any chi > 0, i.e. a spurious "p < .001". Report the p-value as undefined instead,
  # matching the df == 0 handling of CFI/TLI/RMSEA below.
  p_chi <- if (df == 0) NA_real_ else stats::pchisq(chi, df, lower.tail = F)

  # CFI and TLI compare the model and baseline noncentralities, which their definitions
  # (Bentler, 1990; Tucker & Lewis, 1973) require on a single scaling constant. `chi_cfi`
  # and `chi_null_cfi` carry the model and baseline statistics on that common scale (.gof()
  # passes the ML/ULS discrepancy on the (N - 1) noncentrality scale; the scaled sandwich
  # path leaves them at their defaults, since its scaled model and baseline statistics are
  # already comparable). When the baseline is undefined (e.g. a degenerate scaled baseline)
  # CFI/TLI are NA while the model-only RMSEA/AIC/BIC/p remain computable.
  inc <- .incremental_fit(chi_cfi, df, chi_null_cfi, df_null)
  CFI <- inc$CFI
  TLI <- inc$TLI

  if (df != 0) {

    # RMSEA point estimate and bounds, formula 12.6 from Kline (2016, Principles and
    # Practice of Structural Equation Modeling). The noncentrality chi - df is divided by
    # (N - 1) while the model chi-square it draws on carries the Bartlett multiplier; this
    # conventional (N - 1) scaling differs from a strictly self-consistent Browne & Cudeck
    # (1992) interval by well under 1%.
    RMSEA <- sqrt(max(0, chi - df) / (df * (N - 1)))
    if (RMSEA > 1) RMSEA <- 1

    if (isTRUE(ci)) {

      # Analytic 90% RMSEA confidence bounds via the noncentrality solver (Browne &
      # Cudeck, 1992); skipped when ci = FALSE (e.g. per-replicate bootstrap fits).
      lambda_l <- .rmsea_lambda(chi, df, .95)
      lambda_u <- .rmsea_lambda(chi, df, .05)

      RMSEA_LB <- sqrt(lambda_l / (df * (N - 1)))
      RMSEA_UB <- sqrt(lambda_u / (df * (N - 1)))

      if (RMSEA_LB > 1) RMSEA_LB <- 1
      if (RMSEA_UB > 1) RMSEA_UB <- 1

    } else {

      RMSEA_LB <- NA_real_
      RMSEA_UB <- NA_real_

    }

  } else {

    RMSEA <- 0
    RMSEA_LB <- 0
    RMSEA_UB <- 0

  }

  AIC <- chi - 2 * df
  BIC <- chi - log(N) * df
  n_params <- m * (m + 1) / 2 - df
  ECVI <- (chi + 2 * n_params) / (N - 1)

  list(p_chi = p_chi, CFI = CFI, TLI = TLI, RMSEA = RMSEA, RMSEA_LB = RMSEA_LB,
       RMSEA_UB = RMSEA_UB, AIC = AIC, BIC = BIC, ECVI = ECVI, p_null = p_null)
}

# Degrees of freedom of an EFA solution with `m` variables and `q` factors, ((m - q)^2 -
# (m + q)) / 2 (equivalently factanal's 0.5 * ((m - q)^2 - m - q)). Shared so the model
# fit and the identification check use one formula.
.efa_df <- function(m, q) ((m - q)^2 - (m + q)) / 2

.gof <- function(L, # The loading/ pattern matrix
                 R, # The correlation matrix
                 N, # The number of cases
                 method, # The estimation method
                 Fm, # Minimized error
                 ci = TRUE) { # Compute the analytic RMSEA confidence bounds
  m <- nrow(L)
  q <- ncol(L)

  # dfs
  df <- .efa_df(m, q)

  # Model-implied (no-uniqueness) correlations LL'; reused for the residual indices
  # below and the chi-square's Sigma, so form the p x p product once.
  LLt <- L %*% t(L)

  ### compute CAF
  delta_hat <- R - LLt
  diag(delta_hat) <- 1
  CAF <- .compute_caf(delta_hat)

  ### compute RMSR
  RMSR <- .rmsr(delta_hat)

  ### compute SRMR (standardized root mean square residual; Bentler, 1995). The
  ### model-implied diagonal is 1, so only the off-diagonal residuals contribute; the
  ### denominator is the count of non-redundant elements p(p + 1)/2.
  SRMR <- sqrt(sum(delta_hat[upper.tri(delta_hat)] ^ 2) / (m * (m + 1) / 2))

  # Model chi-square: the Bartlett-corrected (Bartlett, 1951) ML discrepancy
  # F = tr(Sigma^-1 R) - log|Sigma^-1 R| - p, evaluated at the model-implied correlation
  # matrix Sigma = LL' (unit diagonal), times (N - 1 - (2p + 5)/6 - (2q)/3). For ML this
  # equals the ML objective times the Bartlett multiplier (matching stats::factanal); for
  # ULS it is the proper chi-square-distributed statistic (matching psych::fa(fm =
  # "minres")), not the raw least-squares residual sum of squares treated as Wishart. NA
  # for PAF, missing N, underidentified df, or a non-PD model-implied matrix (e.g. Heywood
  # cases), where the discrepancy is undefined. DWLS is also excluded: the ML discrepancy is
  # not its fit function -- the appropriate categorical statistic is the mean-and-variance-
  # adjusted chi-square from the full asymptotic covariance -- so the chi-square-derived block
  # is left undefined here rather than reported on an inapplicable scale.
  # The Bartlett multiplier goes non-positive for small N relative to the number of
  # variables, which would turn the discrepancy into a negative (meaningless) chi-square;
  # guard it here so the statistic falls through to the chi NA branch below.
  mult <- .bartlett_mult(N, m, q)
  if (!(method %in% c("PAF", "DWLS")) && !is.na(N) && df >= 0 && mult > 0) {
    Sigma <- LLt
    diag(Sigma) <- 1
    # discrepancy F (>= 0, clamped); the reported chi-square is F * mult.
    Fchi <- tryCatch({
      # F = tr(Sigma^-1 R) - log|Sigma^-1 R| - p, computed via a determinant split
      # (log|Sigma^-1 R| = log|R| - log|Sigma|) to avoid forming the explicit inverse
      # and to stay stable for ill-conditioned Sigma.
      ldSigma <- determinant(Sigma, logarithm = TRUE)
      ldR <- determinant(R, logarithm = TRUE)
      val <- sum(diag(solve(Sigma, R))) +
        as.numeric(ldSigma$modulus) - as.numeric(ldR$modulus) - m
      if (!is.finite(val) || ldSigma$sign <= 0 || ldR$sign <= 0) {
        NA_real_
      } else {
        # clamp the floating-point dust of a (near-)perfect fit to 0
        max(0, val)
      }
    }, error = function(e) NA_real_)
    chi <- Fchi * mult
  } else {
    chi <- NA_real_
    Fchi <- NA_real_
  }

  # null model: reuse the log-determinant already computed for the model chi-square above
  # instead of letting .null_chisq() recompute determinant(R). Undefined (NA) chi-square
  # leaves the whole chi-derived block NA.
  if (is.na(chi)) {
    chi <- NA
    chi_null <- NA
    df_null <- NA
    chi_cfi <- NA
    chi_null_cfi <- NA
  } else {
    chi_null <- .null_chisq(R, N, ld = ldR)
    df_null <- (m**2 - m) / 2
    # CFI/TLI compare the model and baseline noncentralities and so need both on one
    # scaling constant (Bentler, 1990; Tucker & Lewis, 1973). The reported model and
    # baseline chi-squares keep their own Bartlett corrections (the model matching
    # factanal, the baseline being Bartlett's test of sphericity), but the factor-count
    # term (2q)/3 sits only in the model multiplier and would bias the ratio. Put both on
    # the common (N - 1) scale -- the noncentrality scale RMSEA also uses -- for CFI/TLI.
    chi_cfi <- Fchi * (N - 1)
    chi_null_cfi <- -as.numeric(ldR$modulus) * (N - 1)
  }

  # CFI/TLI/RMSEA/AIC/BIC/ECVI and the p-values. Shared with the scaled-chi-square
  # (sandwich-SE) path, which supplies its own scaled model and baseline statistics
  # (already comparable, so it uses the CFI/TLI scale defaults).
  idx <- .chi_fit_indices(chi, df, chi_null, df_null, N, m, ci,
                          chi_cfi = chi_cfi, chi_null_cfi = chi_null_cfi)

  out <- list(
    chi = chi,
    df = df,
    p_chi = idx$p_chi,
    CAF = CAF,
    RMSR = RMSR,
    SRMR = SRMR,
    CFI = idx$CFI,
    TLI = idx$TLI,
    RMSEA = idx$RMSEA,
    RMSEA_LB = idx$RMSEA_LB,
    RMSEA_UB = idx$RMSEA_UB,
    AIC = idx$AIC,
    BIC = idx$BIC,
    ECVI = idx$ECVI,
    Fm = Fm,
    chi_null = chi_null,
    df_null = df_null,
    p_null = idx$p_null
  )

}
