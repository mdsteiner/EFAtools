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

# Independence-model (baseline) chi-square of the null model (model-implied matrix =
# identity), -log|R| * mult. With `corrected = TRUE` (default) mult is the Bartlett
# multiplier N - 1 - (2p + 5)/6: this is Bartlett's test of sphericity, the reported
# baseline statistic shared with BARTLETT() and the CFI/TLI baseline in .gof(). With
# `corrected = FALSE` mult is (N - 1), the uncorrected discrepancy scale on which the
# RMSEA noncentrality is built (the 0-factor reference in HULL() and the null model in SMT()).
.null_chisq <- function(R, N, ld = determinant(R, logarithm = TRUE),
                        corrected = TRUE) {
  p <- ncol(R)
  # Use the log-determinant directly: det(R) underflows to 0 for large,
  # near-singular (but still positive-definite) matrices, which would otherwise
  # turn the statistic into NA and propagate into .gof(), SMT(), and BARTLETT().
  # `ld` defaults to that decomposition; callers that already have it (e.g. .gof()
  # for the model chi-square) pass it in to avoid recomputing determinant(R).
  if (ld$sign <= 0 || !is.finite(ld$modulus)) return(NA_real_)
  # Guard the multiplier: for small N relative to p the Bartlett multiplier turns
  # non-positive, which would flip the baseline chi-square sign and propagate a
  # spurious statistic into .gof(), SMT(), and BARTLETT().
  mult <- if (corrected) .bartlett_mult(N, p) else N - 1
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

# RMSEA point estimate from a model chi-square on the (N - 1) noncentrality scale: the
# self-consistent Browne & Cudeck (1992) quantity sqrt(max(0, chi_cfi - df) / (df * (N - 1))).
# Callers guard df > 0 and N > 1 (a positive divisor); the cap at 1, where wanted, is applied by
# the caller. Shared by .chi_fit_indices(), HULL(), and the EFA_POOLED() D2 pooler so the RMSEA
# convention lives in one place.
.rmsea_point <- function(chi_cfi, df, N) {
  sqrt(max(0, chi_cfi - df) / (df * (N - 1)))
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
# in, so both report the block by exactly the same formulas. `chi_cfi` / `chi_null_cfi` are the
# model and baseline statistics on the common (N - 1) noncentrality scale used for CFI/TLI/RMSEA;
# they are required (not defaulted to `chi`/`chi_null`) so a caller cannot silently report those
# approximation indices off the Bartlett-corrected test scale.
.chi_fit_indices <- function(chi, df, chi_null, df_null, N, m, ci,
                             chi_cfi, chi_null_cfi) {

  if (is.na(chi)) {
    return(list(p_chi = NA_real_, CFI = NA_real_, TLI = NA_real_, RMSEA = NA_real_,
                RMSEA_LB = NA_real_, RMSEA_UB = NA_real_, AIC = NA_real_,
                BIC = NA_real_, ECVI = NA_real_, p_null = NA_real_))
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
  # path passes its scaled model and baseline statistics, which are already comparable).
  # When the baseline is undefined (e.g. a degenerate scaled baseline)
  # CFI/TLI are NA while the model-only RMSEA/AIC/BIC/p remain computable.
  inc <- .incremental_fit(chi_cfi, df, chi_null_cfi, df_null)
  CFI <- inc$CFI
  TLI <- inc$TLI

  if (df != 0) {

    # RMSEA point estimate and bounds (Steiger & Lind, 1980; Browne & Cudeck, 1992).
    # The noncentrality estimate chi_cfi - df and its divisor df * (N - 1) are both on the
    # uncorrected (N - 1) discrepancy scale, so the RMSEA is the self-consistent
    # Browne-Cudeck quantity and shares one noncentrality scale with CFI/TLI. The Bartlett
    # small-sample multiplier is applied only to the reported model chi-square test (matching
    # stats::factanal), not to this approximation index, for which it has no role.
    RMSEA <- .rmsea_point(chi_cfi, df, N)
    if (RMSEA > 1) RMSEA <- 1

    if (isTRUE(ci)) {

      # Analytic 90% RMSEA confidence bounds via the noncentrality solver (Browne &
      # Cudeck, 1992); skipped when ci = FALSE (e.g. per-replicate bootstrap fits).
      lambda_l <- .rmsea_lambda(chi_cfi, df, .95)
      lambda_u <- .rmsea_lambda(chi_cfi, df, .05)

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

# Two-stage / full-information ML likelihood-ratio chi-square for cor_method = "fiml". The EM
# moments `fiml` (the saturated mean `mu`, covariance `sigma`, and saturated observed-data
# log-likelihood `logl`, plus the raw `data`) come from .prepare_cor_input(); the analysed
# correlation R is cov2cor(sigma). The model-implied correlation R_model = LL' + diag(1 - h2)
# has a unit diagonal (diag(LL') = rowSums(L^2) = h2), so it is the model-implied LL' (the `LLt`
# .gof() already formed) with the diagonal set to 1; put it back on the covariance scale with the
# Stage-1 variances d = diag(sigma) as
# Sigma_model = D^(1/2) R_model D^(1/2). The model leaves the mean unrestricted (mu_model = mu),
# so chi = 2 (logl_sat - logl_model) and df = ((p - q)^2 - (p + q))/2 are unchanged (the mean
# and variance parameters are shared with the saturated model and cancel). The independence
# baseline (free means + free variances, zero covariances) factorises under MAR, so its MLE is
# each variable's available-case mean and ML (/n_j) variance; chi_null = 2 (logl_sat - logl_null)
# with df_null = p(p - 1)/2. No Bartlett correction applies to a likelihood-ratio statistic, so
# the CFI/TLI/RMSEA noncentrality scale uses the same statistics (chi_cfi = chi, chi_null_cfi =
# chi_null). This helper returns the plain two-stage likelihood-ratio statistic; `.gof()` reports
# the Satorra-Bentler-corrected two-stage statistic (.fiml_scaled_test()) in its place, and the
# plain LRT here stands only as the fallback when that correction cannot be formed -- referenced to
# chi^2(df) it is only approximate under the two-stage estimator (Yuan, Marshall, & Bentler, 2002,
# Psychometrika 67:95-121). Returns NA throughout for PAF (no discrepancy), missing N,
# underidentified df, or a non-positive-definite Sigma_model (e.g. a Heywood case), matching how
# the ML/ULS discrepancy block leaves those undefined.
.gof_fiml_chisq <- function(LLt, N, method, df, m, fiml) {

  if (!(method %in% c("ML", "ULS")) || is.na(N) || df < 0) {
    return(list(chi = NA_real_, chi_null = NA_real_, df_null = NA_real_,
                chi_cfi = NA_real_, chi_null_cfi = NA_real_))
  }

  # The saturated log-likelihood was accumulated by the EM over rows with at least one observed
  # value; re-apply the same filter and share one missingness-pattern grouping across the model
  # and baseline log-likelihoods so all three use the identical row set. The filter matters on
  # the point-estimate path, where fiml$data is the user's raw data and may keep a fully-missing
  # row the EM dropped when forming logl_sat (the bootstrap resample pool already excludes such
  # rows).
  data <- as.matrix(fiml$data)         # EFA() accepts a data frame; match .fiml_loglik()
  obs <- !is.na(data)
  keep <- rowSums(obs) > 0L
  data <- data[keep, , drop = FALSE]
  patterns <- .fiml_patterns(obs[keep, , drop = FALSE])
  logl_sat <- fiml$logl

  d <- diag(fiml$sigma)
  R_model <- LLt                       # the p x p product .gof() already formed
  diag(R_model) <- 1
  Sigma_model <- R_model * tcrossprod(sqrt(d))
  logl_model <- tryCatch(.fiml_loglik(data, fiml$mu, Sigma_model, patterns = patterns),
                         error = function(e) NA_real_)
  # A non-positive-definite Sigma_model (e.g. a Heywood case) leaves the model deviance
  # undefined; NA the whole chi-square block together, matching how the ML/ULS discrepancy path
  # NAs the model and baseline in lockstep when the discrepancy is undefined.
  if (is.na(logl_model)) {
    return(list(chi = NA_real_, chi_null = NA_real_, df_null = NA_real_,
                chi_cfi = NA_real_, chi_null_cfi = NA_real_))
  }
  chi <- max(0, 2 * (logl_sat - logl_model))

  # Independence baseline: per-variable available-case mean and ML (/n_j) variance. A diagonal
  # Sigma makes the joint observed-data log-likelihood equal the sum of the univariate marginals,
  # so this is the exact independence-model FIML log-likelihood.
  mu_null <- colMeans(data, na.rm = TRUE)
  var_null <- colMeans(sweep(data, 2L, mu_null, "-")^2, na.rm = TRUE)
  logl_null <- tryCatch(
    .fiml_loglik(data, mu_null, diag(var_null, nrow = m), patterns = patterns),
    error = function(e) NA_real_)
  chi_null <- if (is.na(logl_null)) NA_real_ else max(0, 2 * (logl_sat - logl_null))
  df_null <- m * (m - 1) / 2

  list(chi = chi, chi_null = chi_null, df_null = df_null,
       chi_cfi = chi, chi_null_cfi = chi_null)
}

.gof <- function(L, # The loading/ pattern matrix
                 R, # The correlation matrix
                 N, # The number of cases
                 method, # The estimation method
                 Fm, # Minimized error
                 ci = TRUE, # Compute the analytic RMSEA confidence bounds
                 fiml = NULL) { # EM moments + raw data for the FIML LRT (cor_method = "fiml")
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

  # Model and baseline chi-square. For cor_method = "fiml" both are likelihood-ratio statistics
  # computed from the EM moments and the raw data (.gof_fiml_chisq); otherwise the Bartlett-
  # corrected ML/ULS discrepancy below. Either way the residual indices above and the shared
  # .chi_fit_indices() call below stay the single source of the rest of the block.
  if (is.null(fiml)) {

    # Model chi-square: the Bartlett-corrected (Bartlett, 1951) ML discrepancy
    # F = tr(Sigma^-1 R) - log|Sigma^-1 R| - p, evaluated at the model-implied correlation
    # matrix Sigma = LL' (unit diagonal), times (N - 1 - (2p + 5)/6 - (2q)/3). For ML this
    # equals the ML objective times the Bartlett multiplier (matching stats::factanal); for
    # ULS the same ML/Wishart discrepancy is evaluated at the ULS-fitted Sigma (matching
    # psych::fa(fm = "minres")), rather than treating the raw least-squares residual sum of
    # squares as the statistic. Its chi-square reference distribution is asymptotically exact
    # under ML and is used here as the conventional approximation for ULS. NA
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
      # chi is already NA_real_ here; NA out the baseline and common-scale quantities to
      # match, so the whole undefined chi-square block is reported as a numeric (NA_real_) NA.
      chi_null <- NA_real_
      df_null <- NA_real_
      chi_cfi <- NA_real_
      chi_null_cfi <- NA_real_
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

  } else {

    fc <- .gof_fiml_chisq(LLt, N, method, df, m, fiml)
    chi <- fc$chi
    chi_null <- fc$chi_null
    df_null <- fc$df_null
    chi_cfi <- fc$chi_cfi
    chi_null_cfi <- fc$chi_null_cfi

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

  # cor_method = "fiml": replace the plain two-stage likelihood-ratio chi-square (computed above)
  # with the asymptotically-correct Satorra-Bentler-corrected two-stage statistic (Yuan, Marshall,
  # & Bentler, 2002), built on the saturated FIML asymptotic covariance. The plain LRT is referenced
  # to chi^2(df), which is not the two-stage estimator's reference distribution, so the p-value and
  # the noncentrality-based CFI/TLI/RMSEA are biased. Applied to every FIML fit (point estimate and
  # bootstrap replicates), independent of `se`; .fiml_scaled_test() returns NULL (PAF, a just/under-
  # identified model, or a degenerate/non-PD saturated covariance), where the NA / likelihood-ratio
  # block above stands.
  if (!is.null(fiml)) {
    # AIC/BIC/ECVI are likelihood-ratio information criteria with no standard interpretation under
    # the corrected two-stage statistic; leave them NA for every FIML fit. .apply_scaled_test() also
    # NAs them on the scaled path, but it does not run on the plain-LRT and just-identified (df == 0)
    # fallbacks, so NA them here too rather than ship a chi-square-derived value the two-stage
    # estimator cannot support.
    out$AIC <- NA_real_
    out$BIC <- NA_real_
    out$ECVI <- NA_real_
    st <- .fiml_scaled_test(L, R, N, method, df, m, fiml)
    if (!is.null(st)) out <- .apply_scaled_test(out, st, N)
  }

  out

}
