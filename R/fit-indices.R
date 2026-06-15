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
    # compute variance proportions
    if (ncol(L_rot) > 1) {
      vars <- colSums(L_rot^2)
    }
    else {
      vars <- sum(L_rot^2)
    }
  } else {
    # compute variance proportions
    vars <- diag(Phi %*% t(L_rot) %*% L_rot)
  }

  # Compute the explained variances. The code is based on the psych::fac() function
  # total variance (sum of communalities and uniquenesses)
  h2 <- diag(L_unrot %*% t(L_unrot))
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

# Independence-model (baseline) chi-square: the Bartlett-corrected ML discrepancy of the
# null model (model-implied matrix = identity), -log|R| * (N - 1 - (2p + 5)/6). Shared by
# .gof() and HULL so the model and baseline chi-square sit on the same scale.
.null_chisq <- function(R, N) {
  p <- ncol(R)
  # Use the log-determinant directly: det(R) underflows to 0 for large,
  # near-singular (but still positive-definite) matrices, which would otherwise
  # turn the statistic into NA and propagate into .gof(), SMT(), and BARTLETT().
  ld <- determinant(R, logarithm = TRUE)
  if (ld$sign <= 0 || !is.finite(ld$modulus)) return(NA_real_)
  # Guard the Bartlett multiplier: for small N relative to p it turns
  # non-positive, which would flip the baseline chi-square sign and propagate a
  # spurious statistic into .gof(), SMT(), and BARTLETT().
  mult <- N - 1 - (2 * p + 5) / 6
  if (mult <= 0) return(NA_real_)
  -as.numeric(ld$modulus) * mult
}

# Noncentrality parameter for an RMSEA confidence bound: solve pchisq(chi, df, ncp) = goal
# for ncp via stats::uniroot (the 90% CI lower bound uses goal = .95, the upper bound goal =
# .05; Browne & Cudeck, 1992). Returns 0 when the model chi-square already lies below the
# target quantile, so the bound collapses to 0. Shared by .gof() and SMT().
.rmsea_lambda <- function(chi, df, goal) {
  if (stats::pchisq(chi, df = df, ncp = 0) >= goal) {
    p_chi_fun <- function(x, val, df, goal) goal - stats::pchisq(val, df, ncp = x)
    stats::uniroot(f = p_chi_fun, interval = c(1e-10, 10000), val = chi, df = df,
                   goal = goal, extendInt = "upX", maxiter = 100L)$root
  } else {
    0
  }
}

.gof <- function(L, # The loading/ pattern matrix
                 R, # The correlation matrix
                 N, # The number of cases
                 method, # The estimation method
                 Fm) { # Minimized error
  m <- nrow(L)
  q <- ncol(L)

  # dfs
  df <- ((m - q)**2 - (m + q)) / 2

  ### compute CAF
  delta_hat <- R - (L %*% t(L))
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
  # cases), where the discrepancy is undefined.
  if (method != "PAF" && !is.na(N) && df >= 0) {
    Sigma <- L %*% t(L)
    diag(Sigma) <- 1
    chi <- tryCatch({
      # F = tr(Sigma^-1 R) - log|Sigma^-1 R| - p, computed via a determinant split
      # (log|Sigma^-1 R| = log|R| - log|Sigma|) to avoid forming the explicit inverse
      # and to stay stable for ill-conditioned Sigma.
      ldSigma <- determinant(Sigma, logarithm = TRUE)
      ldR <- determinant(R, logarithm = TRUE)
      Fchi <- sum(diag(solve(Sigma, R))) +
        as.numeric(ldSigma$modulus) - as.numeric(ldR$modulus) - m
      if (!is.finite(Fchi) || ldSigma$sign <= 0 || ldR$sign <= 0) {
        NA_real_
      } else {
        # The Bartlett multiplier goes non-positive for small N relative to the
        # number of variables, which would turn the discrepancy into a negative
        # (meaningless) chi-square; return NA so it falls through to the existing
        # chi NA branch instead of masquerading as perfect fit.
        mult <- N - 1 - (2 * m + 5) / 6 - (2 * q) / 3
        if (mult <= 0) {
          NA_real_
        } else {
          # clamp the floating-point dust of a (near-)perfect fit to 0
          max(0, Fchi) * mult
        }
      }
    }, error = function(e) NA_real_)
  } else {
    chi <- NA_real_
  }

  if (!is.na(chi)) {

    ### compute CFI

    # null model
    chi_null <- .null_chisq(R, N)
    df_null <- (m**2 - m) / 2
    p_null <- stats::pchisq(chi_null, df_null, lower.tail = F)

    # current model
    p_chi <- stats::pchisq(chi, df, lower.tail = F)

    ### compute CFI (Bentler, 1990): the model and baseline noncentralities
    # d = chi - df are each floored at 0 before the ratio, so an over-fitting
    # model (d <= 0) cannot deflate the index and CFI is 1 when the baseline
    # itself does not misfit (denominator 0). The flooring keeps CFI in [0, 1]
    # without a post-hoc clamp and matches lavaan::fitMeasures().
    d_m <- max(chi - df, 0)
    d_null <- max(chi_null - df_null, 0)
    CFI <- if (df == 0 || max(d_m, d_null) == 0) 1 else 1 - d_m / max(d_m, d_null)

    ### compute TLI (Tucker-Lewis index / NNFI; Tucker & Lewis, 1973)
    if (df == 0) {
      TLI <- 1
    } else {
      ratio_null <- chi_null / df_null
      # ratio_null == 1 (baseline chi-square equal to its df) leaves the index
      # undefined (0/0); report a perfect 1 as in the just-identified case.
      TLI <- if (ratio_null == 1) 1 else (ratio_null - chi / df) / (ratio_null - 1)
    }


    ### compute RMSEA, incl. 90% confidence intervals if df are not 0
    if(df != 0){

      # formula 12.6 from Kline 2016; Principles and practices of...
      RMSEA <- sqrt(max(0, chi - df) / (df * (N - 1)))

    lambda_l <- .rmsea_lambda(chi, df, .95)
    lambda_u <- .rmsea_lambda(chi, df, .05)

    RMSEA_LB <- sqrt(lambda_l / (df * (N - 1)))
    RMSEA_UB <- sqrt(lambda_u / (df * (N - 1)))

    if(RMSEA > 1) RMSEA <- 1
    if(RMSEA_LB > 1) RMSEA_LB <- 1
    if(RMSEA_UB > 1) RMSEA_UB <- 1

    } else {

      RMSEA <- 0
      RMSEA_LB <- 0
      RMSEA_UB <- 0

    }

    ### compute AIC and BIC based on chi square
    AIC <- chi - 2 * df
    BIC <- chi - log(N) * df

    ### compute ECVI (expected cross-validation index; Browne & Cudeck, 1989)
    n_params <- m * (m + 1) / 2 - df
    ECVI <- (chi + 2 * n_params) / (N - 1)

  } else {
    chi <- NA
    p_chi <- NA
    CFI <- NA
    TLI <- NA
    RMSEA <- NA
    RMSEA_LB <- NA
    RMSEA_UB <- NA
    AIC <- NA
    BIC <- NA
    ECVI <- NA
    chi_null <- NA
    df_null <- NA
    p_null <- NA
  }

  out <- list(
    chi = chi,
    df = df,
    p_chi = p_chi,
    CAF = CAF,
    RMSR = RMSR,
    SRMR = SRMR,
    CFI = CFI,
    TLI = TLI,
    RMSEA = RMSEA,
    RMSEA_LB = RMSEA_LB,
    RMSEA_UB = RMSEA_UB,
    AIC = AIC,
    BIC = BIC,
    ECVI = ECVI,
    Fm = Fm,
    chi_null = chi_null,
    df_null = df_null,
    p_null = p_null
  )

}
