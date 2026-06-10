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

  if (method != "PAF" && !is.na(N) && df >=0) {

    ### compute CFI

    # null model
    chi_null <- sum(R[upper.tri(R)] ^ 2) * (N - 1)
    df_null <- (m**2 - m) / 2
    delta_hat_null <- chi_null - df_null
    p_null <- stats::pchisq(chi_null, df_null, lower.tail = F)

    # current model
    # formula 12.1 from Kline 2016; Principles and practices of...
    # This should also hold for ULS solutions -> Bentler & Bonett (1980)
    # Significance Tests  and Goodness  of Fit in the Analysis of Covariance
    # Structures. Psychological Bulletin, 88(3), 588-606
    chi <- Fm * (N - 1)
    p_chi <- stats::pchisq(chi, df, lower.tail = F)
    delta_hat_m <- chi - df
    CFI <- (delta_hat_null - delta_hat_m) / delta_hat_null
    if (CFI > 1 || df == 0) CFI <- 1
    if (CFI < 0) CFI <- 0


    ### compute RMSEA, incl. 90% confidence intervals if df are not 0
    if(df != 0){

      # formula 12.6 from Kline 2016; Principles and practices of...
      RMSEA <- sqrt(max(0, chi - df) / (df * (N - 1)))

    p_chi_fun <- function(x, val, df, goal){goal - stats::pchisq(val, df, ncp = x)}

    if (stats::pchisq(chi, df = df, ncp = 0) >= .95) {
      lambda_l <- stats::uniroot(f = p_chi_fun, interval = c(1e-10, 10000), val = chi,
                          df = df, goal = .95, extendInt = "upX",
                          maxiter = 100L)$root
    } else {
      lambda_l <- 0
    }

    if (stats::pchisq(chi, df = df, ncp = 0) >= .05) {
      lambda_u <- stats::uniroot(f = p_chi_fun, interval = c(1e-10, 10000),
                          val = chi, df = df, goal = .05,
                          extendInt = "upX", maxiter = 100L)$root
    }
    else {
      lambda_u <- 0
    }

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

  } else {
    chi <- NA
    p_chi <- NA
    CFI <- NA
    RMSEA <- NA
    RMSEA_LB <- NA
    RMSEA_UB <- NA
    AIC <- NA
    BIC <- NA
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
    CFI = CFI,
    RMSEA = RMSEA,
    RMSEA_LB = RMSEA_LB,
    RMSEA_UB = RMSEA_UB,
    AIC = AIC,
    BIC = BIC,
    Fm = Fm,
    chi_null = chi_null,
    df_null = df_null,
    p_null = p_null
  )

}
