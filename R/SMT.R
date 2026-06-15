#' Sequential Chi Square Model Tests, RMSEA lower bound, and AIC
#'
#' Sequential Chi Square Model Tests (SMT) are a factor retention method where
#' multiple
#' EFAs with increasing numbers of factors are fitted and the number of factors
#' for which the Chi Square value first becomes non-significant is taken as the
#' suggested number of factors.
#' Preacher, Zhang, Kim, & Mels (2013) suggested a similar approach with the
#' lower bound of the 90% confidence interval of the Root Mean Square Error of
#' Approximation (RMSEA; Browne & Cudeck, 1992; Steiger & Lind, 1980), and with
#' the Akaike Information Criterion (AIC). For the RMSEA, the
#' number of factors for which this lower bound first falls below .05 is the
#' suggested number of factors to retain. For the AIC, it is the number of factors
#' where the AIC is lowest.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations.
#' @param N numeric. The number of observations. Needs only be specified if a
#' correlation matrix is used.
#' @param use character. Passed to [stats::cor()] if raw
#' data is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. Passed to [stats::cor()].
#'  Default is "pearson".
#'
#' @details
#' As a first step in the procedure, a maximum number of factors to extract is
#' determined for which the model is still over-identified (df > 0).
#'
#' Then, EFAs with increasing numbers of factors from 1 to the maximum number are
#' fitted with maximum likelihood estimation.
#'
#' For the SMT, first the significance of the chi
#' square value for a model with 0 factors is determined. If this value is
#' not significant, 0 factors are suggested to retain. If it is significant,
#' a model with 1 factor is estimated and the significance of its chi square value
#' is determined, and so on, until a non-significant result is obtained. The
#' suggested number of factors is the number of factors for the model where the
#' chi square value first becomes non-significant.
#'
#' Regarding the RMSEA, the suggested number of factors is the number of factors
#' for the model where the lower bound of the 90% confidence interval of the
#' RMSEA first falls below the .05 threshold.
#'
#' Regarding the AIC, the suggested number of factors is the number of factors
#' for the model with the lowest AIC.
#'
#' In comparison with other prominent factor retention criteria, SMT performed
#' well at determining the number of factors to extract in EFA (Auerswald &
#' Moshagen, 2019). The RMSEA lower bound also performed well at determining the true
#' number of factors, while the AIC performed well at determining the
#' most generalizable model (Preacher, Zhang, Kim, & Mels, 2013).
#'
#' The `SMT` function can also be called together with other factor
#' retention criteria in the [N_FACTORS()] function.
#'
#' @returns An object of class `efa_retention` (see [print.efa_retention()] for
#'   the print method). Its main fields are:
#' \item{n_factors}{A named numeric vector (`"chi"`, `"RMSEA"`, `"AIC"`) with the
#'   suggested number of factors from the sequential chi-square model tests, the
#'   RMSEA lower bound, and the AIC.}
#' \item{results}{A list with one record per criterion, each holding the criterion
#'   values for the null model (zero factors) through the maximum number of
#'   factors.}
#' \item{settings}{A list of the settings used.}
#'
#' @source Auerswald, M., & Moshagen, M. (2019). How to determine the number of
#' factors to retain in exploratory factor analysis: A comparison of extraction
#' methods under realistic conditions. Psychological Methods, 24(4), 468–491.
#' https://doi.org/10.1037/met0000200
#' @source Browne, M.W., & Cudeck, R. (1992). Alternative ways of assessing model
#' fit. Sociological Methods and Research, 21, 230–258.
#' @source Preacher, K. J., Zhang G., Kim, C., & Mels, G. (2013). Choosing the
#' Optimal Number of Factors in Exploratory Factor Analysis: A Model Selection
#' Perspective, Multivariate Behavioral Research, 48(1), 28-56,
#' doi:10.108/00273171.2012.710386
#' @source Steiger, J. H., & Lind, J. C. (1980, May). Statistically based tests
#' for the number of common factors. Paper presented at the annual meeting of
#' the Psychometric Society, Iowa City, IA.
#'
#' @seealso Other factor retention criteria: [CD()], [EKC()],
#' [HULL()], [KGC()], [PARALLEL()]
#'
#' [N_FACTORS()] as a wrapper function for this and all the
#' above-mentioned factor retention criteria.
#'
#' @export
#'
#' @examples
#' SMT_base <- SMT(test_models$baseline$cormat, N = 500)
#' SMT_base
#'
SMT <- function(x, N = NA, use = c("pairwise.complete.obs", "all.obs",
                                     "complete.obs", "everything",
                                     "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall")){

  # Perform argument checks
  .assert_cor_input(x)

  checkmate::assert_count(N, na.ok = TRUE)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(
    x, N = N, use = use, cor_method = cor_method, N_policy = "required",
    N_required_msg = c("{.arg N} is {.val NA}.",
                       "i" = "Provide {.arg N} or raw data."))
  R <- prep$R
  N <- prep$N

  # Prepare objects for sequential tests
  max_fac <- .det_max_factors(ncol(R))

  if(max_fac <= 0){
    cli::cli_abort(
      c("The model is underidentified, or just identified with a single factor; SMTs cannot be performed.",
        "i" = "Provide more indicators."),
      class = "efa_smt_underidentified"
    )
  }

  ps <- vector("double", max_fac)
  RMSEA_LB <- vector("double", max_fac)
  AIC <- vector("double", max_fac)
  nfac_chi <- NA
  nfac_RMSEA <- NA
  nfac_AIC <- NA

    # sequentially perform EFAs with 1 to the maximum number of factors
    for (i in seq_len(max_fac)) {

      temp <- suppressWarnings(suppressMessages(EFA(R, n_factors = i,
                                                    method = "ML",
                                                    rotate ="none", N = N)))
      ps[i] <- temp$fit_indices$p_chi
      RMSEA_LB[i] <- temp$fit_indices$RMSEA_LB
      AIC[i] <- temp$fit_indices$AIC

    }

  # With which number of factors does the chi square first become
  # non-significant?

  # Null-model (zero-factor) statistics. These depend only on R and N (the
  # model-implied matrix is the identity), so compute them directly rather than
  # reading them off the last fitted EFA, whose fit_indices come back NA when
  # that model is degenerate (e.g. a Heywood / non-positive-definite case).
  # Mirrors the null-model block in .gof().
  m <- ncol(R)
  chi_null <- .null_chisq(R, N)
  df_null <- (m^2 - m) / 2
  p_null <- stats::pchisq(chi_null, df_null, lower.tail = FALSE)

  # First check if 0 factors already result in nonsignificant chi square
    if(isTRUE(p_null > 0.05)){

      nfac_chi <- 0

    } else if(any(ps > 0.05, na.rm = TRUE)) {

      nfac_chi <- which(ps > 0.05)[1]

    } else {

      nfac_chi <- NA

    }

  # Calculate RMSEA (incl. lower bound of 90% CI) and AIC for the null model
  # (chi_null and df_null were computed above from R and N).
  RMSEA_LB_null <- sqrt(.rmsea_lambda(chi_null, df_null, .95) / (df_null * (N - 1)))

  AIC_null <- chi_null - 2 * df_null

  # With which number of factors does the RMSEA lower bound first fall below .05?
  if(isTRUE(RMSEA_LB_null < .05)){

    nfac_RMSEA <- 0

  } else {

    if(any(RMSEA_LB < .05, na.rm = TRUE)){

      nfac_RMSEA <- which(RMSEA_LB < .05)[1]

    } else {

      nfac_RMSEA <- NA

    }

  }

  # With which number of factors is the AIC lowest? (which.min returns the first
  # minimum, so ties yield a single, well-defined suggestion)
  AIC_all <- c(AIC_null, AIC)
  nfac_AIC <- which.min(AIC_all) - 1

  # one record per criterion (values for the null model through max_fac factors)
  results <- list(
    list(name = "chi", label = "Sequential chi-square model tests",
         n_factors = nfac_chi, plot_type = "none",
         x = 0:max_fac, y = c(p_null, ps)),
    list(name = "RMSEA", label = "Lower bound of RMSEA 90% CI",
         n_factors = nfac_RMSEA, plot_type = "none",
         x = 0:max_fac, y = c(RMSEA_LB_null, RMSEA_LB)),
    list(name = "AIC", label = "Akaike Information Criterion",
         n_factors = nfac_AIC, plot_type = "none",
         x = 0:max_fac, y = AIC_all)
  )

  output <- .new_efa_retention(
    "SMT",
    results = results,
    settings = list(N = N, use = use, cor_method = cor_method)
  )

  return(output)

}
