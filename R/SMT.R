#' Sequential Chi Square Model Tests, RMSEA lower bound, and AIC
#'
#' Sequential Chi Square Model Tests (SMT) are a factor retention method where
#' multiple
#' EFAs with increasing numbers of factors are fitted and the number of factors
#' for which the Chi Square value first becomes non-significant is taken as the
#' suggested number of factors.
#'
#' Preacher, Zhang, Kim, & Mels (2013) suggested a similiar approach with the
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
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
#' data is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. Passed to \code{\link[stats:cor]{stats::cor}}.
#'  Default is "pearson".
#'
#' @details
#' As a first step in the procedure, a maximum number of factors to extract is
#' determined with min(nvar - 2, nvar / 2 + 2), where nvar is the number of
#' variables or indicators.
#'
#' Then, EFAs with increasing numbers of factors from 1 to the maximum number are
#' fitted with maximum likelihood estimation.
#'
#' For the SMT, first the significance of the chi
#' square value for a model with 0 factors is determined. If this value is
#' not significant, 0 factors are suggested to retain. If it is not significant,
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
#' for the model where the AIC is lowest.
#'
#' In comparison with other prominent factor retention criteria, SMT performed
#' well at determining the number of factors to extract in EFA (Auerswald &
#' Moshagen, 2019). The RMSEA lower bound performed well at determining the true
#' number of factors as well, while the AIC performed well at determining the
#' most generalizable model (Preacher, Zhang, Kim, & Mels, 2013).
#'
#' The \code{SMT} function can also be called together with other factor
#' retention criteria in the \code{\link{N_FACTORS}} function.
#'
#' @return A list of class SMT containing
#' \item{nfac_chi}{The number of factors to retain according to the chi square.}
#' \item{nfac_RMSEA}{The number of factors to retain according to the RMSEA lower
#' bound}
#' \item{nfac_AIC}{The number of factors to retain according to the AIC}
#' \item{p_null}{The p-value for the null model (zero factors)}
#' \item{ps_chi}{The p-values for EFA models with increasing numbers of factors,
#' starting with 1 factor}
#' \item{RMSEA_LBs}{The lower bounds of the 90% confidence interval for the RMSEA
#' for EFA models with increasing numbers of factors, starting with 1 factor}
#' \item{AICs}{The AICs for EFA models with increasing numbers of factors,
#' starting with 1 factor}
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
#' @seealso Other factor retention criteria: \code{\link{CD}}, \code{\link{EKC}},
#' \code{\link{HULL}}, \code{\link{KGC}}, \code{\link{PARALLEL}}
#'
#' \code{\link{N_FACTORS}} as a wrapper function for this and all the
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
  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data."))

  }

  checkmate::assert_count(N, na.ok = TRUE)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    R <- x

  } else {

    cli::cli_alert_info(cli::col_cyan("'x' was not a correlation matrix. Correlations are found from entered raw data."))

    R <- stats::cor(x, use = use, method = cor_method)
    colnames(R) <- colnames(x)
    N <- nrow(x)

  }

  if (is.na(N)) {

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Argument 'N' was NA. Either provide N or raw data."))

  }

  # Prepare objects for sequential tests
  max_fac <- floor(ncol(x)/2)
  ps <- vector("double", max_fac)
  RMSEA_LB <- vector("double", max_fac)
  AIC <- vector("double", max_fac)

  # First check if 0 factors already result in nonsignificant chi square
  zeromod <- suppressWarnings(suppressMessages(EFA(x, n_factors = 1,
                                                   method = "ML",
                                                   rotation = "none", N = N)))
  p_null <- stats::pchisq(zeromod$fit_indices$chi_null,
                          zeromod$fit_indices$df_null, lower.tail = F)

  if(p_null > 0.05){

    nfac_chi <- 0

  } else {

    # sequentially perform EFAs with 1 to the maximum number of factors
    for (i in 1:max_fac) {

      temp <- suppressWarnings(suppressMessages(EFA(x, n_factors = i,
                                                    method = "ML",
                                                    rotate ="none", N = N)))
      ps[i] <- stats::pchisq(temp$fit_indices$chi, temp$fit_indices$df,
                             lower.tail = FALSE)
      RMSEA_LB[i] <- temp$fit_indices$RMSEA_LB
      AIC[i] <- temp$fit_indices$AIC

    }

    # With which number of factors does the chi square first become
    # non-significant?
    if(any(ps > 0.05, na.rm = TRUE)) {

      nfac_chi <- which(ps > 0.05)[1]

    } else {

      nfac_chi <- 0

    }

    # With which number of factors does the RMSEA first fall below .05?
    if(any(RMSEA_LB < .05, na.rm = TRUE)){

      nfac_RMSEA <- which(RMSEA_LB < .05)[1]

    } else {

      nfac_RMSEA <- 0

    }

    # With which number of factors is the AIC lowest?
    nfac_AIC <- which(AIC == min(AIC))

  }

  # Prepare the output
  output <- list(nfac_chi = nfac_chi,
                 nfac_RMSEA = nfac_RMSEA,
                 nfac_AIC = nfac_AIC,
                 p_null = p_null,
                 ps_chi = ps,
                 RMSEA_LBs = RMSEA_LB,
                 AICs = AIC,
                 settings = list(N = N,
                                 use = use,
                                 cor_method = cor_method))

  class(output) <- "SMT"

  return(output)

}
