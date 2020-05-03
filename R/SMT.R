#' Sequential Model Tests for Chi Square and RMSEA lower bound
#'
#' Sequential Chi Square Model Tests (SMT) are a factor retention method where
#' multiple
#' EFAs with increasing numbers of factors are fitted and the number of factors
#' for which the Chi Square value first becomes non-significant is taken as the
#' suggested number of factors.
#' Preacher, Zhang, Kim, & Mels (2013) suggested a similiar approach with the
#' lower bound of the 90% confidence interval of the Root Mean Square Error of
#' Approximation (RMSEA; Browne & Cudeck, 1992; Steiger & Lind, 1980), where the
#' number of factors for which this lower bound first falls below .05 is the
#' suggested number of factors to retain.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations.
#' @param N numeric. The number of observations. Needs only be specified if a
#' correlation matrix is used.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
#' data is given as input. Default is "pairwise.complete.obs".
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
#' In comparison with other prominent factor retention criteria, SMT performed
#' well at determining the number of factors to extract in EFA (Auerswald &
#' Moshagen, 2019). The RMSEA lower bound performed well at determining the true
#' number of factors as well (Preacher, Zhang, Kim, & Mels, 2013).
#'
#' @return A list of class SMT containing
#' \item{nfac_chi}{The number of factors to retain according to the chi square.}
#' \item{nfac_RMSEA}{The number of factors to retain according to the RMSEA lowe
#' bound}
#' \item{p_null}{The p-value for the null model (zero factors)}
#' \item{ps_chi}{The p-values for EFA models with increasing numbers of factors,
#' starting with 1 factor}
#' \item{RMSEA.LBs}{The lower bounds of the 90% confidence interval for the RMSEA
#' for EFA models with increasing numbers of factors, starting with 1 factor}
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
#' @export
#'
#' @examples
#' SMT_base <- SMT(test_models$baseline$cormat, N = 500)
#' SMT_base
#'
SMT <- function(x, N = NULL, use = c("pairwise.complete.obs", "all.obs",
                                     "complete.obs", "everything",
                                     "na.or.complete")){

  use <- match.arg(use)

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    if(any(is.na(x))){

      stop("The correlation matrix you entered contains missing values.
           Eigenvalues cannot be computed.")

    }

    R <- x

  } else {

    message("x was not a correlation matrix. Correlations are found from entered
            raw data.")

    R <- stats::cor(x, use = use)
    colnames(R) <- colnames(x)
    N <- nrow(x)

  }

  if (is.null(N)) {

    stop("Argument 'N' was NULL. Either provide N or raw data.")

  }

  # Prepare objects for sequential tests
  max_fac <- min(ncol(x)-2, ncol(x)/2 + 2)
  ps <- vector("double", max_fac)
  RMSEA.LB <- vector("double", max_fac)

  # First check if 0 factors already result in nonsignificant chi square
  zeromod <- EFA(x, n_factors = 1, method = "ML", rotation = "none", N = N)
  p_null <- stats::pchisq(zeromod$fit_indices$chi_null,
                          zeromod$fit_indices$df_null, lower.tail = F)

  if(p_null > 0.05){

    nfac_chi <- 0

  } else {

    # sequentially perform EFAs with 1 to the maximum number of factors
    for (i in 1:max_fac) {

      temp <- EFA(x, n_factors = i, method = "ML", rotate ="none", N = N)
      ps[i] <- stats::pchisq(temp$fit_indices$chi, temp$fit_indices$df,
                             lower.tail = FALSE)
      RMSEA.LB[i] <- temp$fit_indices$RMSEA.LB

    }

    # With which number of factors does the chi square first become
    # non-significant?
    if(any(ps > 0.05, na.rm = TRUE)) {

      nfac_chi <- which(ps > 0.05)[1]

    } else {

      nfac_chi <- 0

    }

    # With which number of factors does the RMSEA first fall below .05?
    if(any(RMSEA.LB < .05, na.rm = TRUE)){

      nfac_RMSEA <- which(RMSEA.LB < .05)[1]

    } else {

      nfac_RMSEA <- 0

    }

  }

  # Prepare the output
  output <- list(nfac_chi = nfac_chi,
                 nfac_RMSEA = nfac_RMSEA,
                 p_null = p_null,
                 ps_chi = ps,
                 RMSEA.LBs = RMSEA.LB)

  class(output) <- "SMT"

  return(output)

}
