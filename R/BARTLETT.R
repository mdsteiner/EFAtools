#' Bartletts test of spericity
#'
#' This function tests whether a correlation matrix is significantly different
#' from an identity matrix (Bartlett, 1951). If the Bartletts test is not
#' significant, the correlation matrix is not suitable for factor analysis
#' because the variables show too little covariance.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations.
#' @param N numeric. The number of observations. Needs only be specified if a
#' correlation matrix is used.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
#' is given as input. Default is "pairwise.complete.obs".
#'
#' @details Bartlett (1951) proposed this statistic to determine a correlation
#' matrix' suitability for factor analysis (see also
#' \code{\link[EFAtools:KMO]{EFAtools:KMO}}). The statistic is approximately
#' chi square distributed with \eqn{df = \fraq{p(p - 1)}} and is given by
#'
#' \deqn{chi^2 = -log(det(R)) (N - 1 - (2 * p + 5)/6)}
#'
#' where \eqn{det(R)} is the determinant of the correlation matrix, \eqn{N} is
#' the sample size, and \eqn{p} is the number of variables.
#'
#' This function was heavily influenced by the
#' \code{\link[psych:cortest.bartlett]{cortest.bartlett}} function
#' from the psych package.
#'
#' @return A list containing
#' \item{chisq}{The chi square statistic.}
#' \item{p_value}{The p value of the chi square statistic.}
#' \item{df}{The degrees of freedom for the chi square statistic.}
#'
#' @source Bartlett, M. S. (1951). The effect of standardization on a Chi-square
#' approximation in factor analysis. Biometrika, 38, 337-344.
#'
#' @export
#'
#' @examples
#' BARTLETT(test_models$baseline$cormat, N = 500)
#'
BARTLETT <- function(x, N = NA, use = c("pairwise.complete.obs", "all.obs",
                                        "complete.obs", "everything",
                                        "na.or.complete")){

  use <- match.arg(use)

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    R <- x

    if (is.na(N)) {

      stop("Argument 'N' was NA, Bartletts test could not be executed.
           Please provide either N or raw data.")

    }

  } else {

    message("x was not a correlation matrix. Correlations are found from entered
            raw data.")

    R <- stats::cor(x, use = use)
    colnames(R) <- colnames(x)
    N <- nrow(x)

  }

  # Check if correlation matrix is invertable, if it is not, stop with message
  R_i <- try(solve(R))

  if (class(R_i) == "try-error") {
    stop("Correlation matrix is singular, Bartletts test cannot be executed")
  }

  # Check if correlation matrix is positive definite
  if(any(eigen(R, symmetric = TRUE, only.values = TRUE)$values <= 0)){

    R <- psych::cor.smooth(R)

  }

  # Calculate test statistic
  p <- nrow(R)
  detR <- det(R)
  statistic <- -log(detR) * (N - 1 - (2 * p + 5)/6)
  df <- p * (p - 1)/2
  pval <- pchisq(statistic, df, lower.tail = FALSE)

  # prepare the output
  output <- list(chisq = statistic, p_value = pval, df = df)

  class(output) <- "BARTLETT"

  return(output)

}
