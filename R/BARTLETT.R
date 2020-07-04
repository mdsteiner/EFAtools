#' Bartlett's test of sphericity
#'
#' This function tests whether a correlation matrix is significantly different
#' from an identity matrix (Bartlett, 1951). If the Bartlett's test is not
#' significant, the correlation matrix is not suitable for factor analysis
#' because the variables show too little covariance.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations.
#' @param N numeric. The number of observations. Needs only be specified if a
#' correlation matrix is used.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
#' is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. Passed to \code{\link[stats:cor]{stats::cor}}.
#' Default is "pearson".
#'
#' @details Bartlett (1951) proposed this statistic to determine a correlation
#' matrix' suitability for factor analysis. The statistic is approximately
#' chi square distributed with \eqn{df = \frac{p(p - 1)}{2}} and is given by
#'
#' \deqn{chi^2 = -log(det(R)) (N - 1 - (2 * p + 5)/6)}
#'
#' where \eqn{det(R)} is the determinant of the correlation matrix, \eqn{N} is
#' the sample size, and \eqn{p} is the number of variables.
#'
#' This tests requires multivariate normality. If this condition is not met,
#' the Kaiser-Meyer-Olkin criterion (\code{\link[EFAtools]{KMO}})
#' can still be used.
#'
#' This function was heavily influenced by the \code{\link[psych:cortest.bartlett]{psych::cortest.bartlett}} function from the psych package.
#'
#' The \code{BARTLETT} function can also be called together with the
#'  (\code{\link[EFAtools]{KMO}}) function and with factor retention criteria
#'  in the \code{\link{N_FACTORS}} function.
#'
#' @return A list containing
#' \item{chisq}{The chi square statistic.}
#' \item{p_value}{The p value of the chi square statistic.}
#' \item{df}{The degrees of freedom for the chi square statistic.}
#' \item{settings}{A list of the settings used.}
#'
#' @source Bartlett, M. S. (1951). The effect of standardization on a Chi-square
#' approximation in factor analysis. Biometrika, 38, 337-344.
#'
#' @seealso \code{\link[EFAtools]{KMO}} for another measure to determine
#'  suitability for factor analysis.
#'
#'  \code{\link{N_FACTORS}} as a wrapper function for this function,
#'  \code{\link[EFAtools]{KMO}} and several factor retention criteria.
#'
#' @export
#'
#' @examples
#' BARTLETT(test_models$baseline$cormat, N = 500)
#'
BARTLETT <- function(x, N = NA, use = c("pairwise.complete.obs", "all.obs",
                                        "complete.obs", "everything",
                                        "na.or.complete"),
                     cor_method = c("pearson", "spearman", "kendall")){

  # Perform argument checks
  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n"))

  }

  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  checkmate::assert_count(N, na.ok = TRUE)

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    R <- x

    if (is.na(N)) {

      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Argument 'N' was NA, Bartlett's test could not be executed. Please provide either N or raw data.\n"))

    }

  } else {

    message(cli::col_cyan(cli::symbol$info, " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n"))

    if (!is.na(N)) {
      warning(crayon::yellow$bold("!"), crayon::yellow(" 'N' was set and data entered. Taking N from data.\n"))
    }

    R <- stats::cor(x, use = use, method = cor_method)
    colnames(R) <- colnames(x)
    N <- nrow(x)

  }

  # Check if correlation matrix is invertable, if it is not, stop with message
  R_i <- try(solve(R), silent = TRUE)

  if (inherits(R_i, "try-error")) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Correlation matrix is singular, Bartlett's test cannot be executed.\n"))
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
  pval <- stats::pchisq(statistic, df, lower.tail = FALSE)

  # prepare the output
  settings <- list(N = N,
                   use = use,
                   cor_method = cor_method)

  output <- list(chisq = statistic,
                 p_value = pval,
                 df = df,
                 settings = settings)

  class(output) <- "BARTLETT"

  return(output)

}
