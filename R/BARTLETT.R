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
#' @param use character. Passed to [stats::cor()] if raw data
#' is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. Passed to [stats::cor()].
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
#' the Kaiser-Meyer-Olkin criterion ([EFAtools::KMO()])
#' can still be used.
#'
#' This function was heavily influenced by the [psych::cortest.bartlett()] function from the psych package.
#'
#' The `BARTLETT` function can also be called together with the
#'  ([EFAtools::KMO()]) function and with factor retention criteria
#'  in the [N_FACTORS()] function.
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
#' @seealso [EFAtools::KMO()] for another measure to determine
#'  suitability for factor analysis.
#'
#'  [N_FACTORS()] as a wrapper function for this function,
#'  [EFAtools::KMO()] and several factor retention criteria.
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
  .assert_cor_input(x)

  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  checkmate::assert_count(N, na.ok = TRUE)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(
    x, N = N, use = use, cor_method = cor_method, N_policy = "required",
    singular_tail = "Bartlett's test cannot be executed",
    N_required_msg = c("{.arg N} is {.val NA}; Bartlett's test could not be executed.",
                       "i" = "Provide {.arg N} or raw data."))
  R <- prep$R
  N <- prep$N

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
