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
#' @param cor_method character. Correlation computed from raw data: `"pearson"`,
#'   `"spearman"`, or `"kendall"` (passed to [stats::cor()]), or `"poly"` /
#'   `"tetra"` for polychoric / tetrachoric correlations of ordinal / binary data
#'   (a two-step estimator).
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
#' the Kaiser-Meyer-Olkin criterion ([efa_kmo()])
#' can still be used.
#'
#' This function was heavily influenced by the [psych::cortest.bartlett()] function from the psych package.
#'
#' The `efa_bartlett` function can also be called together with the
#'  ([efa_kmo()]) function and with factor retention criteria
#'  in the [efa_retain()] function.
#'
#' @return A list containing
#' \item{chisq}{The chi square statistic, or `NA` if `N` is too small for the
#'   Bartlett correction (i.e. \eqn{N - 1 - (2p + 5)/6 \le 0}).}
#' \item{p_value}{The p value of the chi square statistic, or `NA` when `chisq`
#'   is `NA`.}
#' \item{df}{The degrees of freedom for the chi square statistic.}
#' \item{settings}{A list of the settings used.}
#'
#' @source Bartlett, M. S. (1951). The effect of standardization on a Chi-square
#' approximation in factor analysis. Biometrika, 38, 337-344.
#'
#' @family factor analysis suitability
#'
#' @seealso [efa_kmo()] for another measure to determine
#'  suitability for factor analysis.
#'
#'  [efa_retain()] as a wrapper function for this function,
#'  [efa_kmo()] and several factor retention criteria.
#'
#' @export
#'
#' @examples
#' efa_bartlett(test_models$baseline$cormat, N = 500)
#'
efa_bartlett <- function(x, N = NA, use = c("pairwise.complete.obs", "all.obs",
                                            "complete.obs", "everything",
                                            "na.or.complete"),
                         cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")){

  # Perform argument checks
  .assert_cor_input(x)

  use <- .match_arg_ci(use)
  cor_method <- .match_arg_ci(cor_method)
  checkmate::assert_count(N, na.ok = TRUE)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(
    x, N = N, use = use, cor_method = cor_method, N_policy = "required",
    singular_tail = "Bartlett's test cannot be executed",
    N_required_msg = c("{.arg N} is {.val NA}; Bartlett's test could not be executed.",
                       "i" = "Provide {.arg N} or raw data."))
  R <- prep$R
  N <- prep$N

  # Calculate test statistic. Bartlett's test of sphericity is the likelihood-ratio
  # test of the independence model, i.e. the null-model chi-square shared with the
  # CFI baseline in .gof().
  p <- nrow(R)
  statistic <- .null_chisq(R, N)
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

  # the trailing legacy class is load-bearing: it keeps `inherits(x, "BARTLETT")`
  # working in code written against the superseded name. Do not drop it.
  class(output) <- c("efa_bartlett", "BARTLETT")

  return(output)

}
