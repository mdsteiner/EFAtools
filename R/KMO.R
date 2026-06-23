#' Kaiser-Meyer-Olkin criterion
#'
#' This function computes the Kaiser-Meyer-Olkin (KMO) criterion overall and for
#' each variable in a correlation matrix. The KMO represents the degree to
#' which each observed variable is predicted by the other variables in the
#' dataset and with this indicates the suitability for factor analysis.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#'  correlations.
#' @param use character. Passed to [stats::cor()] if raw
#'  data is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. Correlation computed from raw data: `"pearson"`,
#'   `"spearman"`, or `"kendall"` (passed to [stats::cor()]), or `"poly"` /
#'   `"tetra"` for polychoric / tetrachoric correlations of ordinal / binary data
#'   (a two-step estimator with no empty-cell continuity correction).
#' Default is "pearson".
#'
#' @details Kaiser (1970) proposed this index, originally called measure of
#' sampling adequacy (MSA), that indicates how near the inverted correlation
#' matrix \eqn{R^{-1}} is to a diagonal matrix \eqn{S} to determine a given
#' correlation matrix's (\eqn{R}) suitability for factor analysis.
#' The index is
#' \deqn{KMO = \frac{\sum_{i \neq j} r_{ij}^2}{\sum_{i \neq j} r_{ij}^2 + \sum_{i \neq j} q_{ij}^2}}
#' with \eqn{Q = SR^{-1}S} and S = \eqn{(diag R^{-1})^{-1/2}} where
#' \eqn{\sum_{i \neq j} r_{ij}^2} is the sum of squares of the off-diagonal
#' elements of \eqn{R} and \eqn{\sum_{i \neq j} q_{ij}^2} is the sum of squares of
#' the off-diagonal elements of \eqn{Q} (see also Cureton & D'Agostino, 1983).
#'
#' So KMO varies between 0 and 1, with larger values indicating higher suitability
#' for factor analysis. Kaiser and Rice (1974) suggest that KMO should at least
#' exceed .50 for a correlation matrix to be suitable for factor analysis.
#'
#' This function was heavily influenced by the [psych::KMO()]
#' function.
#'
#' See also [BARTLETT()] for another test of suitability for factor
#' analysis.
#'
#' The `KMO` function can also be called together with the
#' [BARTLETT()] function and with factor retention criteria in the
#'  [N_FACTORS()] function.
#'
#' @return A list containing
#' \item{KMO}{Overall KMO.}
#' \item{KMO_i}{KMO for each variable.}
#' \item{settings}{A list of the settings used.}
#'
#' @export
#'
#' @source Kaiser, H. F. (1970). A second generation little jiffy. Psychometrika,
#' 35, 401-415.
#' @source Kaiser, H. F. & Rice, J. (1974). Little jiffy, mark IV. Educational
#' and Psychological Measurement, 34, 111-117.
#' @source Cureton, E. E. & D'Agostino, R. B. (1983). Factor analysis: An
#'  applied approach. Hillsdale, N.J.: Lawrence Erlbaum Associates, Inc.
#'
#' @seealso [BARTLETT()] for another measure to determine
#' suitability for factor analysis.
#'
#' [N_FACTORS()] as a wrapper function for this function,
#' [BARTLETT()] and several factor retention criteria.
#'
#' @examples
#' KMO(test_models$baseline$cormat)
KMO <- function(x, use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                           "everything", "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")) {

  # Perform argument checks
  .assert_cor_input(x)

  use <- match.arg(use)
  cor_method <- match.arg(cor_method)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, use = use, cor_method = cor_method,
                             N_policy = "none")
  R <- prep$R

  KMO_list <- .compute_kmo(R)
  KMO <- KMO_list$KMO
  KMO_i <- KMO_list$KMO_i

  if(!is.null(colnames(R))){

    names(KMO_i) <- colnames(R)

  } else if(!is.null(rownames(R))) {

    names(KMO_i) <- rownames(R)

  } else {

    names(KMO_i) <- paste0("V", seq_len(ncol(R)))

  }

  # Prepare settings
  settings <- list(use = use,
                   cor_method = cor_method)

  output <- list(KMO = KMO,
                 KMO_i = KMO_i,
                 settings = settings)

  class(output) <- "KMO"

  return(output)

}

# outsourced computation of KMO to
# use it when computing CAF which needs
# different input checks
.compute_kmo <- function(R) {

  # Check if correlation matrix is invertable, if it is not, stop with message
  R_i <- try(solve(R), silent = TRUE)

  if (inherits(R_i, "try-error")) {
    cli::cli_abort("The matrix is singular; KMO cannot be computed.",
                   class = "efa_cor_singular")
  }

  # Start computations
  S <- diag(diag(R_i)^(-0.5))
  Q <- S %*% R_i %*% S
  diag(Q) <- 0
  diag(R) <- 0
  sumQ2 <- sum(Q^2)
  sumR2 <- sum(R^2)
  KMO <- sumR2/(sumR2 + sumQ2)
  KMO_i <- colSums(R^2)/(colSums(R^2) + colSums(Q^2))
  list(
    KMO = KMO,
    KMO_i = KMO_i
  )
}
