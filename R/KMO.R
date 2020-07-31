#' Kaiser-Meyer-Olkin criterion
#'
#' This function computes the Kaiser-Meyer-Olkin (KMO) criterion overall and for
#' each variable in a correlation matrix. The KMO represents the degree to
#' which each observed variable is predicted by the other variables in the
#' dataset and with this indicates the suitability for factor analysis.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#'  correlations.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
#'  data is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. Passed to \code{\link[stats:cor]{stats::cor}}.
#' Default is "pearson".
#'
#' @details Kaiser (1970) proposed this index, originally called measure of
#' sampling adequacy (MSA), that indicates how near the inverted correlation
#' matrix \eqn{R^{-1}} is to a diagonal matrix \eqn{S} to determine a given
#' correlation matrix's (\eqn{R}) suitability for factor analysis.
#' The index is
#' \deqn{KMO = \frac{\sum\limits_{i<j}\sum r_{ij}^2}{\sum\limits_{i<j}\sum r_{ij}^2 + \sum\limits_{i<j}\sum q_{ij}^2}}
#' with \eqn{Q = SR^{-1}S} and S = \eqn{(diag R^{-1})^{-1/2}} where
#' \eqn{\sum\limits_{i<j}\sum r_{ij}^2} is the sum of squares of the upper
#' off-diagonal elements of \eqn{R} and \eqn{\sum\limits_{i<j}\sum q_{ij}^2} is the
#' sum of squares of the upper off-diagonal elements of \eqn{Q} (see also Cureton & D'Augustino, 1983).
#'
#' So KMO varies between 0 and 1, with larger values indicating higher suitability
#' for factor analysis. Kaiser and Rice (1974) suggest that KMO should at least
#' exceed .50 for a correlation matrix to be suitable for factor analysis.
#'
#' This function was heavily influenced by the \code{\link[psych:KMO]{psych::KMO}}
#' function.
#'
#' See also \code{\link{BARTLETT}} for another test of suitability for factor
#' analysis.
#'
#' The \code{KMO} function can also be called together with the
#' \code{\link{BARTLETT}} function and with factor retention criteria in the
#'  \code{\link{N_FACTORS}} function.
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
#' @source Cureton, E. E. & D'Augustino, R. B. (1983). Factor analysis: An
#'  applied approach. Hillsdale, N.J.: Lawrence Erlbaum Associates, Inc.
#'
#' @seealso \code{\link{BARTLETT}} for another measure to determine
#' suitability for factor analysis.
#'
#' \code{\link{N_FACTORS}} as a wrapper function for this function,
#' \code{\link{BARTLETT}} and several factor retention criteria.
#'
#' @examples
#' KMO(test_models$baseline$cormat)
KMO <- function(x, use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                           "everything", "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall")) {

  # Perform argument checks
  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n"))

  }

  use <- match.arg(use)
  cor_method <- match.arg(cor_method)

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    R <- x

  } else {

    message(cli::col_cyan(cli::symbol$info, " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n"))

    R <- stats::cor(x, use = use, method = cor_method)
    colnames(R) <- colnames(x)
    N <- nrow(x)

  }

  # Check if correlation matrix is invertable, if it is not, stop with message
  R_i <- try(solve(R), silent = TRUE)

  if (inherits(R_i, "try-error")) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Correlation matrix is singular, no further analyses are performed.\n"))
  }

  # Check if correlation matrix is positive definite
  if(any(eigen(R, symmetric = TRUE, only.values = TRUE)$values <= 0)){

    R <- psych::cor.smooth(R)
    R_i <- try(solve(R), silent = TRUE)

  }

  # Start computations
  S <- diag(diag(R_i)^(-1/2))
  Q <- S %*% R_i %*% S
  diag(Q) <- 0
  diag(R) <- 0
  sumQ2 <- sum(Q^2)
  sumR2 <- sum(R^2)
  KMO <- sumR2/(sumR2 + sumQ2)
  KMO_i <- colSums(R^2)/(colSums(R^2) + colSums(Q^2))

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
