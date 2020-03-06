#' Kaiser-Meyer-Olkin criterion
#'
#' This function computes the Kaiser-Meyer-Olkin (KMO) criterion overall and for
#'  each variable in a correlation matrix. The KMO represents the degree to
#' which each observed variable is predicted by the other variables in the
#' dataset and with this the suitability for factor analysis. A KMO > 0.6 is
#' considered adequate.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#'  correlations.
#' @param cors logical. If \code{TRUE} (default) a correlation matrix is expected
#' in x, otherwise raw data is expected.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
#'  is given as input. Note that in this case \code{cors} must be set to
#'  \code{FALSE}. Default is "pairwise.complete.obs".
#'
#' @details Kaiser (1970) proposed this index, originally called measure of
#' sampling adequacy (MSA), that indicates how near the inverted correlation
#' matrix \eqn{R^{-1}} is to a diagonal matrix \eqn{S} to determine a given
#' correlation matrix \eqn{R}'s suitability for factor analysis.
#' The index is
#' \deqn{KMO = \frac{\underset{i<j}{\sum\sum} r_{ij}^2}{\underset{i<j}{\sum\sum} r_{ij}^2 + \underset{i<j}{\sum\sum} q_{ij}^2}
#' with \eqn{Q = SR^{-1}S} and S = \eqn{(diag R^{-1})^{-1/2}} where
#' \eqn{\underset{i<j}{\sum\sum} r_{ij}^2}} is the sum of squares of the upper
#' off-diagonal elements of \eqn{R} and \eqn{\underset{i<j}{\sum\sum} q_{ij}^2} is the
#' sum of squares of the upper off-diagnoal elements of \eqn{Q} (see also Cureton & DeAugustino, 1983).
#'
#' So KMO varies between 0 and 1, with larger values indicating higher suitability
#' for factor analysis. Kaiser and Rice (1974) suggest that KMO should at least
#' exceed .50 for a correlation matrix to be suitable for factor analysis.
#'
#'
#' @return A list containing
#' \item{KMO}{Overall KMO.}
#' \item{KMO_i}{KMO for each variable.}
#'
#' @export
#'
#' @source Kaiser, H. F. (1970). A second generation little jiffy. Psychometrika,
#' 35, 401-415.
#' @source Kaiser, H. F. & Rice, J. (1974). Little jiffy, mark IV. Educational
#' and Psychological Measurement, 34, 111-117.
#' @source Cureton, E. E. & DeAugustino, R. B. (1983). Factor analysis: An applied
#' approach. Hillsdale, N.J.: Lawrence Erlbaum Associates, Inc.
#'
#' @examples
#'
#' KMO(IDS2_R)
KMO <- function(x, cors = TRUE, use = "pairwise.complete.obs") {

  # create R correlation matrix object, if from data, using
  # pairwise binary correlations
  if (isTRUE(cors)) {
    R <- x

    # test whether a real correlation matrix is used
    if (nrow(R) != ncol(R)) {
      stop("Entered data is no correlation matrix but cors = TRUE. Either set ",
           "cors = FALSE if you entered raw data, or enter a correlation matrix.")
    }

  } else {
    R <- stats::cor(x, use = use)
    colnames(R) <- colnames(x)
  }

  # Check if correlation matrix is invertable, if it is not, stop with message
  R_i <- try(solve(R))

  if (class(R_i) == "try-error") {
    stop("Matrix is singular, KMO is not computed")
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

  output <- list(KMO = KMO, KMO_i = KMO_i)

  return(output)

}
