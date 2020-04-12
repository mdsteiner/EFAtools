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
#' sum of squares of the upper off-diagonal elements of \eqn{Q} (see also Cureton & DeAugustino, 1983).
#'
#' So KMO varies between 0 and 1, with larger values indicating higher suitability
#' for factor analysis. Kaiser and Rice (1974) suggest that KMO should at least
#' exceed .50 for a correlation matrix to be suitable for factor analysis.
#'
#' This function was heavily influenced by the \code{\link[psych:KMO]{KMO}}
#' function from the psych package.
#'
#' See also \code{\link[EFAtools:BARTLETT]{EFAtools:BARTLETT}} for another test
#' of suitability for factor analysis.
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
#' @source Cureton, E. E. & DeAugustino, R. B. (1983). Factor analysis: An
#'  applied approach. Hillsdale, N.J.: Lawrence Erlbaum Associates, Inc.
#'
#' @examples
#' KMO(IDS2_R)
#'
KMO <- function(x, use = c("all.obs", "complete.obs", "pairwise.complete.obs",
                           "everything", "na.or.complete")) {

  use <- match.arg(use)

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    R <- x

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

  if(!is.null(colnames(R))){

    names(KMO_i) <- colnames(R)

  } else if(!is.null(rownames(R))) {

    names(KMO_i) <- rownames(R)

  }

  output <- list(KMO = KMO, KMO_i = KMO_i)

  class(output) <- "KMO"

  return(output)

}
