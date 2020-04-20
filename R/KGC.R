#' Kaiser-Guttman Criterion
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations.
#' @param eigen_type character. On what the eigenvalues should be found. Can be
#'  either "PCA", "SMC", or "EFA". If using "PCA", the diagonal values of the
#'  correlation matrices are left to be 1. If using "SMC", the diagonal of the
#'  correlation matrices is replaced by the squared multiple correlations (SMCs)
#'  of the indicators. If using "EFA", eigenvalues are found on the correlation
#'  matrices with the final communalities of an exploratory factor analysis
#'  solution (default is principal axis factoring extracting 1 factor) as
#'  diagonal.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
#'  data is given as input. Default is "pairwise.complete.obs".
#' @param n_factors numeric. Number of factors to extract if
#' \code{eigen_type = "EFA"}. Default is 1.
#' @param ... Additional arguments passed to \code{\link[EFA]{EFA}}.
#'
#' @details DETAILS
#'
#'
#' @return A list of class KGC containing
#'
#' \item{eigenvalues}{A vector containing the eigenvalues.}
#' \item{n_factors}{The number of factors to retain according to the Kaiser-
#' Guttmann criterion.}
#'
#' @source SOURCE
#'
#' @export
#'
#' @examples
#' KGC(IDS2_R, eigen_type = "PCA")
#'
KGC <- function(x, eigen_type = c("PCA", "SMC", "EFA"),
                use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                        "everything", "na.or.complete"), n_factors = 1, ...){

  eigen_type <- match.arg(eigen_type)
  use <- match.arg(use)

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    if(any(is.na(x))){

      stop("The correlation matrix you entered contains missing values. Eigenvalues
           cannot be computed.")

    }

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
    stop("Correlation matrix is singular, factor analysis is not possible")
  }

  if(eigen_type == "SMC") {

    # Calculate SMCs and replace diagonal of correlation matrix with these
    inverse_R <- solve(R)
    SMCs <- 1 - 1 / diag(inverse_R)
    diag(R) <- SMCs

  } else if(eigen_type == "EFA") {

    # Do an EFA to get final communality estimates and replace diagonal of
    # correlation matrix with these

    EFA_h2 <- EFA(R, n_factors = n_factors, ...)$h2
    diag(R) <- EFA_h2

  }

  # Calculate eigenvalues and determine number of factors
  eigen_R <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  n_fac <- sum(eigen_R >= 1)

  # prepare settings
  settings <- list(eigen_type = eigen_type)

  # Prepare the output
  output <- list(
    eigenvalues = eigen_R,
    n_factors = n_fac,
    settings = settings
  )

  class(output) <- "KGC"

  return(output)

}
