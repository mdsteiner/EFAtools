#' Kaiser-Guttman Criterion
#'
#' Probably the most popular factor retention criterion. Kaiser and Guttman suggested
#' to retain as many factors as there are sample eigenvalues greater than 1.
#' This is why the criterion is also known as eigenvalues-greater-than-one rule.
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
#' @details Originally, the Kaiser-Guttman criterion was intended for the use
#' with prinicpal components, hence with eigenvalues derived from the original
#' correlation matrix. This can be done here be setting \code{eigen_type} to
#' "PCA". However, it is well-known that this criterion is often inaccurate and
#' that it tends to overestimate the number of factors, especially for unidimensional
#' or orthogonal factor structures (e.g., Zwick & Velicer, 1986).
#'
#' The criterion's inaccuracy in these cases is somewhat addressed if it is
#' applied on the correlation matrix with communalities in the diagonal, either
#' initial communalities estimated from SMCs (done setting \code{eigen_type} to
#' "SMC") or final communality estimates from an EFA (done setting \code{eigen_type}
#' to "EFA"; see Auerswald & Moshagen, 2019). However, although this variant
#' of the KGC is more accurate in some cases compared to the traditional KGC, it
#' is at the same time less accurate than the PCA-variant in other cases, and it
#' is still often less accurate than other factor retention methods, for
#' example parallel analysis (\link{\code{PARALLEL}}), the Hull method
#' \link{\code{HULL}}, or sequential \eqn{chi^2} model tests (\link{\code{SMT}};
#' see Auerswald & Moshagen, 2019).
#'
#' @return A list of class KGC containing
#'
#' \item{eigen_PCA}{A vector containing the eigenvalues found with PCA.}
#' \item{eigen_SMC}{A vector containing the eigenvalues found with SMCs.}
#' \item{eigen_EFA}{A vector containing the eigenvalues found with EFA.}
#' \item{n_fac_PCA}{The number of factors to retain according to the Kaiser-
#' Guttmann criterion with PCA eigenvalues type.}
#' \item{n_fac_SMC}{The number of factors to retain according to the Kaiser-
#' Guttmann criterion with SMC eigenvalues type.}
#' \item{n_fac_EFA}{The number of factors to retain according to the Kaiser-
#' Guttmann criterion with EFA eigenvalues type.}
#'
#' @source Auerswald, M., & Moshagen, M. (2019). How to determine the number of
#' factors to retain in exploratory factor analysis: A comparison of extraction
#' methods under realistic conditions. Psychological Methods, 24(4), 468–491.
#' https://doi.org/10.1037/met0000200
#'
#' @source Guttman, L. (1954). Some necessary conditions for common-factor analysis.
#' Psychometrika, 19, 149 –161. http://dx.doi.org/10.1007/BF02289162
#'
#' @source Kaiser, H. F. (1960). The application of electronic computers to factor
#' analysis. Educational and Psychological Measurement, 20, 141–151.
#' http://dx.doi.org/10.1177/001316446002000116
#'
#' @source Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
#' determining the number of components to retain. Psychological Bulletin, 99,
#' 432–442. http://dx.doi.org/10.1037/0033-2909.99.3.432
#'
#' @export
#'
#' @example
#' KGC(test_models$baseline$cormat, eigen_type = c("PCA", "SMC"))
#'
KGC <- function(x, eigen_type = c("PCA", "SMC", "EFA"),
                use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                        "everything", "na.or.complete"), n_factors = 1, ...){

  # Perform argument checks
  if(!inherits(x, c("matrix", "data.frame"))){

    stop("x is neither a matrix nor a dataframe. Either provide a correlation
    matrix or a dataframe or matrix with raw data.")

  }

  eigen_type <- match.arg(eigen_type, several.ok = TRUE)
  use <- match.arg(use)
  checkmate::assert_count(n_factors)

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

  if (inherits(R_i, "try-error")) {
    stop("Correlation matrix is singular, factor analysis is not possible")
  }

  # Store original R for later
  R_orig <- R

  # Prepare objects
  n_fac_PCA <- NA
  n_fac_SMC <- NA
  n_fac_EFA <- NA
  eigen_R_PCA <- NA
  eigen_R_SMC <- NA
  eigen_R_EFA <- NA

  if("PCA" %in% eigen_type) {

    # Calculate eigenvalues and determine number of factors
    eigen_R_PCA <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    n_fac_PCA <- sum(eigen_R_PCA >= 1)

  }

  if("SMC" %in% eigen_type) {

    # Calculate SMCs and replace diagonal of correlation matrix with these
    inverse_R <- solve(R)
    SMCs <- 1 - 1 / diag(inverse_R)
    diag(R) <- SMCs

    # Calculate eigenvalues and determine number of factors
    eigen_R_SMC <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    n_fac_SMC <- sum(eigen_R_SMC >= 1)

  }

  if("EFA" %in% eigen_type) {

    R <- R_orig

    # Do an EFA to get final communality estimates and replace diagonal of
    # correlation matrix with these
    EFA_h2 <- EFA(R, n_factors = n_factors, ...)$h2
    diag(R) <- EFA_h2

    # Calculate eigenvalues and determine number of factors
    eigen_R_EFA <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    n_fac_EFA <- sum(eigen_R_EFA >= 1)

  }

  # prepare settings
  settings <- list(eigen_type = eigen_type)

  # Prepare the output
  output <- list(
    eigen_PCA = eigen_R_PCA,
    eigen_SMC = eigen_R_SMC,
    eigen_EFA = eigen_R_EFA,
    n_fac_PCA = n_fac_PCA,
    n_fac_SMC = n_fac_SMC,
    n_fac_EFA = n_fac_EFA,
    settings = settings
  )

  class(output) <- "KGC"

  return(output)

}
