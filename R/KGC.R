#' Kaiser-Guttman Criterion
#'
#' Probably the most popular factor retention criterion. Kaiser and Guttman suggested
#' to retain as many factors as there are sample eigenvalues greater than 1.
#' This is why the criterion is also known as eigenvalues-greater-than-one rule.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations.
#' @param eigen_type character. On what the eigenvalues should be found. Can be
#'  either "PCA", "SMC", or "EFA", or some combination of them. If using "PCA",
#'  the diagonal values of the correlation matrices are left to be 1. If using
#'  "SMC", the diagonal of the
#'  correlation matrices is replaced by the squared multiple correlations (SMCs)
#'  of the indicators. If using "EFA", eigenvalues are found on the correlation
#'  matrices with the final communalities of an exploratory factor analysis
#'  solution (default is principal axis factoring extracting 1 factor) as
#'  diagonal.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
#'  data is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. Passed to \code{\link[stats:cor]{stats::cor}}.
#' Default is "pearson".
#' @param n_factors numeric. Number of factors to extract if "EFA" is included in
#' \code{eigen_type}. Default is 1.
#' @param ... Additional arguments passed to \code{\link{EFA}}. For example,
#' to change the extraction method (PAF is default).
#'
#' @details Originally, the Kaiser-Guttman criterion was intended for the use
#' with prinicpal components, hence with eigenvalues derived from the original
#' correlation matrix. This can be done here by setting \code{eigen_type} to
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
#' example parallel analysis (\code{\link{PARALLEL}}), the Hull method
#' \code{\link{HULL}}, or sequential \eqn{chi^2} model tests (\code{\link{SMT}};
#' see Auerswald & Moshagen, 2019).
#'
#' The \code{KGC} function can also be called together with other factor
#' retention criteria in the \code{\link{N_FACTORS}} function.
#'
#' @return A list of class KGC containing
#'
#' \item{eigen_PCA}{ A vector containing the eigenvalues found with PCA.}
#' \item{eigen_SMC}{ A vector containing the eigenvalues found with SMCs.}
#' \item{eigen_EFA}{ A vector containing the eigenvalues found with EFA.}
#' \item{n_fac_PCA}{ The number of factors to retain according to the Kaiser-
#' Guttmann criterion with PCA eigenvalues type.}
#' \item{n_fac_SMC}{ The number of factors to retain according to the Kaiser-
#' Guttmann criterion with SMC eigenvalues type.}
#' \item{n_fac_EFA}{ The number of factors to retain according to the Kaiser-
#' Guttmann criterion with EFA eigenvalues type.}
#' \item{settings}{A list of the settings used.}
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
#' @seealso Other factor retention criteria: \code{\link{CD}}, \code{\link{EKC}},
#' \code{\link{HULL}}, \code{\link{PARALLEL}}, \code{\link{SMT}}
#'
#' \code{\link{N_FACTORS}} as a wrapper function for this and all the
#' above-mentioned factor retention criteria.
#'
#' @export
#'
#' @examples
#' KGC(test_models$baseline$cormat, eigen_type = c("PCA", "SMC"))
KGC <- function(x, eigen_type = c("PCA", "SMC", "EFA"),
                use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                        "everything", "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall"), n_factors = 1,
                ...){

  # Perform argument checks
  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n"))

  }

  eigen_type <- match.arg(eigen_type, several.ok = TRUE)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  checkmate::assert_count(n_factors)

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    R <- x

  } else {

    message(cli::col_cyan(cli::symbol$info, " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n"))

    R <- stats::cor(x, use = use, method = cor_method)
    colnames(R) <- colnames(x)
  }

  # Check if correlation matrix is invertable, if it is not, stop with message
  R_i <- try(solve(R), silent = TRUE)

  if (inherits(R_i, "try-error")) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Correlation matrix is singular, no further analyses are performed.\n"))

  }

  # Check if correlation matrix is positive definite, if it is not,
  # smooth the matrix (cor.smooth throws a warning)
  if (any(eigen(R, symmetric = TRUE, only.values = TRUE)$values <= 0)) {

    R <- psych::cor.smooth(R)

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
    EFA_h2 <- suppressMessages(suppressWarnings(EFA(R, n_factors = n_factors, ...)$h2))
    diag(R) <- EFA_h2

    # Calculate eigenvalues and determine number of factors
    eigen_R_EFA <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    n_fac_EFA <- sum(eigen_R_EFA >= 1)

  }

  # prepare settings
  settings <- list(eigen_type = eigen_type,
                   use = use,
                   cor_method = cor_method,
                   n_factors = n_factors)

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
