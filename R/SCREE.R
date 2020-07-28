#' Scree Plot
#'
#' The scree plot was originally introduced by Cattell (1966) to perform the
#' scree test. In a scree plot, the eigenvalues of the factors / components are
#' plotted against the index of the factors / components, ordered from 1 to N
#' factors components, hence from largest to smallest eigenvalue. According to
#' the scree test, the number of factors / components to retain is the number of
#' factors / components to the left of the "elbow" (where the curve starts to
#' level off) in the scree plot.
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
#' @details As the scree test requires visual examination, the test has been
#' especially criticized for its subjectivity and with this low inter-rater
#' reliability. Moreover, a scree plot can be ambiguous if there are either no
#' clear "elbow" or multiple "elbows", making it difficult to judge just where
#' the eigenvalues do level off. Finally, the scree test has also been found to
#' be less accurate than other factor retention criteria. For all these reasons,
#' the scree test has been recommended against, at least for exclusive use as a
#' factor retention criterion (Zwick & Velicer, 1986)
#'
#' The \code{SCREE} function can also be called together with other factor
#' retention criteria in the \code{\link{N_FACTORS}} function.
#'
#' @return A list of class SCREE containing
#'
#' \item{eigen_PCA}{ A vector containing the eigenvalues found with PCA.}
#' \item{eigen_SMC}{ A vector containing the eigenvalues found with SMCs.}
#' \item{eigen_EFA}{ A vector containing the eigenvalues found with EFA.}
#' \item{settings}{A list of the settings used.}
#'
#' @source Cattell, R. B. (1966). The scree test for the number of factors.
#' Multivariate Behavioral Research, 1(2), 245–276.
#' https://doi.org/10.1207/s15327906mbr0102_10
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
#' SCREE(test_models$baseline$cormat, eigen_type = c("PCA", "SMC"))
SCREE <- function(x, eigen_type = c("PCA", "SMC", "EFA"),
                  use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                        "everything", "na.or.complete"),
                  cor_method = c("pearson", "spearman", "kendall"),
                  n_factors = 1, ...){

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
  eigen_R_PCA <- NA
  eigen_R_SMC <- NA
  eigen_R_EFA <- NA

  if("PCA" %in% eigen_type) {

    # Calculate eigenvalues and determine number of factors
    eigen_R_PCA <- eigen(R, symmetric = TRUE, only.values = TRUE)$values

  }

  if("SMC" %in% eigen_type) {

    # Calculate SMCs and replace diagonal of correlation matrix with these
    inverse_R <- solve(R)
    SMCs <- 1 - 1 / diag(inverse_R)
    diag(R) <- SMCs

    # Calculate eigenvalues and determine number of factors
    eigen_R_SMC <- eigen(R, symmetric = TRUE, only.values = TRUE)$values

  }

  if("EFA" %in% eigen_type) {

    R <- R_orig

    # Do an EFA to get final communality estimates and replace diagonal of
    # correlation matrix with these
    EFA_h2 <- suppressMessages(suppressWarnings(EFA(R, n_factors = n_factors, ...)$h2))
    diag(R) <- EFA_h2

    # Calculate eigenvalues and determine number of factors
    eigen_R_EFA <- eigen(R, symmetric = TRUE, only.values = TRUE)$values

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
    settings = settings
  )

  class(output) <- "SCREE"

  return(output)

}
