#' Empirical Kaiser Criterion
#'
#' The empirical Kaiser criterion incorporates random sampling variations of the
#' eigenvalues from the Kaiser-Guttman criterion (\code{\link{KGC}}; see Auerswald & Moshagen
#' , 2019; Braeken & van Assen, 2017). The code is based on Auerswald and Moshagen
#' (2019).
#'
#' @param x data.frame or matrix. data.frame or matrix of raw data or matrix with
#' correlations.
#' @param N numeric. The number of observations. Only needed if x is a correlation
#'  matrix.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
#'  data is given as input. Default is  \code{"pairwise.complete.obs"}.
#' @param cor_method character. Passed to \code{\link[stats:cor]{stats::cor}}.
#'  Default is  \code{"pearson"}.
#'
#' @details The Kaiser-Guttman criterion was defined with the intend that a factor
#'  should only be extracted if it explains at least as much variance as a single
#'  factor (see \code{\link{KGC}}). However, this only applies to population-level
#'  correlation matrices. Due to sampling variation, the KGC strongly overestimates
#'  the number of factors to retrieve (e.g., Zwick & Velicer, 1986). To account
#'  for this and to introduce a factor retention method that performs well with
#'  small number of indicators and correlated factors (cases where the performance
#'  of parallel analysis, see \code{\link{PARALLEL}}, is known to deteriorate)
#'  Braeken and van Assen (2017) introduced the empirical Kaiser criterion in
#'  which a series of reference eigenvalues is created as a function of the
#'  variables-to-sample-size ratio and the observed eigenvalues.
#'
#'  Braeken and van Assen (2017) showed that "(a) EKC performs about as well as
#'  parallel analysis for data arising from the null, 1-factor, or orthogonal
#'  factors model; and (b) clearly outperforms parallel analysis for the specific
#'  case of oblique factors, particularly whenever factor intercorrelation is
#'  moderate to high and the number of variables per factor is small, which is
#'  characteristic of many applications these days" (p.463-464).
#'
#'  The \code{EKC} function can also be called together with other factor
#'   retention criteria in the \code{\link{N_FACTORS}} function.
#'
#' @return A list of class EKC containing
#'
#' \item{eigenvalues}{A vector containing the eigenvalues found on the correlation matrix of the entered data.}
#' \item{n_factors}{The number of factors to retain according to the empirical Kaiser criterion.}
#' \item{references}{The reference eigenvalues.}
#' \item{settings}{A list with the settings used.}
#'
#' @source Auerswald, M., & Moshagen, M. (2019). How to determine the number of
#' factors to retain in exploratory factor analysis: A comparison of extraction
#' methods under realistic conditions. Psychological Methods, 24(4), 468–491.
#' https://doi.org/10.1037/met0000200
#'
#' @source Braeken, J., & van Assen, M. A. (2017). An empirical Kaiser criterion.
#' Psychological Methods, 22, 450 – 466. http://dx.doi.org/10.1037/ met0000074
#'
#' @source Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
#' determining the number of components to retain. Psychological Bulletin, 99,
#' 432–442. http://dx.doi.org/10.1037/0033-2909.99.3.432
#'
#' @seealso Other factor retention criteria: \code{\link{CD}},
#'  \code{\link{HULL}}, \code{\link{KGC}}, \code{\link{PARALLEL}},
#'  \code{\link{SMT}}
#'
#'   \code{\link{N_FACTORS}} as a wrapper function for this and all
#'   the above-mentioned factor retention criteria.
#' @export
#'
#' @examples
#' EKC(test_models$baseline$cormat, N = 500)
EKC <- function(x, N = NA,
                use = c("pairwise.complete.obs", "all.obs",
                           "complete.obs", "everything",
                           "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall")) {

  # Perform argument checks
  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n"))

  }

  checkmate::assert_count(N, na.ok = TRUE)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    R <- x

    if (is.na(N)) {

      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Argument 'N' was NA but correlation matrix was entered. Please either provide N or raw data.\n"))

    }

  } else {

    message(cli::col_cyan(cli::symbol$info, " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n"))

    if (!is.na(N)) {
      warning(crayon::yellow$bold("!"), crayon::yellow(" 'N' was set and data entered. Taking N from data.\n"))
    }

    R <- stats::cor(x, use = use, method = cor_method)
    colnames(R) <- colnames(x)
    N <- nrow(x)

  }

  # Check if correlation matrix is invertable, if it is not, stop with message
  R_i <- try(solve(R), silent = TRUE)

  if (inherits(R_i, "try-error")) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Correlation matrix is singular, no further analyses are performed\n"))
  }

  # Check if correlation matrix is positive definite, if it is not,
  # smooth the matrix (cor.smooth throws a warning)
  if (any(eigen(R, symmetric = TRUE, only.values = TRUE)$values <= 0)) {

    R <- psych::cor.smooth(R)

  }

  p <- ncol(R)

  # eigenvalues
  lambda <- eigen(R, symmetric = TRUE, only.values = TRUE)$values

  # reference values
  refs <- vector("double", p)
  for (i in seq_len(p)) {
    refs[i] <- max( ((1 + sqrt(p / N))^2) * (p - sum(refs))/
                        (p - i + 1), 1)

  }

  out <- list(
    eigenvalues = lambda,
    n_factors = which(lambda <= refs)[1] - 1,
    references = refs,
    settings = list(
      use = use,
      cor_method = cor_method,
      N = N
    )
  )

  class(out) <- "EKC"

  return(out)

}
