#' Estimate factor scores for an EFA model
#'
#' This is a wrapper function for
#' \code{\link[psych:factor.scores]{psych::factor.scores}} to be used directly
#' with an output from \code{\link{EFA}} or by manually specifying the factor
#' loadings and intercorrelations. Calculates factor scores according to the
#' specified methods if raw data are provided, and only factor weights if a
#' correlation matrix is provided.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data (needed to get
#' factor scores) or matrix with correlations.
#' @param f object of class \code{\link{EFA}} or matrix.
#' @param Phi matrix. A matrix of factor intercorrelations. Only needs to be
#' specified if a factor loadings matrix is entered directly into \code{f}.
#' Default is \code{NULL}, in which case all intercorrelations are assumed to be zero.
#' @param method character. The method used to calculate factor scores. One of
#' "Thurstone" (regression-based; default), "tenBerge", "Anderson", "Bartlett",
#' "Harman", or "components".
#' See \code{\link[psych:factor.scores]{psych::factor.scores}} for details.
#' @param impute character. Whether and how missing values in \code{x} should
#' be imputed. One of "none" (default, only complete cases are scored), "median",
#' or "mean".
#'
#' @return A list of class FACTOR_SCORES containing the following:
#'
#' \item{scores}{The factor scores (only if raw data are provided.)}
#' \item{weights}{The factor weights.}
#' \item{r.scores}{The correlations of the factor score estimates.}
#' \item{missing}{A vector of the number of missing observations per subject
#' (only if raw data are provided.}
#' \item{R2}{Multiple R2 of the scores with the factors.}
#' \item{settings}{A list of the settings used.}
#'
#' @export
#'
#' @examples
#' # Example with raw data with method "Bartlett" and no imputation
#' EFA_raw <- EFA(DOSPERT_raw, n_factors = 10, type = "EFAtools", method = "PAF",
#'                rotation = "oblimin")
#' fac_scores_raw <- FACTOR_SCORES(DOSPERT_raw, f = EFA_raw, method = "Bartlett",
#'                                 impute = "none")
#'
#' # Example with a correlation matrix (does not return factor scores)
#' EFA_cor <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                type = "EFAtools", method = "PAF", rotation = "oblimin")
#' fac_scores_cor <- FACTOR_SCORES(test_models$baseline$cormat, f = EFA_cor)
#'
FACTOR_SCORES <- function(x, f, Phi = NULL,
                          method = c("Thurstone", "tenBerge", "Anderson",
                                     "Bartlett", "Harman", "components"),
                          impute = c("none", "means", "median")){

  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n"))

  }

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    message(cli::col_cyan(cli::symbol$info, " 'x' is a correlation matrix, factor scores cannot be computed. Enter raw data to get factor scores.\n"))

  }

method <- match.arg(method)
impute <- match.arg(impute)
checkmate::assert_matrix(Phi, null.ok = TRUE)

if(!inherits(f, c("EFA", "matrix", "LOADINGS"))){

  stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'f' is neither an object of class EFA nor a matrix, nor of class LOADINGS.\n"))

}

if(inherits(f, c("EFA"))){

  Phi <- f$Phi

  if(f$settings$rotation != "none"){
    f <- unclass(f$rot_loadings)
  } else {
    f <- unclass(f$unrot_loadings)
  }

} else {

  f <- unclass(f)

  if(is.null(Phi)){

    message(cli::col_cyan(cli::symbol$info, " Phi argument was left NULL and factor loadings were entered directly in f. Assuming uncorrelated factors.\n"))

  }

}

out_fac_scores <- psych::factor.scores(x = x, f = f, Phi = Phi, method = method,
                                       rho = NULL, impute = impute)

settings <- list(method = method,
                 impute = impute)

output <- c(out_fac_scores,
            settings = list(settings))

class(output) <- "FACTOR_SCORES"

return(output)

}
