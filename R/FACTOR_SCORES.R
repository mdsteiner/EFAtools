#' Estimate factor scores for an EFA model
#'
#' This is a wrapper function for
#' [psych::factor.scores()] to be used directly
#' with an output from [EFA()] or by manually specifying the factor
#' loadings and intercorrelations. Calculates factor scores according to the
#' specified methods if raw data are provided, and only factor weights if a
#' correlation matrix is provided.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data (needed to get
#' factor scores) or matrix with correlations.
#' @param f object of class [EFA()] or matrix.
#' @param Phi matrix. A matrix of factor intercorrelations. Only needs to be
#' specified if a factor loadings matrix is entered directly into `f`.
#' Default is `NULL`, in which case all intercorrelations are assumed to be zero.
#' @param rho matrix. Used when `x` is a matrix of raw data and the user
#' wishes to get factor scores for a correlation matrix other than Pearson's
#' (e.g. polychoric). Defaults to `NULL`, in which case
#' [psych::factor.scores()] uses
#' `cor(x, use = "pairwise")`.
#' @param method character. The method used to calculate factor scores. One of
#' "Thurstone" (regression-based; default), "tenBerge", "Anderson", "Bartlett",
#' "Harman", or "components".
#' See [psych::factor.scores()] for details.
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
#' # Example with raw data with method "Bartlett"
#' EFA_raw <- EFA(DOSPERT_raw, n_factors = 10, type = "EFAtools", method = "PAF",
#'                rotation = "oblimin", randomStarts = 0)
#' fac_scores_raw <- FACTOR_SCORES(DOSPERT_raw, f = EFA_raw, method = "Bartlett")
#'
#' # Same as above, but with raw data AND a correlation matrix
#' cor_pearson <- cor(DOSPERT_raw)
#' EFA_cor_pearson <- EFA(cor_pearson, n_factors = 10, N = nrow(DOSPERT_raw),
#'                        type = "EFAtools", method = "PAF",
#'                        rotation = "oblimin", randomStarts = 0)
#' fac_scores_cor_pearson <- FACTOR_SCORES(DOSPERT_raw, f = EFA_cor_pearson,
#'                                         rho = cor_pearson,
#'                                         method = "Bartlett")
#'
#' # Scores between two alternatives above are identical
#' isTRUE(all.equal(fac_scores_raw$scores, fac_scores_cor_pearson$scores,
#'                  check.attributes = FALSE))
#'
#' # Example with a correlation matrix only (does not return factor scores)
#' EFA_cor <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                type = "EFAtools", method = "PAF", rotation = "oblimin")
#' fac_scores_cor <- FACTOR_SCORES(test_models$baseline$cormat, f = EFA_cor)
#'
FACTOR_SCORES <- function(x, f, Phi = NULL, rho = NULL,
                          method = c("Thurstone", "tenBerge", "Anderson",
                                     "Bartlett", "Harman", "components")){

  .assert_cor_input(x)

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    cli::cli_inform(
      c("i" = "{.arg x} is a correlation matrix; factor scores cannot be computed. Enter raw data to get factor scores."),
      class = "efa_scores_needs_raw"
    )

  }

method <- match.arg(method)
checkmate::assert_matrix(Phi, null.ok = TRUE)

if(!inherits(f, c("EFA", "matrix", "LOADINGS"))){

  cli::cli_abort("{.arg f} must be an {.cls EFA} object, a matrix, or a {.cls LOADINGS} object.",
                 class = "efa_scores_bad_f")

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

    cli::cli_inform(
      c("i" = "{.arg Phi} was left {.code NULL} and loadings were entered directly in {.arg f}; assuming uncorrelated factors."),
      class = "efa_scores_phi_null"
    )

  }

}

out_fac_scores <- psych::factor.scores(x = x, f = f, Phi = Phi, method = method,
                                       rho = rho)

settings <- list(method = method)

output <- c(out_fac_scores,
            settings = list(settings))

class(output) <- "FACTOR_SCORES"

return(output)

}
