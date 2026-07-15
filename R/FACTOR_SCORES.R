#' Estimate factor scores for an EFA model
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `FACTOR_SCORES()` has been superseded by [efa_scores()], which is the
#' recommended interface going forward. It remains available and unchanged so
#' existing code keeps working.
#'
#' A convenience wrapper around [efa_scores()] that returns factor scores and
#' weights in a compact list. Factor scores are calculated according to the
#' specified method if raw data are provided, and only factor weights if a
#' correlation matrix is provided.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data (needed to get
#' factor scores) or matrix with correlations.
#' @param f object of class [efa_fit()] or matrix.
#' @param Phi matrix. A matrix of factor intercorrelations. Only needs to be
#' specified if a factor loadings matrix is entered directly into `f`; for an
#' [efa_fit()] object the intercorrelations are taken from the object. Default is
#' `NULL`, in which case the intercorrelations of a directly supplied loading
#' matrix are assumed to be zero.
#' @param rho matrix. Correlation matrix used to derive the scoring weights.
#' Defaults to `NULL`, in which case the matrix the EFA in `f` was fit on
#' (`f$orig_R`) is used, so the weights stay consistent with the loadings even for
#' a non-Pearson correlation (e.g. polychoric); for a directly supplied loading
#' matrix, `x` itself is used when it is a correlation matrix, otherwise the
#' Pearson correlation of `x`. Pass a matrix here to score against a different
#' correlation.
#' @param method character. The method used to calculate factor scores. One of
#' "Thurstone" (regression-based; default), "tenBerge", "Anderson", "Bartlett",
#' "Harman", or "components".
#'
#' @return A list of class FACTOR_SCORES containing the following:
#'
#' \item{scores}{The factor scores (only if raw data are provided.)}
#' \item{weights}{The factor weights.}
#' \item{r.scores}{The correlations of the factor score estimates.}
#' \item{missing}{Whether the raw data contained missing values (only if raw data
#' are provided).}
#' \item{R2}{The squared factor-score determinacy for each factor: the squared
#' correlation between a factor and its estimated score. For orthogonal factors
#' this equals the squared multiple correlation between the factor and the
#' observed variables; for oblique factors it is the score-specific determinacy.
#' See [efa_scores()] for the underlying score-quality diagnostics.}
#' \item{settings}{A list of the settings used.}
#'
#' @seealso [efa_scores()] for the factor-score weights together with the full
#' set of score-quality diagnostics (determinacy, univocality, and Guttman
#' indeterminacy index) and a print/summary method.
#'
#' @export
#'
#' @examples
#' # Example with raw data with method "Bartlett"
#' EFA_raw <- efa_fit(DOSPERT_raw, n_factors = 10, method = "PAF",
#'                    rotation = "oblimin",
#'                    rotate_control = rotate_control(random_starts = 1))
#' fac_scores_raw <- FACTOR_SCORES(DOSPERT_raw, f = EFA_raw, method = "Bartlett")
#'
#' # Same as above, but with raw data AND a correlation matrix
#' cor_pearson <- cor(DOSPERT_raw)
#' EFA_cor_pearson <- efa_fit(cor_pearson, n_factors = 10, N = nrow(DOSPERT_raw),
#'                            method = "PAF", rotation = "oblimin",
#'                            rotate_control = rotate_control(random_starts = 1))
#' fac_scores_cor_pearson <- FACTOR_SCORES(DOSPERT_raw, f = EFA_cor_pearson,
#'                                         rho = cor_pearson,
#'                                         method = "Bartlett")
#'
#' # Scores between two alternatives above are identical
#' isTRUE(all.equal(fac_scores_raw$scores, fac_scores_cor_pearson$scores,
#'                  check.attributes = FALSE))
#'
#' # Example with a correlation matrix only (does not return factor scores)
#' EFA_cor <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                    method = "PAF", rotation = "oblimin")
#' fac_scores_cor <- FACTOR_SCORES(test_models$baseline$cormat, f = EFA_cor)
#'
FACTOR_SCORES <- function(x, f, Phi = NULL, rho = NULL,
                          method = c("Thurstone", "tenBerge", "Anderson",
                                     "Bartlett", "Harman", "components")){

  method <- match.arg(method)

  # "Thurstone" is the regression method under efa_scores' naming; every other
  # label is shared. Delegate the whole computation (and the classed conditions).
  es <- efa_scores(x, f = f, Phi = Phi, rho = rho,
                   method = if (method == "Thurstone") "regression" else method)

  # Reassemble the legacy FACTOR_SCORES output shape.
  output <- list(
    scores = es$scores,
    weights = es$weights,
    r.scores = es$r.scores
  )

  # Only raw data (with returned scores) carries the `missing` field.
  if (!is.null(es$scores)) {
    output$missing <- anyNA(x)
  }

  # Squared factor-score determinacy, named by factor.
  output$R2 <- stats::setNames(es$determinacy$rho2, rownames(es$determinacy))

  # Preserve the original (user-facing) method label in the settings.
  output$settings <- list(method = method)

  class(output) <- "FACTOR_SCORES"

  output

}
