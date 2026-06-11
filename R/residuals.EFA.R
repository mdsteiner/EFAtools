#' Extract residuals from an EFA object
#'
#' Returns the residual correlation matrix of an [EFA()] or [EFA_POOLED()]
#' solution. Residuals are a pure extractor here; their diagnostics and a
#' formatted display are part of [summary()] of the EFA object.
#'
#' @param object a list of class EFA. Output from [EFA()] or [EFA_POOLED()].
#' @param type character. Which residuals to return. `"raw"` (default) returns
#'   `orig_R - model_implied_R`; `"standardized"` returns the standardized
#'   residuals (residuals divided by their standard errors), available only when
#'   the object was fitted with bootstrap standard errors.
#' @param ... Further arguments (currently unused).
#'
#' @returns A numeric matrix of residual correlations.
#'
#' @export
#' @method residuals EFA
#'
#' @examples
#' efa <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500)
#' residuals(efa)
residuals.EFA <- function(object, type = c("raw", "standardized"), ...) {
  type <- match.arg(type)

  if (identical(type, "standardized")) {
    if (is.null(object$standardized_residuals)) {
      cli::cli_abort(
        c(
          "Standardized residuals are not available for this object.",
          "i" = "They require bootstrap standard errors ({.code se = \"np-boot\"})."
        ),
        class = "efa_no_standardized_residuals"
      )
    }
    return(object$standardized_residuals)
  }

  object$residuals
}
