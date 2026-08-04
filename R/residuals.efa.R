#' Extract residuals from an efa object
#'
#' Returns the residual correlation matrix of an [efa_fit()] or [efa_mi()]
#' solution. Residuals are a pure extractor here; their diagnostics and a
#' formatted display are part of [summary()] of the efa object.
#'
#' @param object a list of class `efa`. Output from [efa_fit()] or [efa_mi()].
#' @param type character. Which residuals to return. `"raw"` (default) returns
#'   `orig_R - model_implied_R`; `"standardized"` returns the standardized
#'   residuals (residuals divided by their standard errors), available only when
#'   the object was fitted with bootstrap standard errors.
#' @param ... Further arguments (currently unused).
#'
#' @returns A numeric matrix of residual correlations.
#'
#' @export
#' @method residuals efa
#'
#' @examples
#' efa <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500)
#' residuals(efa)
residuals.efa <- function(object, type = c("raw", "standardized"), ...) {
  type <- .match_arg_ci(type)

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
