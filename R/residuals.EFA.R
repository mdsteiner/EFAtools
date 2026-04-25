#' Residuals function for EFA objects
#'
#' @param object a list of class EFA. Output from \code{\link{EFA}} or
#'   \code{\link{EFA_POOLED}}.
#' @param print logical. Whether to print residuals. If FALSE, they are just returned.
#' @param digits numeric. The number of digits to round printed residuals.
#' @param ... Further arguments.
#'
#' @export
#' @method residuals EFA
residuals.EFA <- function(object, print = TRUE, digits = 3, ...) {

  out <- list(
    residuals = object$residuals
  )

  if ("standardized_residuals" %in% names(object)) {
    out$standardized_residuals <- object$standardized_residuals
  }

  if (isFALSE(print)) {

    return(out)
  }

  cli::cli_h2("Raw residuals, based on {.code orig_R - model_implied_R}")

  print(round(object$residuals, digits = digits))



  if ("standardized_residuals" %in% names(object)) {

    cli::cli_h2("Standardized residuals, based on {.code residuals / SE}")

    print(round(object$standardized_residuals, digits = digits))

  }

  invisible(out)
}
