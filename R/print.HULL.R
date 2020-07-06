#' Print function for HULL objects
#'
#' @param x a list of class HULL. Output from the \code{\link{HULL}} function.
#' @param plot logical. Whether to plot the results.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print HULL
#'
#' @examples
#' \donttest{
#' HULL(test_models$baseline$cormat, N = 500, method = "ML")
#' }
print.HULL <- function(x, plot = TRUE, ...) {

  method <- x$settings$method
  gof <- x$settings$gof

  cat("Hull Analysis performed testing ", crayon::bold("0"), " to ",
      crayon::bold(x$n_fac_max), " factors.", sep = "")

  if(length(gof) == 1){
  cat("\n")
  cat(crayon::bold(method), " estimation and the ",
      .settings_string(gof), " fit index was used.", sep = "")
  cat("\n")
  } else {
    cat("\n")
    cat(crayon::bold(method), " estimation and the ",
        .settings_string(gof), " fit indices were used.", sep = "")
    cat("\n")
  }

  cat("\n")
  cat(cli::rule(left = "Number of factors suggested by the Hull method",
                col = "blue"))
  cat("\n")

  if("CAF" %in% gof){
    cat("\n")
    cat(crayon::blue(cli::symbol$circle_dotted, "With", crayon::bold("CAF:   ")),
        crayon::bold(x$n_fac_CAF))
    cat("\n")
  }

  if("CFI" %in% gof){
    cat(crayon::blue(cli::symbol$circle_dotted, "With", crayon::bold("CFI:   ")),
        crayon::bold(x$n_fac_CFI))
    cat("\n")
  }

  if("RMSEA" %in% gof){
    cat(crayon::blue(cli::symbol$circle_dotted, "With", crayon::bold("RMSEA: ")),
        crayon::bold(x$n_fac_RMSEA))
    cat("\n")
  }

  if (isTRUE(plot)) {
    graphics::plot(x)
  }

}
