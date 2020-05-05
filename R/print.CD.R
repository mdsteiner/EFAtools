#' Print function for CD objects
#'
#' @param x a list of class CD. Output from \link{CD} function.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print CD
#'
print.CD <- function(x, ...) {

  nfac <- x$n_factors

  cat("\n")

  if(nfac == 1){
    cat("Comparison Data analysis suggests", crayon::bold(nfac), "factor")
    cat("\n")
    cat("\n")
  } else {
    cat("Comparison Data analysis suggests", crayon::bold(nfac), "factors")
    cat("\n")
    cat("\n")
  }

  graphics::plot(x)

}
