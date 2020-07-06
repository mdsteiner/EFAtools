#' Print function for CD objects
#'
#' @param x a list of class CD. Output from \link{CD} function.
#' @param plot logical. Whether to plot the results.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print CD
#'
#' @examples
#' \donttest{
#' # determine n factors of the GRiPS
#' CD(GRiPS_raw)
#' }
#'
print.CD <- function(x, plot = TRUE, ...) {

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

  if (isTRUE(plot)) {
    graphics::plot(x)
  }


}
