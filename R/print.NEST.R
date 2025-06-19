#' Print function for NEST objects
#'
#' @param x a list of class NEST Output from \code{\link{NEST}} function.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print NEST
#'
#' @examples
#' NEST(test_models$baseline$cormat, N = 500)
#'
print.NEST <- function(x, ...) {

  nfac <- x$n_factors

  cat("\n")

  if(nfac == 1){
    cat("NEST suggests", crayon::bold(nfac), "factor")
    cat("\n")
    cat("\n")
  } else {
    cat("NEST suggests", crayon::bold(nfac), "factors")
    cat("\n")
    cat("\n")
  }

}
