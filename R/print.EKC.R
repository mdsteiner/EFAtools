#' Print function for EKC objects
#'
#' @param x a list of class EKC. Output from \link{EKC} function.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print EKC
#'
#' @examples
#' EKC_base <- EKC(test_models$baseline$cormat, N = 500)
#' EKC_base
#'
print.EKC <- function(x, ...) {

  nfac <- x$n_factors

  cat("\n")

  if(nfac == 1){
    cat("Empirical Kaiser criterion suggests", crayon::bold(nfac), "factor")
    cat("\n")
    cat("\n")
  } else {
    cat("Empirical Kaiser criterion suggests", crayon::bold(nfac), "factors")
    cat("\n")
    cat("\n")
  }

  graphics::plot(x)

}
