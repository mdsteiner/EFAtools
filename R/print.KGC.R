#' Print function for KGC objects
#'
#' @param x a list of class KGC. Output from \link{KGC} function.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print KGC
#'
#' @examples
#' KGC_base <- KGC(test_models$baseline$cormat, eigen_type = "PCA")
#' KGC_base
#'
print.KGC <- function(x, ...) {

  nfac <- x$n_factors

  cat("\n")
  cat("Eigenvalues were found using", crayon::bold(x$settings$eigen_type))
  cat("\n")
  cat("\n")

  if(nfac == 1){
  cat("Kaiser-Guttmann criterion suggests", crayon::bold(nfac), "factor")
  cat("\n")
  } else {
    cat("Kaiser-Guttmann criterion suggests", crayon::bold(nfac), "factors")
    cat("\n")
  }

  graphics::plot(x)

}
