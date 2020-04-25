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

  cat("\n")
  cat("Eigenvalues were found using", crayon::bold(x$settings$eigen_type))
  cat("\n")
  cat("\n")

  cat("Kaiser-Guttmann criterion suggests", crayon::bold(x$n_factors), "factors")
  cat("\n")

  graphics::plot(x)

}
