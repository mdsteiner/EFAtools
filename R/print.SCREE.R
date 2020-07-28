#' Print function for SCREE objects
#'
#' @param x a list of class SCREE Output from \link{SCREE} function.
#' @param plot logical. Whether to plot the results.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print SCREE
#'
#' @examples
#' SCREE_base <- SCREE(test_models$baseline$cormat)
#' SCREE_base
#'
print.SCREE <- function(x, plot = TRUE, ...) {

  eigen_type <-x$settings$eigen_type

  cat("\n")
  cat("Eigenvalues were found using ", .settings_string(eigen_type), ".",
      sep = "")
  cat("\n")
  cat("\n")

  cat("\n")

  if (isTRUE(plot)) {
    graphics::plot(x)
  }

}
