#' Print function for MAP objects
#'
#' @param x a list of class MAP. Output from \code{\link{MAP}} function.
#' @param plot logical. Whether to plot the results.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print MAP
#'
#' @examples
#' MAP_base <- MAP(test_models$baseline$cormat, N = 500)
#' MAP_base
#'
print.MAP <- function(x, plot = TRUE, ...) {

  n_factors_TR2 <- x$n_factors_TR2
  n_factors_TR4 <- x$n_factors_TR4

  cat(cli::rule("Number of factors suggested by MAP",
                col = "blue"))
  cat("\n")
  cat("\n")


    cat(crayon::blue(cli::symbol$bullet, "Original implementation (TR2): "),
        crayon::bold(n_factors_TR2))
    cat("\n")

    cat(crayon::blue(cli::symbol$bullet,
                     "Revised implementation (TR4): "),
        crayon::bold(n_factors_TR4))
    cat("\n")


}
