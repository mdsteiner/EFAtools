#' Print function for EKC objects
#'
#' @param x a list of class EKC. Output from \code{\link{EKC}} function.
#' @param plot logical. Whether to plot the results.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print EKC
#'
#' @examples
#' EKC_base <- EKC(test_models$baseline$cormat, N = 500)
#' EKC_base
#'
print.EKC <- function(x, plot = TRUE, ...) {

  n_factors_BvA2017 <- x$n_factors_BvA2017
  n_factors_AM2019 <- x$n_factors_AM2019
  ekc_type <- x$settings$type

  cat(cli::rule("Number of factors suggested by EKC",
                col = "blue"))
  cat("\n")
  cat("\n")

  if("BvA2017" %in% ekc_type){

    cat(crayon::blue(cli::symbol$bullet, "Original implementation (Braeken & van Assen, 2017): "),
        crayon::bold(n_factors_BvA2017))
    cat("\n")

  }

  if("AM2019" %in% ekc_type){

    cat(crayon::blue(cli::symbol$bullet,
                     "Adapted implementation (Auerswald & Moshagen, 2019): "),
        crayon::bold(n_factors_AM2019))
    cat("\n")

  }

  cli::cli_alert_info("Different implementations of EKC exist. Make sure to report which one you relied on (see EKC help page for details)")


  if (isTRUE(plot)) {
    graphics::plot(x)
  }

}
