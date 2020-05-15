#' Print function for HULL objects
#'
#' @param x a list of class HULL Output from the \link{HULL} function.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print HULL
#'
#' @examples
#' \dontrun{
#' HULL(IDS2_R, N = 2000)
#' }
print.HULL <- function(x, ...) {

  cat(paste0("Hull Analysis performed, testing ", crayon::bold("0"), " to ",
             crayon::bold(max(x$solutions[, 1])),
             " factors."))

  cat("\n")
  cat(paste0(crayon::bold(max(x$ctrl$method)), " estimation and the ",
             crayon::bold(max(x$ctrl$gof)), " fit index was used."))
  cat("\n")
  cat("\n")
  cat("Number of factors suggested by the Hull method:")
  cat("\n")

  cat(paste("    ",cli::symbol$bullet, crayon::bold(x$n_factors)))

  graphics::plot(x)

}
