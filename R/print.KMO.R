#' Print KMO object
#'
#' @param x list of class KMO (output from the \link{KMO} function)
#' @param ... additional arguments passed to print
#'
#' @return
#' @method print KMO
#'
#' @export
print.KMO <- function(x, ...) {

  cat("\n")
  cat("Kaiser-Meyer-Olkin criterion")
  cat("\n")
  cat("\n")

  cat("Overall: ", round(x$KMO, 3))
  cat("\n")
  cat("\n")

  cat("For each variable: ")
  cat("\n")
  print(round(x$KMO_i, 3))

}
