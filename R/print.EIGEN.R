#' Print EIGEN object
#'
#' @param x class EIGEN vector
#' @param cutoff numeric. The number above which to print loadings in bold
#'  default is .3.
#' @param digits numeric. Passed to \code{\link[base]{round}}. Number of digits
#'  to round the loadings to (default is 3).
#' @param ... additional arguments passed to print
#'
#' @return A string created using the crayon package
#' @method print EIGEN
#' @export
print.EIGEN <- function(x, cutoff = .3, digits = 3, ...) {

  x <- round(x, digits = digits)
  temp_i <- ifelse(x >= 1, crayon::bold(x), ifelse(x < 0, crayon::red(x),
                                                 crayon::silver(x)))

  for (ss in seq(1, length(x), 7)) {
    if (length(x) > ss + 6) {
      tt <- ss + 6
    } else {
      tt <- length(x)
    }
    if (ss == 1) {
      temp <- stringr::str_c(temp_i[ss:tt], collapse = "  ")
    } else {
      temp <- stringr::str_c(temp, "\n", stringr::str_c(temp_i[ss:tt],
                                                        collapse = "  "))
    }

  }

  temp <- stringr::str_c(temp, "\n")

  # print the results to the console
  cat(temp)
}
