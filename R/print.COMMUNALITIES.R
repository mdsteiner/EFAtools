#' A print method for communalities
#'
#' @param x class COMMUNALITY vector
#' @param cutoff numeric. The number above which to print loadings in bold
#'  default is .3.
#' @param digits numeric. Passed to \code{\link[base]{round}}. Number of digits
#'  to round the loadings to (default is 3).
#' @param ... additional arguments passed to print
#' @method print COMMUNALITIES
#' @return A string created using the crayon package
#' @export
print.COMMUNALITIES <- function(x, cutoff = .3, digits = 3, ...) {

  temp_i <- NULL

  for (ii in seq_along(x)) {
    if (x[ii] >= 1 || x[ii] < 0) {
      temp_i <- c(temp_i, crayon::red$bold(.numformat(round(x[ii], digits = digits),
                                                      digits = digits)))
    } else {
      temp_i <- c(temp_i, .numformat(round(x[ii], digits = digits),
                                     digits = digits))
    }
  }

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
