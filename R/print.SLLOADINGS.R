#' A print method for factor loadings
#'
#' @param x class SLLOADINGS matrix.
#' @param cutoff numeric. The number above which to print loadings in bold
#'  default is .2.
#' @param digits numeric. Passed to \code{\link[base]{round}}. Number of digits
#'  to round the loadings to (default is 3).
#'
#' @return A string created using the crayon package
#' @export
print.SLLOADINGS <- function(x, cutoff = .2, digits = 3) {

  # so the console is not too crouded, round to digits decimals
  x <- round(x, digits = digits)

  n_col <- ncol(x)

  # create the string to paste using the crayon package
  temp <- apply(x, 1, function(i, cutoff, n_col){
    tt <- NULL
    for (kk in 1:n_col) {
      if (kk <= n_col - 2) {
        if (abs(i[kk]) < cutoff) {
          tt <- c(tt, crayon::silver(i[kk]))
        } else if (abs(i[kk]) <= 1) {
          tt <- c(tt, crayon::bold(i[kk]))
        } else {
          tt <- c(tt, crayon::red$bold(i[kk]))
        }
      } else {
        if (i[kk] <= 1 & i[kk] >= 0) {
          tt <- c(tt, crayon::black(i[kk]))
        } else {
          tt <- c(tt, crayon::red$bold(i[kk]))
        }
      }
    }
    paste(tt, collapse = "\t")
  }, cutoff = cutoff, n_col = n_col)
  temp <- paste(temp, collapse = "\n")

  # warn from Heywood cases
  if (sum(x[, -n_col] > 1) == 1) {
    temp <- paste(temp,
                  crayon::red$bold("\nWarning: Results contain a Heywood case!"),
                  collapse = "\n")
  } else if (sum(x[, -n_col] > 1) > 1) {
    temp <- paste(temp, crayon::red$bold("\nWarning: Results contain",
                                         sum(x[, -n_col] > 1), "Heywood cases!"),
                  collapse = "\n")
  }

  # print the results to the console
  cat(temp)
}
