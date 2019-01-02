#' A print method for factor loadings
#'
#' @param x class LOADINGS matrix.
#' @param cutoff numeric. The number above which to print loadings in bold
#'  default is .3.
#' @param digits numeric. Passed to \code{\link[base]{round}}. Number of digits
#'  to round the loadings to (default is 3).
#'
#' @return A string created using the crayon package
#' @export
print.LOADINGS <- function(x, cutoff = .3, digits = 3) {

  # so the console is not too crouded, round to digits decimals
  x <- round(x, digits = digits)

  # create the string to paste using the crayon package
  temp <- apply(x, 1, function(i, cutoff){
    paste(ifelse(abs(i) < cutoff, crayon::silver(i),
                 ifelse(abs(i) <= 1, crayon::bold(i),
                        crayon::red$bold(i))),
          collapse = "\t")
  }, cutoff = cutoff)
  temp <- paste(temp, collapse = "\n")

  # warn from Heywood cases
  if (sum(x > 1) == 1) {
    temp <- paste(temp, crayon::red$bold("\nWarning: Results contain a Heywood case!"), collapse = "\n")
  } else if (sum(x > 1) > 1) {
    temp <- paste(temp, crayon::red$bold("\nWarning: Results contain",
                                         sum(x > 1), "Heywood cases!"),
                  collapse = "\n")
  }

  # print the results to the console
  cat(temp)
}
