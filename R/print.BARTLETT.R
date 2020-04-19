#' Print BARTLETT object
#'
#' @param x list of class BARTLETT (output from the \link{BARTLETT} function)
#' @param ... additional arguments passed to print
#'
#' @return
#' @method print BARTLETT
#'
#' @export
print.BARTLETT <- function(x, ...) {

  pval <- x$p_value

  if(!is.na(pval) && !is.null(pval)){

  if(pval < .05){

    cat("\n")
    cat("The Bartletts test of sphericity was", crayon::green$bold("significant"),
    "at an alpha level of .05.")
    cat("\n")
    cat("These data are probably suitable for factor analysis.")

  } else {

    cat("\n")
    cat("The Bartletts test of sphericity was", crayon::red$bold("not significant"),
    "at an alpha level of .05.")
    cat("\n")
    cat("These data are probably not suitable for factor analysis.")

  }

} else {

    cat("\n")
    cat("The Bartletts test of sphericity did not render a result.")
    cat("\n")

}

  cat("\n")
  cat("\n")
  cat("\U1D712\U00B2(", x$df, ") = ", round(x$chisq, 2), ", ", crayon::italic("p"),
      " = ", round(pval, 4), sep = "")
  cat("\n")

}