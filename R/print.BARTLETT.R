#' Print BARTLETT object
#'
#' @param x list of class BARTLETT (output from the \link{BARTLETT} function)
#' @param ... additional arguments passed to print
#'
#' @method print BARTLETT
#'
#' @export
#'
#' @examples
#' BARTLETT(test_models$baseline$cormat, N = 500)
#'
print.BARTLETT <- function(x, ...) {

  pval <- x$p_value

  if(!is.na(pval) && !is.null(pval)){

  if(pval < .05){

    cat("\n")
    cat(crayon::green$bold(cli::symbol$tick), "The",
        crayon::bold("Bartlett's test of sphericity"), "was",
        crayon::green$bold("significant"),
    "at an alpha level of .05.")
    cat("\n")
    cat(crayon::bold(" "), "These data are probably suitable for factor analysis.")

  } else {

    cat("\n")
    cat(crayon::red$bold(cli::symbol$cross),
        "The Bartlett's test of sphericity was", crayon::red$bold("not significant"),
    "at an alpha level of .05.")
    cat("\n")
    cat(crayon::bold(" "), "These data are probably not suitable for factor analysis.")

  }

} else {

    cat("\n")
    cat(crayon::yellow(crayon::bold("!"), "The Bartlett's test of sphericity did not render a result."))
    cat("\n")

}

  cat("\n")
  cat("\n")
  cat(crayon::bold("  "), "\U1D712\U00B2(", x$df, ") = ", round(x$chisq, 2), ", ",
      crayon::italic("p"), ifelse(pval < .001, " < .001",
                                  paste(" = ", round(pval, 3), sep = "")),
      sep = "")
  cat("\n")

}
