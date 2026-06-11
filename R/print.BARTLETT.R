#' Print BARTLETT object
#'
#' @param x list of class BARTLETT (output from the [BARTLETT] function)
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
    cat(.efa_style(cli::symbol$tick, c("green", "bold")), "The",
        .efa_style("Bartlett's test of sphericity", "bold"), "was",
        .efa_style("significant", c("green", "bold")),
    "at an alpha level of .05.")
    cat("\n")
    cat(.efa_style(" ", "bold"), "These data are probably suitable for factor analysis.")

  } else {

    cat("\n")
    cat(.efa_style(cli::symbol$cross, c("red", "bold")),
        "The Bartlett's test of sphericity was",
        .efa_style("not significant", c("red", "bold")),
    "at an alpha level of .05.")
    cat("\n")
    cat(.efa_style(" ", "bold"), "These data are probably not suitable for factor analysis.")

  }

} else {

    cat("\n")
    cat(.efa_style(paste(.efa_style("!", "bold"),
                         "The Bartlett's test of sphericity did not render a result."),
                   "yellow"))
    cat("\n")

}

  cat("\n")
  cat("\n")
  cat(.efa_style("  ", "bold"), "\U1D712\U00B2(", x$df, ") = ", round(x$chisq, 2), ", ",
      .efa_style("p", "italic"), ifelse(pval < .001, " < .001",
                                  paste(" = ", round(pval, 3), sep = "")),
      sep = "")
  cat("\n")

}
