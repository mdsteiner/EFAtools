#' Print KMO object
#'
#' @param x list of class KMO (output from the [KMO] function)
#' @param ... additional arguments passed to print
#'
#' @method print KMO
#'
#' @export
#'
#' @examples
#' KMO_base <- KMO(test_models$baseline$cormat)
#' KMO_base
#'
print.KMO <- function(x, ...) {

  KMO <- x$KMO

  cat("\n")
  cat(cli::rule(left = "Kaiser-Meyer-Olkin criterion (KMO)"))
  cat("\n")

  if(!is.na(KMO) && !is.null(KMO)){

    if(KMO >= .9){
      symb <- .efa_style(cli::symbol$tick, c("green", "bold"))
      label <- .efa_style("marvellous", c("green", "bold"))
    } else if(KMO >= .8){
      symb <- .efa_style(cli::symbol$tick, c("green", "bold"))
      label <- .efa_style("meritorious", c("green", "bold"))
    } else if(KMO >= .7){
      symb <- .efa_style(cli::symbol$tick, c("green", "bold"))
      label <- .efa_style("middling", c("green", "bold"))
    } else if(KMO >= .6){
      symb <- .efa_style("!", c("yellow", "bold"))
      label <- .efa_style("mediocre", c("yellow", "bold"))
    } else if (KMO >= .5){
      symb <- .efa_style(cli::symbol$cross, c("red", "bold"))
      label <- .efa_style("miserable", c("red", "bold"))
    } else {
      symb <- .efa_style(cli::symbol$cross, c("red", "bold"))
      label <- .efa_style("unacceptable", c("red", "bold"))
    }

    cat("\n")
    cat(symb, " The overall KMO value for your data is ", label, ".", sep = "")
    cat("\n")

    if(KMO < .5){
      cat(.efa_style(" ", "bold"), "These data are not suitable for factor analysis.")
      cat("\n")
    } else if(KMO < .6){
      cat(.efa_style(" ", "bold"), "These data are hardly suitable for factor analysis.")
      cat("\n")
    } else {
      cat(.efa_style(" ", "bold"), "These data are probably suitable for factor analysis.")
      cat("\n")
    }

    cat("\n")
    cat(.efa_style(" ", "bold"), "Overall:", .efa_style(round(KMO, 3), "bold"))
    cat("\n")
    cat("\n")

    cat(.efa_style(" ", "bold"), "For each variable:")
    cat("\n")
    print(round(x$KMO_i, 3))

  } else {

    cat("\n")
    cat(.efa_style(paste(.efa_style("!", "bold"),
                         "Sorry, the KMO value for your data is not available."),
                   "yellow"))
    cat("\n")

  }

}
