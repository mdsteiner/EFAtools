#' Print KMO object
#'
#' @param x list of class KMO (output from the \link{KMO} function)
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
  cat(cli::rule(left = "Kaiser-Meyer-Olkin criterion (KMO)", col = "blue"))
  cat("\n")

  if(!is.na(KMO) && !is.null(KMO)){

    if(KMO >= .9){
      symb <- crayon::green$bold(cli::symbol$tick)
      label <- crayon::green$bold("marvellous")
    } else if(KMO >= .8){
      symb <- crayon::green$bold(cli::symbol$tick)
      label <- crayon::green$bold("meritorious")
    } else if(KMO >= .7){
      symb <- crayon::green$bold(cli::symbol$tick)
      label <- crayon::green$bold("middling")
    } else if(KMO >= .6){
      symb <- crayon::yellow$bold("!")
      label <- crayon::yellow$bold("mediocre")
    } else if (KMO >= .5){
      symb <- crayon::red$bold(cli::symbol$cross)
      label <- crayon::red$bold("miserable")
    } else {
      symb <- crayon::red$bold(cli::symbol$cross)
      label <- crayon::red$bold("unacceptable")
    }

    cat("\n")
    cat(symb, " The overall KMO value for your data is ", label, ".", sep = "")
    cat("\n")

    if(KMO < .5){
      cat(crayon::bold(" "), "These data are not suitable for factor analysis.")
      cat("\n")
    } else if(KMO < .6){
      cat(crayon::bold(" "), "These data are hardly suitable for factor analysis.")
      cat("\n")
    } else {
      cat(crayon::bold(" "), "These data are probably suitable for factor analysis.")
      cat("\n")
    }

    cat("\n")
    cat(crayon::bold(" "), crayon::blue("Overall:"), crayon::bold(round(KMO, 3)))
    cat("\n")
    cat("\n")

    cat(crayon::bold(" "), crayon::blue("For each variable:"))
    cat("\n")
    print(round(x$KMO_i, 3))

  } else {

    cat("\n")
    cat(crayon::yellow(crayon::bold("!"), "Sorry, the KMO value for your data is not available."))
    cat("\n")

  }

}
