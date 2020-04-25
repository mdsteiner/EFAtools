#' Print KMO object
#'
#' @param x list of class KMO (output from the \link{KMO} function)
#' @param ... additional arguments passed to print
#'
#' @return
#' @method print KMO
#'
#' @export
#'
#' @example
#' KMO_base <- KMO(test_models$baseline$cormat)
#' KMO_base
#'
print.KMO <- function(x, ...) {

  KMO <- x$KMO

  cat("\n")
  cat("Kaiser-Meyer-Olkin criterion (KMO)")
  cat("\n")

  if(!is.na(KMO) && !is.null(KMO)){

    if(KMO >= .9){
      label = crayon::green$bold("marvellous")
    } else if(KMO >= .8){
      label = crayon::green$bold("meritorious")
    } else if(KMO >= .7){
      label = crayon::green$bold("middling")
    } else if(KMO >= .6){
      label = crayon::orange$bold("mediocre")
    } else if (KMO >= .5){
      label = crayon::red$bold("miserable")
    } else {
      label = crayon::red$bold("unacceptable")
    }

    cat("\n")
    cat("The overall KMO value for your data is ", label, ".", sep = "")
    cat("\n")

    if(KMO < .5){
      cat("These data are not suitable for factor analysis.")
      cat("\n")
    } else if(KMO < .6){
      cat("These data are hardly suitable for factor analysis.")
      cat("\n")
    }

  } else {

    cat("\n")
    cat("The overall KMO value for your data is not available.")
    cat("\n")

  }

  cat("\n")
  cat("Overall:", crayon::bold(round(KMO, 3)))
  cat("\n")
  cat("\n")

  cat("For each variable:")
  cat("\n")
  print(round(x$KMO_i, 3))

}
