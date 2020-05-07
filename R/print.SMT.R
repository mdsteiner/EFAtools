#' Print SMT object
#'
#' @param x list of class SMT (output from the \link{SMT} function)
#' @param ... additional arguments passed to print
#'
#' @method print SMT
#'
#' @export
#'
#' @examples
#' SMT_base <- SMT(test_models$baseline$cormat, N = 500)
#' SMT_base
#'
print.SMT <- function(x, ...) {

  nfac_chi <- x$nfac_chi
  nfac_RMSEA <- x$nfac_RMSEA
  nfac_AIC <- x$nfac_AIC

  if(nfac_chi == 1){

    cat("\n")
    cat("Sequential \U1D712\U00B2 Model Tests suggest", crayon::bold(nfac_chi),
        "factor")
    cat("\n")
    cat("\n")

  } else {

    cat("\n")
    cat("Sequential \U1D712\U00B2 Model Tests suggest", crayon::bold(nfac_chi),
        "factors")
    cat("\n")
    cat("\n")

  }

  if(nfac_RMSEA == 1){

  cat("RMSEA lower bound suggests", crayon::bold(nfac_RMSEA), "factor")
  cat("\n")
  cat("\n")

  } else {

    cat("RMSEA lower bound suggests", crayon::bold(nfac_RMSEA), "factors")
    cat("\n")
    cat("\n")

  }

  if(nfac_AIC == 1){

    cat("AIC suggests", crayon::bold(nfac_AIC), "factor")
    cat("\n")
    cat("\n")

  } else {

    cat("AIC suggests", crayon::bold(nfac_AIC), "factors")
    cat("\n")
    cat("\n")

  }

}
