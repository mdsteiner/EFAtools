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

  if(!is.na(nfac_chi)){

  cat("\n")
  cat("Sequential \U1D712\U00B2 Model Tests suggest ", crayon::bold(nfac_chi),
      " factor", ifelse(nfac_chi > 1 | nfac_chi == 0 | is.na(nfac_chi), "s.", "."),
      sep = "")
  cat("\n")
  cat("\n")

  }

  if(!is.na(nfac_RMSEA)){

  cat("Lower bound of RMSEA 90% confidence interval suggests ",
      crayon::bold(nfac_RMSEA), " factor",
      ifelse(nfac_RMSEA > 1 | nfac_RMSEA == 0 | is.na(nfac_RMSEA), "s.", "."),
      sep = "")
  cat("\n")
  cat("\n")

  }

  if(!is.na(nfac_AIC)){

  cat("AIC suggests ", crayon::bold(nfac_AIC), " factor",
      ifelse(nfac_AIC > 1 | nfac_AIC == 0 | is.na(nfac_AIC), "s.", "."),
      sep = "")
  cat("\n")
  cat("\n")

  }

}
