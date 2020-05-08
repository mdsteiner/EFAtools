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

  cat("\n")
  cat("Sequential \U1D712\U00B2 Model Tests suggest ", crayon::bold(nfac_chi),
      " factor", ifelse(nfac_chi > 1, "s.", "."), sep = "")
  cat("\n")
  cat("\n")

  cat("RMSEA lower bound suggests ", crayon::bold(nfac_RMSEA), " factor",
      ifelse(nfac_RMSEA > 1, "s.", "."), sep = "")
  cat("\n")
  cat("\n")

  cat("AIC suggests ", crayon::bold(nfac_AIC), " factor",
      ifelse(nfac_AIC > 1, "s.", "."), sep = "")
  cat("\n")
  cat("\n")


}
