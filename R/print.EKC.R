#' Print function for EKC objects
#'
#' @param x a list of class EKC. Output from \code{\link{EKC}} function.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print EKC
#'
#' @examples
#' EKC_base <- EKC(test_models$baseline$cormat, N = 500)
#' EKC_base
#'
print.EKC <- function(x, ...) {

  nfac <- x$n_factors

  cat("\n")
  cat("Empirical Kaiser criterion suggests ", crayon::bold(nfac), " factor",
        ifelse(nfac > 1 | nfac == 0 | is.na(nfac), "s.", "."), sep = "")
  cat("\n")
  cat("\n")

  graphics::plot(x)

}
