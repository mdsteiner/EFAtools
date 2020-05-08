#' Print function for KGC objects
#'
#' @param x a list of class KGC. Output from \link{KGC} function.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print KGC
#'
#' @examples
#' KGC_base <- KGC(test_models$baseline$cormat, eigen_type = "PCA")
#' KGC_base
#'
print.KGC <- function(x, ...) {

  nfac_PCA <- x$n_fac_PCA
  nfac_SMC <- x$n_fac_SMC
  nfac_EFA <- x$n_fac_EFA
  eigen_type <-x$settings$eigen_type

  cat("\n")
  cat("Eigenvalues were found using ", .settings_string(eigen_type), sep = "")
  cat("\n")
  cat("\n")

  cat("Kaiser-Guttmann criterion suggests:")
  cat("\n")
  cat("\n")

  if(!is.na(nfac_PCA)){

    cat("   ", crayon::bold(nfac_PCA), " factor", ifelse(nfac_PCA > 1, "s", ""),
      " with PCA-determined eigenvalues", sep = "")
    cat("\n")

  }

  if(!is.na(nfac_SMC)){

    cat("   ", crayon::bold(nfac_SMC), " factor", ifelse(nfac_SMC > 1, "s", ""),
        " with SMC-determined eigenvalues", sep = "")
    cat("\n")

  }

  if(!is.na(nfac_EFA)){

    cat("   ", crayon::bold(nfac_EFA), " factor", ifelse(nfac_EFA > 1, "s", ""),
        " with EFA-determined eigenvalues", sep = "")
    cat("\n")

  }

  cat("\n")

  graphics::plot(x)

}
