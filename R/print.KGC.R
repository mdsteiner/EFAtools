#' Print function for KGC objects
#'
#' @param x a list of class KGC. Output from \link{KGC} function.
#' @param plot logical. Whether to plot the results.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print KGC
#'
#' @examples
#' KGC_base <- KGC(test_models$baseline$cormat)
#' KGC_base
#'
print.KGC <- function(x, plot = TRUE, ...) {

  nfac_PCA <- x$n_fac_PCA
  nfac_SMC <- x$n_fac_SMC
  nfac_EFA <- x$n_fac_EFA
  eigen_type <-x$settings$eigen_type

  cat("\n")
  cat("Eigenvalues were found using ", .settings_string(eigen_type), ".",
      sep = "")
  cat("\n")
  cat("\n")

  cat(cli::rule("Number of factors suggested by Kaiser-Guttmann criterion",
                col = "blue"))
  cat("\n")
  cat("\n")

  if("PCA" %in% eigen_type){

    cat(crayon::blue(cli::symbol$circle_dotted, "With PCA-determined eigenvalues: "),
                     crayon::bold(nfac_PCA))
    cat("\n")

  }

  if("SMC" %in% eigen_type){

    cat(crayon::blue(cli::symbol$circle_dotted,
                     "With SMC-determined eigenvalues: "),
        crayon::bold(nfac_SMC))
    cat("\n")

  }

  if("EFA" %in% eigen_type){

    cat(crayon::blue(cli::symbol$circle_dotted,
                     "With EFA-determined eigenvalues: "),
        crayon::bold(nfac_EFA))
    cat("\n")

  }

  cat("\n")

  if (isTRUE(plot)) {
    graphics::plot(x)
  }

}
