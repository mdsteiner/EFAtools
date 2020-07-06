#' Print function for PARALLEL objects
#'
#' @param x a list of class PARALLEL. Output from \link{PARALLEL} function.
#' @param plot logical. Whether to plot the results.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print PARALLEL
#'
#' @examples
#' \donttest{
#' # example without real data
#' PARALLEL(N = 500, n_vars = 10)
#'
#' # example with correlation matrix and "ML" estimation
#' PARALLEL(test_models$case_11b$cormat, N = 500, method = "ML")
#' }
print.PARALLEL <- function(x, plot = TRUE, ...) {

  eigen_type <- x$settings$eigen_type

  cat("Parallel Analysis performed using ", crayon::bold(x$settings$n_datasets),
      " simulated random data sets", sep = "")
  cat("\n")

  cat("Eigenvalues were found using ", .settings_string(x$settings$eigen_type),
      "", sep = "")
  cat("\n")

  if (isTRUE(x$settings$x_dat)) {

    cat("\n")
    cat("Decision rule used:", crayon::bold(x$settings$decision_rule))
    cat("\n")
    cat("\n")
    cat(cli::rule(left = crayon::bold("Number of factors to retain according to"),
                  col = "blue"))
    cat("\n")
    cat("\n")

    if("PCA" %in% eigen_type){
      cat(crayon::blue(cli::symbol$circle_dotted, "PCA-determined eigenvalues: "),
          crayon::bold(x$n_fac_PCA))
      cat("\n")
    }

    if("SMC" %in% eigen_type){
      cat(crayon::blue(cli::symbol$circle_dotted, "SMC-determined eigenvalues: "),
          crayon::bold(x$n_fac_SMC))
      cat("\n")
    }

    if("EFA" %in% eigen_type){
      cat(crayon::blue(cli::symbol$circle_dotted, "EFA-determined eigenvalues: "),
          crayon::bold(x$n_fac_EFA))
      cat("\n")
    }

  } else {

    cat("\n")
    cat(cli::col_cyan(cli::symbol$info, " No data were entered to base number of factors on. Plotting simulated eigenvalues."))
    cat("\n")

  }

  if (isTRUE(plot)) {
    graphics::plot(x)
  }

}
