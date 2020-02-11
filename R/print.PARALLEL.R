#' Print function for PARALLEL objects
#'
#' @param x a list of class PARALLEL. Output from \link{PARALLEL} function.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print PARALLEL
#'
#' @examples
#' \dontrun{
#' # example without real data
#' PARALLEL(n_cases = 500, n_vars = 10)
#' # example without correlation matrix
#' PARALLEL(test_models$case_11b$cormat, n_cases = test_models$case_11b$N)
#' }
print.PARALLEL <- function(x, ...) {

  cat(paste("Parallel Analysis performed using", x$ctrl$n_datasets,
               "simulated data sets.\n"))

  cat("Eigenvalues were found using", crayon::bold(x$ctrl$eigen_type))
  cat("\n")


  if (isTRUE(x$ctrl$x_dat)) {

    tt <- paste("Decision rule used:",
                       crayon::bold(x$ctrl$decision_rule))
    cat("\n")
    cat(tt)
    cat("\n")
    cat("\n")
    cat("Number of factors to retain:\n")

    if (x$ctrl$data_type == "sim") {
      ptt <- "Based on simulated data:"
    } else if (x$ctrl$data_type == "resample") {
      ptt <- "Based on resampled data:"
    }
    tt <- paste("    ",cli::symbol$bullet, ptt,
                       crayon::bold(x$n_factors))
    cat(tt)
    cat("\n")

  } else {

    cat("\n")
    cat("No data was entered to base number on factors on. Showing simulated eigenvalues:")
    cat("\n")
    cat("\n")
    print(x$eigenvalues)

  }

}
