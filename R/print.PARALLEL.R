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
#' PARALLEL(ncases = 500, nvars = 10)
#' # example without correlation matrix
#' PARALLEL(test_models$case_11b$cormat, ncases = test_models$case_11b$N)
#' }
print.PARALLEL <- function(x, ...) {

  cat(paste("Parallel Analysis performed using", x$ctrl$ndatasets,
               "simulated data sets.\n"))

  if (isTRUE(x$ctrl$x_dat) && isTRUE(x$ctrl$sim) || isTRUE(x$ctrl$resample)) {

    tt <- paste("Decision rule used:",
                       crayon::bold(x$ctrl$factor_rule))
    cat("\n")
    cat(tt)
    cat("\n")
    cat("\n")
    cat("Number of factors to retain:\n")

    if (isTRUE(x$ctrl$sim)) {

      tt <- paste("    ",cli::symbol$bullet, "Based on simulated data:",
                         crayon::bold(x$n_factors_sim))
      cat(tt)
      cat("\n")

    }

    if (isTRUE(x$ctrl$resample)) {
      tt <- paste("    ",cli::symbol$bullet, "Based on resampled data:",
                         crayon::bold(x$n_factors_resample))
      cat(tt)
      cat("\n")
    }

  } else {

    cat("\n")
    cat("No data was entered to base number on factors on. Showing simulated eigenvalues:")
    cat("\n")
    cat("\n")
    print(x$eigenvalues)

  }

}
