#' Plot PARALLEL object
#'
#' Plot Method showing a summarized output of the \link{PARALLEL} function
#'
#' @param x list of class PARALLEL. An output from the \link{PARALLEL} function.
#' @param ... not used.
#'
#' @export
#' @method plot PARALLEL
#'
#' @examples
#' \dontrun{
#' # example without correlation matrix
#' x <- PARALLEL(test_models$case_11b$cormat, n_cases = test_models$case_11b$N)
#' plot(x)
#' }
plot.PARALLEL <- function(x, ...) {

  graphics::plot.new()
  graphics::plot.window(xlim = c(1, x$ctrl$n_vars), ylim = c(0, max(x$eigenvalues) + .25))
  graphics::axis(1, 1:x$ctrl$n_vars)
  graphics::axis(2, seq(0, max(x$eigenvalues) + .25, round((max(x$eigenvalues) + .25) / 6, 1)),
       las = 1)
  graphics::mtext("Indicators", side = 1, line = 3, cex = 1.5, padj =-.5)
  graphics::mtext("Eigenvalues", side = 2, line = 3, cex = 1.5, padj =.5)

  if (isTRUE(x$ctrl$x_dat)) {
    graphics::lines(1:x$ctrl$n_vars, x$eigenvalues[,"Real Eigenvalues"])
    graphics::points(1:x$ctrl$n_vars, x$eigenvalues[,"Real Eigenvalues"], pch = 16)
  }

  cols <- viridisLite::viridis(ncol(x$eigenvalues) - x$ctrl$x_dat, end = .8)
  names(cols) <- colnames(x$eigenvalues)[(as.numeric(x$ctrl$x_dat) + 1):ncol(x$eigenvalues)]

    graphics::lines(1:x$ctrl$n_vars, x$eigenvalues[,"Means"], lty = 2, lwd = 1.25,
          col = cols[1])

    for (perc_i in x$ctrl$percent) {
      graphics::lines(1:x$ctrl$n_vars, x$eigenvalues[,paste(perc_i, "Percentile")],
            lty = 2, col = cols[paste(perc_i, "Percentile")], lwd = 1.25)
    }

    if (x$ctrl$data_type == "sim") {
      text <- paste0("from simulated data: ", x$n_factors, "; ")
    } else if (x$ctrl$data_type == "sim" && isTRUE(x$ctrl$replace)) {
      text <- paste0("from resampled data with replacement: ", x$n_factors, "; ")
    } else if (x$ctrl$data_type == "sim" && isFALSE(x$ctrl$replace)) {
      text <- paste0("from resampled data without replacement: ", x$n_factors, "; ")
    }

  factors_text <- paste0("N Factors: ", text,
                         "Decision Rule: ", x$ctrl$decision_rule)

  graphics::mtext(factors_text)

  graphics::legend(round(x$ctrl$n_vars / 2), max(x$eigenvalues) - .2,
         colnames(x$eigenvalues), lty = c(1, rep(2, length(cols))),
         col = c("black", cols))


}
