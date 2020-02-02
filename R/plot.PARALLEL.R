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
#' x <- PARALLEL(test_models$case_11b$cormat, ncases = test_models$case_11b$N)
#' plot(x)
#' }
plot.PARALLEL <- function(x, ...) {

  graphics::plot.new()
  graphics::plot.window(xlim = c(1, x$ctrl$nvars), ylim = c(0, max(x$eigenvalues) + .25))
  graphics::axis(1, 1:x$ctrl$nvars)
  graphics::axis(2, seq(0, max(x$eigenvalues) + .25, round((max(x$eigenvalues) + .25) / 6, 1)),
       las = 1)
  graphics::mtext("Indicators", side = 1, line = 3, cex = 1.5, padj =-.5)
  graphics::mtext("Eigenvalues", side = 2, line = 3, cex = 1.5, padj =.5)

  if (isTRUE(x$ctrl$x_dat)) {
    graphics::lines(1:x$ctrl$nvars, x$eigenvalues[,"Real Eigenvalues"])
    graphics::points(1:x$ctrl$nvars, x$eigenvalues[,"Real Eigenvalues"], pch = 16)
  }

  cols <- viridisLite::viridis(ncol(x$eigenvalues) - x$ctrl$x_dat, end = .8)
  names(cols) <- colnames(x$eigenvalues)[(as.numeric(x$ctrl$x_dat) + 1):ncol(x$eigenvalues)]
  text_resample <- NULL
  text_sim <- NULL
  if (isTRUE(x$ctrl$sim)) {

    graphics::lines(1:x$ctrl$nvars, x$eigenvalues[,"Means Sim"], lty = 2, lwd = 1.25,
          col = cols[1])

    for (perc_i in x$ctrl$percent) {
      graphics::lines(1:x$ctrl$nvars, x$eigenvalues[,paste(perc_i, "Percentile Sim")],
            lty = 2, col = cols[paste(perc_i, "Percentile Sim")], lwd = 1.25)
    }

    if (!is.null(x$n_factors_sim)) {
      text_sim <- paste0("from simulated data: ", x$n_factors_sim, "; ")
    }
  }

  if (isTRUE(x$ctrl$resample)) {

    graphics::lines(1:x$ctrl$nvars, x$eigenvalues[,"Means Resample"], lty = 2, lwd = 1.25,
          col = cols[1])

    for (perc_i in x$ctrl$percent) {
      graphics::lines(1:x$ctrl$nvars, x$eigenvalues[,paste(perc_i, "Percentile Resample")],
            lty = 2, col = cols[paste(perc_i, "Percentile Resample")], lwd = 1.25)
    }

    if (!is.null(x$n_factors_resample)) {
      text_resample <- paste0("from resampled data: ", x$n_factors_resample, "; ")
    }
  }

  factors_text <- paste0("N Factors: ", text_sim, text_resample,
                         "Decision Rule: ", x$ctrl$factor_rule)

  graphics::mtext(factors_text)

  graphics::legend(round(x$ctrl$nvars / 2), max(x$eigenvalues) - .2,
         colnames(x$eigenvalues), lty = c(1, rep(2, length(cols))),
         col = c("black", cols))


}
