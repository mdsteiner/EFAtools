#' Plot CD object
#'
#' Plot method showing a summarized output of the \link{CD} function
#'
#' @param x a list of class CD. An output from the \link{CD} function.
#' @param ... not used.
#'
#' @export
#' @method plot CD
plot.CD <- function(x, ...) {

  x_len <- x$n_factors

  ys <- colMeans(x$RMSE_eigenvalues)[seq_len(x$n_factors)]

  graphics::plot.new()
  graphics::plot.window(xlim = c(1, x_len),
                        ylim = c(0, max(pretty(ys))))
  graphics::axis(1, seq_len(x_len))
  graphics::axis(2, pretty(c(0, ys)),
                 las = 1)

  graphics::mtext("N Factors", side = 1, line = 3, cex = 1.5, padj =-.5)
  graphics::mtext("RMSE eigenvalues", side = 2, line = 3, cex = 1.5, padj =.25)

  graphics::lines(seq_len(x_len), ys)
  graphics::points(seq_len(x_len), ys, pch = 16)

  if (!is.na(x$n_factors)) {
    graphics::points(x$n_factors, ys[x$n_factors],
                     pch = 1, cex = 2, col = "red")
    graphics::text(x$n_factors, ys[x$n_factors],
                   x$n_factors, pos = 3, cex = 1.5, col = "red",
                   font = 1, offset = .75)
  }

  factors_text <- paste0("N factors suggested by comparison data analysis: ",
                         x$n_factors)

  graphics::title(factors_text, cex.main = 1.3)

}
