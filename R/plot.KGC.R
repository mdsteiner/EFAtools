#' Plot KGC object
#'
#' Plot method showing a summarized output of the \link{KGC} function
#'
#' @param x a list of class KGC. An output from the \link{KGC} function.
#' @param ... not used.
#'
#' @export
#' @method plot KGC
#'
#' @examples
#' KGC_base <- KGC(test_models$baseline$cormat, eigen_type = "PCA")
#' plot(KGC_base)
#'
plot.KGC <- function(x, ...) {

  eigvls <- x$eigenvalues
  x_len <- length(eigvls)

  graphics::plot.new()
  graphics::plot.window(xlim = c(1, x_len),
                        ylim = c(min(eigvls) - .2,
                                 max(eigvls) + .5))
  graphics::axis(1, 1:x_len)
  graphics::axis(2, round(seq(min(eigvls) - .2,
                              max(eigvls) + .2,
                              round(diff(c(min(eigvls) - .2,
                                           max(eigvls)) + .2) / 6, 1)), 1),
                 las = 1)

  graphics::mtext("Indicators", side = 1, line = 3, cex = 1.5, padj =-.5)
  graphics::mtext("Eigenvalues", side = 2, line = 3, cex = 1.5, padj =.5)

  graphics::lines(1:x_len, eigvls)
  graphics::points(1:x_len, eigvls, pch = 16)
  graphics::abline(h = 1, lty = 2)

  if (!is.na(x$n_factors)) {
      graphics::points(x$n_factors, eigvls[x$n_factors],
                       pch = 1, cex = 2, col = "red")
      graphics::text(x$n_factors, eigvls[x$n_factors],
                     x$n_factors, pos = 3, cex = 1.5, col = "red",
                     font = 1, offset = .75)
    }

  factors_text <- paste0("N factors suggested by Kaiser-Guttman criterion: ",
                         x$n_factors)

  graphics::title(factors_text, cex.main = 1.3)

}
