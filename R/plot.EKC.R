#' Plot EKC object
#'
#' Plot method showing a summarized output of the \link{EKC} function
#'
#' @param x a list of class EKC. An output from the \link{EKC} function.
#' @param ... not used.
#'
#' @export
#' @method plot EKC
#'
#' @examples
#' EKC_base <- EKC(test_models$baseline$cormat, N = 500)
#' plot(EKC_base)
#'
plot.EKC <- function(x, ...) {

  eigvls <- x$eigenvalues
  p_eigvls <- pretty(eigvls)
  x_len <- length(eigvls)

  graphics::plot.new()
  graphics::plot.window(xlim = c(1, x_len),
                        ylim = c(min(p_eigvls),
                                 max(p_eigvls)))
  graphics::axis(1, seq_len(x_len))
  graphics::axis(2, p_eigvls, las = 1)

  graphics::mtext("Indicators", side = 1, line = 3, cex = 1.5, padj =-.5)
  graphics::mtext("Eigenvalues", side = 2, line = 3, cex = 1.5, padj =.5)

  graphics::lines(seq_len(x_len), eigvls)
  graphics::points(seq_len(x_len), eigvls, pch = 16)

  graphics::lines(seq_len(x_len), x$references, lty = 2, lwd = 1.25,
                  col = "darkgray")

  if (!is.na(x$n_factors)) {
    graphics::points(x$n_factors, eigvls[x$n_factors],
                     pch = 1, cex = 2, col = "red")
    graphics::text(x$n_factors, eigvls[x$n_factors],
                   x$n_factors, pos = 3, cex = 1.5, col = "red",
                   font = 1, offset = .75)
  }

  factors_text <- paste0("N factors suggested by EKC: ",
                         x$n_factors)

  graphics::title(factors_text, cex.main = 1.3)

}
