#' Plot HULL object
#'
#' Plot Method showing a summarized output of the \link{HULL} function
#'
#' @param x list of class HULL An output from the \link{HULL} function.
#' @param ... not used.
#'
#' @export
#' @method plot HULL
#'
#' @examples
#' \dontrun{
#' # example without correlation matrix
#' x <- HULL(IDS2_R, n_cases = 2000)
#' plot(x)
#' }
plot.HULL <- function(x, ...) {

  sol <- as.data.frame(x$solutions)

  graphics::plot.new()
  graphics::plot.window(xlim = c(min(sol$df) - min(sol$df) / 100 * 10,
                                 max(sol$df) + max(sol$df) / 100 * 10),
                        ylim = c(0, round(max(sol$f), 1) + .1), xaxs = "i",
                        yaxs = "i")
  graphics::axis(1, c(0, sol$df, round(max(sol$df) + max(sol$df) / 100 * 10)))
  graphics::axis(2, seq(0, round(max(sol$f), 1) + .1, round((round(max(sol$f), 1) + .1) / 6, 1)),
                 las = 1)
  graphics::mtext(expression(italic(df)), side = 1, line = 3, cex = 1.5, padj =-.5)
  graphics::mtext(expression(italic(f)), side = 2, line = 3, cex = 1.5, padj =.5)
  graphics::title("Hull Method with CAF")

  graphics::points(sol$df[is.na(sol$st)], sol$f[is.na(sol$st)], pch = 16,
                   col = "darkgrey")
  graphics::points(sol$df[!is.na(sol$st)], sol$f[!is.na(sol$st)], pch = 16,
                   cex = 1.5, type = "b")
  graphics::points(sol$df[which.max(sol$st)], sol$f[which.max(sol$st)], pch = 16,
                   cex = 1.5)
  graphics::points(sol$df[which.max(sol$st)], sol$f[which.max(sol$st)], pch = 1,
                   cex = 2.5, col = "red")
  graphics::text(sol$df[is.na(sol$st)], sol$f[is.na(sol$st)],
                 sol$`n factors`[is.na(sol$st)], pos = 1, col = "darkgrey")
  graphics::text(sol$df[!is.na(sol$st) & sol$`n factors` != which.max(sol$st) - 1],
                 sol$f[!is.na(sol$st) & sol$`n factors` != which.max(sol$st) - 1],
                 sol$`n factors`[!is.na(sol$st) & sol$`n factors` != which.max(sol$st) - 1],
                 pos = 3, cex = 1.25)
  graphics::text(sol$df[which.max(sol$st)], sol$f[which.max(sol$st)],
                 sol$`n factors`[which.max(sol$st)],
                 pos = 3, cex = 1.5, col = "red", font = 1, offset = .75)

}
