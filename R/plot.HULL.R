#' Plot HULL object
#'
#' Plot method showing a summarized output of the \code{\link{HULL}} function
#'
#' @param x list of class HULL. An output from the \code{\link{HULL}} function.
#' @param ... not used.
#'
#' @export
#' @method plot HULL
#'
#' @examples
#' \donttest{
#' x <- HULL(test_models$baseline$cormat, N = 500, method = "ML")
#' plot(x)
#'}
plot.HULL <- function(x, ...) {

  gof <- x$settings$gof
  method <- x$settings$method
  sol_CAF <- as.data.frame(x$solutions_CAF)
  sol_CFI <- as.data.frame(x$solutions_CFI)
  sol_RMSEA <- as.data.frame(x$solutions_RMSEA)
  n_fac_CAF <- x$n_fac_CAF
  n_fac_CFI <- x$n_fac_CFI
  n_fac_RMSEA <- x$n_fac_RMSEA

  if(!all(is.na(sol_CAF)) && !all(is.na(sol_CAF[, 4]))){

    .plot_HULL_helper(sol_CAF, method = method, gof = "CAF")

  }

  if(!all(is.na(sol_CFI)) && !all(is.na(sol_CFI[, 4]))){

    .plot_HULL_helper(sol_CFI, method = method, gof = "CFI")

  }

  if(!all(is.na(sol_RMSEA)) && !all(is.na(sol_RMSEA[, 4]))){

    .plot_HULL_helper(sol_RMSEA, method = method, gof = "RMSEA")

  }


}

.plot_HULL_helper <- function(sol, method, gof){
  # Sol = solution output for given fit index
  # Method = estimation method of EFA (PAF, ULS, ML)
  # gof = the fit index
  # nfac_suggested = number of factors suggested by HULL with given fit index

  df <- sol$df
  p_df <- pretty(df)

  fit_ind <- sol[[gof]]
  p_fit_ind <- pretty(c(0, fit_ind, max(fit_ind) * 1.25))
  st <- sol$st
  nfac_tested <- sol$nfactors

  graphics::plot.new()
  graphics::plot.window(xlim = c(min(p_df), max(p_df)),
                        ylim = c(min(p_fit_ind), max(p_fit_ind)),
                        xaxs = "i")
  graphics::axis(1, p_df)
  graphics::axis(2, p_fit_ind, las = 1)
  graphics::mtext(expression(italic(df)), side = 1, line = 3, cex = 1.5, padj =-.5)
  graphics::mtext(expression(italic(f)), side = 2, line = 3, cex = 1.5, padj =.5)
  graphics::title(paste("Hull Method with", method, "estimation and", gof))

  graphics::points(df[is.na(st)], fit_ind[is.na(st)], pch = 16,
                   col = "darkgrey")
  graphics::points(df[!is.na(st)], fit_ind[!is.na(st)], pch = 16,
                   cex = 1.5, type = "b")
  graphics::points(df[which.max(st)], fit_ind[which.max(st)], pch = 16,
                   cex = 1.5)
  graphics::points(df[which.max(st)], fit_ind[which.max(st)], pch = 1,
                   cex = 2.5, col = "red")
  if (any(is.na(st))) {
    graphics::text(df[is.na(st)], fit_ind[is.na(st)],
                   nfac_tested[is.na(st)], pos = 1, col = "darkgrey")
  }

  graphics::text(df[!is.na(st) & nfac_tested != which.max(st) - 1],
                 fit_ind[!is.na(st) & nfac_tested != which.max(st) - 1],
                 nfac_tested[!is.na(st) & nfac_tested != which.max(st) - 1],
                 pos = 3, cex = 1.25)
  graphics::text(df[which.max(st)], fit_ind[which.max(st)],
                 nfac_tested[which.max(st)],
                 pos = 3, cex = 1.5, col = "red", font = 1, offset = .75)


}
