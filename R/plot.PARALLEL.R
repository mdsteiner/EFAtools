#' Plot PARALLEL object
#'
#' Plot method showing a summarized output of the \link{PARALLEL} function
#'
#' @param x list of class PARALLEL. An output from the \link{PARALLEL} function.
#' @param ... not used.
#'
#' @export
#' @method plot PARALLEL
#'
#' @examples
#' \donttest{
#' # example with correlation matrix and "ML" estimation
#' x <- PARALLEL(test_models$case_11b$cormat, N = 500, method = "ML")
#' plot(x)
#' }
plot.PARALLEL <- function(x, ...) {

  n_vars <- x$settings$n_vars
  eigen_type <- x$settings$eigen_type
  eigen_PCA <- x$eigenvalues_PCA
  eigen_SMC <- x$eigenvalues_SMC
  eigen_EFA <- x$eigenvalues_EFA
  n_fac_PCA <- x$n_fac_PCA
  n_fac_SMC <- x$n_fac_SMC
  n_fac_EFA <- x$n_fac_EFA
  x_dat <- x$settings$x_dat
  decision_rule <- x$settings$decision_rule
  percent <- x$settings$percent

  # Create plots depending on eigen_type and if real data were entered or not
  if("PCA" %in% eigen_type){

    .plot_PA_helper(eigenvalues = eigen_PCA, n_vars = n_vars, n_fac = n_fac_PCA,
                    x_dat = x_dat, decision_rule = decision_rule,
                    percent = percent, eigen_type = "PCA")

  }

  if("SMC" %in% eigen_type){

    .plot_PA_helper(eigenvalues = eigen_SMC, n_vars = n_vars, n_fac = n_fac_SMC,
                    x_dat = x_dat, decision_rule = decision_rule,
                    percent = percent, eigen_type = "SMC")

  }

  if("EFA" %in% eigen_type){

    .plot_PA_helper(eigenvalues = eigen_EFA, n_vars = n_vars, n_fac = n_fac_EFA,
                    x_dat = x_dat, decision_rule = decision_rule,
                    percent = percent, eigen_type = "EFA")

  }

}

.plot_PA_helper <- function(eigenvalues, n_vars, n_fac, x_dat,
                            decision_rule = decision_rule, percent = percent,
                            eigen_type){

  p_eigen <- pretty(c(min(eigenvalues) * .9, eigenvalues,
                      max(eigenvalues) * 1.15))
  graphics::plot.new()
  graphics::plot.window(xlim = c(1, n_vars),
                        ylim = c(min(p_eigen), max(p_eigen)))
  graphics::axis(1, seq_len(n_vars))
  graphics::axis(2, p_eigen, las = 1)
  graphics::mtext("Indicators", side = 1, line = 3, cex = 1.5, padj =-.5)
  graphics::mtext("Eigenvalues", side = 2, line = 3, cex = 1.5, padj =.5)

  if (isTRUE(x_dat)) {
    graphics::lines(seq_len(n_vars), eigenvalues[,"Real Eigenvalues"])
    graphics::points(seq_len(n_vars), eigenvalues[,"Real Eigenvalues"], pch = 16)
    if (!is.na(n_fac)) {
      graphics::points(n_fac, eigenvalues[n_fac,"Real Eigenvalues"],
                       pch = 1, cex = 2, col = "red")
      graphics::text(n_fac, eigenvalues[n_fac,"Real Eigenvalues"],
                     n_fac, pos = 3, cex = 1.5, col = "red",
                     font = 1, offset = .75)
    }

  }

  cols <- viridisLite::viridis(ncol(eigenvalues) - x_dat, end = .8)
  names(cols) <- colnames(eigenvalues)[(as.numeric(x_dat) + 1):ncol(eigenvalues)]

  graphics::lines(seq_len(n_vars), eigenvalues[,"Means"], lty = 2, lwd = 1.25,
                  col = cols[1])

  for (perc_i in percent) {
    graphics::lines(seq_len(n_vars), eigenvalues[,paste(perc_i, "Percentile")],
                    lty = 2, col = cols[paste(perc_i, "Percentile")], lwd = 1.25)
  }

  factors_text <- paste0("N Factors with Decision Rule '", decision_rule,
                         "' and Eigen Type '", eigen_type, "': ", n_fac)

  graphics::title(factors_text)

  graphics::legend("topright",
                   colnames(eigenvalues), lty = c(1, rep(2, length(cols))),
                   col = c("black", cols))


}
