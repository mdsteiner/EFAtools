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
#' KGC_base <- KGC(test_models$baseline$cormat)
#' plot(KGC_base)
#'
plot.KGC <- function(x, ...) {

  eigen_PCA <- x$eigen_PCA
  eigen_SMC <- x$eigen_SMC
  eigen_EFA <- x$eigen_EFA
  nfac_PCA <- x$n_fac_PCA
  nfac_SMC <- x$n_fac_SMC
  nfac_EFA <- x$n_fac_EFA

  # Create plots
  if(!is.na(nfac_PCA)){

    .plot_KGC_helper(eigen = eigen_PCA, n_fac = nfac_PCA, eigen_name = "PCA")

  }

  if (!is.na(nfac_SMC)) {

    .plot_KGC_helper(eigen = eigen_SMC, n_fac = nfac_SMC, eigen_name = "SMC")

  }

  if (!is.na(nfac_EFA)) {

    .plot_KGC_helper(eigen = eigen_EFA, n_fac = nfac_EFA, eigen_name = "EFA")

  }

}

.plot_KGC_helper <- function(eigen, n_fac, eigen_name){
  # eigen = eigenvalues found with specific type (PCA, SMC or EFA)
  # n_fac = number of factors suggested with this type
  # eigen_name = name of type (one of "PCA", "SMC", "EFA")

x_len <- length(eigen)

p_eigen <- pretty(c(min(eigen) * .9, eigen, max(eigen) * 1.1))

graphics::plot.new()
graphics::plot.window(xlim = c(1, x_len),
                      ylim = c(min(p_eigen), max(p_eigen)))
graphics::axis(1, seq_len(x_len))
graphics::axis(2, p_eigen, las = 1)

graphics::mtext("Indicators", side = 1, line = 3, cex = 1)
graphics::mtext("Eigenvalues", side = 2, line = 3, cex = 1, padj =.5)

graphics::lines(seq_len(x_len), eigen)
graphics::points(seq_len(x_len), eigen, pch = 16)
graphics::abline(h = 1, lty = 2)

graphics::points(n_fac, eigen[n_fac], pch = 1, cex = 2, col = "red")

graphics::text(n_fac, eigen[n_fac], n_fac,
               pos = 4, cex = 1.2, col = "red",
               font = 1, offset = .75)

title <- paste0("N factors suggested by Kaiser-Guttman criterion with ",
                eigen_name, ": ", n_fac)

graphics::title(title, cex.main = 1.2)

}
