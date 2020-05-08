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

  eigen_PCA <- x$eigen_PCA
  eigen_SMC <- x$eigen_SMC
  eigen_EFA <- x$eigen_EFA
  nfac_PCA <- x$n_fac_PCA
  nfac_SMC <- x$n_fac_SMC
  nfac_EFA <- x$n_fac_EFA

  if(!is.na(nfac_PCA)){

  x_len <- length(eigen_PCA)

  graphics::plot.new()
  graphics::plot.window(xlim = c(1, x_len),
                        ylim = c(min(eigen_PCA) - .2,
                                 max(eigen_PCA) + .5))
  graphics::axis(1, 1:x_len)
  graphics::axis(2, round(seq(min(eigen_PCA) - .2,
                              max(eigen_PCA) + .5,
                              round(diff(c(min(eigen_PCA) - .2,
                                           max(eigen_PCA)) + .5) / 6, 1)), 1),
                 las = 1)

  graphics::mtext("Indicators", side = 1, line = 3, cex = 1.5, padj =-.5)
  graphics::mtext("Eigenvalues", side = 2, line = 3, cex = 1.5, padj =.5)

  graphics::lines(1:x_len, eigen_PCA)
  graphics::points(1:x_len, eigen_PCA, pch = 16)
  graphics::abline(h = 1, lty = 2)

  graphics::points(nfac_PCA, eigen_PCA[nfac_PCA], pch = 1, cex = 2, col = "red")

  graphics::text(nfac_PCA, eigen_PCA[nfac_PCA], nfac_PCA,
                 pos = 4, cex = 1.2, col = "red",
                 font = 1, offset = .75)

  title <- paste0("N factors suggested by Kaiser-Guttman criterion with PCA: ",
                  nfac_PCA)

  graphics::title(title, cex.main = 1.2)

  }

  if (!is.na(nfac_SMC)) {

    x_len <- length(eigen_SMC)

    graphics::plot.new()
    graphics::plot.window(xlim = c(1, x_len),
                          ylim = c(min(eigen_SMC) - .2,
                                   max(eigen_SMC) + .5))
    graphics::axis(1, 1:x_len)
    graphics::axis(2, round(seq(min(eigen_SMC) - .2,
                                max(eigen_SMC) + .5,
                                round(diff(c(min(eigen_SMC) - .2,
                                             max(eigen_SMC)) + .5) / 6, 1)), 1),
                   las = 1)

    graphics::mtext("Indicators", side = 1, line = 3, cex = 1.5, padj =-.5)
    graphics::mtext("Eigenvalues", side = 2, line = 3, cex = 1.5, padj =.5)

    graphics::lines(1:x_len, eigen_SMC)
    graphics::points(1:x_len, eigen_SMC, pch = 16)
    graphics::abline(h = 1, lty = 2)

    graphics::points(nfac_SMC, eigen_SMC[nfac_SMC], pch = 1, cex = 2,
                     col = "red")

    graphics::text(nfac_SMC, eigen_SMC[nfac_SMC], nfac_SMC, pos = 4, cex = 1.2,
                   col = "red", font = 1, offset = .75)

    title <- paste0("N factors suggested by Kaiser-Guttman criterion with SMC: ",
                    nfac_SMC)

    graphics::title(title, cex.main = 1.2)

  }

  if (!is.na(nfac_EFA)) {

    x_len <- length(eigen_EFA)

    graphics::plot.new()
    graphics::plot.window(xlim = c(1, x_len),
                          ylim = c(min(eigen_EFA) - .2,
                                   max(eigen_EFA) + .5))
    graphics::axis(1, 1:x_len)
    graphics::axis(2, round(seq(min(eigen_EFA) - .2,
                                max(eigen_EFA) + .5,
                                round(diff(c(min(eigen_EFA) - .2,
                                             max(eigen_EFA)) + .5) / 6, 1)), 1),
                   las = 1)

    graphics::mtext("Indicators", side = 1, line = 3, cex = 1.5, padj =-.5)
    graphics::mtext("Eigenvalues", side = 2, line = 3, cex = 1.5, padj =.5)

    graphics::lines(1:x_len, eigen_EFA)
    graphics::points(1:x_len, eigen_EFA, pch = 16)
    graphics::abline(h = 1, lty = 2)


    graphics::points(nfac_EFA, eigen_EFA[nfac_EFA], pch = 1, cex = 2,
                         col = "red")
    graphics::text(nfac_EFA, eigen_EFA[nfac_EFA], nfac_EFA, pos = 4, cex = 1.2,
                   col = "red", font = 1, offset = .75)

    title <- paste0("N factors suggested by Kaiser-Guttman criterion with EFA: ",
                    nfac_EFA)

    graphics::title(title, cex.main = 1.2)

  }


}
