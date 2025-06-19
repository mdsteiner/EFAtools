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

  ekc_type <- x$settings$type

  if("BvA2017" %in% ekc_type){

    .plot_EKC_helper(eigvls = x$eigenvalues,
                     references = x$references$BvA2017,
                     n_factors = x$n_factors_BvA2017,
                     type = "BvA2017")

  }

  if("AM2019" %in% ekc_type){

    .plot_EKC_helper(eigvls = x$eigenvalues,
                     references = x$references$AM2019,
                     n_factors = x$n_factors_AM2019,
                     type = "AM2019")

  }

}



.plot_EKC_helper <- function(eigvls, references, n_factors, type){
  # eigen = eigenvalues found with specific type (PCA, SMC or EFA)
  # references = reference eigenvalues
  # n_factors = number of factors suggested by EKC
  # type = EKC type used

  p_eigvls <- pretty(eigvls)
  x_len <- length(eigvls)

  graphics::plot.new()
  graphics::plot.window(xlim = c(1, x_len),
                        ylim = c(min(p_eigvls),
                                 max(p_eigvls)))
  graphics::axis(1, seq_len(x_len))
  graphics::axis(2, p_eigvls, las = 1)

  graphics::mtext("Factor", side = 1, line = 3, cex = 1.5, padj =-.5)
  graphics::mtext("Eigenvalues", side = 2, line = 3, cex = 1.5, padj =.5)

  graphics::lines(seq_len(x_len), eigvls)
  graphics::points(seq_len(x_len), eigvls, pch = 16)

  graphics::lines(seq_len(x_len), references, lty = 2, lwd = 1.25,
                  col = "darkgray")

  if (!is.na(n_factors)) {
    graphics::points(n_factors, eigvls[n_factors],
                     pch = 1, cex = 2, col = "red")
    graphics::text(n_factors, eigvls[n_factors],
                   n_factors, pos = 3, cex = 1.5, col = "red",
                   font = 1, offset = .75)
  }

  factors_text <- paste0("N factors suggested by EKC (type '", type, "'): ",
                         n_factors)

  graphics::title(factors_text, cex.main = 1.3)

}

