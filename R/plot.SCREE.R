#' Plot SCREE object
#'
#' Plot method showing a summarized output of the \link{SCREE} function
#'
#' @param x a list of class SCREE An output from the \link{SCREE} function.
#' @param ... not used.
#'
#' @export
#' @method plot SCREE
#'
#' @examples
#' SCREE_base <- SCREE(test_models$baseline$cormat)
#' plot(SCREE_base)
#'
plot.SCREE <- function(x, ...) {

  eigen_PCA <- x$eigen_PCA
  eigen_SMC <- x$eigen_SMC
  eigen_EFA <- x$eigen_EFA

  # Create plots
  if(!all(is.na(eigen_PCA))){

    .plot_SCREE_helper(eigen = eigen_PCA, eigen_name = "PCA")

  }

  if (!all(is.na(eigen_SMC))) {

    .plot_SCREE_helper(eigen = eigen_SMC, eigen_name = "SMC")

  }

  if (!all(is.na(eigen_EFA))) {

    .plot_SCREE_helper(eigen = eigen_EFA, eigen_name = "EFA")

  }

}

.plot_SCREE_helper <- function(eigen, eigen_name){
  # eigen = eigenvalues found with specific type (PCA, SMC or EFA)
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

  title <- paste0("Scree plot with ", eigen_name, "-determined eigenvalues",
                  sep = "")

  graphics::title(title, cex.main = 1.2)

}
