#' Plot COMPARE object
#'
#' Plot method for the [COMPARE()] function showing the distribution of the
#' absolute differences between the two compared objects as a violin plot with
#' jittered points. Differences above the threshold are highlighted.
#'
#' @param x list. An object of class COMPARE (output from the [COMPARE()] function).
#' @param ... not used.
#'
#' @returns A ggplot object showing the absolute differences, with differences
#'  above `plot_red` highlighted in red.
#'
#' @importFrom rlang .data
#' @export
#' @method plot COMPARE
#'
#' @examples
#' # A type SPSS EFA to mimick the SPSS implementation
#' EFA_SPSS_5 <- EFA(IDS2_R, n_factors = 5, type = "SPSS")
#'
#' # A type psych EFA to mimick the psych::fa() implementation
#' EFA_psych_5 <- EFA(IDS2_R, n_factors = 5, type = "psych")
#'
#' # compare the two and plot the differences
#' comp <- COMPARE(EFA_SPSS_5$unrot_loadings, EFA_psych_5$unrot_loadings,
#'                 x_labels = c("SPSS", "psych"))
#' plot(comp)
plot.COMPARE <- function(x, ...) {

  diff <- x$diff
  plot_red <- x$settings$plot_red
  x_labels <- x$settings$x_labels

  if (length(c(diff)) <= 2) {
    cli::cli_abort(
      "{.fun plot.COMPARE} needs more than two differences to plot; this comparison has {length(c(diff))}.",
      class = "efa_compare_too_few_to_plot"
    )
  }

  # prepare variable for plot
  diffs <- as.vector(abs(diff))
  diff_dat <- data.frame(
    diffs = diffs,
    color = ifelse(diffs >= plot_red, "large difference", "acceptable difference"),
    comp = paste(x_labels, collapse = " vs. "),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(diff_dat, ggplot2::aes(.data$comp, .data$diffs,
                                         col = .data$color)) +
    ggplot2::geom_violin(col = "grey20", width = .7, linewidth = .7) +
    ggplot2::geom_hline(yintercept = plot_red, lty = 2, alpha = .5,
                        linewidth = 1.25) +
    ggplot2::geom_jitter(alpha = .5, width = 0.05, height = 0, size = 2) +
    ggplot2::scale_color_manual(values = c("acceptable difference" = "black",
                                           "large difference" = "red")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      subtitle = paste("Threshold for difference coloring:", plot_red),
      x = "Compared Variables",
      y = "Absolute Difference"
    ) +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(size = 11, face = "bold"),
      axis.text = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 13, face = "bold")
    )

}
