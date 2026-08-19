#' Plot efa_compare object
#'
#' Plot method for the [efa_compare()] function showing the distribution of the
#' absolute differences between the two compared objects as a violin plot with
#' jittered points. Differences above the threshold are highlighted.
#'
#' @param x list. An object of class `efa_compare` (output from the
#'  [efa_compare()] function).
#' @param plot_red numeric or `NULL`. Threshold above which to draw the absolute
#'  differences in red, documented in [efa_compare()]. `NULL` (default) uses the
#'  value recorded in `x$settings`; supplying one overrides it for this plot only,
#'  so the comparison need not be recomputed to redraw it at another threshold.
#' @param ... not used.
#'
#' @returns A ggplot object showing the absolute differences, with differences
#'  above `plot_red` highlighted in red.
#'
#' @importFrom rlang .data
#' @export
#' @method plot efa_compare
#'
#' @examples
#' # A type SPSS EFA to mimick the SPSS implementation
#' EFA_SPSS_5 <- efa_fit(IDS2_R, n_factors = 5,
#'                       estimate_control = estimate_control(type = "SPSS"),
#'                       rotate_control = rotate_control(type = "SPSS"))
#'
#' # A type psych EFA to mimick the psych::fa() implementation
#' EFA_psych_5 <- efa_fit(IDS2_R, n_factors = 5,
#'                        estimate_control = estimate_control(type = "psych"),
#'                        rotate_control = rotate_control(type = "psych"))
#'
#' # compare the two and plot the differences
#' comp <- efa_compare(EFA_SPSS_5$unrot_loadings, EFA_psych_5$unrot_loadings,
#'                     x_labels = c("SPSS", "psych"))
#' plot(comp)
plot.efa_compare <- function(x, plot_red = NULL, ...) {

  diff <- x$diff
  # As for the printed report's display controls, the recorded threshold is a drawing
  # setting: an argument supplied here overrides it for this plot, NULL falls back to it.
  if (is.null(plot_red)) {
    plot_red <- x$settings$plot_red
  } else {
    checkmate::assert_number(plot_red)
  }
  x_labels <- x$settings$x_labels

  if (length(c(diff)) <= 2) {
    cli::cli_abort(
      "{.fun plot.efa_compare} needs more than two differences to plot; this comparison has {length(c(diff))}.",
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
    ggplot2::labs(
      subtitle = paste("Threshold for difference coloring:", plot_red),
      x = "Compared Variables",
      y = "Absolute Difference"
    ) +
    .gg_theme() +
    # The subtitle gives the threshold that decides the two point colours, so a legend
    # adds nothing.
    ggplot2::theme(legend.position = "none")

}
