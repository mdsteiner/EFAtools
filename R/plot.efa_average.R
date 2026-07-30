#' Plot efa_average object
#'
#' Plot method showing a summarized output of the [efa_average()] function
#'
#' @param x list. An output from the [efa_average()] function.
#' @param ... not used.
#'
#' @returns A ggplot object showing, for each indicator and factor, the minimum,
#'  maximum, and average (mean or median) loading across the averaged solutions.
#'  Each panel carries a point at the average, a bar spanning the minimum to the
#'  maximum with a tick at each endpoint, and a grey band marking the loadings
#'  that fall below the salience threshold; the caption names the four marks.
#'
#' @importFrom rlang .data
#' @export
#' @method plot efa_average
#'
#' @examples
#' \donttest{
#' EFA_aver <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500)
#' plot(EFA_aver)
#' }
#'
plot.efa_average <- function(x, ...) {

  averaging <- x$settings$averaging

  # Prepare data
  dat <- lapply(x$loadings, function(temp){
    df <- as.data.frame(unclass(temp))

    data.frame(
      rowname  = rep(rownames(df), each = ncol(df)),
      colname  = rep(colnames(df), times = nrow(df)),
      loadings = as.vector(t(as.matrix(df))),
      stringsAsFactors = FALSE
    )
  })

  dat <- do.call(cbind, dat)

  dat <- dat[, c("average.rowname", "average.colname", "average.loadings",
                 "min.loadings", "max.loadings")]
  names(dat) <- c("row_ind", "col_ind", "average", "min", "max")
  dat$row_ind <- factor(dat$row_ind, levels = rownames(x$loadings$average))
  dat$col_ind <- factor(dat$col_ind, levels = colnames(x$loadings$average))

  # Create plot faceted for variables and factors

  # No aesthetic is mapped to a scale (the four marks are drawn by four separate layers
  # with fixed appearance), so ggplot2 renders no legend; name the marks in a caption
  # instead. Without it the panel is uninterpretable: the grey band in particular is the
  # decision-relevant element and nothing else in the figure identifies it.
  centre <- if (averaging == "mean") "mean" else "median"
  caption <- paste0("Point = ", centre, "; bar = min-max range; ticks = endpoints; ",
                    "grey band = |loading| < ", format(x$settings$salience_threshold))

  plot_load <- ggplot2::ggplot(dat) +
    ggplot2::geom_segment(ggplot2::aes(x = min, xend = max, y = 0, yend = 0)) +
    ggplot2::geom_segment(ggplot2::aes(x = min, xend = min, y = -0.5, yend = 0.5)) +
    ggplot2::geom_segment(ggplot2::aes(x = max, xend = max, y = -0.5, yend = 0.5)) +
    ggplot2::geom_rect(xmin = -x$settings$salience_threshold,
                       xmax = x$settings$salience_threshold,
                       ymin = -2, ymax = 2, fill = ggplot2::alpha("grey", 0.3)) +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::geom_point(ggplot2::aes(.data$average, 0), color = "darkred") +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$row_ind),
                        cols = ggplot2::vars(.data$col_ind),
                        switch = "y") +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(paste0("Minimum, Maximum, and ",
                                ifelse(averaging == "mean", "Mean", "Median"),
                                " Loadings")) +
    # The x axis is the loading scale; the y axis is a structural dummy (every mark sits
    # at y = 0, and the vertical extent only gives the ticks and the band their height),
    # so it keeps neither a title nor labels.
    ggplot2::labs(x = "Loading", caption = caption) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.2),
          axis.ticks.x = ggplot2::element_line(color = "black", linewidth = 0.2),
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          panel.grid.minor.y = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.spacing.y = ggplot2::unit(0, "mm"),
          strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 0),
          strip.text.x = ggplot2::element_text(face = "bold"),
          strip.background.x = ggplot2::element_rect(color = "black", linewidth = 0.2),
          panel.border = ggplot2::element_rect(color = "gray", fill = NA,
                                               linewidth = 0.2)
          )

  plot_load

}
