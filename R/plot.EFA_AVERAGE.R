#' Plot EFA_AVERAGE object
#'
#' Plot method showing a summarized output of the \link{EFA_AVERAGE} function
#'
#' @param x list. An output from the \link{EFA_AVERAGE} function.
#' @param ... not used.
#'
#' @importFrom rlang .data
#' @export
#' @method plot EFA_AVERAGE
#'
#' @examples
#' EFA_aver <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500)
#' EFA_aver
#'
plot.EFA_AVERAGE <- function(x, ...) {

  averaging <- x$settings$averaging

  # Prepare data
  dat <- lapply(x$loadings, function(temp){
    temp <- as.data.frame(unclass(temp))

    temp <- temp  %>%
      tibble::rownames_to_column() %>%
      tidyr::pivot_longer(-.data$rowname, names_to = "colname", values_to = "loadings")

    return(temp)
  })

  dat <- do.call(cbind, dat)

  dat <- dat[, c("average.rowname", "average.colname", "average.loadings",
                 "min.loadings", "max.loadings")]
  names(dat) <- c("row_ind", "col_ind", "average", "min", "max")
  dat$row_ind <- factor(dat$row_ind, levels = rownames(x$loadings$average))
  dat$col_ind <- factor(dat$col_ind, levels = colnames(x$loadings$average))

  # Create plot faceted for variables and factors

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
                                " Loadings"), ) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_line(color = "black", size = 0.2),
          axis.ticks.x = ggplot2::element_line(color = "black", size = 0.2),
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          panel.grid.minor.y = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.spacing.y = ggplot2::unit(0, "mm"),
          strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 0),
          strip.text.x = ggplot2::element_text(face = "bold"),
          strip.background.x = ggplot2::element_rect(color = "black", size = 0.2),
          panel.border = ggplot2::element_rect(color = "gray", fill = NA,
                                               size = 0.2)
          )

  print(plot_load)

}
