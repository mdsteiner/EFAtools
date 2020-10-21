#' Plot EFA_AVERAGE object
#'
#' Plot method showing a summarized output of the \link{EFA_AVERAGE} function
#'
#' @param x list. An output from the \link{EFA_AVERAGE} function.
#' @param ... not used.
#'
#' @export
#' @method plot EFA_AVERAGE
#'
#' @examples
#' EFA_aver <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500)
#' EFA_aver
#'
plot.EFA_AVERAGE <- function(x, ...) {

  aggregation <- x$settings$aggregation

  # Prepare data
  dat <- lapply(x$loadings, function(x){
    x <- as.data.frame(unclass(x))

    x <- x  %>%
      tibble::rownames_to_column() %>%
      tidyr::pivot_longer(-rowname, names_to = "colname", values_to = "loadings")

    return(x)
  })

  dat <- do.call(cbind, dat)

  dat <- dat %>%
    dplyr::select(row_ind = aggregate.rowname,
           col_ind = aggregate.colname,
           aggregate = aggregate.loadings,
           min = min.loadings,
           max = max.loadings) %>%
    dplyr::mutate(row_ind = factor(row_ind, levels = rownames(x$loadings$aggregate)),
           col_ind = factor(col_ind, levels = colnames(x$loadings$aggregate)))

  # Create plot faceted for variables and factors

  plot_load <- ggplot2::ggplot(dat) +
    ggplot2::geom_segment(ggplot2::aes(x = min, xend = max, y = 0, yend = 0)) +
    ggplot2::geom_segment(ggplot2::aes(x = min, xend = min, y = -0.5, yend = 0.5)) +
    ggplot2::geom_segment(ggplot2::aes(x = max, xend = max, y = -0.5, yend = 0.5)) +
    ggplot2::geom_rect(xmin = -x$settings$salience_threshold,
                       xmax = x$settings$salience_threshold,
                       ymin = -2, ymax = 2, fill = ggplot2::alpha("grey", 0.3)) +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::geom_point(ggplot2::aes(aggregate, 0), color = "darkred") +
    ggplot2::facet_grid(rows = ggplot2::vars(row_ind),
                        cols = ggplot2::vars(col_ind),
                        switch = "y") +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(paste0("Minimum, Maximum, and ",
                                ifelse(aggregation == "mean", "Mean", "Median"),
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
