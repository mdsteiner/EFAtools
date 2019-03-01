#' Compare loadings or communalities of different implementations
#'
#' The functions takes two objects of the same dimensions containing numeric
#' information (loadings or communalities) and prints the following summary
#' information of the absolute difference of the objects: mean, median, min, max,
#' number above which the max decimals to round to where
#' all corresponding elements of x and y are still equal are displayed, and the
#' difference object itself.
#'
#' @param x matrix, dataframe, or vector. Loadings or communalities of one
#'  implementation.
#' @param y matrix, dataframe, or vector. Loadings or communalities of another
#'  implementation to compare to x.
#' @param reorder logical. Whether factors should be reordered according to their
#'   colnames.
#' @param digits numeric. Number of decimals to print in the output
#'  (default is 4).
#' @param m_red numeric. Number above which the mean and median should be printed
#'  in red (i.e., if .001 is used, the mean will be in red if it is larger than
#'  .001, otherwise it will be displayed in green. Default is .001).
#' @param range_red numeric. Number above which the min and max should be printed
#'  in red (i.e., if .001 is used, min and max will be in red if the max is larger
#'   than .001, otherwise it will be displayed in green. Default is .001). Note that
#'   the color of min also depends on max.
#' @param round_red  numeric. Number above which the max decimals to round to where
#' all corresponding elements of x and y are still equal are displayed in red
#' (i.e., if 3 is used, the number will be in red if it is smaller than
#'  3, otherwise it will be displayed in green. Default is 3).
#' @param print_diff logical. Whether the difference vector or matrix should be
#'  printed or not (default is TRUE).
#' @param na.rm logical. Whether NAs should be removed in the mean, median, min,
#'  and max functions. Default is FALSE.
#' @param x_labels character. A vector of length two containing identifying
#'  labels for the two objects x and y that will be compared. These will be used
#'  as labels on the x-axis of the plot. If left to NULL, "Var 1" and "Var 2" will
#'  be used.
#' @param plot logical. If TRUE (default), a plot illustrating the differences
#'  will be shown.
#' @param plot_red numeric. Threshold above which to plot the absolute differences
#'  in red.
#'
#' @return Print out a comparison summary.
#' @export
compare <- function(x, y, reorder = TRUE, digits = 4, m_red = .001,
                    range_red = .001, round_red = 3, print_diff = TRUE,
                    na.rm = FALSE, x_labels = c("x", "y"), plot = TRUE,
                    plot_red = .001)  {


  # reclass data.frames and tibbles to matrices so the stats functions afterwards
  # work
  if (!(class(x) %in% c("numeric", "COMMUNALITIES", "EIGEN")) &&
      !(class(y) %in% c("numeric", "COMMUNALITIES", "EIGEN"))) {

    if (any(class(x) == "data.frame")) {
      x <- as.matrix(x)
    } else if (any(class(x) %in% c("loadings", "LOADINGS", "SLLOADINGS"))) {
      x <- unclass(x)
    }

    if (any(class(y) == "data.frame")) {
      y <- as.matrix(y)
    } else if (any(class(y) %in% c("loadings", "LOADINGS", "SLLOADINGS"))) {
      y <- unclass(y)
    }

    # check if dimensions match:
    if (any(dim(x) != dim(y))) {

      stop("x and y have different dimensions. Compare only works with identical dimensions")

    }

  } else if (class(x) %in% c("numeric", "COMMUNALITIES", "EIGEN") &&
             class(y) %in% c("numeric", "COMMUNALITIES", "EIGEN")) {

    x <- unclass(x)
    y <- unclass(y)

    if (length(x) != length(y)) {

      stop("x and y have different lengths Compare only works with identical dimensions")

    }

  }

  if (isTRUE(reorder)) {

    if (class(x) == "matrix") {

      if (!is.null(colnames(x)) && !is.null(colnames(y))) {
        ind_x <- order(colnames(x))
        x <- x[, ind_x]

        ind_y <- order(colnames(y))
        y <- y[, ind_y]
      }


    } else if (class(x) == "numeric") {

      if (!is.null(names(x)) && !is.null(names(y))) {
        ind_x <- order(names(x))
        x <- x[ind_x]

        ind_y <- order(names(y))
        y <- y[ind_y]
      }

    }


  }

  # compute differences and statistics
  diff <- x - y
  mean_abs_diff <- round(mean(abs(diff), na.rm = na.rm), digits = digits)
  median_abs_diff <- round(stats::median(abs(diff), na.rm = na.rm),
                           digits = digits)

  min_abs_diff <- round(min(abs(diff), na.rm = na.rm), digits = digits)
  max_abs_diff <- round(max(abs(diff), na.rm = na.rm), digits = digits)

  are_equal_v <- c()

  max_dec <- min(c(.decimals(x), .decimals(y)))
  for (ii in 1:max_dec) {
    are_equal_v[ii] <- all(round(x, digits = ii) == round(y, digits = ii))
  }

  are_equal <- utils::tail(which(are_equal_v), 1)

  # prepare to print statistics

  if (mean_abs_diff <= m_red) {
    mean_out <- crayon::green$bold(.numformat(mean_abs_diff, digits,
                                              TRUE))
  } else {
    mean_out <- crayon::red$bold(.numformat(mean_abs_diff, digits, TRUE))
  }

  if (median_abs_diff <= m_red) {
    median_out <- crayon::green$bold(.numformat(median_abs_diff, digits,
                                                TRUE))
  } else {
    median_out <- crayon::red$bold(.numformat(median_abs_diff, digits,
                                              TRUE))
  }

  if (max_abs_diff <= range_red) {
    max_out <- crayon::green$bold(.numformat(max_abs_diff, digits, TRUE))
    min_out <- crayon::green$bold(.numformat(min_abs_diff, digits, TRUE))
  } else {
    max_out <- crayon::red$bold(.numformat(max_abs_diff, digits, TRUE))
    min_out <- crayon::red$bold(.numformat(min_abs_diff, digits, TRUE))
  }

  if (length(are_equal) == 0) {
    equal_out <- crayon::red$bold("none")
  } else if (are_equal < round_red) {
    equal_out <- crayon::red$bold(are_equal)
  } else {
    equal_out <- crayon::green$bold(are_equal)
  }



  cat("Mean [min, max] absolute difference: ")
  cat(paste0(mean_out, " [", min_out, ", ", max_out, "]"))
  cat("\n")
  cat(paste0("Median absolute difference: ", median_out))
  cat("\n")
  cat(paste0("Max decimals to round to where numbers are still equal: ",
             equal_out))
  cat("\n")
  cat(paste0("Minimum number of decimals provided: ", crayon::blue$bold(max_dec)))

  # create the difference object
  if (isTRUE(print_diff)) {
    cat("\n")
    cat("\n")
    if (class(diff) == "matrix") {

      out_diff <- .get_compare_matrix(diff, digits = digits, r_red = range_red)

    } else {

      out_diff <- .get_compare_vector(diff, digits = digits, r_red = range_red)

    }

  }
  cat(out_diff)

  if (isTRUE(plot)) {

    if (length(x_labels) < 2) {
      warning("Less than two x_labels specified, using 'x' and 'y'.")
      x_labels <- c("x", "y")
    } else if (length(x_labels) > 2) {
      warning("More than two x_labels specified, using only the first two.")
      x_labels <- x_labels[1:2]
    }

    # prepare variable for plot
    diff_dat <- tibble::tibble(diffs = as.vector(abs(diff))) %>%
      dplyr::mutate(color = dplyr::case_when(diffs >= plot_red ~ "large difference",
                                             TRUE ~ "acceptable difference"),
                    comp = paste(x_labels, collapse = " vs. "))

    diff_plot <- ggplot2::ggplot(diff_dat, ggplot2::aes_string("comp", "diffs",
                                                               col = "color")) +
      ggplot2::geom_violin(col = "grey20", width = .7, size = .7) +
      ggplot2::geom_hline(yintercept = plot_red, lty = 2, alpha = .5,
                          size = 1.25) +
      ggplot2::geom_jitter(alpha = .5, width = 0.05, height = 0, size = 2) +
      ggplot2::scale_color_manual(values = c("black", "red")) +
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
        axis.title = ggplot2::element_text(size = 13,face = "bold")
      )

      diff_plot


  }


}
