#' Print COMPARE object
#'
#' Print Method showing a summarized output of the \code{\link{COMPARE}} function.
#'
#' @param x  list. An object of class COMPARE to be printed
#' @param ... Further arguments for print.
#'
#' @export
#' @method print COMPARE
#' @examples
#' # A type SPSS EFA to mimick the SPSS implementation
#' EFA_SPSS_5 <- EFA(IDS2_R, n_factors = 5, type = "SPSS")
#'
#' # A type psych EFA to mimick the psych::fa() implementation
#' EFA_psych_5 <- EFA(IDS2_R, n_factors = 5, type = "psych")
#'
#' # compare the two
#' COMPARE(EFA_SPSS_5$unrot_loadings, EFA_psych_5$unrot_loadings,
#'         x_labels = c("SPSS", "psych"))
print.COMPARE <- function(x, ...) {

  # extract summary statistics
  diff <- x$diff
  mean_abs_diff <- x$mean_abs_diff
  median_abs_diff <- x$median_abs_diff
  min_abs_diff <- x$min_abs_diff
  max_abs_diff <- x$max_abs_diff
  max_dec <- x$max_dec
  are_equal <- x$are_equal

  # extract control settings
  digits <- x$settings$digits
  m_red <- x$settings$m_red
  range_red <- x$settings$range_red
  round_red <- x$settings$round_red
  print_diff <- x$settings$print_diff
  x_labels <- x$settings$x_labels
  plot <- x$settings$plot
  plot_red <- x$settings$plot_red

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

  if (is.na(are_equal)) {
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
  cat(paste0("Max decimals where all numbers are equal: ",
             equal_out))
  cat("\n")
  cat(paste0("Minimum number of decimals provided: ", crayon::blue$bold(max_dec)))

  # create the difference object
  if (isTRUE(print_diff)) {
    cat("\n")
    cat("\n")
    if (inherits(diff, "matrix")) {

      out_diff <- .get_compare_matrix(diff, digits = digits, r_red = range_red)

    } else {

      out_diff <- .get_compare_vector(diff, digits = digits, r_red = range_red)

    }

    cat(out_diff)

  }


  if (isTRUE(plot) && length(c(diff)) > 2) {

    if (length(x_labels) < 2) {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Less than two x_labels specified, using 'x' and 'y'."))
      x_labels <- c("x", "y")
    } else if (length(x_labels) > 2) {
      warning(crayon::yellow$bold("!"), crayon::yellow(" More than two x_labels specified, using only the first two."))
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

    print(diff_plot)

  }


}
