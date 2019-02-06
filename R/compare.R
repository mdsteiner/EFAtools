#' Compare loadings of different implementations
#'
#' @param x matrix or dataframe. Loadings of one implementation.
#' @param y matrix or dataframe. Loadings of another implementation to compare to
#' x.
#' @param digits numeric. Number of decimals to round loadings to and compare
#'   whether they are all equal.
#' @param reorder logical. Whether factors should be reordered according to their
#'   colnames.
#'
#' @return Print out a comparison summary.
#' @export
compare <- function(x, y, round_digits = 4, reorder = TRUE, print_digits = 4) {

  if (isTRUE(reorder) && !(class(x) %in% c("numeric", "COMMUNALITIES")) &&
      !(class(y) %in% c("numeric", "COMMUNALITIES"))) {

    if (any(class(x) == "data.frame")) {
      x <- as.matrix(x)
    }

    if (any(class(y) == "data.frame")) {
      y <- as.matrix(y)
    }

    ind_x <- order(colnames(x))
    x <- x[, ind_x]
    ind_y <- order(colnames(y))
    y <- y[, ind_y]

  }

  # compute differences and statistics
  diff <- x - y
  mean_abs_diff <- mean(abs(diff))
  median_abs_diff <- stats::median(abs(diff))

  min_abs_diff <- min(abs(diff))
  max_abs_diff <- max(abs(diff))

  are_equal <- all(round(x, digits = round_digits) == round(y, digits = round_digits))

  # prepare to print statistics

  if (mean_abs_diff <= .001) {
    mean_out <- crayon::green$bold(.numformat(mean_abs_diff, print_digits,
                                              TRUE))
  } else {
    mean_out <- crayon::red$bold(.numformat(mean_abs_diff, print_digits, TRUE))
  }

  if (median_abs_diff <= .001) {
    median_out <- crayon::green$bold(.numformat(median_abs_diff, print_digits,
                                                TRUE))
  } else {
    median_out <- crayon::red$bold(.numformat(median_abs_diff, print_digits,
                                              TRUE))
  }

  if (min_abs_diff <= .00001) {
    min_out <- crayon::green$bold(.numformat(min_abs_diff, print_digits, TRUE))
  } else {
    min_out <- crayon::red$bold(.numformat(min_abs_diff, print_digits, TRUE))
  }

  if (max_abs_diff <= .001) {
    max_out <- crayon::green$bold(.numformat(max_abs_diff, print_digits, TRUE))
  } else {
    max_out <- crayon::red$bold(.numformat(max_abs_diff, print_digits, TRUE))
  }

  if (isTRUE(are_equal)) {
    equal_out <- crayon::green$bold("All equal")
  } else {
    equal_out <- crayon::red$bold("Not all equal")
  }


  cat("Mean [min, max] absolute difference: ")
  cat(paste0(mean_out, " [", min_out, ", ", max_out, "]"))
  cat("\n")
  cat(paste0("Median absolute difference: ", median_out))
  cat("\n")
  cat(paste0("Loadings equal when rounded to ", round_digits, " decimals: ",
             equal_out))

}
