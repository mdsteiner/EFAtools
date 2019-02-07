#' Compare loadings of different implementations
#'
#' @param x matrix or dataframe. Loadings of one implementation.
#' @param y matrix or dataframe. Loadings of another implementation to compare to
#' x.
#' @param reorder logical. Whether factors should be reordered according to their
#'   colnames.
#' @param print_digits numeric. Number of decimals to print in the output
#'  (default is 4).
#'
#' @return Print out a comparison summary.
#' @export
compare <- function(x, y, reorder = TRUE, print_digits = 4) {

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

  are_equal_v <- c()
  for (ii in 1:20) {
    are_equal_v[ii] <- all(round(x, digits = ii) == round(y, digits = ii))
  }

  are_equal <- utils::tail(which(are_equal_v), 1)

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

  if (length(are_equal) == 0) {
    equal_out <- crayon::red$bold("none")
  } else if (are_equal < 3) {
    equal_out <- crayon::red$bold(are_equal)
  } else {
    equal_out <- crayon::green$bold(are_equal)
  }

  cat("Mean [min, max] absolute difference: ")
  cat(paste0(mean_out, " [", min_out, ", ", max_out, "]"))
  cat("\n")
  cat(paste0("Median absolute difference: ", median_out))
  cat("\n")
  cat(paste0("Max decimals to round to where loadings are equal: ",
             equal_out))

}
