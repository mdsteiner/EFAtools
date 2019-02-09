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
#' @param print_digits numeric. Number of decimals to print in the output
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
#'
#' @return Print out a comparison summary.
#' @export
compare <- function(x, y, reorder = TRUE, print_digits = 4, m_red = .001,
                    range_red = .001, round_red = 3, print_diff = TRUE,
                    na.rm = FALSE)  {


  # reclass data.frames and tibbles to matrices so the stats functions afterwards
  # work
  if (!(class(x) %in% c("numeric", "COMMUNALITIES", "EIGEN")) &&
      !(class(y) %in% c("numeric", "COMMUNALITIES", "EIGEN"))) {

    if (any(class(x) == "data.frame")) {
      x <- as.matrix(x)
    }

    if (any(class(y) == "data.frame")) {
      y <- as.matrix(y)
    }

    # check if dimensions match:
    if (dim(x) != dim(y)) {

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

      if (!is.null(colnames(x))) {
        ind_x <- order(colnames(x))
        x <- x[, ind_x]
      }

      if(!is.null(colnames(y))) {
        ind_y <- order(colnames(y))
        y <- y[, ind_y]
      }


    } else if (class(x) == "numeric") {

      if (!is.null(names(x))) {
        ind_x <- order(names(x))
        x <- x[, ind_x]
      }

      if (!is.null(names(y))) {
        ind_y <- order(names(y))
        y <- y[, ind_y]
      }

    }


  }

  # compute differences and statistics
  diff <- x - y
  mean_abs_diff <- mean(abs(diff), na.rm = na.rm)
  median_abs_diff <- stats::median(abs(diff), na.rm = na.rm)

  min_abs_diff <- min(abs(diff), na.rm = na.rm)
  max_abs_diff <- max(abs(diff), na.rm = na.rm)

  are_equal_v <- c()
  for (ii in 1:20) {
    are_equal_v[ii] <- all(round(x, digits = ii) == round(y, digits = ii))
  }

  are_equal <- utils::tail(which(are_equal_v), 1)

  # prepare to print statistics

  if (mean_abs_diff <= m_red) {
    mean_out <- crayon::green$bold(.numformat(mean_abs_diff, print_digits,
                                              TRUE))
  } else {
    mean_out <- crayon::red$bold(.numformat(mean_abs_diff, print_digits, TRUE))
  }

  if (median_abs_diff <= m_red) {
    median_out <- crayon::green$bold(.numformat(median_abs_diff, print_digits,
                                                TRUE))
  } else {
    median_out <- crayon::red$bold(.numformat(median_abs_diff, print_digits,
                                              TRUE))
  }

  if (max_abs_diff <= range_red) {
    max_out <- crayon::green$bold(.numformat(max_abs_diff, print_digits, TRUE))
    min_out <- crayon::green$bold(.numformat(min_abs_diff, print_digits, TRUE))
  } else {
    max_out <- crayon::red$bold(.numformat(max_abs_diff, print_digits, TRUE))
    min_out <- crayon::red$bold(.numformat(min_abs_diff, print_digits, TRUE))
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
  cat(paste0("Max decimals to round to where loadings are equal: ",
             equal_out))

  # create the difference object
  if (isTRUE(print_diff)) {
    cat("\n")
    cat("\n")
    if (any(class(x) == "LOADINGS") || any(class(y) == "LOADINGS")) {
      class(diff) <- "LOADINGS"
      print.LOADINGS(diff, digits = print_digits)
    } else if (any(class(x) == "SLLOADINGS") || any(class(y) == "SLLOADINGS")) {
      class(diff) <- "SLLOADINGS"
      print.SLLOADINGS(diff, digits = print_digits)
    } else {
      print(round(diff, digits = print_digits))
    }

  }

}
