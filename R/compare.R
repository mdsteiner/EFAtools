#' Compare two vectors or matrices (communalities or loadings)
#'
#' The functions takes two objects of the same dimensions containing numeric
#' information (loadings or communalities) and returns a list of class COMPARE
#' containing summary information of the absolute difference of the objects, as
#' well as the differences themselves.
#'
#' @param x matrix, dataframe, or vector. Loadings or communalities of one
#'  implementation.
#' @param y matrix, dataframe, or vector. Loadings or communalities of another
#'  implementation to compare to x.
#' @param reorder logical. Whether elements / columns should be reordered.
#'  (loading) matrices are reordered according to Tuckers correspondence coefficient,
#'  other objects according to the names.
#' @param digits numeric. Number of decimals to print in the output (default is 4).
#' @param m_red numeric. Number above which the mean and median should be printed
#'  in red (i.e., if .001 is used, the mean will be in red if it is larger than
#'  .001, otherwise it will be displayed in green. Default is .001).
#' @param range_red numeric. Number above which the min and max should be printed
#'  in red (i.e., if .001 is used, min and max will be in red if the max is larger
#'   than .001, otherwise it will be displayed in green. Default is .001). Note that
#'   the color of min also depends on max, that is min will be displayed in the
#'   same color as max.
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
#'  in red (default is .001).
#'
#' @return A list of class COMPARE containing summary statistics on the differences
#'  of x and y.
#'
#' \item{diff}{The vector or matrix containing the differences between x and y.}
#' \item{mean_abs_diff}{The mean absolute difference between x and y.}
#' \item{median_abs_diff}{The median absolute difference between x and y.}
#' \item{min_abs_diff}{The minimum absolute difference between x and y.}
#' \item{max_abs_diff}{The maximum absolute difference between x and y.}
#' \item{max_dec}{The maximum number of decimals to which a comparison makes sense.
#'  For example, if x contains only values up to the third decimals, and y is a
#'  normal double, max_dec will be three.}
#' \item{are_equal}{The maximal number of decimals to which x and y can be rounded
#'  and are still equal.}
#' \item{diff_corres}{The number of differing variable to factor correspondences
#'  between x and y.}
#' \item{ctrl}{List of control settings needed for the print method print.COMPARE.}
#'
#' @export
#'
#' @examples
#' # A type SPSS EFA to mimick the SPSS implementation
#' EFA_SPSS_5 <- EFA(IDS2_R, n_factors = 5, type = "SPSS")
#'
#' # A type psych EFA to mimick the psych::fa() implementation
#' EFA_psych_5 <- EFA(IDS2_R, n_factors = 5, type = "psych")
#'
#' # compare the two
#' compare(EFA_SPSS_5$unrot_loadings, EFA_psych_5$unrot_loadings,
#'         x_labels = c("SPSS", "psych"))
compare <- function(x,
                    y,
                    reorder = TRUE,
                    digits = 4,
                    m_red = .001,
                    range_red = .001,
                    round_red = 3,
                    print_diff = TRUE,
                    na.rm = FALSE,
                    x_labels = c("x", "y"),
                    plot = TRUE,
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

      stop("x and y have different lengths Compare only works with identical
           dimensions")

    }

  }

  if (isTRUE(reorder)) {

    if (class(x) == "matrix") {

      n_factors <- ncol(x)

      if (ncol(x) > 1 && ncol(y) > 1) {

        # get Tucker's conguence coefficients
        congruence <- .factor_congruence(x, y)

        # factor order for y
        factor_order <- apply(congruence, 1, which.max)

        # obtain signs to reflect signs of y if necessary
        factor_sign <- sapply(1:n_factors, function(ll, congruence, factor_order){
          sign(congruence[ll, factor_order[ll]])
        }, congruence = congruence, factor_order = factor_order)

        factor_sign <- rep(factor_sign, each = nrow(x))

        # reorder
        y <- y[, factor_order]

        # reflect signs if necessary
        y <- y * factor_sign

        # factor correspondences
        # factor correspondence x
        x_corres <- apply(x, 1, function(mm) which.max(abs(mm)))

        # factor correspondence y
        y_corres <- apply(y, 1, function(mm) which.max(abs(mm)))

        # different factor correspondences
        diff_corres <- sum(x_corres != y_corres)

      } else if (ncol(x) == 1 && ncol(y) == 1) {
        # different factor correspondences
        diff_corres <- 0
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

  if(length(are_equal) == 0){
    are_equal <- NA
  }

  ctrl <- list(
    digits = digits,
    m_red = m_red,
    range_red = range_red,
    round_red = round_red,
    print_diff = print_diff,
    x_labels = x_labels,
    plot = plot,
    plot_red = plot_red
  )

  # create output list
  out <- list(
    diff = diff,
    mean_abs_diff = mean_abs_diff,
    median_abs_diff = median_abs_diff,
    min_abs_diff = min_abs_diff,
    max_abs_diff = max_abs_diff,
    max_dec = max_dec,
    are_equal = are_equal,
    diff_corres = diff_corres,
    ctrl = ctrl
  )

  class(out) <- "COMPARE"

  out
}
