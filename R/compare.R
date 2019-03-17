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
#' @return A list of class COMPARE containing summary statistics on the differences
#'  of x and y.
#'
#' \item{diff}{The vector or matrix containing the }
#' \item{mean_abs_diff}{Initial communality estimates from PAF.}
#' \item{median_abs_diff}{Final communality estimates from unrotated loadings.}
#' \item{min_abs_diff}{The number of iterations needed for convergence in PAF.}
#' \item{max_abs_diff}{Eigen values of the original correlation matrix.}
#' \item{max_dec}{Initial eigenvalues, obtained from the correlation matrix
#'  with the initial communality estimates as diagonal in PAF.}
#' \item{are_equal}{Eigenvalues of the final iteration in PAF.}
#' \item{ctrl}{List of control settings used in the print method.}
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
    ctrl = ctrl
  )

  class(out) <- "COMPARE"

  out
}
