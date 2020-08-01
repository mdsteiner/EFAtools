#' Compare two vectors or matrices (communalities or loadings)
#'
#' The function takes two objects of the same dimensions containing numeric
#' information (loadings or communalities) and returns a list of class COMPARE
#' containing summary information of the differences of the objects.
#'
#' @param x matrix, or vector. Loadings or communalities of a factor
#'  analysis output.
#' @param y matrix, or vector. Loadings or communalities of another
#'  factor analysis output to compare to x.
#' @param reorder character. Whether and how elements / columns should be
#' reordered. If "congruence" (default), reordering is done according to Tuckers
#' correspondence coefficient, if "names", objects according to their names,
#' if "none", no reordering is done.
#' @param corres logical. Whether factor correspondences should be compared if a
#'  matrix is entered.
#' @param thresh numeric. The threshold to classify a pattern coefficient as substantial. Default is .3.
#' @param digits numeric. Number of decimals to print in the output. Default is 4.
#' @param m_red numeric. Number above which the mean and median should be printed
#'  in red (i.e., if .001 is used, the mean will be in red if it is larger than
#'  .001, otherwise it will be displayed in green.) Default is .001.
#' @param range_red numeric. Number above which the min and max should be printed
#'  in red (i.e., if .001 is used, min and max will be in red if the max is larger
#'   than .001, otherwise it will be displayed in green. Default is .001). Note that
#'   the color of min also depends on max, that is min will be displayed in the
#'   same color as max.
#' @param round_red  numeric. Number above which the max decimals to round to where
#' all corresponding elements of x and y are still equal are displayed in red
#' (i.e., if 3 is used, the number will be in red if it is smaller than
#'  3, otherwise it will be displayed in green). Default is 3.
#' @param print_diff logical. Whether the difference vector or matrix should be
#'  printed or not. Default is TRUE.
#' @param na.rm logical. Whether NAs should be removed in the mean, median, min,
#'  and max functions. Default is FALSE.
#' @param x_labels character. A vector of length two containing identifying
#'  labels for the two objects x and y that will be compared. These will be used
#'  as labels on the x-axis of the plot. Default is "x" and "y".
#' @param plot logical. If TRUE (default), a plot illustrating the differences
#'  will be shown.
#' @param plot_red numeric. Threshold above which to plot the absolute differences
#'  in red. Default is .001.
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
#' \item{are_equal}{The maximal number of decimals to which all elements of x and y
#'  are equal.}
#' \item{diff_corres}{The number of differing variable-to-factor correspondences
#'  between x and y, when only the highest loading is considered.}
#' \item{diff_corres_cross}{The number of differing variable-to-factor correspondences
#'  between x and y when all loadings \code{>= thresh} are considered.}
#' \item{g}{The root mean squared distance (RMSE) between x and y.}
#' \item{settings}{List of the settings used.}
#'
#' @export
#'
#' @examples
#' # A type SPSS EFA to mimick the SPSS implementation
#' EFA_SPSS_6 <- EFA(test_models$case_11b$cormat, n_factors = 6, type = "SPSS")
#'
#' # A type psych EFA to mimick the psych::fa() implementation
#' EFA_psych_6 <- EFA(test_models$case_11b$cormat, n_factors = 6, type = "psych")
#'
#' # compare the two
#' COMPARE(EFA_SPSS_6$unrot_loadings, EFA_psych_6$unrot_loadings,
#'         x_labels = c("SPSS", "psych"))
COMPARE <- function(x,
                    y,
                    reorder = c("congruence", "names", "none"),
                    corres = TRUE,
                    thresh = .3,
                    digits = 4,
                    m_red = .001,
                    range_red = .001,
                    round_red = 3,
                    print_diff = TRUE,
                    na.rm = FALSE,
                    x_labels = c("x", "y"),
                    plot = TRUE,
                    plot_red = .01)  {

  reorder <- match.arg(reorder)
  checkmate::assert_flag(corres)
  checkmate::assert_number(thresh)
  checkmate::assert_count(digits)
  checkmate::assert_number(m_red)
  checkmate::assert_number(range_red)
  checkmate::assert_number(round_red)
  checkmate::assert_flag(print_diff)
  checkmate::assert_flag(na.rm)
  checkmate::assert_character(x_labels, len = 2)
  checkmate::assert_flag(plot)
  checkmate::assert_number(plot_red)

  # reclass data.frames and tibbles to matrices so the stats functions afterwards
  # work
  if ((inherits(x, c("loadings", "LOADINGS", "SLLOADINGS", "matrix"))) &&
      (inherits(y, c("loadings", "LOADINGS", "SLLOADINGS", "matrix")))) {

    if (inherits(x, c("loadings", "LOADINGS", "SLLOADINGS"))) {
      x <- unclass(x)
    }

    if (inherits(y, c("loadings", "LOADINGS", "SLLOADINGS"))) {
      y <- unclass(y)
    }

    # check if dimensions match:
    if (any(dim(x) != dim(y))) {

      stop(crayon::red$bold(cli::symbol$circle_cross),
           crayon::red(" 'x' and 'y' have different dimensions. Can only compare matrices with identical dimensions.\n"))

    }

  } else if (inherits(x, c("numeric", "integer")) &&
             inherits(y, c("numeric", "integer"))) {

    if (length(x) != length(y)) {

      stop(crayon::red$bold(cli::symbol$circle_cross),
           crayon::red(" 'x' and 'y' have different lengths Compare only works with identical dimensions.\n"))

    }

  } else {

    stop(crayon::red$bold(cli::symbol$circle_cross),
         crayon::red(" 'x' is of class", class(x), "and 'y' is of class",
                     class(y), "but must be numeric vectors or matrices\n"))

  }

  if (inherits(x, "matrix")) {

    n_factors <- ncol(x)

    if (reorder == "congruence" && n_factors > 1) {
      # get Tucker's congruence coefficients
      congruence <- .factor_congruence(x, y)

      if (any(is.na(congruence))) {
        stop(crayon::red$bold(cli::symbol$circle_cross),
             crayon::red(" Tucker's congruence coefficients contained NAs, cannot reorder columns based on congruence. Try another reordering method.\n"))
      }

      # factor order for y
      factor_order <- apply(abs(congruence), 1, which.max)

      # obtain signs to reflect signs of y if necessary
      factor_sign <- sapply(seq_len(n_factors),
                            function(ll, congruence, factor_order){
                              sign(congruence[ll, factor_order[ll]])
                              }, congruence = congruence,
                            factor_order = factor_order)

      factor_sign <- rep(factor_sign, each = nrow(x))

      # reorder
      y <- y[, factor_order]

      # reflect signs if necessary
      y <- y * factor_sign
    } else if (reorder == "names" && n_factors > 1) {

      if (!is.null(colnames(x)) && !is.null(colnames(y))) {
        x <- x[, order(colnames(x))]
        y <- y[, order(colnames(y))]

        if(!all(colnames(x) == colnames(y))) {
          warning(crayon::yellow$bold("!"), crayon::yellow(" reorder = 'names' was used but colnames of x and y were not identical. Results might be inaccurate.\n"))
        }
      } else if (is.null(colnames(x)) || is.null(colnames(y))) {
        warning(crayon::yellow$bold("!"), crayon::yellow(" reorder was set to 'names' but at least one of 'x' and 'y' was not named. Proceeding without reordering.\n"))
      }

    }

    if (n_factors > 1 && isTRUE(corres)) {
      # factor correspondences
      corres_list <- .factor_corres(x, y, thresh = thresh)
      diff_corres <- corres_list$diff_corres
      diff_corres_cross <- corres_list$diff_corres_cross
    } else {
      diff_corres <- 0
      diff_corres_cross <- 0
    }


  } else if (inherits(x, c("numeric", "integer"))) {

    if (reorder == "congruence" && !is.null(names(x)) && !is.null(names(y))){

      warning(crayon::yellow$bold("!"), crayon::yellow(" reorder was set to 'congruence', but this only works for matrices. To reorder vectors, set reorder = 'names'. Proceeding without reordering.\n"))

    } else if (reorder == "names") {

      if (!is.null(names(x)) && !is.null(names(y))) {

        x <- x[order(names(x))]
        y <- y[order(names(y))]

        if (!all(names(x) == names(y))) {
          warning(crayon::yellow$bold("!"), crayon::yellow(" reorder = 'names' was used but names of x and y were not identical. Results might be inaccurate.\n"))
        }



      } else if (is.null(names(x)) || is.null(names(y))) {
        warning(crayon::yellow$bold("!"), crayon::yellow(" reorder was set to 'names' but at least one of 'x' and 'y' was not named. Proceeding without reordering.\n"))
      }

      }

      diff_corres <- NA
      diff_corres_cross <- NA

      g <- NA

  }

  # compute differences and statistics
  diff <- x - y

  if(inherits(x, "matrix")) {
    g <- sqrt(sum(diag(t(diff) %*% (diff))) / prod(dim(x)))
  } else {
    g <- sqrt((sum(diff ** 2) / length(diff)))
  }


  mean_abs_diff <- mean(abs(diff), na.rm = na.rm)
  median_abs_diff <- stats::median(abs(diff), na.rm = na.rm)

  min_abs_diff <- min(abs(diff), na.rm = na.rm)
  max_abs_diff <- max(abs(diff), na.rm = na.rm)

  are_equal_v <- c()

  max_dec <- min(c(.decimals(x), .decimals(y)))

  class(x) <- "character"
  class(y) <- "character"
  x <- gsub("-|\\.", "", x)
  y <- gsub("-|\\.", "", y)
  for (ii in seq_len((max_dec + 1))) {
    are_equal_v[ii] <- all(substr(x, 1, ii) == substr(y, 1, ii))
  }

  are_equal <- utils::tail(which(are_equal_v), 1)

  if(length(are_equal) == 0){
    are_equal <- 0
  } else {
     are_equal <- are_equal - 1
  }

  settings <- list(
    reorder = reorder,
    corres = corres,
    digits = digits,
    thresh = thresh,
    m_red = m_red,
    range_red = range_red,
    round_red = round_red,
    print_diff = print_diff,
    na.rm = na.rm,
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
    diff_corres_cross = diff_corres_cross,
    g = g,
    settings = settings
  )

  class(out) <- "COMPARE"

  out
}
