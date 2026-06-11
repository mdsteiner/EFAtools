#' Print COMPARE object
#'
#' Print Method showing a summarized output of the [COMPARE()] function.
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

  # prepare to print statistics

  if (mean_abs_diff <= m_red) {
    mean_out <- .efa_style(.numformat(mean_abs_diff, digits, TRUE), c("green", "bold"))
  } else {
    mean_out <- .efa_style(.numformat(mean_abs_diff, digits, TRUE), c("red", "bold"))
  }

  if (median_abs_diff <= m_red) {
    median_out <- .efa_style(.numformat(median_abs_diff, digits, TRUE), c("green", "bold"))
  } else {
    median_out <- .efa_style(.numformat(median_abs_diff, digits, TRUE), c("red", "bold"))
  }

  if (max_abs_diff <= range_red) {
    max_out <- .efa_style(.numformat(max_abs_diff, digits, TRUE), c("green", "bold"))
    min_out <- .efa_style(.numformat(min_abs_diff, digits, TRUE), c("green", "bold"))
  } else {
    max_out <- .efa_style(.numformat(max_abs_diff, digits, TRUE), c("red", "bold"))
    min_out <- .efa_style(.numformat(min_abs_diff, digits, TRUE), c("red", "bold"))
  }

  if (is.na(are_equal)) {
    equal_out <- .efa_style("none", c("red", "bold"))
  } else if (are_equal < round_red) {
    equal_out <- .efa_style(are_equal, c("red", "bold"))
  } else {
    equal_out <- .efa_style(are_equal, c("green", "bold"))
  }



  cat("Mean [min, max] absolute difference: ")
  cat(paste0(mean_out, " [", min_out, ", ", max_out, "]"))
  cat("\n")
  cat(paste0("Median absolute difference: ", median_out))
  cat("\n")
  cat(paste0("Max decimals where all numbers are equal: ",
             equal_out))
  cat("\n")
  cat(paste0("Minimum number of decimals provided: ", .efa_style(max_dec, "bold")))

  # create the difference object
  if (isTRUE(print_diff)) {
    cat("\n")
    cat("\n")
    .print_compare_diff(diff, digits = digits, r_red = range_red)
  }

}

# Render a COMPARE difference object through the shared matrix renderer, colouring cells
# whose absolute difference exceeds `r_red` (the COMPARE `range_red` threshold). A matrix is
# shown with its variable rows and factor columns; a vector becomes a single column (one
# value per row), with no column header.
.print_compare_diff <- function(diff, digits, r_red) {
  is_vector <- !inherits(diff, "matrix")
  if (is_vector) {
    diff <- matrix(diff, ncol = 1, dimnames = list(names(diff), NULL))
  }

  .print_efa_matrix(diff, role = "compare", digits = digits, cutoff = r_red,
                    header = !is_vector)
}
