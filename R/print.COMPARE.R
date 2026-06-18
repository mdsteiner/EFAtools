#' Print and format a COMPARE object
#'
#' `print()` shows a summarised output of the [COMPARE()] function: the mean
#' (with its range), median, the number of decimals to which all numbers agree,
#' and the minimum number of decimals provided, followed (optionally) by the
#' table of elementwise differences. `format()` assembles the same report and
#' returns it as a character vector; `print()` is `cat(format(x), sep = "\n")`.
#' The lines follow the active console theme, so they are plain when colours are
#' disabled (for example when captured into a file or stripped with
#' [cli::ansi_strip()]).
#'
#' @param x An object of class `COMPARE` (output from [COMPARE()]).
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines (styled to the active console theme;
#'   plain when colours are disabled).
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
#' comp <- COMPARE(EFA_SPSS_5$unrot_loadings, EFA_psych_5$unrot_loadings,
#'                 x_labels = c("SPSS", "psych"))
#' comp
#'
#' # format() returns the same lines as plain text:
#' writeLines(format(comp))
#'
print.COMPARE <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.COMPARE
#' @export
#' @method format COMPARE
format.COMPARE <- function(x, ...) {

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

  # Style each statistic green when it clears its reduction threshold and red otherwise
  # (so a smaller difference reads as a closer match). The values keep their decimal
  # padding, so the lines are emitted verbatim below.
  if (mean_abs_diff <= m_red) {
    mean_out <- .efa_style(.efa_num(mean_abs_diff, digits, TRUE), c("green", "bold"))
  } else {
    mean_out <- .efa_style(.efa_num(mean_abs_diff, digits, TRUE), c("red", "bold"))
  }

  if (median_abs_diff <= m_red) {
    median_out <- .efa_style(.efa_num(median_abs_diff, digits, TRUE), c("green", "bold"))
  } else {
    median_out <- .efa_style(.efa_num(median_abs_diff, digits, TRUE), c("red", "bold"))
  }

  if (max_abs_diff <= range_red) {
    max_out <- .efa_style(.efa_num(max_abs_diff, digits, TRUE), c("green", "bold"))
    min_out <- .efa_style(.efa_num(min_abs_diff, digits, TRUE), c("green", "bold"))
  } else {
    max_out <- .efa_style(.efa_num(max_abs_diff, digits, TRUE), c("red", "bold"))
    min_out <- .efa_style(.efa_num(min_abs_diff, digits, TRUE), c("red", "bold"))
  }

  if (is.na(are_equal)) {
    equal_out <- .efa_style("none", c("red", "bold"))
  } else if (are_equal < round_red) {
    equal_out <- .efa_style(are_equal, c("red", "bold"))
  } else {
    equal_out <- .efa_style(are_equal, c("green", "bold"))
  }

  cli::cli_format_method({
    # Fixed-format statistic lines: the values carry decimal padding and conditional
    # green/red styling, so emit them verbatim (cli never reflows them mid-token).
    cli::cli_verbatim(paste0("Mean [min, max] absolute difference: ",
                             mean_out, " [", min_out, ", ", max_out, "]"))
    cli::cli_verbatim(paste0("Median absolute difference: ", median_out))
    cli::cli_verbatim(paste0("Max decimals where all numbers are equal: ", equal_out))
    cli::cli_verbatim(paste0("Minimum number of decimals provided: ",
                             .efa_style(max_dec, "bold")))

    if (isTRUE(print_diff)) {
      cli::cli_text("")
      .efa_emit_lines(.compare_diff_lines(diff, digits = digits, r_red = range_red))
    }
  })
}

# Render a COMPARE difference object through the shared matrix renderer into verbatim lines,
# colouring cells whose absolute difference exceeds `r_red` (the COMPARE `range_red`
# threshold) red. A matrix is shown with its variable rows and factor columns; a vector
# becomes a single unlabelled column (one value per row), so its column-header row is dropped.
.compare_diff_lines <- function(diff, digits, r_red) {
  is_vector <- !inherits(diff, "matrix")
  if (is_vector) {
    diff <- matrix(diff, ncol = 1, dimnames = list(names(diff), NULL))
  }

  lines <- .efa_format_matrix(
    values = diff,
    row_labels = .efa_variable_names(diff),
    col_labels = .efa_factor_names(diff),
    col_roles = rep("compare", ncol(diff)),
    cutoff = r_red,
    digits = digits
  )

  if (is_vector) lines[-1L] else lines
}
