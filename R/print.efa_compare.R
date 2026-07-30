#' Print and format an efa_compare object
#'
#' `print()` shows a summarised output of the [efa_compare()] function: the mean
#' (with its range), median, and root mean squared distance (RMSE) of the
#' differences, the number of decimals to which all numbers agree, the minimum
#' number of decimals provided, and (for matrices) the number of differing
#' indicator-to-factor correspondences, followed (optionally) by the table of
#' elementwise differences. `format()` assembles the same report and returns it
#' as a character vector; `print()` is `cat(format(x), sep = "\n")`.
#' The lines follow the active console theme, so they are plain when colours are
#' disabled (for example when captured into a file or stripped with
#' [cli::ansi_strip()]).
#'
#' The line reporting the minimum number of decimals provided is shown only when it
#' carries information: two ordinary double matrices carry the full double precision,
#' for which the count is uninformative and the line is omitted.
#'
#' @param x An object of class `efa_compare` (output from [efa_compare()]).
#' @param digits,m_red,range_red,round_red,print_diff Display controls, documented
#'   in [efa_compare()]. Each defaults to `NULL`, meaning the value `efa_compare()`
#'   recorded in `x$settings` is used; supplying one overrides it for this call
#'   only, so the comparison need not be recomputed to change the printed report.
#' @param ... Passed from `print()` to `format()`; not otherwise used.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines (styled to the active console theme;
#'   plain when colours are disabled).
#'
#' @export
#' @method print efa_compare
#' @examples
#' # A type SPSS EFA to mimick the SPSS implementation
#' EFA_SPSS_5 <- efa_fit(IDS2_R, n_factors = 5,
#'                       estimate_control = estimate_control(type = "SPSS"),
#'                       rotate_control = rotate_control(type = "SPSS"))
#'
#' # A type psych EFA to mimick the psych::fa() implementation
#' EFA_psych_5 <- efa_fit(IDS2_R, n_factors = 5,
#'                        estimate_control = estimate_control(type = "psych"),
#'                        rotate_control = rotate_control(type = "psych"))
#'
#' # compare the two
#' comp <- efa_compare(EFA_SPSS_5$unrot_loadings, EFA_psych_5$unrot_loadings,
#'                     x_labels = c("SPSS", "psych"))
#' comp
#'
#' # format() returns the same lines as a character vector:
#' writeLines(format(comp))
#'
#' # the display settings can be changed without recomputing the comparison:
#' print(comp, digits = 2, print_diff = FALSE)
#'
print.efa_compare <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.efa_compare
#' @export
#' @method format efa_compare
format.efa_compare <- function(x, digits = NULL, m_red = NULL, range_red = NULL,
                               round_red = NULL, print_diff = NULL, ...) {

  # extract summary statistics
  diff <- x$diff
  mean_abs_diff <- x$mean_abs_diff
  median_abs_diff <- x$median_abs_diff
  min_abs_diff <- x$min_abs_diff
  max_abs_diff <- x$max_abs_diff
  max_dec <- x$max_dec
  are_equal <- x$are_equal
  g <- x$g
  diff_corres <- x$diff_corres
  diff_corres_cross <- x$diff_corres_cross

  # extract control settings
  corres <- x$settings$corres
  thresh <- x$settings$thresh

  # The display controls are recorded at construction but govern only how the report
  # reads, so an argument supplied here overrides the recorded value for this call and
  # anything left NULL falls back to it. Supplied values are validated exactly as
  # efa_compare() validates its own.
  if (is.null(digits)) digits <- x$settings$digits else checkmate::assert_count(digits)
  if (is.null(m_red)) m_red <- x$settings$m_red else checkmate::assert_number(m_red)
  if (is.null(range_red)) range_red <- x$settings$range_red else checkmate::assert_number(range_red)
  if (is.null(round_red)) round_red <- x$settings$round_red else checkmate::assert_number(round_red)
  if (is.null(print_diff)) print_diff <- x$settings$print_diff else checkmate::assert_flag(print_diff)

  # Style each statistic green when it clears its reduction threshold and red otherwise
  # (so a smaller difference reads as a closer match). The values keep their decimal
  # padding, so the lines are emitted verbatim below. They are differences between
  # loadings, i.e. bounded coefficients, so they follow the package's number convention
  # and drop the leading zero; the two decimal *counts* below are integers and do not.
  if (mean_abs_diff <= m_red) {
    mean_out <- .efa_style(.efa_num(mean_abs_diff, digits), c("green", "bold"))
  } else {
    mean_out <- .efa_style(.efa_num(mean_abs_diff, digits), c("red", "bold"))
  }

  if (median_abs_diff <= m_red) {
    median_out <- .efa_style(.efa_num(median_abs_diff, digits), c("green", "bold"))
  } else {
    median_out <- .efa_style(.efa_num(median_abs_diff, digits), c("red", "bold"))
  }

  if (max_abs_diff <= range_red) {
    max_out <- .efa_style(.efa_num(max_abs_diff, digits), c("green", "bold"))
    min_out <- .efa_style(.efa_num(min_abs_diff, digits), c("green", "bold"))
  } else {
    max_out <- .efa_style(.efa_num(max_abs_diff, digits), c("red", "bold"))
    min_out <- .efa_style(.efa_num(min_abs_diff, digits), c("red", "bold"))
  }

  if (is.na(are_equal)) {
    equal_out <- .efa_style("none", c("red", "bold"))
  } else if (are_equal < round_red) {
    equal_out <- .efa_style(are_equal, c("red", "bold"))
  } else {
    equal_out <- .efa_style(are_equal, c("green", "bold"))
  }

  cli::cli_format_method({
    .print_efa_rule("Summary statistics")

    # Fixed-format statistic lines: the values carry decimal padding and conditional
    # green/red styling, so emit them verbatim (cli never reflows them mid-token).
    cli::cli_verbatim(paste0("Mean [min, max] absolute difference: ",
                             mean_out, " [", min_out, ", ", max_out, "]"))
    cli::cli_verbatim(paste0("Median absolute difference: ", median_out))
    cli::cli_verbatim(paste0("Root mean squared distance (RMSE): ",
                             .efa_style(.efa_num(g, digits), "bold")))
    cli::cli_verbatim(paste0("Max decimals where all numbers agree in absolute value: ",
                             equal_out))

    # The decimals actually carried by the inputs only bound the comparison when one of
    # them was rounded. `.decimals()` renders with `digits = 15`, so an unrounded double
    # always reaches that many decimal places (more for magnitudes below 1), and the count
    # then says nothing about the two solutions; only a smaller count is informative. The
    # threshold therefore tracks the `digits =` argument of `.decimals()` in
    # R/format-helpers.R and must be revisited with it.
    if (max_dec < 15) {
      cli::cli_verbatim(paste0("Minimum number of decimals provided: ",
                               .efa_style(max_dec, "bold")))
    }

    # Differing indicator-to-factor correspondences. Only reported when they were
    # actually compared: they are undefined for vector input (NA), and both
    # `corres = FALSE` and a single-factor matrix skip the comparison and leave a
    # placeholder 0, which must not be printed as if the two solutions had been found to
    # agree. The single-factor case is recognised the same way `.compare_loadings()`
    # decides to skip it, on the number of columns compared. Green when they do agree
    # (0 differing), red otherwise -- mirroring the "smaller difference reads as a closer
    # match" colouring of the statistics above.
    if (isTRUE(corres) && !is.na(diff_corres) && is.matrix(diff) && ncol(diff) > 1L) {
      corres_style <- function(v) {
        .efa_style(v, if (v == 0) c("green", "bold") else c("red", "bold"))
      }
      cli::cli_verbatim(paste0(
        "Differing indicator-to-factor correspondences: ",
        corres_style(diff_corres), " (highest loading), ",
        corres_style(diff_corres_cross), " (all |loadings| >= ", format(thresh), ")"))
    }

    if (isTRUE(print_diff)) {
      .print_efa_rule("Elementwise differences")
      .efa_emit_lines(.compare_diff_lines(diff, digits = digits, r_red = range_red))
    }
  })
}

# Render an efa_compare difference object through the shared matrix renderer into verbatim
# lines, colouring cells whose absolute difference exceeds `r_red` (the efa_compare
# `range_red` threshold) red. A matrix is shown with its variable rows and factor columns; a vector
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
