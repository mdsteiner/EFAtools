#' Print a loading matrix
#'
#' @details
#' The method prints a loading matrix in a compact, console-oriented table.
#' Loadings with absolute value greater than or equal to `cutoff` are emphasized,
#' smaller loadings are de-emphasized, and Heywood-relevant communality/
#' uniqueness values are marked when `h2` is supplied. Long variable names can
#' be truncated, abbreviated, or printed in full. If the matrix has many factor
#' columns, the table is split into column blocks so that the output remains
#' readable in narrower consoles.
#'
#' If `h2` is named and `x` has row names, `h2` is matched to the row names of
#' `x` before any optional row sorting is applied. If `x` has no row names, a
#' named `h2` vector is used in the supplied order.
#'
#' @param x a loading matrix of class `efa_loadings` (or the legacy `LOADINGS`).
#' @param cutoff numeric. The value at or above which loadings are emphasized;
#'  default is .3.
#' @param digits numeric. Passed to \code{\link[base:Round]{round}}. Number of digits
#'  to round the loadings to (default is 3).
#' @param max_name_length numeric. The maximum length of the variable names to
#'  display. Everything beyond this will be cut from the right unless
#'  `name_style = "abbreviate"` or `name_style = "full"` is used. Cutting never
#'  leaves two variables sharing a row label: if it would (as for items with a
#'  long common prefix), the names are abbreviated instead, and numbered if
#'  needed.
#' @param h2 numeric. Vector of communalities to print. If named and `x`
#'  has row names, names are used to align communalities to rows.
#' @param color logical. Whether to apply console styling using \pkg{cli}.
#'  Default is `TRUE`.
#' @param name_style character. How to shorten variable names longer than
#'  `max_name_length`. `"truncate"` cuts names from the right,
#'  `"abbreviate"` uses [base::abbreviate()], and
#'  `"full"` prints full names.
#' @param max_factor_name_length numeric or `NULL`. Optional maximum length
#'  of factor names. If `NULL`, factor names are not shortened.
#' @param max_factors_per_block numeric or `NULL`. Maximum number of factor
#'  columns to print per block. If `NULL`, the number is chosen from the
#'  console width.
#' @param sort_loadings character. Optional row sorting. `"none"` preserves
#'  the input order, `"primary"` groups rows by the factor with the largest
#'  absolute loading, and `"clustered"` additionally sorts within each
#'  factor by the size of the primary loading.
#' @param legend logical. Whether to append a short explanation of the styling.
#'  Default is `FALSE` for standalone loading matrices. The legend is only printed
#'  when the styling it describes is actually rendered (`color = TRUE` and a
#'  colour-capable console); in plain output it is omitted and this argument has
#'  no effect.
#' @param ... additional arguments passed to print or format
#'
#' @returns `print()` returns its argument `x` invisibly; it is
#'   `cat(format(x, ...), sep = "\n")` followed by a blank line for console
#'   spacing. `format()` returns a character vector with the table lines (styled
#'   to the active console theme; plain when colours are disabled).
#'
#' @method print efa_loadings
#' @export
#'
#' @examples
#' EFAtools_PAF <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                         estimator = "PAF", rotation = "promax")
#' EFAtools_PAF
#'
#' # format() returns the same lines as a character vector:
#' writeLines(format(EFAtools_PAF$rot_loadings))
#'
print.efa_loadings <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  # a trailing blank line, so a loading table printed at the console is set off from
  # whatever follows it
  cat("\n")

  invisible(x)
}

#' @rdname print.efa_loadings
#' @method format efa_loadings
#' @export
format.efa_loadings <- function(x, cutoff = .3, digits = 3, max_name_length = 10,
                           h2 = NULL, color = TRUE,
                           name_style = c("truncate", "abbreviate", "full"),
                           max_factor_name_length = NULL,
                           max_factors_per_block = NULL,
                           sort_loadings = c("none", "primary", "clustered"),
                           legend = FALSE, ...) {

  name_style <- .match_arg_ci(name_style)
  sort_loadings <- .match_arg_ci(sort_loadings)

  .validate_loadings_print_args(
    x = x,
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    h2 = h2,
    color = color,
    max_factor_name_length = max_factor_name_length,
    max_factors_per_block = max_factors_per_block,
    legend = legend
  )

  spec <- .loadings_print_spec(
    x = x,
    h2 = h2,
    max_name_length = max_name_length,
    name_style = name_style,
    max_factor_name_length = max_factor_name_length,
    sort_loadings = sort_loadings
  )

  cli::cli_format_method({
    # Through .efa_emit_lines(), not cli_verbatim() directly: a table too wide for the
    # console is rendered as stacked column blocks separated by a blank line, and
    # cli_verbatim() drops interior empty strings, which would run the blocks together.
    .efa_emit_lines(.efa_format_matrix(
      values = spec$values,
      row_labels = spec$var_names,
      col_labels = spec$factor_names,
      col_roles = spec$col_type,
      cutoff = cutoff,
      digits = digits,
      color = color,
      max_factors_per_block = max_factors_per_block
    ))

    # The legend only describes console styling, so it is emitted only when that styling
    # is actually rendered. In plain output the loadings themselves carry the information,
    # and a legend naming marks the reader cannot see would be misleading.
    if (isTRUE(legend) && .efa_styling_visible(color)) {
      cli::cli_verbatim("")
      cli::cli_verbatim(.efa_loadings_legend(cutoff, spec$has_h2, digits, color))
    }
  })
}

.validate_loadings_print_args <- function(x, cutoff, digits, max_name_length,
                                          h2, color, max_factor_name_length,
                                          max_factors_per_block, legend) {
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort("{.arg x} must be a numeric matrix.",
                   class = "efa_print_invalid_x")
  }

  if (nrow(x) < 1L || ncol(x) < 1L) {
    cli::cli_abort("{.arg x} must have at least one row and one column.",
                   class = "efa_print_invalid_x")
  }

  .efa_check_nonneg_number(cutoff, "cutoff", "efa_print_invalid_cutoff")
  .efa_check_nonneg_int(digits, "digits", "efa_print_invalid_digits")
  .efa_check_pos_int(max_name_length, "max_name_length",
                     "efa_print_invalid_max_name_length")
  .efa_check_opt_pos_int(max_factor_name_length, "max_factor_name_length",
                         "efa_print_invalid_max_factor_name_length")
  .efa_check_opt_pos_int(max_factors_per_block, "max_factors_per_block",
                         "efa_print_invalid_max_factors_per_block")
  .efa_check_flag(color, "color", "efa_print_invalid_color")
  .efa_check_flag(legend, "legend", "efa_print_invalid_legend")

  if (!is.null(h2)) {
    if (!is.numeric(h2) || length(h2) != nrow(x)) {
      cli::cli_abort("{.arg h2} must be a numeric vector with one value per row of {.arg x}.",
                     class = "efa_print_invalid_h2")
    }
  }

  invisible(TRUE)
}

.loadings_print_spec <- function(x, h2, max_name_length,
                                 name_style, max_factor_name_length,
                                 sort_loadings) {
  x <- as.matrix(x)
  n_factors <- ncol(x)

  factor_names <- colnames(x)
  if (is.null(factor_names)) {
    factor_names <- paste0("F", seq_len(n_factors))
  }

  original_var_names <- rownames(x)
  has_row_names <- !is.null(original_var_names)
  var_names <- original_var_names
  if (!has_row_names) {
    var_names <- paste0("V", seq_len(nrow(x)))
  }

  has_h2 <- !is.null(h2)
  if (has_h2) {
    h2 <- .align_loadings_h2(
      h2 = h2,
      var_names = var_names,
      has_row_names = has_row_names
    )
  }

  # Shorten before reordering, so a row label is a property of the variable rather than of
  # where it happens to sit: the numbering that resolves a collision is assigned in vector
  # order, and deriving it from the sorted order would label the same item differently under
  # different `sort_loadings`.
  display_var_names <- .shorten_loadings_names(
    var_names,
    max_length = max_name_length,
    name_style = name_style
  )

  row_order <- .efa_loading_row_order(x, sort_loadings)
  x <- x[row_order, , drop = FALSE]
  display_var_names <- display_var_names[row_order]
  if (has_h2) {
    h2 <- h2[row_order]
  }

  col_type <- rep("loading", n_factors)
  values <- x

  if (has_h2) {
    values <- cbind(values, "h2" = h2, "u2" = 1 - h2)
    factor_names <- c(factor_names, "h2", "u2")
    col_type <- c(col_type, "h2", "u2")
  }

  display_factor_names <- .shorten_loadings_factor_names(
    factor_names,
    max_length = max_factor_name_length
  )

  list(
    values = values,
    var_names = display_var_names,
    factor_names = display_factor_names,
    col_type = col_type,
    has_h2 = has_h2
  )
}

.align_loadings_h2 <- function(h2, var_names, has_row_names = TRUE) {
  # Only align by name when the loading matrix itself has row names. Generated
  # fallback names (V1, V2, ...) should not be used to reorder a named h2 vector.
  h2_names <- names(h2)
  has_complete_h2_names <- !is.null(h2_names) && all(nzchar(h2_names))

  if (isTRUE(has_row_names) && has_complete_h2_names) {
    if (!all(var_names %in% h2_names)) {
      cli::cli_abort(
        "If {.arg h2} is named and {.arg x} has row names, {.code names(h2)} must include all row names of {.arg x}.",
        class = "efa_print_invalid_h2"
      )
    }
    h2 <- h2[var_names]
  }

  as.numeric(h2)
}

# Shorten variable names to `max_length` for a loading table. Cutting from the right can map
# distinct items onto one row label -- "premeditation_1" through "premeditation_11" all become
# "premeditat" -- leaving rows the reader cannot tell apart, with nothing on screen to signal
# it. So a truncation that collides falls back to abbreviation and, if that still collides, to
# numbered labels: a row label is never shared by two different variables. Names the caller
# supplied as duplicates are left as they are; those collisions are not the shortening's doing.
.shorten_loadings_names <- function(x, max_length, name_style) {
  x <- as.character(x)

  if (identical(name_style, "full")) {
    return(x)
  }

  if (all(nchar(x) <= max_length)) {
    return(x)
  }

  if (identical(name_style, "abbreviate")) {
    return(.abbreviate_loadings_names(x, max_length))
  }

  out <- substr(x, 1L, max_length)
  if (anyDuplicated(out) < 1L || anyDuplicated(x) > 0L) {
    return(out)
  }

  # `abbreviate()` warns on non-ASCII input, and a table repairing its own row labels must
  # not raise a warning for doing so; such names go straight to the numbering below. The
  # test is `abbreviate()`'s own, so the two agree on what counts as non-ASCII.
  if (all(nchar(x, type = "bytes") == nchar(x, type = "chars"))) {
    abbreviated <- .abbreviate_loadings_names(x, max_length)
    if (anyDuplicated(abbreviated) < 1L) {
      return(abbreviated)
    }
    out <- abbreviated
  }

  .disambiguate_loadings_names(out, max_length)
}

# `abbreviate()` with `strict = TRUE` never exceeds `minlength`, so the result still respects
# the requested width; `method = "both.sides"` also drops characters from a shared prefix,
# which is what separates names that differ only near their end.
.abbreviate_loadings_names <- function(x, max_length) {
  as.character(abbreviate(
    x,
    minlength = max_length,
    strict = TRUE,
    method = "both.sides"
  ))
}

# Last-resort disambiguation of labels that are still duplicated: the numeric suffix
# `make.unique()` appends is budgeted out of `max_length` by trimming the label it is appended
# to, so a disambiguated label keeps the requested width instead of growing past it. Trimming
# can in principle re-collide, and a width too narrow to hold a suffix cannot be budgeted at
# all; both fall back to the untrimmed form, which `make.unique()` guarantees is unique. An
# over-wide label is still readable, an ambiguous one is not.
.disambiguate_loadings_names <- function(x, max_length) {
  suffixed <- make.unique(x, sep = "~")
  suffix <- substr(suffixed, nchar(x) + 1L, nchar(suffixed))

  keep <- max_length - nchar(suffix)
  if (any(keep < 1L)) {
    return(suffixed)
  }

  out <- paste0(substr(x, 1L, keep), suffix)
  if (anyDuplicated(out) > 0L) suffixed else out
}

.shorten_loadings_factor_names <- function(x, max_length = NULL) {
  x <- as.character(x)

  if (is.null(max_length) || all(nchar(x) <= max_length)) {
    return(x)
  }

  substr(x, 1L, max_length)
}

# Build the cli legend lines describing the loading-matrix styling. Styling is dropped when
# `color = FALSE` or colours are off, so the lines embed cleanly in plain output.
.efa_loadings_legend <- function(cutoff, has_h2, digits, color = TRUE) {
  cutoff_str <- .efa_num(cutoff, digits = digits, pad = FALSE)

  bold_word <- if (isTRUE(color)) cli::style_bold("bold") else "bold"
  grey_word <- if (isTRUE(color)) cli::col_grey("grey") else "grey"

  lines <- c(
    "Legend:",
    paste0("  ", bold_word, " = |loading| >= ", cutoff_str),
    paste0("  ", grey_word, " = below cutoff")
  )

  if (isTRUE(has_h2)) {
    red_word <- if (isTRUE(color)) {
      cli::style_bold(cli::col_red("red h2/u2"))
    } else {
      "red h2/u2"
    }
    lines <- c(lines, paste0("  ", red_word, " = Heywood-relevant value"))
  }

  lines
}
