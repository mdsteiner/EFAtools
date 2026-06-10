#' Print LOADINGS object
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
#' @param x class LOADINGS matrix.
#' @param cutoff numeric. The number above which to print loadings in bold
#'  default is .3.
#' @param digits numeric. Passed to \code{\link[base:Round]{round}}. Number of digits
#'  to round the loadings to (default is 3).
#' @param max_name_length numeric. The maximum length of the variable names to
#'  display. Everything beyond this will be cut from the right unless
#'  `name_style = "abbreviate"` or `name_style = "full"` is used.
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
#'  Default is `FALSE` for standalone loading matrices.
#' @param ... additional arguments passed to print or format
#'
#' @method print LOADINGS
#' @export
#'
#' @examples
#' EFAtools_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                     type = "EFAtools", method = "PAF", rotation = "promax")
#' EFAtools_PAF
#'
print.LOADINGS <- function(x, cutoff = .3, digits = 3, max_name_length = 10,
                           h2 = NULL, color = TRUE,
                           name_style = c("truncate", "abbreviate", "full"),
                           max_factor_name_length = NULL,
                           max_factors_per_block = NULL,
                           sort_loadings = c("none", "primary", "clustered"),
                           legend = FALSE, ...) {

  out <- format.LOADINGS(x,
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    h2 = h2,
    color = color,
    name_style = name_style,
    max_factor_name_length = max_factor_name_length,
    max_factors_per_block = max_factors_per_block,
    sort_loadings = sort_loadings,
    legend = legend,
    ...
  )

  cat(out, sep = "\n")
  cat("\n")

  invisible(x)
}

#' @rdname print.LOADINGS
#' @method format LOADINGS
#' @export
format.LOADINGS <- function(x, cutoff = .3, digits = 3, max_name_length = 10,
                            h2 = NULL, color = TRUE,
                            name_style = c("truncate", "abbreviate", "full"),
                            max_factor_name_length = NULL,
                            max_factors_per_block = NULL,
                            sort_loadings = c("none", "primary", "clustered"),
                            legend = FALSE, ...) {

  name_style <- match.arg(name_style)
  sort_loadings <- match.arg(sort_loadings)

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
    cli::cli_verbatim(.efa_format_matrix(
      values = spec$values,
      row_labels = spec$var_names,
      col_labels = spec$factor_names,
      col_roles = spec$col_type,
      cutoff = cutoff,
      digits = digits,
      color = color,
      max_factors_per_block = max_factors_per_block
    ))

    if (isTRUE(legend)) {
      cli::cli_verbatim("")
      cli::cli_verbatim(.efa_loadings_legend(cutoff, spec$has_h2, digits, color))
    }
  })
}

.validate_loadings_print_args <- function(x, cutoff, digits, max_name_length,
                                          h2, color, max_factor_name_length,
                                          max_factors_per_block, legend) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix.", call. = FALSE)
  }

  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop("`x` must have at least one row and one column.", call. = FALSE)
  }

  if (!is.numeric(cutoff) || length(cutoff) != 1L ||
      !is.finite(cutoff) || cutoff < 0) {
    stop("`cutoff` must be a single finite non-negative number.", call. = FALSE)
  }

  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 0 || digits != as.integer(digits)) {
    stop("`digits` must be a single finite non-negative integer.", call. = FALSE)
  }

  if (!is.numeric(max_name_length) || length(max_name_length) != 1L ||
      !is.finite(max_name_length) || max_name_length < 1 ||
      max_name_length != as.integer(max_name_length)) {
    stop("`max_name_length` must be a single finite positive integer.", call. = FALSE)
  }

  if (!is.null(max_factor_name_length) &&
      (!is.numeric(max_factor_name_length) || length(max_factor_name_length) != 1L ||
       !is.finite(max_factor_name_length) || max_factor_name_length < 1 ||
       max_factor_name_length != as.integer(max_factor_name_length))) {
    stop("`max_factor_name_length` must be NULL or a single finite positive integer.", call. = FALSE)
  }

  if (!is.null(max_factors_per_block) &&
      (!is.numeric(max_factors_per_block) || length(max_factors_per_block) != 1L ||
       !is.finite(max_factors_per_block) || max_factors_per_block < 1 ||
       max_factors_per_block != as.integer(max_factors_per_block))) {
    stop("`max_factors_per_block` must be NULL or a single finite positive integer.", call. = FALSE)
  }

  if (!is.logical(color) || length(color) != 1L || is.na(color)) {
    stop("`color` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(legend) || length(legend) != 1L || is.na(legend)) {
    stop("`legend` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(h2)) {
    if (!is.numeric(h2) || length(h2) != nrow(x)) {
      stop("`h2` must be a numeric vector with one value per row of `x`.", call. = FALSE)
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

  row_order <- .loadings_row_order(x, sort_loadings)
  x <- x[row_order, , drop = FALSE]
  var_names <- var_names[row_order]
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

  display_var_names <- .shorten_loadings_names(
    var_names,
    max_length = max_name_length,
    name_style = name_style
  )

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

.loadings_row_order <- function(x, sort_loadings = c("none", "primary", "clustered")) {
  sort_loadings <- match.arg(sort_loadings)

  if (identical(sort_loadings, "none") || nrow(x) < 2L || ncol(x) < 1L) {
    return(seq_len(nrow(x)))
  }

  abs_x <- abs(x)
  abs_x[is.na(abs_x)] <- -Inf
  primary_factor <- max.col(abs_x, ties.method = "first")
  primary_loading <- apply(abs_x, 1L, max, na.rm = TRUE)
  primary_loading[!is.finite(primary_loading)] <- -Inf

  if (identical(sort_loadings, "primary")) {
    return(order(primary_factor, seq_along(primary_factor)))
  }

  order(primary_factor, -primary_loading, seq_along(primary_factor))
}

.align_loadings_h2 <- function(h2, var_names, has_row_names = TRUE) {
  # Only align by name when the loading matrix itself has row names. Generated
  # fallback names (V1, V2, ...) should not be used to reorder a named h2 vector.
  h2_names <- names(h2)
  has_complete_h2_names <- !is.null(h2_names) && all(nzchar(h2_names))

  if (isTRUE(has_row_names) && has_complete_h2_names) {
    if (!all(var_names %in% h2_names)) {
      stop(
        "If `h2` is named and `x` has row names, names(h2) must include all row names of `x`.",
        call. = FALSE
      )
    }
    h2 <- h2[var_names]
  }

  as.numeric(h2)
}

.shorten_loadings_names <- function(x, max_length, name_style) {
  x <- as.character(x)

  if (identical(name_style, "full")) {
    return(x)
  }

  if (all(nchar(x) <= max_length)) {
    return(x)
  }

  if (identical(name_style, "abbreviate")) {
    return(as.character(abbreviate(
      x,
      minlength = max_length,
      strict = TRUE,
      method = "both.sides"
    )))
  }

  substr(x, 1L, max_length)
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
