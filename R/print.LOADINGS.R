#' Print LOADINGS object
#'
#' @param x class LOADINGS matrix.
#' @param cutoff numeric. The number above which to print loadings in bold
#'  default is .3.
#' @param digits numeric. Passed to \code{\link[base:Round]{round}}. Number of digits
#'  to round the loadings to (default is 3).
#' @param max_name_length numeric. The maximum length of the variable names to
#'  display. Everything beyond this will be cut from the right unless
#'  \code{name_style = "abbreviate"} or \code{name_style = "full"} is used.
#' @param h2 numeric. Vector of communalities to print. If named and \code{x}
#'  has row names, names are used to align communalities to rows.
#' @param color logical. Whether to apply console styling using \pkg{crayon}.
#'  Default is \code{TRUE}.
#' @param name_style character. How to shorten variable names longer than
#'  \code{max_name_length}. \code{"truncate"} cuts names from the right,
#'  \code{"abbreviate"} uses \code{\link[base:abbreviate]{abbreviate}}, and
#'  \code{"full"} prints full names.
#' @param max_factor_name_length numeric or \code{NULL}. Optional maximum length
#'  of factor names. If \code{NULL}, factor names are not shortened.
#' @param max_factors_per_block numeric or \code{NULL}. Maximum number of factor
#'  columns to print per block. If \code{NULL}, the number is chosen from the
#'  console width.
#' @param sort_loadings character. Optional row sorting. \code{"none"} preserves
#'  the input order, \code{"primary"} groups rows by the factor with the largest
#'  absolute loading, and \code{"clustered"} additionally sorts within each
#'  factor by the size of the primary loading.
#' @param legend logical. Whether to append a short explanation of the styling.
#'  Default is \code{FALSE} for standalone loading matrices.
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
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    name_style = name_style,
    max_factor_name_length = max_factor_name_length,
    sort_loadings = sort_loadings
  )

  .format_loadings_table(
    spec,
    color = color,
    max_factors_per_block = max_factors_per_block,
    legend = legend
  )
}

.validate_loadings_print_args <- function(x, cutoff, digits, max_name_length,
                                          h2, color, max_factor_name_length,
                                          max_factors_per_block, legend) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix.", call. = FALSE)
  }

  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff) || cutoff < 0) {
    stop("`cutoff` must be a single non-negative number.", call. = FALSE)
  }

  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      digits < 0 || digits != as.integer(digits)) {
    stop("`digits` must be a single non-negative integer.", call. = FALSE)
  }

  if (!is.numeric(max_name_length) || length(max_name_length) != 1L ||
      is.na(max_name_length) || max_name_length < 1 ||
      max_name_length != as.integer(max_name_length)) {
    stop("`max_name_length` must be a single positive integer.", call. = FALSE)
  }

  if (!is.null(max_factor_name_length) &&
      (!is.numeric(max_factor_name_length) || length(max_factor_name_length) != 1L ||
       is.na(max_factor_name_length) || max_factor_name_length < 1 ||
       max_factor_name_length != as.integer(max_factor_name_length))) {
    stop("`max_factor_name_length` must be NULL or a single positive integer.", call. = FALSE)
  }

  if (!is.null(max_factors_per_block) &&
      (!is.numeric(max_factors_per_block) || length(max_factors_per_block) != 1L ||
       is.na(max_factors_per_block) || max_factors_per_block < 1 ||
       max_factors_per_block != as.integer(max_factors_per_block))) {
    stop("`max_factors_per_block` must be NULL or a single positive integer.", call. = FALSE)
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

.loadings_print_spec <- function(x, h2, cutoff, digits, max_name_length,
                                 name_style, max_factor_name_length,
                                 sort_loadings) {
  x <- as.matrix(x)
  n_factors <- ncol(x)

  factor_names <- colnames(x)
  if (is.null(factor_names)) {
    factor_names <- paste0("F", seq_len(n_factors))
  }

  var_names <- rownames(x)
  if (is.null(var_names)) {
    var_names <- paste0("V", seq_len(nrow(x)))
  }

  has_h2 <- !is.null(h2)
  if (has_h2) {
    h2 <- .align_loadings_h2(h2, var_names)
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

  value_strings <- .format_loadings_values(values, digits = digits)

  list(
    values = values,
    value_strings = value_strings,
    var_names = display_var_names,
    factor_names = display_factor_names,
    col_type = col_type,
    cutoff = cutoff,
    digits = digits,
    n_factors = n_factors,
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

.align_loadings_h2 <- function(h2, var_names) {
  if (!is.null(names(h2)) && !is.null(var_names) && all(nzchar(names(h2)))) {
    if (!all(var_names %in% names(h2))) {
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
  if (identical(name_style, "full")) {
    return(x)
  }

  if (all(nchar(x) <= max_length)) {
    return(x)
  }

  if (identical(name_style, "abbreviate")) {
    return(as.character(utils::abbreviate(
      x,
      minlength = max_length,
      strict = TRUE,
      method = "both.sides"
    )))
  }

  substr(x, 1L, max_length)
}

.shorten_loadings_factor_names <- function(x, max_length = NULL) {
  if (is.null(max_length) || all(nchar(x) <= max_length)) {
    return(x)
  }

  substr(x, 1L, max_length)
}

.format_loadings_values <- function(x, digits) {
  out <- matrix("", nrow = nrow(x), ncol = ncol(x))

  for (jj in seq_len(ncol(x))) {
    for (ii in seq_len(nrow(x))) {
      out[ii, jj] <- .format_loadings_number(x[ii, jj], digits = digits)
    }
  }

  out
}

.format_loadings_number <- function(x, digits) {
  if (is.na(x)) {
    return("NA")
  }

  .numformat(round(x, digits = digits), digits = digits)
}

.format_loadings_table <- function(spec, color = TRUE,
                                   max_factors_per_block = NULL,
                                   legend = FALSE) {
  row_width <- max(nchar(spec$var_names))
  col_widths <- pmax(
    nchar(spec$factor_names),
    apply(spec$value_strings, 2L, function(z) max(nchar(z)))
  )

  blocks <- .loadings_column_blocks(
    spec = spec,
    row_width = row_width,
    col_widths = col_widths,
    max_factors_per_block = max_factors_per_block
  )

  out <- character(0)
  multi_block <- length(blocks) > 1L

  for (bb in seq_along(blocks)) {
    cols <- blocks[[bb]]

    if (multi_block) {
      if (bb > 1L) {
        out <- c(out, "")
      }
      out <- c(out, .loadings_block_label(spec, cols, bb, length(blocks), color))
    }

    out <- c(out, .format_loadings_block(
      spec = spec,
      cols = cols,
      row_width = row_width,
      col_widths = col_widths,
      color = color
    ))
  }

  if (isTRUE(legend)) {
    out <- c(out, "", .loadings_legend(spec, color = color))
  }

  out
}

.loadings_column_blocks <- function(spec, row_width, col_widths,
                                    max_factors_per_block = NULL) {
  factor_cols <- seq_len(spec$n_factors)
  aux_cols <- if (isTRUE(spec$has_h2)) {
    (spec$n_factors + 1L):length(spec$col_type)
  } else {
    integer(0)
  }

  if (length(factor_cols) <= 1L) {
    return(list(c(factor_cols, aux_cols)))
  }

  if (!is.null(max_factors_per_block)) {
    split_points <- ceiling(seq_along(factor_cols) / max_factors_per_block)
    return(lapply(split(factor_cols, split_points), function(cols) c(cols, aux_cols)))
  }

  console_width <- .loadings_console_width()
  full_width <- .loadings_table_width(c(factor_cols, aux_cols), row_width, col_widths)
  if (full_width <= console_width) {
    return(list(c(factor_cols, aux_cols)))
  }

  blocks <- list()
  current <- integer(0)

  for (col in factor_cols) {
    candidate <- c(current, col, aux_cols)
    if (length(current) > 0L &&
        .loadings_table_width(candidate, row_width, col_widths) > console_width) {
      blocks[[length(blocks) + 1L]] <- c(current, aux_cols)
      current <- col
    } else {
      current <- c(current, col)
    }
  }

  if (length(current) > 0L) {
    blocks[[length(blocks) + 1L]] <- c(current, aux_cols)
  }

  blocks
}

.loadings_console_width <- function() {
  width <- tryCatch(cli::console_width(), error = function(e) NA_integer_)
  if (!is.finite(width) || width < 40L) {
    return(80L)
  }
  as.integer(width)
}

.loadings_table_width <- function(cols, row_width, col_widths) {
  # Two spaces between the row-label column and every value column.
  row_width + 2L + sum(col_widths[cols]) + 2L * (length(cols) - 1L)
}

.loadings_block_label <- function(spec, cols, block, n_blocks, color) {
  factor_cols <- cols[cols <= spec$n_factors]
  first <- spec$factor_names[min(factor_cols)]
  last <- spec$factor_names[max(factor_cols)]

  label <- if (identical(first, last)) {
    paste0("Factors ", first, " (block ", block, "/", n_blocks, ")")
  } else {
    paste0("Factors ", first, "-", last, " (block ", block, "/", n_blocks, ")")
  }

  .loadings_style(label, style = "label", color = color)
}

.format_loadings_block <- function(spec, cols, row_width, col_widths, color = TRUE) {
  header <- .loadings_style(
    stringr::str_c(
      stringr::str_pad("", row_width),
      "  ",
      stringr::str_c(
        stringr::str_pad(spec$factor_names[cols], col_widths[cols], side = "both"),
        collapse = "  "
      )
    ),
    style = "label",
    color = color
  )

  rows <- vapply(seq_len(nrow(spec$values)), function(ii) {
    cells <- vapply(cols, function(jj) {
      cell <- stringr::str_pad(spec$value_strings[ii, jj], col_widths[jj], side = "left")
      .style_loadings_cell(
        cell = cell,
        value = spec$values[ii, jj],
        col_type = spec$col_type[jj],
        cutoff = spec$cutoff,
        color = color
      )
    }, character(1L))

    stringr::str_c(
      .loadings_style(
        stringr::str_pad(spec$var_names[ii], row_width, side = "right"),
        style = "label",
        color = color
      ),
      "  ",
      stringr::str_c(cells, collapse = "  ")
    )
  }, character(1L))

  c(header, rows)
}

.style_loadings_cell <- function(cell, value, col_type, cutoff, color) {
  if (identical(col_type, "loading")) {
    if (is.na(value)) {
      return(cell)
    }

    if (abs(value) < cutoff) {
      return(.loadings_style(cell, style = "weak", color = color))
    }

    return(.loadings_style(cell, style = "salient", color = color))
  }

  if (identical(col_type, "h2") && is.finite(value) && value > 1) {
    return(.loadings_style(cell, style = "heywood", color = color))
  }

  if (identical(col_type, "u2") && is.finite(value) && value < 0) {
    return(.loadings_style(cell, style = "heywood", color = color))
  }

  cell
}

.loadings_legend <- function(spec, color = TRUE) {
  cutoff <- .format_loadings_number(spec$cutoff, digits = spec$digits)

  if (isTRUE(color)) {
    text <- paste0(
      crayon::italic("Legend: \n"),
      crayon::bold("  bold"), " = |loading| >= ", cutoff,
      "\n", crayon::silver("  grey"), " = below cutoff"
    )
    if (isTRUE(spec$has_h2)) {
      text <- paste0(text, "\n", crayon::red$bold("  red h2/u2"), " = Heywood-relevant value")
    }
    return(text)
  }

  text <- paste0(crayon::italic("Legend: "), "salient loadings use cutoff |loading| >= ", cutoff)
  if (isTRUE(spec$has_h2)) {
    text <- paste0(text, "; h2 > 1 or u2 < 0 indicates a Heywood-relevant value")
  }
  text
}

.loadings_style <- function(x, style = c("label", "weak", "salient", "heywood"),
                            color = TRUE) {
  style <- match.arg(style)

  if (!isTRUE(color)) {
    return(x)
  }

  switch(style,
    label = crayon::blue(x),
    weak = crayon::silver(x),
    salient = crayon::bold(x),
    heywood = crayon::red$bold(x)
  )
}
