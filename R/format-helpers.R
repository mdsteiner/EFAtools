# Number formatting and styled-table builders shared by the print methods: padded decimal
# formatting, the LOADINGS/COMPARE cell renderers, and small settings/CI string helpers.

# Apply named logical styles to a pre-assembled string, mapping each to its cli primitive
# (in order). A no-op when `style` is empty or when colours are off (the cli primitives are
# themselves no-ops then). One home for all print colouring. Blue is
# intentionally absent: rules and labels use cli's default colour.
.efa_style <- function(x, style = NULL) {
  styles <- list(red = cli::col_red, green = cli::col_green, yellow = cli::col_yellow,
                 bold = cli::style_bold, italic = cli::style_italic)
  for (s in style) {
    x <- styles[[s]](x)
  }
  x
}

#' Format numbers for print method
#'
#' Helper function used in the print method for class LOADINGS and SLLOADINGS.
#' Strips the 0 in front of the decimal point of a number if number < 1, only
#' keeps the first `digits` number of digits, and adds an empty space in
#' front of the number if the number is positive. This way all returned strings
#' (except for those > 1, which are exceptions in LOADINGS) have the same number
#' of characters.
#'
#' @param x numeric. Number to be formatted.
#' @param digits numeric. Number of digits after the comma to keep.
#' @param print_zero logical. Whether, if a number is between \[-1, 1\], the
#'  zero should be omitted or printed (default is FALSE, i.e. omit zeros).
#' @param pad logical. Whether, if a number starts with a 0 and the 0 is not printed
#'  a white-space should be added.
#'
#' @return A formated number
.numformat <- function(x, digits = 2, print_zero = FALSE,
                       pad = TRUE) {

  if (isFALSE(print_zero)) {

    ncode <- paste0("%.", digits, "f")
    x <- sub("^(-?)0.", "\\1.", sprintf(ncode, x))
    if (isTRUE(pad)) {
      x <- stringr::str_pad(x, digits + 2, "left")
    }


  } else {

    ncode <- paste0("%.", digits, "f")
    x <- sprintf(ncode, x)

    if (isTRUE(pad)) {
      x <- stringr::str_pad(x, digits + 3, "left")
    }

  }
  x

}

# Canonical number formatter for printed output: rounds to `digits`, drops the leading 0 of
# values in (-1, 1) unless `print_zero`, optionally left-pads to a common width, and renders
# non-finite values as "NA". Accepts a scalar, vector, or matrix and returns the same shape
# (dimnames preserved). Produces the same strings as the older `.numformat`-style helpers.
.efa_num <- function(x, digits = 3, print_zero = FALSE, pad = TRUE) {

  dims <- dim(x)
  dn <- dimnames(x)
  v <- as.numeric(x)
  finite <- is.finite(v)

  out <- sprintf(paste0("%.", digits, "f"), v)
  if (isFALSE(print_zero)) {
    out <- sub("^(-?)0.", "\\1.", out)
  }
  out[!finite] <- "NA"

  if (isTRUE(pad)) {
    width <- if (isFALSE(print_zero)) digits + 2L else digits + 3L
    out <- cli::ansi_align(out, width = width, align = "right")
  }

  if (!is.null(dims)) {
    out <- matrix(out, nrow = dims[1L], ncol = dims[2L], dimnames = dn)
  }
  out
}

# Render a numeric matrix as decimal-aligned console lines. `col_roles` tags each column:
# "loading" columns are split into vertically stacked blocks when the table is wider than the
# console (or `max_factors_per_block` is set), with the h2/u2 columns repeated in every block;
# non-loading roles are auxiliary (kept together in one block, never split): "corr"
# (factor-correlation and variance tables) is rendered plain, and "compare" (COMPARE
# difference tables) colours cells whose |x| exceeds `cutoff` red. Individual rows are never
# wrapped. Salient loadings (|x| >= cutoff) are bold, weaker ones grey, and Heywood-relevant
# cells (|loading| > 1, h2 > 1, u2 < 0) red; styling is dropped when `color = FALSE` or colours
# are off. With `lower_only = TRUE` the strictly-upper triangle is left blank (for symmetric
# matrices such as factor intercorrelations).
.efa_format_matrix <- function(values, row_labels, col_labels, col_roles,
                               cutoff = 0, digits = 3, color = TRUE,
                               max_factors_per_block = NULL, lower_only = FALSE) {

  values <- as.matrix(values)
  n_col <- ncol(values)
  loading_cols <- which(col_roles == "loading")
  aux_cols <- which(col_roles != "loading")

  cell_str <- .efa_num(values, digits = digits, print_zero = FALSE, pad = FALSE)
  if (isTRUE(lower_only)) {
    cell_str[upper.tri(cell_str)] <- ""
  }

  # Display width (not code-point count) so the column maths matches cli::ansi_align below.
  row_width <- max(cli::ansi_nchar(row_labels, type = "width"), 1L)
  col_widths <- vapply(seq_len(n_col), function(j) {
    max(cli::ansi_nchar(col_labels[j], type = "width"),
        max(cli::ansi_nchar(cell_str[, j], type = "width")))
  }, integer(1L))

  blocks <- .efa_matrix_blocks(loading_cols, aux_cols, row_width, col_widths,
                               max_factors_per_block)

  out <- character(0)
  multi_block <- length(blocks) > 1L

  for (bb in seq_along(blocks)) {
    cols <- blocks[[bb]]

    if (multi_block) {
      if (bb > 1L) {
        out <- c(out, "")
      }
      out <- c(out, .efa_matrix_block_label(col_labels, loading_cols, cols,
                                            bb, length(blocks)))
    }

    out <- c(out, .efa_matrix_block(
      cell_str = cell_str,
      values = values,
      col_roles = col_roles,
      col_labels = col_labels,
      row_labels = row_labels,
      cols = cols,
      row_width = row_width,
      col_widths = col_widths,
      cutoff = cutoff,
      color = color
    ))
  }

  # Drop trailing blank padding (e.g. the upper triangle under `lower_only`, or a centred
  # last header cell) so the lines carry no trailing whitespace. Styled cells end in an ANSI
  # reset, so this only ever trims plain padding.
  sub("[ ]+$", "", out)
}

# Assemble the header and body lines for one column block of `.efa_format_matrix()`.
.efa_matrix_block <- function(cell_str, values, col_roles, col_labels, row_labels,
                              cols, row_width, col_widths, cutoff, color) {

  header_cells <- vapply(cols, function(j) {
    # Style the label first, then centre it, so the alignment padding stays outside the
    # ANSI wrap and the trailing-whitespace trim in `.efa_format_matrix()` can reach it.
    label <- if (isTRUE(color)) cli::style_bold(col_labels[j]) else col_labels[j]
    cli::ansi_align(label, width = col_widths[j], align = "center")
  }, character(1L))
  header <- paste0(strrep(" ", row_width), "  ",
                   paste(header_cells, collapse = "  "))

  rows <- vapply(seq_len(nrow(cell_str)), function(ii) {
    cells <- vapply(cols, function(jj) {
      plain <- cli::ansi_align(cell_str[ii, jj], width = col_widths[jj],
                               align = "right")
      .efa_style_cell(plain, values[ii, jj], col_roles[jj], cutoff, color)
    }, character(1L))
    label <- cli::ansi_align(row_labels[ii], width = row_width, align = "left")
    paste0(label, "  ", paste(cells, collapse = "  "))
  }, character(1L))

  c(header, rows)
}

# Style a single padded cell by column role (no-op when colour is off).
.efa_style_cell <- function(cell, value, role, cutoff, color) {
  if (!isTRUE(color)) {
    return(cell)
  }

  if (identical(role, "loading")) {
    if (is.na(value)) {
      return(cell)
    }
    if (abs(value) > 1) {
      return(cli::style_bold(cli::col_red(cell)))
    }
    if (abs(value) < cutoff) {
      return(cli::col_grey(cell))
    }
    return(cli::style_bold(cell))
  }

  if (identical(role, "h2") && is.finite(value) && value > 1) {
    return(cli::col_red(cell))
  }
  if (identical(role, "u2") && is.finite(value) && value < 0) {
    return(cli::col_red(cell))
  }
  # "compare" columns (COMPARE difference tables) flag salient differences: cells whose
  # absolute value exceeds `cutoff` (the COMPARE `range_red` threshold) are coloured red.
  if (identical(role, "compare") && is.finite(value) && abs(value) > cutoff) {
    return(cli::col_red(cell))
  }

  cell
}

# Split the loading columns into console-fitting blocks; the h2/u2 (aux) columns repeat in
# every block. With `max_factors_per_block` the split is fixed; otherwise it follows the
# console width.
.efa_matrix_blocks <- function(loading_cols, aux_cols, row_width, col_widths,
                               max_factors_per_block = NULL) {
  if (length(loading_cols) <= 1L) {
    return(list(c(loading_cols, aux_cols)))
  }

  if (!is.null(max_factors_per_block)) {
    split_points <- ceiling(seq_along(loading_cols) / max_factors_per_block)
    return(lapply(split(loading_cols, split_points),
                  function(cols) c(cols, aux_cols)))
  }

  console_width <- .efa_console_width()
  full_width <- .efa_table_width(c(loading_cols, aux_cols), row_width, col_widths)
  if (full_width <= console_width) {
    return(list(c(loading_cols, aux_cols)))
  }

  blocks <- list()
  current <- integer(0)

  for (col in loading_cols) {
    candidate <- c(current, col, aux_cols)
    if (length(current) > 0L &&
        .efa_table_width(candidate, row_width, col_widths) > console_width) {
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

.efa_console_width <- function() {
  width <- tryCatch(cli::console_width(), error = function(e) NA_integer_)
  if (!is.finite(width) || width < 40L) {
    return(80L)
  }
  as.integer(width)
}

# Printed width of a block: row label + two spaces + value columns + inter-column gaps.
.efa_table_width <- function(cols, row_width, col_widths) {
  row_width + 2L + sum(col_widths[cols]) + 2L * (length(cols) - 1L)
}

.efa_matrix_block_label <- function(col_labels, loading_cols, cols, block, n_blocks) {
  factor_cols <- intersect(cols, loading_cols)
  first <- col_labels[min(factor_cols)]
  last <- col_labels[max(factor_cols)]

  if (identical(first, last)) {
    paste0("Factors ", first, " (block ", block, "/", n_blocks, ")")
  } else {
    paste0("Factors ", first, "-", last, " (block ", block, "/", n_blocks, ")")
  }
}

# Render a numeric matrix through the shared renderer and print it (with a trailing blank
# line), deriving row/column labels from its dimnames. `role` fills every column: "corr"
# (factor-correlation/variance tables) is plain; "compare" (COMPARE difference tables)
# colours cells with |x| > `cutoff` red. `lower_only` blanks the strictly-upper triangle
# (symmetric matrices); `header = FALSE` drops the column-header row (e.g. an unlabelled
# single-column vector difference). This is the one consumer-facing wrapper around
# `.efa_format_matrix()` for non-loading tables.
.print_efa_matrix <- function(values, role = "corr", digits = 3, cutoff = 0,
                              lower_only = FALSE, header = TRUE) {
  values <- as.matrix(values)

  lines <- .efa_format_matrix(
    values = values,
    row_labels = .efa_variable_names(values),
    col_labels = .efa_factor_names(values),
    col_roles = rep(role, ncol(values)),
    cutoff = cutoff,
    digits = digits,
    lower_only = lower_only
  )

  if (isFALSE(header)) {
    lines <- lines[-1L]
  }

  cat(lines, sep = "\n")
  cat("\n")

  invisible(NULL)
}

.decimals <- function(x) {

  if ((is.null(dim(x)) && !(inherits(x, c("numeric", "integer")))) ||
      (!is.null(dim(x)) && !(inherits(x, c("matrix", "loadings", "LOADINGS",
                                          "SLLOADINGS"))))) {
    cli::cli_abort("{.arg x} must be a numeric vector or matrix, not {.cls {class(x)}}.",
                   class = "efa_not_numeric")
  }

  if (!is.null(dim(x))) {

    max(apply(x, 1:2, function(ll) {
      if (abs(ll - round(ll)) > .Machine$double.eps^0.5) {
        nchar(strsplit(sub('0+$', '', as.character(ll)), ".",
                       fixed = TRUE)[[1]][[2]])
      } else {
        return(0)
      }
    }))

  } else if (length(x) > 1) {

    max(sapply(x, function(ll) {
      if (abs(ll - round(ll)) > .Machine$double.eps^0.5) {
        nchar(strsplit(sub('0+$', '', as.character(ll)), ".", fixed = TRUE)[[1]][[2]])
      } else {
        return(0)
      }
    }))


  } else {

    if (abs(x - round(x)) > .Machine$double.eps^0.5) {
      nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
    } else {
      return(0)
    }

  }

}

# Create a string (used to print settings in print functions, to use in cat()
# with sep = "")
.settings_string <- function(x){

  n <- length(x)

if(n == 1){

  c(.efa_style(x, "bold"))

} else if (n == 2){

  c(.efa_style(x[1], "bold"), " and ", .efa_style(x[2], "bold"))

} else if (n > 2){

  c(paste(.efa_style(x[seq_len(n-1)], "bold"), collapse = ", "), ", and ",
    .efa_style(x[n], "bold"))

}

}


.paste_gof_ci <- function(indices, index) {
  paste0(" [", .numformat(indices$lower[index],
                         pad = FALSE), ", ",
         .numformat(indices$upper[index],
                           pad = FALSE),
         "]")
}
