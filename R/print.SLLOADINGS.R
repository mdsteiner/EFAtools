#' Print SLLOADINGS object
#'
#' @details
#' Prints a Schmid-Leiman loading matrix (general factor, group factors, and the
#' communality/uniqueness columns) as a styled, decimal-aligned table. Loadings with
#' absolute value greater than or equal to `cutoff` are emphasised, smaller loadings are
#' de-emphasised, and Heywood-relevant cells (a loading or communality above 1, or a
#' negative uniqueness) are highlighted. If the matrix has many columns or the console is
#' narrow, the table is split into stacked column blocks so the output stays readable.
#'
#' @param x class SLLOADINGS matrix.
#' @param cutoff numeric. The number above which to print loadings in bold
#'  (default is .2).
#' @param digits numeric. Passed to \code{\link[base:Round]{round}}. Number of digits
#'  to round the loadings to (default is 3).
#' @param color logical. Whether to apply console styling using \pkg{cli}.
#'  Default is `TRUE`.
#' @param ... additional arguments passed to print or format.
#'
#' @method print SLLOADINGS
#' @export
#'
#' @examples
#' EFA_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                type = "EFAtools", method = "PAF", rotation = "promax")
#' SL(EFA_mod, type = "EFAtools", method = "PAF")
#'
print.SLLOADINGS <- function(x, cutoff = .2, digits = 3, color = TRUE, ...) {

  mat <- unclass(x)
  n_col <- ncol(mat)

  factor_names <- colnames(mat)
  if (is.null(factor_names)) {
    factor_names <- paste0("F", seq_len(n_col))
  }

  var_names <- rownames(mat)
  if (is.null(var_names)) {
    var_names <- paste0("V", seq_len(nrow(mat)))
  }
  var_names <- .shorten_loadings_names(var_names, max_length = 10,
                                       name_style = "truncate")

  # SLLOADINGS columns are [g, F1..Fn, h2, u2]; the last two are communality/uniqueness.
  col_roles <- rep("loading", n_col)
  col_roles[n_col - 1L] <- "h2"
  col_roles[n_col] <- "u2"

  loading_cols <- seq_len(n_col - 2L)
  # Count items (rows) flagged as Heywood cases: a loading with absolute value
  # above 1, a communality (h2) above 1, or a negative uniqueness. These are
  # deterministically coupled within an item, so each affected item counts once.
  n_heywood <- sum(
    rowSums(abs(mat[, loading_cols, drop = FALSE]) > 1, na.rm = TRUE) > 0 |
      (!is.na(mat[, n_col - 1L]) & mat[, n_col - 1L] > 1) |
      (!is.na(mat[, n_col]) & mat[, n_col] < 0))

  cat(cli::cli_format_method({
    cli::cli_verbatim(.efa_format_matrix(
      values = mat,
      row_labels = var_names,
      col_labels = factor_names,
      col_roles = col_roles,
      cutoff = cutoff,
      digits = digits,
      color = color
    ))

    if (n_heywood >= 1L) {
      cli::cli_verbatim("")
      if (n_heywood == 1L) {
        cli::cli_alert_warning("Results contain a Heywood case!")
      } else {
        cli::cli_alert_warning("Results contain {n_heywood} Heywood cases!")
      }
    }
  }), sep = "\n")
  cat("\n")

  invisible(x)
}

#' @rdname print.SLLOADINGS
#' @method format SLLOADINGS
#' @export
format.SLLOADINGS <- function(x, ...) {
  # `print()` ends with a blank line for console spacing, which capture.output
  # records as a trailing empty element; drop it so format() returns only the
  # rendered table lines (plain, un-styled).
  out <- cli::ansi_strip(utils::capture.output(print(x, ...)))
  if (length(out) > 0L && !nzchar(out[length(out)])) out[-length(out)] else out
}
