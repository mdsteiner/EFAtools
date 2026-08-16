#' Print an efa_sl_loadings object
#'
#' @details
#' Prints a Schmid-Leiman loading matrix (general factor, group factors, and the
#' communality/uniqueness columns) as a styled, decimal-aligned table. Loadings with
#' absolute value greater than or equal to `cutoff` are emphasised, smaller loadings are
#' de-emphasised, and Heywood-relevant cells (a loading or communality above 1, or a
#' negative uniqueness) are highlighted. If the matrix has many columns or the console is
#' narrow, the table is split into stacked column blocks so the output stays readable.
#'
#' @param x class efa_sl_loadings matrix.
#' @param cutoff numeric. The value at or above which loadings are emphasized
#'  (default is .2). The default is lower than the .3 of an ordinary loading
#'  table ([print.efa_loadings()]): the group-factor loadings are residualized,
#'  that is, they carry only the variance left once the general factor has been
#'  partialled out, and are therefore smaller than the corresponding first-order
#'  loadings.
#' @param digits numeric. Passed to \code{\link[base:Round]{round}}. Number of digits
#'  to round the loadings to (default is 3).
#' @param max_name_length numeric. The maximum length of the variable names to
#'  display; see [print.efa_loadings()].
#' @param color logical. Whether to apply console styling using \pkg{cli}.
#'  Default is `TRUE`.
#' @param name_style character. How to shorten variable names longer than
#'  `max_name_length`; see [print.efa_loadings()].
#' @param max_factors_per_block numeric or `NULL`. Maximum number of factor
#'  columns to print per block. If `NULL`, the number is chosen from the
#'  console width.
#' @param sort_loadings character. Optional row sorting; see
#'  [print.efa_loadings()]. The default `"none"` keeps the input order. When
#'  sorting is requested, rows are grouped by their largest *group*-factor
#'  loading: the general factor is left out of the comparison, since it is the
#'  largest loading of almost every item and sorting on it would leave the order
#'  untouched.
#' @param ... additional arguments passed to print or format.
#'
#' @returns `print()` returns its argument `x` invisibly; it is
#'   `cat(format(x, ...), sep = "\n")` followed by a blank line for console
#'   spacing. `format()` returns a character vector with the table lines (styled
#'   to the active console theme; plain when colours are disabled).
#'
#' @method print efa_sl_loadings
#' @export
#'
#' @examples
#' EFA_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                    estimator = "PAF", rotation = "promax")
#' efa_schmid_leiman(EFA_mod, estimator = "PAF")
#'
print.efa_sl_loadings <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  # a trailing blank line, so a loading table printed at the console is set off from
  # whatever follows it
  cat("\n")

  invisible(x)
}

#' @rdname print.efa_sl_loadings
#' @method format efa_sl_loadings
#' @export
format.efa_sl_loadings <- function(x, cutoff = .2, digits = 3, max_name_length = 10,
                                   color = TRUE,
                                   name_style = c("truncate", "abbreviate", "full"),
                                   max_factors_per_block = NULL,
                                   sort_loadings = c("none", "primary", "clustered"),
                                   ...) {

  name_style <- .match_arg_ci(name_style)
  sort_loadings <- .match_arg_ci(sort_loadings)

  mat <- unclass(x)

  # The same validator format.efa_loadings() runs, so the shared display API rejects the
  # same values with the same wording and condition classes in both places. The arguments
  # this class does not offer are passed at their no-op values.
  .validate_loadings_print_args(
    x = mat,
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    h2 = NULL,
    color = color,
    max_factor_name_length = NULL,
    max_factors_per_block = max_factors_per_block,
    legend = FALSE
  )

  n_col <- ncol(mat)

  factor_names <- colnames(mat)
  if (is.null(factor_names)) {
    factor_names <- paste0("F", seq_len(n_col))
  }

  var_names <- rownames(mat)
  if (is.null(var_names)) {
    var_names <- paste0("V", seq_len(nrow(mat)))
  }

  # The columns are [g, F1..Fn, h2, u2]; the last two are communality/uniqueness.
  col_roles <- rep("loading", n_col)
  col_roles[n_col - 1L] <- "h2"
  col_roles[n_col] <- "u2"

  loading_cols <- seq_len(n_col - 2L)

  # Shortened before any reordering, for the same reason as in .loadings_print_spec(): a row
  # label has to be a property of the variable, not of where the sorting puts it.
  var_names <- .shorten_loadings_names(var_names, max_length = max_name_length,
                                       name_style = name_style)

  # Row sorting compares the group factors only. The general factor carries the largest
  # loading of almost every item in a Schmid-Leiman solution, so including it would assign
  # every row to the same primary factor and leave the order as it was.
  group_cols <- setdiff(loading_cols, 1L)
  if (length(group_cols) > 0L) {
    row_order <- .efa_loading_row_order(mat[, group_cols, drop = FALSE], sort_loadings)
    mat <- mat[row_order, , drop = FALSE]
    var_names <- var_names[row_order]
  }

  # Count items (rows) flagged as Heywood cases: a loading with absolute value
  # above 1, a communality (h2) above 1, or a negative uniqueness. These are
  # deterministically coupled within an item, so each affected item counts once.
  n_heywood <- sum(
    rowSums(abs(mat[, loading_cols, drop = FALSE]) > 1, na.rm = TRUE) > 0 |
      (!is.na(mat[, n_col - 1L]) & mat[, n_col - 1L] > 1) |
      (!is.na(mat[, n_col]) & mat[, n_col] < 0))

  cli::cli_format_method({
    # Through .efa_emit_lines(), for the same reason as in format.efa_loadings(): the blank
    # line separating stacked column blocks would not survive cli_verbatim().
    .efa_emit_lines(.efa_format_matrix(
      values = mat,
      row_labels = var_names,
      col_labels = factor_names,
      col_roles = col_roles,
      cutoff = cutoff,
      digits = digits,
      color = color,
      max_factors_per_block = max_factors_per_block
    ))

    if (n_heywood >= 1L) {
      cli::cli_verbatim("")
      if (n_heywood == 1L) {
        cli::cli_alert_warning("Results contain a Heywood case!")
      } else {
        cli::cli_alert_warning("Results contain {n_heywood} Heywood cases!")
      }
    }
  })
}
