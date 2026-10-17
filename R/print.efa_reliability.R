#' Print and format a reliability object
#'
#' `print()` shows the reliability coefficients (omega total, hierarchical, and subscale,
#' standardized alpha, and the H index) of each factor and the common-variance indices (ECV
#' and PUC) of the general factor, for one solution or for each group. `format()` returns the same lines as a character
#' vector. The lines follow the active console theme, and are plain when colours are disabled.
#'
#' @param x An object of class `efa_reliability`.
#' @param digits Integer. The number of decimal places the coefficients are
#'   rounded to. Default is 3.
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines.
#'
#' @family reliability coefficients
#'
#' @export
#' @method print efa_reliability
#'
#' @examples
#' efa_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                    estimator = "PAF", rotation = "promax")
#' rel <- efa_reliability(efa_mod)
#' rel
#'
#' # format() returns the same lines as a character vector:
#' writeLines(format(rel))
#'
print.efa_reliability <- function(x, digits = 3, ...) {
  cat(format(x, digits = digits, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.efa_reliability
#' @export
#' @method format efa_reliability
format.efa_reliability <- function(x, digits = 3, ...) {

  reg <- .reliability_registry()
  groups <- unique(x$group)
  # A single ungrouped solution (group NA) prints without a group header.
  show_group <- !(length(groups) == 1L && is.na(groups))

  # The factor-by-coefficient matrix of one kind of coefficient. Rows are keyed on the level
  # and the factor, because a group factor can have the label of the whole-scale row.
  build_mat <- function(sub, kind_sel) {
    keep <- reg[reg$kind == kind_sel & reg$coefficient %in% sub$coefficient, ,
                drop = FALSE]
    if (nrow(keep) == 0L) return(NULL)
    sk <- sub[sub$coefficient %in% keep$coefficient, , drop = FALSE]
    key <- paste(sk$level, sk$factor, sep = "\r")
    rows <- unique(key)
    facs <- sk$factor[match(rows, key)]
    mat <- matrix(NA_real_, length(rows), nrow(keep),
                  dimnames = list(facs, keep$label))
    for (i in seq_len(nrow(keep))) {
      hit <- sk$coefficient == keep$coefficient[i]
      mat[match(key[hit], rows), keep$label[i]] <- sk$value[hit]
    }
    mat
  }

  cli::cli_format_method({

    for (grp in groups) {
      sub <- x[x$group %in% grp, , drop = FALSE]

      if (isTRUE(show_group)) {
        cli::cli_text("")
        cli::cli_text("Group {.emph {grp}}:")
      }

      # An undefined coefficient prints as a blank cell.
      rel <- build_mat(sub, "reliability")
      if (!is.null(rel)) {
        .print_efa_rule("Reliability coefficients")
        .efa_emit_lines(.efa_corr_lines(rel, digits = digits, na = ""))
      }

      cv <- build_mat(sub, "common_variance")
      if (!is.null(cv)) {
        .print_efa_rule("Common-variance indices")
        .efa_emit_lines(.efa_corr_lines(cv, digits = digits, na = ""))
      }
    }

  })
}
