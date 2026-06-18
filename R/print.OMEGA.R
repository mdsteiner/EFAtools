#' Print and format an OMEGA object
#'
#' `print()` shows the omega coefficients computed by [OMEGA()]: omega total (and,
#' for multi-factor solutions, omega hierarchical, omega subscale, the H index,
#' the explained common variance, and the percent of uncontaminated correlations)
#' for the general factor and the group factors, for a single group or for each
#' group. `format()` assembles the same report and returns it as a character
#' vector; `print()` is `cat(format(x), sep = "\n")`. The lines follow the active
#' console theme, so they are plain when colours are disabled (for example when
#' captured into a file or stripped with [cli::ansi_strip()]).
#'
#' @param x An object of class `OMEGA` (output from [OMEGA()]).
#' @param digits Integer. The number of decimal places the coefficients are
#'   rounded to (passed to [round()]). Default is 3.
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines (styled to the active console theme;
#'   plain when colours are disabled).
#'
#' @method print OMEGA
#'
#' @export
#'
#' @examples
#' efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                type = "EFAtools", method = "PAF", rotation = "promax")
#' sl_mod <- SL(efa_mod, type = "EFAtools", method = "PAF")
#'
#' om <- OMEGA(sl_mod, type = "EFAtools",
#'             factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2)
#' om
#'
#' # format() returns the same lines as plain text:
#' writeLines(format(om))
#'
print.OMEGA <- function(x, digits = 3, ...) {
  cat(format(x, digits = digits, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.OMEGA
#' @export
#' @method format OMEGA
format.OMEGA <- function(x, digits = 3, ...) {

  cli::cli_format_method({

    # Render an omega-coefficient matrix as before (base R's matrix layout with blank NA
    # cells), emitted verbatim so cli does not reflow the aligned columns.
    emit_matrix <- function(mat) {
      cli::cli_verbatim(utils::capture.output(
        print(round(unclass(mat), digits = digits), na.print = "")))
    }

    if (is.list(x)) {

      # In case of multiple groups
      group_names <- names(x)

      if (length(x[[1]]) == 1) {

        # A single factor, omega total only
        cli::cli_text("{.strong Omega total for the single factor for each group:}")
        for (i in seq_along(group_names)) {
          grp <- group_names[i]
          value <- round(x[[i]], digits = digits)
          cli::cli_text("")
          cli::cli_text("Group {.emph {grp}}: {value}")
        }

      } else if (length(x[[1]]) == 2) {

        # A single factor, omega total and H index
        cli::cli_text("{.strong Omega total and H index for the single factor for each group:}")
        for (i in seq_along(group_names)) {
          grp <- group_names[i]
          omega <- round(x[[i]][1], digits = digits)
          h_index <- round(x[[i]][2], digits = digits)
          cli::cli_text("")
          cli::cli_text("Group {.emph {grp}}:")
          cli::cli_text("Omega: {omega}")
          cli::cli_text("H index: {h_index}")
        }

      } else {

        # The full coefficient matrix
        desc <- if (ncol(x[[1]]) == 3) {
          paste0("Omega total, omega hierarchical, and omega subscale for the ",
                 "general factor (top row) and the group factors for each group:")
        } else {
          paste0("Omega total, omega hierarchical, omega subscale, H index, ",
                 "explained common variance (ECV), and percent of uncontaminated ",
                 "correlations (PUC) for the general factor (top row) and omegas and ",
                 "H index for the group factors for each group:")
        }
        cli::cli_text("{.strong {desc}}")

        for (i in seq_along(group_names)) {
          grp <- group_names[i]
          cli::cli_text("")
          cli::cli_text("Group {.emph {grp}}:")
          emit_matrix(x[[i]])
        }

      }

    } else {

      # In case of a single group
      if (length(x) == 1) {

        value <- round(unclass(x), digits = digits)
        cli::cli_text("{.strong Omega total for the single factor:} {value}")

      } else if (length(x) == 2) {

        omega <- round(unclass(x)[1], digits = digits)
        h_index <- round(unclass(x)[2], digits = digits)
        cli::cli_text("{.strong Omega total and H index for the single factor:}")
        cli::cli_text("")
        cli::cli_text("Omega: {omega}")
        cli::cli_text("H index: {h_index}")

      } else {

        desc <- if (ncol(x) == 3) {
          paste0("Omega total, omega hierarchical, and omega subscale for the ",
                 "general factor (top row) and the group factors:")
        } else {
          paste0("Omega total, omega hierarchical, omega subscale, H index, ",
                 "explained common variance (ECV), and percent of uncontaminated ",
                 "correlations (PUC) for the general factor (top row) and omegas and ",
                 "H index for the group factors:")
        }
        cli::cli_text("{.strong {desc}}")
        cli::cli_text("")
        emit_matrix(x)

      }

    }
  })
}
