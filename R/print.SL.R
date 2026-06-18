#' Print and format an SL object
#'
#' `print()` shows a summarised output of the [SL] function: a model header (when
#' the settings are available), the Schmid-Leiman loading matrix, and the
#' variances accounted for. `format()` assembles the same report and returns it as
#' a character vector; `print()` is `cat(format(x), sep = "\n")`. The lines follow
#' the active console theme, so they are plain when colours are disabled (for
#' example when captured into a file or stripped with [cli::ansi_strip()]).
#'
#' @param x An object of class `SL` (output from [SL()]).
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines (styled to the active console theme;
#'   plain when colours are disabled).
#'
#' @export
#' @method print SL
#'
#' @examples
#' EFA_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                type = "EFAtools", method = "PAF", rotation = "promax")
#' sl_mod <- SL(EFA_mod, type = "EFAtools", method = "PAF")
#' sl_mod
#'
#' # format() returns the same lines as plain text:
#' writeLines(format(sl_mod))
#'
print.SL <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.SL
#' @export
#' @method format SL
format.SL <- function(x, ...) {

  cli::cli_format_method({

    if (!all(is.na(x$settings))) {
      # extract the settings for EFA not depending on the method
      method <- x$settings$method
      type <- x$settings$type

      # Emitted verbatim (not via cli_text) so the "setting = 'value'" tokens are never split
      # across a line break; the bold values therefore use style_bold() rather than {.strong}.
      cli::cli_text("")
      cli::cli_verbatim(paste0(
        "EFA for second-order loadings performed with type = '",
        cli::style_bold(type), "' and method = '", cli::style_bold(method), "'"
      ))

      if (!is.null(x$settings$max_iter) && x$iter >= x$settings$max_iter) {
        cli::cli_text("")
        cli::cli_alert_danger(
          "Maximum number of iterations reached without convergence"
        )
      }
    }

    # print the loadings and the variances
    .print_efa_rule("Schmid-Leiman Solution")
    .efa_emit_lines(.efa_capture_loadings(x$sl))

    .print_efa_rule("Variances Accounted for")
    .efa_emit_lines(.efa_corr_lines(x$vars_accounted, digits = 3))
  })
}
