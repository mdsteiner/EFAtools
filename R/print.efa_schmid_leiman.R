#' Print and format an efa_schmid_leiman object
#'
#' `print()` shows a summarised output of the [efa_schmid_leiman()] function: a
#' model header (when the settings are available), the Schmid-Leiman loading
#' matrix, and the variances accounted for. `format()` assembles the same report
#' and returns it as a character vector; `print()` is `cat(format(x), sep = "\n")`.
#' The lines follow the active console theme, so they are plain when colours are
#' disabled (for example when captured into a file or stripped with
#' [cli::ansi_strip()]).
#'
#' @param x An object of class `efa_schmid_leiman` (output from
#'   [efa_schmid_leiman()]).
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines (styled to the active console theme;
#'   plain when colours are disabled).
#'
#' @export
#' @method print efa_schmid_leiman
#'
#' @examples
#' EFA_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                    method = "PAF", rotation = "promax")
#' sl_mod <- efa_schmid_leiman(EFA_mod, method = "PAF")
#' sl_mod
#'
#' # format() returns the same lines as plain text:
#' writeLines(format(sl_mod))
#'
print.efa_schmid_leiman <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.efa_schmid_leiman
#' @export
#' @method format efa_schmid_leiman
format.efa_schmid_leiman <- function(x, ...) {

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

      # Surface non-convergence of the second-order EFA, keyed on its convergence
      # code (with the iteration-count fallback), matching print.efa().
      if (.efa_iteration_nonconvergence(list(convergence = x$convergence,
                                             iter = x$iter,
                                             max_iter = x$settings$max_iter))) {
        cli::cli_text("")
        cli::cli_alert_danger(.efa_nonconvergence_banner(method))
      }
    }

    # print the loadings and the variances
    .print_efa_rule("Schmid-Leiman Solution")
    .efa_emit_lines(.efa_capture_loadings(x$sl))

    .print_efa_rule("Variances Accounted for")
    .efa_emit_lines(.efa_corr_lines(x$vars_accounted, digits = 3))
  })
}
