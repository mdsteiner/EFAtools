#' Print and format an efa_bartlett object
#'
#' `print()` reports the outcome of [efa_bartlett()]'s test of sphericity: a verdict
#' on whether the test was significant (and what that implies for the suitability
#' of the data for factor analysis), followed by the chi-square statistic, its
#' degrees of freedom, and the p-value. `format()` assembles the same report and
#' returns it as a character vector; `print()` is `cat(format(x), sep = "\n")`.
#' The lines follow the active console theme, so they are plain when colours are
#' disabled (for example when captured into a file or stripped with
#' [cli::ansi_strip()]).
#'
#' @param x An object of class `efa_bartlett` (output from [efa_bartlett()]).
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines (styled to the active console theme;
#'   plain when colours are disabled).
#'
#' @export
#'
#' @method print efa_bartlett
#'
#' @examples
#' bart <- efa_bartlett(test_models$baseline$cormat, N = 500)
#' bart
#'
#' # format() returns the same lines as plain text:
#' writeLines(format(bart))
#'
print.efa_bartlett <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.efa_bartlett
#' @export
#' @method format efa_bartlett
format.efa_bartlett <- function(x, ...) {
  pval <- x$p_value

  cli::cli_format_method({
    cli::cli_text("")

    if (!is.null(pval) && !is.na(pval)) {
      if (pval < .05) {
        significant <- cli::col_green(cli::style_bold("significant"))
        cli::cli_alert_success(
          "The {.strong Bartlett's test of sphericity} was {significant} at an alpha level of .05."
        )
        cli::cli_text("These data are probably suitable for factor analysis.")
      } else {
        not_significant <- cli::col_red(cli::style_bold("not significant"))
        cli::cli_alert_danger(
          "The Bartlett's test of sphericity was {not_significant} at an alpha level of .05."
        )
        cli::cli_text("These data are probably not suitable for factor analysis.")
      }
    } else {
      cli::cli_alert_warning("The Bartlett's test of sphericity did not render a result.")
    }

    # The chi-square line is shown for every object, including the "no result" case above, so
    # render the statistic and p-value defensively when they are NA/NULL (.screen_p_str()
    # handles the p-value tail, shared with the screening report). It is a fixed-format
    # statistic line, emitted verbatim so cli never reflows it mid-token.
    chisq_str <- if (is.null(x$chisq) || is.na(x$chisq)) "NA" else round(x$chisq, 2)
    p_str <- .screen_p_str(pval)

    cli::cli_text("")
    cli::cli_verbatim(paste0("\U1D712\U00B2(", x$df, ") = ", chisq_str, ", ",
                             cli::style_italic("p"), p_str))
  })
}
