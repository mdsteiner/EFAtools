#' Print and format a BARTLETT object
#'
#' `print()` reports the outcome of [BARTLETT()]'s test of sphericity: a verdict
#' on whether the test was significant (and what that implies for the suitability
#' of the data for factor analysis), followed by the chi-square statistic, its
#' degrees of freedom, and the p-value. `format()` assembles the same report and
#' returns it as a character vector; `print()` is `cat(format(x), sep = "\n")`.
#' The lines follow the active console theme, so they are plain when colours are
#' disabled (for example when captured into a file or stripped with
#' [cli::ansi_strip()]).
#'
#' @param x An object of class `BARTLETT` (output from [BARTLETT()]).
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines (styled to the active console theme;
#'   plain when colours are disabled).
#'
#' @export
#'
#' @method print BARTLETT
#'
#' @examples
#' bart <- BARTLETT(test_models$baseline$cormat, N = 500)
#' bart
#'
#' # format() returns the same lines as plain text:
#' writeLines(format(bart))
#'
print.BARTLETT <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.BARTLETT
#' @export
#' @method format BARTLETT
format.BARTLETT <- function(x, ...) {
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
    # render the statistic and p-value defensively when they are NA/NULL. It is a fixed-format
    # statistic line, emitted verbatim so cli never reflows it mid-token.
    chisq_str <- if (is.null(x$chisq) || is.na(x$chisq)) "NA" else round(x$chisq, 2)
    p_str <- if (is.null(pval) || is.na(pval)) {
      " = NA"
    } else if (pval < .001) {
      " < .001"
    } else {
      paste0(" = ", round(pval, 3))
    }

    cli::cli_text("")
    cli::cli_verbatim(paste0("\U1D712\U00B2(", x$df, ") = ", chisq_str, ", ",
                             cli::style_italic("p"), p_str))
  })
}
