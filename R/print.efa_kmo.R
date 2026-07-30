#' Print and format an efa_kmo object
#'
#' `print()` shows the Kaiser-Meyer-Olkin (KMO) criterion computed by [efa_kmo()]: a
#' titled section with a verdict on the overall KMO value (and what it implies for
#' the suitability of the data for factor analysis), the overall value, and the
#' per-variable KMO values. `format()` assembles the same report and returns it as
#' a character vector; `print()` is `cat(format(x), sep = "\n")`. The lines follow
#' the active console theme, so they are plain when colours are disabled (for
#' example when captured into a file or stripped with [cli::ansi_strip()]).
#'
#' @param x An object of class `efa_kmo` (output from [efa_kmo()]).
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines (styled to the active console theme;
#'   plain when colours are disabled).
#'
#' @export
#'
#' @method print efa_kmo
#'
#' @examples
#' KMO_base <- efa_kmo(test_models$baseline$cormat)
#' KMO_base
#'
#' # format() returns the same lines as a character vector:
#' writeLines(format(KMO_base))
#'
print.efa_kmo <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.efa_kmo
#' @export
#' @method format efa_kmo
format.efa_kmo <- function(x, ...) {
  KMO <- x$KMO

  cli::cli_format_method({
    cli::cli_text("")
    cli::cli_rule(left = "Kaiser-Meyer-Olkin criterion (KMO)")
    cli::cli_text("")

    if (!is.null(KMO) && !is.na(KMO)) {

      # Kaiser's verbal bands for the KMO value (Kaiser, 1970; Kaiser & Rice, 1974). The band
      # carries its label, its verdict styling (green = suitable, yellow = mediocre, red =
      # unsuitable) and the suitability it implies; .kmo_band() holds the cut-offs, so this
      # report and efa_screen()'s cannot band the same value differently.
      band <- .kmo_band(KMO)

      label <- band$colour(cli::style_bold(band$label))
      band$alert("The overall KMO value for your data is {label}.", wrap = TRUE)
      cli::cli_text("These data are {band$suitability} suitable for factor analysis.")

      cli::cli_text("")
      overall <- cli::style_bold(as.character(round(KMO, 3)))
      cli::cli_text("Overall: {overall}")

      cli::cli_text("")
      cli::cli_text("For each variable:")
      # The per-variable values are kept as base R's named-vector layout, emitted verbatim so
      # cli does not reflow the aligned columns.
      cli::cli_verbatim(utils::capture.output(print(round(x$KMO_i, 3))))

    } else {
      cli::cli_alert_warning("Sorry, the KMO value for your data is not available.",
                             wrap = TRUE)
    }
  })
}
