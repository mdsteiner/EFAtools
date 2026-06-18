#' Print and format a KMO object
#'
#' `print()` shows the Kaiser-Meyer-Olkin (KMO) criterion computed by [KMO()]: a
#' titled section with a verdict on the overall KMO value (and what it implies for
#' the suitability of the data for factor analysis), the overall value, and the
#' per-variable KMO values. `format()` assembles the same report and returns it as
#' a character vector; `print()` is `cat(format(x), sep = "\n")`. The lines follow
#' the active console theme, so they are plain when colours are disabled (for
#' example when captured into a file or stripped with [cli::ansi_strip()]).
#'
#' @param x An object of class `KMO` (output from [KMO()]).
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines (styled to the active console theme;
#'   plain when colours are disabled).
#'
#' @export
#'
#' @method print KMO
#'
#' @examples
#' KMO_base <- KMO(test_models$baseline$cormat)
#' KMO_base
#'
#' # format() returns the same lines as plain text:
#' writeLines(format(KMO_base))
#'
print.KMO <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.KMO
#' @export
#' @method format KMO
format.KMO <- function(x, ...) {
  KMO <- x$KMO

  cli::cli_format_method({
    cli::cli_text("")
    cli::cli_rule(left = "Kaiser-Meyer-Olkin criterion (KMO)")
    cli::cli_text("")

    if (!is.null(KMO) && !is.na(KMO)) {

      # Kaiser's verbal bands for the KMO value. Each band carries its verdict styling (green =
      # suitable, yellow = mediocre, red = unsuitable; Kaiser, 1970; Kaiser & Rice, 1974) and the
      # suitability it implies, so the cutoffs are written once. Ordered high to low; the first
      # band the value clears applies.
      bands <- list(
        list(min = .9,   label = "marvellous",   colour = cli::col_green,
             alert = cli::cli_alert_success, suitability = "probably"),
        list(min = .8,   label = "meritorious",  colour = cli::col_green,
             alert = cli::cli_alert_success, suitability = "probably"),
        list(min = .7,   label = "middling",     colour = cli::col_green,
             alert = cli::cli_alert_success, suitability = "probably"),
        list(min = .6,   label = "mediocre",     colour = cli::col_yellow,
             alert = cli::cli_alert_warning, suitability = "probably"),
        list(min = .5,   label = "miserable",    colour = cli::col_red,
             alert = cli::cli_alert_danger,  suitability = "hardly"),
        list(min = -Inf, label = "unacceptable", colour = cli::col_red,
             alert = cli::cli_alert_danger,  suitability = "not")
      )
      band <- Find(function(b) KMO >= b$min, bands)

      label <- band$colour(cli::style_bold(band$label))
      band$alert("The overall KMO value for your data is {label}.")
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
      cli::cli_alert_warning("Sorry, the KMO value for your data is not available.")
    }
  })
}
