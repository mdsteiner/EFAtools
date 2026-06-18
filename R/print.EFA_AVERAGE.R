#' Print and format an EFA_AVERAGE object
#'
#' `print()` shows a summarised output of the [EFA_AVERAGE] function: the
#' averaging settings, the error/convergence/Heywood/admissibility rates, the
#' indicator-to-factor correspondences, the averaged loadings (and, for oblique
#' solutions, the factor intercorrelations), the variances accounted for, and the
#' model fit. `format()` assembles the same report and returns it as a character
#' vector; `print()` is `cat(format(x), sep = "\n")`. The lines follow the active
#' console theme, so they are plain when colours are disabled (for example when
#' captured into a file or stripped with [cli::ansi_strip()]).
#'
#' @param x An object of class `EFA_AVERAGE` (output from [EFA_AVERAGE()]).
#' @param stat character. A vector with the statistics to print. Possible inputs
#' are "average", "sd", "range", "min", and "max". Default is "average" and
#' "range".
#' @param plot logical. Whether a plot of the average and min- max loadings should
#' be created. Default is FALSE. If more than 10 factors are extracted, no plot is
#' created. Only used by `print()`.
#' @param ...  Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines (styled to the active console theme;
#'   plain when colours are disabled).
#'
#' @export
#'
#' @method print EFA_AVERAGE
#'
#' @examples
#' \dontrun{
#' EFA_aver <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500)
#' EFA_aver
#'
#' # format() returns the same lines as plain text:
#' writeLines(format(EFA_aver))
#' }
print.EFA_AVERAGE <- function(x, stat = c("average", "range"),
                              plot = FALSE, ...) {
  cat(format(x, stat = stat, ...), sep = "\n")

  # Plotting is a genuine side-effect, so it stays in print() (format() is text-only). Skip it
  # when no solution survived: the loadings are all NA and format() has already flagged this.
  if (isTRUE(plot) && !.efa_average_no_solutions(x$implementations_grid)) {

    if (ncol(x$loadings$average) <= 10) {

      print(plot(x))

    } else {

      obj_name <- deparse1(substitute(x))
      cli::cli_inform(c(
        "i" = "The factor solution contained more than 10 factors, no plot was generated.",
        "i" = "If you still want to create the plot, use {.code plot({obj_name})}."
      ))

    }

  }

  invisible(x)
}

#' @rdname print.EFA_AVERAGE
#' @export
#' @method format EFA_AVERAGE
format.EFA_AVERAGE <- function(x, stat = c("average", "range"), ...) {

  checkmate::assert_subset(stat, c("average", "sd", "range", "min", "max"),
                           empty.ok = FALSE)

  # extract settings
  settings <- x$settings
  method <- settings$method
  N <- settings$N
  grid <- x$implementations_grid
  averaging <- settings$averaging

  # settings that were varied: keep only the grid columns that correspond to EFA
  # settings. The per-model outcome columns appended to the grid (errors,
  # convergence, Heywood/admissibility flags, fit indices) are never part of the
  # settings list, so this excludes them without a column to maintain.
  varied_settings <- grid[, intersect(names(grid), names(settings)), drop = FALSE]
  varied_settings <- apply(varied_settings, 2, function(x) unique(x[!is.na(x)]))
  varied_settings <- sapply(varied_settings, length)
  varied_settings <- names(varied_settings[varied_settings > 1])

  # quantities embedded in the summary sentences
  no_efas <- nrow(grid)
  averaging_method <- if (averaging == "median") {
    "median"
  } else {
    paste0("mean (trim = ", settings$trim, ")")
  }
  error_pct <- paste0(round(mean(grid$errors, na.rm = TRUE) * 100), "%")
  converged_pct <- paste0(round(mean(grid$converged == 0, na.rm = TRUE) * 100), "%")
  heywood_pct <- paste0(round(mean(grid$heywood, na.rm = TRUE) * 100), "%")
  admissible_pct <- paste0(round(mean(grid$admissible, na.rm = TRUE) * 100), "%")

  cli::cli_format_method({

    # Prose sentences wrap to the console width; the varied settings collapse into a
    # comma-and-"and" list automatically (each value emphasised individually).
    cli::cli_text("")
    cli::cli_text("Averaging performed with averaging method {.strong {averaging_method}} across {.strong {no_efas}} EFAs, varying the following settings: {.strong {varied_settings}}.")

    cli::cli_text("")
    cli::cli_text("The error rate is at {.strong {error_pct}}. Of the solutions that did not result in an error, {.strong {converged_pct}} converged, {.strong {heywood_pct}} contained Heywood cases, and {.strong {admissible_pct}} were admissible.")

    # If no solutions were achieved across which averaging could be performed, flag it and stop
    # after the summary (there are no loadings or fit indices to show). Otherwise show the full
    # report.
    if (.efa_average_no_solutions(grid)) {

      cli::cli_text("")
      cli::cli_alert_warning("No solutions were achieved across which averaging was possible. Best try again with a different number of factors.")

    } else {

      fit <- x$fit_indices
      rownames(fit) <- fit$index

      # Indicator-to-factor correspondences
      .efa_section_rule("Indicator-to-Factor Correspondences")
      salience_threshold <- settings$salience_threshold
      cli::cli_text("")
      cli::cli_text("For each cell, the proportion of solutions including the respective indicator-to-factor correspondence. A salience threshold of {.strong {salience_threshold}} was used to determine indicator-to-factor correspondences.")
      cli::cli_text("")
      .efa_emit_lines(.efa_capture_loadings(
        structure(x$ind_fac_corres, class = "LOADINGS"),
        cutoff = 1e-4, digits = 2))

      # Loadings
      .efa_section_rule("Loadings")
      .print_average(x, what = "loadings", stat = stat, averaging = averaging)

      # Factor intercorrelations for oblique solutions
      if (!all(is.na(x$Phi))) {
        .efa_section_rule("Factor Intercorrelations from Oblique Solutions")
        .print_average(x, what = "Phi", stat = stat, averaging = averaging)
      }

      # Variances accounted for
      .efa_section_rule("Variances Accounted for")
      .print_average(x, what = "vars_accounted", stat = stat, averaging = averaging)

      # Model fit
      if (fit["df", "average"] == 0) {
        cli::cli_text("")
        cli::cli_alert_warning("The model is just identified (df = 0). Goodness of fit indices may not be interpretable.")
      }

      .efa_section_rule("Model Fit")
      cli::cli_text("")
      cli::cli_verbatim(paste0("       ", if (averaging == "mean") "M" else "Md",
                               " (SD) [Min; Max]"))

      if (all(method == "PAF") || is.na(N)) {

        lines <- .gof_lines(fit, ind = c("caf", "rmsr", "srmr"),
                            ind_name = c("CAF:  ", "RMSR: ", "SRMR: "),
                            print_zero = c(FALSE, FALSE, FALSE),
                            digits = c(2, 2, 2))
        lines <- c(lines, paste0("df: ",
          .efa_num(fit["df", "average"], 0, print_zero = TRUE)))
        cli::cli_verbatim(lines)

      } else {

        chisq <- .gof_lines(fit, ind = "chisq", ind_name = "\U1D712\U00B2: ",
                            print_zero = TRUE, digits = 2)
        df_line <- paste0("df: ",
                          .efa_num(fit["df", "average"], 0, print_zero = TRUE))
        rest <- .gof_lines(fit,
          ind = c("p_chi", "cfi", "tli", "rmsea", "aic", "bic",
                  "ecvi", "caf", "rmsr", "srmr"),
          ind_name = c(.efa_style("p: ", "italic"), "CFI: ", "TLI: ",
                       "RMSEA: ", "AIC: ", "BIC: ", "ECVI: ", "CAF: ",
                       "RMSR: ", "SRMR: "),
          print_zero = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
                         TRUE, FALSE, FALSE, FALSE),
          digits = c(3, 2, 2, 2, 2, 2, 2, 2, 2, 2))
        cli::cli_verbatim(c(chisq, df_line, rest))

      }

    }
  })
}

# TRUE when no fitted solution survived to be averaged (every run errored, did not converge,
# or was a Heywood case), so there are no loadings or fit indices to report. Consumed by both
# format() (to flag it and stop after the summary) and print() (to suppress plotting).
.efa_average_no_solutions <- function(grid) {
  all((grid$converged != 0 | grid$errors | grid$heywood) %in% TRUE)
}

# Emit a double-ruled (==) header for a major EFA_AVERAGE section: a blank line followed by
# the rule, mirroring the original `cat(rule(line = 2))` spacing (the following blank, if any,
# is supplied by the next element). The double rule keeps the major/minor (==/--) visual
# hierarchy, with the minor sub-rules emitted by `.print_efa_rule()`.
.efa_section_rule <- function(title) {
  cli::cli_text("")
  cli::cli_verbatim(cli::rule(left = cli::style_bold(title), line = 2))
  invisible(NULL)
}

.print_average <- function(x, what, stat, averaging) {

  # Emit one statistic's table into the cli report: average/min/max loadings are
  # LOADINGS-classed (styled loading table), captured verbatim with their colour; everything
  # else (sd/range loadings, Phi, variances) is a plain "corr" matrix. Phi is symmetric, so
  # only its lower triangle shows.
  emit <- function(stat_key) {
    if (what == "loadings" && stat_key %in% c("average", "min", "max")) {
      .efa_emit_lines(.efa_capture_loadings(x$loadings[[stat_key]]))
    } else if (what == "Phi") {
      .efa_emit_lines(.efa_corr_lines(x$Phi[[stat_key]], lower_only = TRUE))
    } else if (what == "loadings") {
      .efa_emit_lines(.efa_corr_lines(x$loadings[[stat_key]]))
    } else {
      .efa_emit_lines(.efa_corr_lines(x$vars_accounted[[stat_key]]))
    }
  }

  if ("average" %in% stat) {
    .print_efa_rule(if (averaging == "mean") "Mean" else "Median")
    emit("average")
  }

  if ("sd" %in% stat) {
    .print_efa_rule("Standard Deviation")
    emit("sd")
  }

  if ("range" %in% stat) {
    .print_efa_rule("Range")
    emit("range")
  }

  if ("min" %in% stat) {
    .print_efa_rule("Minimum")
    emit("min")
  }

  if ("max" %in% stat) {
    .print_efa_rule("Maximum")
    emit("max")
  }

  invisible(NULL)
}

# Build the "M (SD) [Min; Max]" report line for each averaged fit index. `.efa_num` drops the
# leading zero of values in (-1, 1) when print_zero is FALSE (and keeps it otherwise) and
# preserves the sign, so it already produces the stripped/unstripped form each index needs.
# The integer-valued indices (print_zero = TRUE, e.g. AIC/BIC) are left-padded to align the
# bracketed range; the decimal-only ones are not.
.gof_lines <- function(fit, ind, ind_name, print_zero, digits) {
  fmt <- function(stat, i) {
    .efa_num(fit[ind[i], stat], digits = digits[i],
             print_zero = print_zero[i], pad = print_zero[i])
  }

  vapply(seq_along(ind), function(i) {
    paste0(ind_name[i], fmt("average", i), " (", fmt("sd", i), ") [",
           fmt("min", i), "; ", fmt("max", i), "]")
  }, character(1L))
}
