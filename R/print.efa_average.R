#' Print and format an efa_average object
#'
#' `print()` shows a summarised output of the [efa_average()] function: the
#' averaging settings, the error/convergence/Heywood/admissibility rates, the
#' indicator-to-factor correspondences, the averaged loadings (and, for oblique
#' solutions, the factor intercorrelations), the variances accounted for, and the
#' model fit. `format()` assembles the same report and returns it as a character
#' vector; `print()` is `cat(format(x), sep = "\n")`. The lines follow the active
#' console theme, so they are plain when colours are disabled (for example when
#' captured into a file or stripped with [cli::ansi_strip()]).
#'
#' @param x An object of class `efa_average` (output from [efa_average()]).
#' @param stat character. A vector with the statistics to print. Possible inputs
#' are "average", "sd", "range", "min", and "max". Default is "average" and
#' "range".
#' @param plot logical. Whether a plot of the average and min- max loadings should
#' be created. Default is FALSE. If more than 10 factors are extracted, no plot is
#' created. Only used by `print()`.
#' @param ...  Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines.
#'
#' @export
#'
#' @method print efa_average
#'
#' @examples
#' \donttest{
#' EFA_aver <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500)
#' EFA_aver
#'
#' # format() returns the same lines as a character vector:
#' writeLines(format(EFA_aver))
#' }
print.efa_average <- function(x, stat = c("average", "range"),
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

#' @rdname print.efa_average
#' @export
#' @method format efa_average
format.efa_average <- function(x, stat = c("average", "range"), ...) {

  checkmate::assert_subset(stat, c("average", "sd", "range", "min", "max"),
                           empty.ok = FALSE)

  # extract settings
  settings <- x$settings
  estimator <- settings$estimator
  N <- settings$N
  grid <- x$implementations_grid
  averaging <- settings$averaging

  # settings that were varied: keep only the grid columns that correspond to EFA
  # settings. The per-model outcome columns appended to the grid (errors,
  # convergence, Heywood/admissibility flags, fit indices) are never part of the
  # settings list, so this excludes them without a column to maintain.
  varied_settings <- grid[, intersect(names(grid), names(settings)), drop = FALSE]
  # vapply keeps one count per column: apply() would simplify to a matrix when every
  # column has the same number of distinct values, and the count would then run over
  # the cells instead of the columns.
  n_unique <- vapply(varied_settings, function(z) length(unique(z[!is.na(z)])),
                     integer(1))
  varied_settings <- names(n_unique[n_unique > 1])
  # The grid and settings keep the historical key `P_type`; the argument is `p_type`, so
  # display it under that name.
  varied_settings[varied_settings == "P_type"] <- "p_type"

  # quantities embedded in the averaging sentence (the rates are built by
  # .efa_emit_average_rates() from the grid alone)
  no_efas <- nrow(grid)
  averaging_method <- if (averaging == "median") {
    "median"
  } else {
    paste0("mean (trim = ", settings$trim, ")")
  }

  cli::cli_format_method({

    # Prose sentences wrap to the console width; the varied settings collapse into a
    # comma-and-"and" list automatically (each value emphasised individually).
    cli::cli_text("")
    cli::cli_text("Averaging performed with averaging method {.strong {averaging_method}} across {.strong {no_efas}} EFAs, varying the following settings: {.strong {varied_settings}}.")

    cli::cli_text("")
    .efa_emit_average_rates(grid)

    # If no solutions were achieved across which averaging could be performed, flag it and stop
    # after the summary (there are no loadings or fit indices to show). Otherwise show the full
    # report.
    if (.efa_average_no_solutions(grid)) {

      cli::cli_text("")
      cli::cli_alert_warning("No solutions were achieved across which averaging was possible. Best try again with a different number of factors.",
                             wrap = TRUE)

    } else {

      fit <- x$fit_indices
      rownames(fit) <- fit$index

      # Indicator-to-factor correspondences
      .efa_section_rule("Indicator-to-Factor Correspondences")
      salience_threshold <- settings$salience_threshold
      cli::cli_text("")
      cli::cli_text("For each cell, the proportion of solutions including the respective indicator-to-factor correspondence. A salience threshold of {.strong {salience_threshold}} was used to determine indicator-to-factor correspondences.")
      cli::cli_text("")
      .efa_emit_lines(format(
        structure(x$ind_fac_corres, class = c("efa_loadings", "LOADINGS")),
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
        cli::cli_alert_warning("The model is just identified (df = 0). Goodness of fit indices may not be interpretable.",
                               wrap = TRUE)
      }

      .efa_section_rule("Model Fit")
      cli::cli_text("")
      # `averaging` is exactly "mean" or "median", so it names the aggregation the block
      # actually reports (the header below abbreviates the same choice as M or Md).
      cli::cli_text("The fit indices are the {averaging} of the per-solution fit indices, not the fit of the averaged loadings printed above, which are a cell-wise summary rather than a fitted solution.")

      # Effective n per index class. .extract_data() records a solution's fit indices only
      # once it converged without a Heywood case, and leaves the chi-square-based ones NA
      # for PAF throughout, so the two classes can be averaged over different numbers of
      # solutions. Each class is counted on the column that defines it (`chisq` for the
      # chi-square-based indices, `caf` for the residual-based ones), which is exactly the
      # non-missing count the averaging itself used.
      n_resid <- sum(!is.na(grid$caf))
      chisq_based <- !(all(estimator == "PAF") || is.na(N))
      cli::cli_text("")
      if (chisq_based) {
        n_chisq <- sum(!is.na(grid$chisq))
        cli::cli_text("Chi-square-based indices averaged over {n_chisq} of {no_efas} solution{?s}.")
      }
      cli::cli_text("CAF, RMSR, and SRMR averaged over {n_resid} of {no_efas} solution{?s}.")
      cli::cli_text("")
      cli::cli_verbatim(paste0("       ", if (averaging == "mean") "M" else "Md",
                               " (SD) [Min; Max]"))

      if (!chisq_based) {

        lines <- .gof_lines(fit, ind = c("caf", "rmsr", "srmr"),
                            ind_name = c("CAF:  ", "RMSR: ", "SRMR: "),
                            print_zero = c(FALSE, FALSE, FALSE),
                            digits = c(2, 2, 2))
        lines <- c(lines, paste0("df: ",
          .efa_num(fit["df", "average"], 0, print_zero = TRUE)))
        cli::cli_verbatim(lines)

      } else {

        chisq <- .gof_lines(fit, ind = "chisq", ind_name = "\U03C7\U00B2: ",
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

# Emit the error/convergence/Heywood/admissibility summary for an implementations grid.
#
# Each conditional rate has its own denominator, because .extract_data() records an outcome
# only once the preceding stage succeeded: `converged` is set for every solution that did not
# error, while `heywood` and `admissible` are set only for those that then converged (the
# unassessed entries stay NA and drop out of the rate). Name the denominator each rate is
# conditioned on, and drop a clause when its denominator is empty: mean(<all NA>, na.rm =
# TRUE) is NaN, and a rate over no solutions is undefined rather than 0%. Each rate is
# therefore computed on the branch that reports it, which leaves an empty denominator no
# variable to sit in.
.efa_emit_average_rates <- function(grid) {

  pct <- function(x) paste0(round(mean(x, na.rm = TRUE) * 100), "%")

  error_pct <- pct(grid$errors)
  clauses <- "The error rate is at {.strong {error_pct}}."

  # `%in% TRUE` keeps an NA flag from propagating into the condition, as in
  # .efa_average_no_solutions() above.
  if (all(grid$errors %in% TRUE)) {

    clauses <- c(clauses, "No solution could be fitted, so convergence, Heywood cases, and admissibility could not be assessed.")

  } else {

    converged_pct <- pct(grid$converged == 0)
    clauses <- c(clauses, "Of the solutions that did not result in an error, {.strong {converged_pct}} converged.")

    if (!any(grid$converged == 0, na.rm = TRUE)) {

      clauses <- c(clauses, "No solution converged, so Heywood cases and admissibility could not be assessed.")

    } else {

      heywood_pct <- pct(grid$heywood)
      admissible_pct <- pct(grid$admissible)
      clauses <- c(clauses, "Of the solutions that converged, {.strong {heywood_pct}} contained Heywood cases and {.strong {admissible_pct}} were admissible.")

    }

  }

  # A single cli_text() so the clauses wrap as one paragraph rather than one block each.
  cli::cli_text(paste(clauses, collapse = " "))

  invisible(NULL)
}

# Emit a double-ruled (==) header for a major efa_average section: a blank line followed by
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
      .efa_emit_lines(format(x$loadings[[stat_key]]))
    } else if (what == "Phi") {
      .efa_emit_lines(.efa_corr_lines(x$Phi[[stat_key]], lower_only = TRUE))
    } else if (what == "loadings") {
      .efa_emit_lines(.efa_corr_lines(x$loadings[[stat_key]]))
    } else {
      .efa_emit_lines(.efa_corr_lines(x$vars_accounted[[stat_key]]))
    }
  }

  # Section title per statistic, in the order they are printed (which is fixed here, not
  # taken from the order of `stat`). The range is labelled with its definition, since the
  # cells are widths, not the [Min; Max] intervals the Model Fit block reports.
  titles <- c(average = if (averaging == "mean") "Mean" else "Median",
              sd = "Standard Deviation",
              range = "Range (max \u2212 min)",
              min = "Minimum",
              max = "Maximum")

  for (stat_key in intersect(names(titles), stat)) {
    .print_efa_rule(titles[[stat_key]])
    emit(stat_key)
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
