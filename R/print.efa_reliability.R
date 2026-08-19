#' Print and format a reliability object
#'
#' `print()` shows the reliability coefficients for the general factor and the
#' group factors, for a single group or for each group: the reliability
#' coefficients (omega total, omega hierarchical, and omega subscale, standardized
#' Cronbach's alpha, and the H index) and the common-variance indices (the explained
#' common variance, ECV, and the percent of uncontaminated correlations, PUC).
#' `format()` assembles the same report and returns it as a character vector;
#' `print()` is `cat(format(x), sep = "\n")`. The lines follow the active console
#' theme, so they are plain when colours are disabled (for example when captured
#' into a file or stripped with [cli::ansi_strip()]).
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
  settings <- attr(x, "settings")
  groups <- unique(x$group)
  # A single unnamed group (group = NA) prints without a group header; any named or
  # multigroup result labels each group.
  show_group <- !(length(groups) == 1L && is.na(groups))

  # Reconstruct the factor-by-coefficient matrix for one kind (reliability /
  # common-variance) from the group's long rows, so the shared matrix renderer can lay it
  # out as OMEGA does. Only factors and coefficients that carry a value for this kind
  # appear (e.g. ECV/PUC contribute a general-factor row only), so no all-NA rows are
  # rendered.
  build_mat <- function(sub, kind_sel) {
    keep <- reg[reg$kind == kind_sel & reg$coefficient %in% sub$coefficient, ,
                drop = FALSE]
    if (nrow(keep) == 0L) return(NULL)
    sk <- sub[sub$coefficient %in% keep$coefficient, , drop = FALSE]
    # Rows are keyed on the level and the factor together, and filled by index rather than
    # by label. The label alone does not identify a row: the whole-scale row carries one of
    # its own ("g" for a hierarchy, "total" for a correlated-factors solution) and a group
    # factor may be named the same, in which case indexing by label would write both into
    # one row and print the group factor's value as the whole scale's. The rendered labels
    # are still the factor names, so a collision shows as two rows of one name rather than
    # as a silently dropped value.
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

    if (!is.null(settings) && !is.null(settings$variance)) {
      vtext <- switch(settings$variance,
                      correlation = "the correlation matrix",
                      sums_load = "the model-implied composite variance",
                      settings$variance)
      cli::cli_text("")
      cli::cli_text("Total variance from {vtext}.")
    }

    # A correlated-factors solution has no general factor, so relate the two omega
    # columns rather than print them side by side unexplained. Which relation holds is a
    # property of the solution: a factor's composite also carries true score variance
    # from the other factors, through its cross-loadings and through any factor
    # correlations, so the two columns differ -- in either direction, since the total is
    # a quadratic form whose cross terms take the sign of the other columns' sums times
    # their correlations. Only under exact simple structure and no factor correlations
    # does the composite draw on nothing else and the two coincide. Decide from the
    # values as printed, so the sentence always describes the table below it. Only shown
    # when both columns are there (a coefficients = subset may drop one).
    if (isTRUE(settings$no_general) &&
        all(c("omega_total", "omega_subscale") %in% x$coefficient)) {
      gr <- x[x$level == "group", , drop = FALSE]
      tot <- gr[gr$coefficient == "omega_total", , drop = FALSE]
      sub_om <- gr[gr$coefficient == "omega_subscale", , drop = FALSE]
      # Paired on the group and the factor rather than on position: an undefined
      # coefficient is dropped row by row, so two equally long columns need not describe
      # the same factors, and a multigroup result repeats every factor label once per
      # group -- pairing on the label alone would match only each factor's first group and
      # leave the counts below unequal, so no multigroup table could ever be found to
      # coincide. Only the factors carrying both can be compared, and a factor carrying one
      # of them alone cannot establish equality, so it takes the general wording. Compared
      # through the printer's own formatter rather than against a tolerance, so "equals"
      # says the two columns show the same number at this `digits` -- a tolerance of one
      # last decimal calls values equal that straddle a rounding boundary and print apart.
      key <- function(d) paste(d$group, d$factor, sep = "\r")
      both <- intersect(key(tot), key(sub_om))
      shown <- function(d) .efa_num(d$value[match(both, key(d))], digits = digits,
                                    pad = FALSE)
      coincide <- length(both) > 0 &&
        length(both) == nrow(tot) && length(both) == nrow(sub_om) &&
        identical(shown(tot), shown(sub_om))
      cli::cli_text("")
      if (isTRUE(coincide)) {
        cli::cli_text("Correlated-factors solution: with no general factor, each factor's subscale omega equals its total omega.")
      } else {
        cli::cli_text("Correlated-factors solution: a factor's total omega counts the true score variance its composite receives from every factor, through its cross-loadings and any factor correlations; its subscale omega counts only that factor's own contribution.")
      }
    }

    for (grp in groups) {
      sub <- x[x$group %in% grp, , drop = FALSE]

      if (isTRUE(show_group)) {
        cli::cli_text("")
        cli::cli_text("Group {.emph {grp}}:")
      }

      # The tidy result carries no NA values (the builder drops them), so any NA hole in
      # the reconstructed factor-by-coefficient grid is a coefficient that is undefined for
      # that factor (e.g. omega subscale / H on the whole-scale "total" row of a
      # correlated-factors solution). Print those blank rather than as "NA", as OMEGA's
      # printer does.
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
