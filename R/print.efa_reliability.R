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
    facs <- unique(sk$factor)
    mat <- matrix(NA_real_, length(facs), nrow(keep),
                  dimnames = list(facs, keep$label))
    for (i in seq_len(nrow(keep))) {
      ri <- sk[sk$coefficient == keep$coefficient[i], , drop = FALSE]
      mat[ri$factor, keep$label[i]] <- ri$value
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

    # A correlated-factors solution has no general factor, so a group factor's omega
    # hierarchical is zero and its omega total and omega subscale coincide exactly.
    # Say so, rather than print two identical columns unexplained. Only when both
    # columns are actually shown (a coefficients = subset may drop one).
    if (isTRUE(settings$no_general) &&
        all(c("omega_total", "omega_subscale") %in% x$coefficient)) {
      cli::cli_text("")
      cli::cli_text("Correlated-factors solution: with no general factor, each factor's subscale omega equals its total omega.")
    }

    for (grp in groups) {
      sub <- x[x$group %in% grp, , drop = FALSE]

      if (isTRUE(show_group)) {
        cli::cli_text("")
        cli::cli_text("Group {.emph {grp}}:")
      }

      # The tidy result carries no NA values (the builder drops them), so any NA hole in
      # the reconstructed factor-by-coefficient grid is a coefficient that is undefined for
      # that factor (e.g. omega subscale / H on the whole-scale g row of a correlated-factors
      # solution). Print those blank rather than as "NA", as OMEGA's printer does.
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
