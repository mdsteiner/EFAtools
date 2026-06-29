#' Print and summarise an EFA object
#'
#' `print()` shows a concise overview of an [EFA()] or [EFA_POOLED()] solution:
#' a model header, the loading matrix (with the factor intercorrelations for
#' oblique solutions), the variances accounted for, and the model fit.
#' [summary()] returns a `summary.EFA` object whose print method adds the full
#' diagnostics: model and simple-structure diagnostics, confidence-interval
#' tables, the structure matrix, multiple-imputation uncertainty (for pooled
#' objects), and residual diagnostics. `format()` assembles the same report and
#' returns it as a character vector; `print()` is `cat(format(x), sep = "\n")`.
#' The lines follow the active console theme, so they are plain when colours are
#' disabled (for example when captured into a file or stripped with
#' [cli::ansi_strip()]).
#'
#' @details
#' The methods are shared by single-imputation `EFA` objects and pooled
#' `EFA_POOLED` objects. For `EFA_POOLED` objects the header reports the number
#' of imputations and the alignment/pooling settings; confidence intervals and a
#' multiple-imputation uncertainty summary are shown by [summary()] when the
#' pooled object carries bootstrap/MI quantities.
#'
#' In [summary()], `ci_filter` controls which loading intervals are shown:
#' `"salient"` reports intervals for loadings whose absolute point estimate is at
#' least `cutoff`, `"nonzero"` reports intervals excluding zero, and `"all"`
#' reports every finite interval.
#'
#' @param x,object An object of class `EFA` or `EFA_POOLED`; for the
#'   `summary.EFA` methods, the object returned by [summary()].
#' @param cutoff numeric. The absolute value at or above which loadings are
#'   emphasised in the loading table. Default is .3.
#' @param digits numeric. Number of decimal places for the printed tables.
#'   Default is 3.
#' @param max_name_length numeric. Maximum length of the variable names to
#'   display; longer names are cut from the right.
#' @param sort_loadings character. Optional row sorting for the loading table.
#'   See [print.LOADINGS()].
#' @param show_loading_legend logical. Whether to print a short legend for the
#'   loading-table styling. Default is `TRUE`.
#' @param max_factors_per_block numeric or `NULL`. Maximum number of factor
#'   columns per loading-table block. If `NULL`, chosen from the console width.
#' @param ci character. Which confidence intervals [summary()] shows, if
#'   available. `"auto"` and `"separate"` print CI sections when CIs were
#'   computed; `"none"` suppresses them. Default is `"auto"`.
#' @param ci_filter character. Which loading CIs [summary()] prints: `"salient"`
#'   (default), `"all"`, or `"nonzero"`; see Details.
#' @param diagnostics_top_n numeric. Maximum number of item-level entries
#'   [summary()] prints per simple-structure diagnostic.
#' @param residual_cutoff numeric. Absolute residual cutoff for the residual
#'   diagnostics in [summary()]. Default is .1.
#' @param residual_top_n numeric. Maximum number of residuals [summary()] prints.
#'   Use `Inf` to print all residuals above `residual_cutoff`.
#' @param show_structure logical. Whether [summary()] prints the structure matrix
#'   for oblique solutions when available. Default is `TRUE`.
#' @param cross_loading_cutoff numeric. Cutoff for counting cross-loadings in the
#'   [summary()] diagnostics. Defaults to `cutoff`.
#' @param min_primary_gap numeric. Minimum desired absolute difference between the
#'   largest and second-largest absolute loading of an item, used in the
#'   [summary()] diagnostics.
#' @param min_salient_per_factor numeric. Minimum number of salient indicators per
#'   factor used in the [summary()] diagnostics. Default is 3.
#' @param show_mi_diagnostics logical or `NULL`. Whether [summary()] prints a
#'   multiple-imputation uncertainty summary for pooled EFAs. `NULL` shows it for
#'   pooled objects.
#' @param ... Further arguments passed to [print.LOADINGS()].
#'
#' @returns `print()` and the print method for `summary.EFA` objects return their
#'   argument invisibly. `format()` returns a character vector with the report
#'   lines (styled to the active console theme; plain when colours are disabled).
#'   `summary()` returns an object of class `summary.EFA`.
#'
#' @export
#'
#' @method print EFA
#'
#' @examples
#' mod <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'            method = "PAF", rotation = "promax")
#' mod
#'
#' # The full diagnostics, CI tables, and residual diagnostics:
#' summary(mod)
#'
#' # format() returns plain text, e.g. for embedding in a report:
#' writeLines(format(mod))
#'
print.EFA <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.EFA
#' @export
#' @method print EFA_POOLED
print.EFA_POOLED <- function(x, ...) {
  print.EFA(x, ...)
}

#' @rdname print.EFA
#' @export
#' @method format EFA
format.EFA <- function(x, cutoff = .3, digits = 3, max_name_length = 10,
                       sort_loadings = c("none", "primary", "clustered"),
                       show_loading_legend = TRUE,
                       max_factors_per_block = NULL, ...) {
  sort_loadings <- match.arg(sort_loadings)

  .render_efa(x,
    view = "brief",
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    sort_loadings = sort_loadings,
    show_loading_legend = show_loading_legend,
    max_factors_per_block = max_factors_per_block,
    ...
  )
}

#' @rdname print.EFA
#' @export
#' @method format EFA_POOLED
format.EFA_POOLED <- function(x, ...) {
  format.EFA(x, ...)
}

#' @rdname print.EFA
#' @export
#' @method summary EFA
summary.EFA <- function(object, cutoff = .3, digits = 3, max_name_length = 10,
                        ci = c("auto", "none", "separate"),
                        ci_filter = c("salient", "all", "nonzero"),
                        diagnostics_top_n = 10,
                        residual_cutoff = .1,
                        residual_top_n = 10,
                        show_structure = TRUE,
                        sort_loadings = c("none", "primary", "clustered"),
                        show_loading_legend = TRUE,
                        cross_loading_cutoff = cutoff,
                        min_primary_gap = .20,
                        min_salient_per_factor = 3,
                        max_factors_per_block = NULL,
                        show_mi_diagnostics = NULL, ...) {
  ci <- match.arg(ci)
  ci_filter <- match.arg(ci_filter)
  sort_loadings <- match.arg(sort_loadings)

  opts <- list(
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    ci = ci,
    ci_filter = ci_filter,
    diagnostics_top_n = diagnostics_top_n,
    residual_cutoff = residual_cutoff,
    residual_top_n = residual_top_n,
    show_structure = show_structure,
    sort_loadings = sort_loadings,
    show_loading_legend = show_loading_legend,
    cross_loading_cutoff = cross_loading_cutoff,
    min_primary_gap = min_primary_gap,
    min_salient_per_factor = min_salient_per_factor,
    max_factors_per_block = max_factors_per_block,
    show_mi_diagnostics = show_mi_diagnostics
  )

  structure(list(efa = object, opts = opts), class = "summary.EFA")
}

#' @rdname print.EFA
#' @export
#' @method summary EFA_POOLED
summary.EFA_POOLED <- function(object, ...) {
  summary.EFA(object, ...)
}

#' @rdname print.EFA
#' @export
#' @method print summary.EFA
print.summary.EFA <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.EFA
#' @export
#' @method format summary.EFA
format.summary.EFA <- function(x, ...) {
  opts <- utils::modifyList(x$opts, list(...))
  do.call(.render_efa, c(list(x$efa, view = "full"), opts))
}

# Render the EFA report at the requested depth: "brief" (print) shows the header,
# loadings (with factor intercorrelations), variances accounted for, and model fit;
# "full" (summary) additionally shows model and simple-structure diagnostics, the
# CI tables, the structure matrix, MI uncertainty, and residual diagnostics.
.render_efa <- function(x, view = c("brief", "full"),
                        cutoff = .3, digits = 3, max_name_length = 10,
                        ci = c("auto", "none", "separate"),
                        ci_filter = c("salient", "all", "nonzero"),
                        diagnostics_top_n = 10,
                        residual_cutoff = .1,
                        residual_top_n = 10,
                        show_structure = TRUE,
                        sort_loadings = c("none", "primary", "clustered"),
                        show_loading_legend = TRUE,
                        cross_loading_cutoff = cutoff,
                        min_primary_gap = .20,
                        min_salient_per_factor = 3,
                        max_factors_per_block = NULL,
                        show_mi_diagnostics = NULL, ...) {
  view <- match.arg(view)
  ci <- match.arg(ci)
  ci_filter <- match.arg(ci_filter)
  sort_loadings <- match.arg(sort_loadings)
  full <- identical(view, "full")

  .efa_validate_print_options(
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    diagnostics_top_n = diagnostics_top_n,
    residual_cutoff = residual_cutoff,
    residual_top_n = residual_top_n,
    show_structure = show_structure,
    show_loading_legend = show_loading_legend,
    cross_loading_cutoff = cross_loading_cutoff,
    min_primary_gap = min_primary_gap,
    min_salient_per_factor = min_salient_per_factor,
    max_factors_per_block = max_factors_per_block,
    show_mi_diagnostics = show_mi_diagnostics
  )

  spec <- .efa_print_spec(x)
  show_mi <- if (!is.null(show_mi_diagnostics)) {
    isTRUE(show_mi_diagnostics)
  } else {
    isTRUE(spec$is_pooled)
  }

  # Assemble the whole report inside a single cli container so prose wraps to the console
  # width; `format()` returns these lines and `print()` is `cat(format(x), sep = "\n")`.
  cli::cli_format_method({
    .print_efa_header(spec)

    if (full) {
      .print_efa_diagnostics_section(x, spec,
        cutoff = cutoff,
        digits = digits,
        cross_loading_cutoff = cross_loading_cutoff,
        min_primary_gap = min_primary_gap,
        min_salient_per_factor = min_salient_per_factor
      )
    }

    .print_efa_loadings_section(x, spec,
      cutoff = cutoff,
      digits = digits,
      max_name_length = max_name_length,
      ci = if (full) ci else "none",
      ci_filter = ci_filter,
      show_structure = isTRUE(show_structure) && full,
      sort_loadings = sort_loadings,
      show_loading_legend = show_loading_legend,
      max_factors_per_block = max_factors_per_block,
      ...
    )

    if (full) {
      .print_efa_simple_structure_section(x, spec,
        cutoff = cutoff,
        digits = digits,
        max_name_length = max_name_length,
        cross_loading_cutoff = cross_loading_cutoff,
        min_primary_gap = min_primary_gap,
        min_salient_per_factor = min_salient_per_factor,
        diagnostics_top_n = diagnostics_top_n,
        sort_loadings = sort_loadings
      )
    }

    .print_efa_variances_section(x, spec, digits = digits)

    # A negative df (underidentified) stops the report after the warning; df == 0 (just
    # identified) emits its caveat but the fit/bootstrap/residual sections still follow.
    if (.print_efa_identification_warning(spec)) {
      .print_efa_model_fit_section(x, spec)
      .print_efa_bootstrap_note(x, spec, ci = if (full) ci else "none")

      if (full && isTRUE(show_mi)) {
        .print_efa_mi_diagnostics_section(x, spec, digits = digits)
      }

      if (full) {
        .print_efa_residuals_section(x,
          residual_cutoff = residual_cutoff,
          residual_top_n = residual_top_n,
          digits = digits
        )
      }
    }
  })
}

.efa_print_spec <- function(x) {
  settings <- x$settings
  if (is.null(settings)) {
    settings <- list()
  }

  se <- settings$se
  is_pooled <- inherits(x, "EFA_POOLED") || isTRUE(settings$pooled)

  list(
    method = settings$method,
    rotation = settings$rotation,
    rotation_diagnostics = settings$rotation_diagnostics,
    type = settings$type,
    cor_method = settings$cor_method,
    se = se,
    N = settings$N,
    fit = x$fit_indices,
    h2 = x$h2,
    heywood = x$heywood,
    np_boot = identical(.efa_setting_text(se), "np-boot"),
    ci = .efa_ci_level(x),
    rmsea_ci_level = .efa_rmsea_ci_level(x),
    b_boot = settings$b_boot,
    max_iter = settings$max_iter,
    iter = x$iter,
    convergence = x$convergence,
    is_pooled = is_pooled,
    n_imputations = .efa_n_imputations(x),
    target_method = settings$target_method,
    align_unrotated = settings$align_unrotated,
    fit_pool_method = settings$fit_pool_method,
    component_se = settings$component_se,
    has_ci = !is.null(x$CI)
  )
}

.print_efa_header <- function(spec) {
  cli::cli_text("")

  if (isTRUE(spec$is_pooled)) {
    .print_efa_pooled_header(spec)
  } else {
    # Emitted verbatim (not via cli_text) so the "setting = 'value'" tokens are never split
    # across a line break; the bold values therefore use style_bold() rather than {.strong}.
    cli::cli_verbatim(paste0(
      "EFA performed with type = '", cli::style_bold(spec$type),
      "', method = '", cli::style_bold(spec$method),
      "', and rotation = '", cli::style_bold(spec$rotation), "'."
    ))
  }

  # Flag two-stage FIML correlations (saturated mean/covariance EM-estimated from raw data
  # with missing values; Yuan, Marshall, & Bentler, 2002) so the analysed matrix is not read
  # as an ordinary complete-case correlation. Emitted verbatim (plain, no ANSI).
  if (identical(.efa_setting_text(spec$cor_method), "fiml")) {
    cli::cli_verbatim("Correlations: FIML (two-stage, missing data)")
  }

  if (.efa_iteration_nonconvergence(spec)) {
    cli::cli_text("")
    cli::cli_alert_danger(.efa_nonconvergence_banner(spec$method))
  }

  invisible(NULL)
}


.print_efa_pooled_header <- function(spec) {
  n_imputations <- .efa_setting_text(spec$n_imputations)

  across <- if (nzchar(n_imputations)) {
    paste0(" across ", cli::style_bold(n_imputations), " imputations")
  } else {
    ""
  }

  # Verbatim, for the same reason as the unpooled header above.
  cli::cli_verbatim(paste0(
    "Pooled EFA", across, " performed with type = '", cli::style_bold(spec$type),
    "', method = '", cli::style_bold(spec$method),
    "', and rotation = '", cli::style_bold(spec$rotation), "'."
  ))

  pooling_settings <- .efa_pooled_settings_text(spec)
  if (length(pooling_settings) > 0L) {
    cli::cli_verbatim(paste0("Pooling settings: ",
                             paste(pooling_settings, collapse = ", "), "."))
  }

  invisible(NULL)
}

.efa_pooled_settings_text <- function(spec) {
  out <- character(0)

  # The MI2S sandwich route fits a single model to the pooled inputs, so the
  # per-imputation alignment and chi-square pooling settings (target_method,
  # align_unrotated, fit_pool_method) do not apply and are not reported.
  if (!identical(spec$component_se, "sandwich")) {
    target_method <- .efa_setting_text(spec$target_method)
    if (nzchar(target_method) && !identical(spec$rotation, "none")) {
      out <- c(out, paste0("target_method = '", target_method, "'"))
    }

    align_unrotated <- .efa_setting_text(spec$align_unrotated)
    if (nzchar(align_unrotated)) {
      out <- c(out, paste0("align_unrotated = '", align_unrotated, "'"))
    }

    fit_pool_method <- .efa_setting_text(spec$fit_pool_method)
    if (nzchar(fit_pool_method)) {
      out <- c(out, paste0("fit_pool_method = '", fit_pool_method, "'"))
    }
  }

  if (is.finite(spec$rmsea_ci_level) &&
      abs(spec$rmsea_ci_level - .90) > sqrt(.Machine$double.eps)) {
    out <- c(out, paste0("rmsea_ci_level = ", .efa_ci_level_text(spec$rmsea_ci_level)))
  }

  out
}

.print_efa_loadings_section <- function(x, spec, cutoff, digits, max_name_length,
                                        ci, ci_filter, show_structure,
                                        sort_loadings, show_loading_legend,
                                        max_factors_per_block, ...) {
  if (identical(spec$rotation, "none")) {
    .print_efa_rule("Unrotated Loadings")

    if (is.null(x$unrot_loadings)) {
      cli::cli_text("No unrotated loading matrix available.")
      return(invisible(NULL))
    }

    .efa_emit_lines(.efa_capture_loadings(x$unrot_loadings,
      cutoff = cutoff,
      digits = digits,
      max_name_length = max_name_length,
      h2 = x$h2,
      sort_loadings = sort_loadings,
      legend = show_loading_legend,
      max_factors_per_block = max_factors_per_block,
      ...
    ))
    .print_efa_heywood_warning(spec$heywood, spec$h2)
    .print_efa_loading_ci_section(x, spec,
      component = "unrot_loadings",
      loadings = x$unrot_loadings,
      loading_label = "unrotated",
      cutoff = cutoff,
      digits = digits,
      max_name_length = max_name_length,
      ci = ci,
      ci_filter = ci_filter,
      sort_loadings = sort_loadings
    )
    return(invisible(NULL))
  }

  .print_efa_rule("Rotated Loadings")

  if (is.null(x$rot_loadings)) {
    cli::cli_text("No rotated loading matrix available.")
    return(invisible(NULL))
  }

  .efa_emit_lines(.efa_capture_loadings(x$rot_loadings,
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    h2 = x$h2,
    sort_loadings = sort_loadings,
    legend = show_loading_legend,
    max_factors_per_block = max_factors_per_block,
    ...
  ))
  .print_efa_heywood_warning(spec$heywood, spec$h2)
  .print_efa_loading_ci_section(x, spec,
    component = "rot_loadings",
    loadings = x$rot_loadings,
    loading_label = "rotated",
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    ci = ci,
    ci_filter = ci_filter,
    sort_loadings = sort_loadings
  )

  if (!is.null(x$Phi)) {
    phi <- as.matrix(x$Phi)
    factor_names <- .efa_factor_names(phi)
    dimnames(phi) <- list(factor_names, factor_names)

    .print_efa_rule("Factor Intercorrelations")
    .efa_emit_lines(.efa_corr_lines(phi, digits = digits, lower_only = TRUE))
    .print_efa_phi_ci_section(x, spec, digits, ci)
  }

  if (!is.null(x$Structure) && isTRUE(show_structure)) {
    .print_efa_structure_section(x,
      cutoff = cutoff,
      digits = digits,
      max_name_length = max_name_length,
      sort_loadings = sort_loadings,
      max_factors_per_block = max_factors_per_block,
      ...
    )
  }

  invisible(NULL)
}

.print_efa_variances_section <- function(x, spec, digits = 3) {
  .print_efa_rule("Variances Accounted for")

  if (identical(spec$rotation, "none")) {
    vars_accounted <- x$vars_accounted
  } else {
    vars_accounted <- x$vars_accounted_rot
  }

  if (is.null(vars_accounted)) {
    cli::cli_text("No variance-accounted table available.")
    return(invisible(NULL))
  }

  .efa_emit_lines(.efa_corr_lines(vars_accounted, digits = digits))

  invisible(NULL)
}

# Format a "corr"-role matrix (factor intercorrelations or variances accounted for) into
# decimal-aligned lines for embedding in the cli report via `cli::cli_verbatim()`. Reuses
# the shared `.efa_format_matrix()` renderer, so these tables align consistently with the
# loading table; rows/columns are labelled from the matrix dimnames and `lower_only = TRUE`
# blanks the strictly-upper triangle (for symmetric matrices such as Phi).
.efa_corr_lines <- function(values, digits = 3, lower_only = FALSE) {
  values <- as.matrix(values)
  .efa_format_matrix(
    values = values,
    row_labels = .efa_variable_names(values),
    col_labels = .efa_factor_names(values),
    col_roles = rep("corr", ncol(values)),
    digits = digits,
    lower_only = lower_only
  )
}

# Capture an already-styled loading/structure table as verbatim lines for embedding in the
# cli report. Reuses all of `print.LOADINGS`' formatting (h2/sorting/legend/column blocks)
# and preserves colour; the trailing blank line `print.LOADINGS` adds for console spacing is
# dropped so the surrounding cli layout controls the spacing.
.efa_capture_loadings <- function(x, ...) {
  out <- utils::capture.output(print(x, ...))
  if (length(out) > 0L && !nzchar(out[length(out)])) out <- out[-length(out)]
  out
}

# Emit pre-formatted table lines into the cli report. `cli::cli_verbatim()` silently drops
# interior empty-string elements, which would erase the blank line a loading table puts
# before its legend and the blank separators between column blocks of a wide table. So emit
# each run of non-empty lines verbatim (no reflow) and each blank line as `cli_text("")`,
# preserving the table's own vertical spacing.
.efa_emit_lines <- function(lines) {
  if (length(lines) < 1L) {
    return(invisible(NULL))
  }

  blank <- !nzchar(lines)
  run <- cumsum(c(TRUE, blank[-1L] != blank[-length(blank)]))
  for (r in unique(run)) {
    idx <- which(run == r)
    if (blank[idx[1L]]) {
      for (i in idx) cli::cli_text("")
    } else {
      cli::cli_verbatim(lines[idx])
    }
  }

  invisible(NULL)
}

# Emit a bullet list whose items embed user-supplied variable/factor names. `cli::cli_ul()`
# treats each item as a cli/glue template, so a name containing `{` or `}` would be evaluated
# as an expression (corrupting output or aborting). Double the braces so they render
# literally.
.efa_emit_bullets <- function(items) {
  items <- gsub("}", "}}", gsub("{", "{{", items, fixed = TRUE), fixed = TRUE)
  cli::cli_ul(items)
  invisible(NULL)
}

.print_efa_identification_warning <- function(spec) {
  df <- .efa_fit_scalar(spec$fit, "df")

  if (!is.finite(df)) {
    return(TRUE)
  }

  if (df == 0) {
    cli::cli_text("")
    cli::cli_alert_warning(
      "The model is just identified (df = 0). Goodness of fit indices may not be interpretable."
    )
    return(TRUE)
  }

  if (df < 0) {
    cli::cli_text("")
    cli::cli_alert_warning(
      "The model is underidentified (df < 0). No goodness of fit indices were calculated."
    )
    return(FALSE)
  }

  TRUE
}

.print_efa_model_fit_section <- function(x, spec) {
  fit <- spec$fit
  fit_ci <- .efa_fit_ci_labels(x, spec)

  .print_efa_rule("Model Fit")

  if (is.null(fit) || length(fit) < 1L) {
    cli::cli_text("No model-fit indices available.")
    return(invisible(NULL))
  }

  # One "name[ CI-tag]: value[ CI]" data line. The CI tag is emitted only for
  # indices whose own CI suffix is available, so partially populated fit-index
  # CI lists do not label plain point estimates as interval estimates.
  fit_line <- function(name, key, ci = fit_ci[[key]], tag = TRUE, ...) {
    if (is.null(ci)) ci <- ""
    label <- if (isTRUE(tag) && nzchar(ci)) fit_ci$label else ""
    paste0(name, label, ": ", .efa_format_fit_value(fit, key, pad = FALSE, ...), ci)
  }
  is_finite_fit <- function(key) is.finite(.efa_fit_scalar(fit, key))

  method <- .efa_setting_text(spec$method)
  df <- .efa_fit_scalar(fit, "df")
  chi <- .efa_fit_scalar(fit, "chi")
  p_chi <- .efa_fit_scalar(fit, "p_chi")

  if (identical(method, "PAF") || .efa_is_missing_number(spec$N) ||
      !is.finite(chi) || !is.finite(df)) {
    lines <- fit_line("CAF", "CAF")
    if (is_finite_fit("SRMR")) {
      lines <- c(lines, fit_line("SRMR", "SRMR"))
    }
    lines <- c(lines, paste0("df: ",
      .efa_format_fit_value(fit, "df", digits = 0, print_zero = TRUE, pad = FALSE)))
    .efa_emit_lines(lines)
    return(invisible(NULL))
  }

  p_text <- if (!is.finite(p_chi)) {
    " = NA"
  } else if (p_chi < .001) {
    " < .001"
  } else {
    paste0(" = ", .efa_format_number(p_chi, digits = 3, pad = FALSE))
  }

  df_text <- .efa_format_fit_value(fit, "df", digits = 0, print_zero = TRUE, pad = FALSE)
  chi_text <- .efa_format_fit_value(fit, "chi", digits = 2, print_zero = TRUE, pad = FALSE)

  rmsea_label <- paste0("RMSEA [", .efa_ci_level_text(spec$rmsea_ci_level), " CI]")
  rmsea_value <- paste0(
    .efa_format_fit_value(fit, "RMSEA", pad = FALSE), " [",
    .efa_format_fit_value(fit, "RMSEA_LB", pad = FALSE), "; ",
    .efa_format_fit_value(fit, "RMSEA_UB", pad = FALSE), "]"
  )

  # When the sandwich path supplies a scaled (Satorra-Bentler) chi-square, mark the line so
  # it is not read as an ordinary chi-square (the CFI/TLI/RMSEA below are derived from it).
  # A D2-pooled fit is flagged the same way: its chi-square is the Li et al. (1991)
  # asymptotic statistic df * F_D2, whose p-value comes from the D2 reference F(df1, df2)
  # rather than the chi-square(df) tail (the note below spells this out when they diverge).
  is_d2 <- isTRUE(spec$is_pooled) && identical(fit[["pool_method"]], "D2")
  chi_prefix <- if (!is.null(fit[["chi_scaled_type"]]) && !is.na(fit[["chi_scaled_type"]])) {
    "scaled "
  } else if (is_d2) {
    "D2-pooled "
  } else {
    ""
  }

  # The chi-square line carries an italic "p"; emit verbatim (with style_italic) so the
  # styling survives and the line is not reflowed.
  lines <- paste0(chi_prefix, "\u03c7\u00b2(", df_text, ") = ", chi_text, ", ",
                  cli::style_italic("p"), p_text)
  lines <- c(lines, fit_line("CFI", "CFI"))
  if (is_finite_fit("TLI")) {
    lines <- c(lines, fit_line("TLI", "TLI"))
  }
  lines <- c(lines, paste0(rmsea_label, fit_ci$label, ": ", rmsea_value, fit_ci$RMSEA))
  lines <- c(lines,
             fit_line("AIC", "AIC", fit_ci$AIC, print_zero = TRUE),
             fit_line("BIC", "BIC", fit_ci$BIC, print_zero = TRUE))
  if (is_finite_fit("ECVI")) {
    lines <- c(lines, fit_line("ECVI", "ECVI", print_zero = TRUE))
  }
  lines <- c(lines, fit_line("CAF", "CAF"))
  if (is_finite_fit("SRMR")) {
    lines <- c(lines, fit_line("SRMR", "SRMR"))
  }

  .efa_emit_lines(lines)

  # The D2-pooled chi-square (Li, Meng, Raghunathan & Rubin, 1991) is the asymptotic
  # statistic df * F_D2 and its p-value comes from the reference F(df1, df2); state that
  # so the reported p is not mistaken for the chi-square(df) tail probability. The two
  # references coincide only as df2 -> Inf (no between-imputation variance), so the note is
  # shown when df2 is finite.
  if (is_d2) {
    d2_df2 <- if (is.null(fit$mi_diagnostics)) NA_real_ else fit$mi_diagnostics$D2_df2
    if (is.finite(d2_df2)) {
      d2_df2_text <- .efa_format_plain_number(d2_df2, digits = 1)
      cli::cli_text(
        "{.emph Note:} the pooled \u03c7\u00b2 is the D2 statistic; its {.emph p} uses the D2 reference F({df_text}, {d2_df2_text}), not the \u03c7\u00b2({df_text}) tail."
      )
    }
  }

  invisible(NULL)
}

.efa_fit_ci_labels <- function(x, spec) {
  out <- list(
    label = "",
    CAF = "",
    RMSR = "",
    SRMR = "",
    CFI = "",
    TLI = "",
    RMSEA = "",
    AIC = "",
    BIC = "",
    ECVI = ""
  )

  if (!isTRUE(spec$has_ci)) {
    return(out)
  }

  fitind_ci <- .efa_get_fit_ci(x, spec)
  if (is.null(fitind_ci)) {
    return(out)
  }

  out$label <- paste0(" [", .efa_ci_header_label(spec), "-CI]")
  out$CAF <- .efa_paste_gof_ci(fitind_ci, "CAF")
  out$RMSR <- .efa_paste_gof_ci(fitind_ci, "RMSR")
  out$SRMR <- .efa_paste_gof_ci(fitind_ci, "SRMR")
  out$CFI <- .efa_paste_gof_ci(fitind_ci, "CFI")
  out$TLI <- .efa_paste_gof_ci(fitind_ci, "TLI")
  out$RMSEA <- .efa_paste_gof_ci(fitind_ci, "RMSEA")
  out$AIC <- .efa_paste_gof_ci(fitind_ci, "AIC")
  out$BIC <- .efa_paste_gof_ci(fitind_ci, "BIC")
  out$ECVI <- .efa_paste_gof_ci(fitind_ci, "ECVI")

  out
}

.print_efa_bootstrap_note <- function(x, spec, ci = NULL) {
  if (!isTRUE(spec$has_ci)) {
    return(invisible(NULL))
  }

  # Analytic SEs produce Wald intervals, not resampled ones; describe their source rather
  # than a bootstrap-sample count. For a NON-pooled analytic fit the Wald CIs are attached
  # only to the unrotated loadings (displayed only for an unrotated solution with CIs
  # enabled), so suppress the note when no such table is shown. Pooled objects, however,
  # attach and display rotated-loading CIs on both the information and MI2S routes, so the
  # note must still appear for a rotated pooled solution. The wording follows the SE method:
  # se = "information" draws Wald CIs from the expected information matrix (Rubin-pooled with
  # plain Rubin (1987) df when pooled), while se = "sandwich" draws them from the robust (Godambe)
  # sandwich covariance (built on the multiple-imputation-pooled asymptotic covariance when
  # pooled), so the pooled MI2S note must not claim information-matrix Rubin pooling.
  if (!isTRUE(spec$np_boot)) {
    if (identical(ci, "none") ||
        (!identical(spec$rotation, "none") && !isTRUE(spec$is_pooled))) {
      return(invisible(NULL))
    }
    cli::cli_text("")
    if (isTRUE(spec$is_pooled) && identical(spec$component_se, "sandwich")) {
      cli::cli_text("{.emph Note:} Wald CIs from the robust sandwich covariance built on the multiple-imputation-pooled asymptotic covariance.")
    } else if (isTRUE(spec$is_pooled)) {
      cli::cli_text("{.emph Note:} Wald CIs from the expected information matrix, Rubin-pooled across imputations (Rubin's 1987 degrees of freedom).")
    } else if (identical(.efa_setting_text(spec$cor_method), "fiml")) {
      # Under cor_method = "fiml" both analytic settings yield the corrected two-stage sandwich, so
      # the note describes that covariance rather than the expected information matrix.
      cli::cli_text("{.emph Note:} Wald CIs from the corrected two-stage (FIML) sandwich covariance.")
    } else if (identical(spec$se, "sandwich")) {
      cli::cli_text("{.emph Note:} Wald CIs from the robust (Godambe) sandwich covariance.")
    } else {
      cli::cli_text("{.emph Note:} Wald CIs from the expected information matrix.")
    }
    return(invisible(NULL))
  }

  note_label <- if (isTRUE(spec$is_pooled)) {
    "Bootstrap/MI CIs"
  } else {
    "Bootstrap CIs"
  }

  note <- note_label

  b_text <- .efa_bootstrap_sample_text(spec)
  if (nzchar(b_text)) {
    note <- paste0(note, " based on ", b_text)
  }

  target_rotation_note <- .efa_target_rotation_note(x, spec)
  if (nzchar(target_rotation_note)) {
    note <- paste0(note, "; ", target_rotation_note)
  }

  note <- paste0(note, ".")

  cli::cli_text("")
  cli::cli_text("{.emph Note:} {note}")

  invisible(NULL)
}

.print_efa_residuals_section <- function(x, residual_cutoff = .1,
                                     residual_top_n = 10, digits = 3) {
  .print_efa_rule("Residual Diagnostics")

  if (is.null(x$residuals)) {
    cli::cli_text("No residual matrix available.")
    return(invisible(NULL))
  }

  residuals <- as.matrix(x$residuals)
  idx <- which(upper.tri(residuals) & is.finite(residuals), arr.ind = TRUE)

  if (nrow(idx) < 1L) {
    cli::cli_text("No finite off-diagonal residuals available.")
    return(invisible(NULL))
  }

  values <- residuals[idx]
  keep <- abs(values) > residual_cutoff
  n_large <- sum(keep, na.rm = TRUE)
  largest_abs <- max(abs(values), na.rm = TRUE)

  cutoff_text <- .efa_format_plain_number(residual_cutoff, digits)
  largest_text <- .efa_format_plain_number(largest_abs, digits)
  cli::cli_text("Residual cutoff: |r| > {cutoff_text}")
  cli::cli_text("Number of large residuals: {n_large}")
  cli::cli_text("Largest absolute residual: {largest_text}")

  if (n_large < 1L) {
    cli::cli_text("")
    cli::cli_text("No absolute residuals > {cutoff_text} occurred.")
    cli::cli_text("")
    cli::cli_text("Inspect the residual matrix for details (e.g., with residuals()).")
    return(invisible(NULL))
  }

  large_idx <- idx[keep, , drop = FALSE]
  large_values <- values[keep]
  ord <- order(abs(large_values), decreasing = TRUE)
  large_idx <- large_idx[ord, , drop = FALSE]
  large_values <- large_values[ord]

  if (is.finite(residual_top_n)) {
    n_print <- min(length(large_values), as.integer(residual_top_n))
  } else {
    n_print <- length(large_values)
  }

  large_idx <- large_idx[seq_len(n_print), , drop = FALSE]
  large_values <- large_values[seq_len(n_print)]

  # The residual matrix is symmetric, so prefer column names but fall back to row
  # names (unlike `.efa_variable_names()`, which only consults row names).
  var_names <- colnames(residuals)
  if (is.null(var_names)) {
    var_names <- rownames(residuals)
  }
  if (is.null(var_names)) {
    var_names <- paste0("V", seq_len(ncol(residuals)))
  }

  heading <- if (n_print < n_large) {
    paste0("Largest residuals (top ", n_print, " of ", n_large, "):")
  } else {
    "Largest residuals:"
  }

  bullets <- vapply(seq_along(large_values), function(i) {
    paste0(var_names[large_idx[i, 1]], " ~~ ", var_names[large_idx[i, 2]],
           ": ", .efa_format_plain_number(large_values[i], digits))
  }, character(1L))

  cli::cli_text("")
  cli::cli_text("{heading}")
  .efa_emit_bullets(bullets)
  cli::cli_text("")
  cli::cli_text("Inspect the residual matrix for details (e.g., with residuals()).")

  invisible(NULL)
}

.print_efa_loading_ci_section <- function(x, spec, component, loadings, loading_label,
                                          cutoff, digits, max_name_length,
                                          ci, ci_filter,
                                          sort_loadings = c("none", "primary", "clustered")) {
  sort_loadings <- match.arg(sort_loadings)

  if (!.efa_should_print_ci(x, ci)) {
    return(invisible(NULL))
  }

  loading_ci <- .efa_get_component_ci(x, component)
  if (is.null(loading_ci)) {
    return(invisible(NULL))
  }

  loadings <- as.matrix(loadings)
  lower <- as.matrix(loading_ci$lower)
  upper <- as.matrix(loading_ci$upper)

  if (!all(dim(loadings) == dim(lower)) || !all(dim(loadings) == dim(upper))) {
    return(invisible(NULL))
  }

  row_order <- .efa_loading_row_order(loadings, sort_loadings)
  loadings <- loadings[row_order, , drop = FALSE]
  lower <- lower[row_order, , drop = FALSE]
  upper <- upper[row_order, , drop = FALSE]

  keep <- .efa_select_loading_cis(loadings, lower, upper, cutoff, ci_filter)
  title <- .efa_loading_ci_title(spec, loading_label, ci_filter)
  .print_efa_rule(title)

  if (!any(keep, na.rm = TRUE)) {
    cli::cli_text("No loading CIs matched ci_filter = '{ci_filter}'.")
    return(invisible(NULL))
  }

  idx <- which(keep, arr.ind = TRUE)
  row_names <- .efa_truncate_names(.efa_variable_names(loadings), max_name_length)
  factor_names <- .efa_factor_names(loadings)

  .print_efa_ci_table(
    key = row_names[idx[, 1]],
    second_key = factor_names[idx[, 2]],
    estimate = loadings[idx],
    lower = lower[idx],
    upper = upper[idx],
    digits = digits,
    key_header = "Variable",
    second_key_header = "Factor",
    highlight = abs(loadings[idx]) >= cutoff
  )

  invisible(NULL)
}

.print_efa_phi_ci_section <- function(x, spec, digits, ci) {
  if (!.efa_should_print_ci(x, ci)) {
    return(invisible(NULL))
  }

  phi_ci <- .efa_get_component_ci(x, "Phi")
  if (is.null(phi_ci) || is.null(x$Phi)) {
    return(invisible(NULL))
  }

  phi <- as.matrix(x$Phi)
  lower <- as.matrix(phi_ci$lower)
  upper <- as.matrix(phi_ci$upper)

  if (!all(dim(phi) == dim(lower)) || !all(dim(phi) == dim(upper))) {
    return(invisible(NULL))
  }

  keep <- lower.tri(phi) & is.finite(phi) & is.finite(lower) & is.finite(upper)
  if (!any(keep, na.rm = TRUE)) {
    return(invisible(NULL))
  }

  idx <- which(keep, arr.ind = TRUE)
  factor_names <- .efa_factor_names(phi)

  .print_efa_rule(.efa_phi_ci_title(spec))
  .print_efa_ci_table(
    key = paste0(factor_names[idx[, 2]], " ~~ ", factor_names[idx[, 1]]),
    second_key = NULL,
    estimate = phi[idx],
    lower = lower[idx],
    upper = upper[idx],
    digits = digits,
    key_header = "Factors",
    second_key_header = NULL,
    highlight = NULL
  )

  invisible(NULL)
}

# Count Heywood cases: prefer the detector's `heywood` field (covers communalities
# >= 1 and ML/ULS boundary uniquenesses); fall back to communalities >= 1 for objects
# that carry no such field (e.g. pooled solutions).
.efa_heywood_count <- function(heywood, h2 = NULL) {
  if (!is.null(heywood)) {
    length(heywood)
  } else if (!is.null(h2)) {
    sum(h2 >= 1 + .Machine$double.eps, na.rm = TRUE)
  } else {
    0L
  }
}

.print_efa_heywood_warning <- function(heywood, h2 = NULL) {
  n_heywood <- .efa_heywood_count(heywood, h2)

  if (n_heywood >= 1) {
    cli::cli_text("")
    cli::cli_alert_danger("{n_heywood} Heywood case{?s} {?was/were} detected!")
  }

  invisible(NULL)
}

.print_efa_rule <- function(title) {
  cli::cli_text("")
  cli::cli_rule(left = "{.strong {title}}")
  cli::cli_text("")

  invisible(NULL)
}

.efa_should_print_ci <- function(x, ci) {
  !identical(ci, "none") && !is.null(x$CI)
}

.efa_get_component_ci <- function(x, component) {
  if (is.null(x$CI)) {
    return(NULL)
  }

  ci <- x$CI[[component]]
  if (.efa_is_ci_pair(ci)) {
    return(ci)
  }

  NULL
}

.efa_get_fit_ci <- function(x, spec = .efa_print_spec(x)) {
  if (is.null(x$CI)) {
    return(NULL)
  }

  components <- if (isTRUE(spec$is_pooled)) {
    c("fit_indices_descriptive", "fit_indices")
  } else {
    c("fit_indices", "fit_indices_descriptive")
  }

  for (component in components) {
    ci <- x$CI[[component]]
    if (.efa_is_ci_pair(ci)) {
      return(ci)
    }
  }

  NULL
}

.efa_is_ci_pair <- function(x) {
  is.list(x) && !is.null(x$lower) && !is.null(x$upper)
}

.efa_select_loading_cis <- function(loadings, lower, upper, cutoff, ci_filter) {
  finite <- is.finite(loadings) & is.finite(lower) & is.finite(upper)

  if (identical(ci_filter, "all")) {
    return(finite)
  }

  if (identical(ci_filter, "nonzero")) {
    return(finite & (lower > 0 | upper < 0))
  }

  finite & abs(loadings) >= cutoff
}

.efa_n_imputations <- function(x) {
  n_imputations <- x$settings$n_imputations

  if ((is.null(n_imputations) || length(n_imputations) < 1L || is.na(n_imputations[1])) &&
      !is.null(x$fits)) {
    n_imputations <- length(x$fits)
  }

  if (is.null(n_imputations) || length(n_imputations) < 1L || is.na(n_imputations[1])) {
    return(NA_integer_)
  }

  as.integer(n_imputations[1])
}

.efa_setting_text <- function(x) {
  if (is.null(x) || length(x) < 1L || is.na(x[1])) {
    return("")
  }

  as.character(x[1])
}

.efa_ci_level <- function(x) {
  val <- x$settings$ci
  if (is.null(val) || length(val) < 1L || is.na(val[1])) {
    return(NA_real_)
  }

  out <- suppressWarnings(as.numeric(val[1]))
  if (length(out) < 1L || is.na(out)) {
    return(NA_real_)
  }

  out
}

.efa_rmsea_ci_level <- function(x) {
  val <- x$settings$rmsea_ci_level
  if (is.null(val) || length(val) < 1L || is.na(val[1])) {
    return(.90)
  }

  out <- suppressWarnings(as.numeric(val[1]))
  if (length(out) < 1L || !is.finite(out) || out <= 0 || out >= 1) {
    return(.90)
  }

  out
}

.efa_ci_level_text <- function(level) {
  if (!is.finite(level)) {
    return("")
  }

  paste0(formatC(level * 100, format = "fg", digits = 4, width = -1), "%")
}

.efa_ci_source_text <- function(spec) {
  if (isTRUE(spec$is_pooled)) {
    # Label the pooled CI source by route. The information route Rubin-pools the
    # per-imputation Wald intervals across imputations ("Wald/MI"). The sandwich
    # (MI2S) route fits a single model to the pooled inputs and carries the
    # imputation uncertainty in the pooled asymptotic covariance rather than by
    # per-parameter Rubin pooling, so its intervals are Wald intervals from that
    # single sandwich fit ("sandwich/MI") -- not Rubin-pooled. The bootstrap
    # route's within-imputation variance is bootstrap-based ("bootstrap/MI").
    if (identical(spec$component_se, "information")) {
      return("Wald/MI")
    }
    if (identical(spec$component_se, "sandwich")) {
      return("sandwich/MI")
    }
    return("bootstrap/MI")
  }

  # Analytic intervals are Wald (estimate +/- z * SE); only the bootstrap CIs are percentile.
  if (!isTRUE(spec$np_boot)) {
    return("Wald")
  }

  "bootstrap"
}

.efa_target_rotation_note <- function(x, spec) {
  counts <- .efa_valid_target_rotations(x)
  if (is.null(counts) || length(counts) < 1L || all(is.na(counts))) {
    return("")
  }

  counts <- as.integer(stats::na.omit(counts))
  if (length(counts) < 1L) {
    return("")
  }

  total <- .efa_total_bootstrap_samples(x, spec)
  if (!is.na(total) && !any(counts < total)) {
    return("")
  }

  component_text <- .efa_target_rotation_component_text(x)
  count_text <- .efa_valid_target_rotation_count_text(counts, spec)
  total_text <- if (!is.na(total)) {
    paste0(" out of ", total)
  } else {
    ""
  }
  scope_text <- if (isTRUE(spec$is_pooled)) {
    " per imputation"
  } else {
    ""
  }

  paste0(
    component_text,
    " based on ", count_text,
    " valid target-rotated bootstrap samples",
    total_text,
    scope_text
  )
}

.efa_valid_target_rotations <- function(x) {
  if (!is.null(x$SE$valid_target_rotations)) {
    return(x$SE$valid_target_rotations)
  }

  if (!is.null(x$MI$bootstrap_rotation_valid)) {
    return(x$MI$bootstrap_rotation_valid)
  }

  NULL
}

.efa_total_bootstrap_samples <- function(x, spec) {
  if (!is.null(spec$b_boot) && length(spec$b_boot) >= 1L &&
      is.finite(spec$b_boot[1])) {
    return(as.integer(spec$b_boot[1]))
  }

  if (!is.null(x$replicates$unrot_loadings)) {
    arr <- x$replicates$unrot_loadings
    if (length(dim(arr)) >= 3L) {
      return(dim(arr)[3L])
    }
  }

  NA_integer_
}

.efa_target_rotation_component_text <- function(x) {
  if (!is.null(x$Phi)) {
    return("rotated-loading, factor-intercorrelation, and structure-coefficient CIs")
  }

  "rotated-loading CIs"
}

.efa_valid_target_rotation_count_text <- function(counts, spec) {
  counts <- as.integer(counts)

  if (!isTRUE(spec$is_pooled) || length(unique(counts)) == 1L) {
    return(as.character(counts[1L]))
  }

  paste0("between ", min(counts), " and ", max(counts))
}

.efa_ci_header_label <- function(spec) {
  level <- .efa_ci_level_text(spec$ci)
  source <- .efa_ci_source_text(spec)

  if (nzchar(level)) {
    paste(level, source)
  } else {
    source
  }
}

.efa_loading_ci_title <- function(spec, loading_label, ci_filter) {
  filter_text <- switch(ci_filter,
    salient = "salient",
    all = "all",
    nonzero = "nonzero",
    ci_filter
  )

  paste0(
    .efa_ci_header_label(spec),
    " CIs for ", filter_text, " ", loading_label, " loadings"
  )
}

.efa_phi_ci_title <- function(spec) {
  paste0(.efa_ci_header_label(spec), " CIs for factor intercorrelations")
}

.efa_truncate_names <- function(x, max_name_length) {
  x <- as.character(x)
  too_long <- nchar(x) > max_name_length
  x[too_long] <- substr(x[too_long], 1, max_name_length)
  x
}

# Format a vector of CI numbers for the CI tables. `.efa_num` is vectorized, so one
# call formats the whole vector; the result feeds data.frame columns that drop names.
.efa_format_ci_num <- function(x, digits) {
  .efa_num(round(x, digits = digits), digits = digits)
}

.print_efa_ci_table <- function(key, second_key = NULL, estimate, lower, upper,
                                digits, key_header, second_key_header = NULL,
                                highlight = NULL) {
  has_second_key <- !is.null(second_key)

  # Right-pad each cell to its column width (base equivalent of str_pad side = "right").
  pad <- function(s, w) formatC(s, width = w, flag = "-")

  estimate_chr <- .efa_format_ci_num(estimate, digits)
  lower_chr <- .efa_format_ci_num(lower, digits)
  upper_chr <- .efa_format_ci_num(upper, digits)

  if (has_second_key) {
    plain_rows <- data.frame(
      key = as.character(key),
      second_key = as.character(second_key),
      estimate = estimate_chr,
      lower = lower_chr,
      upper = upper_chr,
      stringsAsFactors = FALSE
    )
    headers <- c(key_header, second_key_header, "est", "lower", "upper")
  } else {
    plain_rows <- data.frame(
      key = as.character(key),
      estimate = estimate_chr,
      lower = lower_chr,
      upper = upper_chr,
      stringsAsFactors = FALSE
    )
    headers <- c(key_header, "est", "lower", "upper")
  }

  widths <- vapply(seq_along(headers), function(j) {
    max(nchar(c(headers[j], plain_rows[[j]])), na.rm = TRUE)
  }, numeric(1))

  header <- paste(
    mapply(pad, headers, widths, USE.NAMES = FALSE),
    collapse = "  "
  )

  rows <- vapply(seq_len(nrow(plain_rows)), function(i) {
    cells <- as.character(unlist(plain_rows[i, , drop = FALSE], use.names = FALSE))
    padded <- mapply(pad, cells, widths, USE.NAMES = FALSE)

    if (has_second_key && !is.null(highlight) && isTRUE(highlight[i])) {
      padded[3] <- cli::style_bold(padded[3])
    }

    paste(padded, collapse = "  ")
  }, character(1L))

  # A fixed-width table; emit verbatim so cli does not reflow the aligned columns.
  .efa_emit_lines(c(header, rows))

  invisible(NULL)
}


.efa_validate_print_options <- function(cutoff, digits, max_name_length,
                                        diagnostics_top_n,
                                        residual_cutoff, residual_top_n,
                                        show_structure, show_loading_legend,
                                        cross_loading_cutoff, min_primary_gap,
                                        min_salient_per_factor,
                                        max_factors_per_block,
                                        show_mi_diagnostics) {
  if (!is.numeric(cutoff) || length(cutoff) != 1L ||
      !is.finite(cutoff) || cutoff < 0) {
    cli::cli_abort("{.arg cutoff} must be a single finite non-negative number.",
                   class = "efa_print_invalid_cutoff")
  }

  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits < 0 || digits != as.integer(digits)) {
    cli::cli_abort("{.arg digits} must be a single finite non-negative integer.",
                   class = "efa_print_invalid_digits")
  }

  if (!is.numeric(max_name_length) || length(max_name_length) != 1L ||
      !is.finite(max_name_length) || max_name_length < 1 ||
      max_name_length != as.integer(max_name_length)) {
    cli::cli_abort("{.arg max_name_length} must be a single finite positive integer.",
                   class = "efa_print_invalid_max_name_length")
  }

  if (!.efa_is_top_n(diagnostics_top_n)) {
    cli::cli_abort("{.arg diagnostics_top_n} must be a positive integer or {.val Inf}.",
                   class = "efa_print_invalid_diagnostics_top_n")
  }

  if (!is.numeric(residual_cutoff) || length(residual_cutoff) != 1L ||
      !is.finite(residual_cutoff) || residual_cutoff < 0) {
    cli::cli_abort("{.arg residual_cutoff} must be a single finite non-negative number.",
                   class = "efa_print_invalid_residual_cutoff")
  }

  if (!.efa_is_top_n(residual_top_n)) {
    cli::cli_abort("{.arg residual_top_n} must be a positive integer or {.val Inf}.",
                   class = "efa_print_invalid_residual_top_n")
  }

  if (!is.logical(show_structure) || length(show_structure) != 1L ||
      is.na(show_structure)) {
    cli::cli_abort("{.arg show_structure} must be {.val TRUE} or {.val FALSE}.",
                   class = "efa_print_invalid_show_structure")
  }

  if (!is.logical(show_loading_legend) || length(show_loading_legend) != 1L ||
      is.na(show_loading_legend)) {
    cli::cli_abort("{.arg show_loading_legend} must be {.val TRUE} or {.val FALSE}.",
                   class = "efa_print_invalid_show_loading_legend")
  }

  if (!is.numeric(cross_loading_cutoff) || length(cross_loading_cutoff) != 1L ||
      !is.finite(cross_loading_cutoff) || cross_loading_cutoff < 0) {
    cli::cli_abort("{.arg cross_loading_cutoff} must be a single finite non-negative number.",
                   class = "efa_print_invalid_cross_loading_cutoff")
  }

  if (!is.numeric(min_primary_gap) || length(min_primary_gap) != 1L ||
      !is.finite(min_primary_gap) || min_primary_gap < 0) {
    cli::cli_abort("{.arg min_primary_gap} must be a single finite non-negative number.",
                   class = "efa_print_invalid_min_primary_gap")
  }

  if (!is.numeric(min_salient_per_factor) || length(min_salient_per_factor) != 1L ||
      !is.finite(min_salient_per_factor) || min_salient_per_factor < 1 ||
      min_salient_per_factor != as.integer(min_salient_per_factor)) {
    cli::cli_abort("{.arg min_salient_per_factor} must be a single finite positive integer.",
                   class = "efa_print_invalid_min_salient_per_factor")
  }

  if (!is.null(max_factors_per_block) &&
      (!is.numeric(max_factors_per_block) || length(max_factors_per_block) != 1L ||
       !is.finite(max_factors_per_block) || max_factors_per_block < 1 ||
       max_factors_per_block != as.integer(max_factors_per_block))) {
    cli::cli_abort("{.arg max_factors_per_block} must be {.val NULL} or a single finite positive integer.",
                   class = "efa_print_invalid_max_factors_per_block")
  }

  if (!is.null(show_mi_diagnostics) &&
      (!is.logical(show_mi_diagnostics) || length(show_mi_diagnostics) != 1L ||
       is.na(show_mi_diagnostics))) {
    cli::cli_abort("{.arg show_mi_diagnostics} must be {.val TRUE}, {.val FALSE}, or {.val NULL}.",
                   class = "efa_print_invalid_show_mi_diagnostics")
  }

  invisible(TRUE)
}

.efa_iteration_nonconvergence <- function(spec) {
  # All current fits report a convergence code (0 = converged); a non-zero code
  # means the estimator stopped before meeting its tolerance. The iteration-count
  # fallback below covers objects that carry no convergence code but do report iter
  # and max_iter (e.g. solutions saved by an older package version).
  convergence <- .efa_as_scalar_number(spec$convergence)
  if (is.finite(convergence) && convergence != 0) {
    return(TRUE)
  }

  iter <- .efa_as_scalar_number(spec$iter)
  max_iter <- .efa_as_scalar_number(spec$max_iter)

  is.finite(iter) && is.finite(max_iter) && iter >= max_iter
}

# Message for the non-convergence banner (shared by print.EFA() and print.SL()).
# PAF non-convergence always means the iteration cap was hit; the optimiser-based
# estimators report a convergence code that need not be an iteration limit.
.efa_nonconvergence_banner <- function(method) {
  if (identical(method, "PAF")) {
    "Maximum number of iterations reached without convergence"
  } else {
    "The optimiser did not converge"
  }
}

.efa_as_scalar_number <- function(x) {
  if (is.null(x) || length(x) < 1L) {
    return(NA_real_)
  }

  out <- tryCatch(
    suppressWarnings(as.numeric(x[1])),
    error = function(e) NA_real_
  )
  if (length(out) < 1L) {
    return(NA_real_)
  }

  out
}

.efa_is_missing_number <- function(x) {
  val <- .efa_as_scalar_number(x)
  !is.finite(val)
}

.efa_fit_scalar <- function(fit, name) {
  if (is.null(fit) || is.null(fit[[name]]) || length(fit[[name]]) < 1L) {
    return(NA_real_)
  }

  out <- tryCatch(
    suppressWarnings(as.numeric(fit[[name]][1])),
    error = function(e) NA_real_
  )
  if (length(out) < 1L) {
    return(NA_real_)
  }

  out
}

.efa_format_number <- function(x, digits = 2, print_zero = FALSE, pad = TRUE) {
  if (length(x) != 1L || !is.finite(x)) {
    return("NA")
  }

  .efa_num(round(x, digits = digits), digits = digits,
           print_zero = print_zero, pad = pad)
}

.efa_format_fit_value <- function(fit, name, digits = 2,
                                  print_zero = FALSE, pad = TRUE) {
  .efa_format_number(
    .efa_fit_scalar(fit, name),
    digits = digits,
    print_zero = print_zero,
    pad = pad
  )
}

.efa_bootstrap_sample_text <- function(spec) {
  b_boot <- .efa_as_scalar_number(spec$b_boot)
  if (!is.finite(b_boot)) {
    return("")
  }

  label <- if (isTRUE(spec$is_pooled)) {
    "bootstrap samples per imputation"
  } else {
    "bootstrap samples"
  }

  paste(as.integer(b_boot), label)
}

.efa_alignment_summary <- function(x) {
  alignment <- x$alignment
  if (is.null(alignment)) {
    return("")
  }

  out <- character(0)

  method <- .efa_setting_text(alignment$method)
  if (!nzchar(method)) {
    method <- .efa_setting_text(x$settings$target_method)
  }
  if (nzchar(method)) {
    out <- c(out, paste0("method = '", method, "'"))
  }

  if (!is.null(alignment$converged)) {
    converged <- isTRUE(alignment$converged)
    out <- c(out, if (converged) "converged" else "not converged")
  }

  failures <- alignment$point_rotation_failures
  if (!is.null(failures) && length(failures) > 0L) {
    failures <- failures[!is.na(failures)]
    if (length(failures) > 0L) {
      out <- c(out, paste0("point-alignment failures = ", length(failures)))
    }
  }

  paste(out, collapse = ", ")
}

.efa_is_top_n <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) &&
    (is.infinite(x) || (x >= 1 && x == as.integer(x)))
}

.efa_main_loadings <- function(x, spec) {
  if (identical(spec$rotation, "none")) {
    if (is.null(x$unrot_loadings)) {
      return(NULL)
    }
    return(as.matrix(x$unrot_loadings))
  }

  if (is.null(x$rot_loadings)) {
    return(NULL)
  }

  as.matrix(x$rot_loadings)
}

.efa_variable_names <- function(x) {
  rn <- rownames(x)
  if (is.null(rn)) {
    rn <- paste0("V", seq_len(nrow(x)))
  }
  rn
}

.efa_factor_names <- function(x) {
  cn <- colnames(x)
  if (is.null(cn)) {
    cn <- paste0("F", seq_len(ncol(x)))
  }
  cn
}

.print_efa_diagnostics_section <- function(x, spec, cutoff, digits,
                                           cross_loading_cutoff,
                                           min_primary_gap,
                                           min_salient_per_factor) {
  loadings <- .efa_main_loadings(x, spec)

  title <- if (isTRUE(spec$is_pooled)) {
    "Pooled Model Diagnostics"
  } else {
    "Model Diagnostics"
  }
  .print_efa_rule(title)

  if (is.null(loadings)) {
    cli::cli_text("No loading matrix available for diagnostics.")
    return(invisible(NULL))
  }

  .efa_print_key_value("Factors", ncol(loadings))
  .efa_print_key_value("Variables", nrow(loadings))

  n_text <- .efa_setting_text(spec$N)
  if (nzchar(n_text) && !identical(n_text, "NA")) {
    .efa_print_key_value("N", n_text)
  }

  if (isTRUE(spec$is_pooled)) {
    n_imp <- .efa_setting_text(spec$n_imputations)
    if (nzchar(n_imp)) {
      .efa_print_key_value("Imputations", n_imp)
    }
    pooling <- .efa_pooled_settings_text(spec)
    if (length(pooling) > 0L) {
      .efa_print_key_value("Pooling", paste(pooling, collapse = ", "))
    }

    alignment_text <- .efa_alignment_summary(x)
    if (nzchar(alignment_text)) {
      .efa_print_key_value("Alignment", alignment_text)
    }
  }

  if (isTRUE(spec$np_boot) && !is.null(spec$b_boot)) {
    sample_label <- if (isTRUE(spec$is_pooled)) {
      "Bootstrap samples per imputation"
    } else {
      "Bootstrap samples"
    }
    .efa_print_key_value(sample_label, spec$b_boot)

    valid_text <- .efa_valid_target_rotation_summary(x, spec)
    if (nzchar(valid_text)) {
      .efa_print_key_value("Valid target-rotated samples", valid_text)
    }
  }

  # Distinct local optima the rotation found across its random starts, counted over the
  # starts that converged (only present for the gradient-projection rotations that use
  # random starts; omitted when none converged).
  rot_diag <- spec$rotation_diagnostics
  if (!is.null(rot_diag) && isTRUE(rot_diag$n_converged >= 1L)) {
    .efa_print_key_value(
      "Rotation local optima",
      paste0(rot_diag$n_distinct_minima, " distinct from ",
             rot_diag$n_starts, " random starts")
    )
  }

  # Count Heywood cases the same way as the printed Heywood note (see
  # `.efa_heywood_count()`).
  heywood <- .efa_heywood_count(spec$heywood, spec$h2)
  .efa_print_key_value("Heywood cases", heywood)

  cross_salient <- abs(loadings) >= cross_loading_cutoff
  cross_salient[is.na(cross_salient)] <- FALSE
  cross_loading_items <- sum(rowSums(cross_salient) > 1L)

  salient <- abs(loadings) >= cutoff
  salient[is.na(salient)] <- FALSE
  no_salient_items <- sum(rowSums(salient) < 1L)
  .efa_print_key_value(
    paste0("Cross-loading items (|loading| >= ",
           .efa_format_plain_number(cross_loading_cutoff, digits), ")"),
    cross_loading_items
  )
  .efa_print_key_value(
    paste0("Items without salient loading (|loading| >= ",
           .efa_format_plain_number(cutoff, digits), ")"),
    no_salient_items
  )

  weak_factors <- .efa_factors_below_indicator_count(
    loadings,
    cutoff = cutoff,
    min_salient_per_factor = min_salient_per_factor
  )
  .efa_print_key_value(
    paste0("Factors with fewer than ", min_salient_per_factor,
           " salient indicators"),
    length(weak_factors)
  )

  gap_count <- .efa_small_primary_gap_count(loadings,
    cutoff = cutoff,
    min_primary_gap = min_primary_gap
  )
  .efa_print_key_value(
    paste0("Items with primary-loading gap < ",
           .efa_format_plain_number(min_primary_gap, digits)),
    gap_count
  )

  if (!is.null(x$residuals)) {
    largest_resid <- .efa_largest_abs_residual(x$residuals)
    if (is.finite(largest_resid)) {
      .efa_print_key_value("Largest |residual|",
        .efa_format_plain_number(largest_resid, digits)
      )
    }
  }

  if (!is.null(x$Phi)) {
    high_phi <- .efa_high_factor_correlation_count(x$Phi, cutoff = .85)
    value <- if (high_phi > 0L) {
      as.character(high_phi)
    } else {
      "none"
    }
    .efa_print_key_value("Factor intercorrelations > .85", value)
  }

  invisible(NULL)
}

.efa_print_key_value <- function(key, value) {
  # A "label: value" data line; emit verbatim so a long key or a settings value (e.g. the
  # pooled "Pooling: ..." line) is never reflowed mid-token.
  cli::cli_verbatim(paste0(key, ": ", as.character(value)))
  invisible(NULL)
}

.efa_valid_target_rotation_summary <- function(x, spec) {
  counts <- .efa_valid_target_rotations(x)
  if (is.null(counts) || length(counts) < 1L || all(is.na(counts))) {
    return("")
  }

  counts <- as.integer(stats::na.omit(counts))
  if (length(counts) < 1L) {
    return("")
  }

  total <- .efa_total_bootstrap_samples(x, spec)
  count_text <- .efa_valid_target_rotation_count_text(counts, spec)

  if (!is.na(total)) {
    count_text <- paste0(count_text, " out of ", total)
  }

  if (isTRUE(spec$is_pooled)) {
    count_text <- paste0(count_text, " per imputation")
  }

  count_text
}

.efa_factors_below_indicator_count <- function(loadings, cutoff,
                                               min_salient_per_factor) {
  salient <- abs(loadings) >= cutoff
  salient[is.na(salient)] <- FALSE
  counts <- colSums(salient)
  which(counts < min_salient_per_factor)
}

.efa_small_primary_gap_count <- function(loadings, cutoff, min_primary_gap) {
  length(.efa_small_primary_gap_rows(loadings, cutoff = cutoff,
                                     min_primary_gap = min_primary_gap))
}

.efa_largest_abs_residual <- function(residuals) {
  residuals <- as.matrix(residuals)
  values <- residuals[upper.tri(residuals)]
  values <- values[is.finite(values)]
  if (length(values) < 1L) {
    return(NA_real_)
  }
  max(abs(values), na.rm = TRUE)
}

.efa_high_factor_correlation_count <- function(phi, cutoff = .85) {
  phi <- as.matrix(phi)
  values <- phi[lower.tri(phi)]
  values <- values[is.finite(values)]
  sum(abs(values) > cutoff, na.rm = TRUE)
}

.print_efa_simple_structure_section <- function(x, spec, cutoff, digits,
                                                max_name_length,
                                                cross_loading_cutoff,
                                                min_primary_gap,
                                                min_salient_per_factor,
                                                diagnostics_top_n,
                                                sort_loadings) {
  loadings <- .efa_main_loadings(x, spec)
  row_order <- .efa_loading_row_order(loadings, sort_loadings)
  loadings <- loadings[row_order, , drop = FALSE]

  summary <- .efa_simple_structure_summary(
    loadings = loadings,
    cutoff = cutoff,
    cross_loading_cutoff = cross_loading_cutoff,
    min_primary_gap = min_primary_gap,
    min_salient_per_factor = min_salient_per_factor,
    digits = digits,
    max_name_length = max_name_length
  )

  if (!summary$has_output) {
    return(invisible(NULL))
  }

  .print_efa_rule("Simple Structure Diagnostics")

  if (length(summary$no_salient) > 0L) {
    .efa_print_limited_bullets(
      heading = paste0("Items with no salient loading, |loading| >= ",
                       .efa_format_plain_number(cutoff, digits), ":"),
      values = summary$no_salient,
      top_n = diagnostics_top_n
    )
  }

  if (length(summary$cross_loadings) > 0L) {
    .efa_print_limited_bullets(
      heading = paste0("Items with cross-loadings, |loading| >= ",
                       .efa_format_plain_number(cross_loading_cutoff, digits),
                       " on multiple factors:"),
      values = summary$cross_loadings,
      top_n = diagnostics_top_n
    )
  }

  if (length(summary$small_gaps) > 0L) {
    .efa_print_limited_bullets(
      heading = paste0("Items with primary-loading gap < ",
                       .efa_format_plain_number(min_primary_gap, digits), ":"),
      values = summary$small_gaps,
      top_n = diagnostics_top_n
    )
  }

  if (length(summary$weak_factors) > 0L) {
    .efa_print_limited_bullets(
      heading = paste0("Factors with fewer than ", min_salient_per_factor,
                       " salient indicators:"),
      values = summary$weak_factors,
      top_n = Inf
    )
  }

  invisible(NULL)
}

.efa_simple_structure_summary <- function(loadings, cutoff, cross_loading_cutoff,
                                          min_primary_gap,
                                          min_salient_per_factor,
                                          digits, max_name_length) {
  var_names <- .efa_truncate_names(.efa_variable_names(loadings), max_name_length)
  factor_names <- .efa_factor_names(loadings)

  abs_loadings <- abs(loadings)
  salient_for_items <- abs_loadings >= cutoff
  salient_for_items[is.na(salient_for_items)] <- FALSE

  cross_salient <- abs_loadings >= cross_loading_cutoff
  cross_salient[is.na(cross_salient)] <- FALSE

  no_salient <- var_names[rowSums(salient_for_items) < 1L]

  cross_rows <- which(rowSums(cross_salient) > 1L)
  cross_loadings <- vapply(cross_rows, function(i) {
    idx <- which(cross_salient[i, ])
    paste0(var_names[i], ": ",
      .efa_format_loading_pairs(loadings[i, idx], factor_names[idx], digits)
    )
  }, character(1L))

  small_gap_rows <- .efa_small_primary_gap_rows(loadings,
    cutoff = cutoff,
    min_primary_gap = min_primary_gap
  )
  small_gaps <- vapply(small_gap_rows, function(i) {
    ord <- order(abs(loadings[i, ]), decreasing = TRUE)
    idx <- ord[seq_len(min(2L, length(ord)))]
    paste0(var_names[i], ": ",
      .efa_format_loading_pairs(loadings[i, idx], factor_names[idx], digits)
    )
  }, character(1L))

  factor_counts <- colSums(salient_for_items)
  weak_factor_idx <- which(factor_counts < min_salient_per_factor)
  if (length(weak_factor_idx) > 0L) {
    weak_factors <- paste0(
      factor_names[weak_factor_idx],
      ": ",
      factor_counts[weak_factor_idx],
      " salient indicator",
      ifelse(factor_counts[weak_factor_idx] == 1L, "", "s")
    )
  } else {
    weak_factors <- NULL
  }


  list(
    no_salient = no_salient,
    cross_loadings = cross_loadings,
    small_gaps = small_gaps,
    weak_factors = weak_factors,
    has_output = length(no_salient) > 0L ||
      length(cross_loadings) > 0L ||
      length(small_gaps) > 0L ||
      length(weak_factors) > 0L
  )
}

.efa_small_primary_gap_rows <- function(loadings, cutoff, min_primary_gap) {
  if (ncol(loadings) < 2L) {
    return(integer(0))
  }

  abs_loadings <- abs(loadings)
  abs_loadings[is.na(abs_loadings)] <- -Inf
  sorted <- t(apply(abs_loadings, 1L, sort, decreasing = TRUE))
  primary <- sorted[, 1L]
  secondary <- sorted[, 2L]

  which(is.finite(primary) & primary >= cutoff &
          is.finite(secondary) & (primary - secondary) < min_primary_gap)
}

.efa_format_loading_pairs <- function(values, factor_names, digits) {
  ord <- order(abs(values), decreasing = TRUE)
  values <- values[ord]
  factor_names <- factor_names[ord]

  paste0(
    factor_names,
    " = ",
    vapply(values, .efa_format_plain_number, character(1L), digits = digits),
    collapse = ", "
  )
}

.efa_print_limited_bullets <- function(heading, values, top_n) {
  cli::cli_text("{heading}")

  if (length(values) < 1L) {
    return(invisible(NULL))
  }

  n_print <- if (is.finite(top_n)) {
    min(length(values), as.integer(top_n))
  } else {
    length(values)
  }

  bullets <- values[seq_len(n_print)]
  if (n_print < length(values)) {
    bullets <- c(bullets, paste0("... ", length(values) - n_print, " more not shown"))
  }

  .efa_emit_bullets(bullets)
  cli::cli_text("")
  invisible(NULL)
}

.print_efa_structure_section <- function(x, cutoff, digits, max_name_length,
                                         sort_loadings, max_factors_per_block,
                                         ...) {
  .print_efa_rule("Structure Matrix")
  .efa_emit_lines(.efa_capture_loadings(x$Structure,
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    h2 = NULL,
    sort_loadings = sort_loadings,
    legend = FALSE,
    max_factors_per_block = max_factors_per_block,
    ...
  ))
  invisible(NULL)
}

.print_efa_mi_diagnostics_section <- function(x, spec, digits = 3) {
  if (!isTRUE(spec$is_pooled) || is.null(x$MI)) {
    return(invisible(NULL))
  }

  fmi_values <- .efa_collect_mi_values(x$MI, pattern = "fmi")
  riv_values <- .efa_collect_mi_values(x$MI, pattern = "riv")

  if (length(fmi_values) < 1L && length(riv_values) < 1L) {
    return(invisible(NULL))
  }

  .print_efa_rule("MI Uncertainty Summary")

  if (length(fmi_values) > 0L) {
    .efa_print_key_value("Largest available FMI",
      .efa_format_plain_number(max(fmi_values, na.rm = TRUE), digits)
    )
    .efa_print_key_value("Median available FMI",
      .efa_format_plain_number(stats::median(fmi_values, na.rm = TRUE), digits)
    )
  }

  if (length(riv_values) > 0L) {
    .efa_print_key_value("Largest available RIV",
      .efa_format_plain_number(max(riv_values, na.rm = TRUE), digits)
    )
    .efa_print_key_value("Median available RIV",
      .efa_format_plain_number(stats::median(riv_values, na.rm = TRUE), digits)
    )
  }

  invisible(NULL)
}

.efa_collect_mi_values <- function(x, pattern, path = "") {
  pieces <- .efa_collect_mi_values_impl(x, pattern = pattern, path = path)
  if (length(pieces) < 1L) {
    return(numeric(0))
  }

  values <- unlist(pieces, use.names = FALSE)
  values[is.finite(values)]
}

.efa_collect_mi_values_impl <- function(x, pattern, path = "") {
  if (is.list(x)) {
    nm <- names(x)
    if (is.null(nm)) {
      nm <- rep("", length(x))
    }

    pieces <- vector("list", length(x))
    for (i in seq_along(x)) {
      child_path <- paste(c(path, nm[i]), collapse = "/")
      pieces[[i]] <- .efa_collect_mi_values_impl(x[[i]], pattern, child_path)
    }
    return(unlist(pieces, recursive = FALSE, use.names = FALSE))
  }

  if (!is.numeric(x) || !grepl(pattern, path, ignore.case = TRUE)) {
    return(list())
  }

  list(as.numeric(x))
}

.efa_loading_row_order <- function(x, sort_loadings = c("none", "primary", "clustered")) {
  sort_loadings <- match.arg(sort_loadings)

  if (identical(sort_loadings, "none") || nrow(x) < 2L || ncol(x) < 1L) {
    return(seq_len(nrow(x)))
  }

  abs_x <- abs(x)
  abs_x[is.na(abs_x)] <- -Inf
  primary_factor <- max.col(abs_x, ties.method = "first")
  primary_loading <- apply(abs_x, 1L, max, na.rm = TRUE)
  primary_loading[!is.finite(primary_loading)] <- -Inf

  if (identical(sort_loadings, "primary")) {
    return(order(primary_factor, seq_along(primary_factor)))
  }

  order(primary_factor, -primary_loading, seq_along(primary_factor))
}

.efa_format_plain_number <- function(x, digits = 3, pad = FALSE) {
  if (length(x) != 1L || is.na(x)) {
    return("NA")
  }

  .efa_num(round(x, digits = digits), digits = digits, pad = pad)
}


.efa_paste_gof_ci <- function(fitind_ci, name) {
  if (!.efa_is_ci_pair(fitind_ci)) {
    return("")
  }

  lower_names <- names(fitind_ci$lower)
  upper_names <- names(fitind_ci$upper)
  if (is.null(lower_names) || is.null(upper_names) ||
      !(name %in% lower_names) || !(name %in% upper_names)) {
    return("")
  }

  tryCatch(.paste_gof_ci(fitind_ci, name), error = function(e) "")
}
