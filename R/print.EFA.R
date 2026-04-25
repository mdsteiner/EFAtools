#' Print EFA object
#'
#' Print Method showing a summarized output of the \link{EFA} function
#'
#' @param x list. An object of class EFA to be printed
#' @param cutoff numeric. Passed to \code{\link[EFAtools:print.LOADINGS]{print.LOADINGS}}.
#' The number above which to print loadings in bold. Default is .3.
#' @param digits numeric. Passed to \code{\link[EFAtools:print.LOADINGS]{print.LOADINGS}}
#' Number of digits to round the loadings to (default is 3).
#' @param max_name_length numeric. Passed to \code{\link[EFAtools:print.LOADINGS]{print.LOADINGS}}.
#' The maximum length of the variable names to display. Everything beyond this
#' will be cut from the right.
#' @param ci character. Whether to print confidence intervals for loadings and
#' factor intercorrelations, if available. \code{"auto"} and \code{"separate"}
#' print separate CI sections when CIs were computed; \code{"none"} suppresses
#' these sections. Default is \code{"auto"}.
#' @param ci_filter character. Which loading CIs to print. \code{"salient"}
#' prints CIs for loadings with absolute values greater than or equal to
#' \code{cutoff}; \code{"all"} prints all loading CIs; \code{"nonzero"} prints
#' loading CIs that exclude zero. Default is \code{"salient"}.
#' @param details character. Amount of information to print. \code{"standard"}
#' prints the regular output, \code{"compact"} suppresses longer diagnostic
#' sections, and \code{"full"} additionally prints available optional sections.
#' @param diagnostics logical. Whether to print model and simple-structure
#' diagnostics. Default is \code{TRUE}.
#' @param diagnostics_top_n numeric. Maximum number of item-level diagnostic
#' entries to print per diagnostic type.
#' @param residual_cutoff numeric. Absolute residual cutoff used in the residual
#' diagnostics section. Default is .1.
#' @param residual_top_n numeric. Maximum number of residuals to print. Use
#' \code{Inf} to print all residuals above \code{residual_cutoff}.
#' @param show_structure logical. Whether to print the structure matrix for
#' oblique solutions when available. Default is \code{FALSE}.
#' @param sort_loadings character. Optional row sorting for loading tables.
#' See \code{\link[EFAtools:print.LOADINGS]{print.LOADINGS}}.
#' @param show_loading_legend logical. Whether to print a compact legend for
#' loading-table styling. Default is \code{TRUE}.
#' @param cross_loading_cutoff numeric. Cutoff for counting cross-loadings in
#' the diagnostics. Defaults to \code{cutoff}.
#' @param min_primary_gap numeric. Minimum desired absolute difference between
#' the largest and second-largest absolute loading of an item. Used only for
#' descriptive diagnostics.
#' @param min_salient_per_factor numeric. Minimum number of salient indicators
#' per factor used in the diagnostics. Default is 3.
#' @param max_factors_per_block numeric or \code{NULL}. Maximum number of factor
#' columns per loading-table block. If \code{NULL}, this is chosen from the
#' console width.
#' @param show_mi_diagnostics logical or \code{NULL}. Whether to print a compact
#' MI uncertainty summary for pooled EFAs when available. \code{NULL} prints it
#' only for \code{details = "full"}.
#' @param ... Further arguments passed to \code{\link[EFAtools:print.LOADINGS]{print.LOADINGS}}.
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
#' # With SEs and CIs, needs raw-data
#' mod <- EFA(GRiPS_raw, n_factors = 1, method = "ML", se = "np-boot")
#' mod
#' # with additional infos
#' print(mod, details = "full")
#' # omit some diagnostic parts
#' print(mod, details = "compact")
#'
print.EFA <- function(x, cutoff = .3, digits = 3, max_name_length = 10,
                      ci = c("auto", "none", "separate"),
                      ci_filter = c("salient", "all", "nonzero"),
                      details = c("standard", "compact", "full"),
                      diagnostics = TRUE,
                      diagnostics_top_n = 10,
                      residual_cutoff = .1,
                      residual_top_n = 10,
                      show_structure = FALSE,
                      sort_loadings = c("none", "primary", "clustered"),
                      show_loading_legend = TRUE,
                      cross_loading_cutoff = cutoff,
                      min_primary_gap = .20,
                      min_salient_per_factor = 3,
                      max_factors_per_block = NULL,
                      show_mi_diagnostics = NULL, ...) {
  ci <- match.arg(ci)
  ci_filter <- match.arg(ci_filter)
  details <- match.arg(details)
  sort_loadings <- match.arg(sort_loadings)

  .print_EFA(x,
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    ci = ci,
    ci_filter = ci_filter,
    details = details,
    diagnostics = diagnostics,
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
    show_mi_diagnostics = show_mi_diagnostics,
    ...
  )

  invisible(x)
}

#' @rdname print.EFA
#' @export
#' @method print EFA_POOLED
print.EFA_POOLED <- function(x, cutoff = .3, digits = 3, max_name_length = 10,
                             ci = c("auto", "none", "separate"),
                             ci_filter = c("salient", "all", "nonzero"),
                             details = c("standard", "compact", "full"),
                             diagnostics = TRUE,
                             diagnostics_top_n = 10,
                             residual_cutoff = .1,
                             residual_top_n = 10,
                             show_structure = FALSE,
                             sort_loadings = c("none", "primary", "clustered"),
                             show_loading_legend = TRUE,
                             cross_loading_cutoff = cutoff,
                             min_primary_gap = .20,
                             min_salient_per_factor = 3,
                             max_factors_per_block = NULL,
                             show_mi_diagnostics = NULL, ...) {
  ci <- match.arg(ci)
  ci_filter <- match.arg(ci_filter)
  details <- match.arg(details)
  sort_loadings <- match.arg(sort_loadings)

  .print_EFA(x,
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    ci = ci,
    ci_filter = ci_filter,
    details = details,
    diagnostics = diagnostics,
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
    show_mi_diagnostics = show_mi_diagnostics,
    ...
  )

  invisible(x)
}

#' @rdname print.EFA
#' @export
#' @method format EFA
format.EFA <- function(x, cutoff = .3, digits = 3, max_name_length = 10,
                       ci = c("auto", "none", "separate"),
                       ci_filter = c("salient", "all", "nonzero"),
                       details = c("standard", "compact", "full"),
                       diagnostics = TRUE,
                       diagnostics_top_n = 10,
                       residual_cutoff = .1,
                       residual_top_n = 10,
                       show_structure = FALSE,
                       sort_loadings = c("none", "primary", "clustered"),
                       show_loading_legend = TRUE,
                       cross_loading_cutoff = cutoff,
                       min_primary_gap = .20,
                       min_salient_per_factor = 3,
                       max_factors_per_block = NULL,
                       show_mi_diagnostics = NULL, ...) {
  ci <- match.arg(ci)
  ci_filter <- match.arg(ci_filter)
  details <- match.arg(details)
  sort_loadings <- match.arg(sort_loadings)

  utils::capture.output(
    .print_EFA(x,
      cutoff = cutoff,
      digits = digits,
      max_name_length = max_name_length,
      ci = ci,
      ci_filter = ci_filter,
      details = details,
      diagnostics = diagnostics,
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
      show_mi_diagnostics = show_mi_diagnostics,
      ...
    )
  )
}

#' @rdname print.EFA
#' @export
#' @method format EFA_POOLED
format.EFA_POOLED <- function(x, cutoff = .3, digits = 3, max_name_length = 10,
                              ci = c("auto", "none", "separate"),
                              ci_filter = c("salient", "all", "nonzero"),
                              details = c("standard", "compact", "full"),
                              diagnostics = TRUE,
                              diagnostics_top_n = 10,
                              residual_cutoff = .1,
                              residual_top_n = 10,
                              show_structure = FALSE,
                              sort_loadings = c("none", "primary", "clustered"),
                              show_loading_legend = TRUE,
                              cross_loading_cutoff = cutoff,
                              min_primary_gap = .20,
                              min_salient_per_factor = 3,
                              max_factors_per_block = NULL,
                              show_mi_diagnostics = NULL, ...) {
  ci <- match.arg(ci)
  ci_filter <- match.arg(ci_filter)
  details <- match.arg(details)
  sort_loadings <- match.arg(sort_loadings)

  utils::capture.output(
    .print_EFA(x,
      cutoff = cutoff,
      digits = digits,
      max_name_length = max_name_length,
      ci = ci,
      ci_filter = ci_filter,
      details = details,
      diagnostics = diagnostics,
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
      show_mi_diagnostics = show_mi_diagnostics,
      ...
    )
  )
}

.print_EFA <- function(x, cutoff = .3, digits = 3, max_name_length = 10,
                       ci = c("auto", "none", "separate"),
                       ci_filter = c("salient", "all", "nonzero"),
                       details = c("standard", "compact", "full"),
                       diagnostics = TRUE,
                       diagnostics_top_n = 10,
                       residual_cutoff = .1,
                       residual_top_n = 10,
                       show_structure = FALSE,
                       sort_loadings = c("none", "primary", "clustered"),
                       show_loading_legend = TRUE,
                       cross_loading_cutoff = cutoff,
                       min_primary_gap = .20,
                       min_salient_per_factor = 3,
                       max_factors_per_block = NULL,
                       show_mi_diagnostics = NULL, ...) {
  ci <- match.arg(ci)
  ci_filter <- match.arg(ci_filter)
  details <- match.arg(details)
  sort_loadings <- match.arg(sort_loadings)

  .efa_validate_print_options(
    diagnostics = diagnostics,
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
  show_structure <- .efa_resolve_show_structure(show_structure, details)
  show_mi_diagnostics <- .efa_resolve_show_mi_diagnostics(
    show_mi_diagnostics,
    details = details,
    spec = spec
  )

  .print_efa_header(spec)

  if (isTRUE(diagnostics) && !identical(details, "compact")) {
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
    ci = ci,
    ci_filter = ci_filter,
    show_structure = show_structure,
    sort_loadings = sort_loadings,
    show_loading_legend = show_loading_legend,
    max_factors_per_block = max_factors_per_block,
    ...
  )

  if (isTRUE(diagnostics) && !identical(details, "compact")) {
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

  if (!identical(details, "compact")) {
    .print_efa_variances_section(x, spec)
  }

  if (!.print_efa_identification_warning(spec)) {
    return(invisible(x))
  }

  .print_efa_model_fit_section(x, spec)
  .print_efa_bootstrap_note(x, spec)

  if (isTRUE(show_mi_diagnostics)) {
    .print_efa_mi_diagnostics_section(x, spec, digits = digits)
  }

  if (!identical(details, "compact")) {
    .print_efa_residuals_section(x,
      residual_cutoff = residual_cutoff,
      residual_top_n = residual_top_n,
      digits = digits
    )
  }

  invisible(x)
}

.efa_print_spec <- function(x) {
  se <- x$settings$se
  is_pooled <- inherits(x, "EFA_POOLED") || isTRUE(x$settings$pooled)

  list(
    method = x$settings$method,
    rotation = x$settings$rotation,
    type = x$settings$type,
    N = x$settings$N,
    fit = x$fit_indices,
    h2 = x$h2,
    np_boot = isTRUE(se == "np-boot"),
    ci = .efa_ci_level(x),
    b_boot = x$settings$b_boot,
    max_iter = x$settings$max_iter,
    iter = x$iter,
    is_pooled = is_pooled,
    n_imputations = .efa_n_imputations(x),
    target_method = x$settings$target_method,
    align_unrotated = x$settings$align_unrotated,
    fit_pool_method = x$settings$fit_pool_method,
    has_boot_ci = !is.null(x$boot.CI)
  )
}

.print_efa_header <- function(spec) {
  cat("\n")

  if (isTRUE(spec$is_pooled)) {
    .print_efa_pooled_header(spec)
  } else {
    cat("EFA performed with type = '", crayon::bold(spec$type),
      "', method = '", crayon::bold(spec$method),
      "', and rotation = '", crayon::bold(spec$rotation), "'.",
      sep = ""
    )
    cat("\n")
  }

  if (!is.null(spec$max_iter) && spec$iter > spec$max_iter) {
    cat("\n")
    cat(crayon::red$bold(
      cli::symbol$cross,
      "Maximum number of iterations reached",
      "without convergence"
    ))
    cat("\n")
  }
}


.print_efa_pooled_header <- function(spec) {
  cat("Pooled EFA")

  n_imputations <- .efa_setting_text(spec$n_imputations)
  if (nzchar(n_imputations)) {
    cat(" across ", crayon::bold(n_imputations), " imputations", sep = "")
  }

  cat(" performed with type = '", crayon::bold(spec$type),
    "', method = '", crayon::bold(spec$method),
    "', and rotation = '", crayon::bold(spec$rotation), "'.",
    sep = ""
  )
  cat("\n")

  pooling_settings <- .efa_pooled_settings_text(spec)
  if (length(pooling_settings) > 0L) {
    cat("Pooling settings: ", paste(pooling_settings, collapse = ", "), ".", sep = "")
    cat("\n")
  }

  invisible(NULL)
}

.efa_pooled_settings_text <- function(spec) {
  out <- character(0)

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

  out
}

.print_efa_loadings_section <- function(x, spec, cutoff, digits, max_name_length,
                                        ci, ci_filter, show_structure,
                                        sort_loadings, show_loading_legend,
                                        max_factors_per_block, ...) {
  if (identical(spec$rotation, "none")) {
    .print_efa_rule("Unrotated Loadings")
    print(x$unrot_loadings,
      cutoff = cutoff,
      digits = digits,
      max_name_length = max_name_length,
      h2 = x$h2,
      sort_loadings = sort_loadings,
      legend = show_loading_legend,
      max_factors_per_block = max_factors_per_block,
      ...
    )
    .print_efa_heywood_warning(spec$h2)
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
  print(x$rot_loadings,
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    h2 = x$h2,
    sort_loadings = sort_loadings,
    legend = show_loading_legend,
    max_factors_per_block = max_factors_per_block,
    ...
  )
  .print_efa_heywood_warning(spec$h2)
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
    .print_efa_rule("Factor Intercorrelations")
    cat(.get_compare_matrix(x$Phi,
      r_red = Inf,
      n_char = 17,
      var_names = paste0("F", seq_len(ncol(x$Phi)))
    ))
    .print_efa_phi_ci_section(x, spec, digits, ci)
  }

  if (!is.null(x$Structure)) {
    if (isTRUE(show_structure)) {
      .print_efa_structure_section(x,
        cutoff = cutoff,
        digits = digits,
        max_name_length = max_name_length,
        sort_loadings = sort_loadings,
        max_factors_per_block = max_factors_per_block,
        ...
      )
    }
    # else {
    #   .print_efa_structure_note()
    # }
  }

  invisible(NULL)
}

.print_efa_variances_section <- function(x, spec) {
  .print_efa_rule("Variances Accounted for")

  if (identical(spec$rotation, "none")) {
    vars_accounted <- x$vars_accounted
  } else {
    vars_accounted <- x$vars_accounted_rot
  }

  cat(.get_compare_matrix(vars_accounted, r_red = Inf, n_char = 17))

  invisible(NULL)
}

.print_efa_identification_warning <- function(spec) {
  fit <- spec$fit

  if (fit$df == 0) {
    cat("\n")
    cat(
      crayon::yellow$bold("!"),
      crayon::yellow(
        " The model is just identified (df = 0). Goodness of fit indices may not be interpretable."
      )
    )
    cat("\n")
    return(TRUE)
  }

  if (fit$df < 0) {
    cat("\n")
    cat(
      crayon::yellow$bold("!"),
      crayon::yellow(
        " The model is underidentified (df < 0). No goodness of fit indices were calculated."
      )
    )
    cat("\n")
    return(FALSE)
  }

  TRUE
}

.print_efa_model_fit_section <- function(x, spec) {
  fit <- spec$fit
  fit_ci <- .efa_fit_ci_labels(x, spec)

  .print_efa_rule("Model Fit")

  if (identical(spec$method, "PAF") || is.na(spec$N)) {
    cat(crayon::blue("CAF", fit_ci$label, ":"),
      .numformat(fit$CAF), fit_ci$CAF, "\n",
      sep = ""
    )
    cat(crayon::blue("RMSR", fit_ci$label, ":"),
      .numformat(fit$RMSR), fit_ci$RMSR, "\n",
      sep = ""
    )
    cat(crayon::blue("df: "),
      .numformat(fit$df, 0, print_zero = TRUE), "\n",
      sep = ""
    )
    return(invisible(NULL))
  }

  cat(crayon::blue("\U1D712\U00B2(", sep = ""), fit$df,
    crayon::blue(") = ", sep = ""),
    .numformat(fit$chi, 2, print_zero = TRUE), ", ",
    crayon::blue(crayon::italic("p")),
    ifelse(
      fit$p_chi < .001,
      " < .001",
      paste(
        crayon::blue(ifelse(fit$p_chi < 1, " =", " = ")),
        .numformat(fit$p_chi, 3),
        sep = ""
      )
    ),
    "\n",
    sep = ""
  )

  cat(crayon::blue("CFI", fit_ci$label, ": "),
    .numformat(fit$CFI, pad = FALSE), fit_ci$CFI, "\n",
    sep = ""
  )
  cat(crayon::blue("RMSEA [90% CI]", fit_ci$label, ": "),
    paste0(
      .numformat(fit$RMSEA), " [",
      ifelse(fit$RMSEA_LB < 1, substr(.numformat(fit$RMSEA_LB), 2, 4),
        .numformat(fit$RMSEA_LB)
      ),
      ifelse(fit$RMSEA_UB < 1, ";", "; "),
      .numformat(fit$RMSEA_UB), "]",
      sep = ""
    ),
    fit_ci$RMSEA, "\n",
    sep = ""
  )
  cat(crayon::blue("AIC", fit_ci$label, ": "),
    .numformat(fit$AIC, print_zero = TRUE), fit_ci$AIC, "\n",
    sep = ""
  )
  cat(crayon::blue("BIC", fit_ci$label, ": "),
    .numformat(fit$BIC, print_zero = TRUE), fit_ci$BIC, "\n",
    sep = ""
  )
  cat(crayon::blue("CAF", fit_ci$label, ":"),
    .numformat(fit$CAF), fit_ci$CAF, "\n",
    sep = ""
  )
  cat(crayon::blue("RMSR", fit_ci$label, ":"),
    .numformat(fit$RMSR), fit_ci$RMSR, "\n",
    sep = ""
  )

  invisible(NULL)
}

.efa_fit_ci_labels <- function(x, spec) {
  out <- list(
    label = "",
    CAF = "",
    RMSR = "",
    CFI = "",
    RMSEA = "",
    AIC = "",
    BIC = ""
  )

  if (!isTRUE(spec$np_boot)) {
    return(out)
  }

  fitind_ci <- .efa_get_fit_ci(x)
  if (is.null(fitind_ci)) {
    return(out)
  }

  out$label <- paste0(" [", .efa_ci_header_label(spec), "-CI]")
  out$CAF <- .efa_paste_gof_ci(fitind_ci, "CAF")
  out$RMSR <- .efa_paste_gof_ci(fitind_ci, "RMSR")
  out$CFI <- .efa_paste_gof_ci(fitind_ci, "CFI")
  out$RMSEA <- .efa_paste_gof_ci(fitind_ci, "RMSEA")
  out$AIC <- .efa_paste_gof_ci(fitind_ci, "AIC")
  out$BIC <- .efa_paste_gof_ci(fitind_ci, "BIC")

  out
}

.print_efa_bootstrap_note <- function(x, spec) {
  if (!isTRUE(spec$np_boot)) {
    return(invisible(NULL))
  }

  note_label <- if (isTRUE(spec$is_pooled)) {
    "Bootstrap/MI CIs"
  } else {
    "Bootstrap CIs"
  }

  cat("\n", crayon::italic("Note: "), note_label, " based on ",
    spec$b_boot, " Bootstrap samples",
    sep = ""
  )

  target_rotation_note <- .efa_target_rotation_note(x, spec)
  if (nzchar(target_rotation_note)) {
    cat("; ", target_rotation_note, sep = "")
  }

  cat(".\n")

  invisible(NULL)
}

.print_efa_residuals_section <- function(x, residual_cutoff = .1,
                                     residual_top_n = 10, digits = 3) {
  .print_efa_rule("Residual Diagnostics")

  if (is.null(x$residuals)) {
    cat("No residual matrix available.\n")
    return(invisible(NULL))
  }

  residuals <- as.matrix(x$residuals)
  idx <- which(upper.tri(residuals) & is.finite(residuals), arr.ind = TRUE)

  if (nrow(idx) < 1L) {
    cat("No finite off-diagonal residuals available.\n")
    return(invisible(NULL))
  }

  values <- residuals[idx]
  keep <- abs(values) > residual_cutoff
  n_large <- sum(keep, na.rm = TRUE)
  largest_abs <- max(abs(values), na.rm = TRUE)

  cat(crayon::blue("Residual cutoff: "), "|r| > ",
    .efa_format_plain_number(residual_cutoff, digits), "\n",
    sep = ""
  )
  cat(crayon::blue("Number of large residuals: "), n_large, "\n", sep = "")
  cat(crayon::blue("Largest absolute residual: "),
    .efa_format_plain_number(largest_abs, digits), "\n",
    sep = ""
  )

  if (n_large < 1L) {
    cat("\nNo absolute residuals > ",
      .efa_format_plain_number(residual_cutoff, digits),
      " occurred.\n",
      sep = ""
    )
    cat("\nInspect the residual matrix for details (e.g., with residuals()).\n")
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

  var_names <- colnames(residuals)
  if (is.null(var_names)) {
    var_names <- paste0("V", seq_len(ncol(residuals)))
  }

  cat("\n")
  heading <- if (n_print < n_large) {
    paste0("Largest residuals (top ", n_print, " of ", n_large, "):")
  } else {
    "Largest residuals:"
  }
  cat(crayon::blue(heading), "\n", sep = "")

  for (i in seq_along(large_values)) {
    cat(cli::symbol$bullet, " ",
      var_names[large_idx[i, 1]], " ~~ ", var_names[large_idx[i, 2]],
      ": ", .efa_format_plain_number(large_values[i], digits), "\n",
      sep = ""
    )
  }

  cat("\nInspect the residual matrix for details (e.g., with residuals()).\n")

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
    cat("No loading CIs matched ci_filter = '", ci_filter, "'.\n", sep = "")
    return(invisible(NULL))
  }

  idx <- which(keep, arr.ind = TRUE)
  row_names <- rownames(loadings)
  if (is.null(row_names)) {
    row_names <- paste0("V", seq_len(nrow(loadings)))
  }
  row_names <- .efa_truncate_names(row_names, max_name_length)

  factor_names <- colnames(loadings)
  if (is.null(factor_names)) {
    factor_names <- paste0("F", seq_len(ncol(loadings)))
  }

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
  factor_names <- colnames(phi)
  if (is.null(factor_names)) {
    factor_names <- paste0("F", seq_len(ncol(phi)))
  }

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

.print_efa_heywood_warning <- function(h2) {
  n_heywood <- sum(h2 >= 1 + .Machine$double.eps)

  if (n_heywood == 1) {
    cat(crayon::red$bold("\nWarning: A Heywood case was detected!"))
    cat("\n")
  } else if (n_heywood > 1) {
    cat(crayon::red$bold("\nWarning:", n_heywood, "Heywood cases were detected!"))
    cat("\n")
  }

  invisible(NULL)
}

.print_efa_rule <- function(title) {
  cat("\n")
  cat(cli::rule(left = crayon::bold(title), col = "blue"))
  cat("\n")
  cat("\n")

  invisible(NULL)
}

.efa_should_print_ci <- function(x, ci) {
  !identical(ci, "none") && !is.null(x$boot.CI)
}

.efa_get_component_ci <- function(x, component) {
  if (is.null(x$boot.CI)) {
    return(NULL)
  }

  ci <- x$boot.CI[[component]]
  if (.efa_is_ci_pair(ci)) {
    return(ci)
  }

  NULL
}

.efa_get_fit_ci <- function(x) {
  if (is.null(x$boot.CI)) {
    return(NULL)
  }

  for (component in c("fit_indices", "fit_indices_pooled_algorithm", "fit_indices_descriptive")) {
    ci <- x$boot.CI[[component]]
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

  as.numeric(val[1])
}

.efa_ci_level_text <- function(level) {
  if (!is.finite(level)) {
    return("")
  }

  paste0(formatC(level * 100, format = "fg", digits = 4, width = -1), "%")
}

.efa_ci_source_text <- function(spec) {
  if (isTRUE(spec$is_pooled)) {
    return("bootstrap/MI")
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
  if (!is.null(x$boot.SE$valid_target_rotations)) {
    return(x$boot.SE$valid_target_rotations)
  }

  if (!is.null(x$boot.MI$bootstrap_rotation_valid)) {
    return(x$boot.MI$bootstrap_rotation_valid)
  }

  NULL
}

.efa_total_bootstrap_samples <- function(x, spec) {
  if (!is.null(spec$b_boot) && length(spec$b_boot) >= 1L &&
      is.finite(spec$b_boot[1])) {
    return(as.integer(spec$b_boot[1]))
  }

  if (!is.null(x$boot.arrays$unrot_loadings)) {
    arr <- x$boot.arrays$unrot_loadings
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

.efa_format_ci_num <- function(x, digits) {
  vapply(x, function(z) {
    .numformat(round(z, digits = digits), digits = digits)
  }, character(1))
}

.print_efa_ci_table <- function(key, second_key = NULL, estimate, lower, upper,
                                digits, key_header, second_key_header = NULL,
                                highlight = NULL) {
  has_second_key <- !is.null(second_key)

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
    stringr::str_pad(headers, widths, side = "right"),
    collapse = "  "
  )
  cat(crayon::blue(header), "\n", sep = "")

  for (i in seq_len(nrow(plain_rows))) {
    cells <- as.character(unlist(plain_rows[i, , drop = FALSE], use.names = FALSE))
    padded <- stringr::str_pad(cells, widths, side = "right")

    if (has_second_key) {
      padded[1] <- crayon::blue(padded[1])
      padded[2] <- crayon::blue(padded[2])
      if (!is.null(highlight) && isTRUE(highlight[i])) {
        padded[3] <- crayon::bold(padded[3])
      }
    } else {
      padded[1] <- crayon::blue(padded[1])
    }

    cat(paste(padded, collapse = "  "), "\n", sep = "")
  }

  invisible(NULL)
}


.efa_validate_print_options <- function(diagnostics, diagnostics_top_n,
                                        residual_cutoff, residual_top_n,
                                        show_structure, show_loading_legend,
                                        cross_loading_cutoff, min_primary_gap,
                                        min_salient_per_factor,
                                        max_factors_per_block,
                                        show_mi_diagnostics) {
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) {
    stop("`diagnostics` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!.efa_is_top_n(diagnostics_top_n)) {
    stop("`diagnostics_top_n` must be a positive integer or Inf.", call. = FALSE)
  }

  if (!is.numeric(residual_cutoff) || length(residual_cutoff) != 1L ||
      is.na(residual_cutoff) || residual_cutoff < 0) {
    stop("`residual_cutoff` must be a single non-negative number.", call. = FALSE)
  }

  if (!.efa_is_top_n(residual_top_n)) {
    stop("`residual_top_n` must be a positive integer or Inf.", call. = FALSE)
  }

  if (!is.logical(show_structure) || length(show_structure) != 1L ||
      is.na(show_structure)) {
    stop("`show_structure` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(show_loading_legend) || length(show_loading_legend) != 1L ||
      is.na(show_loading_legend)) {
    stop("`show_loading_legend` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(cross_loading_cutoff) || length(cross_loading_cutoff) != 1L ||
      is.na(cross_loading_cutoff) || cross_loading_cutoff < 0) {
    stop("`cross_loading_cutoff` must be a single non-negative number.", call. = FALSE)
  }

  if (!is.numeric(min_primary_gap) || length(min_primary_gap) != 1L ||
      is.na(min_primary_gap) || min_primary_gap < 0) {
    stop("`min_primary_gap` must be a single non-negative number.", call. = FALSE)
  }

  if (!is.numeric(min_salient_per_factor) || length(min_salient_per_factor) != 1L ||
      is.na(min_salient_per_factor) || min_salient_per_factor < 1 ||
      min_salient_per_factor != as.integer(min_salient_per_factor)) {
    stop("`min_salient_per_factor` must be a single positive integer.", call. = FALSE)
  }

  if (!is.null(max_factors_per_block) &&
      (!is.numeric(max_factors_per_block) || length(max_factors_per_block) != 1L ||
       is.na(max_factors_per_block) || max_factors_per_block < 1 ||
       max_factors_per_block != as.integer(max_factors_per_block))) {
    stop("`max_factors_per_block` must be NULL or a single positive integer.", call. = FALSE)
  }

  if (!is.null(show_mi_diagnostics) &&
      (!is.logical(show_mi_diagnostics) || length(show_mi_diagnostics) != 1L ||
       is.na(show_mi_diagnostics))) {
    stop("`show_mi_diagnostics` must be TRUE, FALSE, or NULL.", call. = FALSE)
  }

  invisible(TRUE)
}

.efa_is_top_n <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) &&
    (is.infinite(x) || (x >= 1 && x == as.integer(x)))
}

.efa_resolve_show_structure <- function(show_structure, details) {
  isTRUE(show_structure) || identical(details, "full")
}

.efa_resolve_show_mi_diagnostics <- function(show_mi_diagnostics, details, spec) {
  if (!is.null(show_mi_diagnostics)) {
    return(isTRUE(show_mi_diagnostics))
  }

  isTRUE(spec$is_pooled) && identical(details, "full")
}

.efa_main_loadings <- function(x, spec) {
  if (identical(spec$rotation, "none")) {
    return(as.matrix(x$unrot_loadings))
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

  h2 <- spec$h2
  heywood <- if (is.null(h2)) {
    0L
  } else {
    sum(h2 >= 1 + .Machine$double.eps, na.rm = TRUE)
  }
  .efa_print_key_value("Heywood cases", heywood)

  sal <- abs(loadings) >= cross_loading_cutoff
  sal[is.na(sal)] <- FALSE
  salient_per_item <- rowSums(sal)
  cross_loading_items <- sum(salient_per_item > 1L)
  no_salient_items <- sum(salient_per_item < 1L)
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
  cat(crayon::blue(paste0(key, ": ")), value, "\n", sep = "")
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
  if (ncol(loadings) < 2L) {
    return(0L)
  }

  abs_loadings <- abs(loadings)
  abs_loadings[is.na(abs_loadings)] <- -Inf
  sorted <- t(apply(abs_loadings, 1L, sort, decreasing = TRUE))
  primary <- sorted[, 1L]
  secondary <- sorted[, 2L]
  sum(is.finite(primary) & primary >= cutoff &
        is.finite(secondary) & (primary - secondary) < min_primary_gap)
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
  cat(crayon::blue(heading), "\n", sep = "")

  if (length(values) < 1L) {
    return(invisible(NULL))
  }

  n_print <- if (is.finite(top_n)) {
    min(length(values), as.integer(top_n))
  } else {
    length(values)
  }

  for (i in seq_len(n_print)) {
    cat(cli::symbol$bullet, " ", values[i], "\n", sep = "")
  }

  if (n_print < length(values)) {
    cat(cli::symbol$bullet, " ... ", length(values) - n_print,
      " more not shown\n",
      sep = ""
    )
  }

  cat("\n")
  invisible(NULL)
}

.print_efa_structure_section <- function(x, cutoff, digits, max_name_length,
                                         sort_loadings, max_factors_per_block,
                                         ...) {
  .print_efa_rule("Structure Matrix")
  print(x$Structure,
    cutoff = cutoff,
    digits = digits,
    max_name_length = max_name_length,
    h2 = NULL,
    sort_loadings = sort_loadings,
    legend = FALSE,
    max_factors_per_block = max_factors_per_block,
    ...
  )
  invisible(NULL)
}

.print_efa_structure_note <- function() {
  cat("\n")
  cat(crayon::italic("Note: "),
    "Structure coefficients are available; use show_structure = TRUE to print them.",
    "\n",
    sep = ""
  )
  invisible(NULL)
}

.print_efa_mi_diagnostics_section <- function(x, spec, digits = 3) {
  if (!isTRUE(spec$is_pooled) || is.null(x$boot.MI)) {
    return(invisible(NULL))
  }

  fmi_values <- .efa_collect_mi_values(x$boot.MI, pattern = "fmi")
  riv_values <- .efa_collect_mi_values(x$boot.MI, pattern = "riv")

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
  out <- numeric(0)

  if (is.list(x)) {
    nm <- names(x)
    if (is.null(nm)) {
      nm <- rep("", length(x))
    }

    for (i in seq_along(x)) {
      child_path <- paste(c(path, nm[i]), collapse = "/")
      out <- c(out, .efa_collect_mi_values(x[[i]], pattern, child_path))
    }
    return(out)
  }

  if (!is.numeric(x) || !grepl(pattern, path, ignore.case = TRUE)) {
    return(out)
  }

  values <- as.numeric(x)
  values <- values[is.finite(values)]
  values
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

  .numformat(round(x, digits = digits), digits = digits, pad = pad)
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
