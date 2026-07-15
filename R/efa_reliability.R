#' Reliability and common-variance coefficients for a factor solution
#'
#' Computes model-based reliability coefficients (McDonald's omega total,
#' hierarchical, and subscale; standardized Cronbach's alpha; and the H index)
#' together with the bifactor common-variance indices (explained common variance,
#' ECV, and percent of uncontaminated correlations, PUC) for the general factor
#' and each group factor, and returns them as a tidy, long-format table. The
#' coefficients can be obtained from a Schmid-Leiman solution
#' ([efa_schmid_leiman()] or [psych::schmid()]), an oblique [efa_fit()]
#' (correlated-factors) solution, a `lavaan` single-factor, second-order, or
#' bifactor fit, a raw bifactor loading matrix, or manually supplied components.
#'
#' @details ## Coefficients
#'
#' The reliability coefficients are McDonald's omegas (McDonald, 1978, 1985,
#' 1999), standardized Cronbach's alpha (Cronbach, 1951), and the H index
#' (construct replicability; Hancock & Mueller, 2001). The omegas give the share
#' of true score variance in a unit-weighted composite: omega total the share due
#' to all factors together, omega hierarchical the share due to the general
#' factor only, and omega subscale the share due to the group factors (for the
#' whole scale) or to the specific group factor (for a subscale composite).
#' Alpha is computed from the standardized correlation matrix. The H index is
#' the reliability of an optimally weighted composite; low values indicate a
#' factor that is not well defined by its indicators.
#'
#' The common-variance indices ECV and PUC (Bonifay et al., 2015; Reise et al.,
#' 2013; Rodriguez et al., 2016a, 2016b) describe the general factor and so are
#' reported for the general factor only. The ECV is the share of the common
#' variance that is explained by the general factor. The PUC is the proportion
#' of correlations that reflect general-factor variance alone (those between
#' indicators of different group factors); the higher it is, the more similar
#' the general factor is to the single factor of a unidimensional model.
#'
#' ## Input
#'
#' An oblique [efa_fit()] object is scored as the correlated-factors model it is
#' (having no general factor, it omits the bifactor indices -- omega hierarchical,
#' ECV, and PUC), and a bare loading matrix is read as a raw bifactor solution
#' (general factor in the first column). For a correlated-factors [efa_fit()] solution `variance` is always
#' `"correlation"`. The indicator-to-factor correspondences come from `factor_map`
#' when it is supplied; otherwise each variable is assigned to the group factor on
#' which it loads most strongly. For `lavaan` input the composite variances are
#' model-implied (`variance` is not used), and the coefficients are computed per
#' group.
#'
#' @param model an [efa_schmid_leiman()], `schmid` ([psych::schmid()]), [efa_fit()]
#'   (oblique), or `lavaan` object; a raw bifactor loading matrix (general factor
#'   first); or `NULL` to supply the components manually via `g_load`, `s_load`,
#'   `u2`, and `var_names`.
#' @param coefficients character. An optional subset of the coefficients to
#'   return, any of `"omega_total"`, `"omega_hierarchical"`, `"omega_subscale"`,
#'   `"alpha"`, `"H"`, `"ECV"`, and `"PUC"`. Default `NULL` returns all of them.
#' @param g_name character. The name of the general factor in the `lavaan`
#'   solution. Only needed for a `lavaan` second-order or bifactor fit. Default is
#'   `"g"`.
#' @param group_names character. An optional vector of group names for a `lavaan`
#'   multiple-group fit. Its length must match the number of groups.
#' @param factor_map matrix. A logical or 0/1 matrix indicating which variable
#'   corresponds to which group factor, with the same dimensions as the group
#'   loading matrix (cross-loadings are allowed). If `NULL` (default), each
#'   variable is assigned to the group factor on which it loads most strongly. Not
#'   used for `lavaan` input.
#' @param variance character. The total-variance denominator for the coefficients:
#'   `"correlation"` (default) takes it from the correlation matrix (the observed-score
#'   omega, as in [psych::omega()]); `"sums_load"` uses the model-implied composite
#'   variance from the squared loading sums and the uniquenesses, which needs no
#'   correlation matrix and so is the way to score a bare loading matrix or
#'   manual components given without one. Some inputs fix the convention: `lavaan` is
#'   always model-implied and an oblique [efa_fit()] is always correlation-based.
#' @param var_names character. Subtest names in the row order of the loadings.
#'   Only needed when `model` is `NULL`.
#' @param fac_names character. An optional vector of group-factor names in the
#'   column order of the loadings. Taken from the input if `NULL`.
#' @param g_load numeric. General-factor loadings. Only needed when `model` is
#'   `NULL`.
#' @param s_load matrix. Group-factor loadings. Only needed when `model` is
#'   `NULL`.
#' @param u2 numeric. Uniquenesses. Only needed when `model` is `NULL` (or to
#'   override the communality-based default for a raw bifactor matrix).
#' @param cormat matrix. A correlation matrix used when `variance =
#'   "correlation"`. If `NULL`, it is taken from the input object or reconstructed
#'   from `pattern` and `Phi` where possible.
#' @param pattern matrix. Pattern coefficients from an oblique solution, used with
#'   `Phi` to reconstruct `cormat` when `model` is `NULL` and no `cormat` is given.
#' @param Phi matrix. Factor intercorrelations from an oblique solution, used with
#'   `pattern`.
#'
#' @returns An object of class `efa_reliability`: a long-format data frame with
#'   one row per computed coefficient, with columns
#'   \item{coefficient}{the coefficient name (e.g. `"omega_total"`).}
#'   \item{level}{`"general"` for the general-factor row, `"group"` otherwise.}
#'   \item{factor}{the factor label (`"g"` for the general factor).}
#'   \item{group}{the group label, or `NA` for a single ungrouped solution.}
#'   \item{value}{the coefficient value.}
#'   Structurally undefined cells (for example ECV and PUC on a group factor) are
#'   omitted. The object carries a `settings` attribute (the total-variance
#'   convention used) and a `kind` attribute tagging each coefficient as a
#'   reliability coefficient or a common-variance index, and has a [print()]
#'   method.
#'
#' @source McDonald, R. P. (1978). Generalizability in factorable domains: Domain
#'   validity and generalizability. Educational and Psychological Measurement, 38,
#'   75-79.
#' @source McDonald, R. P. (1985). Factor analysis and related methods. Hillsdale,
#'   NJ: Erlbaum.
#' @source McDonald, R. P. (1999). Test theory: A unified treatment. Mahwah,
#'   NJ: Erlbaum.
#' @source Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
#'   tests. Psychometrika, 16, 297-334.
#' @source Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
#'   reliability within latent variable systems. In R. Cudeck, S. du Toit, & D.
#'   Sörbom (Eds.), Structural equation modeling: Present and future (pp.
#'   195-216). Lincolnwood, IL: Scientific Software International.
#' @source Bonifay, W. E., Reise, S. P., Scheines, R., & Meijer, R. R. (2015). When
#'   are multidimensional data unidimensional enough for structural equation
#'   modeling? An evaluation of the DETECT multidimensionality index. Structural
#'   Equation Modeling, 22, 504-516.
#' @source Reise, S. P., Scheines, R., Widaman, K. F., & Haviland, M. G. (2013).
#'   Multidimensionality and structural coefficient bias in structural equation
#'   modeling: A bifactor perspective. Educational and Psychological Measurement,
#'   73, 5-26.
#' @source Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016a). Applying bifactor
#'   statistical indices in the evaluation of psychological measures. Journal of
#'   Personality Assessment, 98, 223-237.
#' @source Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016b). Evaluating
#'   bifactor models: Calculating and interpreting statistical indices.
#'   Psychological Methods, 21, 137-150.
#'
#' @family reliability coefficients
#' @seealso [OMEGA()], the superseded function that returns these same
#'   coefficients in a wide, per-factor layout.
#'
#' @export
#'
#' @examples
#' ## From an oblique EFA (correlated-factors) solution. With no factor_map, each
#' ## item is auto-assigned to its highest-loading factor.
#' efa_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                    method = "PAF", rotation = "promax")
#' efa_reliability(efa_mod)
#'
#' ## From a Schmid-Leiman solution, with an explicit indicator-to-factor map.
#' sl_mod <- efa_schmid_leiman(efa_mod, method = "PAF")
#' fc <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2
#' efa_reliability(sl_mod, factor_map = fc)
#'
#' ## Request a subset of the coefficients only.
#' efa_reliability(sl_mod, factor_map = fc,
#'                 coefficients = c("omega_total", "alpha"))
#'
#' ## From a lavaan bifactor solution.
#' \donttest{
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#' mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
#'         F2 =~ V7 + V8 + V9 + V10 + V11 + V12
#'         F3 =~ V13 + V14 + V15 + V16 + V17 + V18
#'         g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
#'              V13 + V14 + V15 + V16 + V17 + V18'
#' fit <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
#'                    sample.nobs = 500, estimator = "ml", orthogonal = TRUE)
#' efa_reliability(fit, g_name = "g")
#' }
#' }
#'
efa_reliability <- function(model = NULL,
                            coefficients = NULL,
                            g_name = "g",
                            group_names = NULL,
                            factor_map = NULL,
                            variance = c("correlation", "sums_load"),
                            var_names = NULL,
                            fac_names = NULL,
                            g_load = NULL,
                            s_load = NULL,
                            u2 = NULL,
                            cormat = NULL,
                            pattern = NULL,
                            Phi = NULL) {

  variance <- match.arg(variance)
  checkmate::assert_string(g_name)
  checkmate::assert_character(group_names, null.ok = TRUE)

  # Validate the coefficient selector against the surfaced menu.
  menu <- .reliability_registry()$coefficient
  if (!is.null(coefficients)) {
    checkmate::assert_character(coefficients, any.missing = FALSE, min.len = 1L)
    bad <- setdiff(coefficients, menu)
    if (length(bad) > 0) {
      cli::cli_abort(
        c("Unknown {.arg coefficients}: {.val {bad}}.",
          "i" = "Choose from {.val {menu}}."),
        class = "efa_reliability_bad_coefficients"
      )
    }
  }

  # Dispatch to the right adapter, then score each factor solution through the
  # shared reliability core (add_rel = TRUE so standardized alpha is available;
  # the registry drops the core's CR/AVE columns). lavaan input is scored per
  # group with the model-implied composite variance.
  if (inherits(model, "lavaan")) {

    .require_lavaan()
    adapt <- .rel_adapt_lavaan(model, g_name = g_name, group_names = group_names)
    used_variance <- adapt$variance

    # Distinct group labels: the group column keys each block of the tidy result,
    # so duplicate labels would silently merge groups.
    if (length(adapt$groups) > 1L && anyDuplicated(adapt$group_names)) {
      cli::cli_abort(
        c("{.arg group_names} must be unique.",
          "i" = "Duplicate labels would merge groups in the result."),
        class = "efa_reliability_group_names"
      )
    }

    if (isTRUE(adapt$higher_order)) {
      cli::cli_inform(
        c("i" = "The specified general factor is a second-order factor; coefficients are computed on the Schmid-Leiman transformed second-order solution."),
        class = "efa_reliability_g_second_order"
      )
    }

    if (isTRUE(adapt$few_loadings)) {
      cli::cli_inform(
        c("i" = "Some variables have fewer than two loadings; did you enter a bifactor model? Provide a bifactor model, a second-order model, or a single-factor model."),
        class = "efa_reliability_few_loadings"
      )
    }

    mats <- vector("list", length(adapt$groups))
    informed_single <- FALSE
    for (i in seq_along(adapt$groups)) {
      grp <- adapt$groups[[i]]
      if (isTRUE(grp$single)) {
        if (!informed_single) {
          cli::cli_inform(
            c("i" = "The model contained a single factor; only omega total and the H index are returned."),
            class = "efa_reliability_single_factor"
          )
          informed_single <- TRUE
        }
        # A single factor has no group factors: score it directly into a
        # general-factor-only omega-total / H matrix the result builder accepts.
        sf <- .rel_single_factor(grp$g_load, grp$u2)
        mats[[i]] <- matrix(c(sf[["omega"]], sf[["H"]]), nrow = 1,
                            dimnames = list("g", c("tot", "H")))
      } else {
        mats[[i]] <- .reliability_core(grp, used_variance, add_ind = TRUE,
                                       add_rel = TRUE)
      }
    }

    # One matrix (single ungrouped fit) or a named list (multiple groups).
    x <- if (length(mats) == 1L) mats[[1]] else
      stats::setNames(mats, adapt$group_names)

  } else {

    spec <- .rel_dispatch_spec(model, factor_map = factor_map,
                               var_names = var_names, fac_names = fac_names,
                               g_load = g_load, s_load = s_load, u2 = u2,
                               cormat = cormat, pattern = pattern, Phi = Phi)

    # An all-zero general factor marks a correlated-factors (no-general-factor) solution
    # -- an oblique EFA, or a bifactor/manual input whose general column is zero.
    no_general <- isTRUE(all(spec$g_load == 0))

    # A correlated-factors solution has correlated group factors, so the model-implied
    # "sums_load" composite variance -- which sums the squared loading-column sums as if the
    # factors were orthogonal -- omits the factor correlations and understates the true
    # composite variance. Only the correlation-based total is correct, so such a solution is
    # always scored in correlation mode and always needs a correlation matrix (using
    # "sums_load" without one would silently return the wrong omega total). A solution with a
    # general factor needs the correlation matrix only in correlation mode. Either way, abort
    # with an actionable message when correlation mode is required but no correlation matrix
    # is available (the adapters supply one when they can; an SL object from a bare loading
    # matrix, or manual components without cormat/pattern/Phi, leave it NULL) rather than
    # divide the omega denominators by sum(NULL) = 0 and return non-finite values.
    if (no_general) {
      if (is.null(spec$cormat)) {
        cli::cli_abort(
          c("A correlation matrix is required for a correlated-factors solution but none is available.",
            "i" = "Supply {.arg cormat}, or {.arg pattern} and {.arg Phi}."),
          class = "efa_reliability_need_cormat"
        )
      }
      if (variance == "sums_load") {
        # The default is "correlation", so this warns only when "sums_load" was explicit.
        cli::cli_warn(
          c("{.code variance = \"sums_load\"} is not available for a correlated-factors solution.",
            "i" = "Its composite variance assumes uncorrelated factors; using {.code variance = \"correlation\"} instead."),
          class = "efa_reliability_variance_override"
        )
        variance <- "correlation"
      }
    } else if (variance == "correlation" && is.null(spec$cormat)) {
      cli::cli_abort(
        c("A correlation matrix is required for {.code variance = \"correlation\"} but none is available.",
          "i" = "Supply {.arg cormat} (or {.arg pattern} and {.arg Phi}), or use {.code variance = \"sums_load\"}."),
        class = "efa_reliability_need_cormat"
      )
    }

    used_variance <- variance
    x <- .reliability_core(spec, variance, add_ind = TRUE, add_rel = TRUE)

    # A correlated-factors solution (no general factor) does not identify the
    # general-factor decomposition, so the bifactor-specific indices are dropped: omega
    # hierarchical and ECV are structurally zero, and PUC ("percent of uncontaminated
    # correlations") presupposes a general factor for the cross-factor correlations to be
    # uncontaminated of. On the whole-scale (g) row the omega subscale and H index are
    # further artifacts of the synthetic all-zero general-factor column rather than
    # coefficients (the H of an all-zero loading vector is 0, and the whole-scale subscale
    # omega does not partition the composite without a general factor). Drop all of these
    # so the result surfaces only the well-defined coefficients -- whole-scale omega total
    # and alpha, and each group factor's congeneric omega, H, and alpha. NA-ing the cells
    # lets the result builder omit them as any other undefined coefficient (the
    # group-factor omega subscale and H stay on their own rows).
    if (no_general) {
      x[, c("hier", "ECV", "PUC")] <- NA
      x["g", c("sub", "H")] <- NA
    }

  }

  res <- .reliability_result(x, settings = list(variance = used_variance))

  if (!is.null(coefficients)) {
    # Note any requested coefficient this solution does not define (e.g. omega
    # hierarchical for a correlated-factors EFA, or anything but omega total and H for a
    # single-factor lavaan fit), so a reduced or empty selection is not silent.
    undefined <- setdiff(coefficients, unique(res$coefficient))
    if (length(undefined) > 0) {
      cli::cli_inform(
        c("i" = "{cli::qty(undefined)}The coefficient{?s} {.val {undefined}} {?is/are} not defined for this solution and {?is/are} omitted."),
        class = "efa_reliability_coef_undefined"
      )
    }
    res <- .reliability_select(res, coefficients)
  }

  res

}


# Dispatch a non-lavaan input to its front-end adapter, returning the normalized
# reliability spec. Mirrors OMEGA()'s model dispatch, adding the oblique EFA and
# raw bifactor-matrix paths that OMEGA does not expose. The internal adapters take a
# `type` that gates how the item-to-factor map is derived when `factor_map` is
# absent; efa_reliability always auto-maps to the highest-loading factor (adapter
# type "psych"), so a supplied map is used and an omitted one is derived.
.rel_dispatch_spec <- function(model, factor_map, var_names, fac_names,
                               g_load, s_load, u2, cormat, pattern, Phi) {

  if (inherits(model, "SL")) {
    .rel_adapt_SL(model, factor_corres = factor_map, type = "psych",
                  cormat = cormat, fac_names = fac_names)
  } else if (inherits(model, "schmid")) {
    .rel_adapt_schmid(model, factor_corres = factor_map, type = "psych",
                      cormat = cormat, fac_names = fac_names)
  } else if (inherits(model, "EFA")) {
    .rel_adapt_efa(model, factor_corres = factor_map, type = "psych",
                   cormat = cormat, fac_names = fac_names)
  } else if (is.matrix(model)) {
    if (!is.numeric(model)) {
      cli::cli_abort(
        c("A matrix {.arg model} must be a numeric bifactor loading matrix (general factor in the first column).",
          "i" = "Enter numeric loadings, or pass a fitted model object instead."),
        class = "efa_reliability_bad_matrix"
      )
    }
    .rel_adapt_bifactor(model, factor_corres = factor_map, u2 = u2,
                        cormat = cormat, fac_names = fac_names)
  } else if (is.null(model)) {
    if (is.null(var_names) || is.null(g_load) || is.null(s_load) || is.null(u2)) {
      cli::cli_abort(
        "Specify all of {.arg var_names}, {.arg g_load}, {.arg s_load}, and {.arg u2}.",
        class = "efa_reliability_missing_args"
      )
    }
    # Validate the manually supplied components, as OMEGA() does: a non-matrix
    # s_load would otherwise be coerced (e.g. a vector to a one-column matrix) and
    # silently yield wrong coefficients instead of an error.
    if (!inherits(s_load, c("matrix", "SLLOADINGS"))) {
      cli::cli_abort(
        c("Invalid {.arg s_load}.",
          "i" = "Supply a Schmid-Leiman or bifactor group-loading matrix of class {.cls matrix} or {.cls efa_sl_loadings}."),
        class = "efa_reliability_bad_s_load"
      )
    }
    checkmate::assert_numeric(g_load)
    checkmate::assert_numeric(u2)
    # The components must agree in length so they combine into one item-by-factor spec;
    # otherwise they only fail deep in the algebra with an opaque base-R error.
    n_items <- nrow(as.matrix(s_load))
    if (length(g_load) != n_items || length(u2) != n_items ||
        length(var_names) != n_items) {
      cli::cli_abort(
        c("{.arg g_load}, {.arg u2}, and {.arg var_names} must each have one entry per variable.",
          "i" = "{.arg s_load} has {n_items} row{?s}, but {.arg g_load} has {length(g_load)}, {.arg u2} has {length(u2)}, and {.arg var_names} has {length(var_names)}."),
        class = "efa_reliability_length_mismatch"
      )
    }
    .rel_adapt_manual(g_load = g_load, s_load = s_load, u2 = u2,
                      var_names = var_names, factor_corres = factor_map,
                      type = "psych", cormat = cormat, pattern = pattern, Phi = Phi,
                      fac_names = fac_names)
  } else {
    cli::cli_abort(
      c("Invalid {.arg model}.",
        "i" = "Enter a {.cls lavaan}, {.fn psych::schmid}, {.cls SL}, or {.cls EFA} object, a bifactor loading matrix, or specify {.arg var_names}, {.arg g_load}, {.arg s_load}, and {.arg u2}."),
      class = "efa_reliability_invalid_model"
    )
  }

}


# Keep only the requested coefficients in an efa_reliability result, preserving the
# object's schema, settings, and kind tagging (base data.frame subsetting drops the
# custom attributes).
.reliability_select <- function(res, coefficients) {
  settings <- attr(res, "settings")
  kind <- attr(res, "kind")
  out <- res[res$coefficient %in% coefficients, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "settings") <- settings
  attr(out, "kind") <- kind[names(kind) %in% coefficients]
  class(out) <- c("efa_reliability", "data.frame")
  out
}
