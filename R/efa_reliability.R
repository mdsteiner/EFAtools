#' Reliability coefficients for a factor solution
#'
#' Computes model-based reliability coefficients for a factor solution: McDonald's omega
#' (total, hierarchical, and subscale), standardized Cronbach's alpha, and the H index. For a
#' solution with a general factor, it also computes the common-variance indices ECV and PUC.
#' The result is a long-format table with one row per coefficient.
#'
#' @details ## Coefficients
#'
#' Each row of the result describes one composite, that is, the unit-weighted sum of a set of
#' variables: the whole scale, or the variables assigned to one group factor. Each omega is the
#' share of the composite's variance that the model gives to the factors (McDonald, 1978,
#' 1985, 1999; Zinbarg et al., 2005, 2006):
#'
#' - **Omega total** (`omega_total`): the share due to all factors. On a group-factor row, it
#'   is the omega S of Rodriguez et al. (2016a, 2016b).
#' - **Omega hierarchical** (`omega_hierarchical`): the share due to the general factor. On the
#'   whole-scale row, it tells whether the total score measures one construct.
#' - **Omega subscale** (`omega_subscale`): the share due to the group factors; on a
#'   group-factor row, the share due to its own group factor. There, it is the omega HS of
#'   Rodriguez et al. (2016a, 2016b): what the subscale score measures in addition to the
#'   general factor.
#'
#' On the whole-scale row, omega total is omega hierarchical plus omega subscale. On a
#' group-factor row, omega total can be larger, because cross-loadings on other group factors
#' add true-score variance.
#'
#' The other coefficients are:
#'
#' - **Alpha** (`alpha`): standardized Cronbach's alpha (Cronbach, 1951). It assumes
#'   essentially tau-equivalent variables and is a lower bound otherwise. For a
#'   multidimensional scale, it is rarely the coefficient to report (Gignac, 2014).
#' - **H index** (`H`): the reliability of an optimally weighted composite (Hancock & Mueller,
#'   2001). A low value shows that the variables do not define the factor well.
#' - **ECV** and **PUC** (general factor only): the share of the common variance due to the
#'   general factor, and the proportion of correlations between variables of different group
#'   factors (Bonifay et al., 2015; Reise et al., 2013; Rodriguez et al., 2016a, 2016b).
#'   Together, they tell whether a unidimensional model is acceptable.
#'
#' A solution without a general factor (a correlated-factors solution) gives omega total,
#' alpha, and the H index of each factor, and omega total and alpha of the whole scale, whose
#' row is labelled `"total"`. Omega subscale is not given: it is the share in addition to a
#' general factor, and when the factors correlate, it is not a share of the variance. A
#' single-factor solution gives omega total, alpha, and the H index. Use
#' `semTools::compRelSEM()` for composite reliability and `semTools::AVE()` for the average
#' variance extracted of a `lavaan` fit.
#'
#' ## Total variance
#'
#' With `variance = "correlation"`, the variance of a composite is the sum of its block of the
#' correlation matrix, so residual correlations that the model does not reproduce change
#' omega total. With `"sums_load"`, it is the model-implied variance. `lavaan` input always uses the
#' model-implied variance, which includes freed residual covariances. [psych::omega()] gives
#' other values: it computes the whole-scale omega total from the residual variances, and it
#' leaves cross-loadings out of the whole-scale omega subscale and of the omega total of a
#' subscale.
#'
#' ## Reverse-keyed and ordinal variables
#'
#' The coefficients describe the raw sum of the variables as fitted. Variables keyed in the
#' opposite direction decrease the true-score variance of that sum, and the function warns
#' when it finds such variables. Reverse-code them before you fit the solution (Flora, 2020).
#'
#' With polychoric or tetrachoric correlations, or a `lavaan` fit with `ordered` variables, the
#' coefficients describe the sum of the latent continuous responses, not the sum of the ordinal
#' answers. The two can be very different when the answers have few categories or are skewed. Green
#' and Yang (2009) give an omega for the ordinal sum score; this function does not compute it.
#'
#' ## Input
#'
#' - **An [efa_fit()] solution.** An oblique solution is a correlated-factors solution. An
#'   orthogonal solution is scored in the same way, with factor correlations of 0. A
#'   `"bifactorT"` solution is a bifactor solution, whose general factor is the factor that
#'   loads on every variable (the column whose smallest absolute loading is largest). A
#'   `"bifactorQ"` solution, whose group factors are correlated, and an unrotated solution with
#'   two or more factors are refused.
#' - **An [efa_schmid_leiman()] or [psych::schmid()] solution**, or the loading table of an
#'   [efa_schmid_leiman()] solution: a bifactor solution.
#' - **A loading matrix.** Without `Phi`, the first column is the general factor and the other
#'   columns are the group factors. With `Phi`, the matrix is the pattern matrix of a
#'   correlated-factors solution (not the structure matrix). Give the variable names as row
#'   names. A pattern matrix needs `cormat` for `variance = "correlation"`.
#' - **A `lavaan` fit**: single-factor, correlated-factors (cross-loadings are permitted),
#'   second-order, or bifactor. A fit is a bifactor solution if variables load on a factor named
#'   `g_name` and on other factors. Fit it with `orthogonal = TRUE`,
#'   because the coefficients need uncorrelated factors. A second-order solution is
#'   Schmid-Leiman transformed. A multiple-group fit is scored per group.
#'
#' Without `factor_map`, each variable is assigned to the group factor with its largest
#' absolute loading, and a variable without a group loading to none. A loading matrix or a
#' Schmid-Leiman solution in which most variables have an exact zero group loading (a target
#' matrix, or a transformed `lavaan` solution) uses its non-zero pattern instead. So does a
#' `lavaan` fit, but in an ESEM block, each variable is assigned to the block's factor with its
#' largest loading.
#'
#' @param model An [efa_fit()], [efa_schmid_leiman()], [psych::schmid()], or `lavaan`
#'   solution, or a loading matrix (see Input).
#' @param coefficients character. A subset of the coefficients to return: any of
#'   `"omega_total"`, `"omega_hierarchical"`, `"omega_subscale"`, `"alpha"`, `"H"`, `"ECV"`,
#'   and `"PUC"`. The default `NULL` returns all that the solution defines.
#' @param g_name character. The name of the general factor of a `lavaan` bifactor or
#'   second-order fit. Default is `"g"`.
#' @param group_names character. Labels for the groups of a `lavaan` multiple-group fit, one
#'   per group.
#' @param factor_map matrix. A logical or 0/1 matrix, variables by group factors, that assigns
#'   the variables to the group factors. A variable can belong to more than one factor. The
#'   columns must be in the column order of the solution. If `NULL` (default), see Input. Not
#'   used for `lavaan` input.
#' @param variance character. The variance of each composite: `"correlation"` (default) takes
#'   it from the correlation matrix; `"sums_load"` uses the model-implied variance and needs no
#'   correlation matrix. Not used for `lavaan` input (see Details).
#' @param fac_names character. Labels for the group factors, in the column order of the
#'   loadings; for a single-factor solution, the label of its factor. Not used for `lavaan`
#'   input.
#' @param u2 numeric. Uniquenesses, one per variable, to replace the ones derived from a loading
#'   matrix in `model`.
#' @param cormat matrix or data.frame. The correlation matrix of the variables, on the metric
#'   of the loadings, for `variance = "correlation"`. Named variables are matched by name, else
#'   by position. If `NULL`, it is taken from the solution where available.
#' @param Phi matrix. The factor intercorrelations of a pattern matrix in `model`. Not used for
#'   a fitted solution, and refused for the loading table of an [efa_schmid_leiman()] solution.
#' @param var_names,g_load,s_load,pattern `r lifecycle::badge("deprecated")` Give the loadings
#'   as one matrix in `model` instead.
#'
#' @returns An object of class `efa_reliability`: a data frame with one row per coefficient
#'   and the columns
#'   \item{coefficient}{the coefficient name.}
#'   \item{level}{`"general"` for the general factor (or the factor of a single-factor
#'     solution), `"total"` for the whole scale of a correlated-factors solution, and
#'     `"group"` for a group factor.}
#'   \item{factor}{the factor label: `"g"` for the general factor, `"total"` for the whole
#'     scale of a correlated-factors solution.}
#'   \item{group}{the group of a `lavaan` multiple-group fit, else `NA`.}
#'   \item{value}{the coefficient value.}
#'   Coefficients that the solution does not define are omitted. The attributes are
#'   `settings` (the variance used, and whether the solution has no general factor), `kind`
#'   (reliability coefficient or common-variance index), and `factor_map` (the assignment of
#'   the variables to the group factors that was used: a logical matrix, a list of them for a
#'   multiple-group fit, or `NULL` for a single factor).
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
#' @source Gignac, G. E. (2014). On the inappropriateness of using items to
#'   calculate total scale score reliability via coefficient alpha for
#'   multidimensional scales. European Journal of Psychological Assessment, 30,
#'   130-139.
#' @source Flora, D. B. (2020). Your coefficient alpha is probably wrong, but which
#'   coefficient omega is right? A tutorial on using R to obtain better reliability
#'   estimates. Advances in Methods and Practices in Psychological Science, 3,
#'   484-501.
#' @source Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
#'   structural equation modeling: An alternative to coefficient alpha. Psychometrika,
#'   74, 155-167.
#' @source Zinbarg, R. E., Revelle, W., Yovel, I., & Li, W. (2005). Cronbach's
#'   alpha, Revelle's beta, and McDonald's omega H: Their relations with each
#'   other and two alternative conceptualizations of reliability. Psychometrika,
#'   70, 123-133.
#' @source Zinbarg, R. E., Yovel, I., Revelle, W., & McDonald, R. P. (2006).
#'   Estimating generalizability to a latent variable common to all of a scale's
#'   indicators: A comparison of estimators for omega H. Applied Psychological
#'   Measurement, 30, 121-144.
#' @source Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
#'   reliability within latent variable systems. In R. Cudeck, S. du Toit, & D.
#'   Sörbom (Eds.), Structural equation modeling: Present and future - A
#'   Festschrift in honor of Karl Jöreskog (pp. 195-216). Lincolnwood, IL:
#'   Scientific Software International.
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
#' @seealso [efa_fit()] for the solution these are computed from, and [OMEGA()], the
#'   superseded function that returns these same coefficients in a wide, per-factor layout.
#'
#' @export
#'
#' @examples
#' ## From an oblique EFA (correlated-factors) solution. Each variable is assigned to the
#' ## factor with its largest loading.
#' efa_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                    estimator = "PAF", rotation = "promax")
#' efa_reliability(efa_mod)
#'
#' ## From a Schmid-Leiman solution, with an explicit variable-to-factor map.
#' sl_mod <- efa_schmid_leiman(efa_mod, estimator = "PAF")
#' fc <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2
#' efa_reliability(sl_mod, factor_map = fc)
#'
#' ## Request a subset of the coefficients only.
#' efa_reliability(sl_mod, factor_map = fc,
#'                 coefficients = c("omega_total", "alpha"))
#'
#' ## From an oblique pattern matrix and its factor intercorrelations. This is
#' ## the same correlated-factors solution, and gives the same coefficients.
#' efa_reliability(efa_mod$rot_loadings, Phi = efa_mod$Phi,
#'                 cormat = test_models$baseline$cormat)
#'
#' ## From lavaan fits: a bifactor solution, and a correlated-factors one.
#' \donttest{
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#' mod_cf <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
#'            F2 =~ V7 + V8 + V9 + V10 + V11 + V12
#'            F3 =~ V13 + V14 + V15 + V16 + V17 + V18'
#' mod <- paste(mod_cf, 'g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
#'                            V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18',
#'              sep = "\n")
#' fit <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
#'                    sample.nobs = 500, estimator = "ml", orthogonal = TRUE)
#' efa_reliability(fit, g_name = "g")
#'
#' ## No general factor: omega hierarchical, omega subscale, ECV, and PUC are omitted.
#' fit_cf <- lavaan::cfa(mod_cf, sample.cov = test_models$baseline$cormat,
#'                       sample.nobs = 500, estimator = "ml")
#' efa_reliability(fit_cf)
#' }
#' }
#'
efa_reliability <- function(model = NULL,
                            coefficients = NULL,
                            g_name = "g",
                            group_names = NULL,
                            factor_map = NULL,
                            variance = c("correlation", "sums_load"),
                            var_names = lifecycle::deprecated(),
                            fac_names = NULL,
                            g_load = lifecycle::deprecated(),
                            s_load = lifecycle::deprecated(),
                            u2 = NULL,
                            cormat = NULL,
                            pattern = lifecycle::deprecated(),
                            Phi = NULL) {

  variance <- .match_arg_ci(variance)
  checkmate::assert_string(g_name)
  checkmate::assert_character(group_names, null.ok = TRUE)

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

  # The component arguments are deprecated; from here on, an absent one is NULL.
  components <- c(var_names = lifecycle::is_present(var_names),
                  g_load = lifecycle::is_present(g_load),
                  s_load = lifecycle::is_present(s_load),
                  pattern = lifecycle::is_present(pattern))
  if (!components[["var_names"]]) var_names <- NULL
  if (!components[["g_load"]]) g_load <- NULL
  if (!components[["s_load"]]) s_load <- NULL
  if (!components[["pattern"]]) pattern <- NULL

  if (is.null(model)) {
    if (!any(components)) {
      cli::cli_abort(
        c("Supply a factor solution in {.arg model}.",
          "i" = "For example, an {.fn efa_fit}, {.fn efa_schmid_leiman}, or {.cls lavaan} solution, or a loading matrix."),
        class = "efa_reliability_missing_args"
      )
    }
    lifecycle::deprecate_warn(
      when = "1.2.0",
      what = I("Supplying the components in `g_load`, `s_load`, `var_names`, and `pattern`"),
      details = c(
        "i" = paste("Give the loadings as one matrix in `model` instead: `cbind(g_load, s_load)`",
                    "for a bifactor solution, or the pattern matrix together with `Phi` for a",
                    "correlated-factors solution. Give the variable names as its row names."),
        "i" = "Instead of `pattern`, supply `cormat = psych::factor.model(pattern, Phi, U2 = FALSE)`."
      )
    )
  }

  lavaan_fit <- inherits(model, "lavaan")
  if (lavaan_fit) .require_lavaan()

  # Report arguments that the input does not read, so that they are not dropped silently.
  if (!is.null(model)) {
    unused <- c("g_load", "s_load", "var_names", "pattern")
    if (!is.matrix(model)) unused <- c(unused, "Phi", "u2")
    if (lavaan_fit) unused <- c(unused, "cormat", "factor_map", "fac_names")
    given <- unused[!vapply(mget(unused), is.null, logical(1))]
    if (length(given) > 0) {
      cli::cli_warn(
        c("{.arg {given}} {?is/are} not used for this input.",
          "i" = if (lavaan_fit) {
            "A {.cls lavaan} fit is scored from its own estimates, and its factor labels come from the model syntax."
          } else {
            "{.arg model} carries its own loadings, uniquenesses, and factor correlations."
          }),
        class = "efa_reliability_args_ignored"
      )
    }
  }

  if (!is.null(group_names)) {
    multigroup <- lavaan_fit && lavaan::lavInspect(model, "ngroups") > 1L
    if (!multigroup) {
      cli::cli_warn(
        c("{.arg group_names} is only used for a {.cls lavaan} multiple-group fit.",
          "i" = "Every other input is scored as a single ungrouped solution, whose {.field group} column is {.code NA}."),
        class = "efa_reliability_arg_not_applicable"
      )
    }
  }

  if (lavaan_fit) {

    adapt <- .rel_adapt_lavaan(model, g_name = g_name, group_names = group_names)

    # The group column keys each block of the result, so labels must be unique.
    if (length(adapt$groups) > 1L && anyDuplicated(adapt$group_names)) {
      cli::cli_abort(
        c("{.arg group_names} must be unique.",
          "i" = "Duplicate labels would merge groups in the result."),
        class = "efa_reliability_group_names"
      )
    }

    specs <- adapt$groups
    no_general <- isTRUE(adapt$correlated)
    used_variance <- adapt$variance

  } else {

    spec <- .rel_dispatch_spec(model, factor_map = factor_map,
                               var_names = var_names, fac_names = fac_names,
                               g_load = g_load, s_load = s_load, u2 = u2,
                               cormat = cormat, pattern = pattern, Phi = Phi)

    # A solution with one factor is scored as such, whichever slot the factor came in.
    single_spec <- .rel_single_factor_spec(spec)
    single <- !is.null(single_spec)

    n_fac <- if (single) 1L else ncol(as.matrix(spec$s_load))
    if (!is.null(fac_names) && length(fac_names) != n_fac) {
      cli::cli_abort(
        c("{.arg fac_names} must have one name per factor.",
          "x" = "It has {length(fac_names)} name{?s}, but the solution has {n_fac} factor{?s}.",
          "i" = if (single) {
            "A single-factor solution has no group factors; {.arg fac_names} labels its one factor."
          } else {
            "Name the group factors in the column order of the loadings; the general factor is labelled {.val g}."
          }),
        class = "efa_reliability_fac_names_length"
      )
    }

    if (single) spec <- c(list(single = TRUE, fac_label = spec$fac_names), single_spec)
    no_general <- !single && isTRUE(all(spec$g_load == 0))

    if (variance == "correlation" && is.null(spec$cormat)) {
      cli::cli_abort(
        c("A correlation matrix is required for {.code variance = \"correlation\"} but none is available.",
          "i" = "Supply {.arg cormat}, or use {.code variance = \"sums_load\"}."),
        class = "efa_reliability_need_cormat"
      )
    }

    specs <- list(spec)
    used_variance <- variance

  }

  mats <- lapply(specs, function(s) {
    x <- .reliability_core(s, used_variance, add_ind = TRUE, add_rel = TRUE,
                           arg = "factor_map")
    if (isTRUE(s$single)) {
      .rel_drop_single_factor(x, fac_names = s$fac_label)
    } else if (no_general) {
      .rel_drop_general(x)
    } else {
      x
    }
  })

  has_map <- !all(vapply(specs, function(s) isTRUE(s$single), logical(1)))
  maps <- lapply(specs, function(s) {
    if (isTRUE(s$single)) return(NULL)
    m <- as.matrix(s$map) == 1
    dimnames(m) <- list(s$var_names, s$fac_names)
    m
  })

  if (length(mats) > 1L) {
    mats <- stats::setNames(mats, adapt$group_names)
    maps <- stats::setNames(maps, adapt$group_names)
  } else {
    mats <- mats[[1]]
    maps <- maps[[1]]
  }

  res <- .reliability_result(mats, settings = list(variance = used_variance,
                                                   no_general = no_general))

  if (!is.null(coefficients)) {
    undefined <- setdiff(coefficients, unique(res$coefficient))
    if (length(undefined) > 0) {
      cli::cli_inform(
        c("i" = "{cli::qty(undefined)}The coefficient{?s} {.val {undefined}} {?is/are} not defined for this solution and {?is/are} omitted."),
        class = "efa_reliability_coef_undefined"
      )
    }
    res <- .reliability_select(res, coefficients)
  }

  if (has_map) attr(res, "factor_map") <- maps
  res

}


# Read a loading matrix and its `Phi` as one correlated-factors solution. Shared by both matrix
# routes, so a pattern matrix gives the same coefficients with or without its loading class.
.rel_read_pattern <- function(model, Phi, u2, factor_map, cormat, fac_names) {
  .assert_args({
    checkmate::assert_numeric(model, any.missing = FALSE, finite = TRUE)
    if (!is.null(u2)) {
      checkmate::assert_numeric(u2, any.missing = FALSE, finite = TRUE,
                                len = nrow(model))
    }
  })
  L <- as.matrix(unclass(model))
  .rel_check_phi(Phi, ncol(L))
  # The communalities that the pattern implies under Phi, as on the efa_fit() route.
  if (is.null(u2)) u2 <- 1 - diag(L %*% Phi %*% t(L))
  var_names <- rownames(L)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(nrow(L)))
  if (is.null(fac_names)) {
    fac_names <- colnames(L)
    if (is.null(fac_names)) fac_names <- paste0("F", seq_len(ncol(L)))
  }
  .rel_adapt_manual(g_load = rep(0, nrow(L)), s_load = L, u2 = u2,
                    var_names = var_names, factor_corres = factor_map,
                    type = "psych", cormat = cormat, Phi = Phi,
                    fac_names = fac_names)
}


# Send a non-lavaan input to its adapter and return the normalized spec. Without `factor_map`,
# the fitted-object adapters assign each variable to its largest group loading (type "psych").
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
    # Recorded before the factor columns replace `model`, which drops the class.
    declared_hierarchy <- inherits(model, c("efa_sl_loadings", "SLLOADINGS"))
    if (declared_hierarchy) {
      # A Schmid-Leiman table is a bifactor matrix followed by its h2 and u2 columns. Drop
      # those two (by name, else by position) and take the uniquenesses from the table.
      sl <- unclass(model)
      nms <- colnames(sl)
      if (is.null(nms) && ncol(sl) > 3L) {
        display <- seq_len(ncol(sl)) > ncol(sl) - 2L
        u2_col <- ncol(sl)
      } else {
        if (is.null(nms)) nms <- character(ncol(sl))
        display <- nms %in% c("h2", "u2")
        u2_col <- match("u2", nms)
      }
      if (is.null(u2) && !is.na(u2_col)) u2 <- sl[, u2_col]
      model <- sl[, !display, drop = FALSE]
      # Mapped as the Schmid-Leiman object is, so both routes give the same coefficients.
      if (is.null(factor_map)) factor_map <- .rel_default_map(model[, -1, drop = FALSE])
    } else if (inherits(model, c("efa_loadings", "LOADINGS"))) {
      # The loading table of an ordinary solution: its first column is not a general factor.
      # With `Phi`, it is the pattern of a correlated-factors solution. A single column is one
      # factor and needs neither, so it continues to the bifactor route below.
      if (!is.null(Phi)) {
        return(.rel_read_pattern(model, Phi = Phi, u2 = u2, factor_map = factor_map,
                                 cormat = cormat, fac_names = fac_names))
      }
      if (ncol(model) > 1L) {
        cli::cli_abort(
          c("{.arg model} is a loading matrix of an ordinary factor solution, not a bifactor loading matrix.",
            "i" = "Pass the {.fn efa_fit} object itself for a correlated-factors solution, or, if this is its pattern matrix, supply the matching factor intercorrelations in {.arg Phi} to score it as one.",
            "i" = "To read it as a hierarchy instead, transform it with {.fn efa_schmid_leiman} first.",
            "i" = "If it really is a bifactor solution, pass {.code unclass(x)} to read its first column as the general factor."),
          class = "efa_reliability_pattern_matrix"
        )
      }
    } else if (!anyNA(model) && .is_cormat(model)) {
      # A correlation matrix read as loadings gives plausible but meaningless coefficients.
      # The anyNA() guard lets the loading check below report a missing value.
      cli::cli_abort(
        c("{.arg model} is a correlation matrix, not a bifactor loading matrix.",
          "i" = "Fit a solution first, e.g. with {.fn efa_fit} or {.fn efa_schmid_leiman}, and pass that.",
          "i" = "To score a solution against these correlations, pass the matrix as {.arg cormat}."),
        class = "efa_reliability_cormat_as_model"
      )
    }
    # Without `Phi`, the matrix is a hierarchy. With `Phi`, a Schmid-Leiman table contradicts
    # itself and is refused; a matrix without a loading class (which `[` and rbind() drop) is
    # read as a pattern matrix, and a warning states this reading. One column is one factor, for
    # which `Phi` changes nothing.
    if (!is.null(Phi)) {
      if (ncol(model) > 1L && declared_hierarchy) {
        cli::cli_abort(
          c("{.arg Phi} describes correlated group factors, which the loading table of a Schmid-Leiman solution does not have.",
            "i" = "To score the pattern matrix of an oblique solution, pass the {.fn efa_fit} object, or the pattern matrix in {.arg model} with this {.arg Phi}.",
            "i" = "Otherwise drop {.arg Phi}: the group factors of a hierarchy are uncorrelated."),
          class = "efa_reliability_phi_with_general"
        )
      }
      if (ncol(model) > 1L) {
        # Read before the warning, so a reading that fails is reported on its own.
        spec <- .rel_read_pattern(model, Phi = Phi, u2 = u2, factor_map = factor_map,
                                  cormat = cormat, fac_names = fac_names)
        cli::cli_warn(
          c("{.arg model} is read as the pattern matrix of a correlated-factors solution, because {.arg Phi} gives its factor correlations.",
            "i" = "Every column is a group factor, so omega hierarchical, omega subscale, ECV, and PUC are not reported.",
            "i" = "Drop {.arg Phi} to read the first column as a general factor and score the matrix as a bifactor one."),
          class = "efa_reliability_phi_as_pattern"
        )
        return(spec)
      }
      .rel_check_phi(Phi, ncol(model))
    }
    # A missing or infinite loading or uniqueness would give an ordinary-looking table.
    .assert_args({
      checkmate::assert_numeric(model, any.missing = FALSE, finite = TRUE)
      if (!is.null(u2)) {
        checkmate::assert_numeric(u2, any.missing = FALSE, finite = TRUE,
                                  len = nrow(model))
      }
    })

    .rel_adapt_bifactor(model, factor_corres = factor_map, u2 = u2,
                        cormat = cormat, fac_names = fac_names)
  } else if (is.null(model)) {
    if (is.null(var_names) || is.null(g_load) || is.null(s_load) || is.null(u2)) {
      cli::cli_abort(
        "Specify all of {.arg var_names}, {.arg g_load}, {.arg s_load}, and {.arg u2}.",
        class = "efa_reliability_missing_args"
      )
    }
    if (!inherits(s_load, c("matrix", "SLLOADINGS"))) {
      cli::cli_abort(
        c("Invalid {.arg s_load}.",
          "i" = "Supply a Schmid-Leiman or bifactor group-loading matrix of class {.cls matrix} or {.cls efa_sl_loadings}."),
        class = "efa_reliability_bad_s_load"
      )
    }
    # A missing or infinite component would give an ordinary-looking table, which no later
    # check can detect.
    .assert_args({
      checkmate::assert_numeric(s_load, any.missing = FALSE, finite = TRUE)
      checkmate::assert_numeric(g_load, any.missing = FALSE, finite = TRUE)
      checkmate::assert_numeric(u2, any.missing = FALSE, finite = TRUE)
    })
    if (.is_cormat(s_load)) {
      cli::cli_abort(
        c("{.arg s_load} is a correlation matrix, not a group-factor loading matrix.",
          "i" = "Supply the group-factor loadings of a Schmid-Leiman or bifactor solution.",
          "i" = "To score the components against these correlations, pass the matrix as {.arg cormat}."),
        class = "efa_reliability_cormat_as_s_load"
      )
    }
    # The h2 and u2 columns of a whole Schmid-Leiman table are not group factors.
    display_cols <- intersect(c("h2", "u2"), colnames(s_load))
    if (length(display_cols) > 0) {
      cli::cli_abort(
        c("{.arg s_load} carries the {.val {display_cols}} display column{?s} of a Schmid-Leiman loading table, which {?is/are} not group factors.",
          "i" = "Pass the group-factor columns only, e.g. {.code s_load = x$sl[, c(\"F1\", \"F2\")]}.",
          "i" = "To score the whole table as one solution, pass it as {.arg model} instead."),
        class = "efa_reliability_s_load_display_cols"
      )
    }
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
        "i" = "Enter an {.fn efa_fit}, {.fn efa_schmid_leiman}, {.fn psych::schmid}, or {.cls lavaan} solution, or a loading matrix.",
        "i" = "From an {.cls efa_average} object, use the averaged loading matrix in
               {.code $loadings$average}, with {.code $Phi$average} as {.arg Phi} and
               {.code $orig_R} as {.arg cormat}."),
      class = "efa_reliability_invalid_model"
    )
  }

}


# Keep only the requested coefficients, with the attributes that data frame subsetting drops.
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
