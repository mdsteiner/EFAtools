# Superseded functions
# =====================
#
# Home for the UPPERCASE names that have been superseded by lowercase `efa_*`
# equivalents. Each old name is a thin, silent bridge to its `efa_*` implementation:
# no message, no warning, no runtime signal. Most simply forward their arguments; the
# ones whose dots reach a fit repack the old flat tuning knobs into the control objects
# first (see the discipline below). Either way scripts that predate the rename keep
# working. The one deliberate exception is an ABBREVIATED tuning-knob name in a forwarded
# `...` (e.g. `max_it` for `max_iter`, `typ` for `type`): the old wrappers reached EFA()
# through do.call(), where R partial-matched such a name against the formals ahead of
# `...`, and it now aborts with `efa_unused_dots` instead. Silently running a different
# fit is worse than a loud error, so .drop_unknown_frozen_dots() lets it travel on to
# fail (see the comment there).
#
# SUPERSEDED DISCIPLINE (freeze the CONTRACT, not the implementation) -- every
# later wrapper added here MUST follow it, so the `superseded` badge stays honest:
#
#   * Never extend an old name's signature; keep each old argument list frozen
#     forever. New capabilities arrive only as new `efa_*` arguments (or
#     `estimate_control()` / `rotate_control()` knobs), which the frozen list
#     cannot name. Recorded exception: EFA() itself carries the rename of
#     `P_type`/`randomStarts` to `p_type`/`random_starts`, with silent
#     lifecycle::deprecated() catchers keeping every named legacy call working
#     unchanged. The wrappers that forward `...` (PARALLEL, KGC, HULL, SCREE,
#     NEST, N_FACTORS, SL, EFA_POOLED) are the exception: a new named `efa_*` formal IS reachable
#     through their dots. When adding one, decide deliberately whether the old
#     name should reach it, and intercept it here if it should not.
#   * A `...` that reaches a fit must be repacked with `.repack_flat_dots()` before it is
#     forwarded. The old flat interface let a caller tune the fit through these dots
#     (`type`, `max_iter`, `k`, ... -- "additional arguments passed to EFA()"), but those
#     knobs now live in the control objects and no fitting function has a formal for any of
#     them, so a bare knob forwarded on would be dropped and the fit would quietly run the
#     default preset. Repacking keeps the old call tuning the fit exactly as it used to.
#     PARALLEL, KGC, HULL, SCREE, NEST, SL, N_FACTORS and EFA_POOLED therefore translate rather
#     than forward verbatim. SL also repacks its frozen `type` argument, which the successor no
#     longer takes as a bare formal.
#   * A repacking wrapper must carry its OWN `@param ...`; never let `@inheritParams` supply
#     it. The successors tell their readers that the estimation knobs do not go in `...` --
#     true for them, and the exact opposite of what a repacking wrapper does with a legacy
#     flat call. Inheriting it documents the wrapper as breaking the code it exists to keep
#     working. Check this whenever a wrapper's dots reach a fit.
#   * Critical bug fixes in shared internals propagate through the wrapper -- this
#     is allowed for superseded functions and is why the wrapper forwards rather
#     than copies. Consciously admitted under this rule: the estimation control a
#     `N_FACTORS()` call repacks now reaches the NEST and SMT criteria, which silently dropped
#     it before, so a legacy `N_FACTORS(x, criteria = "NEST", type = "psych")` starts honouring
#     the knob it always claimed to pass on.
#   * For any DELIBERATE behaviour change to an existing pathway, consciously decide
#     whether to let it reach the old name or to pin the old behaviour inside the
#     wrapper, and document that decision.

# Reject the successor-only `estimator` name when it lands in a wrapper's `...`: the
# wrappers with a frozen `method` formal splice their own `estimator = method` translation
# into the do.call, so a caller's `estimator` dot would collide as an opaque base
# "matched by multiple actual arguments" error instead of naming the mistake.
.reject_estimator_dot <- function(nms, fn, successor) {
  if (!"estimator" %in% nms) return(invisible(NULL))
  cli::cli_abort(
    c("{.arg estimator} is not an argument of {.fn {fn}}.",
      "i" = "{.fn {fn}} selects the estimator with {.arg method}; {.arg estimator} is the
             {.fn {successor}} argument."),
    class = "efa_renamed_arg"
  )
}

#' Bartlett's test of sphericity
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `BARTLETT()` has been superseded by [efa_bartlett()], which is the recommended
#' interface going forward. It remains available and unchanged so existing code
#' keeps working.
#'
#' @inheritParams efa_bartlett
#'
#' @returns A list of class `c("efa_bartlett", "BARTLETT")`, identical to the value
#'   of [efa_bartlett()]; see there for the components.
#'
#' @seealso [efa_bartlett()]
#'
#' @export
BARTLETT <- function(x, N = NA, use = c("pairwise.complete.obs", "all.obs",
                                        "complete.obs", "everything",
                                        "na.or.complete"),
                     cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")) {
  efa_bartlett(x, N = N, use = use, cor_method = cor_method)
}

#' Kaiser-Meyer-Olkin criterion
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `KMO()` has been superseded by [efa_kmo()], which is the recommended interface
#' going forward. It remains available and unchanged so existing code keeps working.
#'
#' @inheritParams efa_kmo
#'
#' @returns A list of class `c("efa_kmo", "KMO")`, identical to the value of
#'   [efa_kmo()]; see there for the components.
#'
#' @seealso [efa_kmo()]
#'
#' @export
KMO <- function(x, use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                           "everything", "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")) {
  efa_kmo(x, use = use, cor_method = cor_method)
}

#' Parallel analysis
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `PARALLEL()` has been superseded by [efa_parallel()], which is the recommended
#' interface going forward. It remains available and unchanged so existing code
#' keeps working.
#'
#' @inheritParams efa_parallel
#' @param ... Further arguments passed on to the [efa_fit()] fits. For example,
#'  `estimator`, to change the estimator (default is "PAF"; PAF is more robust, but it
#'  will take longer compared to "ML" and "ULS"), or one of the estimation tuning knobs
#'  (`type`, `init_comm`, `criterion`, `criterion_type`, `max_iter`, `abs_eigen`,
#'  `start_method`), which are repacked into an [estimate_control()] object so that they
#'  tune the fits exactly as they always did.
#'
#' @returns An object of class `efa_retention`, identical to the value of
#'   [efa_parallel()]; see there for the components.
#'
#' @seealso [efa_parallel()]
#'
#' @export
PARALLEL <- function(x = NULL,
                     N = NA,
                     n_vars = NA,
                     n_datasets = 1000,
                     percent = 95,
                     eigen_type = c("PCA", "SMC", "EFA"),
                     use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                             "everything", "na.or.complete"),
                     cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                     decision_rule = c("means", "percentile", "crawford"),
                     n_factors = 1,
                     ...) {
  do.call(efa_parallel,
          c(list(x, N = N, n_vars = n_vars, n_datasets = n_datasets,
                 percent = percent, eigen_type = eigen_type, use = use,
                 cor_method = cor_method, decision_rule = decision_rule,
                 n_factors = n_factors),
            .repack_flat_dots(.drop_unknown_frozen_dots(list(...)))))
}

#' Empirical Kaiser criterion
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `EKC()` has been superseded by [efa_ekc()], which is the recommended interface
#' going forward. It remains available and unchanged so existing code keeps working.
#'
#' @inheritParams efa_ekc
#'
#' @returns An object of class `efa_retention`, identical to the value of
#'   [efa_ekc()]; see there for the components.
#'
#' @seealso [efa_ekc()]
#'
#' @export
EKC <- function(x, N = NA,
                use = c("pairwise.complete.obs", "all.obs",
                           "complete.obs", "everything",
                           "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                type = "BvA2017") {
  efa_ekc(x, N = N, use = use, cor_method = cor_method, type = type)
}

#' Kaiser-Guttman criterion
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `KGC()` has been superseded by [efa_kgc()], which is the recommended interface
#' going forward. It remains available and unchanged so existing code keeps working.
#'
#' @inheritParams efa_kgc
#' @param ... Further arguments passed on to the [efa_fit()] fit. For example,
#'  `estimator`, to change the estimator (PAF is default), or one of the estimation
#'  tuning knobs (`type`, `init_comm`, `criterion`, `criterion_type`, `max_iter`,
#'  `abs_eigen`, `start_method`), which are repacked into an [estimate_control()]
#'  object so that they tune the fit exactly as they always did.
#'
#' @returns An object of class `efa_retention`, identical to the value of
#'   [efa_kgc()]; see there for the components.
#'
#' @seealso [efa_kgc()]
#'
#' @export
KGC <- function(x, eigen_type = c("PCA", "SMC", "EFA"),
                use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                        "everything", "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                n_factors = 1, ...) {
  do.call(efa_kgc,
          c(list(x, eigen_type = eigen_type, use = use, cor_method = cor_method,
                 n_factors = n_factors),
            .repack_flat_dots(.drop_unknown_frozen_dots(list(...)))))
}

#' Hull method
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `HULL()` has been superseded by [efa_hull()], which is the recommended
#' interface going forward. It remains available and unchanged so existing code
#' keeps working.
#'
#' @inheritParams efa_hull
#' @param method character. The estimator to use; passed to [efa_hull()] as its
#'   `estimator` argument. One of `"PAF"`, `"ULS"`, or `"ML"`.
#' @param ... Further arguments passed on to the [efa_fit()] fits, including the
#'  estimation tuning knobs (`type`, `init_comm`, `criterion`, `criterion_type`,
#'  `max_iter`, `abs_eigen`, `start_method`), which are repacked into an
#'  [estimate_control()] object so that they tune the fits exactly as they always did.
#'  The estimator is selected with `method`.
#'
#' @returns An object of class `efa_retention`, identical to the value of
#'   [efa_hull()]; see there for the components.
#'
#' @seealso [efa_hull()]
#'
#' @export
HULL <- function(x, N = NA, n_fac_theor = NA,
                 method = c("PAF", "ULS", "ML"), gof = c("CAF", "CFI", "RMSEA"),
                 eigen_type = c("SMC", "PCA", "EFA"),
                 use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                         "everything", "na.or.complete"),
                 cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                 n_datasets = 1000, percent = 95,
                 decision_rule = c("means", "percentile", "crawford"),
                 n_factors = 1, ...) {
  .reject_estimator_dot(...names(), "HULL", "efa_hull")
  do.call(efa_hull,
          c(list(x, N = N, n_fac_theor = n_fac_theor, estimator = method, gof = gof,
                 eigen_type = eigen_type, use = use, cor_method = cor_method,
                 n_datasets = n_datasets, percent = percent,
                 decision_rule = decision_rule, n_factors = n_factors),
            .repack_flat_dots(.drop_unknown_frozen_dots(list(...)))))
}

#' Scree plot
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `SCREE()` has been superseded by [efa_scree()], which is the recommended
#' interface going forward. It remains available and unchanged so existing code
#' keeps working.
#'
#' @inheritParams efa_scree
#' @param ... Further arguments passed on to the [efa_fit()] fit. For example,
#'  `estimator`, to change the estimator (PAF is default), or one of the estimation
#'  tuning knobs (`type`, `init_comm`, `criterion`, `criterion_type`, `max_iter`,
#'  `abs_eigen`, `start_method`), which are repacked into an [estimate_control()]
#'  object so that they tune the fit exactly as they always did.
#'
#' @returns An object of class `efa_retention`, identical to the value of
#'   [efa_scree()]; see there for the components.
#'
#' @seealso [efa_scree()]
#'
#' @export
SCREE <- function(x, eigen_type = c("PCA", "SMC", "EFA"),
                  use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                        "everything", "na.or.complete"),
                  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                  n_factors = 1, ...) {
  do.call(efa_scree,
          c(list(x, eigen_type = eigen_type, use = use, cor_method = cor_method,
                 n_factors = n_factors),
            .repack_flat_dots(.drop_unknown_frozen_dots(list(...)))))
}

#' Minimum average partial
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `MAP()` has been superseded by [efa_map()], which is the recommended interface
#' going forward. It remains available and unchanged so existing code keeps working.
#'
#' @inheritParams efa_map
#'
#' @returns An object of class `efa_retention`, identical to the value of
#'   [efa_map()]; see there for the components.
#'
#' @seealso [efa_map()]
#'
#' @export
MAP <- function(x,
                use = c("pairwise.complete.obs", "all.obs",
                        "complete.obs", "everything",
                        "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")) {
  efa_map(x, use = use, cor_method = cor_method)
}

#' Next eigenvalue sufficiency test
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `NEST()` has been superseded by [efa_nest()], which is the recommended
#' interface going forward. It remains available and unchanged so existing code
#' keeps working.
#'
#' @inheritParams efa_nest
#' @param ... Further arguments passed on to the [efa_fit()] fits. For example,
#'  `estimator`, to change the estimator (PAF is default), or one of the estimation
#'  tuning knobs (`type`, `init_comm`, `criterion`, `criterion_type`, `max_iter`,
#'  `abs_eigen`, `start_method`), which are repacked into an [estimate_control()]
#'  object so that they tune the fits exactly as they always did.
#'
#' @returns An object of class `efa_retention`, identical to the value of
#'   [efa_nest()]; see there for the components.
#'
#' @seealso [efa_nest()]
#'
#' @export
NEST <- function(x, N = NA,
                 alpha = .05,
                 use = c("pairwise.complete.obs", "all.obs",
                         "complete.obs", "everything",
                         "na.or.complete"),
                 cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                 n_datasets = 1000,
                 ...) {
  do.call(efa_nest,
          c(list(x, N = N, alpha = alpha, use = use, cor_method = cor_method,
                 n_datasets = n_datasets),
            .repack_flat_dots(.drop_unknown_frozen_dots(list(...)))))
}

#' Sequential model tests
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `SMT()` has been superseded by [efa_smt()], which is the recommended interface
#' going forward. It remains available and unchanged so existing code keeps working.
#'
#' @inheritParams efa_smt
#'
#' @returns An object of class `efa_retention`, identical to the value of
#'   [efa_smt()]; see there for the components.
#'
#' @seealso [efa_smt()]
#'
#' @export
SMT <- function(x, N = NA, use = c("pairwise.complete.obs", "all.obs",
                                   "complete.obs", "everything",
                                   "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")) {
  efa_smt(x, N = N, use = use, cor_method = cor_method)
}

#' Comparison data
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `CD()` has been superseded by [efa_cd()], which is the recommended interface
#' going forward. It remains available and unchanged so existing code keeps working.
#'
#' @inheritParams efa_cd
#'
#' @returns An object of class `efa_retention`, identical to the value of
#'   [efa_cd()]; see there for the components.
#'
#' @seealso [efa_cd()]
#'
#' @export
CD <- function(x, n_factors_max = NA, N_pop = 10000, N_samples = 500, alpha = .30,
               cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
               max_iter = 50) {
  efa_cd(x, n_factors_max = n_factors_max, N_pop = N_pop, N_samples = N_samples,
         alpha = alpha, cor_method = cor_method, max_iter = max_iter)
}

#' Various factor retention criteria
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `N_FACTORS()` has been superseded by [efa_retain()], which is the recommended
#' interface going forward. It remains available and unchanged so existing code
#' keeps working.
#'
#' @inheritParams efa_retain
#' @param criteria character. A vector with the factor retention methods to
#'   perform. Possible inputs are: `"CD"`, `"EKC"`, `"HULL"`, `"KGC"`, `"MAP"`,
#'   `"NEST"`, `"PARALLEL"`, `"SCREE"`, and `"SMT"` (see the details in
#'   [efa_retain()]). By default, a subset of often used, well-performing
#'   methods are performed.
#' @param method character. The estimator to use in the criteria that fit EFA models;
#'   passed to [efa_retain()] as its `estimator` argument. One of `"ML"`, `"PAF"`, or
#'   `"ULS"`.
#' @param ... Further arguments passed on to the [efa_fit()] fits, including the
#'  estimation tuning knobs (`type`, `init_comm`, `criterion`, `criterion_type`,
#'  `abs_eigen`, `start_method`), which are repacked into an [estimate_control()] object
#'  so that they tune the fits exactly as they always did. The estimator is selected with
#'  `method`; `max_iter` is taken by the `max_iter_CD` argument (R matches an abbreviated
#'  name against the arguments before `...`) and so does not reach the fits.
#'
#' @returns A list of class `c("efa_retain", "N_FACTORS")`, identical to the
#'   value of [efa_retain()]; see there for the components.
#'
#' @seealso [efa_retain()]
#'
#' @export
N_FACTORS <- function(x, criteria = c("CD", "EKC", "HULL", "MAP", "NEST", "PARALLEL"),
                      suitability = TRUE, N = NA,
                      use = c("pairwise.complete.obs", "all.obs",
                              "complete.obs", "everything", "na.or.complete"),
                      cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                      n_factors_max = NA, N_pop = 10000, N_samples = 500,
                      alpha = .30, max_iter_CD = 50, n_fac_theor = NA,
                      method = c("ML", "PAF", "ULS"),
                      gof = c("CAF", "CFI", "RMSEA"),
                      eigen_type_HULL = c("SMC", "PCA", "EFA"),
                      eigen_type_other = c("SMC"),
                      n_factors = 1, n_datasets = 1000,
                      percent = 95,
                      decision_rule = c("means", "percentile", "crawford"),
                      ekc_type = c("BvA2017"),
                      n_datasets_nest = 1000, alpha_nest = .05,
                      show_progress = FALSE,
                      ...) {
  .reject_estimator_dot(...names(), "N_FACTORS", "efa_retain")
  do.call(efa_retain,
          c(list(x, criteria = criteria, suitability = suitability, N = N,
                 use = use, cor_method = cor_method, n_factors_max = n_factors_max,
                 N_pop = N_pop, N_samples = N_samples, alpha = alpha,
                 max_iter_CD = max_iter_CD, n_fac_theor = n_fac_theor,
                 estimator = method, gof = gof, eigen_type_HULL = eigen_type_HULL,
                 eigen_type_other = eigen_type_other, n_factors = n_factors,
                 n_datasets = n_datasets, percent = percent,
                 decision_rule = decision_rule, ekc_type = ekc_type,
                 n_datasets_nest = n_datasets_nest, alpha_nest = alpha_nest,
                 show_progress = show_progress),
            .repack_flat_dots(.drop_unknown_frozen_dots(list(...)))))
}

#' Schmid-Leiman transformation
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `SL()` has been superseded by [efa_schmid_leiman()], which is the recommended
#' interface going forward. It remains available and unchanged so existing code
#' keeps working.
#'
#' @inheritParams efa_schmid_leiman
#' @param type character. One of "EFAtools" (default), "psych", "SPSS", or "none". This is
#'  used to control the procedure of the second-order factor analysis. In
#'  [efa_schmid_leiman()] it is set through the `type` of the [estimate_control()] object.
#' @param method character. The estimator for the second-order factor analysis; passed to
#'  [efa_schmid_leiman()] as its `estimator` argument. One of `"PAF"`, `"ML"`, `"ULS"`, or
#'  `"MINRES"`.
#' @param ... Further arguments passed on to the second-order [efa_fit()], including the
#'  estimation tuning knobs (`init_comm`, `criterion`, `criterion_type`, `max_iter`,
#'  `abs_eigen`, `start_method`), which are repacked, together with `type`, into an
#'  [estimate_control()] object so that they tune that fit exactly as they always did.
#'  The estimator is selected with `method`.
#'
#' @returns A list of class `c("efa_schmid_leiman", "SL")`, identical to the value
#'   of [efa_schmid_leiman()]; see there for the components.
#'
#' @seealso [efa_schmid_leiman()]
#'
#' @export
SL <- function(x, Phi = NULL, type = c("EFAtools", "psych", "SPSS", "none"),
               method = c("PAF", "ML", "ULS", "MINRES"), g_name = "g", ...) {
  .reject_estimator_dot(...names(), "SL", "efa_schmid_leiman")
  # The successor-only control objects must not ride through the frozen dots: with a control
  # present, .repack_flat_dots() translates nothing, so the wrapper's own `type` argument --
  # and any flat knob alongside it -- would be silently discarded and the second-order fit
  # would quietly run the control's preset instead of the requested one. No pre-rename code
  # could have passed these names meaningfully, so rejecting them breaks nothing frozen (the
  # same rule EFA() applies).
  new_only <- intersect(...names(), c("estimate_control", "rotate_control"))
  if (length(new_only) > 0L) {
    cli::cli_abort(
      c("{.arg {new_only}} {?is/are} not {?an argument/arguments} of {.fn SL}.",
        "i" = "{.fn SL} takes the tuning knobs flat and selects the preset with
               {.arg type}; the control objects belong to {.fn efa_schmid_leiman}."),
      class = "efa_renamed_arg"
    )
  }
  # Resolve `type` once so the frozen preset argument, and any flat knob the dots carry, reach
  # the second-order fit through the control objects it now takes. The second-order fit is
  # always unrotated, so the frozen filter's retention defaults apply.
  type <- match.arg(type)
  do.call(efa_schmid_leiman,
          c(list(x, Phi = Phi, estimator = method, g_name = g_name),
            .repack_flat_dots(.drop_unknown_frozen_dots(list(...)), type = type)))
}

#' Compare two vectors or matrices (communalities or loadings)
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `COMPARE()` has been superseded by [efa_compare()], which is the recommended
#' interface going forward. It remains available and unchanged so existing code
#' keeps working.
#'
#' @inheritParams efa_compare
#'
#' @returns A list of class `c("efa_compare", "COMPARE")`, identical to the value
#'   of [efa_compare()]; see there for the components.
#'
#' @seealso [efa_compare()]
#'
#' @export
COMPARE <- function(x,
                    y,
                    reorder = c("congruence", "names", "none"),
                    corres = TRUE,
                    thresh = .3,
                    digits = 4,
                    m_red = .001,
                    range_red = .001,
                    round_red = 3,
                    print_diff = TRUE,
                    na.rm = FALSE,
                    x_labels = c("x", "y"),
                    plot = TRUE,
                    plot_red = .01) {
  efa_compare(x, y = y, reorder = reorder, corres = corres, thresh = thresh,
              digits = digits, m_red = m_red, range_red = range_red,
              round_red = round_red, print_diff = print_diff, na.rm = na.rm,
              x_labels = x_labels, plot = plot, plot_red = plot_red)
}

#' Rotate a loading matrix to a target using Procrustes alignment
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `PROCRUSTES()` has been superseded by [efa_procrustes()], which is the
#' recommended interface going forward. It remains available and unchanged so
#' existing code keeps working.
#'
#' @inheritParams efa_procrustes
#'
#' @returns A list identical to the value of [efa_procrustes()]; see there for the
#'   components.
#'
#' @seealso [efa_procrustes()]
#'
#' @export
PROCRUSTES <- function(A,
                       Target,
                       rotation = c("orthogonal", "oblique"),
                       S = NULL,
                       T_init = NULL,
                       oblique_eps = 1e-5,
                       oblique_maxit = 1000,
                       oblique_max_line_search = 10,
                       oblique_step0 = 1,
                       oblique_normalize = FALSE,
                       oblique_random_starts = 0,
                       oblique_screen_keep = 2,
                       oblique_triage_maxit = 25,
                       oblique_triage_improve_tol = 0) {
  efa_procrustes(A, Target = Target, rotation = rotation, S = S,
                 T_init = T_init, oblique_eps = oblique_eps,
                 oblique_maxit = oblique_maxit,
                 oblique_max_line_search = oblique_max_line_search,
                 oblique_step0 = oblique_step0,
                 oblique_normalize = oblique_normalize,
                 oblique_random_starts = oblique_random_starts,
                 oblique_screen_keep = oblique_screen_keep,
                 oblique_triage_maxit = oblique_triage_maxit,
                 oblique_triage_improve_tol = oblique_triage_improve_tol)
}

# `EFA_AVERAGE` and `efa_average` differ only in case, so the two topics cannot both be
# written to `man/`: file names that differ only in case are not portable. `@rdname` moves
# this wrapper's topic aside; `?EFA_AVERAGE` still resolves through the alias.

#' Model averaging across different EFA methods and types
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `EFA_AVERAGE()` has been superseded by [efa_average()], which is the recommended
#' interface going forward. It remains available and unchanged so existing code
#' keeps working.
#'
#' @inheritParams efa_average
#'
#' @param method character vector. Any combination of `"PAF"`, `"ML"`, and `"ULS"`,
#'   the estimators to average across; passed to [efa_average()] as its `estimator`
#'   argument. Default is `"PAF"`.
#' @param P_type character vector. Any combination of `"norm"` and `"unnorm"`.
#'   How the promax target matrix P is computed if `"none"` is among the specified
#'   types and `"promax"` or `"oblique"` is among the specified rotations:
#'   `"unnorm"` uses the unnormalized target matrix of Hendrickson and White
#'   (1964), `"norm"` a normalized one. This frozen argument keeps its original
#'   name; [efa_average()] takes the same setting as `p_type`. Default is
#'   `c("norm", "unnorm")`.
#'
#' @returns The value of [efa_average()], normally a list of class
#'   `c("efa_average", "EFA_AVERAGE")`; see there for the components.
#'
#' @seealso [efa_average()]
#'
#' @rdname EFA_AVERAGE-superseded
#' @export
EFA_AVERAGE <- function(x, n_factors, N = NA, method = "PAF", rotation = "promax",
                        type = "none", averaging = c("mean", "median"), trim = 0,
                        salience_threshold = .3,
                        max_iter = 1e4,
                        init_comm = c("smc", "mac", "unity"),
                        criterion = c(1e-3),
                        criterion_type = c("sum", "max_individual"),
                        abs_eigen = c(TRUE),
                        varimax_type = c("svd", "kaiser"),
                        normalize = TRUE,
                        k_promax = 2:4, k_simplimax = ncol(x),
                        P_type = c("norm", "unnorm"), precision = 1e-5,
                        start_method = c("psych", "factanal"),
                        use = c("pairwise.complete.obs", "all.obs",
                                "complete.obs", "everything", "na.or.complete"),
                        cor_method = c("pearson", "spearman", "kendall", "poly",
                                       "tetra", "fiml"),
                        show_progress = TRUE) {
  efa_average(x, n_factors = n_factors, N = N, estimator = method,
              rotation = rotation, type = type, averaging = averaging,
              trim = trim, salience_threshold = salience_threshold,
              max_iter = max_iter, init_comm = init_comm, criterion = criterion,
              criterion_type = criterion_type, abs_eigen = abs_eigen,
              varimax_type = varimax_type, normalize = normalize,
              k_promax = k_promax, k_simplimax = k_simplimax, p_type = P_type,
              precision = precision, start_method = start_method, use = use,
              cor_method = cor_method, show_progress = show_progress)
}

# `EFA()` is a translating wrapper, not a pure forwarder: `efa_fit()` collects the
# estimation and rotation tuning knobs into two control objects, so the frozen flat
# argument list is repacked into `estimate_control()` / `rotate_control()` here. This is
# why `EFA()` is exempt from the static one-line-forwarder contract the plain forwarders meet;
# its behaviour is pinned by a byte-identical `EFA()` == `efa_fit()` regression test instead.

#' Exploratory factor analysis (EFA)
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `EFA()` has been superseded by [efa_fit()], which is the recommended interface going
#' forward. `efa_fit()` keeps the primary choices (data, factors, estimator, rotation,
#' standard errors) as top-level arguments and collects the estimation and rotation tuning
#' knobs into two control objects, [estimate_control()] and [rotate_control()]. `EFA()`
#' remains available and unchanged -- its full flat argument list still works exactly as
#' before -- so existing code keeps running.
#'
#' @inheritParams efa_fit
#' @param method character. The estimator used to fit the EFA; passed to [efa_fit()]
#'  as its `estimator` argument. One of "PAF", "ML", "ULS", "MINRES" (an accepted
#'  alias of "ULS"), or "DWLS"; see the [efa_fit()] documentation for their
#'  properties and data requirements.
#' @param type character. If one of "EFAtools" (default), "psych", or "SPSS" is
#'  used, and the following arguments with default NA are left with
#'  NA, these implementations are executed according to the respective program
#'  ("psych" and "SPSS") or according to the best solution found in Grieder &
#'  Steiner (2022; "EFAtools"). Individual properties can be adapted using one of
#'  the three types and specifying some of the following arguments. If set to
#'  "none" additional arguments must be specified depending on the `method`
#'  and `rotation` used (see details).
#' @param max_iter numeric. The maximum number of iterations to perform after which
#' the iterative PAF procedure is halted with a warning. If `type` is one of
#' "EFAtools", "SPSS", or "psych", this is automatically specified if `max_iter` is
#' left to be `NA`, but can be overridden by entering a number. Default is
#' `NA`.
#' @param init_comm character. The method to estimate the initial communalities
#' in `PAF`. "smc" will use squared multiple correlations, "mac" will use
#' maximum absolute correlations, "unity" will use 1s (see details).
#' Default is `NA`.
#' @param criterion numeric. The convergence criterion used for PAF.
#' If the change in communalities from one iteration to the next is smaller than
#' this criterion the solution is accepted and the procedure ends.
#' Default is `NA`.
#' @param criterion_type character. Type of convergence criterion used for
#' PAF. "max_individual" selects the maximum change in any of the
#' communalities from one iteration to the next and tests it against the
#' specified criterion. This is also used by SPSS. "sum" takes the difference of
#' the sum of all communalities in one iteration and the sum of all communalities
#' in the next iteration and tests this against the criterion. This procedure is
#' used by the [psych::fa()] function. Default is `NA`.
#' @param abs_eigen logical. Which algorithm to use in the PAF
#' iterations. If FALSE, the loadings are computed from the eigenvalues. This is
#' also used by the [psych::fa()] function. If TRUE the
#' loadings are computed with the absolute eigenvalues as done by SPSS.
#' Default is `NA`.
#' @param k numeric. Either the power used for computing the target matrix P in
#' the promax rotation or the number of 'close to zero loadings' for the simplimax
#' rotation. If left to `NA` (default), the value for promax depends on the specified type.
#' For simplimax, `nrow(L)`, where L is the matrix of unrotated loadings,
#' is used by default.
#' @param normalize logical. If `TRUE`, a kaiser normalization is
#' performed before the specified rotation. Default is `TRUE`.
#' @param p_type character. This specifies how the target
#' matrix P is computed in promax rotation. If "unnorm" it will use the
#' unnormalized target matrix as originally done in Hendrickson and White (1964).
#' This is also used in the psych and stats packages. If "norm" it will use the
#' normalized target matrix as used in SPSS. Default is `NA`.
#' @param precision numeric. The tolerance for stopping in the rotation
#' procedure. Default is 10^-5 for all rotation methods.
#' @param varimax_type character. The type of the varimax rotation performed.
#' If "svd", singular value decomposition is used, as [stats::varimax()] does. If
#' "kaiser", the varimax procedure performed in SPSS is used, following the original
#' procedure from Kaiser (1958) (see details). Default is `NA`.
#' @param order_type character. How to order the factors. "eigen" reorders the
#' factors by descending explained variance; "ss_factors" reorders the factors by
#' descending (unweighted) sum of squared factor loadings per factor. Default is `NA`.
#' @param start_method character. How to specify the starting values for the
#' optimization procedure for ML. Default is "psych" which takes the
#' starting values specified in [psych::fa()]. "factanal" takes the
#' starting values specified in the [stats::factanal()] function.
#' @param random_starts numeric. The number of random starts to use in the
#'  rotation to guard against local minima. Default is 100.
#' @param P_type,randomStarts `r lifecycle::badge("superseded")` Former names of
#'  `p_type` and `random_starts`. Still accepted (silently) for backwards
#'  compatibility; please use the new names.
#' @param ... Additional arguments passed to the rotation procedure (e.g., `maxit` for
#'  the maximum number of iterations).
#'
#' @returns The value of [efa_fit()], a list of class `c("efa", "EFA")`; see there for the
#'   components.
#'
#' @seealso [efa_fit()], [estimate_control()], [rotate_control()]
#'
#' @export
EFA <- function(x, n_factors, N = NA, method = c("PAF", "ML", "ULS", "MINRES", "DWLS"),
                rotation = c("none", "varimax", "equamax", "quartimax", "geominT",
                             "bentlerT", "bifactorT", "promax", "oblimin",
                             "quartimin", "simplimax", "bentlerQ", "geominQ",
                             "bifactorQ"),
                se = c("none", "information", "sandwich", "np-boot"),
                type = c("EFAtools", "psych", "SPSS", "none"), max_iter = NA,
                init_comm = NA, criterion = NA, criterion_type = NA,
                abs_eigen = NA, use = c("pairwise.complete.obs", "all.obs",
                                          "complete.obs", "everything",
                                          "na.or.complete"),
                varimax_type = NA,
                k = NA, normalize = TRUE, p_type = NA, precision = 1e-5,
                order_type = NA, start_method = "psych",
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra",
                               "fiml"),
                b_boot = 1000, ci = .95,
                random_starts = 100, seed = NULL,
                P_type = lifecycle::deprecated(),
                randomStarts = lifecycle::deprecated(),
                ...) {

  # Resolve the former argument spellings before translating to the control objects.
  if (lifecycle::is_present(P_type)) p_type <- P_type
  if (lifecycle::is_present(randomStarts)) random_starts <- randomStarts

  # Resolve `type` once so the estimation and rotation controls share the single preset the
  # old flat interface always applied to both. efa_fit() match.args the remaining choice
  # arguments and runs every input guard, so nothing else is resolved here.
  type <- match.arg(type)

  # The flat interface always ignored a `...` argument its rotation engine could not consume;
  # efa_fit() rejects such an argument. Keep the frozen contract by dropping the unconsumable
  # dots -- by the exact names the engines read -- before efa_fit() sees them. With
  # rotation = "none" no engine runs and the accepted set is empty, so the dots are dropped
  # in full and the call stays silent, as it always has been.
  #
  # `rotation` is resolved here so the lookup gets a canonical name, and the frozen
  # match.arg() is deliberate: efa_fit() matches its choices case-insensitively
  # (`.match_arg_ci()`), a successor-only feature, so the wrapper keeps the case-sensitive
  # matching its callers have always had. Forwarding the canonical value below also keeps
  # the successor's wider matching from reaching the frozen argument list.
  rotation <- match.arg(rotation)
  dots <- list(...)
  # The successor-only arguments must not ride through the frozen dots: none of the three is
  # a name a rotation engine reads, so the ignore-unknown-dots contract would silently drop
  # all of them and the fit would quietly run the flat defaults -- `estimator` losing out to
  # the wrapper's own `estimator = method` translation, the control objects to the ones the
  # wrapper builds below. No pre-rename code could have passed these names meaningfully, so
  # rejecting them breaks nothing frozen.
  new_only <- intersect(names(dots), c("estimator", "estimate_control", "rotate_control"))
  if (length(new_only) > 0L) {
    cli::cli_abort(
      c("{.arg {new_only}} {?is/are} not {?an argument/arguments} of {.fn EFA}.",
        "i" = "{.fn EFA} selects the estimator with {.arg method} and takes the tuning
               knobs flat; {.arg estimator} and the control objects belong to
               {.fn efa_fit}."),
      class = "efa_renamed_arg"
    )
  }
  dots <- dots[names(dots) %in% .rotation_dot_extras(rotation)]

  args <- c(list(x = x, N = N, estimator = method, rotation = rotation,
                 se = se, cor_method = cor_method, use = use, b_boot = b_boot, ci = ci,
                 seed = seed,
                 estimate_control = estimate_control(type = type, init_comm = init_comm,
                                                     criterion = criterion,
                                                     criterion_type = criterion_type,
                                                     max_iter = max_iter, abs_eigen = abs_eigen,
                                                     start_method = start_method),
                 rotate_control = rotate_control(type = type, normalize = normalize,
                                                 precision = precision, order_type = order_type,
                                                 varimax_type = varimax_type, p_type = p_type,
                                                 k = k, random_starts = random_starts)),
            dots)
  # `n_factors` has no default; do.call() would force it here and turn a call that omitted it
  # (e.g. `EFA(1:5)`) into an eager missing-argument error before efa_fit()'s input checks
  # run, so it is added only when supplied.
  if (!missing(n_factors)) args$n_factors <- n_factors
  do.call(efa_fit, args)
}

#' Exploratory factor analysis on multiple data imputations
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `EFA_POOLED()` has been superseded by [efa_mi()], which is the recommended
#' interface going forward. It remains available and unchanged so existing code
#' keeps working.
#'
#' @inheritParams efa_mi
#'
#' @returns The value of [efa_mi()], normally a list of class
#'   `c("efa_mi", "EFA_POOLED", "efa", "EFA")`; see there for the components.
#'
#' @seealso [efa_mi()]
#'
#' @export
EFA_POOLED <- function(data_list,
                       p = 0.05,
                       target_method = c("first_target", "consensus"),
                       align_unrotated = c("signed_tucker_congruence", "none", "procrustes"),
                       fit_pool_method = c("D2"),
                       consensus_args = list(),
                       procrustes_args = list(),
                       rmsea_ci_level = .90,
                       rmsr_upper = TRUE,
                       ...) {
  # The per-imputation fits ARE rotated here, so the frozen filter keeps `rotation` and the
  # extras the selected rotation's engine reads -- the same per-rotation contract EFA()
  # froze. The rotation is read off the dots verbatim: pre-rename code always passed the
  # canonical spelling (the flat interface matched case-sensitively), and an invalid value
  # is rejected downstream by efa_mi() itself.
  dots <- list(...)
  rot <- if (is.character(dots$rotation) && length(dots$rotation) == 1L) {
    dots$rotation
  } else {
    "none"
  }
  # `rmsr_upper` stays a formal (the frozen signature) but is not forwarded: it is deprecated
  # on efa_mi(), and passing the frozen default on would fire that deprecation for every
  # legacy call. It never selected anything -- see the argument's documentation -- so dropping
  # it here leaves the result identical.
  do.call(efa_mi,
          c(list(data_list, p = p, target_method = target_method,
                 align_unrotated = align_unrotated, fit_pool_method = fit_pool_method,
                 consensus_args = consensus_args, procrustes_args = procrustes_args,
                 rmsea_ci_level = rmsea_ci_level),
            .repack_flat_dots(.drop_unknown_frozen_dots(
              dots, extras = .rotation_dot_extras(rot), unrotated = FALSE))))
}
