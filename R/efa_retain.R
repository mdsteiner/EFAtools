#' Various factor retention criteria
#'
#' Choosing the number of factors to retain is one of the most important
#' decisions in an exploratory factor analysis (EFA). Many criteria exist to
#' help with this choice. This function runs several of them together, and can
#' also check whether the data are suitable for factor analysis.
#'
#' @param x data.frame or matrix. Raw data, or a correlation matrix. `"CD"` needs raw
#'   data.
#' @param criteria character. The factor retention criteria to run: one or more of
#'   `"CD"`, `"EKC"`, `"HULL"`, `"KGC"`, `"MAP"`, `"NEST"`, `"PARALLEL"`, `"SCREE"`, and
#'   `"SMT"`. The default runs a subset of commonly used, well-performing criteria (see
#'   details).
#' @param suitability logical. Whether the data should be checked for suitability
#' for factor analysis using Bartlett's test of sphericity and the
#' Kaiser-Meyer-Olkin criterion (see details). Default is `TRUE`.
#' @param N  numeric. The number of observations. Only needed if x is a
#' correlation matrix.
#' @param use character. Passed to [stats::cor()] if raw data is given as input. Default
#'   is `"pairwise.complete.obs"`. `"CD"` always uses only the complete cases: it removes
#'   the incomplete rows with a warning.
#' @param cor_method character. Correlation computed from raw data: `"pearson"`,
#'   `"spearman"`, or `"kendall"` (passed to [stats::cor()]), or `"poly"` / `"tetra"`
#'   for polychoric / tetrachoric correlations (a two-step estimator). With `"poly"` /
#'   `"tetra"`, `"HULL"`, `"NEST"`, `"PARALLEL"`, and `"SMT"` are skipped, and `"CD"` and
#'   Bartlett's test use the Pearson correlations of the raw data (see details). Default
#'   is `"pearson"`.
#' @param n_factors_max numeric. Passed to [efa_cd()]. The maximum number
#' of factors to test against.
#' Larger numbers will increase the duration the procedure takes, but test more
#' possible solutions. If left NA (default), the maximum number of factors for
#' which the model is still over-identified (df > 0) is used.
#' @param N_pop numeric. Passed to [efa_cd()]. Size of finite populations
#' of comparison data. Default is 10000.
#' @param N_samples numeric. Passed to [efa_cd()]. Number of samples drawn
#'  from each population. Default is 500.
#' @param alpha `r lifecycle::badge("deprecated")` Use `alpha_CD` instead.
#' @param alpha_CD numeric. Passed to [efa_cd()] as `alpha`. The alpha level used to
#'   test the significance of the improvement added by an additional factor. Default is
#'   .30.
#' @param max_iter_CD numeric. Passed to [efa_cd()]. The maximum number of
#'  iterations to perform after which the iterative PAF procedure is halted.
#'   Default is 50.
#' @param n_fac_theor numeric. Passed to [efa_hull()]. Theoretical number
#'  of factors to retain. The Hull method uses one plus the larger of this
#'  number and the number of factors suggested by [efa_parallel()] as its
#'  upper bound.
#' @param estimator character. The estimator for the criteria that fit a model: `"ML"`
#'   (default), `"PAF"`, or `"ULS"`. Some criteria have a different default when called
#'   on their own (for example, [efa_nest()] uses `"PAF"`). [efa_kgc()], [efa_scree()],
#'   and [efa_parallel()] use it only for `"EFA"` eigenvalues. [efa_smt()] always uses
#'   maximum likelihood.
#' @param gof character. Passed to [efa_hull()]. The goodness of fit index to use:
#'   `"CAF"` (default), `"CFI"`, or `"RMSEA"`, or any combination of them. With the
#'   `"PAF"` estimator, only the CAF is available. In Lorenzo-Seva, Timmerman, and Kiers
#'   (2011), with factor correlations up to .40, the CFI was the most accurate index. With
#'   strongly correlated factors (e.g., .70), CFI and RMSEA can make the Hull method
#'   retain too few factors. For details on the CAF, see Lorenzo-Seva et al. (2011).
#' @param eigen_type_HULL character. Passed to [efa_parallel()] in [efa_hull()]. The
#'   eigenvalues of the parallel analysis: `"PCA"` (default) uses the correlation
#'   matrix, `"SMC"` puts the squared multiple correlations on its diagonal, and `"EFA"`
#'   puts the communalities of an EFA solution on its diagonal (see details).
#' @param decision_rule_HULL character. Passed to [efa_hull()] as `decision_rule`. The
#'   rule of the parallel analysis that sets the upper bound of the Hull method, with the
#'   same choices as `decision_rule`. Default is `"percentile"` (see details).
#' @param eigen_type_other character. Passed to [efa_kgc()], [efa_scree()], and
#'   [efa_parallel()]. The same types as `eigen_type_HULL`, but one or more of them.
#'   Default is `"PCA"`.
#' @param n_factors numeric. Passed to [efa_parallel()] (also within
#' [efa_hull()]), [efa_kgc()], and [efa_scree()]. Number of
#' factors to extract if `"EFA"` is included in `eigen_type_HULL` or
#'  `eigen_type_other`. Default is 1.
#' @param n_datasets numeric. Passed to [efa_parallel()] (also within
#' [efa_hull()]). The number of reference correlation matrices to simulate. Default is
#' 1000.
#' @param percent numeric. Passed to [efa_parallel()] (also within
#' [efa_hull()]). The percentile to take from the simulated eigenvalues.
#'  Default is 95.
#' @param decision_rule character. Passed to [efa_parallel()] for the `"PARALLEL"`
#'  criterion only (the Hull method uses `decision_rule_HULL`). Which rule to use to
#'  determine the number of factors to retain. Default is `"means"`, which uses the
#'  average simulated eigenvalues. `"percentile"` uses the percentiles specified in
#'  `percent`. `"crawford"` uses the 95th percentile for the first factor and the mean
#'  afterwards (based on Crawford et al., 2010).
#' @param ekc_type `r lifecycle::badge("deprecated")` Accepted and ignored. It
#'   used to select between two ways to compute the [efa_ekc()] reference values.
#'   The `"AM2019"` reference values do not depend on the observed eigenvalues.
#'   They therefore skip the empirical correction that defines the criterion, so
#'   they are no longer computed.
#' @param n_datasets_nest numeric. Passed to [efa_nest()]. The number of reference
#'   correlation matrices to simulate. Default is 1000.
#' @param alpha_nest numeric. Passed to [efa_nest()]. The alpha level to use.
#'   The reference values are the eigenvalues at the (1 - alpha_nest)
#'   percentile. Default is .05.
#' @param show_progress logical. Whether a progress bar should be shown in the
#'   console. Default is FALSE.
#' @param estimate_control an [estimate_control()] object with the estimation settings
#'   for the unrotated [efa_fit()] fits of the criteria that fit a model (all except
#'   [efa_cd()], [efa_ekc()], and [efa_map()]). `NULL` (default) uses the [efa_fit()]
#'   defaults. [efa_smt()] always uses maximum likelihood, so only `start_method` applies
#'   there.
#' @param ... Further arguments passed to [efa_fit()] in the criteria that fit a model,
#'   except [efa_smt()]. Unknown arguments, rotation settings, estimation settings (use
#'   `estimate_control`), and the standard-error arguments (`se`, `b_boot`, `ci`, `seed`)
#'   cause an error. Arguments after `...` must be given by their full name.
#'
#' @details
#' By default, `efa_retain()` checks if the data are suitable for factor analysis. It
#' uses Bartlett's test of sphericity ([efa_bartlett()]), which needs `N`, and the
#' Kaiser-Meyer-Olkin criterion ([efa_kmo()]). A significant Bartlett's test only shows
#' that the correlation matrix is not an identity matrix. This is almost always true for
#' realistic sample sizes. With `cor_method = "poly"` or `"tetra"` and raw data,
#' Bartlett's test uses the Pearson correlations of the raw data. Its chi-square
#' distribution is not valid for polychoric correlations, because these have a larger
#' sampling variance. Zero polychoric correlations give zero Pearson correlations, so the
#' test has the same null hypothesis. The KMO uses the analysed correlations.
#'
#' The available factor retention criteria are the following (see the
#'  respective documentation for details):
#'  \itemize{
#'   \item Comparison data (see [efa_cd()])
#'   \item Empirical Kaiser criterion (see [efa_ekc()])
#'   \item Hull method (see [efa_hull()])
#'   \item Kaiser-Guttman criterion (see [efa_kgc()])
#'   \item Velicer's minimum average partial, MAP (see [efa_map()])
#'   \item Next Eigenvalue Sufficiency Test, NEST (see [efa_nest()])
#'   \item Parallel analysis (see [efa_parallel()])
#'   \item Scree plot (see [efa_scree()])
#'   \item Sequential chi-square model tests, RMSEA lower bound, and AIC
#'     (see [efa_smt()])
#' }
#'
#' The default `criteria` are comparison data, the empirical Kaiser criterion, the
#' Hull method, MAP, NEST, and parallel analysis. No single criterion is the most
#' accurate in all conditions. `efa_retain()` therefore runs several criteria
#' together, and the printed summary gives the range of their suggestions and the
#' most common one (one vote per criterion). Auerswald and Moshagen (2019) compare the
#' criteria and give guidance on the selection.
#'
#' Parallel analysis, also in the Hull method, uses PCA eigenvalues by default. With
#' the eigenvalues of the reduced correlation matrix (`"SMC"`), parallel analysis often
#' retains too many factors, particularly when the model does not fit the data exactly.
#' In Lim and Jahng (2019), parallel analysis with PCA eigenvalues was the most accurate
#' in most conditions. The `"PARALLEL"` criterion compares the observed eigenvalues with
#' the mean simulated eigenvalues by default. To set the upper bound of the Hull method,
#' parallel analysis uses the 95th percentile instead, as Lorenzo-Seva et al. (2011)
#' suggest. This helps to prevent overfactoring (see [efa_hull()]). The Hull method uses
#' the CAF by default (see `gof`).
#'
#' A criterion that cannot use the input is skipped: `"CD"` with a correlation matrix,
#' and `"HULL"`, `"NEST"`, `"PARALLEL"`, and `"SMT"` with `cor_method = "poly"` or
#' `"tetra"`. If you do not set `criteria`, this occurs without a warning. With a
#' correlation matrix and no `N`, the criteria that need `N` (`"EKC"`, `"HULL"`,
#' `"NEST"`, `"PARALLEL"`, and `"SMT"`) are skipped with a warning. The printed result
#' shows all criteria that did not run. With `cor_method = "poly"` or `"tetra"`,
#' comparison data use the Pearson correlations of the raw data (see [efa_cd()]).
#'
#' Comparison data, parallel analysis, NEST, and the Hull method use simulated reference
#' data, so their suggestions can change slightly from run to run. The Hull method uses
#' simulation only for the parallel analysis that sets its upper bound. If parallel
#' analysis also runs and `eigen_type_other` includes `eigen_type_HULL`, the Hull method
#' uses its simulated eigenvalues with `decision_rule_HULL`, and does not simulate again.
#' To make the results reproducible, call [base::set.seed()] first.
#'
#' @returns A list of class `c("efa_retain", "N_FACTORS")`, the trailing class
#'   keeping `inherits(x, "N_FACTORS")` working for code written against the
#'   superseded name. It contains
#' \item{suitability}{A list with the results from [efa_bartlett()] and
#'   [efa_kmo()] (`bartlett` and `kmo`), or `NULL` if
#'   `suitability = FALSE`. `bartlett` is `NULL` if `N` is missing.}
#' \item{outputs}{A named list with one `efa_retention` object per factor
#'   retention criterion that was run, in alphabetical order (see, e.g., [efa_ekc()]).}
#' \item{n_factors}{A named numeric vector with the suggested number of factors per
#'   criterion and variant (e.g. `EKC_BvA2017` or `PARALLEL_PCA`), without the scree plot.}
#' \item{not_run}{A named character vector with the criteria that were skipped
#'   or failed and the reason, or `NULL` if all requested criteria ran.}
#' \item{settings}{A list of the settings used. `gof_used` records the Hull goodness of
#'   fit indices that were computed (`NA` if the Hull method did not run).}
#'
#' @source Auerswald, M., & Moshagen, M. (2019). How to determine the number of
#' factors to retain in exploratory factor analysis: A comparison of extraction
#' methods under realistic conditions. Psychological Methods, 24(4), 468–491.
#' https://doi.org/10.1037/met0000200
#'
#' @source Lim, S., & Jahng, S. (2019). Determining the number of factors using
#' parallel analysis and its recent variants. Psychological Methods, 24(4), 452–467.
#' https://doi.org/10.1037/met0000230
#'
#' @source Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011). The Hull method
#' for selecting the number of common factors. Multivariate Behavioral Research, 46(2),
#' 340-364. https://doi.org/10.1080/00273171.2011.564527
#'
#' @source Ruscio, J., & Roche, B. (2012). Determining the number of factors to
#' retain in an exploratory factor analysis using comparison data of known
#' factorial structure. Psychological Assessment, 24, 282–292.
#' https://doi.org/10.1037/a0025697
#'
#' @seealso [efa_screen()] for data screening before retention, and [efa_fit()] to extract
#'  the chosen number of factors.
#'
#' @family factor retention criteria
#'
#' @export
#'
#' @examples
#' \donttest{
#' # The simulation-based criteria are seeded to make the runs reproducible
#' set.seed(42)
#'
#' # Correlation matrix: the default criteria except "CD", which needs raw data
#' efa_retain(test_models$baseline$cormat, N = 500, n_datasets = 100,
#'            n_datasets_nest = 100)
#'
#' # Raw data: "CD" runs too
#' efa_retain(GRiPS_raw, N_pop = 500, N_samples = 20, n_datasets = 100,
#'            n_datasets_nest = 100)
#' }
efa_retain <- function(x, criteria = c("CD", "EKC", "HULL", "MAP", "NEST", "PARALLEL"),
                       suitability = TRUE, N = NA,
                       use = c("pairwise.complete.obs", "all.obs",
                               "complete.obs", "everything", "na.or.complete"),
                       cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                       n_factors_max = NA, N_pop = 10000, N_samples = 500,
                       alpha = lifecycle::deprecated(),
                       # arguments after `...` must be given by their full name
                       ...,
                       alpha_CD = .30, max_iter_CD = 50, n_fac_theor = NA,
                       estimator = c("ML", "PAF", "ULS"),
                       gof = "CAF",
                       eigen_type_HULL = c("PCA", "SMC", "EFA"),
                       decision_rule_HULL = c("percentile", "means", "crawford"),
                       eigen_type_other = "PCA",
                       n_factors = 1, n_datasets = 1000,
                       percent = 95,
                       decision_rule = c("means", "percentile", "crawford"),
                       ekc_type = lifecycle::deprecated(),
                       n_datasets_nest = 1000, alpha_nest = .05,
                       show_progress = FALSE,
                       estimate_control = NULL){

  # default criteria that cannot use the input are skipped without a warning
  default_criteria <- missing(criteria) || is.null(criteria)

  # refuse here the dots that the criterion fits cannot take (see R/control.R)
  .reject_flat_knobs(...names(), fn = "efa_retain")
  .reject_unknown_fit_dots(...names(), fn = "efa_retain", unrotated = TRUE)
  .reject_rotation_dots(list(...), fn = "efa_retain")
  .assert_cor_input(x)
  .assert_estimate_control(estimate_control)

  if (lifecycle::is_present(alpha)) {
    if (!missing(alpha_CD)) {
      cli::cli_abort(
        c("{.arg alpha} and {.arg alpha_CD} cannot both be supplied.",
          "i" = "{.arg alpha} is the former name of {.arg alpha_CD}; pass only one."),
        class = "efa_renamed_arg"
      )
    }
    lifecycle::deprecate_warn(when = "1.2.0", what = "efa_retain(alpha)",
                              with = "efa_retain(alpha_CD)")
    alpha_CD <- alpha
  }

  criteria <- .match_arg_ci(criteria, several.ok = TRUE,
                            choices = c("CD", "EKC", "HULL", "KGC", "PARALLEL",
                                        "SCREE", "SMT", "NEST", "MAP"))
  eigen_type_HULL <- .match_arg_ci(eigen_type_HULL)
  eigen_type_other <- .match_arg_ci(eigen_type_other, several.ok = TRUE,
                                    choices = c("PCA", "SMC", "EFA"))
  gof <- .match_arg_ci(gof, several.ok = TRUE, choices = c("CAF", "CFI", "RMSEA"))
  cor_method <- .match_arg_ci(cor_method)
  use <- .match_arg_ci(use)
  estimator <- .match_arg_ci(estimator)
  decision_rule <- .match_arg_ci(decision_rule)
  decision_rule_HULL <- .match_arg_ci(decision_rule_HULL)
  if (lifecycle::is_present(ekc_type)) .deprecate_ekc_type("efa_retain(ekc_type)")

  # the checks of the criteria, done before any of them runs
  .assert_args({
    checkmate::assert_flag(suitability)
    checkmate::assert_count(N, na.ok = TRUE, positive = TRUE)
    checkmate::assert_count(n_factors_max, na.ok = TRUE, positive = TRUE)
    checkmate::assert_int(N_pop, lower = 2)
    checkmate::assert_count(N_samples, positive = TRUE)
    checkmate::assert_number(alpha_CD, lower = 0, upper = 1)
    checkmate::assert_count(max_iter_CD, positive = TRUE)
    checkmate::assert_count(n_fac_theor, na.ok = TRUE)
    checkmate::assert_count(n_factors, positive = TRUE)
    checkmate::assert_count(n_datasets, positive = TRUE)
    checkmate::assert_number(percent, lower = 0, upper = 100)
    checkmate::assert_count(n_datasets_nest, positive = TRUE)
    checkmate::assert_number(alpha_nest, lower = 0, upper = 1)
    checkmate::assert_flag(show_progress)
  })

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, N = N, use = use, cor_method = cor_method,
                             N_policy = "optional", inform_from_data = FALSE)
  R <- prep$R
  N <- prep$N

  ## Tests for suitability of factor analysis
  suitability_out <- NULL
  if (suitability) {
    # Bartlett's chi-square is not valid for polychoric correlations, so with raw data
    # it tests the Pearson correlations, which are zero if the polychoric ones are
    # (data.matrix() codes ordered factors by level, as .polychoric() does)
    bartlett_out <- if (is.na(N)) {
      NULL
    } else if (.is_poly_cor(cor_method) && !prep$is_cormat) {
      # a smoothing of this auxiliary matrix is not a smoothing of the analysed one
      pearson <- withCallingHandlers(
        .prepare_cor_input(data.matrix(x), use = use, inform_from_data = FALSE),
        efa_cor_smoothed = function(w) invokeRestart("muffleWarning"))
      efa_bartlett(pearson$R, N = pearson$N, use = use, cor_method = "pearson")
    } else {
      efa_bartlett(R, N = N, use = use, cor_method = cor_method)
    }
    suitability_out <- list(
      bartlett = bartlett_out,
      kmo = efa_kmo(R, use = use, cor_method = cor_method)
    )
  }

  ## Factor retention criteria, driven by the registry
  ctl <- .n_factors_ctl(N = N, use = use, cor_method = cor_method,
                        n_factors_max = n_factors_max, N_pop = N_pop,
                        N_samples = N_samples, alpha_CD = alpha_CD,
                        max_iter_CD = max_iter_CD, n_fac_theor = n_fac_theor,
                        estimator = estimator, gof = gof,
                        eigen_type_HULL = eigen_type_HULL,
                        decision_rule_HULL = decision_rule_HULL,
                        eigen_type_other = eigen_type_other,
                        n_factors = n_factors, n_datasets = n_datasets,
                        percent = percent, decision_rule = decision_rule,
                        n_datasets_nest = n_datasets_nest,
                        alpha_nest = alpha_nest,
                        estimate_control = estimate_control, dots = list(...))

  run <- intersect(names(.retention_registry), criteria)
  # the Hull method reuses the parallel analysis, so that one runs first
  hull_from_pa <- all(c("HULL", "PARALLEL") %in% run) && eigen_type_HULL %in% eigen_type_other
  if (hull_from_pa) {
    run <- setdiff(run, "PARALLEL")
    run <- append(run, "PARALLEL", after = match("HULL", run) - 1L)
  }

  outputs <- list()
  not_run <- character(0)   # id -> reason it was skipped or failed
  kind <- character(0)      # id -> "silent" or "skipped" (both not run), or "failed"
  first_cnd <- NULL         # the first criterion error, chained into the abort below
  no_n <- "needs N, which was not supplied"

  for (id in run) {

    if (show_progress) {
      # cli evaluates the label when the step ends, after `id` changed, so build it now
      cli::cli_progress_step(paste0("Running ", id))
    }

    entry <- .retention_registry[[id]]

    skip <- if (isTRUE(entry$needs_raw) && prep$is_cormat) {
      "needs raw data, but a correlation matrix was supplied"
    } else if (.is_poly_cor(cor_method) && isFALSE(entry$poly_ok)) {
      paste0("does not support cor_method = \"", cor_method, "\"")
    } else if (isTRUE(entry$needs_N) && is.na(N)) {
      no_n
    }
    if (!is.null(skip)) {
      not_run[[id]] <- skip
      # a missing N is likely a mistake, so it is reported also for a default criterion
      kind[[id]] <- if (default_criteria && skip != no_n) "silent" else "skipped"
      next
    }

    if (id == "HULL" && hull_from_pa && !is.null(outputs[["PARALLEL"]])) {
      pa <- outputs[["PARALLEL"]]
      ctl$n_fac_PA <- if (decision_rule_HULL == decision_rule) {
        unname(pa$n_factors[[eigen_type_HULL]])
      } else {
        # the same simulated eigenvalues, compared with the rule of the Hull method
        rec <- pa$results[[match(eigen_type_HULL, names(pa$n_factors))]]
        .determine_factors(decision_rule_HULL, rec$y, do.call(cbind, rec$references),
                           percent)
      }
    }

    out_id <- tryCatch(entry$fun(if (isTRUE(entry$needs_raw)) x else R, ctl),
                       error = function(e) e)

    if (inherits(out_id, "error")) {
      # the unwrapped cli headline, or the formatted message if the headline is empty
      reason <- trimws(paste(out_id$message, collapse = " "))
      if (!nzchar(reason)) reason <- conditionMessage(out_id)
      not_run[[id]] <- gsub("\\s+", " ", trimws(reason))
      kind[[id]] <- "failed"
      if (is.null(first_cnd)) first_cnd <- out_id
      next
    }

    outputs[[id]] <- out_id

  }

  if (show_progress) {
    cli::cli_progress_done()
  }

  # report in registry order, whatever order the criteria ran in
  outputs <- outputs[intersect(names(.retention_registry), names(outputs))]
  ord <- intersect(names(.retention_registry), names(not_run))
  not_run <- not_run[ord]
  kind <- kind[ord]

  # nothing ran: a result with no criteria is not meaningful
  if (length(outputs) == 0) {
    cli::cli_abort(
      c("None of the requested factor retention criteria could be run.",
        .retain_reason_bullets(not_run)),
      class = "efa_no_criteria",
      parent = first_cnd
    )
  }

  skipped <- not_run[kind == "skipped"]
  if (length(skipped) > 0) {
    n_skip <- length(skipped)
    cli::cli_warn(
      c("{n_skip} factor retention criteri{?on was/a were} skipped:",
        .retain_reason_bullets(skipped),
        if (any(skipped == no_n)) {
          c("i" = "Supply {.arg N}, or raw data, to run the criteria that need it.")
        }),
      class = "efa_criterion_skipped"
    )
  }

  failed <- not_run[kind == "failed"]
  if (length(failed) > 0) {
    n_fail <- length(failed)
    cli::cli_warn(
      c("{n_fail} factor retention criteri{?on/a} failed:",
        .retain_reason_bullets(failed)),
      class = "efa_criterion_failed"
    )
  }

  settings <- list(criteria = criteria,
                   suitability = suitability,
                   N = N,
                   use = use,
                   n_factors_max = n_factors_max,
                   N_pop = N_pop,
                   N_samples = N_samples,
                   alpha_CD = alpha_CD,
                   # back-compat alias of alpha_CD
                   alpha = alpha_CD,
                   cor_method = cor_method,
                   max_iter_CD = max_iter_CD,
                   n_fac_theor = n_fac_theor,
                   estimator = estimator,
                   # back-compat alias, as for the frozen P_type key
                   method = estimator,
                   gof = gof,
                   # the Hull indices actually computed (PAF reduces them to the CAF)
                   gof_used = if (is.null(outputs[["HULL"]])) NA_character_
                              else outputs[["HULL"]]$settings$gof,
                   eigen_type_HULL = eigen_type_HULL,
                   decision_rule_HULL = decision_rule_HULL,
                   eigen_type_other = eigen_type_other,
                   n_factors = n_factors,
                   n_datasets = n_datasets,
                   percent = percent,
                   decision_rule = decision_rule,
                   n_datasets_nest = n_datasets_nest,
                   alpha_nest = alpha_nest,
                   estimate_control = estimate_control)

  # Aggregate the suggested numbers of factors as "<id>" or "<id>_<variant>".
  # Visual criteria (the scree plot) make no numeric suggestion and are omitted;
  # for the others NA suggestions are kept (named), so a criterion that ran but
  # could not determine a number stays visible.
  n_factors_out <- unlist(lapply(names(outputs), function(id) {
    if (isTRUE(.retention_registry[[id]]$visual)) return(NULL)
    .retention_key(id, outputs[[id]]$n_factors)
  }))
  if (is.null(n_factors_out)) n_factors_out <- numeric(0)

  output <- list(suitability = suitability_out,
                 outputs = outputs,
                 n_factors = n_factors_out,
                 not_run = if (length(not_run) > 0) not_run else NULL,
                 settings = settings)

  # the trailing legacy class is load-bearing: it keeps `inherits(x, "N_FACTORS")`
  # working in code written against the superseded name. Do not drop it.
  class(output) <- c("efa_retain", "N_FACTORS")

  return(output)

}

# One cli bullet "<id>: <reason>" per element of a named character vector, with the
# braces doubled so that cli shows a reason (often an error message) as it is.
.retain_reason_bullets <- function(reasons) {
  txt <- paste0(names(reasons), ": ", reasons)
  txt <- gsub("}", "}}", gsub("{", "{{", txt, fixed = TRUE), fixed = TRUE)
  stats::setNames(txt, rep("*", length(txt)))
}

#' Print method for efa_retain objects
#'
#' @param x an object of class efa_retain, returned by [efa_retain()].
#' @param ... not used.
#'
#' @returns `print()` returns its argument `x` invisibly; it is
#'   `cat(format(x), sep = "\n")`.
#'
#' @export
#' @method print efa_retain
#'
#' @examples
#' \donttest{
#' efa_retain(test_models$baseline$cormat, criteria = c("EKC", "SMT"), N = 500)
#' }
print.efa_retain <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' Format method for efa_retain objects
#'
#' @param x an object of class efa_retain, returned by [efa_retain()].
#' @param ... not used.
#'
#' @returns A character vector with the report lines (styled to the active
#'   console theme; plain when colours are disabled).
#'
#' @export
#' @method format efa_retain
#'
#' @examples
#' \donttest{
#' nf <- efa_retain(test_models$baseline$cormat, criteria = c("EKC", "SMT"),
#'                  N = 500)
#' writeLines(format(nf))
#' }
format.efa_retain <- function(x, ...) {
  cli::cli_format_method({

    if (!is.null(x$suitability)) {

      cli::cli_rule(left = "Sampling adequacy and sphericity")
      cli::cli_text("")

      # the same lines as in format.efa_screen()
      bart <- x$suitability$bartlett
      pearson <- !is.null(bart) && .is_poly_cor(x$settings$cor_method) &&
        identical(bart$settings$cor_method, "pearson")
      .print_kmo_line(x$suitability$kmo$KMO, digits = 3)
      .print_bartlett_line(bart)
      if (pearson) {
        cli::cli_alert_info(paste("Bartlett's test uses Pearson correlations: its chi-square",
                                  "is not valid for polychoric correlations."), wrap = TRUE)
      }

      cli::cli_text("")
    }

    is_visual <- vapply(names(x$outputs),
                        function(id) isTRUE(.retention_registry[[id]]$visual),
                        logical(1))

    cli::cli_rule(left = "Suggested number of factors")

    # the range and most common suggestion, interpolated as a variable so that cli does
    # not parse it as a template
    summary_line <- .retention_summary(x$outputs[!is_visual])
    if (!is.null(summary_line)) {
      cli::cli_text("")
      cli::cli_text("{summary_line}")
    }

    # one group per criterion with a numeric suggestion: label, context line, bullets
    for (out in x$outputs[!is_visual]) {
      cli::cli_text("")
      cli::cli_text("{out$criterion[['label']]}")
      if (!is.null(out$subtitle)) {
        cli::cli_text("{out$subtitle}")
      }
      cli::cli_ul(.retention_bullets(out$results))
    }

    # visual criteria (e.g. the scree plot) point to the plot
    for (out in x$outputs[is_visual]) {
      cli::cli_text("")
      label <- out$criterion[["label"]]
      cli::cli_bullets(c("i" =
        "{label} provides no numeric suggestion; inspect the plot."))
    }

    # criteria that were skipped or failed
    if (!is.null(x$not_run)) {
      cli::cli_text("")
      cli::cli_rule(left = "Criteria not run")
      cli::cli_text("")
      cli::cli_bullets(.retain_reason_bullets(x$not_run))
    }

  })
}

# One-line summary of the suggestions for format.efa_retain(): their number, range, and, if
# they disagree, the most common one. Each criterion has ONE vote, its own modal variant (none
# if its variants tie), so several variants cannot outvote one. NULL if no criterion
# determined a number, or if a single suggestion would only repeat the bullet below.
.retention_summary <- function(outputs) {

  # per criterion, the numbers it determined (a variant that came out NA determined none)
  per_crit <- lapply(outputs, function(o) o$n_factors[!is.na(o$n_factors)])
  per_crit <- per_crit[lengths(per_crit) > 0]
  if (length(per_crit) == 0) return(NULL)

  nf <- unlist(per_crit, use.names = FALSE)
  n_sug <- length(nf)
  n_crit <- length(per_crit)
  # there is nothing to summarise about a lone suggestion; the bullet below already
  # carries the number
  if (n_sug == 1L) return(NULL)

  counts <- paste0(n_sug, " ", .screen_plural(n_sug, "suggestion", "suggestions"),
                   " from ", n_crit, " ",
                   .screen_plural(n_crit, "criterion", "criteria"))

  if (min(nf) == max(nf)) {
    return(paste0(counts, ", all suggesting ", .retention_count(nf[1]), " ",
                  .screen_plural(nf[1], "factor", "factors"), "."))
  }

  # the most frequent value(s) of a vector, with how often the winner occurred; table()
  # names its cells in ascending numeric order, so tied values come out sorted
  modal <- function(v) {
    tab <- table(v)
    list(values = as.numeric(names(tab))[tab == max(tab)], n = max(tab))
  }
  votes <- unlist(lapply(per_crit, function(v) {
    m <- modal(v)
    if (length(m$values) == 1L) m$values else NULL
  }))

  spread <- paste0(counts, ", ranging from ", .retention_count(min(nf)), " to ",
                   .retention_count(max(nf)), " factors")
  if (is.null(votes)) return(paste0(spread, "."))

  # a most common value needs a value that more than one criterion chose
  best <- modal(votes)
  if (best$n == 1L) return(paste0(spread, "."))

  paste0(spread, " (most common: ",
         .screen_and_list(.retention_count(best$values)), ").")

}

#' Plot method for efa_retain objects
#'
#' Plots every factor-retention criterion in the [efa_retain()] result that has
#' a plottable outcome (see [plot.efa_retention()]). Criteria without a plot
#' (e.g. [efa_map()] or [efa_smt()]) are skipped.
#'
#' @param x an object of class efa_retain, returned by [efa_retain()].
#' @param ... not used.
#'
#' @returns A named list of [ggplot2::ggplot] objects, one per criterion with a
#'   plottable result, or invisibly `NULL` if there is none.
#'
#' @export
#' @method plot efa_retain
#'
#' @examples
#' \donttest{
#' nf <- efa_retain(test_models$baseline$cormat, criteria = c("EKC", "SMT"),
#'                  N = 500)
#' plot(nf)
#' }
plot.efa_retain <- function(x, ...) {

  # plot.efa_retention returns NULL (with a message) for criteria with no plot
  # (e.g. MAP/SMT) or a degenerate one (e.g. CD suggesting 0 factors); drop those
  plots <- lapply(x$outputs, function(o) suppressMessages(plot(o)))
  plots <- Filter(Negate(is.null), plots)

  if (length(plots) == 0) {
    cli::cli_inform("No plot is available for the criteria that were run.",
                    class = "efa_no_plot")
    return(invisible(NULL))
  }

  plots

}
