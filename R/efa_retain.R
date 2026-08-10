#' Various factor retention criteria
#'
#' Among the most important decisions for an exploratory factor analysis (EFA) is
#' the choice of the number of factors to retain. Several factor retention
#' criteria have been developed for this. With this function, various factor
#'  retention criteria can be performed simultaneously. Additionally, the data
#'  can be checked for their suitability for factor analysis.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations. If `"CD"` is included as a criterion, x must be raw
#'  data.
#' @param criteria character. A vector with the factor retention methods to
#' perform. Possible inputs are: `"CD"`, `"EKC"`, `"HULL"`,
#' `"KGC"`, `"MAP"`, `"NEST"`,`"PARALLEL"`, `"SCREE"`, and `"SMT"`
#' (see details). By default, a subset of often used, well-performing methods
#' are performed.
#' @param suitability logical. Whether the data should be checked for suitability
#' for factor analysis using the Bartlett's test of sphericity and the
#' Kaiser-Meyer-Olkin criterion (see details). Default is `TRUE`.
#' @param N  numeric. The number of observations. Only needed if x is a
#' correlation matrix.
#' @param use character. Passed to [stats::cor()] if raw
#' data is given as input. Default is `"pairwise.complete.obs"`.
#' @param cor_method character. Correlation computed from raw data: `"pearson"`,
#'   `"spearman"`, or `"kendall"` (passed to [stats::cor()]), or `"poly"` /
#'   `"tetra"` for polychoric / tetrachoric correlations (a two-step estimator).
#'   `CD`, `PARALLEL`, `NEST`, and
#'   `HULL` compare against simulated continuous data, and `SMT` relies on a
#'   normal-theory chi-square test; none of these support `"poly"` / `"tetra"`,
#'   so they are skipped in that case.
#' Default is  `"pearson"`.
#' @param n_factors_max numeric. Passed to [efa_cd()]. The maximum number
#' of factors to test against.
#' Larger numbers will increase the duration the procedure takes, but test more
#' possible solutions. If left NA (default), the maximum number of factors for
#' which the model is still over-identified (df > 0) is used.
#' @param N_pop numeric. Passed to [efa_cd()]. Size of finite populations
#' of comparison data. Default is 10000.
#' @param N_samples numeric. Passed to [efa_cd()]. Number of samples drawn
#'  from each population. Default is 500.
#' @param alpha numeric. Passed to [efa_cd()]. The alpha level used to test
#'  the significance of the improvement added by an additional factor.
#'  Default is .30.
#' @param max_iter_CD numeric. Passed to [efa_cd()]. The maximum number of
#'  iterations to perform after which the iterative PAF procedure is halted.
#'   Default is 50.
#' @param n_fac_theor numeric. Passed to [efa_hull()]. Theoretical number
#'  of factors to retain. The maximum of this number and the number of factors
#'  suggested by [efa_parallel()] plus one will be used in the Hull method.
#' @param estimator character. Passed to [efa_fit()] in [efa_hull()],
#' [efa_kgc()], [efa_scree()], [efa_parallel()], and [efa_nest()]. The
#' estimator to use. One of  `"PAF"`, `"ULS"`, or  `"ML"`,
#' for principal axis factoring, unweighted least squares, and maximum
#' likelihood, respectively. The default here is `"ML"`, which is not
#' the default of every criterion called standalone (for example [efa_hull()] defaults to
#' `"PAF"`), so a criterion run through `efa_retain()` can differ from the same criterion
#' called directly unless `estimator` is set to match. In [efa_kgc()], [efa_scree()], and
#' [efa_parallel()] it only takes effect when the respective `eigen_type` includes `"EFA"`.
#' @param gof character. Passed to [efa_hull()]. The goodness of fit index
#' to use. Either `"CAF"`, `"CFI"`, or `"RMSEA"`, or any
#' combination of them. With the `"PAF"` estimator, only
#' the CAF can be used as goodness of fit index. For details on the CAF, see
#' Lorenzo-Seva, Timmerman, and Kiers (2011).
#' @param eigen_type_HULL character. Passed to  [efa_parallel()] in
#' [efa_hull()]. On what the
#' eigenvalues should be found in the parallel analysis. Can be one of
#' `"SMC"`, `"PCA"`, or `"EFA"`. If using  `"SMC"` (default),
#' the diagonal of the correlation matrices is
#' replaced by the squared multiple correlations (SMCs) of the indicators. If
#' using  `"PCA"`, the diagonal values of the correlation
#' matrices are left to be 1. If using  `"EFA"`, eigenvalues are found on the
#' correlation  matrices with the final communalities of an EFA solution as
#' diagonal.
#' @param eigen_type_other character. Passed to [efa_kgc()],
#' [efa_scree()], and [efa_parallel()]. The same as eigen_type_HULL,
#' but multiple inputs are possible here (any combination of `"PCA"`, `"SMC"`,
#' and `"EFA"`). Default is `"SMC"`.
#' @param n_factors numeric. Passed to [efa_parallel()] (also within
#' [efa_hull()]), [efa_kgc()], and [efa_scree()]. Number of
#' factors to extract if `"EFA"` is included in `eigen_type_HULL` or
#'  `eigen_type_other`. Default is 1.
#' @param n_datasets numeric. Passed to [efa_parallel()] (also within
#' [efa_hull()]). The number of datasets to simulate. Default is 1000.
#' @param percent numeric. Passed to [efa_parallel()] (also within
#' [efa_hull()]). The percentile to take from the simulated eigenvalues.
#'  Default is 95.
#' @param decision_rule character. Passed to [efa_parallel()] (also within
#'  [efa_hull()]). Which rule to use to determine the number of
#'  factors to retain. Default is `"means"`, which will use the average
#'  simulated eigenvalues. `"percentile"`, uses the percentiles specified
#'  in percent. `"crawford"` uses the 95th percentile for the first factor
#'  and the mean afterwards (based on Crawford et al, 2010).
#' @param ekc_type character. Passed to the `type` argument of [efa_ekc()].
#'   Either `"BvA2017"` for the original implementation by Braeken and van Assen
#'   (2017), or `"AM2019"` for the adapted implementation by Auerswald and Moshagen
#'   (2019).
#' @param n_datasets_nest numeric. The number of datasets to simulate in [efa_nest()]. Default is 1000.
#' @param alpha_nest numeric. The alpha level to use in [efa_nest()] (i.e., 1-alpha percentile of eigenvalues is used for reference values).
#' @param show_progress logical. Whether a progress bar should be shown in the
#'   console. Default is FALSE.
#' @param estimate_control an [estimate_control()] object with the estimation settings for the
#'  [efa_fit()] fits run by the criteria that fit a model ([efa_hull()], [efa_kgc()],
#'  [efa_scree()], [efa_parallel()], [efa_nest()], and [efa_smt()]). `NULL` (default) uses the
#'  [efa_fit()] defaults. [efa_cd()], [efa_ekc()], and [efa_map()] run no [efa_fit()] model, so
#'  it does not apply to them. In [efa_kgc()], [efa_scree()], and [efa_parallel()] it only takes
#'  effect when the respective `eigen_type` includes `"EFA"`, since no model is fitted otherwise,
#'  and [efa_smt()] fits with maximum likelihood by definition, so only `start_method` takes
#'  effect there. All fits are unrotated, so no rotation settings apply.
#' @param ... Further arguments passed to [efa_fit()] in
#' [efa_parallel()], [efa_kgc()], [efa_scree()], and [efa_nest()]. They also reach
#' [efa_hull()], both through the parallel analysis it runs to set its upper bound and
#' through its own candidate fits, so an argument that [efa_fit()] rejects stops the Hull
#' method even when its `eigen_type_HULL` fits no model.
#' The estimation tuning knobs are not passed here; they live in `estimate_control`, and the
#' standard-error arguments (`se`, `b_boot`, `ci`, `seed`) are not accepted because the
#' criterion fits are internal steps whose standard errors are not reported. Note
#' that the arguments listed after `...` must be given by their full name (R matches an abbreviated
#' name only against the arguments before `...`), so that a tuning knob such as `max_iter` cannot
#' be mistaken for `max_iter_CD`.
#'
#' @details
#' By default, the entered data are checked for suitability for factor analysis
#' using the following methods (see respective documentations for details):
#' \itemize{
#'   \item Bartlett's test of sphericity (see [efa_bartlett()])
#'   \item Kaiser-Meyer-Olkin criterion (see [efa_kmo()])
#' }
#'
#' The available factor retention criteria are the following (see respective
#'  documentations for details):
#'  \itemize{
#'   \item Comparison data (see [efa_cd()])
#'   \item Empirical Kaiser criterion (see [efa_ekc()])
#'   \item Hull method (see [efa_hull()])
#'   \item Kaiser-Guttman criterion (see [efa_kgc()])
#'   \item Parallel analysis (see [efa_parallel()])
#'   \item Next Eigenvalue Sufficiency Test, NEST (see [efa_nest()])
#'   \item Scree plot (see [efa_scree()])
#'   \item Sequential chi-square model tests, RMSEA lower bound, and AIC
#'     (see [efa_smt()])
#' }
#'
#' The comparison data, parallel analysis, and NEST criteria compare the data against
#' simulated reference data, so their suggested numbers of factors vary slightly from run
#' to run; the Hull method does too, because it calls [efa_parallel()] to set its upper
#' bound. Call [set.seed()] before `efa_retain()` to make them reproducible.
#'
#' @returns A list of class `c("efa_retain", "N_FACTORS")`, the trailing class
#'   keeping `inherits(x, "N_FACTORS")` working for code written against the
#'   superseded name. It contains
#' \item{suitability}{A list with the results from [efa_bartlett()] and
#'   [efa_kmo()] (`bartlett` and `kmo`), or `NULL` if
#'   `suitability = FALSE`.}
#' \item{outputs}{A named list with one `efa_retention` object per factor
#'   retention criterion that was run (see, e.g., [efa_ekc()]).}
#' \item{n_factors}{A named numeric vector with the suggested number of factors
#'   per criterion and, where a criterion has several variants, per variant
#'   (e.g. `EKC_BvA2017` or `PARALLEL_SMC`). Criteria without a numeric
#'   suggestion (the scree plot) are not included. The "most common" value the printed
#'   summary line reports is counted with one vote per *criterion* -- each criterion's own
#'   modal variant -- so that a criterion with several variants cannot outvote a
#'   single-variant one. A criterion whose variants tie has no modal value and abstains
#'   rather than casting a vote for each of them, and a value that no two criteria agree on
#'   is not reported as most common. Tallying this vector directly counts one vote per
#'   *variant* and can therefore have a different mode.}
#' \item{not_run}{A named character vector with the criteria that were skipped
#'   or failed and the reason, or `NULL` if all requested criteria ran.}
#' \item{settings}{A list of the settings used. Its `criteria` element records the
#'   requested criteria, in the order they were given, while `outputs` and `n_factors`
#'   are in the order in which the criteria were run. `gof` records the requested Hull
#'   goodness-of-fit indices and `gof_used` the ones the Hull method actually computed
#'   (it reduces them to `"CAF"` for the PAF estimator); `gof_used` is `NA` when HULL was
#'   not requested, was skipped, or failed.}
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
#' # Default criteria, with correlation matrix and estimator "ML" (where needed)
#' # This will throw a warning for CD, as no raw data were specified
#' # The simulation-based criteria are seeded to make the run reproducible
#' set.seed(42)
#' nfac_all <- efa_retain(test_models$baseline$cormat, N = 500, estimator = "ML",
#'                        n_datasets = 100, n_datasets_nest = 100)
#'
#' # The same as above, but without "CD"
#' nfac_wo_CD <- efa_retain(test_models$baseline$cormat, criteria = c("EKC",
#'                          "HULL", "PARALLEL", "NEST"), N = 500,
#'                          estimator = "ML", n_datasets = 100,
#'                          n_datasets_nest = 100)
#'
#' # Use PAF instead of ML (this will take longer). For this, gof has
#' # to be set to "CAF" for the Hull method.
#' nfac_PAF <- efa_retain(test_models$baseline$cormat, criteria = c("EKC",
#'                        "HULL", "PARALLEL", "NEST"), N = 500,
#'                        estimator = "PAF", gof = "CAF", n_datasets = 100,
#'                        n_datasets_nest = 100)
#'
#' # Do KGC and PARALLEL with only "PCA" type of eigenvalues
#' nfac_PCA <- efa_retain(test_models$baseline$cormat, criteria = c("EKC",
#'                        "HULL", "PARALLEL", "NEST"), N = 500,
#'                        estimator = "ML", eigen_type_other = "PCA",
#'                        n_datasets = 100, n_datasets_nest = 100)
#'
#' # Use raw data, such that CD can also be performed
#' nfac_raw <- efa_retain(GRiPS_raw, estimator = "ML", N_pop = 500,
#'                         N_samples = 20, n_datasets = 100,
#'                         n_datasets_nest = 100)
#'}
efa_retain <- function(x, criteria = c("CD", "EKC", "HULL", "MAP", "NEST", "PARALLEL"),
                       suitability = TRUE, N = NA,
                       use = c("pairwise.complete.obs", "all.obs",
                               "complete.obs", "everything", "na.or.complete"),
                       cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                       n_factors_max = NA, N_pop = 10000, N_samples = 500,
                       alpha = .30,
                       # `...` sits here deliberately, and must stay ahead of `max_iter_CD`: R
                       # matches an argument to a formal by unique prefix, but only for the
                       # formals BEFORE `...`. With `max_iter_CD` behind it, a bare `max_iter`
                       # can no longer be taken for it -- silently capping the comparison-data
                       # iterations instead of tuning the fits -- and lands in the dots, where
                       # the flat-knob guard below names the control object it belongs to. Every
                       # argument below must therefore be passed by its full name. Keep any new
                       # formal here too unless it can never be confused with a tuning knob.
                       ...,
                       max_iter_CD = 50, n_fac_theor = NA,
                       estimator = c("ML", "PAF", "ULS"),
                       gof = c("CAF", "CFI", "RMSEA"),
                       eigen_type_HULL = c("SMC", "PCA", "EFA"),
                       eigen_type_other = c("SMC"),
                       n_factors = 1, n_datasets = 1000,
                       percent = 95,
                       decision_rule = c("means", "percentile", "crawford"),
                       ekc_type = c("BvA2017"),
                       n_datasets_nest = 1000, alpha_nest = .05,
                       show_progress = FALSE,
                       estimate_control = NULL){

  # Perform argument checks
  # Reject a flat tuning knob here rather than letting it travel into the criteria that fit a
  # model: each criterion is run under a tryCatch, so the knob would surface as an opaque
  # "criterion could not be run" instead of naming the control object it belongs to.
  .reject_flat_knobs(...names(), fn = "efa_retain")
  # A dot the criteria cannot forward to a fit is a misspelling; reject it here rather than
  # letting it be silently dropped (when the selected criteria run no fit) or surface as an
  # opaque "criterion could not be run". The criterion fits are always unrotated, so a
  # rotation setting is refused too (the N_FACTORS() wrapper's repacked rotate_control() object
  # is exempt -- see .reject_rotation_dots()), as are efa_fit()'s inference arguments
  # (`se`/`b_boot`/`ci`/`seed`), which a criterion fit computes and then discards -- see
  # .reject_inference_dots().
  .reject_unknown_fit_dots(...names(), fn = "efa_retain", unrotated = TRUE)
  .reject_rotation_dots(list(...), fn = "efa_retain")
  .assert_cor_input(x)
  .assert_estimate_control(estimate_control)

  ## Perform argument checks and prepare input
  criteria <- .match_arg_ci(criteria, several.ok = TRUE,
                            choices = c("CD", "EKC", "HULL", "KGC", "PARALLEL",
                                        "SCREE", "SMT", "NEST", "MAP"))
  suitability <- checkmate::assert_flag(suitability)
  eigen_type_HULL <- .match_arg_ci(eigen_type_HULL)
  eigen_type_other <- .match_arg_ci(eigen_type_other, several.ok = TRUE,
                                    choices = c("PCA", "SMC", "EFA"))
  gof <- .match_arg_ci(gof, several.ok = TRUE, choices = c("CAF", "CFI", "RMSEA"))
  cor_method <- .match_arg_ci(cor_method)
  use <- .match_arg_ci(use)
  estimator <- .match_arg_ci(estimator)
  decision_rule <- .match_arg_ci(decision_rule)
  ekc_type <- .match_arg_ci(ekc_type, c("BvA2017", "AM2019"), several.ok = TRUE)
  checkmate::assert_number(alpha_nest, lower = 0, upper = 1)
  checkmate::assert_count(n_datasets_nest, na.ok = FALSE, positive = TRUE)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, N = N, use = use, cor_method = cor_method,
                             N_policy = "optional", inform_from_data = FALSE)
  R <- prep$R
  N <- prep$N

  ## Tests for suitability of factor analysis
  suitability_out <- NULL
  if (isTRUE(suitability)) {
    # Bartlett's test of sphericity needs N; the KMO criterion does not. When N
    # is unavailable (a correlation matrix without N), skip Bartlett's test with a
    # note and still report KMO, instead of aborting and thereby blocking the
    # N-free retention criteria too.
    bartlett_out <- if (is.na(N)) {
      cli::cli_warn(
        c("{.arg N} is {.val NA}; Bartlett's test of sphericity was skipped.",
          "i" = "Provide {.arg N} or raw data to include it."),
        class = "efa_suitability_no_n"
      )
      NULL
    } else {
      efa_bartlett(R, N = N, use = use, cor_method = cor_method)
    }
    suitability_out <- list(
      bartlett = bartlett_out,
      kmo = efa_kmo(R, use = use, cor_method = cor_method)
    )
  }

  ## Factor retention criteria, driven by the registry (run in registry order)
  ctl <- .n_factors_ctl(N = N, use = use, cor_method = cor_method,
                        n_factors_max = n_factors_max, N_pop = N_pop,
                        N_samples = N_samples, alpha = alpha,
                        max_iter_CD = max_iter_CD, n_fac_theor = n_fac_theor,
                        estimator = estimator, gof = gof,
                        eigen_type_HULL = eigen_type_HULL,
                        eigen_type_other = eigen_type_other,
                        n_factors = n_factors, n_datasets = n_datasets,
                        percent = percent, decision_rule = decision_rule,
                        ekc_type = ekc_type, n_datasets_nest = n_datasets_nest,
                        alpha_nest = alpha_nest,
                        estimate_control = estimate_control, dots = list(...))

  run <- intersect(names(.retention_registry), criteria)
  outputs <- list()
  not_run <- character(0)   # id -> reason it was skipped or failed
  # the first criterion error, chained into the all-failed abort below so the underlying
  # cause travels with it (a criterion that was skipped raises no condition to chain)
  first_cnd <- NULL

  for (id in run) {

    if (isTRUE(show_progress)) {
      # the message is a glue template evaluated when the bar is terminated -- by which time
      # `id` has advanced to the next criterion -- so build the string eagerly and freeze the
      # criterion name in it, or every completed step is labelled with its successor's name
      cli::cli_progress_step(paste0("Running ", id))
    }

    entry <- .retention_registry[[id]]

    if (isTRUE(entry$needs_raw) && .is_cormat(x)) {
      cli::cli_warn(
        c("{.arg x} is a correlation matrix, but {.val {id}} needs raw data.",
          "i" = "Skipping {.val {id}}."),
        class = "efa_criterion_skipped"
      )
      not_run[[id]] <- "needs raw data, but a correlation matrix was supplied"
      next
    }

    # Some criteria do not support polychoric/tetrachoric correlations (they
    # compare against continuous reference data, or their normal-theory test is
    # invalid for such correlations); skip them with a clear note rather than
    # letting them run on an inappropriate matrix.
    if (.is_poly_cor(cor_method) && isFALSE(entry$poly_ok)) {
      cli::cli_warn(
        c("{.val {id}} does not support {.code cor_method = {.val {cor_method}}}.",
          "i" = "Skipping {.val {id}}."),
        class = "efa_criterion_skipped"
      )
      not_run[[id]] <- paste0("does not support cor_method = \"", cor_method, "\"")
      next
    }

    # a failing criterion is excluded with a warning; the others still run
    out_id <- try(entry$fun(if (isTRUE(entry$needs_raw)) x else R, ctl),
                  silent = TRUE)

    if (inherits(out_id, "try-error")) {
      # `conditionMessage()` returns the message already wrapped to the console width, so its
      # first line is only the first *wrapped* line -- truncating the reason mid-sentence at
      # narrow widths. A cli condition keeps the unwrapped headline in `$message` (the bullets
      # live in `$body`), so take that and squish the newlines the source string carries.
      # A condition whose headline is absent or empty carries its text in `$body` instead, so
      # fall back to the formatted message there rather than reporting an empty reason.
      cnd <- attr(out_id, "condition")
      reason <- trimws(paste(cnd$message, collapse = " "))
      if (!nzchar(reason)) reason <- conditionMessage(cnd)
      reason <- gsub("\\s+", " ", trimws(reason))
      if (is.null(first_cnd)) first_cnd <- cnd
      cli::cli_warn(
        c("{.val {id}} could not be run and is excluded from the results.",
          "i" = "Error: {reason}"),
        class = "efa_criterion_failed"
      )
      not_run[[id]] <- reason
      next
    }

    outputs[[id]] <- out_id

  }

  if (isTRUE(show_progress)) {
    cli::cli_progress_done()
  }

  # nothing ran: a result with no criteria is not meaningful
  if (length(outputs) == 0) {
    cli::cli_abort(
      c("None of the requested factor retention criteria could be run.",
        "x" = "Could not run: {.val {names(not_run)}}.",
        "i" = "See the warnings above for the reason in each case."),
      class = "efa_no_criteria",
      # chains the first criterion error onto the abort, so the cause is reported with it
      # rather than only in the warning that has already scrolled past
      parent = first_cnd
    )
  }

  # Prepare settings here
  settings <- list(criteria = criteria,
                   suitability = suitability,
                   N = N,
                   use = use,
                   n_factors_max = n_factors_max,
                   N_pop = N_pop,
                   N_samples = N_samples,
                   alpha = alpha,
                   cor_method = cor_method,
                   max_iter_CD = max_iter_CD,
                   n_fac_theor = n_fac_theor,
                   estimator = estimator,
                   # back-compat alias, as for the frozen P_type key
                   method = estimator,
                   gof = gof,
                   # `gof` above is what was requested; the Hull method reduces it to CAF
                   # whenever the PAF estimator is used, so record what actually ran
                   # alongside it. NA when HULL was not among the criteria, was skipped, or
                   # failed -- no Hull goodness-of-fit index was computed in those cases.
                   gof_used = if (is.null(outputs[["HULL"]])) NA_character_
                              else outputs[["HULL"]]$settings$gof,
                   eigen_type_HULL = eigen_type_HULL,
                   eigen_type_other = eigen_type_other,
                   n_factors = n_factors,
                   n_datasets = n_datasets,
                   percent = percent,
                   decision_rule = decision_rule,
                   ekc_type = ekc_type,
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

      cli::cli_rule(left = "Tests for the suitability of the data for factor analysis")
      cli::cli_text("")

      bart <- x$suitability$bartlett
      pval <- bart$p_value
      if (!is.null(pval) && !is.na(pval)) {
        stats_text <- paste0("\u03c7\u00b2(", bart$df, ") = ",
                             round(bart$chisq, 2), ", p", .screen_p_str(pval), ".")
        if (.screen_is_sig(pval)) {
          cli::cli_bullets(c("v" = paste0(
            "The Bartlett's test of sphericity was significant at an alpha level of .05: ",
            stats_text,
            " These data are probably suitable for factor analysis.")))
        } else {
          cli::cli_bullets(c("x" = paste0(
            "The Bartlett's test of sphericity was not significant at an alpha level of .05: ",
            stats_text,
            " These data are probably not suitable for factor analysis.")))
        }
      } else {
        cli::cli_bullets(c("!" =
          "The Bartlett's test of sphericity did not render a result."))
      }

      kmo <- x$suitability$kmo$KMO
      if (!is.null(kmo) && !is.na(kmo)) {
        # Kaiser's verbal bands, shared with efa_kmo() and efa_screen() so that the same
        # value is never reported at two different severities
        band <- .kmo_band(kmo)
        kmo_text <- paste0("The Kaiser-Meyer-Olkin criterion is ", band$label,
                           " (KMO = ", round(kmo, 3), "). These data are ",
                           band$suitability, " suitable for factor analysis.")
        cli::cli_bullets(stats::setNames(kmo_text, band$symbol))
      } else {
        cli::cli_bullets(c("!" =
          "The overall KMO value for your data is not available."))
      }

      cli::cli_text("")
    }

    is_visual <- vapply(names(x$outputs),
                        function(id) isTRUE(.retention_registry[[id]]$visual),
                        logical(1))

    cli::cli_rule(left = "Suggested number of factors")

    # one line answering "how many factors, and how much do the criteria disagree?"
    # before the per-criterion detail; interpolated as a variable so that cli
    # substitutes the value instead of re-parsing it as a template
    summary_line <- .retention_summary(x$outputs[!is_visual])
    if (!is.null(summary_line)) {
      cli::cli_text("")
      cli::cli_text("{summary_line}")
    }

    # criteria with a numeric suggestion: one group with one bullet per variant
    # (a variant the criterion could not determine shows as "not applicable"),
    # separated by a blank line for readability
    for (out in x$outputs[!is_visual]) {
      cli::cli_text("")
      cli::cli_text("{out$criterion[['label']]}")
      cli::cli_ul(.retention_bullets(out$results))
    }

    # visual criteria (e.g. the scree plot) are pointed to the plot
    for (out in x$outputs[is_visual]) {
      cli::cli_text("")
      label <- out$criterion[["label"]]
      cli::cli_bullets(c("i" =
        "{label} provides no numeric suggestion; inspect the plot."))
    }

    # criteria that were skipped or failed
    if (!is.null(x$not_run)) {
      cli::cli_text("")
      cli::cli_rule(left = "Criteria that could not be run")
      cli::cli_text("")
      for (id in names(x$not_run)) {
        reason <- x$not_run[[id]]
        cli::cli_bullets(c("!" = "{id}: {reason}"))
      }
    }

  })
}

# One-line summary of the criteria's suggestions for format.efa_retain(): how many
# suggestions how many criteria made, the range they span, and, when they disagree, the
# most common suggestion. The most common value is counted with ONE VOTE PER CRITERION --
# each criterion's own modal variant -- and never per variant, so that a criterion with
# several variants (e.g. the Hull method with three goodness-of-fit indices) does not
# outvote a single-variant one. A criterion whose variants have no single modal value is
# undecided and casts no vote. `outputs` are the criterion results with a numeric
# suggestion (the visual ones are excluded by the caller); returns NULL when none of them
# determined a number, or when a single suggestion would only restate the bullet below.
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
    return(paste0(counts, ", all suggesting ", nf[1], " ",
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

  spread <- paste0(counts, ", ranging from ", min(nf), " to ", max(nf), " factors")
  if (is.null(votes)) return(paste0(spread, "."))

  # a most common value only means something once some value was chosen by more than one
  # criterion. Where every deciding criterion picked a different number there is none, and
  # listing them all under "most common" would read as agreement where there is the
  # maximum possible disagreement, so the range is left to speak for itself.
  best <- modal(votes)
  if (best$n == 1L) return(paste0(spread, "."))

  paste0(spread, " (most common: ",
         .screen_and_list(as.character(best$values)), ").")

}

#' Plot method for efa_retain objects
#'
#' Plots every factor-retention criterion in the [efa_retain()] result that has
#' a plottable outcome (see [plot.efa_retention()]); criteria without a plot
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
