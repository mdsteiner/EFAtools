#' Various Factor Retention Criteria
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
#' (see details). By default, a subset of often used, well-performing methods are performed.
#' @param suitability logical. Whether the data should be checked for suitability
#' for factor analysis using the Bartlett's test of sphericity and the
#' Kaiser-Meyer-Olkin criterion (see details). Default is `TRUE`.
#' @param N  numeric. The number of observations. Only needed if x is a
#' correlation matrix.
#' @param use character. Passed to [stats::cor()] if raw
#' data is given as input. Default is `"pairwise.complete.obs"`.
#' @param cor_method character. Correlation computed from raw data: `"pearson"`,
#'   `"spearman"`, or `"kendall"` (passed to [stats::cor()]), or `"poly"` /
#'   `"tetra"` for polychoric / tetrachoric correlations (a two-step estimator
#'   with no empty-cell continuity correction). Criteria that compare against
#'   simulated continuous data (`CD`, `PARALLEL`, `NEST`, `HULL`) do not support
#'   `"poly"` / `"tetra"` and are skipped in that case.
#' Default is  `"pearson"`.
#' @param n_factors_max numeric. Passed to [CD()]. The maximum number
#' of factors to test against.
#' Larger numbers will increase the duration the procedure takes, but test more
#' possible solutions. Maximum possible is number of variables / 2. Default is
#' NA. If not specified, number of variables / 2 is used.
#' @param N_pop numeric. Passed to [CD()]. Size of finite populations
#' of comparison data. Default is 10000.
#' @param N_samples numeric. Passed to [CD()]. Number of samples drawn
#'  from each population. Default is 500.
#' @param alpha numeric. Passed to [CD()]. The alpha level used to test
#'  the significance of the improvement added by an additional factor.
#'  Default is .30.
#' @param max_iter_CD numeric. Passed to [CD()]. The maximum number of
#'  iterations to perform after which the iterative PAF procedure is halted.
#'   Default is 50.
#' @param n_fac_theor numeric. Passed to [HULL()]. Theoretical number
#'  of factors to retain. The maximum of this number and the number of factors
#'  suggested by [PARALLEL] plus one will be used in the Hull method.
#' @param method character. Passed to [EFA()] in [HULL()],
#' [KGC()], [SCREE()], [PARALLEL()], and [NEST()]. The
#' estimation method to use. One of  `"PAF"`, `"ULS"`, or  `"ML"`,
#' for principal axis factoring, unweighted least squares, and maximum
#' likelihood, respectively.
#' @param gof character. Passed to [HULL()]. The goodness of fit index
#' to use. Either `"CAF"`, `"CFI"`, or `"RMSEA"`, or any
#' combination of them. If `method = "PAF"` is used, only
#' the CAF can be used as goodness of fit index. For details on the CAF, see
#' Lorenzo-Seva, Timmerman, and Kiers (2011).
#' @param eigen_type_HULL character. Passed to  [PARALLEL()] in
#' [HULL()]. On what the
#' eigenvalues should be found in the parallel analysis. Can be one of
#' `"SMC"`, `"PCA"`, or `"EFA"`. If using  `"SMC"` (default),
#' the diagonal of the correlation matrices is
#' replaced by the squared multiple correlations (SMCs) of the indicators. If
#' using  `"PCA"`, the diagonal values of the correlation
#' matrices are left to be 1. If using  `"EFA"`, eigenvalues are found on the
#' correlation  matrices with the final communalities of an EFA solution as
#' diagonal.
#' @param eigen_type_other character. Passed to [KGC()],
#' [SCREE()], and [PARALLEL()]. The same as eigen_type_HULL,
#' but multiple inputs are possible here (any combination of `"PCA"`, `"SMC"`,
#' and `"EFA"`). Default is `"SMC"`.
#' @param n_factors numeric. Passed to [PARALLEL()] (also within
#' [HULL()]), [KGC()], and [SCREE()]. Number of
#' factors to extract if `"EFA"` is included in `eigen_type_HULL` or
#'  `eigen_type_other`. Default is 1.
#' @param n_datasets numeric. Passed to [PARALLEL()] (also within
#' [HULL()]). The number of datasets to simulate. Default is 1000.
#' @param percent numeric. Passed to [PARALLEL()] (also within
#' [HULL()]). The percentile to take from the simulated eigenvalues.
#'  Default is 95.
#' @param decision_rule character. Passed to [PARALLEL()] (also within
#'  [HULL()]). Which rule to use to determine the number of
#'  factors to retain. Default is `"means"`, which will use the average
#'  simulated eigenvalues. `"percentile"`, uses the percentiles specified
#'  in percent. `"crawford"` uses the 95th percentile for the first factor
#'  and the mean afterwards (based on Crawford et al, 2010).
#' @param ekc_type character. Passed to the `type` argument of [EKC()].
#'   Either `"BvA2017"` for the original implementation by Braeken and van Assen
#'   (2017), or `"AM2019"` for the adapted implementation by Auerswald and Moshagen
#'   (2019).
#' @param n_datasets_nest numeric. The number of datasets to simulate in [NEST()]. Default is 1000.
#' @param alpha_nest numeric. The alpha level to use in [NEST()] (i.e., 1-alpha percentile of eigenvalues is used for reference values).
#' @param show_progress logical. Whether a progress bar should be shown in the
#'   console. Default is FALSE.
#' @param ... Further arguments passed to [EFA()] in
#' [PARALLEL()] (also within [HULL()]) and [KGC()].
#'
#' @details
#' By default, the entered data are checked for suitability for factor analysis
#' using the following methods (see respective documentations for details):
#' \itemize{
#' \item{Bartlett's test of sphericity (see [BARTLETT()])}
#' \item{Kaiser-Meyer-Olkin criterion (see [EFAtools::KMO()])}}
#'
#' The available factor retention criteria are the following (see respective
#'  documentations for details):
#'  \itemize{
#' \item{Comparison data (see [CD()])}
#' \item{Empirical Kaiser criterion (see [EKC()])}
#' \item{Hull method (see [HULL()])}
#' \item{Kaiser-Guttman criterion (see [KGC()])}
#' \item{Parallel analysis (see [PARALLEL()])}
#' \item{Next Eigenvalue Sufficiency Test, NEST (see [NEST()])}
#' \item{Scree plot (see [SCREE()])}
#' \item{Sequential chi-square model tests, RMSEA lower bound, and AIC
#' (see [SMT()])}
#' }
#'
#' @returns A list of class N_FACTORS containing
#' \item{suitability}{A list with the results from [BARTLETT()] and
#'   [EFAtools::KMO()] (`bartlett` and `kmo`), or `NULL` if
#'   `suitability = FALSE`.}
#' \item{outputs}{A named list with one `efa_retention` object per factor
#'   retention criterion that was run (see, e.g., [EKC()]).}
#' \item{n_factors}{A named numeric vector with the suggested number of factors
#'   per criterion and, where a criterion has several variants, per variant
#'   (e.g. `EKC_BvA2017` or `PARALLEL_SMC`). Criteria without a numeric
#'   suggestion (the scree plot) are not included.}
#' \item{not_run}{A named character vector with the criteria that were skipped
#'   or failed and the reason, or `NULL` if all requested criteria ran.}
#' \item{settings}{A list of the settings used.}
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Default criteria, with correlation matrix and fit method "ML" (where needed)
#' # This will throw a warning for CD, as no raw data were specified
#' nfac_all <- N_FACTORS(test_models$baseline$cormat, N = 500, method = "ML")
#'
#' # The same as above, but without "CD"
#' nfac_wo_CD <- N_FACTORS(test_models$baseline$cormat, criteria = c("EKC",
#'                         "HULL", "PARALLEL", "NEST"), N = 500,
#'                         method = "ML")
#'
#' # Use PAF instead of ML (this will take longer). For this, gof has
#' # to be set to "CAF" for the Hull method.
#' nfac_PAF <- N_FACTORS(test_models$baseline$cormat, criteria = c("EKC",
#'                         "HULL", "PARALLEL", "NEST"), N = 500,
#'                       method = "PAF", gof = "CAF")
#'
#' # Do KGC and PARALLEL with only "PCA" type of eigenvalues
#' nfac_PCA <- N_FACTORS(test_models$baseline$cormat, criteria = c("EKC",
#'                       "HULL", "PARALLEL", "NEST"), N = 500,
#'                       method = "ML", eigen_type_other = "PCA")
#'
#' # Use raw data, such that CD can also be performed
#' nfac_raw <- N_FACTORS(GRiPS_raw, method = "ML")
#'}
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
                      ...){

  # Perform argument checks
  .assert_cor_input(x)

  ## Perform argument checks and prepare input
  criteria <- match.arg(criteria, several.ok = TRUE,
                        choices = c("CD", "EKC", "HULL", "KGC", "PARALLEL",
                                    "SCREE", "SMT", "NEST", "MAP"))
  suitability <- checkmate::assert_flag(suitability)
  eigen_type_HULL <- match.arg(eigen_type_HULL)
  eigen_type_other <- match.arg(eigen_type_other, several.ok = TRUE,
                                choices = c("PCA", "SMC", "EFA"))
  gof <- match.arg(gof, several.ok = TRUE, choices = c("CAF", "CFI", "RMSEA"))
  cor_method <- match.arg(cor_method)
  use <- match.arg(use)
  method <- match.arg(method)
  decision_rule <- match.arg(decision_rule)
  ekc_type <- match.arg(ekc_type, c("BvA2017", "AM2019"), several.ok = TRUE)
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
      BARTLETT(R, N = N, use = use, cor_method = cor_method)
    }
    suitability_out <- list(
      bartlett = bartlett_out,
      kmo = KMO(R, use = use, cor_method = cor_method)
    )
  }

  ## Factor retention criteria, driven by the registry (run in registry order)
  ctl <- list(N = N, use = use, cor_method = cor_method,
              n_factors_max = n_factors_max, N_pop = N_pop,
              N_samples = N_samples, alpha = alpha, max_iter_CD = max_iter_CD,
              n_fac_theor = n_fac_theor, method = method, gof = gof,
              eigen_type_HULL = eigen_type_HULL,
              eigen_type_other = eigen_type_other, n_factors = n_factors,
              n_datasets = n_datasets, percent = percent,
              decision_rule = decision_rule, ekc_type = ekc_type,
              n_datasets_nest = n_datasets_nest, alpha_nest = alpha_nest,
              dots = list(...))

  run <- intersect(names(.retention_registry), criteria)
  outputs <- list()
  not_run <- character(0)   # id -> reason it was skipped or failed

  for (id in run) {

    if (isTRUE(show_progress)) {
      cli::cli_progress_step("Running {id}")
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

    # Criteria that compare against continuous reference data cannot use
    # polychoric/tetrachoric correlations; skip them with a clear note rather than
    # letting them fail and be reported as a generic error.
    if (.is_poly_cor(cor_method) && isFALSE(entry$poly_ok)) {
      cli::cli_warn(
        c("{.val {id}} does not support {.code cor_method = {.val {cor_method}}}.",
          "i" = "Its reference data are continuous; skipping {.val {id}}."),
        class = "efa_criterion_skipped"
      )
      not_run[[id]] <- paste0("does not support cor_method = \"", cor_method, "\"")
      next
    }

    # a failing criterion is excluded with a warning; the others still run
    out_id <- try(entry$fun(if (isTRUE(entry$needs_raw)) x else R, ctl),
                  silent = TRUE)

    if (inherits(out_id, "try-error")) {
      reason <- conditionMessage(attr(out_id, "condition"))
      # keep the headline only; cli error bodies span several lines
      reason <- strsplit(reason, "\n", fixed = TRUE)[[1]][1]
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
      class = "efa_no_criteria"
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
                   method = method,
                   gof = gof,
                   eigen_type_HULL = eigen_type_HULL,
                   eigen_type_other = eigen_type_other,
                   n_factors = n_factors,
                   n_datasets = n_datasets,
                   percent = percent,
                   decision_rule = decision_rule,
                   ekc_type = ekc_type,
                   n_datasets_nest = n_datasets_nest,
                   alpha_nest = alpha_nest)

  # Aggregate the suggested numbers of factors as "<id>" or "<id>_<variant>".
  # Visual criteria (the scree plot) make no numeric suggestion and are omitted;
  # for the others NA suggestions are kept (named), so a criterion that ran but
  # could not determine a number stays visible.
  n_factors <- unlist(lapply(names(outputs), function(id) {
    if (isTRUE(.retention_registry[[id]]$visual)) return(NULL)
    nf <- outputs[[id]]$n_factors
    names(nf) <- ifelse(names(nf) == id, id, paste(id, names(nf), sep = "_"))
    nf
  }))
  if (is.null(n_factors)) n_factors <- numeric(0)

  output <- list(suitability = suitability_out,
                 outputs = outputs,
                 n_factors = n_factors,
                 not_run = if (length(not_run) > 0) not_run else NULL,
                 settings = settings)

  class(output) <- "N_FACTORS"

  return(output)

}

#' Print method for N_FACTORS objects
#'
#' @param x an object of class N_FACTORS, returned by [N_FACTORS()].
#' @param ... not used.
#'
#' @returns `print()` returns its argument `x` invisibly; it is
#'   `cat(format(x), sep = "\n")`.
#'
#' @export
#' @method print N_FACTORS
#'
#' @examples
#' \donttest{
#' N_FACTORS(test_models$baseline$cormat, criteria = c("EKC", "SMT"), N = 500)
#' }
print.N_FACTORS <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' Format method for N_FACTORS objects
#'
#' @param x an object of class N_FACTORS, returned by [N_FACTORS()].
#' @param ... not used.
#'
#' @returns A character vector with the report lines (styled to the active
#'   console theme; plain when colours are disabled).
#'
#' @export
#' @method format N_FACTORS
#'
#' @examples
#' \donttest{
#' nf <- N_FACTORS(test_models$baseline$cormat, criteria = c("EKC", "SMT"),
#'                 N = 500)
#' writeLines(format(nf))
#' }
format.N_FACTORS <- function(x, ...) {
  cli::cli_format_method({

    if (!is.null(x$suitability)) {

      cli::cli_rule(left = "Tests for the suitability of the data for factor analysis")
      cli::cli_text("")

      bart <- x$suitability$bartlett
      pval <- bart$p_value
      if (!is.null(pval) && !is.na(pval)) {
        p_text <- if (pval < .001) "p < .001" else paste0("p = ", round(pval, 3))
        stats_text <- paste0("\u03c7\u00b2(", bart$df, ") = ",
                             round(bart$chisq, 2), ", ", p_text, ".")
        if (pval < .05) {
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
        # Kaiser's verbal labels for the KMO ranges
        kmo_label <- if (kmo >= .9) {
          "marvellous"
        } else if (kmo >= .8) {
          "meritorious"
        } else if (kmo >= .7) {
          "middling"
        } else if (kmo >= .6) {
          "mediocre"
        } else if (kmo >= .5) {
          "miserable"
        } else {
          "unacceptable"
        }
        verdict <- if (kmo < .5) {
          "These data are not suitable for factor analysis."
        } else if (kmo < .6) {
          "These data are hardly suitable for factor analysis."
        } else {
          "These data are probably suitable for factor analysis."
        }
        kmo_text <- paste0("The Kaiser-Meyer-Olkin criterion is ", kmo_label,
                           " (KMO = ", round(kmo, 3), "). ", verdict)
        kmo_symbol <- if (kmo >= .7) "v" else if (kmo >= .6) "!" else "x"
        cli::cli_bullets(stats::setNames(kmo_text, kmo_symbol))
      } else {
        cli::cli_bullets(c("!" =
          "The overall KMO value for your data is not available."))
      }

      cli::cli_text("")
    }

    cli::cli_rule(left = "Number of factors suggested by the factor retention criteria")

    is_visual <- vapply(names(x$outputs),
                        function(id) isTRUE(.retention_registry[[id]]$visual),
                        logical(1))

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

#' Plot method for N_FACTORS objects
#'
#' Plots every factor-retention criterion in the [N_FACTORS()] result that has
#' a plottable outcome (see [plot.efa_retention()]); criteria without a plot
#' (e.g. [MAP()] or [SMT()]) are skipped.
#'
#' @param x an object of class N_FACTORS, returned by [N_FACTORS()].
#' @param ... not used.
#'
#' @returns A named list of [ggplot2::ggplot] objects, one per criterion with a
#'   plottable result, or invisibly `NULL` if there is none.
#'
#' @export
#' @method plot N_FACTORS
#'
#' @examples
#' \donttest{
#' nf <- N_FACTORS(test_models$baseline$cormat, criteria = c("EKC", "SMT"),
#'                 N = 500)
#' plot(nf)
#' }
plot.N_FACTORS <- function(x, ...) {

  # plot.efa_retention returns NULL (with a message) for criteria with no plot
  # (e.g. MAP/SMT) or a degenerate one (e.g. CD suggesting 0 factors); drop those
  plots <- lapply(x$outputs, function(o) suppressMessages(plot(o)))
  plots <- Filter(Negate(is.null), plots)

  if (length(plots) == 0) {
    cli::cli_inform("No plot is available for the criteria that were run.")
    return(invisible(NULL))
  }

  plots

}
