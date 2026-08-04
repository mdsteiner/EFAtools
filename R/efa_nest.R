#' Next eigenvalue sufficiency test (NEST)
#'
#' NEST uses many synthetic datasets to generate reference eigenvalues against
#' which to compare the empirical eigenvalues. This is similar to parallel
#' analysis, but other than parallel analysis, NEST does not just rely on
#' synthetic eigenvalues based on an identity matrix as null model.
#' It was introduced by Achim (2017), see also Brandenburg and Papenberg (2024) and
#' Caron (2025) for further simulation studies including NEST.
#'
#' @inheritParams efa_kgc
#' @param N numeric. The number of observations. Only needed if x is a correlation
#'  matrix.
#' @param alpha numeric. The alpha level to use (i.e., 1-alpha percentile of eigenvalues is used for reference values).
#' @param cor_method character. One of `"pearson"`, `"spearman"`, or `"kendall"`,
#'   passed to [stats::cor()]. `"poly"` and `"tetra"` are not supported because
#'   `NEST` compares the data against simulated continuous reference data.
#'  Default is  `"pearson"`.
#' @param n_datasets numeric. The number of datasets to simulate. Default is 1000.
#' @param estimate_control an [estimate_control()] object with the estimation settings for the
#'  [efa_fit()] reference-model fits. `NULL` (default) uses the [efa_fit()] defaults. The
#'  reference models are unrotated, so no rotation settings apply.
#' @param ... Additional arguments passed to [efa_fit()]. For example,
#' `estimator`, to change the estimator (default is "PAF"). PAF is more
#' robust, but it will take longer compared to the other estimators
#' available ("ML" and "ULS"). The estimation tuning knobs are not passed here; they
#' live in `estimate_control`, and the standard-error arguments (`se`, `b_boot`, `ci`,
#' `seed`) are not accepted because the reference-model fits are internal steps.
#'
#' @details NEST compares the first empirical eigenvalue against the first eigenvalues
#' of `n_dataset` synthetic datasets based on a null model  (i.e.,
#'  with uncorrelated variables; same as in parallel analysis, see [efa_parallel()]).
#'  The following eigenvalues are compared against synthetic datasets based on an EFA-model with one fewer factors
#'  than the position of the respective empirical eigenvalue. E.g, the second
#'  empirical eigenvalue is compared against synthetic data based on a one-factor
#'  model. In each comparison the \eqn{k}-th empirical eigenvalue is tested against
#'  the \eqn{k}-th largest eigenvalue of the synthetic datasets. The `alpha`-level
#'  defines against which percentile of the synthetic
#'  eigenvalue distribution to compare the empirical eigenvalues against, i.e., an
#'  alpha of .05 (the default) uses the 95th percentile as reference value.
#'
#'  The number of factors tested is capped at \eqn{\lfloor 0.8 \times p \rfloor}
#'  (with \eqn{p} the number of variables; Achim, 2017) and additionally limited so
#'  that the \eqn{(k - 1)}-factor reference model used at each step stays
#'  over-identified. If no empirical eigenvalue falls at or below its reference
#'  within this range, every tested factor is accepted and this capped number is
#'  returned.
#'
#'  Because each reference model carries the factors already retained, NEST does not
#'  lose accuracy for the strongly correlated factor structures where parallel
#'  analysis tends to under-extract, and it was among the more accurate criteria in
#'  the simulation studies of Brandenburg and Papenberg (2024) and Caron (2025). The
#'  price is runtime: a fresh set of `n_datasets` reference datasets is drawn and
#'  eigen-decomposed at every candidate factor count, which makes NEST one of the
#'  slowest criteria available here.
#'
#'  The reference models are fitted without inequality constraints. A Heywood case in
#'  one of them leaves no unique variance to simulate the reference data from, so
#'  NEST aborts rather than continuing from an inadmissible reference.
#'
#'  The reference eigenvalues are obtained from simulated data, so the suggested number
#'  of factors varies slightly from run to run. Call [set.seed()] beforehand to make a
#'  run reproducible.
#'
#'  For details on the method, including simulation studies, see Achim (2017),
#'  Brandenburg and Papenberg (2024), and Caron (2025).
#'
#'  The `efa_nest` function can also be called together with other factor
#'   retention criteria in the [efa_retain()] function.
#'
#'
#' @returns An object of class `efa_retention` (see [print.efa_retention()] and
#'   [plot.efa_retention()] for the print and plot methods). Its main fields are:
#' \item{n_factors}{A named numeric vector (`"NEST"`) with the suggested number of
#'   factors according to the NEST procedure.}
#' \item{results}{A list with a single record holding the empirical eigenvalues
#'   and the reference eigenvalues. Only the positions the search actually tested
#'   carry a reference value; beyond the position at which it stopped the series is
#'   `NA`.}
#' \item{settings}{A list of control settings used.}
#'
#' @source Achim, A. (2017). Testing the number of required dimensions in exploratory factor analysis. The Quantitative Methods for Psychology, 13(1), 64–74. https://doi.org/10.20982/tqmp.13.1.p064
#' @source Brandenburg, N., & Papenberg, M. (2024). Reassessment of innovative methods to determine the number of factors: A simulation-based comparison of exploratory graph analysis and Next Eigenvalue Sufficiency Test. Psychological Methods, 29(1), 21–47. https://doi.org/10.1037/met0000527
#' @source Caron, P.-O. (2025). A Comparison of the Next Eigenvalue Sufficiency Test to Other Stopping Rules for the Number of Factors in Factor Analysis.
#' Educational and Psychological Measurement, Online-first publication. https://doi.org/10.1177/00131644241308528
#'
#' @family factor retention criteria
#'
#' @seealso [efa_retain()] as a wrapper function for this and the other factor
#'   retention criteria.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # with correlation matrix
#' efa_nest(test_models$baseline$cormat, N = 500)
#'
#' # with raw data
#' efa_nest(GRiPS_raw)
#' }
efa_nest <- function(x, N = NA,
                 alpha = .05,
                 use = c("pairwise.complete.obs", "all.obs",
                         "complete.obs", "everything",
                         "na.or.complete"),
                 cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                 n_datasets = 1000,
                 estimate_control = NULL,
                 ...) {


  # Perform argument checks
  .reject_flat_knobs(...names(), fn = "efa_nest")
  .reject_unknown_fit_dots(...names(), fn = "efa_nest", unrotated = TRUE)
  .reject_rotation_dots(list(...), fn = "efa_nest")
  .assert_cor_input(x)

  checkmate::assert_count(N, na.ok = TRUE)
  checkmate::assert_number(alpha, lower = 0, upper = 1)
  use <- .match_arg_ci(use)
  cor_method <- .match_arg_ci(cor_method)
  .assert_estimate_control(estimate_control)
  .reject_poly_reference(cor_method, "efa_nest")
  checkmate::assert_count(n_datasets, na.ok = FALSE,
                          positive = TRUE)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, N = N, use = use, cor_method = cor_method,
                             N_policy = "required")
  R <- prep$R
  N <- prep$N

  emp_eigen <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  nvar <- ncol(R)

  # Number of factors to test. Achim's (2017) rule caps this at floor(.8 * nvar);
  # it is additionally bounded so the (nf - 1)-factor reference model fit at each
  # step stays over-identified (df > 0). The reference data for the nf-th empirical
  # eigenvalue come from an (nf - 1)-factor EFA, whose loadings are only
  # trustworthy when that model has positive degrees of freedom
  # (df = ((nvar - k)^2 - (nvar + k)) / 2 for k factors); the largest reference
  # factor count with df > 0 caps the search at one factor beyond it.
  ref_k <- seq_len(nvar) - 1L
  ref_df <- .efa_df(nvar, ref_k)
  max_fac <- min(floor(.8 * nvar), max(ref_k[ref_df > 0]) + 1L)

  references <- vector(mode = "double", length = max_fac)

  # `stopped` records whether the search ended because an empirical eigenvalue
  # fell at/below its reference (i.e. a factor was rejected). If none does so
  # within the tested range, every tested factor was accepted and the last tested
  # model is retained (see the boundary handling after the loop).
  stopped <- FALSE
  for (nf in 1:max_fac) {

    if (nf == 1) {

      # For the first factor the reference is the null (identity) model: no loadings
      # and unit uniquenesses.
      Lambda <- matrix(numeric(0), nvar, 0)
      Psi <- rep(1, nvar)

    } else {


      # For further factors, use a model with nf - 1 factors as reference. The fit is an
      # internal step run once per tested factor; suppress its warnings so a forwarded
      # estimator (e.g. estimator = "ML") cannot turn per-model Heywood or non-convergence
      # warnings into one warning per loop iteration. A Heywood case is handled
      # explicitly below.
      mm <- suppressWarnings(efa_fit(R, N = N, n_factors = nf - 1,
                                     estimate_control = estimate_control, ...))
      Lambda <- mm$unrot_loadings

      # A Heywood case in the reference model (communality at or above 1) makes
      # the uniqueness negative, so the reference data cannot be simulated; abort
      # with an actionable message instead of passing a negative variance to the
      # C++ simulator.
      Psi <- 1 - mm$h2
      if (any(!is.finite(Psi) | Psi < 0)) {
        cli::cli_abort(
          c("A Heywood case occurred in the {nf - 1}-factor reference model used by NEST.",
            "x" = "A communality at or above 1 leaves no unique variance to simulate the reference data from.",
            "i" = "Try a different {.arg estimator} (e.g. {.val ML} or {.val ULS}) or fewer indicators."),
          class = "efa_nest_heywood"
        )
      }

    }

    # Simulate the reference datasets and return the nf-th largest eigenvalue of each.
    # The shared kernel draws via the factor-score rule X = F*L' + e from the reference
    # model (loadings Lambda, uniquenesses Psi), which is faster than drawing from the
    # model-implied correlation matrix.
    ref_values <- .simulate_cfm_eigen(nf, N, Lambda, Psi, n_datasets)

    references[nf] <- stats::quantile(ref_values, probs = 1 - alpha)
    if (emp_eigen[nf] <= references[nf]) {
      stopped <- TRUE
      break
    }


  }

  # On a rejection the retained model is the one with one fewer factor than the
  # rejected position. At the no-stop boundary (no rejection within the tested
  # range) every tested factor was accepted, so the last tested model
  # (`nf == max_fac`) is retained rather than the just-past index.
  n_factors <- if (stopped) nf - 1L else nf
  # Only the tested positions have a reference; pad the rest with NA so the series
  # aligns with the empirical eigenvalues (the shared eigenvalue plot binds the two
  # into one data frame and drops the NAs).
  references <- c(references[seq_len(nf)],
                  rep(NA_real_, length(emp_eigen) - nf))

  results <- list(list(
    name = "NEST",
    label = "Suggested number of factors",
    n_factors = n_factors,
    plot_type = "eigen",
    x = seq_along(emp_eigen),
    y = emp_eigen,
    reference = references,
    # mark the retained solution in the plot, as the other eigenvalue-based
    # criteria do; there is nothing to mark when no factor is retained
    highlight = if (n_factors >= 1) n_factors else NULL
  ))

  out <- .new_efa_retention(
    "NEST",
    results = results,
    settings = list(alpha = alpha, N = N, n_datasets = n_datasets,
                    use = use, cor_method = cor_method)
  )

  out

}
