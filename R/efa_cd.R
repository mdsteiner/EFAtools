#' Comparison data
#'
#' Factor retention method introduced by Ruscio and Roche (2012). The code was
#' adapted from the CD code by Auerswald and Moshagen (2019) available at
#' <https://osf.io/x5cz2/?view_only=d03efba1fd0f4c849a87db82e6705668>
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data.
#' @param n_factors_max numeric. The maximum number of factors to test against.
#'  Larger numbers will increase the duration the procedure takes, but test more
#'  possible solutions. If left NA (default) the maximum number of factors for
#'  which the model is still over-identified (df > 0) is used.
#' @param N_pop numeric. Size of finite populations of comparison data. Default
#'  is 10000.
#' @param N_samples numeric. Number of samples drawn from each population.
#'  Default is 500.
#' @param alpha numeric. The alpha level used to test the significance of the
#'  improvement added by an additional factor. Default is .30.
#' @param cor_method character. One of `"pearson"`, `"spearman"`, or `"kendall"`,
#'   passed to [stats::cor()]. `"poly"` and `"tetra"` are not supported because
#'   `CD` compares the data against simulated continuous reference data.
#' Default is "pearson".
#' @param max_iter numeric. The maximum number of iterations after which the
#'  iterative PAF procedure inside the comparison-data generation is halted; it
#'  does not cap an EFA of `x`. Default is 50.
#'
#' @details Comparison data (CD) extends parallel analysis by reproducing the
#' observed correlation matrix rather than generating random data: datasets with a
#' known factor structure are generated with an increasing number of factors, and
#' the smallest number for which adding a further factor no longer significantly
#' improves the reproduction of the observed eigenvalues is retained (Ruscio &
#' Roche, 2012).
#'
#' Because it reproduces the observed correlation structure instead of a null model,
#' CD was among the more accurate criteria across a broad range of conditions in
#' Ruscio and Roche (2012). It is, however, the only criterion in this family that
#' requires raw data, and by some margin the most computationally intensive one: a
#' finite population of `N_pop` cases is generated and `N_samples` samples are drawn
#' from it at every candidate factor count. It is therefore a good choice when the
#' raw data are at hand and the runtime is acceptable, and a poor one for a quick
#' look at a correlation matrix.
#'
#' Note that if the data contains missing values, these will be removed for the
#' comparison data procedure using [`stats::na.omit()`][stats::na.fail]. If
#' missing data should be treated differently, e.g., by imputation, do this outside
#' `efa_cd()` and then pass the complete data.
#'
#' @returns An object of class `efa_retention` (see [print.efa_retention()] and
#'   [plot.efa_retention()] for the print and plot methods). Its main fields are:
#' \item{n_factors}{A named numeric vector (`"CD"`) with the suggested number of
#'   factors according to comparison data results.}
#' \item{results}{A list with a single record holding the mean RMSE between the
#'   eigenvalues of the generated and the entered data per number of factors
#'   (used for the plot) and, in `rmse_eigenvalues`, the per-sample RMSE matrix
#'   (rows are samples, columns are factor counts; columns beyond the last tested
#'   factor count are left as zero).}
#' \item{settings}{A list of the settings used.}
#'
#' @source Auerswald, M., & Moshagen, M. (2019). How to determine the number of
#' factors to retain in exploratory factor analysis: A comparison of extraction
#' methods under realistic conditions. Psychological Methods, 24(4), 468–491.
#' https://doi.org/10.1037/met0000200
#'
#' @source Ruscio, J., & Roche, B. (2012). Determining the number of factors to
#' retain in an exploratory factor analysis using comparison data of known
#' factorial structure. Psychological Assessment, 24, 282–292.
#' doi: 10.1037/a0025697
#'
#' @family factor retention criteria
#'
#' @seealso [efa_retain()] as a wrapper function for this and the other factor
#'   retention criteria.
#'
#' @export
#'
#'
#' @examples
#' \donttest{
#' # determine n factors of the GRiPS
#' efa_cd(GRiPS_raw)
#'
#' # determine n factors of the DOSPERT risk subscale
#' efa_cd(DOSPERT_raw)
#'}
efa_cd <- function(x, n_factors_max = NA, N_pop = 10000, N_samples = 500, alpha = .30,
               cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
               max_iter = 50) {

  # Perform argument checks
  .assert_cor_input(x, raw_only = TRUE)

  if (.is_cormat(x)) {
    cli::cli_abort("{.arg x} is a correlation matrix, but CD only works with raw data.",
                   class = "efa_cd_needs_raw")
  }

  if (inherits(x, c("tbl_df", "tbl"))) {
    x <- as.data.frame(x)
  }

  cor_method <- match.arg(cor_method)
  .reject_poly_reference(cor_method, "efa_cd")

  checkmate::assert_count(n_factors_max, na.ok = TRUE)
  checkmate::assert_count(N_pop)
  checkmate::assert_count(N_samples)
  checkmate::assert_number(alpha, lower = 0, upper = 1)
  checkmate::assert_count(max_iter)

  if (any(is.na(x))) {
    n_row_complete <- nrow(x)
    x <- stats::na.omit(x)
    n_row_new <- nrow(x)
    n_rows_removed <- n_row_complete - n_row_new

    cli::cli_warn(
      c("The data contained missing values, removed with {.fn stats::na.omit}.",
        "i" = "{n_rows_removed} row{?s} removed."),
      class = "efa_cd_missing_removed"
    )
  }
  n_cases <- nrow(x)
  k <- ncol(x)

  m_possible <- .det_max_factors(k)

  # Comparison data needs at least one over-identified model (df > 0) to compare
  # against; with too few indicators no such model exists and the search range is
  # empty, so abort rather than silently returning zero factors.
  if (m_possible < 1) {
    cli::cli_abort(
      c("Comparison data cannot be run because no factor model with {k} indicator{?s} is over-identified.",
        "i" = "Provide more indicators."),
      class = "efa_cd_min_indicators"
    )
  }

  if (is.na(n_factors_max) || n_factors_max > m_possible) {

    if (!is.na(n_factors_max) & n_factors_max > m_possible) {
      cli::cli_warn(
        c("{.arg n_factors_max} was set to {n_factors_max}, but at most {m_possible} factor{?s} can be extracted.",
          "i" = "Setting {.arg n_factors_max} to {m_possible}."),
        class = "efa_cd_max_factors"
      )
    }

    n_factors_max <- m_possible

  }

  # Create correlation matrix (x has no missing values at this point: incomplete
  # rows were removed above)
  R <- stats::cor(x, method = cor_method)

  eigvals_real <- eigen(R, symmetric = TRUE, only.values = TRUE)$values

  # initialize objects for iterative procedures
  RMSE_eigvals <- matrix(0, nrow = N_samples, ncol = n_factors_max)
  sig <- TRUE
  n_factors <- 1



  while (n_factors <= n_factors_max && isTRUE(sig)) {

    pop <- .simulate_cfm_empirical(R = R, x = x, n_factors = n_factors, N = N_pop,
                                   cor_method = cor_method, max_iter = max_iter)

    for (j in 1:N_samples) {

      samp <- pop[sample(1:N_pop, size = n_cases, replace = TRUE),]
      R_samp <- stats::cor(samp, method = cor_method)
      eigvals_samp <- eigen(R_samp, symmetric = TRUE, only.values = TRUE)$values
      RMSE_eigvals[j,n_factors] <- sqrt(sum((eigvals_samp - eigvals_real) *
                                         (eigvals_samp - eigvals_real)) / k)
    }

    if (n_factors > 1) {

      sig <- (stats::wilcox.test(RMSE_eigvals[,n_factors],
                          RMSE_eigvals[,(n_factors - 1)], "less")$p.value < alpha)
    }
    if (isTRUE(sig)) {
      n_factors <- n_factors + 1
    }

  }

  n_factors <- n_factors - 1

  settings <- list(
    n_factors_max = n_factors_max,
    N_pop = N_pop,
    N_samples = N_samples,
    alpha = alpha,
    cor_method = cor_method,
    max_iter = max_iter
  )

  # The RMSE curve is plotted over every tested factor count, including the first
  # count whose lack of significant improvement stopped the search (one beyond the
  # retained number), so the flattening that drives the decision stays visible.
  # When the search instead ran to n_factors_max, all columns are populated.
  n_tested <- min(n_factors + 1, n_factors_max)

  # single record: the mean RMSE curve over candidate factor counts (the full
  # per-sample RMSE matrix is kept in rmse_eigenvalues)
  results <- list(list(
    name = "CD",
    label = "Suggested number of factors",
    n_factors = n_factors,
    plot_type = "eigen",
    x = seq_len(n_tested),
    y = colMeans(RMSE_eigvals)[seq_len(n_tested)],
    highlight = if (n_factors >= 1) n_factors else NULL,
    y_label = "RMSE eigenvalues",
    rmse_eigenvalues = RMSE_eigvals
  ))

  out <- .new_efa_retention(
    "CD",
    results = results,
    settings = settings
  )

  return(out)

}
