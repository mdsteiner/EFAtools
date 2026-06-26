#' Comparison Data
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
#' @param max_iter numeric. The maximum number of iterations to perform after
#'  which the iterative PAF procedure is halted. Default is 50.
#'
#' @details "Parallel analysis (PA) is an effective stopping rule that compares
#' the eigenvalues of randomly generated data with those for the actual data.
#' PA takes into account sampling error, and at present it is widely considered
#' the best available method. We introduce a variant of PA that goes even further
#' by reproducing the observed correlation matrix rather than generating random
#' data. Comparison data (CD) with known factorial structure are first generated
#' using 1 factor, and then the number of factors is increased until the
#' reproduction of the observed eigenvalues fails to improve significantly"
#' (Ruscio & Roche, 2012, p. 282).
#'
#' The CD implementation here is based on the code by Ruscio and Roche (2012), but
#' is slightly adapted to increase speed by performing the principal axis factoring
#' using a C++ based function.
#'
#' Note that if the data contains missing values, these will be removed for the
#' comparison data procedure using [`stats::na.omit()`][stats::na.fail]. If
#' missing data should be treated differently, e.g., by imputation, do this outside
#' `CD` and then pass the complete data.
#'
#' The `CD` function can also be called together with other factor retention
#' criteria in the [N_FACTORS()] function.
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
#' @seealso [N_FACTORS()] as a wrapper function for this and the other factor
#'   retention criteria.
#'
#' @export
#'
#'
#' @examples
#' \donttest{
#' # determine n factors of the GRiPS
#' CD(GRiPS_raw)
#'
#' # determine n factors of the DOSPERT risk subscale
#' CD(DOSPERT_raw)
#'}
CD <- function(x, n_factors_max = NA, N_pop = 10000, N_samples = 500, alpha = .30,
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
  .reject_poly_reference(cor_method, "CD")

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

    pop <- .gen_data(x, cor_method = cor_method, n_factors, N_pop,
                     max_iter = max_iter)

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


.gen_data <- function(x, cor_method, n_factors, N, max_trials = 5,
                      initial_multiplier = 1, max_iter = 100) {
  # Steps refer to description in the following article:
  # Ruscio, J., & Kaczetow, W. (2008). Simulating multivariate nonnormal data using an iterative algorithm.
  # Multivariate Behavioral Research, 43(3), 355-381.

  # Initialize variables and (if applicable) set random number seed (step 1) -------------------------------------

  k <- ncol(x)
  sim_dat <- matrix(0, nrow = N, ncol = k)         # Matrix to store the simulated data
  dists <- matrix(0, nrow = N, ncol = k)   # Matrix to store each variable's score distribution
  best_RMSE <- 1                                   # Lowest RMSE correlation
  t_no_impr <- 0                  # Trial counter

  # Generate distribution for each variable (step 2) -------------------------------------------------------------

  for (i in 1:k) {
    dists[,i] <- sort(sample(x[,i], size = N, replace = TRUE))
  }


  # Calculate and store a copy of the target correlation matrix (step 3) -----------------------------------------

  R <- stats::cor(x, method = cor_method)
  R_inter <- R
  # Seed the best-so-far matrices so the no-improvement branch below is safe even
  # if the first iteration does not lower the initial RMSE (best_RMSE = 1).
  R_best <- R_inter
  res_best <- matrix(0, nrow = k, ncol = k)

  # Generate random normal data for shared and unique components, initialize factor loadings (steps 5, 6) --------

  shared_comp <- matrix(stats::rnorm(N * n_factors, 0, 1), nrow = N,
                        ncol = n_factors)
  unique_comp <- matrix(stats::rnorm(N * k, 0, 1), nrow = N, ncol = k)
  shared_load <- matrix(0, nrow = k, ncol = n_factors)
  unique_load <- matrix(0, nrow = k, ncol = 1)

  # Begin loop that ends when specified number of iterations pass without improvement in RMSE correlation --------

  while (t_no_impr < max_trials) {

    # Calculate factor loadings and apply to reproduce desired correlations (steps 7, 8) ---------------------------

    L <- suppressWarnings(.paf_iter(rep(1, k), criterion = .001, R = R_inter,
                   n_fac = n_factors, abs_eig = TRUE, crit_type = 1,
                   max_iter = max_iter)$L)

    shared_load[,1:n_factors] <- L

    # get rid of Heywood cases
    shared_load[shared_load > 1] <- 1
    shared_load[shared_load < -1] <- -1

    if (shared_load[1, 1] < 0) {
      shared_load <- shared_load * -1
    }

    for (i in seq_len(k)) {
      if (sum(shared_load[i,] * shared_load[i,]) < 1) {
        unique_load[i, 1] <-
          (1 - sum(shared_load[i,] * shared_load[i,]))
      } else {
        unique_load[i, 1] <- 0
      }
    }

    unique_load <- sqrt(unique_load)

    for (i in seq_len(k)) {
      sim_dat[, i] <- (shared_comp %*% t(shared_load))[, i] +
        unique_comp[, i] * unique_load[i, 1]
    }


      # Replace normal with nonnormal distributions (step 9) ---------------------------------------------------------

      for (i in seq_len(k)) {
        sim_dat <- sim_dat[sort.list(sim_dat[, i]),]
        sim_dat[,i] <- dists[, i]
      }

      # Calculate RMSE correlation, compare to lowest value, take appropriate action (steps 10, 11, 12) --------------

      R_rep <- stats::cor(sim_dat, method = cor_method)
      R_res <- R - R_rep
      RMSE <- sqrt(sum(R_res[lower.tri(R_res)] * R_res[lower.tri(R_res)]) /
                     (.5 * (k * k - k)))

      if (RMSE < best_RMSE) {
        best_RMSE <- RMSE
        R_best <- R_inter
        res_best <- R_res
        R_inter <- R_inter + initial_multiplier * R_res
        t_no_impr <- 0
      } else {
        t_no_impr <- t_no_impr + 1
        current_multiplier <- initial_multiplier * .5 ^ t_no_impr
        R_inter <- R_best + current_multiplier * res_best
      }
  }

  # Construct the data set with the lowest RMSE correlation (step 13) --------------------------------------------

  L <- suppressWarnings(.paf_iter(rep(1, k), criterion = .001, R = R_best,
                n_fac = n_factors, abs_eig = TRUE, crit_type = 1,
                max_iter = max_iter)$L)
  shared_load[, seq_len(n_factors)] <- L

  shared_load[shared_load > 1] <- 1
  shared_load[shared_load < -1] <- -1
  if (shared_load[1, 1] < 0) {
    shared_load <- shared_load * -1
  }

  for (i in seq_len(k)) {
    if (sum(shared_load[i,] * shared_load[i,]) < 1) {
      unique_load[i, 1] <-
        (1 - sum(shared_load[i,] * shared_load[i,]))
    } else {
      unique_load[i, 1] <- 0
    }
  }

  unique_load <- sqrt(unique_load)

  for (i in seq_len(k)) {
    sim_dat[, i] <- (shared_comp %*% t(shared_load))[, i] +
      unique_comp[, i] * unique_load[i, 1]
  }

  for (i in seq_len(k)) {
    sim_dat <- sim_dat[sort.list(sim_dat[, i]),]
    sim_dat[,i] <- dists[, i]
  }

  # Return the simulated data set (step 14) ----------------------------------------------------------------------

  return(sim_dat)
}
