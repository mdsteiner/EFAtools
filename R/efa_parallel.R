#' Parallel analysis
#'
#' Various methods for performing parallel analysis. This function uses
#' [future_lapply()][future.apply::future_lapply] for which a parallel processing plan can
#' be selected. To do so, call `library(future)` and, for example,
#'  `plan(multisession)`; see examples.
#'
#' @param x matrix or data.frame. The real data to compare the simulated eigenvalues
#'  against. Must not contain variables of classes other than numeric. Can be a
#'  correlation matrix or raw data.
#' @param N numeric. The number of cases / observations to simulate. Only has to
#'  be specified if `x` is either a correlation matrix or `NULL`. If
#'  x contains raw data, `N` is found from the dimensions of `x`.
#' @param n_vars numeric. The number of variables / indicators to simulate.
#' Only has to be specified if `x` is left as `NULL` as otherwise the
#' dimensions are taken from `x`.
#' @param n_datasets numeric. The number of datasets to simulate. Default is 1000.
#' @param percent numeric. The percentile to take from the simulated eigenvalues.
#'  Default is 95.
#' @param eigen_type character. On what the eigenvalues should be found. Can be
#'  either "SMC", "PCA", or "EFA". If using "SMC", the diagonal of the correlation
#'  matrix is replaced by the squared multiple correlations (SMCs) of the
#'  indicators. If using "PCA", the diagonal values of the correlation matrices
#'  are left to be 1. If using "EFA", eigenvalues are found on the correlation
#'  matrices with the final communalities of an EFA solution as diagonal. Default
#'  is `c("PCA", "SMC", "EFA")`, i.e. all three, which costs roughly six times a
#'  single non-EFA type: `"EFA"` fits an EFA to every simulated dataset and
#'  dominates that total. Pass a single type if the run is time-critical.
#' @param use character. Passed to [stats::cor()] if raw data
#' is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. One of `"pearson"`, `"spearman"`, or `"kendall"`,
#'   passed to [stats::cor()]. `"poly"` and `"tetra"` are not supported because
#'   `PARALLEL` compares the data against simulated continuous reference data.
#' Default is "pearson".
#' @param decision_rule character. Which rule to use to determine the number of
#'  factors to retain. Default is `"means"`, which will use the average
#'  simulated eigenvalues. `"percentile"`, uses the percentiles specified
#'  in percent. `"crawford"` uses the 95th percentile for the first factor
#'  and the mean afterwards (based on Crawford et al, 2010). The `"means"` rule
#'  retains a factor whenever its real eigenvalue exceeds the average simulated
#'  one and thus tends to retain more factors than the more conservative
#'  `"percentile"` rule (Glorfeld, 1995).
#' @param n_factors numeric. Number of factors to extract if "EFA" is included in
#' `eigen_type`. Default is 1.
#' @param estimate_control an [estimate_control()] object with the estimation settings for the
#'  [efa_fit()] fits (of both the real and the simulated data) when `"EFA"` is included in
#'  `eigen_type`. `NULL` (default) uses the [efa_fit()] defaults. The fits are unrotated, so no
#'  rotation settings apply.
#' @param ... Additional arguments passed to [efa_fit()]. For example,
#' `estimator`, to change the estimator (default is "PAF"). PAF is more
#' robust, but it will take longer compared to the other estimators
#' available ("ML" and "ULS"). The estimation tuning knobs are not passed here; they live in
#' `estimate_control`, and the standard-error arguments (`se`, `b_boot`, `ci`, `seed`) are
#' not accepted because the fits are internal steps that keep only their eigenvalues.
#'
#' @details Parallel analysis (Horn, 1965) compares the eigenvalues obtained from
#' the sample
#'  correlation matrix against those of null model correlation matrices (i.e.,
#'  with uncorrelated variables) of the same sample size. This way, it accounts
#'  for the variation in eigenvalues introduced by sampling error and thus
#'  eliminates the main problem inherent in the Kaiser-Guttman criterion
#'  ([efa_kgc()]).
#'
#'  Parallel analysis is often argued to be one of the most accurate factor
#'  retention criteria. However, for highly correlated
#'  factor structures it has been shown to underestimate the correct number of
#'  factors. The reason for this is that a null model (uncorrelated variables)
#'  is used as reference. However, when factors are highly correlated, the first
#'  eigenvalue will be much larger compared to the following ones, as
#'  later eigenvalues are conditional on the earlier ones in the sequence and thus
#'  the shared variance is already accounted in the first eigenvalue (e.g.,
#'  Braeken & van Assen, 2017).
#'
#'  The reference eigenvalues are obtained from simulated data, so the suggested number
#'  of factors varies slightly from run to run. Call [set.seed()] beforehand to make a
#'  run reproducible; the result is then also independent of the parallel plan set via
#'  [future::plan()], so it can be reproduced on a machine with a different number of
#'  cores.
#'
#'  The `efa_parallel` function can also be called together with other factor
#'  retention criteria in the [efa_retain()] function.
#'
#' @returns An object of class `efa_retention` (see [print.efa_retention()] and
#'   [plot.efa_retention()] for the print and plot methods). Its main fields are:
#' \item{n_factors}{A named numeric vector with the suggested number of factors for
#'   each requested eigenvalue type (`"PCA"`, `"SMC"`, and/or `"EFA"`). These are
#'   `NA` when no real data are supplied (i.e. only `N` and `n_vars` are given). When
#'   every real eigenvalue exceeds its reference (no crossing is found), all `n_vars`
#'   components are retained and a warning is issued.}
#' \item{results}{A list with one record per eigenvalue type, each holding the real
#'   eigenvalues (when real data were supplied) and the simulated reference
#'   eigenvalues (means and percentiles) used for printing and plotting.}
#' \item{settings}{A list of the settings used.}
#'
#' @source Braeken, J., & van Assen, M. A. (2017). An empirical Kaiser criterion.
#' Psychological Methods, 22, 450 – 466. http://dx.doi.org/10.1037/ met0000074
#'
#' @source Crawford, A. V., Green, S. B., Levy, R., Lo, W. J., Scott, L.,
#' Svetina, D., & Thompson, M. S. (2010). Evaluation of parallel analysis methods
#' for determining the number of factors. Educational and Psychological
#' Measurement, 70(6), 885-901.
#'
#' @source Glorfeld, L. W. (1995). An improvement on Horn's parallel analysis
#' methodology for selecting the correct number of factors to retain. Educational
#' and Psychological Measurement, 55(3), 377-393.
#'
#' @source Horn, J. L. (1965). A rationale and test for the number of factors in
#' factor analysis. Psychometrika, 30(2), 179–185. doi: 10.1007/BF02289447
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
#' # example without real data
#' pa_unreal <- efa_parallel(N = 500, n_vars = 10)
#'
#' # example with correlation matrix with all eigen_types and PAF estimation
#' pa_paf <- efa_parallel(test_models$case_11b$cormat, N = 500)
#'
#' # example with correlation matrix with all eigen_types and ML estimation
#' # this will be faster than the above with PAF)
#' pa_ml <- efa_parallel(test_models$case_11b$cormat, N = 500, estimator = "ML")
#'}
#'
#'\dontrun{
#' # for parallel computation
#' future::plan(future::multisession)
#' pa_faster <- efa_parallel(test_models$case_11b$cormat, N = 500)
#' }

efa_parallel <- function(x = NULL,
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
                     estimate_control = NULL,
                     ...) {

  .reject_flat_knobs(...names(), fn = "efa_parallel")
  .reject_unknown_fit_dots(...names(), fn = "efa_parallel", unrotated = TRUE)
  .reject_rotation_dots(list(...), fn = "efa_parallel")

  if(!is.null(x) && !inherits(x, c("matrix", "data.frame"))){

    cli::cli_abort(
      c("{.arg x} must be {.code NULL}, a correlation matrix, or a data frame/matrix of raw data.",
        "x" = "You supplied {.obj_type_friendly {x}}."),
      class = "efa_input_not_matrix"
    )

  }
  eigen_type <- .match_arg_ci(eigen_type, several.ok = TRUE)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  .reject_poly_reference(cor_method, "efa_parallel")
  decision_rule <- match.arg(decision_rule)
  .assert_estimate_control(estimate_control)
  checkmate::assert_count(n_factors)
  checkmate::assert_count(N, na.ok = TRUE)
  checkmate::assert_count(n_vars, na.ok = TRUE)
  checkmate::assert_count(n_datasets)
  checkmate::assert_number(percent, lower = 0, upper = 100)

  # The simulated datasets are drawn in chunks, one future per chunk, under
  # future.seed = TRUE -- which assigns one L'Ecuyer stream per element of the chunk
  # vector. The number of chunks therefore has to be independent of the number of
  # workers: deriving it from nbrOfWorkers() would make the per-chunk streams, and hence
  # the reference eigenvalues, differ between a sequential and a multisession plan for the
  # same set.seed(). A fixed chunk count keeps a seeded run reproducible on any plan. The
  # trade-off is that the granularity no longer adapts to the pool: a plan with more than
  # 20 workers leaves the surplus idle, so a very wide pool is slower than it would be with
  # worker-matched chunking. 20 is the compromise -- enough chunks to keep a typical pool
  # busy, few enough that the per-chunk dispatch stays negligible. The chunk count never
  # exceeds n_datasets (no empty chunks) and never drops below one, which would divide by
  # zero for the degenerate n_datasets = 0.
  size_vec <- .parallel_chunks(n_datasets, max(1L, min(n_datasets, 20L)))

  # Prepare objects
  results_PCA <- NA
  results_SMC <- NA
  results_EFA <- NA
  eigvals_real_PCA <- NA
  eigvals_real_SMC <- NA
  eigvals_real_EFA <- NA
  n_fac_PCA <- NA
  n_fac_SMC <- NA
  n_fac_EFA <- NA
  x_dat <- FALSE

  if (!is.null(x)){

      .assert_cor_input(x)

      if (!is.na(n_vars)) {
        cli::cli_warn(
          c("Both {.arg n_vars} and raw data were supplied.",
            "i" = "Taking {.arg n_vars} from the data."),
          class = "efa_nvars_from_data"
        )
      }
      n_vars <- ncol(x)
      x_dat <- TRUE

      # Detect or compute the correlation matrix, check it, and smooth it if needed
      prep <- .prepare_cor_input(x, N = N, use = use, cor_method = cor_method,
                                 N_policy = "optional",
                                 singular_tail = "parallel analysis is not possible")
      R <- prep$R
      N <- prep$N
      eigvals_R <- eigen(R, symmetric = TRUE, only.values = TRUE)$values

      if ("PCA" %in% eigen_type) {
        eigvals_real_PCA <- matrix(eigvals_R, ncol = 1)
        colnames(eigvals_real_PCA) <- "Real Eigenvalues"
      }

      if ("SMC" %in% eigen_type) {
        # compute smcs
        R_SMC <- R
        diag(R_SMC) <- .smc_start(R)
        eigvals_real_SMC <- matrix(eigen(R_SMC, symmetric = TRUE,
                                     only.values = TRUE)$values, ncol = 1)
        colnames(eigvals_real_SMC) <- "Real Eigenvalues"
      }

      if ("EFA" %in% eigen_type) {
        # Internal fit used only for its eigenvalues; suppress its warnings so a
        # forwarded estimator that does not converge does not raise a warning from
        # inside efa_parallel().
        eigvals_real_EFA <- matrix(suppressWarnings(
          efa_fit(R, n_factors = n_factors, N = N,
                  estimate_control = estimate_control, ...)$final_eigen),  ncol = 1)
        colnames(eigvals_real_EFA) <- "Real Eigenvalues"
      }

  }

  if (is.na(n_vars)) {
    cli::cli_abort(
      c("{.arg n_vars} was not set and could not be taken from the data.",
        "i" = "Specify {.arg n_vars} and try again."),
      class = "efa_nvars_required"
    )
  }

  if (is.na(N)) {

    cli::cli_abort(
      c("{.arg N} was not set and could not be taken from the data.",
        "i" = "Specify {.arg N} and try again."),
      class = "efa_n_required"
    )

  }

  if (N <= n_vars) {
    cli::cli_abort("{.arg N} must be larger than the number of variables.",
                   class = "efa_n_too_small")
  }

    if ("PCA" %in% eigen_type) {

      eigvals_PCA <- try(future.apply::future_lapply(size_vec, .parallel_sim_eig, N = N,
                                             n_vars = n_vars, eigen_type = 1,
                                             cor_method = cor_method,
                                             future.seed = TRUE),
                         silent = TRUE)

      it_i <- 1
      while (inherits(eigvals_PCA, "try-error") && it_i < 25) {
        eigvals_PCA <- try(future.apply::future_lapply(size_vec, .parallel_sim_eig,
                                                       N = N,
                                                       n_vars = n_vars,
                                                       eigen_type = 1,
                                                       cor_method = cor_method,
                                                       future.seed = TRUE),
                           silent = TRUE)
        it_i <- it_i + 1
      }

      if (inherits(eigvals_PCA, "try-error")) {
        cli::cli_abort(
          c("Eigenvalues from simulated data via {.val PCA} could not be found in 25 tries.",
            "i" = "This is likely due to singular matrices."),
          class = "efa_parallel_sim_failed"
        )
      }

      eigvals_PCA <- do.call(rbind, eigvals_PCA)

      results_PCA <- .parallel_summarise(eigvals_PCA, percent = percent,
                                        n_vars = n_vars)

      colnames(results_PCA) <- c("Means", paste(percent, "Percentile"))

      if (isTRUE(x_dat)) {
        n_fac_PCA <- .determine_factors(decision_rule = decision_rule,
                                        eigvals_real = eigvals_real_PCA,
                                        results = results_PCA,
                                        percent = percent)
      }

    }

    if ("SMC" %in% eigen_type) {

      eigvals_SMC <- try(future.apply::future_lapply(size_vec, .parallel_sim_eig,
                                                     N = N, n_vars = n_vars,
                                                     eigen_type = 2,
                                                     cor_method = cor_method,
                                                     maxit = n_datasets * 10,
                                                     future.seed = TRUE),
                         silent = TRUE)
      it_i <- 1
      while (inherits(eigvals_SMC, "try-error") && it_i < 25) {
        eigvals_SMC <- try(future.apply::future_lapply(size_vec, .parallel_sim_eig,
                                                       N = N,
                                                       n_vars = n_vars,
                                                       eigen_type = 2,
                                                       cor_method = cor_method,
                                                       maxit = n_datasets * 10,
                                                       future.seed = TRUE),
                           silent = TRUE)
        it_i <- it_i + 1
      }

      if (inherits(eigvals_SMC, "try-error")) {
        cli::cli_abort(
          c("Eigenvalues from simulated data via {.val SMCs} could not be found in 25 tries.",
            "i" = "This is likely due to singular matrices."),
          class = "efa_parallel_sim_failed"
        )
      }
      eigvals_SMC <- do.call(rbind, eigvals_SMC)

      results_SMC <- .parallel_summarise(eigvals_SMC, percent = percent,
                                        n_vars = n_vars)

      colnames(results_SMC) <- c("Means", paste(percent, "Percentile"))

      if (isTRUE(x_dat)) {
      n_fac_SMC <- .determine_factors(decision_rule = decision_rule,
                                      eigvals_real = eigvals_real_SMC,
                                      results = results_SMC,
                                      percent = percent)
      }

    }

    if ("EFA" %in% eigen_type) {

      eigvals_EFA <- future.apply::future_lapply(size_vec, .parallel_EFA_sim,
                                             n_vars = n_vars, N = N,
                                             n_factors = n_factors,
                                             cor_method = cor_method,
                                             estimate_control = estimate_control, ...,
                                             future.seed = TRUE)
      eigvals_EFA <- do.call(rbind, eigvals_EFA)

      results_EFA <- .parallel_summarise(eigvals_EFA, percent = percent,
                                        n_vars = n_vars)

      colnames(results_EFA) <- c("Means", paste(percent, "Percentile"))

      if (isTRUE(x_dat)) {
      n_fac_EFA <- .determine_factors(decision_rule = decision_rule,
                                      eigvals_real = eigvals_real_EFA,
                                      results = results_EFA,
                                      percent = percent)
      }

    }

  settings <- list(
    x_dat = x_dat,
    N = N,
    n_vars = n_vars,
    n_datasets = n_datasets,
    percent = percent,
    eigen_type = eigen_type,
    use = use,
    cor_method = cor_method,
    decision_rule = decision_rule,
    n_factors = n_factors
  )

  # one record per requested eigenvalue type: the real eigenvalues (the solid
  # line, absent when no real data are given) plus the simulated reference series
  # (means and percentile) drawn as dashed lines
  sim_list <- list(PCA = results_PCA, SMC = results_SMC, EFA = results_EFA)
  real_list <- list(PCA = eigvals_real_PCA, SMC = eigvals_real_SMC,
                    EFA = eigvals_real_EFA)
  nfac_list <- list(PCA = n_fac_PCA, SMC = n_fac_SMC, EFA = n_fac_EFA)

  results <- list()
  for (et in c("PCA", "SMC", "EFA")) {
    if (!(et %in% eigen_type)) next
    sim <- sim_list[[et]]
    refs <- stats::setNames(lapply(seq_len(ncol(sim)), function(j) sim[, j]),
                            colnames(sim))
    if (isTRUE(x_dat)) {
      n_fac <- nfac_list[[et]]
      y <- as.numeric(real_list[[et]])
      highlight <- if (!is.na(n_fac) && n_fac >= 1) n_fac else NULL
    } else {
      # no real data: no real-eigenvalue series and no suggestion
      n_fac <- NA_real_
      y <- NULL
      highlight <- NULL
    }
    results[[et]] <- list(
      name = et,
      label = paste0(et, " eigenvalues"),
      n_factors = n_fac,
      plot_type = "eigen",
      x = seq_len(n_vars),
      y = y,
      references = refs,
      highlight = highlight
    )
  }

  out <- .new_efa_retention(
    "PARALLEL",
    results = unname(results),
    settings = settings,
    subtitle = paste0("Eigenvalues found using ", cli::ansi_collapse(eigen_type),
                      "; ", n_datasets, " simulated datasets."),
    note = if (isTRUE(x_dat)) {
      paste0("Number of factors retained using the \"", decision_rule,
             "\" decision rule.")
    } else {
      "No data were entered; showing the simulated eigenvalues only. No number of factors is suggested."
    }
  )

  return(out)

}


.parallel_EFA_sim <- function(n_datasets, n_vars, N, n_factors, cor_method,
                              estimate_control = NULL, ...){

  eigvals <- matrix(nrow = n_datasets, ncol = n_vars)

  # The null-model reference data are uncorrelated variables, drawn with the shared
  # multivariate-normal kernel from the identity correlation (the no-factor case).
  R_null <- diag(n_vars)

  for(i in seq_len(n_datasets)){

    x <- .simulate_cfm_mvn(R_null, N)
    R <- stats::cor(x, method = cor_method)
    eigvals_i <- try(suppressWarnings(suppressMessages(
      efa_fit(R, n_factors = n_factors, N = N,
              estimate_control = estimate_control, ...)$final_eigen)), silent = TRUE)
    it_i <- 1
    while (inherits(eigvals_i, "try-error") && it_i < 25) {
      x <- .simulate_cfm_mvn(R_null, N)
      R <- stats::cor(x, method = cor_method)
      eigvals_i <- try(suppressWarnings(suppressMessages(
        efa_fit(R, n_factors = n_factors, N = N,
                estimate_control = estimate_control, ...)$final_eigen)), silent = TRUE)
      it_i <- it_i + 1
    }

    if (inherits(eigvals_i, "try-error")) {
      cli::cli_abort(
        c("Eigenvalues from simulated data via {.val EFA} could not be found in 25 tries.",
          "i" = "This is likely due to singular matrices."),
        class = "efa_parallel_sim_failed"
      )
    }

    eigvals[i,] <- eigvals_i

  }

  return(eigvals)
}

# Reference eigenvalues for one simulation chunk. Horn's (1965) parallel analysis
# requires the reference matrices to be built with the same correlation estimator
# as the observed data. The default Pearson case uses the compiled .parallel_sim()
# for its allocation-light, eigenvalue-only loop; rank-based estimators ("spearman",
# "kendall") draw the null-model data with the shared multivariate-normal kernel
# (.simulate_cfm_mvn on the identity correlation -- uncorrelated variables) and then
# correlate them with stats::cor(), so the reference matches the observed correlation
# estimator rather than always being Pearson.
.parallel_sim_eig <- function(n_datasets, n_vars, N, eigen_type, cor_method,
                              maxit = 10000) {

  if (cor_method == "pearson") {
    return(.parallel_sim(n_datasets, n_vars, N, eigen_type, maxit))
  }

  eig_vals <- matrix(NA_real_, nrow = n_datasets, ncol = n_vars)
  success <- 0L
  iter <- 0L

  # The null-model reference data are uncorrelated variables, drawn with the shared
  # kernel from the identity correlation (the no-factor case).
  R_null <- diag(n_vars)

  # Only the SMC path can reject a draw (a singular simulated matrix) and so needs
  # the maxit retry bound. The PCA path never rejects, so it must run all n_datasets
  # draws regardless of maxit, matching the compiled .parallel_sim() whose PCA loop
  # ignores maxit.
  while (success < n_datasets && (eigen_type != 2 || iter < maxit)) {
    iter <- iter + 1L
    x <- .simulate_cfm_mvn(R_null, N)
    R <- stats::cor(x, method = cor_method)

    if (eigen_type == 2) { # SMC: replace the diagonal with squared multiple
      # correlations; skip the draw if the simulated matrix is singular.
      smc <- try(.smc_start(R), silent = TRUE)
      if (inherits(smc, "try-error")) next
      diag(R) <- smc
    }

    success <- success + 1L
    eig_vals[success, ] <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  }

  if (success < n_datasets) {
    cli::cli_abort("Could not generate enough non-singular matrices.",
                   class = "efa_parallel_sim_failed")
  }

  eig_vals
}

.determine_factors <- function(decision_rule, eigvals_real, results, percent){

# determine the number of factors to retain
if (decision_rule == "crawford") {
  # n factors from resampling
  if ("95 Percentile" %in% colnames(results)) {
    crawford <- c(results[1, "95 Percentile"],
                  results[-1, "Means"])
    n_fac <- which(!(eigvals_real > crawford))[1] - 1
  } else {
    cli::cli_warn(
      c("{.code decision_rule = \"crawford\"} was specified, but the 95th percentile was not used; using means instead.",
        "i" = "To use {.val crawford}, set {.code percent = 95}."),
      class = "efa_parallel_crawford"
    )
    n_fac <- which(!(eigvals_real > results[, "Means"]))[1] - 1
  }

} else if (decision_rule == "means") {
  n_fac <- which(!(eigvals_real > results[, "Means"]))[1] - 1

} else if (decision_rule == "percentile") {

  pp <- paste(percent, "Percentile")
  n_fac <- which(!(eigvals_real > results[, pp]))[1] - 1

}

  # When every real eigenvalue exceeds its reference the rule finds no crossing
  # (`which()` is empty, so the index is NA). Every dimension then sits above the
  # noise reference, so retain all tested components (the same "all-exceed"
  # convention as the empirical Kaiser criterion in [efa_ekc()]) and flag the boundary
  # with a classed warning rather than returning a silent NA.
  if (is.na(n_fac)) {
    n_fac <- length(eigvals_real)
    cli::cli_warn(
      c("All real eigenvalues exceeded the parallel analysis reference; no crossing was found.",
        "i" = "Retaining all {n_fac} component{?s}. This often indicates a near-singular or highly collinear correlation matrix; interpret the suggestion with caution."),
      class = "efa_parallel_no_crossing"
    )
  }

  return(n_fac)

}

.parallel_summarise <- function(eig_vals, percent, n_vars) {

  results <- matrix(NA, nrow = n_vars, ncol = length(percent) + 1)
  results[, 1] <- colMeans(eig_vals)

  # percentile reference series via stats::quantile (type 7, matching
  # psych::fa.parallel) rather than a manual order statistic
  for (root in seq_len(n_vars)) {
    results[root, -1] <- stats::quantile(eig_vals[, root], probs = percent / 100,
                                         names = FALSE, na.rm = TRUE)
  }

  return(results)
}

# Split n_datasets into n_chunks non-negative integer chunks that sum to n_datasets,
# distributing the remainder one per chunk. Avoids a negative final chunk when
# n_datasets is not much larger than the number of chunks.
.parallel_chunks <- function(n_datasets, n_chunks) {
  size_vec <- rep(n_datasets %/% n_chunks, n_chunks)
  rem <- n_datasets %% n_chunks
  if (rem > 0) size_vec[seq_len(rem)] <- size_vec[seq_len(rem)] + 1
  size_vec
}
