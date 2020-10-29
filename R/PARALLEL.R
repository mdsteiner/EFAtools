#' Parallel analysis
#'
#' Various methods for performing parallel analysis. This function uses
#' \link[future.apply]{future_lapply} for which a parallel processing plan can
#' be selected. To do so, call \code{library(future)} and, for example,
#'  \code{plan(multisession)}; see examples.
#'
#' @param x matrix or data.frame. The real data to compare the simulated eigenvalues
#'  against. Must not contain variables of classes other than numeric. Can be a
#'  correlation matrix or raw data.
#' @param N numeric. The number of cases / observations to simulate. Only has to
#'  be specified if \code{x} is either a correlation matrix or \code{NULL}. If
#'  x contains raw data, \code{N} is found from the dimensions of \code{x}.
#' @param n_vars numeric. The number of variables / indicators to simulate.
#' Only has to be specified if \code{x} is left as \code{NULL} as otherwise the
#' dimensions are taken from \code{x}.
#' @param n_datasets numeric. The number of datasets to simulate. Default is 1000.
#' @param percent numeric. The percentile to take from the simulated eigenvalues.
#'  Default is 95.
#' @param eigen_type character. On what the eigenvalues should be found. Can be
#'  either "SMC", "PCA", or "EFA". If using "SMC", the diagonal of the correlation
#'  matrix is replaced by the squared multiple correlations (SMCs) of the
#'  indicators. If using "PCA", the diagonal values of the correlation matrices
#'  are left to be 1. If using "EFA", eigenvalues are found on the correlation
#'  matrices with the final communalities of an EFA solution as diagonal.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
#' is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. Passed to \code{\link[stats:cor]{stats::cor}}
#' Default is "pearson".
#' @param decision_rule character. Which rule to use to determine the number of
#'  factors to retain. Default is \code{"means"}, which will use the average
#'  simulated eigenvalues. \code{"percentile"}, uses the percentiles specified
#'  in percent. \code{"crawford"} uses the 95th percentile for the first factor
#'  and the mean afterwards (based on Crawford et al, 2010).
#' @param n_factors numeric. Number of factors to extract if "EFA" is included in
#' \code{eigen_type}. Default is 1.
#' @param ... Additional arguments passed to \code{\link{EFA}}. For example,
#' the extraction method can be changed here (default is "PAF"). PAF is more
#' robust, but it will take longer compared to the other estimation methods
#' available ("ML" and "ULS").
#'
#' @details Parallel analysis (Horn, 1965) compares the eigenvalues obtained from
#' the sample
#'  correlation matrix against those of null model correlation matrices (i.e.,
#'  with uncorrelated variables) of the same sample size. This way, it accounts
#'  for the variation in eigenvalues introduced by sampling error and thus
#'  eliminates the main problem inherent in the Kaiser-Guttman criterion
#'  (\code{\link{KGC}}).
#'
#'  Three different ways of finding the eigenvalues under the factor model are
#'  implemented, namely "SMC", "PCA", and "EFA". PCA leaves the diagonal elements
#'  of the correlation matrix as they are and is thus equivalent to what is done
#'  in PCA. SMC uses squared multiple correlations as communality estimates with
#'  which the diagonal of the correlation matrix is replaced. Finally, EFA performs
#'  an \code{\link{EFA}} with one factor (can be adapted to more factors) to estimate
#'  the communalities and based on the correlation matrix with these as diagonal
#'  elements, finds the eigenvalues.
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
#'  The \code{PARALLEL} function can also be called together with other factor
#'  retention criteria in the \code{\link{N_FACTORS}} function.
#'
#' @return A list of class PARALLEL containing the following objects
#' \item{eigenvalues_PCA}{A matrix containing the eigenvalues of the real and the simulated data found with eigen_type = "PCA"}
#' \item{eigenvalues_SMC}{A matrix containing the eigenvalues of the real and the simulated data found with eigen_type = "SMC"}
#' \item{eigenvalues_EFA}{A matrix containing the eigenvalues of the real and the simulated data found with eigen_type = "EFA"}
#' \item{n_fac_PCA}{The number of factors to retain according to the parallel procedure with eigen_type = "PCA".}
#' \item{n_fac_SMC}{The number of factors to retain according to the parallel procedure with eigen_type = "SMC".}
#' \item{n_fac_EFA}{The number of factors to retain according to the parallel procedure with eigen_type = "EFA".}
#' \item{settings}{A list of control settings used in the print function.}
#'
#' @source Braeken, J., & van Assen, M. A. (2017). An empirical Kaiser criterion.
#' Psychological Methods, 22, 450 – 466. http://dx.doi.org/10.1037/ met0000074
#'
#' @source Crawford, A. V., Green, S. B., Levy, R., Lo, W. J., Scott, L.,
#' Svetina, D., & Thompson, M. S. (2010). Evaluation of parallel analysis methods
#' for determining the number of factors. Educational and Psychological
#' Measurement, 70(6), 885-901.
#'
#' @source Horn, J. L. (1965). A rationale and test for the number of factors in
#' factor analysis. Psychometrika, 30(2), 179–185. doi: 10.1007/BF02289447
#'
#' @seealso Other factor retention criteria: \code{\link{CD}}, \code{\link{EKC}},
#' \code{\link{HULL}}, \code{\link{KGC}}, \code{\link{SMT}}
#'
#' \code{\link{N_FACTORS}} as a wrapper function for this and all the
#' above-mentioned factor retention criteria.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # example without real data
#' pa_unreal <- PARALLEL(N = 500, n_vars = 10)
#'
#' # example with correlation matrix with all eigen_types and PAF estimation
#' pa_paf <- PARALLEL(test_models$case_11b$cormat, N = 500)
#'
#' # example with correlation matrix with all eigen_types and ML estimation
#' # this will be faster than the above with PAF)
#' pa_ml <- PARALLEL(test_models$case_11b$cormat, N = 500, method = "ML")
#'}
#'
#'\dontrun{
#' # for parallel computation
#' future::plan(future::multisession)
#' pa_faster <- PARALLEL(test_models$case_11b$cormat, N = 500)
#' }

PARALLEL <- function(x = NULL,
                     N = NA,
                     n_vars = NA,
                     n_datasets = 1000,
                     percent = 95,
                     eigen_type = c("PCA", "SMC", "EFA"),
                     use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                             "everything", "na.or.complete"),
                     cor_method = c("pearson", "spearman", "kendall"),
                     decision_rule = c("means", "percentile", "crawford"),
                     n_factors = 1,
                     ...) {


  if(!is.null(x) && !inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither NULL, nor a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data or leave x at NULL.\n"))

  }
  eigen_type <- match.arg(eigen_type, several.ok = TRUE)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  decision_rule <- match.arg(decision_rule)
  checkmate::assert_count(n_factors)
  checkmate::assert_count(N, na.ok = TRUE)
  checkmate::assert_count(n_vars, na.ok = TRUE)
  checkmate::assert_count(n_datasets)
  checkmate::assert_number(percent, lower = 0, upper = 100)

  n_cores <- future::nbrOfWorkers()
  size_vec <- rep(round(n_datasets / n_cores), n_cores - 1)
  size_vec[n_cores] <- n_datasets - sum(size_vec)

  eigenvalues_PCA <- NA
  eigenvalues_SMC <- NA
  eigenvalues_EFA <- NA
  n_fac_PCA <- NA
  n_fac_SMC <- NA
  n_fac_EFA <- NA
  x_dat <- FALSE

  if (!is.null(x)){

      if (!is.na(n_vars)) {
        warning(crayon::yellow$bold("!"), crayon::yellow(" n_vars was set and data entered. Taking n_vars from data\n"))
      }
      n_vars <- ncol(x)
      x_dat <- TRUE

      # Check if it is a correlation matrix
      if(.is_cormat(x)){

        R <- x

      } else {

        message(cli::col_cyan(cli::symbol$info, " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n"))

        if (!is.na(N)) {
          warning(crayon::yellow$bold("!"), crayon::yellow(" 'N' was set and data entered. Taking N from data.\n"))
        }

        R <- stats::cor(x, use = use, method = cor_method)
        colnames(R) <- colnames(x)
        N <- nrow(x)

      }

      # Check if correlation matrix is invertable, if it is not, stop with message
      R_i <- try(solve(R), silent = TRUE)

      if (inherits(R_i, "try-error")) {
        stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Correlation matrix is singular, parallel analysis is not possible.\n"))
      }

      # Check if correlation matrix is positive definite
      if(any(eigen(R)$values <= 0)){

        R <- psych::cor.smooth(R)

      }

      if ("PCA" %in% eigen_type) {
        eigvals_real_PCA <- matrix(eigen(R, symmetric = TRUE,
                                     only.values = TRUE)$values, ncol = 1)
        colnames(eigvals_real_PCA) <- "Real Eigenvalues"
      }

      if ("SMC" %in% eigen_type) {
        # compute smcs
        R_SMC <- R
        diag(R_SMC) <- 1 - (1 / diag(solve(R_SMC)))
        eigvals_real_SMC <- matrix(eigen(R_SMC, symmetric = TRUE,
                                     only.values = TRUE)$values, ncol = 1)
        colnames(eigvals_real_SMC) <- "Real Eigenvalues"
      }

      if ("EFA" %in% eigen_type) {
        eigvals_real_EFA <- matrix(EFA(R, n_factors = n_factors, N = N,
                                   ...)$final_eigen,  ncol = 1)
        colnames(eigvals_real_EFA) <- "Real Eigenvalues"
      }

  }

  if (is.na(n_vars)) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(' "n_vars" was not set and could not be taken from data. Please specify n_vars and try again.\n'))
  }

  if (is.na(N)) {

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(' "N" was not set and could not be taken from data. Please specify N and try again.\n'))

  }

    if ("PCA" %in% eigen_type) {

      eigvals_PCA <- try(future.apply::future_lapply(size_vec, .parallel_sim, N = N,
                                             n_vars = n_vars, eigen_type = 1,
                                             future.seed = TRUE),
                         silent = TRUE)

      it_i <- 1
      while (inherits(eigvals_PCA, "try-error") && it_i < 25) {
        eigvals_PCA <- try(future.apply::future_lapply(size_vec, .parallel_sim,
                                                       N = N,
                                                       n_vars = n_vars,
                                                       eigen_type = 1,
                                                       future.seed = TRUE),
                           silent = TRUE)
        it_i <- it_i + 1
      }

      if (inherits(eigvals_PCA, "try-error")) {
        stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Eigenvalues based on simulated data and PCA could not be found in 25 tries; likely due to occurrences of singular matrices. Aborting.\n"))
      }

      eigvals_PCA <- do.call(rbind, eigvals_PCA)

      results_PCA <- .parallel_summarise(eigvals_PCA, percent = percent,
                                        n_datasets = n_datasets, n_vars = n_vars)

      colnames(results_PCA) <- c("Means", paste(percent, "Percentile"))

      if (isTRUE(x_dat)) {
        n_fac_PCA <- .determine_factors(decision_rule = decision_rule,
                                        eigvals_real = eigvals_real_PCA,
                                        results = results_PCA,
                                        percent = percent)

        eigenvalues_PCA <- cbind(eigvals_real_PCA, results_PCA)

      } else {

        eigenvalues_PCA <- results_PCA

      }

    }

    if ("SMC" %in% eigen_type) {

      eigvals_SMC <- try(future.apply::future_lapply(size_vec, .parallel_sim,
                                                     N = N, n_vars = n_vars,
                                                     eigen_type = 2,
                                                     maxit = n_datasets * 10,
                                                     future.seed = TRUE),
                         silent = TRUE)
      it_i <- 1
      while (inherits(eigvals_SMC, "try-error") && it_i < 25) {
        eigvals_SMC <- try(future.apply::future_lapply(size_vec, .parallel_sim,
                                                       N = N,
                                                       n_vars = n_vars,
                                                       eigen_type = 2,
                                                       maxit = n_datasets * 10,
                                                       future.seed = TRUE),
                           silent = TRUE)
        it_i <- it_i + 1
      }

      if (inherits(eigvals_SMC, "try-error")) {
        stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Eigenvalues based on simulated data and SMCs could not be found in 25 tries; likely due to occurrences of singular matrices. Aborting.\n"))
      }
      eigvals_SMC <- do.call(rbind, eigvals_SMC)

      results_SMC <- .parallel_summarise(eigvals_SMC, percent = percent,
                                        n_datasets = n_datasets, n_vars = n_vars)

      colnames(results_SMC) <- c("Means", paste(percent, "Percentile"))

      if (isTRUE(x_dat)) {
      n_fac_SMC <- .determine_factors(decision_rule = decision_rule,
                                      eigvals_real = eigvals_real_SMC,
                                      results = results_SMC,
                                      percent = percent)

      eigenvalues_SMC <- cbind(eigvals_real_SMC, results_SMC)

      } else {

        eigenvalues_SMC <- results_SMC

      }

    }

    if ("EFA" %in% eigen_type) {

      eigvals_EFA <- future.apply::future_lapply(size_vec, .parallel_EFA_sim,
                                             n_vars = n_vars, N = N,
                                             n_factors = n_factors, ...,
                                             future.seed = TRUE)
      eigvals_EFA <- do.call(rbind, eigvals_EFA)

      results_EFA <- .parallel_summarise(eigvals_EFA, percent = percent,
                                        n_datasets = n_datasets, n_vars = n_vars)

      colnames(results_EFA) <- c("Means", paste(percent, "Percentile"))

      if (isTRUE(x_dat)) {
      n_fac_EFA <- .determine_factors(decision_rule = decision_rule,
                                      eigvals_real = eigvals_real_EFA,
                                      results = results_EFA,
                                      percent = percent)

      eigenvalues_EFA <- cbind(eigvals_real_EFA, results_EFA)

      } else {

        eigenvalues_EFA <- results_EFA

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

  out <- list(
    eigenvalues_PCA = eigenvalues_PCA,
    eigenvalues_SMC = eigenvalues_SMC,
    eigenvalues_EFA = eigenvalues_EFA,
    n_fac_PCA = n_fac_PCA,
    n_fac_SMC = n_fac_SMC,
    n_fac_EFA = n_fac_EFA,
    settings = settings
  )

  class(out) <- "PARALLEL"

  return(out)

}


.parallel_EFA_sim <- function(n_datasets, n_vars, N, n_factors, ...){

  eigvals <- matrix(nrow = n_datasets, ncol = n_vars)

  for(i in seq_len(n_datasets)){

    x <- matrix(stats::rnorm(N * n_vars), nrow = N, ncol = n_vars)
    R <- stats::cor(x)
    eigvals_i <- try(suppressWarnings(suppressMessages(EFA(R, n_factors = n_factors, N = N,
                                          ...)$final_eigen)), silent = TRUE)
    it_i <- 1
    while (inherits(eigvals_i, "try-error") && it_i < 25) {
      x <- matrix(stats::rnorm(N * n_vars), nrow = N, ncol = n_vars)
      R <- stats::cor(x)
      eigvals_i <- try(suppressWarnings(suppressMessages(EFA(R, n_factors = n_factors, N = N,
                                            ...)$final_eigen)), silent = TRUE)
      it_i <- it_i + 1
    }

    if (inherits(eigvals_i, "try-error")) {
      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Eigenvalues based on simulated data and EFA could not be found in 25 tries; likely due to occurrences of singular matrices. Aborting.\n"))
    }

    eigvals[i,] <- eigvals_i

  }

  return(eigvals)
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
    warning(crayon::yellow$bold("!"), crayon::yellow(" decision_rule == 'crawford' is specified, but 95 percentile was not used. Using means instead. To use 'crawford', make sure to specify percent = 95.\n"))
    n_fac <- which(!(eigvals_real > results[, "Means"]))[1] - 1
    decision_rule <- "means"
  }

} else if (decision_rule == "means") {
  n_fac <- which(!(eigvals_real > results[, "Means"]))[1] - 1

} else if (decision_rule == "percentile") {

  pp <- paste(percent, "Percentile")
  n_fac <- which(!(eigvals_real > results[, pp]))[1] - 1

}

  return(n_fac)

}

.parallel_summarise <- function(eig_vals, percent, n_datasets, n_vars) {

  results <- matrix(NA, nrow = n_vars, ncol = length(percent) + 1)
  results[, 1] <- colMeans(eig_vals)

  for (root in seq_len(n_vars)) {
    for (perc_i in seq_along(percent)) {
      ind <- round((percent[perc_i] * n_datasets) / 100)
      results[root, 1 + perc_i] <- sort(eig_vals[, root])[ind]
    }
  }

  return(results)
}
