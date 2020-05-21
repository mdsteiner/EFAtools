#' Parallel analysis
#'
#' Various methods for performing parallel analysis. This function uses
#' \link{future.apply}{future_lapply} for which a parallel processing plan can
#' be selected. To do so, call \code{library(future)} and, for example,
#'  \code{plan(multisession)}; see examples.
#'
#' @param x matrix or data.frame. The real data to compare the simulated eigenvalues
#'  against. Must not contain variables of classes other than numeric. Can be a
#'  correlation matrix or raw data.
#' @param N numeric. The number of cases / observations to simulate. Should only
#'  be specified if \code{x} is left as \code{NULL} as otherwise the dimensions are taken from \code{x}.
#' @param n_vars numeric. The number of variables / indicators to simulate.
#' Should only be specified if \code{x} is left as \code{NULL} as otherwise the
#' dimensions are taken from \code{x}.
#' @param n_datasets numeric. The number of datasets to simulate. Default is 1000.
#' @param percent numeric. A vector of percentiles to take the simulated eigenvalues from.
#'  Default is 95.
#' @param eigen_type character. On what the eigenvalues should be found. Can be
#'  either "SMC", "PCA", or "EFA". If using "SMC", the diagonal of the correlation
#'  matrices is replaced by the squared multiple correlations (SMCs) of the
#'  indicators. If using "PCA", the diagonal values of the correlation matrices
#'  are left to be 1. If using "EFA", eigenvalues are found on the correlation
#'  matrices with the final communalities of an EFA solution as diagonal.
#' @param data_type character. Currently only "sim" is implemented. Finds eigenvalues
#'  for parallel analysis on simulated data.
#' @param replace logical. Currently ignored. Whether, if \code{data_type = "resample"}, the
#'  resampling should be done with replacements or not. Default is \code{TRUE}.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
#' is given as input. Default is "pairwise.complete.obs".
#' @param decision_rule character. Which rule to use to determine the number of
#'  factors to retain. Default is \code{"mean"}, which will use the average
#'  simulated eigenvalues. \code{"Percentile"}, uses the percentiles specified
#'  in percent. \code{"Crawford"} uses the 95th percentile for the first factor
#'  and the mean afterwards (based on Crawford et al, 2010).
#'  @param n_factors numeric. Number of factors to extract if "EFA" is included in
#' \code{eigen_type}. Default is 1.
#'  @param ... Additional arguments passed to \code{\link[EFA]{EFA}}. For example,
#'  the extraction method can be changed here (default is PAF).
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
#'  which the diagonal of the correlation matrix is replaced. Finally, EFA performes
#'  an \code{\link{EFA}} with one factor to estimate
#'  the communalities and based on the correlation matrix with these as diagonal
#'  elements, finds the eigenvalues.
#'
#'  Parallel analysis is often argued to be one of the most accurate factor
#'  retention criteria. However, for highly correlated
#'  factor structures it has been shown to underestimate the correct number of
#'  factors. The reason for this is that a null model (uncorrelated variables)
#'  is used as reference. However, when factors are highly correlated, the first
#'  eigenvalue will relatively be much larger than the following ones, as
#'  later eigenvalues are conditional on the earlier ones in the sequence and thus
#'  the shared variance is already accounted in the first eigenvalue (e.g.,
#'  Braeken & van Assen, 2017).
#'
#' @return A list of class PARALLEL containing the following objects
#' \item{eigenvalues}{A matrix containing the eigenvalues of the real and the simulated data.}
#' \item{n_fac}{The number of factors to retain according to the parallel procedure.}
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
#' @export
#'
#' @examples
#' \dontrun{
#' # example without real data
#' PARALLEL(N = 500, n_vars = 10)
#' # example without correlation matrix
#' PARALLEL(test_models$case_11b$cormat, N = test_models$case_11b$N)
#'
#' # for parallel computation
#' future::plan(future::multisession)
#' PARALLEL(test_models$case_11b$cormat, N = test_models$case_11b$N)
#' }
PARALLEL <- function(x = NULL,
                     N = NA,
                     n_vars = NA,
                     n_datasets = 1000,
                     percent = 95,
                     eigen_type = c("PCA", "SMC", "EFA"),
                     data_type = c("sim"), # , "resample"
                     replace = TRUE,
                     use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                             "everything", "na.or.complete"),
                     decision_rule = c("Means", "Percentile", "Crawford"),
                     n_factors = 1,
                     ...) {

  eigen_type <- match.arg(eigen_type) # SEVERAL OK HERE!
  data_type <- match.arg(data_type)
  use <- match.arg(use)
  decision_rule <- match.arg(decision_rule)
  checkmate::assert_count(n_factors)
  checkmate::assert_count(N, na.ok = TRUE)
  checkmate::assert_count(n_vars, na.ok = TRUE)
  checkmate::assert_count(n_datasets)
  checkmate::assert_number(percent, lower = 0, upper = 100)
  checkmate::assert_flag(replace)

  n_cores <- future::nbrOfWorkers()
  size_vec <- rep(round(n_datasets / n_cores), n_cores - 1)
  size_vec[n_cores] <- n_datasets - sum(size_vec)

  eigvals_real <- NULL
  results <- NULL
  n_fac <- NA
  x_dat <- FALSE

  if (!is.null(x)) {

    if (any(!apply(x, 2, inherits, "numeric"))) {
      warning("Data contained non-numeric columns. Eigenvalues of actual data
              were not computed and no PA was done using resampling.")
    } else {

      if (!is.na(n_vars)) {
        warning("n_vars was set and data entered. Taking n_vars from data")
      }
      n_vars <- ncol(x)
      x_dat <- TRUE


      # Check if it is a correlation matrix
      if(.is_cormat(x)){

        if(any(is.na(x))){

          stop("The correlation matrix you entered contains missing values.
           Analyses are not possible.")

        }

        # if (data_type == "resample") {
        #   warning("data_type was set to resample, but correlation matrix was
        # entered. Resampling can only be done on raw data. Setting data_type to sim")
        #   data_type <- "sim"
        # }

        R <- x

      } else {

        message("x was not a correlation matrix. Correlations are found from entered
            raw data.")

        if (!is.na(N)) {
          warning("N was set and data entered. Taking N from data")
        }

        R <- stats::cor(x, use = use)
        colnames(R) <- colnames(x)
        N <- nrow(x)

      }

      # Check if correlation matrix is invertable, if it is not, stop with message
      R_i <- try(solve(R))

      if (inherits(R_i, "try-error")) {
        stop("Correlation matrix is singular, parallel analysis is not possible")
      }

      # Check if correlation matrix is positive definite
      if(any(eigen(R)$values <= 0)){

        R <- psych::cor.smooth(R)

        }


      if (eigen_type == "SMC") {
        # compute smcs
        diag(R) <- 1 - (1 / diag(solve(R)))
        eigvals_real <- matrix(eigen(R, symmetric = TRUE,
                                     only.values = TRUE)$values, ncol = 1)
      } else if (eigen_type == "PCA") {
        eigvals_real <- matrix(eigen(R, symmetric = TRUE,
                                     only.values = TRUE)$values, ncol = 1)
      } else if (eigen_type == "EFA") {
        eigvals_real <- matrix(EFA(R, n_factors = n_factors, ...)$final_eigen,
                               ncol = 1)
      }


      colnames(eigvals_real) <- "Real Eigenvalues"

      # if (data_type == "resample") {
      #
      #
      #   if (eigen_type == "PCA") {
      #
      #     eigvals <- parallel::mclapply(size_vec, parallel_resample, data = x,
      #                                   eigen_type = 1, replace = replace,
      #                                   mc.cores = n_cores)
      #
      #     eigvals <- do.call(rbind, eigvals)
      #
      #   } else if (eigen_type == "SMC") {
      #
      #     eigvals <- parallel::mclapply(size_vec, parallel_resample, data = x,
      #                                   eigen_type = 2, replace = replace,
      #                                   mc.cores = n_cores)
      #     eigvals <- do.call(rbind, eigvals)
      #
      #   } else if (eigen_type == "EFA") {
      #
      #     eigvals <- parallel::mclapply(size_vec, parallel_paf_resample, data = x,
      #                                   replace = replace, criterion = criterion,
      #                                   crit_type = ifelse(criterion_type == "sums",
      #                                                      2, 1),
      #                                   max_iter = max_iter,
      #                                   mc.cores = n_cores)
      #     eigvals <- do.call(rbind, eigvals)
      #
      #   }
      #
      #   results <- parallel_summarise(eigvals, percent = percent,
      #                                 n_datasets = n_datasets,
      #                                 n_vars = n_vars)
      #
      #   colnames(results) <- c("Means", paste(percent, "Percentile"))
      #
      # }

    }

  }

  if (data_type == "sim") {
    if (is.na(N)) {
      stop('"N" was not set and could not be taken from data. Please specify N
           and try again.')
    }

    if (is.na(n_vars)) {
      stop('"n_vars" was not set and could not be taken from data. Please specify
           n_vars and try again.')
    }

    if (eigen_type == "PCA") {

      eigvals <- future.apply::future_lapply(size_vec, parallel_sim, N = N,
                                             n_vars = n_vars, eigen_type = 1)

      eigvals <- do.call(rbind, eigvals)

    } else if (eigen_type == "SMC") {

      eigvals <- future.apply::future_lapply(size_vec, parallel_sim, N = N,
                                             n_vars = n_vars, eigen_type = 2)
      eigvals <- do.call(rbind, eigvals)

    } else if (eigen_type == "EFA") {

      eigvals <- future.apply::future_lapply(size_vec, .parallel_EFA_sim,
                                             n_vars = n_vars, N = N,
                                             n_factors = n_factors, ...)
      eigvals <- do.call(rbind, eigvals)

    }

    results <- parallel_summarise(eigvals, percent = percent,
                                  n_datasets = n_datasets, n_vars = n_vars)

    colnames(results) <- c("Means", paste(percent, "Percentile"))

  }

  # determine the number of factors to retain
  if (isTRUE(x_dat)) {

    if (decision_rule == "Crawford") {
      # n factors from resampling
        if ("95 Percentile" %in% colnames(results)) {
          crawford <- c(results[1, "95 Percentile"],
                        results[-1, "Means"])
          n_fac <- which(!(eigvals_real > crawford))[1] - 1
        } else {
          warning("decision_rule == 'Crawford' is specified, but 95 percentile
                  was not used. Using Means instead. To use 'Crawford', make sure
                  to specify percent = 95.")
          n_fac <- which(!(eigvals_real > results[, "Means"]))[1] - 1
          decision_rule <- "Means"
        }

    } else if (decision_rule == "Means") {
      n_fac <- which(!(eigvals_real > results[, "Means"]))[1] - 1

    } else if (decision_rule == "Percentile") {

      for (perc_i in percent) {
        pp <- paste(perc_i, "Percentile")
        n_fac[pp] <- which(!(eigvals_real > results[, pp]))[1] - 1
      }

    }
  }

  settings <- list(
    x_dat = x_dat,
    N = N,
    n_vars = n_vars,
    n_datasets = n_datasets,
    percent = percent,
    eigen_type = eigen_type,
    data_type = data_type,
    replace = replace,
    use = use,
    decision_rule = decision_rule
  )

  eigenvalues <- cbind(eigvals_real, results)

  out <- list(
    eigenvalues = eigenvalues,
    n_fac = n_fac,
    settings = settings
  )

  class(out) <- "PARALLEL"

  return(out)

}


.parallel_EFA_sim <- function(n_datasets, n_vars, N, n_factors, ...){

  eigvals <- matrix(nrow = n_vars, ncol = n_datasets)

  for(i in 1:n_datasets){

    x <- matrix(rnorm(N * n_vars), nrow = N, ncol = n_vars)
    R <- stats::cor(x)
    eigvals[, i] <- suppressWarnings(EFA(R, n_factors = n_factors,
                                         ...)$final_eigen)

  }

  return(eigvals)
}
