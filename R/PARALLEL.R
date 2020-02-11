#' Parallel analysis
#'
#' Various methods for performing parallel analysis.
#'
#' @param x matrix or data.frame. The real data to compare the simulated eigenvalues
#'  against. Must not contain variables of classes other than numeric. Can be a
#'  correlation matrix or raw data.
#' @param cors logical. Whether x is a correlation matrix (TRUE), or raw data.
#'  against. Must not contain variables of classes other than numeric.
#' @param n_cases numeric. The number of cases / observations to simulate. Should only
#'  be specified if \code{x} is left as \code{NULL} as otherwise the dimensions are taken from \code{x}.
#' @param n_vars numeric. The number of variables / indicators to simulate. Should only
#'  be specified if \code{x} is left as \code{NULL} as otherwise the dimensions are taken from \code{x}.
#' @param n_datasets numeric. The number of datasets to simulate. Default is 1000.
#' @param percent numeric. A vector of percentiles to take the simulated eigenvalues from.
#'  Default is 95.
#' @param eigen_type character. On what the eigenvalues should be found. Can be
#'  either "PAF", "PCA", or "FA". If using "PAF", the diagonal of the correlation
#'  matrices is replaced by the squared multiple correlations (SMCs) of the
#'  indicators. If using "PCA", the diagonal values of the correlation matrices
#'  are left to be 1. If using "FA", eigenvalues are found on the correlation
#'  matrices with the final communalities of a 1 factor principal axis factoring
#'  solution as diagonal.
#' @param data_type character. Currently only "sim" is implemented. Finds eigenvalues
#'  for parallel analysis on simulated data.
#' @param replace logical. Currently ignored. Whether, if \code{data_type = "resample"}, the
#'  resampling should be done with replacements or not. Default is \code{TRUE}.
#' @param na_action character. What to do with \code{NA}s occurring in \code{x}. \code{"rm"}
#'  removes them list-wise. \code{"none"} returns a warning if \code{NA}s occur and
#'  in this case does not perform resampling, nor compares the the eigen values
#'  of \code{x} against the simulated ones.
#' @param n_cores numeric. Number of cores to perform the simulations and resampling
#'  on. By default this is found using \code{parallel::detectCores()}.
#' @param decision_rule character. Which rule to use to determine the number of
#'  factors to retain. Default is \code{"mean"}, which will use the average
#'  simulated eigenvalues. \code{"Percentile"}, uses the percentiles specified
#'  in percent. \code{"Crawford"} uses the 95th percentile for the first factor
#'  and the mean afterwards (based on Crawford et al, 2010).
#' @param criterion numeric. The convergence criterion.
#'  If the change in communalities from one iteration to the next is smaller than
#'  this criterion the solution is accepted and the procedure ends. Details
#'  depend on criterion_type. Default is \code{.001}.
#' @param criterion_type character. "max_individual" selects the
#'  maximum change in any of the communalities from one iteration to the next
#'  and tests it against the specified criterion. This is also used by SPSS.
#'  "sums" takes difference of the sum of all communalities in one iteration and
#'  the sum of all communalities in the next iteration and tests it against the
#'  criterion. This procedure is used by the \code{\link[psych:fa]{psych::fa}} function.
#'  Default is \code{"sums"}.
#' @param max_iter numeric. The maximum number of iterations to
#'  perform after which the iterative PAF procedure is halted with a warning.
#'  Default is \code{1000}.
#'
#' @return A list of class PARALLEL containing the following objects
#' \item{eigenvalues}{A matrix containing the eigenvalues of the real and the simulated data.}
#' \item{ctrl}{A list of control settings used in the print function.}
#'
#' @source Crawford, A. V., Green, S. B., Levy, R., Lo, W. J., Scott, L., Svetina, D., & Thompson, M. S. (2010). Evaluation of parallel analysis methods for determining the number of factors. Educational and Psychological Measurement, 70(6), 885-901.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # example without real data
#' PARALLEL(n_cases = 500, n_vars = 10)
#' # example without correlation matrix
#' PARALLEL(test_models$case_11b$cormat, n_cases = test_models$case_11b$N)
#' }
PARALLEL <- function(x = NULL,
                     cors = TRUE,
                     n_cases = NA,
                     n_vars = NA,
                     n_datasets = 1000,
                     percent = 95,
                     eigen_type = c("PAF", "PCA", "FA"),
                     data_type = c("sim"), # , "resample"
                     replace = TRUE,
                     na_action = c("pairwise.complete.obs", "everything", "all.obs", "complete.obs", "na.or.complete"),
                     n_cores = parallel::detectCores(),
                     decision_rule = c("Means", "Percentile", "Crawford"),
                     criterion = .001,
                     criterion_type = c("sums", "max_individual"),
                     max_iter = 1000) {

  eigen_type <- match.arg(eigen_type)
  data_type <- match.arg(data_type)
  na_action <- match.arg(na_action)
  decision_rule <- match.arg(decision_rule)
  criterion_type <- match.arg(criterion_type)

  size_vec <- rep(round(n_datasets / n_cores), n_cores - 1)
  size_vec[n_cores] <- n_datasets - sum(size_vec)

  eigvals_real <- NULL
  results <- NULL
  n_factors <- NA
  x_dat <- FALSE

  if (!is.null(x)) {

    if (any(apply(x, 2, class) != "numeric")) {
      warning("Data contained non-numeric columns. Eigenvalues of actual data were not computed and no PA was done using resampling.")
    } else {

      # if (any(is.na(x))) {
      #
      #   if(na_action == "rm" && isFALSE(cors)) {
      #
      #     x <- x[stats::complete.cases(x), ]
      #
      #   }
      #
      # }


      if (!is.na(n_vars)) {
        warning("n_vars was set and data entered. Taking n_vars from data")
      }
      n_vars <- ncol(x)
      x_dat <- TRUE

      if (isTRUE(cors)) {
        R <- x
        # if (data_type == "resample") {
        #   warning("data_type was set to resample, but correlation matrix was entered. Resampling can only be done on raw data. Setting data_type to sim")
        #   data_type <- "sim"
        # }

      } else {

        if (!is.na(n_cases)) {
          warning("n_cases was set and data entered. Taking n_cases from data")
        }

        n_cases <- nrow(x)
        R <- stats::cor(x, use = na_action)

      }

      if (eigen_type == "PAF") {
        # compute smcs
        diag(R) <- 1 - (1 / diag(solve(R)))
        eigvals_real <- matrix(eigen(R, symmetric = TRUE)$values, ncol = 1)
      } else if (eigen_type == "PCA") {
        eigvals_real <- matrix(eigen(R, symmetric = TRUE)$values, ncol = 1)
      } else if (eigen_type == "FA") {
        eigvals_real <- matrix(parallel_paf(R, criterion, ifelse(criterion_type == "sums",
                                                          2, 1), max_iter), ncol = 1)
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
      #   } else if (eigen_type == "PAF") {
      #
      #     eigvals <- parallel::mclapply(size_vec, parallel_resample, data = x,
      #                                   eigen_type = 2, replace = replace,
      #                                   mc.cores = n_cores)
      #     eigvals <- do.call(rbind, eigvals)
      #
      #   } else if (eigen_type == "FA") {
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

    if (eigen_type == "PCA") {

      eigvals <- parallel::mclapply(size_vec, parallel_sim, n_cases = n_cases,
                                    n_vars = n_vars, eigen_type = 1,
                                    mc.cores = n_cores)

      eigvals <- do.call(rbind, eigvals)

    } else if (eigen_type == "PAF") {

      eigvals <- parallel::mclapply(size_vec, parallel_sim, n_cases = n_cases,
                                    n_vars = n_vars, eigen_type = 2,
                                    mc.cores = n_cores)
      eigvals <- do.call(rbind, eigvals)

    } else if (eigen_type == "FA") {

      eigvals <- parallel::mclapply(size_vec, parallel_paf_sim, n_cases = n_cases,
                                    n_vars = n_vars, criterion = criterion,
                                    crit_type = ifelse(criterion_type == "sums",
                                                       2, 1), max_iter = max_iter,
                                    mc.cores = n_cores)
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
          n_factors <- which(!(eigvals_real > crawford))[1] - 1
        } else {
          warning("decision_rule == 'Crawford' is specified, but 95 percentile was not used. Using Means instead. To use 'Crawford', make sure to specify percent = 95.")
          n_factors <- which(!(eigvals_real > results[, "Means"]))[1] - 1
          decision_rule <- "Means"
        }

    } else if (decision_rule == "Means") {
      n_factors <- which(!(eigvals_real > results[, "Means"]))[1] - 1

    } else if (decision_rule == "Percentile") {

      for (perc_i in percent) {
        pp <- paste(perc_i, "Percentile")
        n_factors[pp] <- which(!(eigvals_real > results[, pp]))[1] - 1
      }

    }
  }

  ctrl <- list(
    x_dat = x_dat,
    n_cases = n_cases,
    n_vars = n_vars,
    n_datasets = n_datasets,
    percent = percent,
    eigen_type = eigen_type,
    data_type = data_type,
    replace = replace,
    na_action = na_action,
    decision_rule = decision_rule
  )

  eigenvalues <- cbind(eigvals_real, results)

  out <- list(
    eigenvalues = eigenvalues,
    n_factors = n_factors,
    ctrl = ctrl
  )

  class(out) <- "PARALLEL"

  return(out)

}
