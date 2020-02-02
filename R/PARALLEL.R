#' Parallel analysis
#'
#' Various methods for performing parallel analysis.
#'
#' @param x matrix or data.frame. The real data to compare the simulated eigenvalues
#'  against. Must not contain variables of classes other than numeric. Can be a
#'  correlation matrix or raw data.
#' @param cors logical. Whether x is a correlation matrix (TRUE), or raw data.
#'  against. Must not contain variables of classes other than numeric.
#' @param ncases numeric. The number of cases / observations to simulate. Should only
#'  be specified if \code{x} is left as \code{NULL} as otherwise the dimensions are taken from \code{x}.
#' @param nvars numeric. The number of variables / indicators to simulate. Should only
#'  be specified if \code{x} is left as \code{NULL} as otherwise the dimensions are taken from \code{x}.
#' @param ndatasets numeric. The number of datasets to simulate. Default is 1000.
#' @param percent numeric. A vector of percentiles to take the simulated eigenvalues from.
#'  Default is 95.
#' @param kind character. What kind of parallel analysis to perform. Can be either
#'  "PAF" or "PCA". If using "PAF", the diagonal of the correlation matrices is
#'  replaced by the squared multiple correlations of the indicators. If using "PCA",
#'  the diagonal values of the correlation matrices  are left to be 1.
#' @param sim logical. Whether eigenvalues for parallel analysis should be found
#'  on simulated data. Default is \code{TRUE}.
#' @param resample logical. Whether eigenvalues for parallel analysis should be
#'  found on data sets created by resampling from \code{x}. Default is \code{TRUE}.
#' @param replace logical. Whether, if \code{resample = TRUE}, the resampling should
#'  be done with replacements or not. Default is \code{TRUE}.
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
#' PARALLEL(ncases = 500, nvars = 10)
#' # example without correlation matrix
#' PARALLEL(test_models$case_11b$cormat, ncases = test_models$case_11b$N)
#' }
PARALLEL <- function(x = NULL,
                     cors = TRUE,
                     ncases = NA,
                     nvars = NA,
                     ndatasets = 1000,
                     percent = 95,
                     kind = c("PAF", "PCA"),
                     sim = TRUE,
                     resample = FALSE,
                     replace = TRUE,
                     na_action = c("rm", "none"),
                     n_cores = parallel::detectCores(),
                     decision_rule = c("Means", "Percentile", "Crawford")) {

  kind <- match.arg(kind)
  na_action <- match.arg(na_action)
  decision_rule <- match.arg(decision_rule)

  size_vec <- rep(round(ndatasets / n_cores), n_cores - 1)
  size_vec[n_cores] <- ndatasets - sum(size_vec)

  results_sim <- NULL
  eigvals_real <- NULL
  results_resample <- NULL
  n_factors_resample <- NULL
  n_factors_sim <- NULL
  x_dat <- FALSE

  if (!is.null(x)) {

    if (any(apply(x, 2, class) != "numeric")) {
      warning("Data contained non-numeric columns. Eigenvalues of actual data were not computed and no PA was done using resampling.")
    } else {

      if (any(is.na(x))) {

        if(na_action == "rm") {

          x <- x[stats::complete.cases(x), ]

        }

      }


      if (!is.na(nvars)) {
        warning("nvars was set and data entered. Taking nvars from data")
      }
      nvars <- ncol(x)
      x_dat <- TRUE

      if (isTRUE(cors)) {
        R <- x
        if (isTRUE(resample)) {
          warning("resample was set to TRUE, but correlation matrix was entered. Resampling can only be done on raw data. Setting resample to FALSE")
          resample <- FALSE
        }

      } else {

        if (!is.na(ncases)) {
          warning("ncases was set and data entered. Taking ncases from data")
        }

        ncases <- nrow(x)
        R <- stats::cor(x)

      }



      if (kind == "PAF") {
        # compute smcs
        diag(R) <- 1 - (1 / diag(solve(R)))
      }

      eigvals_real <- matrix(eigen(R, symmetric = TRUE)$values, ncol = 1)
      colnames(eigvals_real) <- "Real Eigenvalues"

      if (isTRUE(resample)) {


        if (kind == "PCA") {

          eigvals <- parallel::mclapply(size_vec, parallel_resample, data = x,
                                        kind = 1, replace = replace,
                                        mc.cores = n_cores)

          eigvals <- do.call(rbind, eigvals)

        } else if (kind == "PAF") {

          eigvals <- parallel::mclapply(size_vec, parallel_resample, data = x,
                                        kind = 2, replace = replace,
                                        mc.cores = n_cores)
          eigvals <- do.call(rbind, eigvals)

        }

        results_resample <- parallel_summarise(eigvals, percent = percent,
                                               ndatasets = ndatasets, nvars = nvars)

        colnames(results_resample) <- c("Means Resample", paste(percent, "Percentile Resample"))

      }

    }

  }

  if (isTRUE(sim)) {

    if (kind == "PCA") {

      eigvals <- parallel::mclapply(size_vec, parallel_sim, ncases = ncases,
                                    nvars = nvars, kind = 1, mc.cores = n_cores)

      eigvals <- do.call(rbind, eigvals)

    } else if (kind == "PAF") {

      eigvals <- parallel::mclapply(size_vec, parallel_sim, ncases = ncases,
                                    nvars = nvars, kind = 2, mc.cores = n_cores)
      eigvals <- do.call(rbind, eigvals)

    }

    results_sim <- parallel_summarise(eigvals, percent = percent,
                                      ndatasets = ndatasets, nvars = nvars)

    colnames(results_sim) <- c("Means Sim", paste(percent, "Percentile Sim"))

  }

  # determine the number of factors to retain
  if (isTRUE(x_dat)) {


    if (decision_rule == "Crawford") {
      # n factors from resampling
      if (!is.null(results_resample)) {
        if ("95 Percentile Resample" %in% colnames(results_resample)) {
          crawford <- c(results_resample[1, "95 Percentile Resample"],
                        results_resample[-1, "Means Resample"])
          n_factors_resample <- which(!(eigvals_real > crawford))[1] - 1
        } else {
          warning("decision_rule == 'Crawford' is specified, but 95 percentile was not used. Using Means instead. To use 'Crawford', make sure to specify percent = 95.")
          n_factors_resample <- which(!(eigvals_real >
                                          results_resample[, "Means Resample"]))[1] - 1
        }
      }

      # n factors from simulated data
      if (!is.null(results_sim)) {
        if ("95 Percentile Sim" %in% colnames(results_sim)) {
          crawford <- c(results_sim[1, "95 Percentile Sim"],
                        results_sim[-1, "Means Sim"])
          n_factors_sim <- which(!(eigvals_real > crawford))[1] - 1
        } else {
          warning("decision_rule == 'Crawford' is specified, but 95 percentile was not used. Using Means instead. To use 'Crawford', make sure to specify percent = 95.")
          n_factors_sim <- which(!(eigvals_real > results_sim[, "Means Sim"]))[1] - 1
        }
      }

    } else if (decision_rule == "Means") {

      # n factors from resampling
      if (!is.null(results_resample)) {
          n_factors_resample <- which(!(eigvals_real >
                                          results_resample[, "Means Resample"]))[1] - 1
      }

      # n factors from simulated data
      if (!is.null(results_sim)) {
          n_factors_sim <- which(!(eigvals_real > results_sim[, "Means Sim"]))[1] - 1
      }

    } else if (decision_rule == "Percentile") {


      # n factors from resampling
      if (!is.null(results_resample)) {
        for (perc_i in percent) {

          pp <- paste(perc_i, "Percentile")

          n_factors_resample[pp] <- which(!(eigvals_real >
                                          results_resample[, paste(pp, "Resample")]))[1] - 1
        }

      }

      # n factors from simulated data
      if (!is.null(results_sim)) {
        for (perc_i in percent) {

          pp <- paste(perc_i, "Percentile")

          n_factors_sim[pp] <- which(!(eigvals_real >
                                         results_sim[, paste(pp, "Sim")]))[1] - 1
        }
      }


    }



  }

  ctrl <- list(
    x_dat = x_dat,
    ncases = ncases,
    nvars = nvars,
    ndatasets = ndatasets,
    percent = percent,
    kind = kind,
    sim = sim,
    resample = resample,
    replace = replace,
    na_action = na_action,
    decision_rule = decision_rule
  )

  eigenvalues <- cbind(eigvals_real, results_sim, results_resample)

  out <- list(
    eigenvalues = eigenvalues,
    n_factors_resample = n_factors_resample,
    n_factors_sim = n_factors_sim,
    ctrl = ctrl
  )

  class(out) <- "PARALLEL"

  return(out)

}
