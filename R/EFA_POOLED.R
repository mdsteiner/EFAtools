#' Exploratory factor analysis on multiple data imputations
#'
#' @author Andreas Soteriades
#'
#' @description
#' Consider a dataset with missing values that has been imputed \eqn{m} times.
#' This function runs \code{\link{EFA}} on each imputed dataset, averages the
#' estimated parameters (loadings) and calculates confidence intervals (CIs) on
#' the pooled (mean) loadings.
#'
#' @details
#' Multiple imputation (van Buuren, 2018) was developed by Rubin (1976). The
#' idea stems from the fact that that imputing one value (single imputation) for
#' the missing value could not be correct in general. Instead, impute the data
#' \eqn{m} times, i.e. generate multiple versions of the data table, where each
#' version has been imputed. Then, fit \eqn{m} models (e.g. \eqn{m} regressions
#' or \eqn{m} EFA models) and average (pool) the \eqn{m} estimated parameters
#' (e.g. the \eqn{m} regression coefficients of a given predictor or the \eqn{m}
#' loadings of a given variable and factor). This average is your estimate.
#'
#' In \code{\link{EFA_POOLED}}, the pooled loading \eqn{\bar{\lambda_{ij}}}for
#' the \eqn{i}th variable and \eqn{j}th factor is calculated as
#' \eqn{\frac{1}{m}\sum_{k=1}^{m}{\lambda_{ij}^k}}, where \eqn{m} is the number
#' of imputations
#'
#' Similarly, for the \eqn{i}th variable, the communality \code{h2} returned by
#' \code{\link{EFA_POOLED}} is calculated from the pooled loadings as
#' \eqn{\sum_{j=1}^{p}{\bar{\lambda_{ij}}^2}}, where \eqn{p} is the number of
#' factors.
#'
#' @source Rubin, D.B. (1976). Inference and Missing Data. Biometrika 63 (3):
#'   581–90.
#' @source van Buuren, S. (2018). Flexible Imputation of Missing Data. Second
#'   Edition. CRC Press. \url{https://stefvanbuuren.name/fimd/}
#'
#' @param data_list A list of length \eqn{m}, where \eqn{m} is the number of
#' imputations. Each list element \eqn{k} is a data frame or matrix of raw data
#' or a matrix with correlations for data imputation \eqn{k}. See argument
#' \code{x} in \code{\link{EFA}}.
#' @param p Numeric. One minus the confidence level of the CIs. Defaults to 0.05
#' for 95\% CIs.
#' @param ... Additional arguments passed to \code{\link{EFA}}.
#'
#' @return A list of class EFA containing (a subset of) the following:
#'
#' \item{h2}{Final communality estimates from the unrotated solution. This is
#' based on the squared pooled unrotated loadings, \code{unrot_loadings}.}
#' \item{unrot_loadings}{Loading matrix containing the final, pooled unrotated
#' loadings.}
#' \item{vars_accounted}{Matrix of explained variances and sums of squared
#' loadings. Based on the pooled unrotated loadings.}
#' \item{fit_indices}{For ML and ULS: Pooled fit indices derived from the unrotated
#' factor loadings: Chi Square, including significance level, degrees of freedom
#' (df), Comparative Fit Index (CFI), Root Mean Square Error of Approximation
#' (RMSEA), including its 90\% confidence interval, Akaike Information Criterion
#' (AIC), Bayesian Information Criterion (BIC), and the common part accounted
#' for (CAF) index as proposed by Lorenzo-Seva, Timmerman, & Kiers (2011).
#' For PAF, only the CAF and dfs are returned.}
#' \item{rot_loadings}{Loading matrix containing the final pooled rotated
#' loadings (pattern matrix).}
#' \item{Phi}{The pooled factor intercorrelations (only for oblique rotations).}
#' \item{Structure}{The pooled structure matrix (only for oblique rotations).}
#' \item{vars_accounted_rot}{Matrix of explained variances and sums of squared
#' loadings. Based on pooled rotated loadings and, for oblique rotations, the
#' pooled factor intercorrelations.}
#' \item{settings}{A list of the settings used.}
#' \item{ci_unrot_loadings}{CIs for the pooled unrotated loadings.}
#' \item{ci_rot_loadings}{CIs for the pooled rotated loadings.}
#' \item{ci_structure_loadings}{CIs for the pooled structure loadings (only for
#' oblique rotations).}
#' \item{ci_phis}{CIs for the pooled factor intercorrelations (only for oblique
#' rotations).}
#' \item{fits}{List of length \eqn{m} with the results from the \code{\link{EFA}}
#' runs for each imputed dataset.}
#'
#' @export
#'
#' @examples
#' dataset <- psych::bfi[1:500, 1:25] # Create dataset
#' apply(dataset, 2, function(x) sum(is.na(x))) # Any NAs in the columns?
#'
#' # Function to impute columns
#' impute_vector <- function(x, func = mean) {
#'   replace(
#'     x,
#'     is.na(x),
#'     func(x, na.rm = TRUE)
#'   )
#' }
#'
#' # Impute data using mean
#' dataset_imp_mean <- apply(
#'   dataset,
#'   2,
#'   impute_vector,
#'   func = mean
#' )
#'
#' # Impute data using mean
#' dataset_imp_median <- apply(
#'   dataset,
#'   2,
#'   impute_vector,
#'   func = median
#' )
#'
#' # List of imputed datasets, one using the mean, the other using the median
#' data_list <- list(
#'   dataset_imp_mean,
#'   dataset_imp_median
#' )
#'
#' # Run EFA on each imputed dataset, pool the results
#' efa_pooled <- EFA_POOLED(data_list, p = 0.05, n_factors = 3,
#'                          type = "EFAtools", method = "ML", rotation = "oblimin")
#'
EFA_POOLED <- function(data_list, p = 0.05, ...) {
  efa_args <- list(...) # Arguments of EFA()

  checkmate::assert_list(data_list)
  lapply(
    data_list,
    checkmate::assert_multi_class,
    c('matrix', 'data frame')
  )
  checkmate::assert_number(p, na.ok = FALSE)

  # Fit EFA for each dataset
  fits <- lapply(
    data_list,
    function(data_list_subset, ...) {
      # data_list_subset is the i-th element of data_list, and is a data frame
      # or matrix of raw data or matrix with correlations, i.e. it is argument
      # x of EFA(). In each lapply() repetition, we modify x to be the i-th
      # element of data_list. We then run EFA() for this x
      efa_args$x <- data_list_subset
      do.call(EFAtools::EFA, efa_args)
    }
  )

  # List of structure loadings matrices
  structure_loadings <- .extract_list_object(fits, 'Structure')

  # For orthogonal rotations, structure_loadings will be empty
  check_empty <- structure_loadings %>%
    Filter(Negate(is.null), .) %>%
    rlang::is_empty()

  # List of unrotated loadings matrices
  unrot_loadings <- fits %>%
    .extract_list_object('unrot_loadings') %>%
    lapply(.change_class, 'matrix')

  # List of rotated loadings matrices
  rot_loadings <- fits %>%
    .extract_list_object('rot_loadings') %>%
    lapply(.change_class, 'matrix')

  # Target-rotate the rotated loadings matrices. This is necessary for averaging
  # parameters estimated using different imputations
  target_rotations <- lapply(
    rot_loadings[-1], # First one is the target- redundant to self-target-rotate
    psych::target.rot,
    keys = rot_loadings[[1]]
  )

  # List of rotated loadings matrices, now target-rotated
  rot_loadings[-1] <- lapply(
    target_rotations,
    function(x) {
      .change_class(x$loadings, 'matrix')
    }
  )

  # Factor correlations
  if (!check_empty) {
    phis <- .extract_list_object(fits, 'Phi')
  }

  # Means of estimated parameters
  mean_unrot_loadings <- .stat_over_list(unrot_loadings, mean)
  mean_rot_loadings <- .stat_over_list(rot_loadings, mean)
  if (!check_empty) {
    mean_structure_loadings <- .stat_over_list(structure_loadings, mean)
    mean_phis <- .stat_over_list(phis, mean)
  }

  # Standard deviations of estimated parameters
  sd_unrot_loadings <- .stat_over_list(unrot_loadings, sd)
  sd_rot_loadings <- .stat_over_list(rot_loadings, sd)
  if (!check_empty) {
    sd_structure_loadings <- .stat_over_list(structure_loadings, sd)
    sd_phis <- .stat_over_list(phis, sd)
  }

  # CIs of estimated parameters
  n_perm <- length(fits)
  ci_unrot_loadings <- .calc_cis(
    mean_unrot_loadings,
    sd_unrot_loadings,
    p,
    n = n_perm
  )

  ci_rot_loadings <- .calc_cis(
    mean_rot_loadings,
    sd_rot_loadings,
    p,
    n = n_perm
  )

  if (!check_empty) {
    ci_structure_loadings <- .calc_cis(
      mean_structure_loadings,
      sd_structure_loadings,
      p,
      n = n_perm
    )

    ci_phis <- .calc_cis(
      mean_phis,
      sd_phis,
      p,
      n = n_perm
    )
  }

  # Communalities from mean unrotated loadings
  mean_h2 <- rowSums(mean_unrot_loadings ^ 2)

  # Variances, mean unrotated loadings
  mean_vars_accounted <- EFAtools:::.compute_vars(
    L_unrot = mean_unrot_loadings,
    L_rot = mean_unrot_loadings,
    Phi = NULL
  )

  # Variances, mean rotated loadings
  mean_vars_accounted_rot <- EFAtools:::.compute_vars(
    L_unrot = mean_rot_loadings,
    L_rot = mean_rot_loadings,
    Phi = mean_phis
  )

  mean_fit_indices <- fits %>%
    .extract_list_object('fit_indices') %>%
    lapply(unlist) %>%
    simplify2array() %>%
    rowMeans() %>%
    as.list()

  # Ensure loadings tables class is consistent with EFA() output
  mean_unrot_loadings <- .change_class(mean_unrot_loadings, "LOADINGS")
  mean_rot_loadings <- .change_class(mean_rot_loadings, "LOADINGS")

  if (!check_empty) {
    results <- list(
      h2 = mean_h2,
      unrot_loadings = mean_unrot_loadings,
      vars_accounted = mean_vars_accounted,
      fit_indices = mean_fit_indices,
      rot_loadings = mean_rot_loadings,
      Phi = mean_phis,
      Structure = mean_structure_loadings,
      vars_accounted_rot = mean_vars_accounted_rot,
      settings = fits[[1]]$settings,
      ci_unrot_loadings = ci_unrot_loadings,
      ci_rot_loadings = ci_rot_loadings,
      ci_structure_loadings = ci_structure_loadings,
      ci_phis = ci_phis,
      fits = fits
    )
  } else {
    results <- list(
      h2 = mean_h2,
      unrot_loadings = mean_unrot_loadings,
      vars_accounted = mean_vars_accounted,
      fit_indices = mean_fit_indices,
      rot_loadings = mean_rot_loadings,
      vars_accounted_rot = mean_vars_accounted_rot,
      settings = fits[[1]]$settings,
      ci_unrot_loadings = ci_unrot_loadings,
      ci_rot_loadings = ci_rot_loadings,
      fits = fits
    )
  }

  # class(results) <- 'EFA'

  return(results)
}
