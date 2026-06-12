#' Extract a list object by its name
#'
#' @author Andreas Soteriades
#'
#' Consider a list of named sub-lists. This function extracts, for each sub-list,
#' the sub-list element that is specified by the user. This function is useful
#' for extracting results from [EFA()] for each permutation run in
#' [EFA_POOLED()].
#'
#' @param alist A list of sub-lists, typically a list of \eqn{m} objects of class
#' `"EFA"`, where \eqn{m} is the number of imputations passed to
#' [EFA_POOLED()].
#' @param object String of length 1. The name of the object to extract e.g.
#' `"h2"` or `"vars_accounted"`.
#'
#' @return A list of length \eqn{m}, with each element containing the extracted
#' `object` for the \eqn{k}th element (\eqn{k = 1,..., m}).
.extract_list_object <- function(alist, object) {
  lapply(
    alist,
    function(x) {
      x[[object]]
    }
  )
}

# Abort (classed) when the suggested lavaan package is unavailable; guards the lavaan
# input paths of OMEGA() and SL(). `call` is forwarded so the error points at the caller.
.require_lavaan <- function(call = rlang::caller_env()) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    cli::cli_abort(
      c("The {.pkg lavaan} package must be installed to use a {.cls lavaan} object as input.",
        "i" = 'Install it with {.code install.packages("lavaan")}.'),
      class = "efa_lavaan_not_installed",
      call = call
    )
  }
}

#' Calculate statistics for a list of matrices
#'
#' @author Andreas Soteriades
#'
#' Given a list of matrices, this function calculates user-supplied statistics
#' (e.g. mean, median) over the matrices. This function is useful for averaging
#' results from [EFA()] (e.g. loadings) over all permutation runs in
#' [EFA_POOLED()].
#'
#' @param alist A list of sub-lists, typically a list of \eqn{m} matrices from
#' `"EFA"`, where \eqn{m} is the number of imputations passed to
#' [EFA_POOLED()].
#' @param stat A function, e.g. `mean` or `sd`.
#'
#' @return A matrix with the aggregated results
#'
.stat_over_list <- function(alist, stat) {
  # In what follows:
  # n: length(alist)
  # n_row: nrow(alist[[i]]), where i is any list element
  # n_col: ncol(alist[[i]]), where i is any list element

  # Convert list to array
  list_to_array <- simplify2array(alist)

  # Apply stat over the array. The 1:2 will apply stat on an array with n_col
  # elements (run apply(list_to_array, 1:2, function(x) x) to see said array).
  # Each element is a n x n_row matrix (run
  # dim(apply(list_to_array, 1:2, function(x) x)) to confirm said array's
  # dimensions are n x n_row x n_col, i.e. n_col matrices with dimensions
  # n x n_row). In each matrix, the i-th row is the i-th column of the matrix in
  # the j-th element of list_to_array.
  # The apply() operation takes the column sum in each element of this new array
  # that has n_col elements, and returns the means in a n_row x n_col matrix.
  apply(
    list_to_array,
    1:2,
    stat
  )
}

#' Covert a `"LOADINGS"` table to matrix or a matrix to `"LOADINGS"`
#'
#' @author Andreas Soteriades
#'
#' The loadings tables returned by [EFA()] are of class
#' `"LOADINGS"`, which prevents applying functions on them. This function
#' allows to change their class to `"matrix"`, and to change back to
#' `"LOADINGS"` when done.
#'
#' @param x A table of class `"matrix"` or `"LOADINGS"`.
#' @param cl A string with the class to change the table to. Should be
#' `"LOADINGS"` or `"matrix"`.
#'
#' @return A table with the loadings, of class either `"LOADINGS"` or
#' `"matrix"`.
.change_class <- function(x, cl = 'matrix') {
  class(x) <- cl
  return(x)
}

#' Confidence intervals around mean
#'
#' @author Andreas Soteriades
#'
#' This function is used internally by [EFA_POOLED()] to calculate
#' confidence intervals (CIs) around the pooled loadings and pooled interfactor
#' correlations.
#'
#' @details
#' The standard error (SE) to use in the creation of the CIs is calculated
#' according to Rubin's (1987) formula (eq. 11 in Hayes & Enders, 2023):
#'
#' \eqn{\sqrt{mean(SE^2) + var(\hat{\theta_{m}}) + var(\hat{\theta_{m}}) / M}}
#'
#' According to Hayes & Enders (2023) p. 42:
#'
#' *\[T\]he first term under the radical represents the average squared
#' standard error (the within imputation sampling variance \[...\]), the second
#' term depends on the variance of the M parameter estimates around their
#' average (the between imputation variance \[...\]), and the final term
#' represents the squared standard error of the pooled estimate \[...\].
#' Conceptually, the first term estimates the sampling error of a complete-data
#' analysis, and the next two terms are essentially correction factors that
#' inflate the standard error to compensate for uncertainty due to the
#' imputations–that is, additional uncertainty (sampling variability) in the
#' parameter estimates caused by missing data.*
#'
#' Currently, it is not possible to calculate the first term, because
#' [EFA()] does not calculate SEs for the loadings. Only the second
#' and third terms are used in the calculation of SE for the CIs.
#'
#' The CI is generally calculated as:
#'
#' CI = Point estimate ± Margin of error,
#'
#' where
#'
#' Margin of error = Critical value × SE of point estimate.
#'
#' To account for situations where the sample size is small, instead of using
#' z-values, the critical value is derived from the *t* distribution with
#' `n - 1` degrees of freedom (Hazra, 2017).
#
#' @param means A matrix with the pooled loadings or pooled interfactor
#' correlations
#' @param sds A matrix with the standard deviations for the loadings or
#' interfactor correlations
#' @param p Numeric. One minus the confidence level of the CIs. Defaults to 0.05
#' for 95% CIs.
#' @param n Numeric. Typically, the number of permutations \eqn{m}. See
#' [EFA_POOLED()].
#'
#' @return A list with the lower and upper CIs.
#'
#' @references
#' Hayes, T. & Enders, C. K. (2023). Maximum likelihood and multiple imputation
#' missing data handling: how they work, and how to make them work in practice.
#' In *APA Handbook of Research Methods in Psychology, Second Edition* Vol.
#' 3. Data Analysis and Research Publication, H. Cooper (Editor-in-Chief).
#'
#' Hazra, A. (2017). Using the confidence interval confidently. *Journal of
#' Thoracic Disease* 9(10), 4125--4130.
#'
#' Rubin, D. B. (1987). Multiple imputation for nonresponse in surveys. Wiley.
#' https://doi.org/10.1002/9780470316696
.calc_cis <- function(means, sds, p = 0.05, n) {

  # Calculate standard error
  rubin_term_1 <- 0 # Currently not calculated- see Details
  rubin_term_2 <- sds ^ 2 # The variances
  rubin_term_3 <- rubin_term_2 / n
  se <- sqrt(rubin_term_1 + rubin_term_2 + rubin_term_3)

  error <- stats::qt(1 - p / 2, n - 1) * se

  lower_ci <- means - error
  upper_ci <-  means + error

  cis <- list(
    lower = lower_ci,
    upper = upper_ci
  )

  return(cis)
}

# varimax criterion for SPSS varimax implementation
.SV <- function(lambda) {

  n <- nrow(lambda)

  # the SPSS manual (ftp://public.dhe.ibm.com/software/analytics/spss/documentation/statistics/23.0/en/client/Manuals/IBM_SPSS_Statistics_Algorithms.pdf)
  # suggests the following formula:
  # sum(n*colSums(lambda**4) - colSums(lambda ** 2) ** 2) / n**2
  # however, the formula below produces results more in line with SPSS
  sum(n*colSums(abs(lambda)) - colSums(lambda ** 4) ** 2) / n**2

}


.onUnload <- function (libpath) {
  library.dynam.unload("EFAtools", libpath)
}

.det_max_factors <- function(m) {
  q <- floor((2*m + 1 - sqrt(8 * m + 9)) / 2)
  if(q < 0) q <- 0
  return(q)
}
