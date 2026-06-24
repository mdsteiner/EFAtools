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

# Second-order Schmid-Leiman group-factor loadings: scale the first-order factor
# loadings by the residual standard deviations of the first-order factors. Only
# the psi diagonal is used (off-diagonal first-order disturbance covariances are
# not residual standard deviations); `nrow` keeps diag() a matrix when there is a
# single first-order factor. Shared by SL() and OMEGA() for lavaan second-order
# input. Schmid & Leiman (1957, Psychometrika).
.sl_group_loadings <- function(loadings, psi, col_names) {
  loadings %*% diag(sqrt(diag(psi)[col_names]), nrow = length(col_names))
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

# Varimax simplicity criterion monitored for convergence by .VARIMAX_SPSS(); the
# rotation stops once it stabilises. This is the criterion SV given in the SPSS
# Statistics Algorithms manual (FACTOR, "Orthogonal Rotations (Harman, 1976)"),
# evaluated on the Kaiser-normalized loadings: for each factor, n times the sum of
# its fourth-power loadings minus the squared sum of its squared loadings, summed
# over factors and divided by n^2 (Kaiser, 1958).
.SV <- function(lambda) {

  n <- nrow(lambda)

  sum(n * colSums(lambda ** 4) - colSums(lambda ** 2) ** 2) / n ** 2

}


.onUnload <- function (libpath) {
  library.dynam.unload("EFAtools", libpath)
}

.det_max_factors <- function(m) {
  q <- floor((2*m + 1 - sqrt(8 * m + 9)) / 2)
  if(q < 0) q <- 0
  return(q)
}

# Eigenvalues of the correlation matrix under each requested diagonal convention used by the
# eigenvalue-based retention criteria: "PCA" keeps the unit diagonal, "SMC" substitutes the
# squared multiple correlations, and "EFA" substitutes the final communalities of an EFA
# solution (default principal axis factoring extracting `n_factors`). Returns a named list
# (PCA/SMC/EFA); an entry is NA for any convention not requested. `...` is forwarded to EFA()
# for the "EFA" convention.
.three_eigen <- function(R, eigen_type, n_factors = 1, ...) {

  eigen_R_PCA <- NA
  eigen_R_SMC <- NA
  eigen_R_EFA <- NA

  if ("PCA" %in% eigen_type) {
    eigen_R_PCA <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  }

  if ("SMC" %in% eigen_type) {
    R_SMC <- R
    diag(R_SMC) <- .smc_start(R)
    eigen_R_SMC <- eigen(R_SMC, symmetric = TRUE, only.values = TRUE)$values
  }

  if ("EFA" %in% eigen_type) {
    R_EFA <- R
    EFA_h2 <- suppressMessages(suppressWarnings(EFA(R, n_factors = n_factors, ...)$h2))
    diag(R_EFA) <- EFA_h2
    eigen_R_EFA <- eigen(R_EFA, symmetric = TRUE, only.values = TRUE)$values
  }

  list(PCA = eigen_R_PCA, SMC = eigen_R_SMC, EFA = eigen_R_EFA)
}
