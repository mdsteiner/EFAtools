#' Extract a list object by its name
#'
#' Consider a list of named sub-lists. This function extracts, for each sub-list,
#' the sub-list element that is specified by the user. This function is useful
#' for extracting results from [efa_fit()] for each imputation in
#' [efa_mi()].
#'
#' @param alist A list of sub-lists, typically a list of \eqn{m} objects of class
#' `"efa"`, where \eqn{m} is the number of imputations passed to
#' [efa_mi()].
#' @param object String of length 1. The name of the object to extract e.g.
#' `"h2"` or `"vars_accounted"`.
#'
#' @return A list of length \eqn{m}, with each element containing the extracted
#' `object` for the \eqn{k}th element (\eqn{k = 1,..., m}).
#'
#' @author Andreas Soteriades
#'
#' @keywords internal
.extract_list_object <- function(alist, object) {
  lapply(
    alist,
    function(x) {
      x[[object]]
    }
  )
}

# Set `seed` for the duration of the *calling* function and restore the caller's
# global RNG state when that function returns, so `.Random.seed` carries no lasting
# side effect: the caller's stream is saved and reinstated on exit or -- if none
# existed yet -- the state set.seed() creates is removed again. The restoring
# on.exit() is registered in the caller's frame (`envir`, via do.call) so it fires
# when the public function returns, not when this helper does; the saved state is
# inlined into the handler with bquote() so it resolves in that frame. `add = TRUE`
# composes with any other on.exit() handlers the caller has registered. A NULL
# `seed` is a no-op. Used to make seeded, worker-count-independent bootstraps and
# simulations reproducible without disturbing the caller's RNG stream.
.set_local_seed <- function(seed, envir = parent.frame()) {
  if (is.null(seed)) {
    return(invisible(NULL))
  }
  seed_existed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  saved_seed <- if (seed_existed) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
  do.call(
    on.exit,
    list(
      bquote({
        if (.(seed_existed)) {
          assign(".Random.seed", .(saved_seed), envir = globalenv())
        } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
          rm(".Random.seed", envir = globalenv())
        }
      }),
      add = TRUE
    ),
    envir = envir
  )
  set.seed(seed)
  invisible(NULL)
}

# Abort (classed) when the suggested lavaan package is unavailable; guards the lavaan
# input paths of OMEGA() and efa_schmid_leiman(). `call` is forwarded so the error points
# at the caller.
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

# Evaluate a block of argument checks and re-raise checkmate's assertion failures under
# `efa_invalid_argument`. checkmate's assert_*() functions signal a plain `simpleError`,
# which leaves an invalid argument as the one condition in the package that cannot be
# caught or tested by class; this keeps checkmate's wording but gives it the class the
# rest of the package guarantees. cli_abort() raises an rlang condition rather than a
# `simpleError`, so classed checks inside the block pass through untouched. Only messages
# carrying checkmate's `mstop()` prefix are relabelled: an "object not found" or any other
# unclassed error from the block is an internal fault, not a user input problem, and is
# rethrown unchanged rather than disguised as one. `call` is forwarded so the error points
# at the user's call, not at this helper.
.assert_args <- function(expr, call = rlang::caller_env()) {
  tryCatch(
    expr,
    simpleError = function(e) {
      if (!startsWith(conditionMessage(e), "Assertion on ")) stop(e)
      cli::cli_abort("{conditionMessage(e)}", class = "efa_invalid_argument", call = call)
    }
  )
}

# Second-order Schmid-Leiman group-factor loadings: scale the first-order factor
# loadings by the residual standard deviations of the first-order factors. Only
# the psi diagonal is used (off-diagonal first-order disturbance covariances are
# not residual standard deviations); `nrow` keeps diag() a matrix when there is a
# single first-order factor. Shared by efa_schmid_leiman() and OMEGA() for lavaan
# second-order input. Schmid & Leiman (1957, Psychometrika).
.sl_group_loadings <- function(loadings, psi, col_names) {
  loadings %*% diag(sqrt(diag(psi)[col_names]), nrow = length(col_names))
}

#' Convert an `"efa_loadings"` table to matrix or a matrix to `"efa_loadings"`
#'
#' The loadings tables returned by [efa_fit()] are of class
#' `c("efa_loadings", "LOADINGS")`, which prevents applying functions on them.
#' This function changes their class to `"matrix"`, and changes it back to
#' `"efa_loadings"` when done.
#'
#' @param x A table of class `"matrix"` or `"efa_loadings"`.
#' @param cl A character vector with the class to change the table to. Should be
#' `c("efa_loadings", "LOADINGS")` or `"matrix"`.
#'
#' @return A table with the loadings, of class either `"efa_loadings"` or
#' `"matrix"`.
#'
#' @author Andreas Soteriades
#'
#' @keywords internal
.change_class <- function(x, cl = 'matrix') {
  class(x) <- cl
  return(x)
}

# Vectorize the lower triangle of a symmetric matrix (diagonal included) in column-major
# order, dropping the duplicated upper triangle. This is the standard vech() ordering (the
# one MBESS::vech() returns), which callers rely on when they pair the vectorized entries
# with duplication-matrix weights.
.vech <- function(M) {
  M <- as.matrix(M)
  M[lower.tri(M, diag = TRUE)]
}

# Reconstruct a symmetric matrix from .vech() output.
.unvech <- function(v, k) {
  M <- matrix(NA_real_, k, k)
  M[lower.tri(M, diag = TRUE)] <- v
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  M
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
# (PCA/SMC/EFA); an entry is NA for any convention not requested. `...` is forwarded to efa_fit()
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
    EFA_h2 <- suppressMessages(suppressWarnings(efa_fit(R, n_factors = n_factors, ...)$h2))
    diag(R_EFA) <- EFA_h2
    eigen_R_EFA <- eigen(R_EFA, symmetric = TRUE, only.values = TRUE)$values
  }

  list(PCA = eigen_R_PCA, SMC = eigen_R_SMC, EFA = eigen_R_EFA)
}

# Fall back to `y` when `x` is NULL. Defined here because the base R version was
# only added in R 4.4.0, below the package's minimum.
`%||%` <- function(x, y) if (is.null(x)) y else x
