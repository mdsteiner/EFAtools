# Factor congruence metrics and solution alignment: Tucker congruence between two loading
# matrices, and reordering/sign-matching a solution to a target.

.factor_congruence <- function(x, y, na.rm = TRUE, skip_checks = FALSE) {

  if (isFALSE(skip_checks)) {

    if (any(is.na(x) | any(is.na(y)))) {
      if (isTRUE(na.rm)) {
        cli::cli_warn(
          c("The input contained missing values.",
            "i" = "The analysis is performed on complete cases."),
          class = "efa_missing_complete"
        )
        if (any(is.na(x))) {
          xc <- x[stats::complete.cases(x), ]
          y <- y[stats::complete.cases(x), ]
          x <- xc
        }
        if (any(is.na(y))) {
          yc <- y[stats::complete.cases(y), ]
          x <- x[stats::complete.cases(y), ]
          y <- yc
        }
      } else {
        cli::cli_warn(
          c("The input contained missing values.",
            "i" = "Check your data or rerun with {.code na.rm = TRUE}."),
          class = "efa_missing_check"
        )
      }
    }

  }


  nx <- dim(x)[2]
  ny <- dim(y)[2]
  cross <- t(y) %*% x
  sumsx <- sqrt(1/diag(t(x) %*% x))
  sumsy <- sqrt(1/diag(t(y) %*% y))
  result <- sumsy * (cross * rep(sumsx, each = ny))
  return(t(result))

}


### Align solution according to factor congruence

.align_solution <- function(L_target, L, Phi = NULL) {

  L_target <- as.matrix(L_target)
  L <- as.matrix(L)
  if (!identical(dim(L_target), dim(L))) {
    cli::cli_abort("{.arg L_target} and {.arg L} must have identical dimensions.",
                   class = "efa_dim_mismatch")
  }

  m <- ncol(L_target)
  target_colnames <- colnames(L_target)
  target_rownames <- rownames(L_target)

  congruence <- .tucker_congruence(L_target, L)

  if (m > 1) {

    cost <- max(abs(congruence), na.rm = TRUE) - abs(congruence)

    factor_order <- as.integer(clue::solve_LSAP(cost))

    # Reorder columns of L so column j is matched to target column j.
    L <- L[, factor_order, drop = FALSE]

    if (!is.null(Phi)) {
      Phi <- Phi[factor_order, factor_order, drop = FALSE]
    }

    # Get signs after reordering.
    factor_sign <- sign(diag(crossprod(L, L_target)))

    # Avoid accidental zero columns if an exact dot product is zero.
    factor_sign[factor_sign == 0] <- 1

    factor_sign_matrix <- diag(factor_sign, nrow = m, ncol = m)

    # Apply signs.
    L <- L %*% factor_sign_matrix

    if (!is.null(Phi)) {
      Phi <- factor_sign_matrix %*% Phi %*% factor_sign_matrix
    }

  } else {

    factor_order <- 1L
    factor_sign <- as.numeric(sign(congruence[1L, 1L]))
    if (factor_sign == 0) factor_sign <- 1

    L <- L * factor_sign

    # Phi unchanged for one factor, because sign^2 = 1.
  }

  dimnames(L) <- list(target_rownames, target_colnames)
  if (!is.null(Phi)) {
    dimnames(Phi) <- list(target_colnames, target_colnames)
  }

  list(
    loadings = L,
    Phi = Phi,
    factor_order = factor_order,
    factor_sign = factor_sign
  )
}

#' Tucker congruence between factors
#'
#' Compute the Tucker congruence matrix between the columns of two loading
#' matrices.
#'
#' @param L1 Numeric matrix.
#' @param L2 Numeric matrix with the same dimensions as `L1`.
#'
#' @returns A square matrix whose `(i, j)` entry is the Tucker congruence
#'   between column `i` of `L1` and column `j` of `L2`.
#'
#' @references
#' Lorenzo-Seva, U., and ten Berge, J. M. F. (2006). Tucker's congruence
#' coefficient as a meaningful index of factor similarity. *Methodology*, 2,
#' 57-64.
.tucker_congruence <- function(L1, L2) {
  L1 <- .procrustes_as_matrix(L1, "L1")
  L2 <- .procrustes_as_matrix(L2, "L2")

  if (!identical(dim(L1), dim(L2))) {
    cli::cli_abort("{.arg L1} and {.arg L2} must have identical dimensions.",
                   class = "efa_dim_mismatch")
  }

  # Column Euclidean norms.
  n1 <- sqrt(colSums(L1 * L1))
  n2 <- sqrt(colSums(L2 * L2))

  if (any(n1 <= sqrt(.Machine$double.eps))) {
    cli::cli_abort("{.arg L1} contains at least one zero or near-zero column.",
                   class = "efa_zero_column")
  }
  if (any(n2 <= sqrt(.Machine$double.eps))) {
    cli::cli_abort("{.arg L2} contains at least one zero or near-zero column.",
                   class = "efa_zero_column")
  }

  # Full pairwise congruence matrix.
  cong <- crossprod(L1, L2) / tcrossprod(n1, n2)

  rownames(cong) <- colnames(L1)
  colnames(cong) <- colnames(L2)
  cong
}
