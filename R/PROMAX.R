## Promax rotation engine: build a varimax base, fit the oblique promax target,
## and reflect/order the solution. Dispatched from .rotate_model(), which resolves
## the `type` preset and guards the single-factor case. Promax reflects and orders
## its varimax base before the fit (ss_factors) or its pattern after the fit
## (eigen), so it finalizes its own solution rather than using .reflect_and_order().
.rotate_promax <- function(L, normalize, P_type, precision, order_type,
                           varimax_type, k) {

  dim_names <- dimnames(L)

  # perform the varimax rotation
  if (varimax_type == "svd") {
    AV <- stats::varimax(L, normalize = normalize, eps = precision)
  } else if (varimax_type == "kaiser") {
    AV <- .VARIMAX_SPSS(L, normalize = normalize, precision = precision)
  }


  if (order_type == "ss_factors") {

    # reflect factors with negative sums
    signs <- .reflect_signs(AV$loadings)
    AV$loadings <- AV$loadings %*% diag(signs, nrow = length(signs))

    # reorder the factors according to largest sums of squares
    ss <- colSums(AV$loadings ^2)
    ss_order <- order(ss, decreasing = TRUE)

    AV$loadings <- AV$loadings[, ss_order]

    AV$rotmat <- (AV$rotmat %*% diag(signs))[, ss_order, drop = FALSE]

    dim_names[[2]] <- dim_names[[2]][ss_order]

  }

  # store the loading matrix in a separate object for use
  A <- AV$loadings

  if (P_type == "unnorm") {

    # Promax target of Hendrickson and White (1964): P = sign(A) |A|^k. Written
    # without dividing by A -- the algebraically identical abs(A^(k + 1)) / A form
    # yields 0/0 = NaN for an exact-zero loading -- as stats::promax also does.
    P <- sign(A) * abs(A)^k

  } else if (P_type == "norm") {

    # SPSS-normalized target: row-normalize A by its row norm before raising to the
    # power, P = sign(A) |A / rowNorm|^k. The row norm is floored so a zero-
    # communality row (norm 0) stays 0 instead of producing 0/0 = NaN, mirroring the
    # Kaiser-weight floor in .VARIMAX_SPSS; an isolated zero loading in a non-zero
    # row already maps to 0 without the floor.
    rn <- sqrt(rowSums(A^2))
    rn[!is.finite(rn) | rn < 1e-15] <- 1
    P <- sign(A) * abs(A / rn)^k
  }

  # run the least squares fit
  fit <- stats::lm.fit(A, P)
  if (fit$rank < ncol(A)) {
    cli::cli_abort(
      c("The varimax base passed to promax is rank-deficient, so the oblique target fit is not identified.",
        "i" = "Inspect the unrotated loadings or extract fewer factors."),
      class = "efa_promax_rank_deficient"
    )
  }
  U <- fit$coefficients

  # rescale the transformation matrix
  D <- diag(solve(t(U) %*% U))
  U <- U %*% diag(sqrt(D))
  dimnames(U) <- NULL

  # compute the factor pattern matrix
  AP <- A %*% U

  # get rotation matrix and prepare output
  U <- AV$rotmat %*% U
  Ui <- solve(U)
  Phi <- Ui %*% t(Ui)

  if (order_type == "eigen") {
    # reflect factors with negative sums
    signs <- .reflect_signs(AP)
    AP <- AP %*% diag(signs, nrow = length(signs))
    Phi <- diag(signs) %*% Phi %*% diag(signs)

    # order the factors by their explained variance, largest first: the
    # Phi-weighted sum of squares diag(Phi L'L), which is the quantity reported as
    # the "SS loadings", so the reported variances decrease monotonically (as in
    # psych)
    ss_rotated <- diag(Phi %*% crossprod(AP))
    eig_order <- order(ss_rotated, decreasing = TRUE)
    AP <- AP[, eig_order, drop = FALSE]
    Phi <- Phi[eig_order, eig_order]

    # the rotation matrix follows the same sign reflection and factor ordering as
    # the pattern, so L_unrot %*% rotmat reproduces the rotated loadings (the
    # ss_factors branch reorders its varimax base before the fit, so U is already
    # in order there)
    U <- (U %*% diag(signs))[, eig_order, drop = FALSE]

    dim_names[[2]] <- dim_names[[2]][eig_order]

  }

  dimnames(AP) <- dim_names

  vars_accounted_rot <- .compute_vars(L_unrot = L, L_rot = AP, Phi = Phi)

  colnames(vars_accounted_rot) <- colnames(AP)

  # get structure matrix
  Structure <- AP %*% Phi
  dimnames(Structure) <- list(rownames(AP), colnames(AP))

  # prepare and return output list
  class(AP) <- "LOADINGS"
  class(Structure) <- "LOADINGS"

  list(rot_loadings = AP,
       Phi = Phi,
       Structure = Structure,
       rotmat = U,
       vars_accounted_rot = vars_accounted_rot)
}
