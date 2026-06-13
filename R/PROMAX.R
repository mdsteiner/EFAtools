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
    signs <- sign(colSums(AV$loadings))
    signs[signs == 0] <- 1
    AV$loadings <- AV$loadings %*% diag(signs)

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

    # this is the formula for P as given by Hendricks and White (1964)
    P <- abs(A^(k + 1)) / A

  } else if (P_type == "norm") {

    # this is the formula as used by SPSS Version 23
    P <- abs(A / sqrt(rowSums(A^2))) ^(k + 1) * (sqrt(rowSums(A^2)) / A)
  }

  # run the least squares fit
  U <- stats::lm.fit(A, P)$coefficients

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
    signs <- sign(colSums(AP))
    signs[signs == 0] <- 1
    AP <- AP %*% diag(signs)

    # order according to communalities
    eig_rotated <- diag(t(AP) %*% AP)
    eig_order <- order(eig_rotated, decreasing = TRUE)
    AP <- AP[, eig_order]

    Phi <- diag(signs) %*% Phi %*% diag(signs)
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
