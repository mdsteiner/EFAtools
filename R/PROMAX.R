#' Title
#'
#' @param L
#' @param k
#' @param do_varimax
#' @param type
#' @param P_type
#' @param kaiser
#' @param precision
#' @param order_type
#'
#' @return
#' @export
#'
#' @examples
PROMAX <- function (x, k = 4, do_varimax = TRUE, type = "EFAdiff",
                    P_type = "HW", kaiser = TRUE, precision = 1e-10,
                    order_type = "SPSS") {
  dim_names <- dimnames(x)

  # if Loadings are not yet rotated orthogonally
  if (do_varimax) {

    # perform the varimax rotation
    AV <- stats::varimax(x, normalize = kaiser, eps = precision)

    if (order_type == "SPSS") {

      # reflect factors with negative sums
      signs <- sign(colSums(AV$loadings))
      signs[signs == 0] <- 1
      AV$loadings <- AV$loadings %*% diag(signs)

      # reorder the factors according to largest sums of squares
      ss <- colSums(AV$loadings ^2)
      AV$loadings <- AV$loadings[, order(ss, decreasing = TRUE)]
    }

    # store the loading matrix in a separate object for use
    A <- AV$loadings
  } else {

    # if the entered data was already varimax rotated, bring it in the
    # appropriate form for further use
    A <- x$loadings
    AV <- x
  }

  if (P_type == "HW") {

    # this is the formula for P as given by Hendricks and White (1964)
    P <- abs(A^(k + 1)) / A
  } else if (P_type == "SPSS") {

    # this is the formula as given by SPSS Version 23
    P <- abs(A / sqrt(rowSums(A^2))) ^(k + 1) * (sqrt(rowSums(A^2)) / A)
  }

  # run the least squares fit
  U <- lm.fit(A, P)$coefficients

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

  if (order_type == "psych") {
    # reflect factors with negative sums
    signs <- sign(colSums(AP))
    signs[signs == 0] <- 1
    AP <- AP %*% diag(signs)

    # order according to eigenvalues
    ev_rotated <- diag(t(AP) %*% AP)
    ev_order <- order(ev_rotated, decreasing = TRUE)
    AP <- AP[, ev_order]

    Phi <- diag(signs) %*% Phi %*% diag(signs)
    Phi <- Phi[ev_order, ev_order]
  }


  # get structure matrix
  Structure <- AP %*% Phi

  dimnames(AP) <- dim_names
  list(loadings = AP, rotmat = U, Phi = Phi, Structure = Structure)
}
