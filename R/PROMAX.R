#' Promax rotation
#'
#' This function implements the promax rotation procedure. It can
#' reproduce the results from \code{\link[psych:fa]{psych::fa}} and the
#' SPSS FACTOR algorithm. To reproduce psych or SPSS promax results,
#' only the \code{type} argument has to be specified additional to \code{x}.
#' The other arguments can be used to control the procedure more flexibly by overriding
#' the default settings. If type = "none" is specified, all arguments with default
#' \code{NULL} have to be specified.
#'
#' @param x matrix or class \code{\link{PAF}} object. Either a matrix containing
#' an unrotated factor solution, or a \code{\link{PAF}} output object.
#' @param type character. If one of "EFAtools" (default), "psych", or "SPSS" is
#'  used, and the following arguments with default NULL are left with
#'  NULL, these implementations are executed as reported in Grieder and Steiner
#'  (2019). Individual properties can be adapted using one of the three types and
#'  specifying some of the following arguments. If set to "none" all
#'  following arguments must be specified.
#' @param kaiser logical. If \code{TRUE} (default), kaiser normalization is
#' performed in the varimax rotation.
#' @param P_type character. This specifies how the target
#'  matrix P is computed in a promax rotation. If "unnorm" it will use the
#'  unnormalized target matrix as originally done in Hendrickson and White (1964).
#'  This is also used done the psych and stats packages. If "norm" it will use the
#'  normalized target matrix as used in SPSS (see below or
#'  Grieder and Steiner, 2019 for details). Default is \code{NULL}.
#' @param precision numeric. The tolerance for stopping in the varimax procecdure.
#'  This is passed to the "eps" argument of the
#'  \code{\link[stats:varimax]{stats::varimax}} function. Default is \code{NULL}.
#' @param order_type character. How to order the factors and when to reflect
#'  their signs. "eigen" will order the factors according,
#'  to their eigenvalues. "ss_factors" will order them according to the sum of
#'  squared loadings. In both cases signs are reflected. See below for details.
#'   Default is \code{NULL}.
#' @param k numeric. The power used for computing the target matrix P in the
#'  promax rotation. Default is \code{NULL}.
#'
#' @details \code{type = "EFAtools"} will use the following argument specification:
#' \code{P_type = "unnorm", precision = 1e-10, order_type = "eigen", k = 3}.
#' \code{type = "psych"} will use the following argument specification:
#' \code{P_type = "unnorm", precision = 1e-5, order_type = "eigen", k = 4}.
#' \code{type = "SPSS"} will use the following argument specification:
#' \code{P_type = "norm", precision = 1e-10, order_type = "ss_factors", k = 4}.
#'
#' The \code{P_type} argument can take two values, "unnorm" and "norm". It controls
#' which formula is used to compute the target matrix P in the promax rotation.
#' "unnorm" uses the formula from Hendrickson and White (1964), specifically:
#' \code{P <- abs(A^(k + 1)) / A},
#' where A is the unnormalized matrix containing varimax rotated loadings.
#' "SPSS" uses the normalized varimax rotated loadings. Specifically it used the
#' following formula, which can be found in the SPSS 23 Algorithms
#' manual:
#' \code{P <- abs(A / sqrt(rowSums(A^2))) ^(k + 1) * (sqrt(rowSums(A^2)) / A)}
#'
#'#### REMOVE HERE AND IN ALL OTHER ROTATIONS
#' The \code{order_type} argument can take two arguments "eigen" or "ss_factors". The
#' order of factors is then determined using the ordered eigenvalues of the
#' promax pattern matrix, if \code{order_type = "eigen"}, and using the ordered
#' sums of squares of factor loadings per factor if \code{order_type = "ss_factors"}.
#'
#' @return A list containing the following
#'
#' \item{rot_loadings}{The promax rotated loadings (the pattern matrix).}
#' \item{rotmat}{The rotation matrix.}
#' \item{Phi}{The factor intercorrelations.}
#' \item{Structure}{The structure matrix.}
#' \item{vars_accounted_rot}{Matrix of explained variances and sums of squared
#' loadings. Based on rotated loadings and factor intercorrelations.}
#' \item{settings}{list. The settings (arguments) used in promax.}
#'
#' @export
#' @examples
#' # call within EFA function:
#' EFA(IDS2_R, n_factors = 5, type = "EFAtools", method = "PAF",
#'    rotation = "promax")
PROMAX <- function (x, type = c("EFAtools", "psych", "SPSS", "none"),
                    kaiser = TRUE, P_type = NULL, precision = NULL,
                    order_type = NULL, k = NULL) {

  if (type == "none") {
    # if type is not one of the three valid inputs, throw an error if not
    # all the other necessary arguments are specified.

    if (is.null(precision) || is.null(order_type) || is.null(k)) {
      stop('One of "P_type", "precision", "order_type", or "k" was NULL and no valid
           "type" was specified. Either use one of "EFAtools", "psych", or "SPSS"
            for type, or specify all other arguments')
    }
  } else if (type == "EFAtools") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (isFALSE(kaiser)) {

      warning("Type and kaiser is specified. kaiser is used with value '",
              kaiser, "'. Results may differ from the specified type")
    }

    if (is.null(P_type)) {
      P_type <- "unnorm"
    } else {
      warning("Type and P_type is specified. P_type is used with value '",
              P_type, "'. Results may differ from the specified type")
    }

    if (is.null(precision)) {
      precision <- 1e-5
    } else {
      warning("Type and precision is specified. precision is used with value '",
              precision, "'. Results may differ from the specified type")
    }

    if (is.null(order_type)) {
      order_type <- "eigen"
    } else {
      warning("Type and order_type is specified. order_type is used with value '",
              order_type, "'. Results may differ from the specified type")
    }

    if (is.null(k)) {
      k <- 3
    } else {
      warning("Type and k is specified. k is used with value '",
              k, "'. Results may differ from the specified type")
    }


  } else if (type == "psych") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (isFALSE(kaiser)) {

      warning("Type and kaiser is specified. kaiser is used with value '",
              kaiser, "'. Results may differ from the specified type")
    }

    if (is.null(P_type)) {
      P_type <- "unnorm"
    } else {
      warning("Type and P_type is specified. P_type is used with value '",
              P_type, "'. Results may differ from the specified type")
    }

    if (is.null(precision)) {
      precision <- 1e-5
    } else {
      warning("Type and precision is specified. precision is used with value '",
              precision, "'. Results may differ from the specified type")
    }

    if (is.null(order_type)) {
      order_type <- "eigen"
    } else {
      warning("Type and order_type is specified. order_type is used with value '",
              order_type, "'. Results may differ from the specified type")
    }

    if (is.null(k)) {
      k <- 4
    } else {
      warning("Type and k is specified. k is used with value '",
              k, "'. Results may differ from the specified type")
    }


  } else if (type == "SPSS") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (isFALSE(kaiser)) {

      warning("Type and kaiser is specified. kaiser is used with value '",
              kaiser, "'. Results may differ from the specified type")
    }

    if (is.null(P_type)) {
      P_type <- "norm"
    } else {
      warning("Type and P_type is specified. P_type is used with value '",
              P_type, "'. Results may differ from the specified type")
    }

    if (is.null(precision)) {
      precision <- 1e-10
    } else {
      warning("Type and precision is specified. precision is used with value '",
              precision, "'. Results may differ from the specified type")
    }

    if (is.null(order_type)) {
      order_type <- "ss_factors"
    } else {
      warning("Type and order_type is specified. order_type is used with value '",
              order_type, "'. Results may differ from the specified type")
    }

    if (is.null(k)) {
      k <- 4
    } else {
      warning("Type and k is specified. k is used with value '",
              k, "'. Results may differ from the specified type")
    }

  }

  # extract loadings and dim names
    L <- x$unrot_loadings
    dim_names <- dimnames(L)

  # store settings used
    settings <- list(kaiser = kaiser, P_type = P_type,
                     precision = precision, order_type = order_type, k = k)

  if (ncol(L) < 2) {
    # prepare and return output list
    output <- list(rot_loadings = L,
                   rotmat = NA,
                   Phi = NA,
                   Structure = NA,
                   vars_accounted_rot = NA,
                   settings = settings)

    warning("Cannot rotate single factor. Unrotated loadings returned.")
    return(output)
  }

  # perform the varimax rotation
  AV <- stats::varimax(L, normalize = kaiser, eps = precision)

  if (order_type == "ss_factors") {

    # reflect factors with negative sums
    signs <- sign(colSums(AV$loadings))
    signs[signs == 0] <- 1
    AV$loadings <- AV$loadings %*% diag(signs)

    # reorder the factors according to largest sums of squares
    ss <- colSums(AV$loadings ^2)
    ss_order <- order(ss, decreasing = TRUE)

    AV$loadings <- AV$loadings[, ss_order]

    AV$rotmat <- AV$rotmat[ss_order, ss_order]

    if (!is.null(dim_names[[2]])) {
      dim_names[[2]] <- dim_names[[2]][ss_order]
    } else {
      dim_names[[2]] <- paste0("F", 1:ncol(AV$loadings))[ss_order]
    }

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

    if (!is.null(dim_names[[2]])) {
      dim_names[[2]] <- dim_names[[2]][eig_order]
    } else {
      dim_names[[2]] <- paste0("F", 1:ncol(AV$loadings))[eig_order]
    }

  }

  dimnames(AP) <- dim_names

  vars_accounted_rot <- .compute_vars(L_unrot = L, L_rot = AP, Phi = Phi)

  colnames(vars_accounted_rot) <- colnames(AP)

  # get structure matrix
  Structure <- AP %*% Phi

  # prepare and return output list
  class(AP) <- "LOADINGS"

  output <- list(rot_loadings = AP,
                 rotmat = U,
                 Phi = Phi,
                 Structure = Structure,
                 vars_accounted_rot = vars_accounted_rot,
                 settings = settings)
  output
}
