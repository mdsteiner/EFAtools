#' Promax rotation
#'
#' This function implements the promax rotation procedure. It can
#' reproduce the results from \code{\link[psych:fa]{psych::fa}} and the
#' SPSS FACTOR algorithm. To reproduce psych or SPSS promax results,
#' only the \code{type} argument has to be specified additional to \code{x}.
#' The other arguments can be used to control the procedure more flexibly.
#'
#' @param x matrix or class \code{\link{PAF}} object. Either a matrix containing
#' an unrotated factor solution, or a \code{\link{PAF}} output object.
#' @param k numeric. The power used for computing the target matrix P in the
#'  promax rotation.
#' @param type character. If one of "EFAdiff" (default), "psych", or "SPSS" is
#'  used, and the following arguments (except kaiser) are left with \code{NULL},
#'  these implementations
#'  are executed as reported in Steiner and Grieder (2019; see details).
#'  Individual properties can be adapted using one of the three types and
#'  specifying some of the following
#'  arguments. If set to another value than one of the three specified above, all
#'  following arguments must be specified.
#' @param kaiser logical. If \code{TRUE} (default), kaiser normalization is
#' performed in the varimax rotation.
#' @param P_type character. Default is \code{NULL} This specifies how the target
#'  matrix P is computed. If "HW" it will use the implementation of Hendrickson and
#'  White (1964), which is also used by the psych and stats packages. If "SPSS"
#'  it will use the SPSS version (see below or Steiner and Grieder, 2019 for details).
#' @param precision numeric. The tolerance for stopping in the varimax procecdure.
#'  This is passed to the "eps" argument of the
#'  \code{\link[stats:varimax]{stats::varimax}} function. Default is \code{NULL}
#' @param order_type character. How to order the factors and when to reflect
#'  their signs. Default is \code{NULL}. "psych" will use the psych method, "SPSS" the
#'  SPSS method. See below for details.
#'
#' @details \code{type = "EFAdiff"} will use the following argument specification:
#' \code{P_type = "HW", precision = 1e-10, order_type = "psych"}.
#' \code{type = "psych"} will use the following argument specification:
#' \code{P_type = "HW", precision = 1e-5, order_type = "psych"}.
#' \code{type = "SPSS"} will use the following argument specification:
#' \code{P_type = "SPSS", precision = 1e-10, order_type = "SPSS"}.
#'
#' The \code{P_type} argument can take two values, "HW" and "SPSS". It controlls
#' which formula is used to compute the target matrix P in the promax rotation.
#' "HW" uses the formula from Hendrickson and White (1964), specifically:
#' \code{P <- abs(A^(k + 1)) / A},
#' where A is the matrix containing varimax rotated loadings.
#' "SPSS" uses a different formula, which can be found in the SPSS 23 Algorithms
#' manual:
#' \code{P <- abs(A / sqrt(rowSums(A^2))) ^(k + 1) * (sqrt(rowSums(A^2)) / A)}
#'
#' The \code{order_type} argument can take two arguments "psych" or "SPSS". The
#' order of factors is then determined using the ordered communalities of the
#' promax pattern matrix, if \code{order_type = "psych"}, and using the ordered
#' sums of squares of factor loadings per factor if \code{order_type = "SPSS"}.
#'loadings = AP, rotmat = U, Phi = Phi, Structure = Structure
#' @return A list of class PROMAX containing the following
#'
#' \item{loadings}{The promax rotated loadings (the pattern matrix).}
#' \item{rotmat}{The rotation matrix.}
#' \item{Phi}{The factor intercorrelations.}
#' \item{Structure}{The structure matrix.}
#'
#' @export
PROMAX <- function (x, k = 4, type = "EFAdiff", kaiser = TRUE, P_type = NULL,
                    precision = NULL, order_type = NULL) {

  if (is.null(type) || !(type %in% c("EFAdiff", "psych", "SPSS"))) {
    # if type is not one of the three valid inputs, throw an error if not
    # all the other necessary arguments are specified.

    if (is.null(P_type) || is.null(precision) || is.null(order_type)) {
      stop('One of "P_type", "precision", or "order_type" was NULL and no valid
           "type" was specified. Either use one of "EFAdiff", "psych", or "SPSS"
            for type, or specify all other arguments')
    }
  } else if (type == "EFAdiff") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (isFALSE(kaiser)) {

      warning("Type and kaiser is specified. kaiser is used with value '",
              kaiser, "'. Results may differ from the specified type")
    }

    if (is.null(P_type)) {
      P_type <- "HW"
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
      order_type <- "psych"
    } else {
      warning("Type and order_type is specified. order_type is used with value '",
              order_type, "'. Results may differ from the specified type")
    }


  } else if (type == "psych") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (isFALSE(kaiser)) {

      warning("Type and kaiser is specified. kaiser is used with value '",
              kaiser, "'. Results may differ from the specified type")
    }

    if (is.null(P_type)) {
      P_type <- "HW"
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
      order_type <- "psych"
    } else {
      warning("Type and order_type is specified. order_type is used with value '",
              order_type, "'. Results may differ from the specified type")
    }


  } else if (type == "SPSS") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (isFALSE(kaiser)) {

      warning("Type and kaiser is specified. kaiser is used with value '",
              kaiser, "'. Results may differ from the specified type")
    }

    if (is.null(P_type)) {
      P_type <- "SPSS"
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
      order_type <- "SPSS"
    } else {
      warning("Type and order_type is specified. order_type is used with value '",
              order_type, "'. Results may differ from the specified type")
    }

  }

  # extract loadings and dim names
  if (all(class(x) == "PAF")) {
    L <- x$loadings
    dim_names <- dimnames(L)
  } else if (all(class(x) == "matrix") |
             any(class(x) %in% c("loadings", "LOADINGS"))) {
    L <- x
    dim_names <- dimnames(L)
  } else {
    stop("x is not of class PAF and not a matrix. Either provide a PAF output
         object, or a matrix containing unrotated factor loadings")
  }

  # perform the varimax rotation
  AV <- stats::varimax(L, normalize = kaiser, eps = precision)

  if (order_type == "SPSS") {

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

  if (P_type == "HW") {

    # this is the formula for P as given by Hendricks and White (1964)
    P <- abs(A^(k + 1)) / A
  } else if (P_type == "SPSS") {

    # this is the formula as given by SPSS Version 23
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

  if (order_type == "psych") {
    # reflect factors with negative sums
    signs <- sign(colSums(AP))
    signs[signs == 0] <- 1
    AP <- AP %*% diag(signs)

    # order according to communalities
    comm_rotated <- diag(t(AP) %*% AP)
    comm_order <- order(comm_rotated, decreasing = TRUE)
    AP <- AP[, comm_order]

    Phi <- diag(signs) %*% Phi %*% diag(signs)
    Phi <- Phi[comm_order, comm_order]

    if (!is.null(dim_names[[2]])) {
      dim_names[[2]] <- dim_names[[2]][comm_order]
    } else {
      dim_names[[2]] <- paste0("F", 1:ncol(AV$loadings))[comm_order]
    }

  }

  dimnames(AP) <- dim_names


  # compute variance proportions
  vars <- diag(Phi %*% t(AP) %*% AP)

  # Compute the explained variances. The code is based on the psych::fac() function
  # total variance (sum of communalities and uniquenesses)
  h2 <- diag(L %*% t(L))
  var_total <- sum(h2 + (1 - h2))
  vars_explained <- rbind(`SS loadings` = vars)
  vars_explained <- rbind(vars_explained, `Proportion Var` = vars / var_total)

  if (ncol(L) > 1) {
    vars_explained <- rbind(vars_explained,
                            `Cumulative Var` = cumsum(vars / var_total))
    vars_explained <- rbind(vars_explained,
                            `Prop of Explained Var` = vars / sum(vars))
    vars_explained <- rbind(vars_explained,
                            `Cum Prop of Explained Var` = cumsum(vars / sum(vars)))
  }
  vars_accounted <- vars_explained

  colnames(vars_accounted) <- colnames(AP)

  # get structure matrix
  Structure <- AP %*% Phi

  # prepare and return output list
  class(AP) <- "LOADINGS"
  output <- list(loadings = AP, rotmat = U, Phi = Phi, Structure = Structure,
                 h2 = h2, vars_accounted = vars_accounted)
  class(output) <- "PROMAX"
  output
}
