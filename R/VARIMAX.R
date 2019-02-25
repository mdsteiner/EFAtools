#' Varimax rotation
#'
#' This function implements the varimax rotation procedure. It can
#' reproduce the results from \code{\link[psych:fa]{psych::fa}} and the
#' SPSS FACTOR algorithm. To reproduce psych or SPSS varimax rotation,
#' only the \code{type} argument has to be specified additional to \code{x}.
#' The other arguments can be used to control the procedure more flexibly.
#'
#' @param x matrix or class \code{\link{PAF}} object. Either a matrix containing
#' an unrotated factor solution, or a \code{\link{PAF}} output object.
#' @param type character. If one of "SG" (default), "psych", or "SPSS" is
#'  used, and the following arguments (except kaiser) are left with \code{NULL},
#'  these implementations
#'  are executed as reported in Steiner and Grieder (2019; see details).
#'  Individual properties can be adapted using one of the three types and
#'  specifying some of the following
#'  arguments. If set to another value than one of the three specified above, all
#'  following arguments must be specified.
#' @param kaiser logical. If \code{TRUE} (default), kaiser normalization is
#'  performed in the varimax rotation.
#' @param precision numeric. The tolerance for stopping in the varimax procecdure.
#'  This is passed to the "eps" argument of the
#'  \code{\link[stats:varimax]{stats::varimax}} function. Default is \code{NULL}
#' @param order_type character. How to order the factors and when to reflect
#'  their signs. Default is \code{NULL}. "eigen" will order the factors according,
#'  to their eigenvalues. "ss_factors" will order them according to the sum of
#'  squared loadings. In both cases signs are reflected. See below for details.
#'
#' @details \code{type = "SG"} will use the following argument specification:
#' \code{precision = 1e-10, order_type = "eigen"}.
#' \code{type = "psych"} will use the following argument specification:
#' \code{precision = 1e-5, order_type = "eigen"}.
#' \code{type = "SPSS"} will use the following argument specification:
#' \code{precision = 1e-10, order_type = "ss_factors"}.
#'
#' The \code{order_type} argument can take two arguments "eigen" or "ss_factors". The
#' order of factors is then determined using the ordered communalities of the
#' promax pattern matrix, if \code{order_type = "eigen"}, and using the ordered
#' sums of squares of factor loadings per factor if \code{order_type = "ss_factors"}.
#'loadings = AP, rotmat = U, Phi = Phi, Structure = Structure
#' @return A list of class VARIMAX containing the following
#'
#' \item{rot_loadings}{The varimax rotated loadings (the pattern matrix).}
#' \item{rotmat}{The rotation matrix.}
#' \item{h2}{The communalities from the unrotated solution.}
#' \item{vars_accounted}{Matrix of explained variances and sums of squared loadings}
#' \item{fit_indices}{Fit indices as returned by
#'  \code{\link[psych:factor.stats]{psych::factor.stats}}}
#' \item{settings}{list. The settings (arguments) used in the varimax.}
#'
#' @export
#' @examples
#' # call within EFA function:
#' EFA(IDS2_R, n_factors = 5, type = "SG", rotation = "varimax")
#'
#' # call as single function
#' unrot <- EFA(IDS2_R, n_factors = 5, type = "SG")
#' VARIMAX(unrot$unrot_loadings, type = "SG")
VARIMAX <- function (x, type = "SG", kaiser = TRUE,
                    precision = NULL, order_type = NULL) {

  if (is.null(type) || !(type %in% c("SG", "psych", "SPSS"))) {
    # if type is not one of the three valid inputs, throw an error if not
    # all the other necessary arguments are specified.

    if (is.null(precision) || is.null(order_type)) {
      stop('One of "precision", or "order_type" was NULL and no valid
           "type" was specified. Either use one of "SG", "psych", or "SPSS"
           for type, or specify all other arguments')
      }
    } else if (type == "SG") {

      # if not specified, set PAF properties. If specified, throw warning that
      # results may not exactly match the specified type

      if (isFALSE(kaiser)) {

        warning("Type and kaiser is specified. kaiser is used with value '",
                kaiser, "'. Results may differ from the specified type")
      }

      if (is.null(precision)) {
        precision <- 1e-10
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


    } else if (type == "psych") {

      # if not specified, set PAF properties. If specified, throw warning that
      # results may not exactly match the specified type

      if (isFALSE(kaiser)) {

        warning("Type and kaiser is specified. kaiser is used with value '",
                kaiser, "'. Results may differ from the specified type")
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

    } else if (type == "SPSS") {

      # if not specified, set PAF properties. If specified, throw warning that
      # results may not exactly match the specified type

      if (isFALSE(kaiser)) {

        warning("Type and kaiser is specified. kaiser is used with value '",
                kaiser, "'. Results may differ from the specified type")
      }

      if (is.null(precision)) {
        precision <- 1e-10
      } else {
        warning("Type and precision is specified. precision is used with value '",
                precision, "'. Results may differ from the specified type")
      }

      if (is.null(order_type)) {
        order_type <- "ss_loadings"
      } else {
        warning("Type and order_type is specified. order_type is used with value '",
                order_type, "'. Results may differ from the specified type")
      }

    }

  # extract loadings and dim names
  if (all(class(x) == "PAF")) {
    L <- x$unrot_loadings
    dim_names <- dimnames(L)

    # N <- x$fit_indices$n.obs

  } else if (all(class(x) == "matrix") |
             any(class(x) %in% c("loadings", "LOADINGS"))) {
    L <- x
    dim_names <- dimnames(L)

  } else {
    stop("x is not of class PAF and not a matrix. Either provide a PAF output
         object, or a matrix containing unrotated factor loadings")
  }

  if (ncol(L) < 2) {
    if (any(class(x) == "PAF")) {
      vars_accounted <- x$vars_accounted
      fit_ind <- x$fit_indices
    } else {
      vars_accounted <- NA
      fit_ind <- NA
    }
    # prepare and return output list
    output <- list(rot_loadings = L, rotmat = NA,
                   h2 = diag(L %*% t(L)),
                   vars_accounted = vars_accounted,
                   fit_indices = fit_ind)
    class(output$h2) <- "COMMUNALITIES"
    class(output) <- "VARIMAX"

    warning("Cannot rotate single factor. Unrotated loadings returned.")
    return(output)
  }

  # perform the varimax rotation
  if (type != "SPSS") {
    AV <- stats::varimax(L, normalize = kaiser, eps = precision)
  } else {
    AV <- VARIMAX_SPSS(L, normalize = kaiser, precision = precision)
  }


  # reflect factors with negative sums
  signs <- sign(colSums(AV$loadings))
  signs[signs == 0] <- 1
  AV$loadings <- AV$loadings %*% diag(signs)

  if (order_type == "ss_factors") {

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

  } else if (order_type == "eigen") {

    # order according to communalities
    eig_rotated <- diag(t(AV$loadings) %*% AV$loadings)
    eig_order <- order(eig_rotated, decreasing = TRUE)
    AV$loadings <- AV$loadings[, eig_order]
    AV$rotmat <- AV$rotmat[eig_order, eig_order]

    if (!is.null(dim_names[[2]])) {
      dim_names[[2]] <- dim_names[[2]][eig_order]
    } else {
      dim_names[[2]] <- paste0("F", 1:ncol(AV$loadings))[eig_order]
    }

  }

  # prepare and return output list
  load_mat <- AV$loadings
  dimnames(load_mat) <- dim_names

  vars_accounted <- .compute_vars(L_unrot = L, L_rot = load_mat)
  colnames(vars_accounted) <- colnames(load_mat)

  if (any(class(x) == "PAF") && all(class(x$fit_indices) == c("psych", "stats"))) {
    # compute fit indices
    fit_ind <- psych::factor.stats(f = load_mat, phi = NULL, r = x$orig_R,
                                   n.obs = x$fit_indices$n.obs)

    if(all(class(x) == "try_error")) {
      fit_ind <- NA
    }

  } else {
    fit_ind <- NA
  }

  # prepare output
  class(load_mat) <- "LOADINGS"

  # prepare settings
  settings <- list(type = type, kaiser = kaiser, precision = precision,
                   order_type = order_type)
  output <- list(rot_loadings = load_mat,
                 rotmat = AV$rotmat,
                 h2 = diag(L %*% t(L)),
                 vars_accounted = vars_accounted,
                 fit_indices = fit_ind,
                 settings = settings)
  class(output$h2) <- "COMMUNALITIES"
  class(output) <- "VARIMAX"
  output
}
