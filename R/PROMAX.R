## Promax rotation
.PROMAX <- function (x, type = c("EFAtools", "psych", "SPSS", "none"),
                    normalize = TRUE, P_type = NA, precision = 1e-5,
                    order_type = NA,  varimax_type = NA, k = NA) {

  if (type == "none") {
    # if type is none, throw an error if not
    # all the other necessary arguments are specified.

    if (is.na(P_type) || is.na(order_type) || is.na(varimax_type) || is.na(k)) {
      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(' One of "P_type", "order_type", "varimax_type", or "k" was NA and no valid "type" was specified. Either use one of "EFAtools", "psych", or "SPSS" for type, or specify all other arguments\n'))
    }
  } else if (type == "EFAtools") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (isFALSE(normalize)) {

      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and normalize is specified. normalize is used with value '", normalize, "'. Results may differ from the specified type\n"))
    }

    if (is.na(P_type)) {
      P_type <- "norm"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and P_type is specified. P_type is used with value '", P_type, "'. Results may differ from the specified type\n"))
    }

    if (is.na(order_type)) {
      order_type <- "eigen"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and order_type is specified. order_type is used with value '", order_type, "'. Results may differ from the specified type\n"))
    }

    if (is.na(varimax_type)) {
      varimax_type <- "svd"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and varimax_type is specified. varimax_type is used with value '", varimax_type, "'. Results may differ from the specified type\n"))
    }

    if (is.na(k)) {
      k <- 4
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and k is specified. k is used with value '", k, "'. Results may differ from the specified type\n"))
    }


  } else if (type == "psych") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (isFALSE(normalize)) {

      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and normalize is specified. normalize is used with value '", normalize, "'. Results may differ from the specified type.\n"))
    }

    if (is.na(P_type)) {
      P_type <- "unnorm"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and P_type is specified. P_type is used with value '", P_type, "'. Results may differ from the specified type\n"))
    }

    if (is.na(order_type)) {
      order_type <- "eigen"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and order_type is specified. order_type is used with value '", order_type, "'. Results may differ from the specified type\n"))
    }

    if (is.na(varimax_type)) {
      varimax_type <- "svd"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and varimax_type is specified. varimax_type is used with value '", varimax_type, "'. Results may differ from the specified type\n"))
    }

    if (is.na(k)) {
      k <- 4
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and k is specified. k is used with value '", k, "'. Results may differ from the specified type\n"))
    }


  } else if (type == "SPSS") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (isFALSE(normalize)) {

      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and normalize is specified. normalize is used with value '", normalize, "'. Results may differ from the specified type.\n"))
    }

    if (is.na(P_type)) {
      P_type <- "norm"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and P_type is specified. P_type is used with value '", P_type, "'. Results may differ from the specified type\n"))
    }

    if (is.na(order_type)) {
      order_type <- "ss_factors"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and order_type is specified. order_type is used with value '", order_type, "'. Results may differ from the specified type\n"))
    }

    if (is.na(varimax_type)) {
      varimax_type <- "kaiser"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and varimax_type is specified. varimax_type is used with value '", varimax_type, "'. Results may differ from the specified type\n"))
    }

    if (is.na(k)) {
      k <- 4
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and k is specified. k is used with value '", k, "'. Results may differ from the specified type\n"))
    }

  }

  # extract loadings and dim names
    L <- x$unrot_loadings
    dim_names <- dimnames(L)

  # store settings used
    settings <- list(normalize = normalize, P_type = P_type,
                     precision = precision, order_type = order_type,
                     varimax_type = varimax_type, k = k)

  if (ncol(L) < 2) {
    # prepare and return output list
    output <- list(rot_loadings = L,
                   Phi = NA,
                   Structure = NA,
                   rotmat = NA,
                   vars_accounted_rot = NA,
                   settings = settings)

    warning(crayon::yellow$bold("!"), crayon::yellow(" Cannot rotate single factor. Unrotated loadings returned.\n"))
    return(output)
  }

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

    AV$rotmat <- AV$rotmat[ss_order, ss_order]

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

    dim_names[[2]] <- dim_names[[2]][eig_order]

  }

  dimnames(AP) <- dim_names

  vars_accounted_rot <- .compute_vars(L_unrot = L, L_rot = AP, Phi = Phi)

  colnames(vars_accounted_rot) <- colnames(AP)

  # get structure matrix
  Structure <- AP %*% Phi

  # prepare and return output list
  class(AP) <- "LOADINGS"

  output <- list(rot_loadings = AP,
                 Phi = Phi,
                 Structure = Structure,
                 rotmat = U,
                 vars_accounted_rot = vars_accounted_rot,
                 settings = settings)
  output
}
