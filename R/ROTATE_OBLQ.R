# Oblique factor rotations with GPArotation package
.ROTATE_OBLQ <- function(x, rotation = c("oblimin", "quartimin", "simplimax",
                                        "bentlerQ", "geominQ", "bifactorQ"),
                        type = c("EFAtools", "psych", "SPSS", "none"),
                        kaiser = TRUE, precision = NULL, order_type = NULL,
                        k = NULL, ...){

  if(is.null(precision)){

    precision <- 1e-5

  }

  if (type == "none") {

    if (is.null(order_type)) {

      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(' "order_type" was NULL and no valid "type" was specified. Either use one of "EFAtools", "psych", or "SPSS" for type, or specify the "order_type" argument\n'))

    }

  } else if (type == "EFAtools") {

    if (isFALSE(kaiser)) {

      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and kaiser is specified. kaiser is used with value '", kaiser, "'. Results may differ from the specified type\n"))
    }

    if (is.null(order_type)) {
      order_type <- "eigen"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and order_type is specified. order_type is used with value '", order_type, "'. Results may differ from the specified type\n"))
    }

  } else if (type == "psych") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (isFALSE(kaiser)) {

      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and kaiser is specified. kaiser is used with value '", kaiser, "'. Results may differ from the specified type\n"))
    }

    if (is.null(order_type)) {
      order_type <- "eigen"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and order_type is specified. order_type is used with value '", order_type, "'. Results may differ from the specified type\n"))
    }

  } else if (type == "SPSS") {

    if (isFALSE(kaiser)) {

      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and kaiser is specified. kaiser is used with value '", kaiser, "'. Results may differ from the specified type\n"))
    }

    if (is.null(order_type)) {
      order_type <- "ss_factors"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and order_type is specified. order_type is used with value '", order_type, "'. Results may differ from the specified type\n"))
    }

  }

  # extract loadings and dim names
  L <- x$unrot_loadings
  dim_names <- dimnames(L)

  # prepare settings
  settings <- list(kaiser = kaiser, precision = precision,
                   order_type = order_type, k = k)

  if (ncol(L) < 2) {

    # prepare and return output list
    output <- list(rot_loadings = L,
                   rotmat = NA,
                   Phi = NA,
                   Structure = NA,
                   vars_accounted_rot = NA,
                   settings = settings)

    warning(crayon::yellow$bold("!"), crayon::yellow(" Cannot rotate single factor. Unrotated loadings returned.\n"))
    return(output)
  }

  # perform the requested rotation
  if(rotation == "bentlerQ"){
    AV <- GPArotation::bentlerQ(L, eps = precision, normalize = kaiser, ...)

  } else if (rotation == "oblimin"){
    AV <- GPArotation::oblimin(L, eps = precision, normalize = kaiser, ...)

  } else if (rotation == "quartimin"){
    AV <- GPArotation::quartimin(L, eps = precision, normalize = kaiser, ...)

  } else if (rotation == "geominQ"){
    AV <- GPArotation::geominQ(L, eps = precision, normalize = kaiser, ...)

  } else if (rotation == "bifactorQ"){
    AV <- GPArotation::bifactorQ(L, eps = precision, normalize = kaiser, ...)

  } else if (rotation == "simplimax"){

    if(is.null(k)){

      k <- nrow(L)
      # update settings
      settings$k <- k

    }

    AV <- GPArotation::simplimax(L, eps = precision, normalize = kaiser,
                                 k = k, ...)

  }

  # get pattern and rotation matrix and factor intercorrelations
  patt_mat <- AV$loadings
  dimnames(patt_mat) <- dim_names
  rotmat <- AV$Th
  Phi <- AV$Phi

  # reflect factors with negative sums
  signs <- sign(colSums(patt_mat))
  signs[signs == 0] <- 1
  patt_mat <- patt_mat %*% diag(signs)

  if (order_type == "ss_factors") {

    # reorder the factors according to largest sums of squares
    ss <- colSums(patt_mat^2)
    ss_order <- order(ss, decreasing = TRUE)

    patt_mat <- patt_mat[, ss_order]

    rotmat <- rotmat[ss_order, ss_order]

    if (!is.null(dim_names[[2]])) {
      dim_names[[2]] <- dim_names[[2]][ss_order]
    } else {
      dim_names[[2]] <- paste0("F", 1:ncol(patt_mat))[ss_order]
    }

  }

    if (order_type == "eigen") {
      # reflect factors with negative sums
      signs <- sign(colSums(patt_mat))
      signs[signs == 0] <- 1
      patt_mat <- patt_mat %*% diag(signs)

      # order according to communalities
      eig_rotated <- diag(t(patt_mat) %*% patt_mat)
      eig_order <- order(eig_rotated, decreasing = TRUE)
      patt_mat <- patt_mat[, eig_order]

      Phi <- diag(signs) %*% Phi %*% diag(signs)
      Phi <- Phi[eig_order, eig_order]

      if (!is.null(dim_names[[2]])) {
        dim_names[[2]] <- dim_names[[2]][eig_order]
      } else {
        dim_names[[2]] <- paste0("F", 1:ncol(patt_mat))[eig_order]
      }

    }

  # get structure matrix
  Structure <- patt_mat %*% Phi

  # compute explained variances
  vars_accounted_rot <- .compute_vars(L_unrot = L, L_rot = patt_mat, Phi = Phi)
  colnames(vars_accounted_rot) <- colnames(patt_mat)

  # prepare and return output list
  class(patt_mat) <- "LOADINGS"

  output <- list(rot_loadings = patt_mat,
                 Phi = Phi,
                 Structure = Structure,
                 rotmat = rotmat,
                 vars_accounted_rot = vars_accounted_rot,
                 settings = settings)
  output
}
