# Oblique factor rotations with GPArotation package
.ROTATE_OBLQ <- function(x, rotation = c("oblimin", "quartimin", "simplimax",
                                        "bentlerQ", "geominQ", "bifactorQ"),
                        type = c("EFAtools", "psych", "SPSS", "none"),
                        normalize = TRUE, precision = 1e-5, order_type = NA,
                        k = NA, ...){

  if (type == "none") {

    if (is.na(order_type)) {

      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(' "order_type" was NA and no valid "type" was specified. Either use one of "EFAtools", "psych", or "SPSS" for type, or specify the "order_type" argument\n'))

    }

  } else if (type == "EFAtools") {

    if (isFALSE(normalize)) {

      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and normalize is specified. normalize is used with value '", normalize, "'. Results may differ from the specified type\n"))
    }

    if (is.na(order_type)) {
      order_type <- "eigen"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and order_type is specified. order_type is used with value '", order_type, "'. Results may differ from the specified type\n"))
    }

  } else if (type == "psych") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (isFALSE(normalize)) {

      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and normalize is specified. normalize is used with value '", normalize, "'. Results may differ from the specified type\n"))
    }

    if (is.na(order_type)) {
      order_type <- "eigen"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and order_type is specified. order_type is used with value '", order_type, "'. Results may differ from the specified type\n"))
    }

  } else if (type == "SPSS") {

    if (isFALSE(normalize)) {

      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and normalize is specified. normalize is used with value '", normalize, "'. Results may differ from the specified type\n"))
    }

    if (is.na(order_type)) {
      order_type <- "ss_factors"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and order_type is specified. order_type is used with value '", order_type, "'. Results may differ from the specified type\n"))
    }

  }

  # extract loadings and dim names
  L <- x$unrot_loadings
  dim_names <- dimnames(L)

  # prepare settings
  settings <- list(normalize = normalize, precision = precision,
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
    AV <- GPArotation::bentlerQ(L, eps = precision, normalize = normalize, ...)

  } else if (rotation == "oblimin"){
    AV <- GPArotation::oblimin(L, eps = precision, normalize = normalize, ...)

  } else if (rotation == "quartimin"){
    AV <- GPArotation::quartimin(L, eps = precision, normalize = normalize, ...)

  } else if (rotation == "geominQ"){
    AV <- GPArotation::geominQ(L, eps = precision, normalize = normalize, ...)

  } else if (rotation == "bifactorQ"){
    AV <- GPArotation::bifactorQ(L, eps = precision, normalize = normalize, ...)

  } else if (rotation == "simplimax"){

    if(is.na(k)){

      k <- nrow(L)
      # update settings
      settings$k <- k

    }

    AV <- GPArotation::simplimax(L, eps = precision, normalize = normalize,
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

    dim_names[[2]] <- dim_names[[2]][ss_order]

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

      dim_names[[2]] <- dim_names[[2]][eig_order]

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
