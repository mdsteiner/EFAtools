# Orthogonal factor rotations with GPArotation package
.ROTATE_ORTH <- function(x, rotation = c("equamax", "quartimax", "geominT",
                                        "bentlerT", "bifactorT"),
                        type = c("EFAtools", "psych", "SPSS", "none"),
                        normalize = TRUE, precision = 1e-5, order_type = NA,
                        ...){

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
                   order_type = order_type)

  if (ncol(L) < 2) {

    # prepare and return output list
    output <- list(rot_loadings = L,
                   rotmat = NA,
                   vars_accounted_rot = NA,
                   settings = settings)

    warning(crayon::yellow$bold("!"), crayon::yellow(" Cannot rotate single factor. Unrotated loadings returned.\n"))
    return(output)
  }

  # perform the requested rotation
  if(rotation == "equamax"){
  AV <- GPArotation::cfT(L, eps = precision, kappa = ncol(L)/(2 * nrow(L)),
                         normalize = normalize, ...)

  } else if (rotation == "bentlerT"){
    AV <- GPArotation::bentlerT(L, eps = precision, normalize = normalize, ...)

  } else if (rotation == "quartimax"){
    AV <- GPArotation::bentlerT(L, eps = precision, normalize = normalize, ...)

  } else if (rotation == "geominT"){
    AV <- GPArotation::geominT(L, eps = precision, normalize = normalize, ...)

  } else if (rotation == "bifactorT"){
    AV <- GPArotation::bifactorT(L, eps = precision, normalize = normalize, ...)

  }

  # Extract loading and rotation matrix
  load_mat <- AV$loadings
  dimnames(load_mat) <- dim_names
  rotmat <- AV$Th

  # reflect factors with negative sums
  signs <- sign(colSums(load_mat))
  signs[signs == 0] <- 1
  load_mat <- load_mat %*% diag(signs)

  if (order_type == "ss_factors") {

    # reorder the factors according to largest sums of squares
    ss <- colSums(load_mat^2)
    ss_order <- order(ss, decreasing = TRUE)

    load_mat <- load_mat[, ss_order]

    rotmat <- rotmat[ss_order, ss_order]

    dim_names[[2]] <- dim_names[[2]][ss_order]

  } else if (order_type == "eigen") {

    # order according to communalities
    eig_rotated <- diag(t(load_mat) %*% load_mat)
    eig_order <- order(eig_rotated, decreasing = TRUE)
    load_mat <- load_mat[, eig_order]
    rotmat <- rotmat[eig_order, eig_order]

    dim_names[[2]] <- dim_names[[2]][eig_order]

  }

  # compute explained variances
  vars_accounted_rot <- .compute_vars(L_unrot = L, L_rot = load_mat)
  colnames(vars_accounted_rot) <- colnames(load_mat)

  # prepare output and return output list
  class(load_mat) <- "LOADINGS"

  output <- list(rot_loadings = load_mat,
                 rotmat = rotmat,
                 vars_accounted_rot = vars_accounted_rot,
                 settings = settings)
  output
}
