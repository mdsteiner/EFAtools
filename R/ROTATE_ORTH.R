# Orthogonal factor rotations with GPArotation package
.ROTATE_ORTH <- function(x, rotation = c("equamax", "quartimax", "geominT",
                                        "bentlerT", "bifactorT"),
                        type = c("EFAtools", "psych", "SPSS", "none"),
                        kaiser = TRUE, precision = NULL, order_type = NULL,
                        ...){

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
                         normalize = kaiser, ...)

  } else if (rotation == "bentlerT"){
    AV <- GPArotation::bentlerT(L, eps = precision, normalize = kaiser, ...)

  } else if (rotation == "quartimax"){
    AV <- GPArotation::bentlerT(L, eps = precision, normalize = kaiser, ...)

  } else if (rotation == "geominT"){
    AV <- GPArotation::geominT(L, eps = precision, normalize = kaiser, ...)

  } else if (rotation == "bifactorT"){
    AV <- GPArotation::bifactorT(L, eps = precision, normalize = kaiser, ...)

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

    if (!is.null(dim_names[[2]])) {
      dim_names[[2]] <- dim_names[[2]][ss_order]
    } else {
      dim_names[[2]] <- paste0("F", 1:ncol(load_mat))[ss_order]
    }

  } else if (order_type == "eigen") {

    # order according to communalities
    eig_rotated <- diag(t(load_mat) %*% load_mat)
    eig_order <- order(eig_rotated, decreasing = TRUE)
    load_mat <- load_mat[, eig_order]
    rotmat <- rotmat[eig_order, eig_order]

    if (!is.null(dim_names[[2]])) {
      dim_names[[2]] <- dim_names[[2]][eig_order]
    } else {
      dim_names[[2]] <- paste0("F", 1:ncol(load_mat))[eig_order]
    }

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
