#' Orthogonal rotations with GPArotation package
#'
#' @param x
#' @param rotation
#' @param kaiser
#' @param precision
#' @param order_type
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ROTATE_ORTH <- function(x, rotation = c("equamax", "quartimax", "geominT",
                                        "bentlerT", "bifactorT"),
                        kaiser = TRUE, precision = NULL, order_type = NULL,
                        ...){

  if (!requireNamespace("GPArotation")){

    stop("To perform the specified rotation, the GPArotation package is needed.
         Please install the GPArotation package or use varimax or promax
         rotation instead.")

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

    warning("Cannot rotate single factor. Unrotated loadings returned.")
    return(output)
  }

  # perform the requested rotation
  if(rotation == "equamax"){
  AV <- GPArotation::cfT(L, eps = precision, kappa = ncol(L)/(2 * nrow(L)),
                         normalize = kaiser, ...)

  } else if (rotation == "bentlerT"){
    AV <- GPArotation::bentlerT(L, eps = precision, normalize = kaiser, ...)

  } else if (rotation == "geominT"){
    AV <- GPArotation::geominT(L, eps = precision, normalize = kaiser, ...)

  } else if (rotation == "bifactorT"){
    AV <- GPArotation::bifactorT(L, eps = precision, normalize = kaiser, ...)

  } else {
    AV <- GPArotation::GPForth(L, method = rotation,
                             normalize = kaiser, eps = precision, ...)

  }

  # reflect factors with negative sums
  signs <- sign(colSums(AV$loadings))
  signs[signs == 0] <- 1
  AV$loadings <- AV$loadings %*% diag(signs)

  if (order_type == "ss_factors") {

    # reorder the factors according to largest sums of squares
    ss <- colSums(AV$loadings^2)
    ss_order <- order(ss, decreasing = TRUE)

    AV$loadings <- AV$loadings[, ss_order]

    AV$Th <- AV$Th[ss_order, ss_order]

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
    AV$Th <- AV$Th[eig_order, eig_order]

    if (!is.null(dim_names[[2]])) {
      dim_names[[2]] <- dim_names[[2]][eig_order]
    } else {
      dim_names[[2]] <- paste0("F", 1:ncol(AV$loadings))[eig_order]
    }

  }

  # prepare and return output list
  load_mat <- AV$loadings
  dimnames(load_mat) <- dim_names

  vars_accounted_rot <- .compute_vars(L_unrot = L, L_rot = load_mat)
  colnames(vars_accounted_rot) <- colnames(load_mat)

  # prepare output
  class(load_mat) <- "LOADINGS"

  output <- list(rot_loadings = load_mat,
                 rotmat = AV$Th,
                 vars_accounted_rot = vars_accounted_rot,
                 settings = settings)
  output
}
