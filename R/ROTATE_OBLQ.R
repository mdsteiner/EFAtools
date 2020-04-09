ROTATE_OBLQ <- function(x, rotation = c("oblimin", "quartimin", "bentlerQ",
                                        "geominQ", "bifactorQ"),
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
                   Phi = NA,
                   vars_accounted_rot = NA,
                   settings = settings)

    warning("Cannot rotate single factor. Unrotated loadings returned.")
    return(output)
  }

  # perform the requested rotation
  if(rotation == "bentlerQ"){
    AV <- GPArotation::bentlerQ(L, eps = precision,
                                kappa = ncol(L)/(2 * nrow(L)),
                                normalize = kaiser, ...)

  } else if (rotation == "geominQ"){
    AV <- GPArotation::geominQ(L, eps = precision, normalize = kaiser, ...)

  } else if (rotation == "bifactorQ"){
    AV <- GPArotation::bifactorQ(L, eps = precision, normalize = kaiser, ...)

  } else {
    AV <- GPArotation::GPForth(L, method = rotation,
                               normalize = kaiser, eps = precision, ...)

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
                 rotmat = rotmat,
                 Phi = Phi,
                 Structure = Structure,
                 vars_accounted_rot = vars_accounted_rot,
                 settings = settings)
  output
}
