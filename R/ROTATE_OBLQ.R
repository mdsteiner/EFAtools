# Oblique factor rotations with GPArotation package
.ROTATE_OBLQ <- function(x, rotation = c("oblimin", "quartimin", "simplimax",
                                        "bentlerQ", "geominQ", "bifactorQ"),
                        type = c("EFAtools", "psych", "SPSS", "none"),
                        normalize = TRUE, precision = 1e-5, order_type = NA,
                        k = NA, randomStarts = 100, ...){

  # Fill the type preset defaults and warn about any pinned rotation arguments.
  resolved <- .resolve_settings(
    type = type,
    user = list(normalize = normalize, order_type = order_type),
    preset = .efa_presets$ROTATE_OBLQ
  )
  normalize <- resolved$normalize
  order_type <- resolved$order_type

  # extract loadings and dim names
  L <- x$unrot_loadings
  dim_names <- dimnames(L)

  # prepare settings
  settings <- list(normalize = normalize, precision = precision,
                   order_type = order_type, k = k,
                   randomStarts = randomStarts)

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
    AV <- GPArotation::bentlerQ(L, eps = precision, normalize = normalize,
                                randomStarts = randomStarts, ...)

  } else if (rotation == "oblimin"){
    AV <- GPArotation::oblimin(L, eps = precision, normalize = normalize,
                               randomStarts = randomStarts, ...)

  } else if (rotation == "quartimin"){
    AV <- GPArotation::quartimin(L, eps = precision, normalize = normalize,
                                 randomStarts = randomStarts, ...)

  } else if (rotation == "geominQ"){
    AV <- GPArotation::geominQ(L, eps = precision, normalize = normalize,
                               randomStarts = randomStarts, ...)

  } else if (rotation == "bifactorQ"){
    AV <- GPArotation::bifactorQ(L, eps = precision, normalize = normalize,
                                 randomStarts = randomStarts, ...)

  } else if (rotation == "simplimax"){

    if(is.na(k)){

      k <- nrow(L)
      # update settings
      settings$k <- k

    }

    AV <- GPArotation::simplimax(L, eps = precision, normalize = normalize,
                                 k = k, randomStarts = randomStarts, ...)

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
  dimnames(Structure) <- list(rownames(patt_mat), colnames(patt_mat))

  # compute explained variances
  vars_accounted_rot <- .compute_vars(L_unrot = L, L_rot = patt_mat, Phi = Phi)
  colnames(vars_accounted_rot) <- colnames(patt_mat)

  # prepare and return output list
  class(patt_mat) <- "LOADINGS"
  class(Structure) <- "LOADINGS"

  output <- list(rot_loadings = patt_mat,
                 Phi = Phi,
                 Structure = Structure,
                 rotmat = rotmat,
                 vars_accounted_rot = vars_accounted_rot,
                 settings = settings)
  output
}
