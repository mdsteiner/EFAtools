## Varimax rotation
.VARIMAX <- function (x, type = c("EFAtools", "psych", "SPSS", "none"),
                     normalize = TRUE, precision = 1e-5, order_type = NA,
                     varimax_type = NA) {

  if (type == "none") {
    # if type is none, throw an error if not
    # all the other necessary arguments are specified.

    if (is.na(order_type) || is.na(varimax_type)) {
      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(' "order_type" or "varimax_type" was NA and no valid "type" was specified. Either use one of "EFAtools", "psych", or "SPSS" for type, or specify all other arguments\n'))
    }

    } else if (type == "EFAtools") {

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

      if (is.na(varimax_type)) {
        varimax_type <- "svd"
      } else {
        warning(crayon::yellow$bold("!"), crayon::yellow(" Type and varimax_type is specified. varimax_type is used with value '", varimax_type, "'. Results may differ from the specified type\n"))
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

      if (is.na(varimax_type)) {
        varimax_type <- "svd"
      } else {
        warning(crayon::yellow$bold("!"), crayon::yellow(" Type and varimax_type is specified. varimax_type is used with value '", varimax_type, "'. Results may differ from the specified type\n"))
      }

    } else if (type == "SPSS") {

      # if not specified, set PAF properties. If specified, throw warning that
      # results may not exactly match the specified type

      if (isFALSE(normalize)) {

        warning(crayon::yellow$bold("!"), crayon::yellow(" Type and normalize is specified. normalize is used with value '", normalize, "'. Results may differ from the specified type\n"))
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

    }

    # extract loadings and dim names
    L <- x$unrot_loadings
    dim_names <- dimnames(L)

    # prepare settings
    settings <- list(normalize = normalize, precision = precision,
                     order_type = order_type, varimax_type = varimax_type)

  if (ncol(L) < 2) {

    # prepare and return output list
    output <- list(rot_loadings = L,
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

    dim_names[[2]] <- dim_names[[2]][ss_order]

  } else if (order_type == "eigen") {

    # order according to communalities
    eig_rotated <- diag(t(AV$loadings) %*% AV$loadings)
    eig_order <- order(eig_rotated, decreasing = TRUE)
    AV$loadings <- AV$loadings[, eig_order]
    AV$rotmat <- AV$rotmat[eig_order, eig_order]

    dim_names[[2]] <- dim_names[[2]][eig_order]

  }

  # prepare and return output list
  load_mat <- AV$loadings
  dimnames(load_mat) <- dim_names

  vars_accounted_rot <- .compute_vars(L_unrot = L, L_rot = load_mat)
  colnames(vars_accounted_rot) <- colnames(load_mat)

  # prepare output
  class(load_mat) <- "LOADINGS"

  output <- list(rot_loadings = load_mat,
                 rotmat = AV$rotmat,
                 vars_accounted_rot = vars_accounted_rot,
                 settings = settings)
  output
}


# SPSS varimax implementation as described in the SPSS manual: ftp://public.dhe.ibm.com/software/analytics/spss/documentation/statistics/23.0/en/client/Manuals/IBM_SPSS_Statistics_Algorithms.pdf

.VARIMAX_SPSS <- function(x, normalize = TRUE, precision = 1e-5) {

  # get dimensions
  n_cols <- ncol(x)
  n_rows <- nrow(x)

  if (isTRUE(normalize)) {
    h2 <- diag(x %*% t(x))
    x <- diag(1/sqrt(h2)) %*% x
  }


  # initialize rotation matrix as identity matrix
  rotmat <- diag(1, n_cols)

  SV_now <- .SV(x)


  for (it in 1:1000) {

    for (col_j in 1:(n_cols - 1)) {
      for (col_k in (col_j + 1):n_cols) {

        u_p <- x[,col_j] ** 2 - x[,col_k] ** 2
        v_p <- 2 * x[,col_j] * x[,col_k]

        A <- sum(u_p)
        B <- sum(v_p)
        C <- sum(u_p ** 2 - v_p ** 2)
        D <- sum(2 * u_p * v_p)

        X <- D - 2 * A * B / n_rows
        Y <- C - (A ** 2 - B ** 2) / n_rows

        P <- 1/4 * atan2(X, Y)

        if (abs(sin(P)) > 1e-15) {
          sub_rot <- matrix(c(cos(P), -sin(P),
                              sin(P), cos(P)),
                            ncol = 2, byrow = TRUE)

          x[, c(col_j, col_k)] <- x[, c(col_j, col_k)] %*% sub_rot
          rotmat[, c(col_j, col_k)] <- rotmat[, c(col_j, col_k)] %*% sub_rot
        }

      }
    }

    SV_old <- SV_now
    SV_now <- .SV(x)

    if (abs(SV_now - SV_old) <= precision) {
      break
    }

  }

  if (isTRUE(normalize)) {
    x <- diag(sqrt(h2)) %*% x
  }

  # Create output list
  list(rotmat = rotmat, loadings = x, iter = it)

}
