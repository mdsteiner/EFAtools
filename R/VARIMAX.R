## Varimax rotation
.VARIMAX <- function (x, type = c("EFAtools", "psych", "SPSS", "none"),
                     kaiser = TRUE, precision = NULL, order_type = NULL) {

  if (type == "none") {
    # if type is none, throw an error if not
    # all the other necessary arguments are specified.

    if (is.null(precision) || is.null(order_type)) {
      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(' One of "precision", or "order_type" was NULL and no valid "type" was specified. Either use one of "EFAtools", "psych", or "SPSS" for type, or specify all other arguments\n'))
    }

    } else if (type == "EFAtools") {

      # if not specified, set PAF properties. If specified, throw warning that
      # results may not exactly match the specified type

      if (isFALSE(kaiser)) {

        warning(crayon::yellow$bold("!"), crayon::yellow(" Type and kaiser is specified. kaiser is used with value '", kaiser, "'. Results may differ from the specified type\n"))
      }

      if (is.null(precision)) {
        precision <- 1e-5
      } else {
        warning(crayon::yellow$bold("!"), crayon::yellow(" Type and precision is specified. precision is used with value '", precision, "'. Results may differ from the specified type\n"))
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

      if (is.null(precision)) {
        precision <- 1e-5
      } else {
        warning(crayon::yellow$bold("!"), crayon::yellow(" Type and precision is specified. precision is used with value '", precision, "'. Results may differ from the specified type\n"))
      }

      if (is.null(order_type)) {
        order_type <- "eigen"
      } else {
        warning(crayon::yellow$bold("!"), crayon::yellow(" Type and order_type is specified. order_type is used with value '", order_type, "'. Results may differ from the specified type\n"))
      }

    } else if (type == "SPSS") {

      # if not specified, set PAF properties. If specified, throw warning that
      # results may not exactly match the specified type

      if (isFALSE(kaiser)) {

        warning(crayon::yellow$bold("!"), crayon::yellow(" Type and kaiser is specified. kaiser is used with value '", kaiser, "'. Results may differ from the specified type\n"))
      }

      if (is.null(precision)) {
        precision <- 1e-10
      } else {
        warning(crayon::yellow$bold("!"), crayon::yellow(" Type and precision is specified. precision is used with value '", precision, "'. Results may differ from the specified type\n"))
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

  # perform the varimax rotation
  AV <- stats::varimax(L, normalize = kaiser, eps = precision)

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
