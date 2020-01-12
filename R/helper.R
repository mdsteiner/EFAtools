#' Format numbers for print method
#'
#' Helper function used in the print method for class LOADINGS and SLLOADINGS.
#' Strips the 0 in front of the decimal point of a number if number < 1, only
#' keeps the first \code{digits} number of digits, and adds an empty space in
#' front of the number if the number is positive. This way all returned strings
#' (except for those > 1, which are exceptions in LOADINGS) have the same number
#' of characters.
#'
#' @param x numeric. Number to be formated.
#' @param digits numeric. Number of digits after the comma to keep.
#' @param print_zero logical. Whether, if a number is between ]-1, 1[, the
#'  zero should be omitted or printed (default is FALSE, i.e. omit zeros).
#'
#' @return A formated number
.numformat <- function(x, digits = 2, print_zero = FALSE) {

  if (isFALSE(print_zero)) {

    ncode <- paste0("%.", digits, "f")
    x <- sub("^(-?)0.", "\\1.", sprintf(ncode, x))
    x <- stringr::str_pad(x, digits + 2, "left")

  } else {

    ncode <- paste0("%.", digits, "f")
    x <- sprintf(ncode, x)
    x <- stringr::str_pad(x, digits + 3, "left")

  }

}


#' Compute explained variances from loadings
#'
#' From unrotated loadings compute the communalities and uniquenesses for total
#' variance. Compute explained variances per factor from rotated loadings (and
#' factor intercorrelations Phi if oblique rotation was used).
#'
#' @param L_unrot matrix. Unrotated factor loadings.
#' @param L_rot matrix. Rotated factor loadings.
#' @param Phi matrix. Factor intercorrelations. Provide only if oblique rotation
#'  is used.
#'
#' @return A matrix with sum of squared loadings, proportion explained variance
#'  from total variance per factor, same as previous but cumulative, Proportion
#'  of explained variance from total explained variance, and same as previous but
#'  cumulative.
.compute_vars <- function(L_unrot, L_rot, Phi = NULL) {

  if (is.null(Phi)) {
    # compute variance proportions
    if (ncol(L_rot) > 1) {
      vars <- colSums(L_rot^2)
    }
    else {
      vars <- sum(L_rot^2)
    }
  } else {
    # compute variance proportions
    vars <- diag(Phi %*% t(L_rot) %*% L_rot)
  }

  # Compute the explained variances. The code is based on the psych::fac() function
  # total variance (sum of communalities and uniquenesses)
  h2 <- diag(L_unrot %*% t(L_unrot))
  var_total <- sum(h2 + (1 - h2))
  vars_explained <- rbind(`SS loadings` = vars)
  vars_explained <- rbind(vars_explained, `Prop Tot Var` = vars / var_total)

  if (ncol(L_rot) > 1) {
    vars_explained <- rbind(vars_explained,
                            `Cum Prop Tot Var` = cumsum(vars / var_total))
    vars_explained <- rbind(vars_explained,
                            `Prop Comm Var` = vars / sum(vars))
    vars_explained <- rbind(vars_explained,
                            `Cum Prop Comm Var` = cumsum(vars / sum(vars)))
  }

  vars_explained
}



.SV <- function(lambda) {

  n <- nrow(lambda)

  sum((n * colSums(lambda ** 4) - colSums(lambda ** 2) ** 2)/ n ** 2)

}


.get_compare_matrix <- function(x, digits = 3, r_red = .001, n_char = 10) {

  # create factor names to display
  factor_names <- colnames(x)
  if (is.null(factor_names)) {
    factor_names <- paste0("F", 1:ncol(x))
  }

  # for equal spacing, fill the factor names such that they match the columns
  fn_nchar <- sapply(factor_names, nchar)
  factor_names[which(fn_nchar > digits + 2)] <- substr(
    factor_names[which(fn_nchar > digits + 2)] , 1, digits + 2)
  factor_names <- stringr::str_pad(factor_names, digits + 2, side = "both")

  var_names <- rownames(x)
  if (is.null(var_names)) {
    var_names <- paste0("V", 1:nrow(x))
  }

  max_char <- max(sapply(var_names, nchar))

  if (max_char > n_char) {
    vn_nchar <- sapply(var_names, nchar)
    var_names[which(vn_nchar > n_char)] <- substr(var_names[which(n_char)] ,
                                              1, n_char)
    max_char <- n_char
  }

  var_names <- stringr::str_pad(var_names, max_char, side = "right")

  n_col <- ncol(x)

  # create the string to paste using the crayon package
  temp <- apply(matrix(1:nrow(x), ncol = 1), 1,
                function(ind, x, cutoff, n_col, vn, digits){
                  i <- x[ind,]

                  tt <- crayon::blue(vn[ind])

                  for (kk in 1:n_col) {
                    if (abs(i[kk]) <= cutoff) {
                      tt <- c(tt, .numformat(round(i[kk], digits = digits),
                                                          digits = digits,
                                             print_zero = TRUE))
                    } else {
                      tt <- c(tt,
                              crayon::red(.numformat(round(i[kk],
                                                                digits = digits),
                                                              digits = digits,
                                                              print_zero = TRUE)))
                    }
                  }
                  stringr::str_c(tt, collapse = "\t")
                }, cutoff = r_red, n_col = n_col, digits = digits, x = x,
                vn = var_names)

  factor_names <- stringr::str_c(factor_names,
                                 collapse = "\t")
  factor_names <- crayon::blue(stringr::str_c( stringr::str_pad(" ", max_char),
                                               "\t", factor_names))


  temp <- stringr::str_c(temp, collapse = "\n")

  temp <- stringr::str_c(factor_names, "\n", temp)


  temp <- stringr::str_c(temp, "\n")
  # print the results to the console

  temp
}


.get_compare_vector <- function(x, digits = 3, r_red = .001) {

  temp_i <- NULL

  for (ii in seq_along(x)) {
    if (abs(x[ii]) > r_red) {
      temp_i <- c(temp_i, crayon::red(.numformat(round(x[ii], digits = digits),
                                                      digits = digits,
                                                      print_zero = TRUE)))
    } else {
      temp_i <- c(temp_i, .numformat(round(x[ii], digits = digits),
                                     digits = digits,
                                     print_zero = TRUE))
    }
  }

  for (ss in seq(1, length(x), 7)) {
    if (length(x) > ss + 6) {
      tt <- ss + 6
    } else {
      tt <- length(x)
    }
    if (ss == 1) {
      temp <- stringr::str_c(temp_i[ss:tt], collapse = "  ")
    } else {
      temp <- stringr::str_c(temp, "\n", stringr::str_c(temp_i[ss:tt],
                                                        collapse = "  "))
    }

  }

  temp <- stringr::str_c(temp, "\n")

  # print the results to the console
  return(temp)
}



.decimals <- function(x) {

  if ((is.null(dim(x)) && !(class(x) %in% c("numeric", "COMMUNALITIES", "EIGEN"))) ||
      (!is.null(dim(x)) && !(class(x) %in% c("matrix", "loadings", "LOADINGS",
                                             "SLLOADINGS")))) {
    stop("x is of class ", class(x), " but must be a numeric vector or matrix")
  }

  if (!is.null(dim(x))) {

    max(apply(x, 1:2, function(ll) {
      if (abs(ll - round(ll)) > .Machine$double.eps^0.5) {
        nchar(strsplit(sub('0+$', '', as.character(ll)), ".", fixed = TRUE)[[1]][[2]])
      } else {
        return(0)
      }
    }))

  } else if (length(x) > 1) {

    max(sapply(x, function(ll) {
      if (abs(ll - round(ll)) > .Machine$double.eps^0.5) {
        nchar(strsplit(sub('0+$', '', as.character(ll)), ".", fixed = TRUE)[[1]][[2]])
      } else {
        return(0)
      }
    }))


  } else {

    if (abs(x - round(x)) > .Machine$double.eps^0.5) {
      nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
    } else {
      return(0)
    }

  }

}


.factor_congruence <- function(x, y, digits = 3, na.rm = TRUE) {

  if (any(is.na(x) | any(is.na(y)))) {
    warning("Some loadings were missing.")
    if (isTRUE(na.rm)) {
      message("Analysis is performed on complete cases")
      if (any(is.na(x))) {
        xc <- x[stats::complete.cases(x), ]
        y <- y[stats::complete.cases(x), ]
        x <- xc
      }
      if (any(is.na(y))) {
        yc <- y[stats::complete.cases(y), ]
        x <- x[stats::complete.cases(y), ]
        y <- yc
      }
    } else {
      warning("Check your data or rerun with na.rm = TRUE")
    }
  }

  nx <- dim(x)[2]
  ny <- dim(y)[2]
  cross <- t(y) %*% x
  sumsx <- sqrt(1/diag(t(x) %*% x))
  sumsy <- sqrt(1/diag(t(y) %*% y))
  result <- matrix(rep(0, nx * ny), ncol = nx)
  result <- round(sumsy * (cross * rep(sumsx, each = ny)),
                  digits)
  return(t(result))

}


.onUnload <- function (libpath) {
  library.dynam.unload("EFAdiff", libpath)
}
