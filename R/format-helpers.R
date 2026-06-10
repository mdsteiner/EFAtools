# Number formatting and styled-table builders shared by the print methods: padded decimal
# formatting, the LOADINGS/COMPARE cell renderers, and small settings/CI string helpers.

#' Format numbers for print method
#'
#' Helper function used in the print method for class LOADINGS and SLLOADINGS.
#' Strips the 0 in front of the decimal point of a number if number < 1, only
#' keeps the first `digits` number of digits, and adds an empty space in
#' front of the number if the number is positive. This way all returned strings
#' (except for those > 1, which are exceptions in LOADINGS) have the same number
#' of characters.
#'
#' @param x numeric. Number to be formatted.
#' @param digits numeric. Number of digits after the comma to keep.
#' @param print_zero logical. Whether, if a number is between \[-1, 1\], the
#'  zero should be omitted or printed (default is FALSE, i.e. omit zeros).
#' @param pad logical. Whether, if a number starts with a 0 and the 0 is not printed
#'  a white-space should be added.
#'
#' @return A formated number
.numformat <- function(x, digits = 2, print_zero = FALSE,
                       pad = TRUE) {

  if (isFALSE(print_zero)) {

    ncode <- paste0("%.", digits, "f")
    x <- sub("^(-?)0.", "\\1.", sprintf(ncode, x))
    if (isTRUE(pad)) {
      x <- stringr::str_pad(x, digits + 2, "left")
    }


  } else {

    ncode <- paste0("%.", digits, "f")
    x <- sprintf(ncode, x)

    if (isTRUE(pad)) {
      x <- stringr::str_pad(x, digits + 3, "left")
    }

  }
  x

}

.get_compare_matrix <- function(x, digits = 3, r_red = .001, n_char = 10,
                                var_names = NULL, factor_names = NULL,
                                gof = FALSE) {

  # create factor names to display
  if (is.null(factor_names)) {
    if(is.null(colnames(x))){
      factor_names <- paste0("F", seq_len(ncol(x)))
    } else {
      factor_names <- colnames(x)
    }
  }

  # for equal spacing, fill the factor names such that they match the columns
  fn_nchar <- sapply(factor_names, nchar)
  factor_names[which(fn_nchar > digits + 2)] <- substr(
    factor_names[which(fn_nchar > digits + 2)] , 1, digits + 2)
  factor_names <- stringr::str_pad(factor_names, digits + 2, side = "both")

  if(gof == FALSE){

  if(is.null(var_names)) {
    if(is.null(rownames(x))){
      var_names <- paste0("V", seq_len(nrow(x)))
    } else {
    var_names <- rownames(x)
    }
  }

  max_char <- max(sapply(var_names, nchar))

  if (max_char > n_char) {
    vn_nchar <- sapply(var_names, nchar)
    var_names[which(vn_nchar > n_char)] <- substr(var_names[which(vn_nchar > n_char)],
                                              1, n_char)
    max_char <- n_char
  }

  var_names <- stringr::str_pad(var_names, max_char, side = "right")

  }

  n_col <- ncol(x)

  # create the string to paste using the crayon package
  temp <- apply(matrix(seq_len(nrow(x)), ncol = 1), 1,
                function(ind, x, cutoff, n_col, vn, digits){
                  i <- x[ind,]

                  tt <- crayon::blue(vn[ind])

                  for (kk in seq_len(n_col)) {
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

  factor_names <- stringr::str_c(factor_names, collapse = "\t")

  if(gof == TRUE){

    factor_names <- crayon::blue(stringr::str_c(factor_names))

  } else {

    factor_names <- crayon::blue(stringr::str_c( stringr::str_pad(" ", max_char),
                                                   "\t", factor_names))
  }


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

  if ((is.null(dim(x)) && !(inherits(x, c("numeric", "integer")))) ||
      (!is.null(dim(x)) && !(inherits(x, c("matrix", "loadings", "LOADINGS",
                                          "SLLOADINGS"))))) {
    cli::cli_abort("{.arg x} must be a numeric vector or matrix, not {.cls {class(x)}}.",
                   class = "efa_not_numeric")
  }

  if (!is.null(dim(x))) {

    max(apply(x, 1:2, function(ll) {
      if (abs(ll - round(ll)) > .Machine$double.eps^0.5) {
        nchar(strsplit(sub('0+$', '', as.character(ll)), ".",
                       fixed = TRUE)[[1]][[2]])
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

# Create a string (used to print settings in print functions, to use in cat()
# with sep = "")
.settings_string <- function(x){

  n <- length(x)

if(n == 1){

  c(crayon::bold(x))

} else if (n == 2){

  c(crayon::bold(x[1]), " and ", crayon::bold(x[2]))

} else if (n > 2){

  c(paste(crayon::bold(x[seq_len(n-1)]), collapse = ", "), ", and ",
    crayon::bold(x[n]))

}

}


.paste_gof_ci <- function(indices, index) {
  paste0(" [", .numformat(indices$lower[index],
                         pad = FALSE), ", ",
         .numformat(indices$upper[index],
                           pad = FALSE),
         "]")
}
