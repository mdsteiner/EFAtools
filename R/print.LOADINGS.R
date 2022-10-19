#' Print LOADINGS object
#'
#' @param x class LOADINGS matrix.
#' @param cutoff numeric. The number above which to print loadings in bold
#'  default is .3.
#' @param digits numeric. Passed to \code{\link[base:Round]{round}}. Number of digits
#'  to round the loadings to (default is 3).
#' @param max_name_length numeric. The maximum length of the variable names to
#'  display. Everything beyond this will be cut from the right.
#' @param ... additional arguments passed to print
#'
#' @method print LOADINGS
#' @export
#'
#' @examples
#' EFAtools_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                     type = "EFAtools", method = "PAF", rotation = "promax")
#' EFAtools_PAF
#'
print.LOADINGS <- function(x, cutoff = .3, digits = 3, max_name_length = 10,
                           ...) {

  # create factor names to display
  factor_names <- colnames(x)
  if (is.null(factor_names)) {
    factor_names <- paste0("F", seq_len(ncol(x)))
  }

  # for equal spacing, fill the factor names such that they match the columns
  fn_nchar <- sapply(factor_names, nchar)
  factor_names[which(fn_nchar > digits + 2)] <- substr(
    factor_names[which(fn_nchar > digits + 2)] , 1, digits + 2)
  factor_names <- stringr::str_pad(factor_names, digits + 2, side = "both")

  var_names <- rownames(x)
  if (is.null(var_names)) {
    var_names <- paste0("V", seq_len(nrow(x)))
  }

  max_char <- max(sapply(var_names, nchar))

  if (max_char > max_name_length) {
    vn_nchar <- sapply(var_names, nchar)
    ind <- which(vn_nchar > max_name_length)
    var_names[ind] <- substr(var_names[ind] , 1, max_name_length)
    max_char <- max_name_length
  }

  var_names <- stringr::str_pad(var_names, max_char, side = "right")

  n_col <- ncol(x)

  # create the string to paste using the crayon package
  temp <- apply(matrix(seq_len(nrow(x)), ncol = 1), 1,
                function(ind, x, cutoff, n_col, vn, digits){
    i <- x[ind,]

    tt <- crayon::blue(vn[ind])

    for (kk in seq_len(n_col)) {
      if (abs(i[kk]) < cutoff) {
        tt <- c(tt, crayon::silver(.numformat(round(i[kk], digits = digits),
                                              digits = digits)))
      } else {
        tt <- c(tt, crayon::bold(.numformat(round(i[kk], digits = digits),
                                            digits = digits)))
      }
    }
    stringr::str_c(tt, collapse = "\t")
  }, cutoff = cutoff, n_col = n_col, digits = digits, x = x, vn = var_names)

  factor_names <- stringr::str_c(factor_names,
                                 collapse = "\t")
  factor_names <- crayon::blue(stringr::str_c( stringr::str_pad(" ", max_char),
                                               "\t", factor_names))

  temp <- stringr::str_c(temp, collapse = "\n")

  temp <- stringr::str_c(factor_names, "\n", temp)

  temp <- stringr::str_c(temp, "\n")
  # print the results to the console
  cat(temp)
}
