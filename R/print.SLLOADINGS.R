#' Print SLLOADINGS object
#'
#' @param x class SLLOADINGS matrix.
#' @param cutoff numeric. The number above which to print loadings in bold
#'  (default is .2).
#' @param digits numeric. Passed to \code{\link[base:Round]{round}}. Number of digits
#'  to round the loadings to (default is 3).
#' @param ... additional arguments passed to print
#'
#' @method print SLLOADINGS
#' @export
#'
#' @examples
#' EFA_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                type = "EFAtools", method = "PAF", rotation = "promax")
#' SL(EFA_mod, type = "EFAtools", method = "PAF")
#'
print.SLLOADINGS <- function(x, cutoff = .2, digits = 3, ...) {

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

  if (max_char > 10) {
    vn_nchar <- sapply(var_names, nchar)
    var_names[which(vn_nchar > 10)] <- substr(var_names[which(vn_nchar > 10)] ,
                                                      1, 10)
    max_char <- 10
  }

  var_names <- stringr::str_pad(var_names, max_char, side = "right")

  n_col <- ncol(x)

  # create the string to paste using the crayon package
  temp <- apply(matrix(seq_len(nrow(x)), ncol = 1), 1,
                function(ind, x, cutoff, n_col, vn, digits){
    i <- x[ind,]

    tt <- crayon::blue(vn[ind])
    for (kk in seq_len(n_col)) {
      if (kk <= n_col - 2) {
        if (abs(i[kk]) < cutoff) {
          tt <- c(tt, crayon::silver(.numformat(round(i[kk], digits = digits),
                                                digits = digits)))
        } else if (abs(i[kk]) <= 1) {
          tt <- c(tt, crayon::bold(.numformat(round(i[kk], digits = digits),
                                              digits = digits)))
        } else {
          tt <- c(tt, crayon::red$bold(.numformat(round(i[kk], digits = digits),
                                                  digits = digits)))
        }
      } else {
        if (i[kk] <= 1 & i[kk] >= 0) {
          tt <- c(tt, .numformat(round(i[kk], digits = digits), digits = digits))
        } else {
          tt <- c(tt, crayon::red$bold(.numformat(round(i[kk], digits = digits),
                                                  digits = digits)))
        }
      }
    }
    stringr::str_c(tt, collapse = "\t")
  }, cutoff = cutoff, n_col = n_col, digits = digits, x = x, vn = var_names)

  # add header
  factor_names <- stringr::str_c(factor_names,
                                 collapse = "\t")
  factor_names <- crayon::blue(stringr::str_c( stringr::str_pad(" ", max_char),
                                               "\t", factor_names))

  temp <- stringr::str_c(temp, collapse = "\n")
  temp <- stringr::str_c(factor_names, "\n", temp)

  # warn from Heywood cases
  if (sum(x[, -n_col] > 1) == 1) {
    temp <- paste(temp,
                  crayon::red$bold("\nWarning: Results contain a Heywood case!"),
                  collapse = "\n")
  } else if (sum(x[, -(n_col-1):-n_col] > 1) > 1) {
    temp <- paste(temp, crayon::red$bold("\nWarning: Results contain",
                                         sum(x[, -(n_col-1):-n_col] > 1),
                                         "Heywood cases!"),
                  collapse = "\n")
  }

  temp <- stringr::str_c(temp, "\n")
  # print the results to the console
  cat(temp)
}
