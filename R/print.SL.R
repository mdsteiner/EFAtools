#' Print SL object
#'
#' Print Method showing a summarized output of the \link{SL} function.
#'
#' @param x list. An object of class SL to be printed
#' @param ...  Further arguments for print.
#'
#' @export
#' @method print SL
#'
#' @examples
#' EFA_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                type = "EFAtools", method = "PAF", rotation = "promax")
#' SL(EFA_mod, type = "EFAtools", method = "PAF")
#'
print.SL <- function(x, ...) {

  if(!all(is.na(x$settings))){
  # extract the settings for EFA not depending on the method
  method <- x$settings$method
  type <- x$settings$type

  # Settings intro message
  cat("\n")
  cat("EFA for second-order loadings performed with type = '",
      crayon::bold(type), "' and method = '", crayon::bold(method),
      "'", sep = "")
  cat("\n")

  if (!is.null(x$settings$max_iter) && x$iter > x$settings$max_iter) {
    cat("\n")
    cat(crayon::red$bold(cli::symbol$cross,
                         "Maximum number of iterations reached",
                         "without convergence"))
    cat("\n")
  }

  }

  # print the loadings and the variances
  cat("\n")
  cat(cli::rule(left = crayon::bold("Schmid-Leiman Solution"), col = "blue"))
  cat("\n")
  cat("\n")
  print(x$sl)
  cat("\n")
  cat("\n")
  cat(cli::rule(left = crayon::bold("Variances Accounted for"), col = "blue"))
  cat("\n")
  cat("\n")
  cat(.get_compare_matrix(x$vars_accounted, r_red = Inf, n_char = 17))

}
