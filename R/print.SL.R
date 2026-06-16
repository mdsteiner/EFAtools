#' Print SL object
#'
#' Print Method showing a summarized output of the [SL] function.
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
      .efa_style(type, "bold"), "' and method = '", .efa_style(method, "bold"),
      "'", sep = "")
  cat("\n")

  if (!is.null(x$settings$max_iter) && x$iter >= x$settings$max_iter) {
    cat("\n")
    cat(.efa_style(paste(cli::symbol$cross,
                         "Maximum number of iterations reached",
                         "without convergence"), c("red", "bold")))
    cat("\n")
  }

  }

  # print the loadings and the variances
  cat("\n")
  cat(cli::rule(left = .efa_style("Schmid-Leiman Solution", "bold")))
  cat("\n")
  cat("\n")
  print(x$sl)
  cat("\n")
  cat("\n")
  cat(cli::rule(left = .efa_style("Variances Accounted for", "bold")))
  cat("\n")
  cat("\n")
  .print_efa_corr_matrix(x$vars_accounted, digits = 3)

}
