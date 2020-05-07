#' Print EFA object
#'
#' Print Method showing a summarized output of the EFA procedure.
#'
#' @param x list. An object of class EFA to be printed
#' @param ...  Further arguments for print.
#'
#' @export
#'
#' @method print EFA
#'
#' @examples
#' EFAtools_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                     type = "EFAtools", method = "PAF", rotation = "promax")
#' EFAtools_PAF
#'
print.EFA <- function(x, ...) {

  # extract the settings not depending on the method or rotation
  method <- x$settings$method
  rotation <- x$settings$rotation
  type <- x$settings$type

  # Settings intro message
  cat(crayon::blue("EFA performed with type = '", crayon::bold(type),
               "', method = '", crayon::bold(method),
               "', and rotation = '", crayon::bold(rotation),
               "'.", sep = ""))
  cat("\n")

  if (!is.null(x$settings$max_iter) && x$iter > x$settings$max_iter) {
    cat("\n")
    cat(crayon::red$bold(cli::symbol$cross,
                         "Maximum number of iterations reached",
                         "without convergence"))
    cat("\n")
  }

  # If no rotation was used, print the unrotated loadings, otherwise print the
  # rotated loadings
  if(rotation == "none"){

    cat("\n")
    cat(crayon::blue$bold("Unrotated Loadings:"))
    cat("\n")
    print(x$unrot_loadings)

  } else {

    cat("\n")
    cat(crayon::blue$bold("Rotated Loadings:"))
    cat("\n")
    print(x$rot_loadings)

    ## Print Phi for oblique solutions
    if(!is.null(x$Phi)){

      cat("\n")
      cat(crayon::blue$bold("Factor Intercorrelations:"))
      cat("\n")
      cat(.get_compare_matrix(x$Phi, r_red = Inf, n_char = 17,
                              var_names = paste0("F", 1:ncol(x$Phi))))

    }
  }

  cat("\n")
  cat(crayon::blue$bold("Variances Accounted for:"))
  cat("\n")
  cat(.get_compare_matrix(x$vars_accounted, r_red = Inf, n_char = 17))

  if(method == "PAF"){

    fitind <- matrix(c(x$fit_indices$CAF, x$fit_indices$df), nrow = 1)
    rownames(fitind) <- ""
    colnames(fitind) <- c("CAF", "df")

  cat("\n")
  cat(crayon::blue$bold("Model Fit:"))
  cat("\n")
  cat(.get_compare_matrix(fitind, r_red = Inf, n_char = 17, gof = TRUE))

  } else {

    fitind <- matrix(c(x$fit_indices$chi, x$fit_indices$df, x$fit_indices$CFI,
                       x$fit_indices$RMSEA, x$fit_indices$RMSEA_LB,
                       x$fit_indices$RMSEA_UB, x$fit_indices$AIC,
                       x$fit_indices$BIC, x$fit_indices$CAF), nrow = 1)
    colnames(fitind) <- c("\U1D712\U00B2", "df", "CFI", "RMSEA", "lower",
                          "upper", "AIC", "BIC", "CAF")

    cat("\n")
    cat(crayon::blue$bold("Model Fit:"))
    cat("\n")
    cat(.get_compare_matrix(fitind, r_red = Inf, n_char = 30, gof = TRUE))

  }

}



