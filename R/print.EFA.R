#' Print EFA object
#'
#' Print Method showing a summarized output of the \link{EFA} function
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
  N <- x$settings$N
  fit <- x$fit_indices

  cat("\n")
  cat("EFA performed with type = '", crayon::bold(type), "', method = '",
       crayon::bold(method), "', and rotation = '", crayon::bold(rotation),
       "'.", sep = "")
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
    cat(cli::rule(left = crayon::bold("Unrotated Loadings"), col = "blue"))
    cat("\n")
    cat("\n")
    print(x$unrot_loadings)

  } else {

    cat("\n")
    cat(cli::rule(left = crayon::bold("Rotated Loadings"), col = "blue"))
    cat("\n")
    cat("\n")
    print(x$rot_loadings)

    ## Print Phi for oblique solutions
    if(!is.null(x$Phi)){

      cat("\n")
      cat(cli::rule(left = crayon::bold("Factor Intercorrelations"), col = "blue"))
      cat("\n")
      cat("\n")
      cat(.get_compare_matrix(x$Phi, r_red = Inf, n_char = 17,
                              var_names = paste0("F", seq_len(ncol(x$Phi)))))

    }
  }

  cat("\n")
  cat(cli::rule(left = crayon::bold("Variances Accounted for"), col = "blue"))
  cat("\n")
  cat("\n")
  cat(.get_compare_matrix(x$vars_accounted, r_red = Inf, n_char = 17))

  if (fit$df == 0) {
    cat("\n")
    cat(crayon::yellow$bold("!"), crayon::yellow(" The model is just identified (df = 0). Goodness of fit indices may not be interpretable."))
    cat("\n")
  }

  if(method == "PAF" || is.na(N)){

    cat("\n")
    cat(cli::rule(left = crayon::bold("Model Fit"), col = "blue"))
    cat("\n")
    cat("\n")
    cat(crayon::blue("CAF:"),
        .numformat(fit$CAF), "\n", sep = "")
    cat(crayon::blue("df: "),
        .numformat(fit$df, 0, print_zero = TRUE), "\n", sep = "")


  } else {

    cat("\n")
    cat(cli::rule(left = crayon::bold("Model Fit"), col = "blue"))
    cat("\n")
    cat("\n")
    cat(crayon::blue("\U1D712\U00B2(", sep = ""), fit$df,
        crayon::blue(") = ", sep = ""),
        .numformat(fit$chi, 2, print_zero = TRUE), ", ",
        crayon::blue(crayon::italic("p")),
        ifelse(fit$p_chi < .001, " < .001",
               paste(crayon::blue(ifelse(fit$p_chi < 1, " =", " = ")),
                     .numformat(fit$p_chi, 3), sep = "")),
               "\n", sep = "")
    cat(crayon::blue(ifelse(fit$CFI < 1, "CFI =", "CFI = ")),
        .numformat(fit$CFI), "\n", sep = "")
    cat(crayon::blue(ifelse(fit$RMSEA < 1, "RMSEA [90% CI] =",
                            "RMSEA [90% CI] = ")),
        paste0(.numformat(fit$RMSEA), " [",
               ifelse(fit$RMSEA_LB < 1, substr(.numformat(fit$RMSEA_LB), 2, 4),
                      .numformat(fit$RMSEA_LB)),
               ifelse(fit$RMSEA_UB < 1, ";", "; "),
               .numformat(fit$RMSEA_UB), "]", sep = ""), "\n", sep = "")
    cat(crayon::blue("AIC = "),
        .numformat(fit$AIC, print_zero = TRUE), "\n", sep = "")
    cat(crayon::blue("BIC = "),
        .numformat(fit$BIC, print_zero = TRUE), "\n", sep = "")
    cat(crayon::blue("CAF ="),
        .numformat(fit$CAF), "\n", sep = "")

  }

}
