#' Print EFA_AVERAGE object
#'
#' Print Method showing a summarized output of the \link{EFA_AVERAGE} function
#'
#' @param x list. An object of class EFA_AVERAGE to be printed
#' @param ...  Further arguments for print.
#'
#' @export
#'
#' @method print EFA_AVERAGE
#'
#' @examples
#'
#' EFA_aver <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500)
#' EFA_aver
#'
print.EFA_AVERAGE <- function(x, ...) {

  # extract settings
  settings <- x$settings
  method <- settings$method
  rotation <- settings$rotation
  N <- settings$N
  grid <- x$implementations_grid

  # settings that were varied
  varied_settings <- grid
  varied_settings[, c("errors", "error_m", "converged", "heywood", "chisq",
                      "p_chi", "cfi", "caf", "rmsea", "aic", "bic")] <- NULL
  varied_settings <- apply(varied_settings, 2, function(x)unique(x[!is.na(x)]))
  varied_settings <- sapply(varied_settings, length)
  varied_settings <- names(varied_settings[varied_settings > 1])

  # extract other stuff
  no_efas <- nrow(grid)
  fit <- x$fit_indices

  cat("\n")
  cat("Averaging performed with aggregation method ",
      crayon::bold(ifelse(settings$aggregation == "median", "median",
                          paste0("mean (trim = ", settings$trim, ")", sep = ""))),
      " across ", crayon::bold(no_efas), " EFAs, ",
      "varying the following settings: ",
      .settings_string(varied_settings), ".", sep = "")
  cat("\n")

  cat("\n")
  cat("Convergence rate is at ",
      crayon::bold(round((1 - sum(grid$converged) / no_efas) * 100), "%", sep = ""),
      " and Heywood cases occurred in ",
      crayon::bold(round((sum(grid$heywood) / no_efas) * 100), "%", sep = ""),
      " of the solutions.", sep = "")
  cat("\n")

  # If no rotation was used, print the unrotated loadings, otherwise print the
  # rotated loadings
  if(rotation == "none"){

    cat("\n")
    cat(cli::rule(left = crayon::bold("Unrotated Loadings"), col = "blue"))
    cat("\n")
    cat("\n")
    print(x$loadings$aggregate)

  } else {

    cat("\n")
    cat(cli::rule(left = crayon::bold("Rotated Loadings"), col = "blue"))
    cat("\n")
    cat("\n")
    print(x$loadings$aggregate)

    ## Print Phi for oblique solutions
    if(!is.null(x$Phi)){

      cat("\n")
      cat(cli::rule(left = crayon::bold("Factor Intercorrelations"), col = "blue"))
      cat("\n")
      cat("\n")
      cat(.get_compare_matrix(x$Phi$aggregate, r_red = Inf, n_char = 17,
                              var_names = paste0("F",
                                                 seq_len(ncol(x$Phi$aggregate)))))

    }
  }

  # Indicator-to-factor correspondences
  cat("\n")
  cat(cli::rule(left = crayon::bold("Indicator-to-Factor Correspondences"),
                col = "blue"))
  cat("\n")
  cat("\n")
  cat("A cutoff of", crayon::bold(settings$salience_threshold),
      "was used to determine indicator-to-factor correspondences.")
  cat("\n")
  cat("\n")
  print(x$ind_fac_corres, cutoff = 0, digits = 1)

  # Print fit indices

#   if (fit$df == 0) {
#     cat("\n")
#     cat(crayon::yellow$bold("!"), crayon::yellow(" The model is just identified (df = 0). Goodness of fit indices may not be interpretable."))
#     cat("\n")
#   }
#
#   if(all(method == "PAF") || is.na(N)){
#
#     cat("\n")
#     cat(cli::rule(left = crayon::bold("Model Fit"), col = "blue"))
#     cat("\n")
#     cat("\n")
#     cat(crayon::blue("CAF:"),
#         .numformat(fit$CAF$aggregate), "\n", sep = "")
#     cat(crayon::blue("df: "),
#         .numformat(fit$df, 0, print_zero = TRUE), "\n", sep = "")
#
#
#   } else {
#
#     cat("\n")
#     cat(cli::rule(left = crayon::bold("Model Fit"), col = "blue"))
#     cat("\n")
#     cat("\n")
#     cat(crayon::blue("\U1D712\U00B2(", sep = ""), fit$df,
#         crayon::blue(") = ", sep = ""),
#         .numformat(fit$chi$aggregate, 2, print_zero = TRUE), ", ",
#         crayon::blue(crayon::italic("p")),
#         ifelse(fit$p_chi$aggregate < .001, " < .001",
#                paste(crayon::blue(ifelse(fit$p_chi$aggregate < 1, " =", " = ")),
#                      .numformat(fit$p_chi$aggregate, 3), sep = "")),
#                "\n", sep = "")
#     cat(crayon::blue(ifelse(fit$CFI$aggregate < 1, "CFI =", "CFI = ")),
#         .numformat(fit$CFI$aggregate), "\n", sep = "")
#     cat(crayon::blue(ifelse(fit$RMSEA$aggregate < 1, "RMSEA =", "RMSEA  = ")),
#         .numformat(fit$RMSEA$aggregate), "\n", sep = "")
#     cat(crayon::blue("AIC = "),
#         .numformat(fit$AIC$aggregate, print_zero = TRUE), "\n", sep = "")
#     cat(crayon::blue("BIC = "),
#         .numformat(fit$BIC$aggregate, print_zero = TRUE), "\n", sep = "")
#     cat(crayon::blue("CAF ="),
#         .numformat(fit$CAF$aggregate), "\n", sep = "")
#
#   }

 }
