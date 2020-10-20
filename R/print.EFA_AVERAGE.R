#' Print EFA_AVERAGE object
#'
#' Print Method showing a summarized output of the \link{EFA_AVERAGE} function
#'
#' @param x list. An object of class EFA_AVERAGE to be printed
#' @param stat. character. A vector with the statistics to print. Possible inputs
#' are "aggregate", "sd", "range", "min", and "max". Default is "aggregate" and
#' "range".
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
print.EFA_AVERAGE <- function(x, stat = c("aggregate", "range"),
                              plot = TRUE, ...) {

  checkmate::assert_subset(stat, c("aggregate", "sd", "range", "min", "max"),
                           empty.ok = FALSE)

  # extract settings
  settings <- x$settings
  method <- settings$method
  rotation <- settings$rotation
  N <- settings$N
  grid <- x$implementations_grid
  aggregation <- settings$aggregation

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
      crayon::bold(ifelse(aggregation == "median", "median",
                          paste0("mean (trim = ", settings$trim, ")", sep = ""))),
      " across ", crayon::bold(no_efas), " EFAs, ",
      "varying the following settings: ",
      .settings_string(varied_settings), ".", sep = "")
  cat("\n")

  cat("\n")
  cat("The error rate is at ",
      crayon::bold(round(mean(grid$errors, na.rm = TRUE) * 100), "%", sep = ""),
                   ". Of the solutions that did not result in an error, ",
      crayon::bold(round(mean(grid$converged == 0, na.rm = TRUE) * 100), "%",
                   sep = ""),
      " converged, ",
      crayon::bold(round(mean(grid$heywood, na.rm = TRUE) * 100), "%", sep = ""),
      " contained Heywood cases, and ",
      crayon::bold(round(mean(grid$admissible, na.rm = TRUE) * 100), "%", sep = ""),
      " were admissible.", sep = "")
  cat("\n")
  cat("\n")

  # Indicator-to-factor correspondences
  cat("\n")
  cat(cli::rule(left = crayon::bold("Indicator-to-Factor Correspondences"),
                col = "blue", line = 2))
  cat("\n")
  cat("\n")
  cat("For each cell, the proportion of solutions including the respective indicator-to-factor correspondence. A salience threshold of",
      crayon::bold(settings$salience_threshold),
      "was used to determine indicator-to-factor correspondences.")
  cat("\n")
  cat("\n")
  cat(print.LOADINGS(x$ind_fac_corres, cutoff = 1e-4, digits = 2))
  cat("\n")

  # Print the loadings

    cat("\n")
    cat(cli::rule(left = crayon::bold("Loadings"), col = "blue", line = 2))
    cat("\n")
    .print_average(x, what = c("loadings"), stat = stat, aggregation = aggregation)
    cat("\n")

  ## Print Phi for oblique solutions
  if(!all(is.na(x$Phi))){

    cat("\n")
    cat(cli::rule(left = crayon::bold("Factor Intercorrelations from Oblique Solutions"), col = "blue", line = 2))
    cat("\n")
    .print_average(x, what = c("Phi"), stat = stat, aggregation = aggregation)
    cat("\n")
  }

  # Variances accounted for
    cat("\n")
    cat(cli::rule(left = crayon::bold("Variances Accounted for"), col = "blue",
                  line = 2))
    cat("\n")
    .print_average(x, what = c("vars_accounted"), stat = stat,
                   aggregation = aggregation)
    cat("\n")


  # Print fit indices

  if (fit$df == 0) {
    cat("\n")
    cat(crayon::yellow$bold("!"), crayon::yellow(" The model is just identified (df = 0). Goodness of fit indices may not be interpretable."))
    cat("\n")
  }

  if(all(method == "PAF") || is.na(N)){

    cat("\n")
    cat(cli::rule(left = crayon::bold("Model Fit"), col = "blue", line = 2))
    cat("\n")
    cat("\n")

#     cat(crayon::blue("CAF:"),
#         .numformat(fit$CAF$aggregate), "\n", sep = "")
#     cat(crayon::blue("df: "),
#         .numformat(fit$df, 0, print_zero = TRUE), "\n", sep = "")
#
#
  } else {

    cat("\n")
    cat(cli::rule(left = crayon::bold("Model Fit"), col = "blue", line = 2))
    cat("\n")
    cat("\n")
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
  }

    # Plot loadings
    if(isTRUE(plot)){

    if(ncol(x$loadings$aggregate) <= 10){

      plot(x)

    } else {

      message(cli::col_cyan(cli::symbol$info, " The factor solution contained more than 10 factors, no plot was generated. If you still want to create the plot, use 'plot(", substitute(x) ,")'.\n"))

    }

    }

}

.print_average <- function(x, what, stat, aggregation){

  if("aggregate" %in% stat){

    if(aggregation == "mean"){

      cat("\n")
      cat(cli::rule(left = crayon::bold("Mean"), col = "blue"))
      cat("\n")
      cat("\n")

    } else {

      cat("\n")
      cat(cli::rule(left = crayon::bold("Median"), col = "blue"))
      cat("\n")
      cat("\n")

    }

    if(what == "loadings"){

      print(x$loadings$aggregate)

      } else if(what == "Phi"){

        cat(.get_compare_matrix(x$Phi$aggregate, r_red = Inf, n_char = 17,
                                var_names = paste0("F",
                                                   seq_len(ncol(x$Phi$aggregate)))))

      } else {

        cat(.get_compare_matrix(x$vars_accounted$aggregate, r_red = Inf,
                                n_char = 17))

      }

    }

  if("sd" %in% stat){

    cat("\n")
      cat(cli::rule(left = crayon::bold("Standard Deviation"), col = "blue"))
      cat("\n")
      cat("\n")

      if(what == "loadings"){

        cat(.get_compare_matrix(x$loadings$sd, r_red = Inf, n_char = 17))

      } else if(what == "Phi"){

        cat(.get_compare_matrix(x$Phi$sd, r_red = Inf, n_char = 17,
                                var_names = paste0("F",
                                                   seq_len(ncol(x$Phi$sd)))))

      } else {

        cat(.get_compare_matrix(x$vars_accounted$sd, r_red = Inf,
                                n_char = 17))

      }

  }

  if("range" %in% stat){

    cat("\n")
    cat(cli::rule(left = crayon::bold("Range"), col = "blue"))
    cat("\n")
    cat("\n")

    if(what == "loadings"){

      cat(.get_compare_matrix(x$loadings$range, r_red = Inf, n_char = 17))

    } else if(what == "Phi"){

      cat(.get_compare_matrix(x$Phi$range, r_red = Inf, n_char = 17,
                              var_names = paste0("F",
                                                 seq_len(ncol(x$Phi$range)))))

    } else {

      cat(.get_compare_matrix(x$vars_accounted$range, r_red = Inf,
                              n_char = 17))

    }

  }

  if("min" %in% stat){

    cat("\n")
    cat(cli::rule(left = crayon::bold("Minimum"), col = "blue"))
    cat("\n")
    cat("\n")

    if(what == "loadings"){

      print(x$loadings$min)

    } else if(what == "Phi"){

      cat(.get_compare_matrix(x$Phi$min, r_red = Inf, n_char = 17,
                              var_names = paste0("F",
                                                 seq_len(ncol(x$Phi$min)))))

    } else {

      cat(.get_compare_matrix(x$vars_accounted$min, r_red = Inf,
                              n_char = 17))

    }

  }

  if("max" %in% stat){

    cat("\n")
    cat(cli::rule(left = crayon::bold("Maximum"), col = "blue"))
    cat("\n")
    cat("\n")

    if(what == "loadings"){

      print(x$loadings$max)

    } else if(what == "Phi"){

      cat(.get_compare_matrix(x$Phi$max, r_red = Inf, n_char = 17,
                              var_names = paste0("F",
                                                 seq_len(ncol(x$Phi$max)))))

    } else {

      cat(.get_compare_matrix(x$vars_accounted$max, r_red = Inf,
                              n_char = 17))

    }

  }

  }

