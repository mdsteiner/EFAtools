#' Print EFA_AVERAGE object
#'
#' Print Method showing a summarized output of the \link{EFA_AVERAGE} function
#'
#' @param x list. An object of class EFA_AVERAGE to be printed
#' @param stat character. A vector with the statistics to print. Possible inputs
#' are "average", "sd", "range", "min", and "max". Default is "average" and
#' "range".
#' @param plot logical. Whether a plot of the average and min- max loadings should
#' be created. Default is TRUE. If more than 10 factors are extracted, no plot is
#' created.
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
print.EFA_AVERAGE <- function(x, stat = c("average", "range"),
                              plot = TRUE, ...) {

  checkmate::assert_subset(stat, c("average", "sd", "range", "min", "max"),
                           empty.ok = FALSE)

  # extract settings
  settings <- x$settings
  method <- settings$method
  rotation <- settings$rotation
  N <- settings$N
  grid <- x$implementations_grid
  averaging <- settings$averaging

  # settings that were varied
  varied_settings <- grid
  varied_settings[, c("errors", "error_m", "converged", "heywood", "chisq",
                      "p_chi", "cfi", "caf", "rmsea", "aic", "bic")] <- NULL
  varied_settings <- apply(varied_settings, 2, function(x)unique(x[!is.na(x)]))
  varied_settings <- sapply(varied_settings, length)
  varied_settings <- names(varied_settings[varied_settings > 1])

  # extract other stuff
  no_efas <- nrow(grid)

  cat("\n")
  cat("Averaging performed with averaging method ",
      crayon::bold(ifelse(averaging == "median", "median",
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

  # If no solutions were achieved across which averaging could be performed,
  # stop here with a message. Else, continue printing loadings etc.
  if(all(grid$converged != 0 | grid$errors | grid$heywood)){

    warning(crayon::yellow$bold("!"), crayon::yellow(" No solutions were achieved across which averaging was possible. Best try again with a different number of factors.\n"))

  } else {

    fit <- x$fit_indices
    rownames(fit) <- fit$index

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
    .print_average(x, what = c("loadings"), stat = stat, averaging = averaging)
    cat("\n")

  ## Print Phi for oblique solutions
  if(!all(is.na(x$Phi))){

    cat("\n")
    cat(cli::rule(left = crayon::bold("Factor Intercorrelations from Oblique Solutions"), col = "blue", line = 2))
    cat("\n")
    .print_average(x, what = c("Phi"), stat = stat, averaging = averaging)
    cat("\n")
  }

  # Variances accounted for
    cat("\n")
    cat(cli::rule(left = crayon::bold("Variances Accounted for"), col = "blue",
                  line = 2))
    cat("\n")
    .print_average(x, what = c("vars_accounted"), stat = stat,
                   averaging = averaging)
    cat("\n")


  # Print fit indices
  if (fit["df", "average"] == 0) {
    cat("\n")
    cat(crayon::yellow$bold("!"), crayon::yellow(" The model is just identified (df = 0). Goodness of fit indices may not be interpretable."))
    cat("\n")
  }

    cat("\n")
    cat(cli::rule(left = crayon::bold("Model Fit"), col = "blue", line = 2))
    cat("\n")
    cat("\n")
    cat(crayon::blue("       ", ifelse(averaging == "mean", "M", "Md"),
                     " (SD) [Min; Max]", sep = ""))
    cat("\n")

  if(all(method == "PAF") || is.na(N)){

    .print_gof(fit, ind = "caf", ind_name = "CAF:  ", print_zero = FALSE, digits = 2)
    cat(crayon::blue("df: "),
        .numformat(fit["df", "average"], 0, print_zero = TRUE), "\n", sep = "")

  } else {

    .print_gof(fit, ind = c("chisq"), ind_name = "\U1D712\U00B2: ",
               print_zero = TRUE, digits = 2)
    cat(crayon::blue("df: "),
        .numformat(fit["df", "average"], 0, print_zero = TRUE), "\n", sep = "")
    .print_gof(fit, ind = c("p_chi", "cfi", "rmsea", "aic", "bic", "caf"),
               ind_name = c(crayon::italic("p: "), "CFI: ", "RMSEA: ",
                            "AIC: ", "BIC: ", "CAF: "),
               print_zero = c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE),
               digits = c(3, 2, 2, 2, 2, 2))

  }

    # Plot loadings
    if(isTRUE(plot)){

    if(ncol(x$loadings$average) <= 10){

      plot(x)

    } else {

      message(cli::col_cyan(cli::symbol$info, " The factor solution contained more than 10 factors, no plot was generated. If you still want to create the plot, use 'plot(", substitute(x) ,")'.\n"))

    }

    }

  }

}

.print_average <- function(x, what, stat, averaging){

  if("average" %in% stat){

    if(averaging == "mean"){

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

      print(x$loadings$average)

      } else if(what == "Phi"){

        cat(.get_compare_matrix(x$Phi$average, r_red = Inf, n_char = 17,
                                var_names = paste0("F",
                                                   seq_len(ncol(x$Phi$average)))))

      } else {

        cat(.get_compare_matrix(x$vars_accounted$average, r_red = Inf,
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

.print_gof <- function(fit, ind, ind_name, print_zero, digits){

    for(i in seq_along(ind)){

      if(ind[i] %in% c("p_chi", "cfi", "rmsea", "caf")){

  cat(crayon::blue(ind_name[i], sep = ""),
      ifelse(round(fit[ind[i], "average"], digits[i]) < 1,
             substr(.numformat(fit[ind[i], "average"], digits = digits[i],
                               print_zero = print_zero[i]),
                    2, digits + 2),
             .numformat(fit[ind[i], "average"], digits = digits[i],
                        print_zero = print_zero[i])), " (",
      ifelse(round(fit[ind[i], "sd"], digits[i]) < 1,
             substr(.numformat(fit[ind[i], "sd"], digits = digits[i],
                               print_zero = print_zero[i]),
                    2, digits + 2),
             .numformat(fit[ind[i], "sd"], digits = digits[i],
                        print_zero = print_zero[i])),
      ") [",
      ifelse(round(fit[ind[i], "min"], digits[i]) < 1,
             substr(.numformat(fit[ind[i], "min"], digits = digits[i],
                               print_zero = print_zero[i]),
                    2, digits + 2),
             .numformat(fit[ind[i], "min"], digits = digits[i],
                        print_zero = print_zero[i])), "; ",
      ifelse(round(fit[ind[i], "max"], digits[i]) < 1,
             substr(.numformat(fit[ind[i], "max"], digits = digits[i],
                               print_zero = print_zero[i]),
                    2, digits + 2),
             .numformat(fit[ind[i], "max"], digits = digits[i],
                        print_zero = print_zero[i])),
      "]\n", sep = "")

      } else {

        cat(crayon::blue(ind_name[i], sep = ""),
            .numformat(fit[ind[i], "average"], digits = digits[i],
                                     print_zero = print_zero[i]), " (",
            .numformat(fit[ind[i], "sd"], digits = digits[i],
                       print_zero = print_zero[i]),
            ") [",
            .numformat(fit[ind[i], "min"], digits = digits[i],
                       print_zero = print_zero[i]),
            "; ",
            .numformat(fit[ind[i], "max"], digits = digits[i],
                       print_zero = print_zero[i]),
            "]\n", sep = "")


      }

    }

  }


