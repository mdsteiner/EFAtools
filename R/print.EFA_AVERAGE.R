#' Print EFA_AVERAGE object
#'
#' Print Method showing a summarized output of the [EFA_AVERAGE] function
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
#' \dontrun{
#' EFA_aver <- EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500)
#' EFA_aver
#' }
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
      .efa_style(ifelse(averaging == "median", "median",
                        paste0("mean (trim = ", settings$trim, ")")), "bold"),
      " across ", .efa_style(no_efas, "bold"), " EFAs, ",
      "varying the following settings: ",
      .settings_string(varied_settings), ".", sep = "")
  cat("\n")

  cat("\n")
  cat("The error rate is at ",
      .efa_style(paste0(round(mean(grid$errors, na.rm = TRUE) * 100), "%"), "bold"),
                   ". Of the solutions that did not result in an error, ",
      .efa_style(paste0(round(mean(grid$converged == 0, na.rm = TRUE) * 100), "%"),
                 "bold"),
      " converged, ",
      .efa_style(paste0(round(mean(grid$heywood, na.rm = TRUE) * 100), "%"), "bold"),
      " contained Heywood cases, and ",
      .efa_style(paste0(round(mean(grid$admissible, na.rm = TRUE) * 100), "%"), "bold"),
      " were admissible.", sep = "")
  cat("\n")
  cat("\n")

  # If no solutions were achieved across which averaging could be performed,
  # stop here with a message. Else, continue printing loadings etc.
  if(all((grid$converged != 0 | grid$errors | grid$heywood) %in% TRUE)){

    cli::cli_warn("No solutions were achieved across which averaging was possible. Best try again with a different number of factors.",
                  class = "efa_average_no_solutions")

  } else {

    fit <- x$fit_indices
    rownames(fit) <- fit$index

  # Indicator-to-factor correspondences
  cat("\n")
  cat(cli::rule(left = .efa_style("Indicator-to-Factor Correspondences", "bold"),
                line = 2))
  cat("\n")
  cat("\n")
  cat("For each cell, the proportion of solutions including the respective indicator-to-factor correspondence. A salience threshold of",
      .efa_style(settings$salience_threshold, "bold"),
      "was used to determine indicator-to-factor correspondences.")
  cat("\n")
  cat("\n")
  cat(format(structure(x$ind_fac_corres, class = "LOADINGS"), cutoff = 1e-4,
             digits = 2), sep = "\n")
  cat("\n")

  # Print the loadings
    cat("\n")
    cat(cli::rule(left = .efa_style("Loadings", "bold"), line = 2))
    cat("\n")
    .print_average(x, what = c("loadings"), stat = stat, averaging = averaging)
    cat("\n")

  ## Print Phi for oblique solutions
  if(!all(is.na(x$Phi))){

    cat("\n")
    cat(cli::rule(left = .efa_style("Factor Intercorrelations from Oblique Solutions", "bold"), line = 2))
    cat("\n")
    .print_average(x, what = c("Phi"), stat = stat, averaging = averaging)
    cat("\n")
  }

  # Variances accounted for
    cat("\n")
    cat(cli::rule(left = .efa_style("Variances Accounted for", "bold"),
                  line = 2))
    cat("\n")
    .print_average(x, what = c("vars_accounted"), stat = stat,
                   averaging = averaging)
    cat("\n")


  # Print fit indices
  if (fit["df", "average"] == 0) {
    cat("\n")
    cat(.efa_style("!", c("yellow", "bold")),
        .efa_style(" The model is just identified (df = 0). Goodness of fit indices may not be interpretable.", "yellow"))
    cat("\n")
  }

    cat("\n")
    cat(cli::rule(left = .efa_style("Model Fit", "bold"), line = 2))
    cat("\n")
    cat("\n")
    cat(paste0("       ", ifelse(averaging == "mean", "M", "Md"),
               " (SD) [Min; Max]"))
    cat("\n")

  if(all(method == "PAF") || is.na(N)){

    .print_gof(fit, ind = "caf", ind_name = "CAF:  ", print_zero = FALSE, digits = 2)
    cat("df: ",
        .numformat(fit["df", "average"], 0, print_zero = TRUE), "\n", sep = "")

  } else {

    .print_gof(fit, ind = c("chisq"), ind_name = "\U1D712\U00B2: ",
               print_zero = TRUE, digits = 2)
    cat("df: ",
        .numformat(fit["df", "average"], 0, print_zero = TRUE), "\n", sep = "")
    .print_gof(fit, ind = c("p_chi", "cfi", "rmsea", "aic", "bic", "caf"),
               ind_name = c(.efa_style("p: ", "italic"), "CFI: ", "RMSEA: ",
                            "AIC: ", "BIC: ", "CAF: "),
               print_zero = c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE),
               digits = c(3, 2, 2, 2, 2, 2))

  }

    # Plot loadings
    if(isTRUE(plot)){

    if(ncol(x$loadings$average) <= 10){

      print(plot(x))

    } else {

      message(cli::col_cyan(cli::symbol$info, " The factor solution contained more than 10 factors, no plot was generated. If you still want to create the plot, use 'plot(", substitute(x) ,")'.\n"))

    }

    }

  }

}

.print_average <- function(x, what, stat, averaging){

  # Render one statistic's table: average/min/max loadings are LOADINGS-classed (styled
  # loading table); everything else (sd/range loadings, Phi, variances) is a plain matrix
  # rendered through the shared renderer. Phi is symmetric, so only its lower triangle shows.
  render <- function(stat_key) {
    if (what == "loadings" && stat_key %in% c("average", "min", "max")) {
      print(x$loadings[[stat_key]])
    } else if (what == "Phi") {
      .print_efa_matrix(x$Phi[[stat_key]], role = "corr", lower_only = TRUE)
    } else if (what == "loadings") {
      .print_efa_matrix(x$loadings[[stat_key]], role = "corr")
    } else {
      .print_efa_matrix(x$vars_accounted[[stat_key]], role = "corr")
    }
  }

  if("average" %in% stat){

    cat("\n")
    cat(cli::rule(left = .efa_style(if (averaging == "mean") "Mean" else "Median", "bold")))
    cat("\n")
    cat("\n")

    render("average")

    }

  if("sd" %in% stat){

    cat("\n")
    cat(cli::rule(left = .efa_style("Standard Deviation", "bold")))
    cat("\n")
    cat("\n")

    render("sd")

  }

  if("range" %in% stat){

    cat("\n")
    cat(cli::rule(left = .efa_style("Range", "bold")))
    cat("\n")
    cat("\n")

    render("range")

  }

  if("min" %in% stat){

    cat("\n")
    cat(cli::rule(left = .efa_style("Minimum", "bold")))
    cat("\n")
    cat("\n")

    render("min")

  }

  if("max" %in% stat){

    cat("\n")
    cat(cli::rule(left = .efa_style("Maximum", "bold")))
    cat("\n")
    cat("\n")

    render("max")

  }
}

.print_gof <- function(fit, ind, ind_name, print_zero, digits){

    for(i in seq_along(ind)){

      if(ind[i] %in% c("p_chi", "cfi", "rmsea", "caf")){

  cat(ind_name[i],
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

        cat(ind_name[i],
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


