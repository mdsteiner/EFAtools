#' Print EFA_AVERAGE object
#'
#' Print Method showing a summarized output of the [EFA_AVERAGE] function
#'
#' @param x list. An object of class EFA_AVERAGE to be printed
#' @param stat character. A vector with the statistics to print. Possible inputs
#' are "average", "sd", "range", "min", and "max". Default is "average" and
#' "range".
#' @param plot logical. Whether a plot of the average and min- max loadings should
#' be created. Default is FALSE. If more than 10 factors are extracted, no plot is
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
                              plot = FALSE, ...) {

  checkmate::assert_subset(stat, c("average", "sd", "range", "min", "max"),
                           empty.ok = FALSE)

  # extract settings
  settings <- x$settings
  method <- settings$method
  rotation <- settings$rotation
  N <- settings$N
  grid <- x$implementations_grid
  averaging <- settings$averaging

  # settings that were varied: keep only the grid columns that correspond to EFA
  # settings. The per-model outcome columns appended to the grid (errors,
  # convergence, Heywood/admissibility flags, fit indices) are never part of the
  # settings list, so this excludes them without a column to maintain.
  varied_settings <- grid[, intersect(names(grid), names(settings)), drop = FALSE]
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
  print(structure(x$ind_fac_corres, class = "LOADINGS"), cutoff = 1e-4,
        digits = 2)

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

    .print_gof(fit, ind = c("caf", "rmsr", "srmr"),
               ind_name = c("CAF:  ", "RMSR: ", "SRMR: "),
               print_zero = c(FALSE, FALSE, FALSE), digits = c(2, 2, 2))
    cat("df: ",
        .numformat(fit["df", "average"], 0, print_zero = TRUE), "\n", sep = "")

  } else {

    .print_gof(fit, ind = c("chisq"), ind_name = "\U1D712\U00B2: ",
               print_zero = TRUE, digits = 2)
    cat("df: ",
        .numformat(fit["df", "average"], 0, print_zero = TRUE), "\n", sep = "")
    .print_gof(fit, ind = c("p_chi", "cfi", "tli", "rmsea", "aic", "bic",
                            "ecvi", "caf", "rmsr", "srmr"),
               ind_name = c(.efa_style("p: ", "italic"), "CFI: ", "TLI: ",
                            "RMSEA: ", "AIC: ", "BIC: ", "ECVI: ", "CAF: ",
                            "RMSR: ", "SRMR: "),
               print_zero = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
                              TRUE, FALSE, FALSE, FALSE),
               digits = c(3, 2, 2, 2, 2, 2, 2, 2, 2, 2))

  }

    # Plot loadings
    if(isTRUE(plot)){

    if(ncol(x$loadings$average) <= 10){

      print(plot(x))

    } else {

      obj_name <- deparse1(substitute(x))
      cli::cli_inform(c(
        "i" = "The factor solution contained more than 10 factors, no plot was generated.",
        "i" = "If you still want to create the plot, use {.code plot({obj_name})}."
      ))

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

  # Render "M (SD) [Min; Max]" for each index. .numformat drops the leading zero
  # of values in (-1, 1) when print_zero is FALSE (and keeps it otherwise) and
  # preserves the sign, so it already produces the stripped/unstripped form each
  # index needs. The integer-valued indices (print_zero = TRUE, e.g. AIC/BIC) are
  # left-padded to align the bracketed range; the decimal-only ones are not.
  fmt <- function(stat, i) {
    .numformat(fit[ind[i], stat], digits = digits[i],
               print_zero = print_zero[i], pad = print_zero[i])
  }

  for (i in seq_along(ind)) {
    cat(ind_name[i],
        fmt("average", i), " (", fmt("sd", i), ") [",
        fmt("min", i), "; ", fmt("max", i), "]\n", sep = "")
  }

}


