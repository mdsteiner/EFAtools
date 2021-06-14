#' Print OMEGA object
#'
#' @param x output of class OMEGA (output from the \link{OMEGA} function)
#' @param digits numeric. Passed to \code{\link[base:Round]{round}}. Number of digits
#' to round to (default is 3).
#' @param ... additional arguments passed to print
#'
#' @method print OMEGA
#'
#' @export
#'
#' @examples
#' efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                type = "EFAtools", method = "PAF", rotation = "promax")
#' sl_mod <- SL(efa_mod, type = "EFAtools", method = "PAF")
#'
#' OMEGA(sl_mod, type = "EFAtools",
#' factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2)
#'
print.OMEGA <- function(x, digits = 3, ...) {

  # In case of multiple groups
  if(is.list(x)){

    group_names <- names(x)

    # If there is only a single factor and only omega was computed
    if(length(x[[1]]) == 1){

      cat(crayon::blue$bold("Omega total for the single factor for each group:"))
      cat("\n")

      for(i in seq_along(group_names)){

        cat("\n")
        cat(crayon::blue("Group ", crayon::italic(group_names[i]), ": ", sep = ""),
            round(x[[i]], digits = digits))

    }

      # If there is only a single factor and omega and H index were computed
    } else if(length(x[[1]]) == 2){

      cat(crayon::blue$bold("Omega total and H index for the single factor for each group:"))

      for(i in seq_along(group_names)){

        cat("\n")
        cat("\n")
        cat(crayon::blue("Group ", crayon::italic(group_names[i]), ": ", sep = ""))
        cat("\n")
        cat("Omega: ", round(x[[i]][1], digits = digits), sep = "")
        cat("\n")
        cat("H index: ", round(x[[i]][2], digits = digits), sep = "")

      }

    } else {

      if(ncol(x[[1]] == 3)){

      cat(crayon::blue$bold("Omega total, omega hierarchical, and omega subscale",
      "for the general factor (top row) and the group factors for each group:"))
      cat("\n")

      } else {

        cat(crayon::blue$bold("Omega total, omega hierarchical, omega subscale,",
        "H index, explained common variance (ECV), and percent of uncontaminated",
        "correlations (PUC) for the general factor (top row) and omegas and",
        "H index for the group factors for each group:"))
        cat("\n")

      }

      for(i in seq_along(group_names)){

        cat("\n")
        cat(crayon::blue("Group ", crayon::italic(group_names[i]), ":", sep = ""))
        cat("\n")
        print(round(unclass(x[[i]]), digits = digits), na.print = "")

      }

    }

  } else { # In case of a single group

    if(length(x) == 1){

      cat(crayon::blue$bold("Omega total for the single factor:"),
                            round(x, digits = digits))

    } else if(length(x) == 2){

      cat(crayon::blue$bold("Omega total and H index for the single factor:"))
      cat("\n")
      cat("\n")
      cat("Omega: ", round(x[1], digits = digits), sep = "")
      cat("\n")
      cat("H index: ", round(x[2], digits = digits), sep = "")

      } else {

        if(ncol(x) == 3){

      cat(crayon::blue$bold("Omega total, omega hierarchical, and omega subscale",
      "for the general factor (top row) and the group factors:"))
      cat("\n")
      cat("\n")

        } else {

          cat(crayon::blue$bold("Omega total, omega hierarchical, omega subscale,",
                                "H index, explained common variance (ECV), and",
                                "percent of uncontaminated correlations (PUC)",
                                "for the general factor (top row) and omegas",
                                "and H index for the group factors:"))
          cat("\n")
          cat("\n")

        }

        print(round(unclass(x), digits = digits), na.print = "")

    }

}

}

