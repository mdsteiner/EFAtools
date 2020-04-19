#' Print OMEGA object
#'
#' @param x output of class OMEGA (output from the \link{OMEGA} function)
#' @param digits numeric. Passed to \code{\link[base]{round}}. Number of digits
#' to round the loadings to (default is 3).
#' @param ... additional arguments passed to print
#'
#' @return
#' @method print OMEGA
#'
#' @export
print.OMEGA <- function(x, digits = 3, ...) {

  # In case of multiple groups
  if(is.list(x)){

    group_names <- names(x)

    if(length(x[[1]]) == 1){

      cat("Omega total for the single factor for each group:")
      cat("\n")

      for(i in 1:length(group_names)){

        cat("\n")
        cat("Group ", crayon::italic(group_names[i]), ": ",
            round(x[[i]], digits = digits), sep = "")

    }

    } else {

      cat("Omega total, omega hierarchical, and omega subscale for the general",
          "factor (top row) and the group factors for each group:")
      cat("\n")

      for(i in 1:length(group_names)){

        cat("\n")
        cat("Group ", crayon::italic(group_names[i]), ":", sep = "")
        cat("\n")
        print(round(x[[i]], digits = digits))

      }

    }

  } else { # In case of a single group

    if(length(x) == 1){

      cat("Omega total for the single factor:", round(x, digits = digits))

    } else {

      cat("Omega total, omega hierarchical, and omega subscale for the general",
      "factor (top row) and the group factors:")
      cat("\n")
      cat("\n")
      print(round(unclass(x), digits = digits))

    }

}

}

