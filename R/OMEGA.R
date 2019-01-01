# Input:
# var_names = a vector with subtest names in the same order as g_load
# fac_names = optional vector of group factor names in the order of the S-L solution
# factor_corres = vector that indicates which factor corresponds to which
#                 variable. The numbers must correspond to the columns of
#                 s_load. E.g. c(3, 3, 3, 1, 1, 2, 2) if the first three
#                 variables load on factor three, the next two on factor
#                 one and the last two on factor two.
# g_load = a vector of g loadings from an S-L solution
# s_load = a matrix of group factor loadings from an S-L solution
# u2 = a vector of uniquenesses (1 - h2) from an S-L solution
# Phi = factor intercorrelations from the oblique solution used for S-L
# pattern = pattern coefficients from the oblique solution used for S-L
# type = one of c("EFAdiff", "psych", "Watkins")

# Input:
# model = a psych Schmid-Leiman output

# Input:
# model: A fit object from lavaan bifactor cfa

# Input:
# model = an output of class "SL"
# fac_names = optional vector of group factor names in the order of the S-L solution
# factor_corres = vector that indicates which factor corresponds to which
#                 variable. The numbers must correspond to the columns of
#                 s_load. E.g. c(3, 3, 3, 1, 1, 2, 2) if the first three
#                 variables load on factor three, the next two on factor
#                 one and the last two on factor two.

#' MacDonald's omega
#'
#' This function finds omega total, omega hierarchical, and omega subscale
#' from a Schmid-Leiman solution or lavaan bifactor solution. The Schmid-Leiman
#' based omegas can either be found using the outputs from
#' \code{\link[psych::schmid]{psych::schmid}} or from \code{\link{SL}}, or, in a
#' more flexible way, by leaving model empty and specifying additional arguments.
#' This last option allows to reproduce results from \code{\link[psych::omega]{psych::omega}}
#' or Watkin's omega program (Watkins, 2013).
#'
#' @param model class \code{\link{SL}}, class \code{\link{schmid}}, or class
#' \code{\link{lavaan}} object. An output object from \code{\link{SL}}, from
#' \code{\link[psych::schmid]{psych::schmid}}, or from a \code{\link{lavaan}}
#' bifactor model. If of class \code{\link{lavaan}}, no other arguments need to
#' be specified. If of class \code{\link{SL}} or \code{\link{schmid}}, only
#' the argument \code{factor_corres} needs to be specified additionally.
#' @param var_names numeric.
#' @param fac_names character.
#' @param factor_corres numeric.
#' @param g_load numeric.
#' @param s_load numeric.
#' @param u2 numeric.
#' @param Phi matrix.
#' @param pattern matrix.
#' @param type character.
#'
#' @details If \code{model} is specified and of class \code{\link{lavaan}},
#' no other arguments need to be specified. If \code{model} is of class
#' \code{\link{SL}} or \code{\link{schmid}}, only the argument \code{factor_corres}
#' needs to be specified additionally.
#' If \code{model = NULL} and \code{type = "EFAdiff"}(default), the arguments
#' \code{var_names}, \code{factor_corres}, \code{g_load}, \code{s_load}, and \code{u2}
#' need to be specified.
#' If \code{model = NULL} and \code{type = "psych"}, the additional arguments
#' \code{Phi} and \code{pattern} must be specified. Additionally, the argument
#' \code{factor_corres} should be left NULL to replicate \code{\link[psych::omega]{psych::omega}}
#' results, where variable-to-factor correspondences are found by taking the highest
#' group factor loading for each variable as the relevant group factor loading.
#' If, however, the argument \code{factor_corres} is specified, the specified
#' variable-to-factor correspondences are taken, with a warning.
#' If \code{model = NULL} and \code{type = "Watkins"}, the \code{u2} argument
#' should be left NULL to replicate results from Watkins' omega program, where
#' uniquenesses are found based on the g loadings and relevant group factor
#' loadings only. If, however, the argument \code{u2} is specified, the specified
#' uniquenesses are taken, with a warning.
#'
#' @return

OMEGA <- function(model = NULL, var_names = NULL, fac_names = NULL,
                  factor_corres = NULL, g_load = NULL,
                  s_load = NULL, u2 = NULL, Phi = NULL, pattern = NULL,
                  type = "EFAdiff"){

  if(all(class(model) == "lavaan")){

    OMEGA_LAVAAN(model = model)

  } else if(all(class(model) == c("psych", "schmid"))) {

    if(is.null(factor_corres)){

      stop("Please specify the argument 'factor_corres'")

    } else {

      OMEGA_PSYCH(model = model, fac_names = fac_names,
                  factor_corres = factor_corres)
    }

  } else if(class(model) == "SL"){

    OMEGA_EFADIFF(model = model, fac_names = fac_names,
                  factor_corres = factor_corres)

  } else {

    if(is.null(var_names) | is.null(g_load) | is.null(s_load)){

      stop("Please specify all of the following arguments: 'var_names', 'g_load',
           's_load'")

    } else {

      OMEGA_FLEX(var_names = var_names, fac_names = fac_names,
                 factor_corres = factor_corres, g_load = g_load, s_load = s_load,
                 u2 = u2, Phi = Phi, pattern = pattern, type = type)

    }
  }
}
