#' MacDonald's omega
#'
#' This function finds omega total, omega hierarchical, and omega subscale
#' from a Schmid-Leiman (S-L) solution or lavaan bifactor solution. The S-L-based
#' omegas can either be found using the outputs from
#' \code{\link[psych:schmid]{psych::schmid}} or from \code{\link{SL}}, or, in a
#' more flexible way, by leaving model empty and specifying additional arguments.
#' By setting the \code{type} argument, results from \code{\link[psych:omega]{psych::omega}}
#' or Watkins' omega program (Watkins, 2013) can be reproduced.
#'
#' @param model class \code{\link{SL}}, class \code{\link{schmid}}, or class
#' \code{\link{lavaan}} object. An output object from \code{\link{SL}}, from
#' \code{\link[psych:schmid]{psych::schmid}}, or from a \code{\link{lavaan}}
#' bifactor model. If of class \code{\link{lavaan}}, no other arguments need to
#' be specified. If of class \code{\link{SL}} or \code{\link{schmid}}, only
#' the argument \code{factor_corres} needs to be specified additionally.
#' @param var_names numeric. A vector with subtest names in the same order as
#' \code{g_load}.
#' @param fac_names character. An optional vector of group factor names in the order
#' of the S-L solution.
#' @param factor_corres numeric. A vector that indicates which variable corresponds
#' to which group factor. Must be in the same order as the S-L solution. E.g.
#' c(3, 3, 3, 1, 1, 2, 2) if the first three variables load on the third group
#' factor of the S-L solution, the next two on the first group factor and the
#' last two on the second group factor. If a variable should not be assigned to
#' any group factor, insert a zero at its position (e.g. c(3, 3, 0, 1, 1, 2, 2),
#' the third variable has no corresponding group factor).
#' @param g_load numeric. A vector of g loadings from an S-L solution.
#' @param s_load matrix. A matrix of group factor loadings from an S-L solution.
#' @param u2 numeric. A vector of uniquenesses from an S-L solution.
#' @param cormat matrix. A correlation matrix to be used when type = "psych" or
#' type = "EFAdiff". If left NULL, the correlation matrix is found based on the
#' pattern matrix and Phi using \code{\link[psych:factor.model]{psych::factor.model}}.
#' @param Phi matrix. Factor intercorrelations from an oblique factor solution.
#' @param pattern matrix. Pattern coefficients from an oblique factor solution.
#' @param variance character. If "correlation", then total variances for the whole
#' scale as well as the group factors are calculated
#' based on the correlation matrix. If "sums_load", then total variances are
#' calculated using the squared sums of g loadings and group factor loadings and
#' the sum of uniquenesses (see details).
#' @param type character. Either "EFAdiff" (default), "psych", or "Watkins".
#'
#' @details If \code{model} is specified and of class \code{\link{lavaan}},
#' no other arguments need to be specified.
#' If \code{model} is of class
#' \code{\link{SL}} or \code{\link{schmid}}, only the argument \code{factor_corres}
#' needs to be specified additionally. There is, however, the option to reproduce
#' Watkins' or \code{\link[psych:omega]{psych::omega}} results by setting the
#' \code{type} argument to "Watkins" or "psych".
#' If \code{model = NULL} and \code{type = "EFAdiff"}(default), the arguments
#' \code{var_names}, \code{factor_corres}, \code{g_load}, \code{s_load}, and \code{u2}
#' need to be specified.
#' If \code{type = "psych"} or \code{type = "EFAdiff"}, either \code{cormat}
#' (recommended) or \code{Phi} and \code{pattern} must be specified.
#' Additionally, the argument \code{factor_corres} should be left NULL to
#' replicate \code{\link[psych:omega]{psych::omega}}
#' results, where variable-to-factor correspondences are found by taking the highest
#' group factor loading for each variable as the relevant group factor loading.
#' If \code{type = "Watkins"}, the \code{u2} argument
#' should be left NULL to replicate results from Watkins' omega program, where
#' uniquenesses are found based on the g loadings and relevant group factor
#' loadings only. If, however, the argument \code{u2} is specified, the specified
#' uniquenesses are taken, with a warning.
#'
#' Explain omegas (incl. formula) for all type settings. 1 - omega tot as
#' variance in the composite that is not explained by g or the group factor.
#'
#' @return
#' A matrix with Omegas for g / the whole scale and for the subscales.
#' \item{omega_tot}{Omega total.}
#' \item{omega_h}{Omega hierarchical.}
#' \item{omega_sub}{Omega subscale.}
#' \item{1 - omega_tot}{1 - Omega total.}
#' @export
OMEGA <- function(model = NULL, var_names = NULL, fac_names = NULL,
                  factor_corres = NULL, g_load = NULL,
                  s_load = NULL, u2 = NULL, Phi = NULL, pattern = NULL,
                  cormat = NULL, variance = NULL, type = "EFAdiff"){

  if(!is.null(model) & (!is.null(var_names) || !is.null(g_load) || !is.null(s_load)
                        || !is.null(u2))){

    warning("You entered a model and specified at least one of the arguments
            var_names, g_load, s_load, or u2. These arguments are ignored.
            To use specific values for these, leave model = NULL and specify
            all arguments separately.")

  }

  if(all(class(model) == "lavaan")){

    OMEGA_LAVAAN(model = model)

  } else if(all(class(model) == c("psych", "schmid")) || all(class(model) == "SL")) {

    if(is.null(factor_corres & all(class(model) == c("psych", "schmid")))){

      stop("Please specify the argument 'factor_corres'")

    } else {

      OMEGA_FLEX(model = model, var_names = var_names, fac_names = fac_names,
                 factor_corres = factor_corres, g_load = g_load, s_load = s_load,
                 u2 = u2, Phi = Phi, pattern = pattern, cormat = cormat,
                 variance = variance, type = type)
    }

    } else {

      if(!is.null(model)){

        stop("Invalid input for model. Either enter a lavaan, psych::schmid or
             SL object or specify the arguments var_names, g_load, and s_load.")

      }

      if(is.null(var_names) | is.null(g_load) | is.null(s_load)){

      stop("Please specify all of the following arguments: 'var_names', 'g_load',
           's_load'")

        } else {

          OMEGA_FLEX(model = model, var_names = var_names, fac_names = fac_names,
                 factor_corres = factor_corres, g_load = g_load, s_load = s_load,
                 u2 = u2, Phi = Phi, pattern = pattern, cormat = cormat,
                 variance = variance, type = type)

        }

      }
}
