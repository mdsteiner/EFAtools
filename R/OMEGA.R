#' MacDonald's omega
#'
#'ALSO FOR LAVAAN SINGLE FACTOR SOLUTION
#'
#' This function finds omega total, omega hierarchical, and omega subscale
#' from a Schmid-Leiman (SL) solution or lavaan bifactor solution. The SL-based
#' omegas can either be found using the outputs from
#' \code{\link[psych:schmid]{psych::schmid}} or from \code{\link{SL}}, or, in a
#' more flexible way, by leaving \code{model = NULL} and specifying additional arguments.
#' By setting the \code{type} argument, results from
#' \code{\link[psych:omega]{psych::omega}} can be reproduced.
#'
#' @param model class \code{\link{SL}}, class \code{\link{schmid}}, or class
#' \code{\link{lavaan}} object. That is, an output object from \code{\link{SL}}, from
#' \code{\link[psych:schmid]{psych::schmid}}, or from a \code{\link{lavaan}}
#' bifactor model. If of class \code{\link{lavaan}}, only \code{g_name} needs to
#' be specified additionally. If of class \code{\link{SL}} or \code{\link{schmid}}, only
#' the arguments \code{factor_corres} and \code{cormat} need to be specified
#' additionally.
#' @param g_name character. The name of the general factor from the lavaan bifactor
#' solution. This needs only be specified if \code{model} is a class \code{\link{lavaan}}
#' object.
#' @param var_names character. A vector with subtest names in the order
#' of the rows from the SL solution. This needs only be specified if \code{model}
#' is left \code{NULL}
#' @param fac_names character. An optional vector of group factor names in the
#' order of the columns of the SL solution.
#' @param factor_corres numeric. A vector that indicates which variable corresponds
#' to which group factor. Must be in the same order as the SL solution. For example
#' c(3, 3, 3, 1, 1, 2, 2) if the first three variables load on the third group
#' factor of the SL solution, the next two on the first group factor and the
#' last two on the second group factor. If a variable should not be assigned to
#' any group factor, insert a zero at its position (e.g. c(3, 3, 0, 1, 1, 2, 2),
#' the third variable has no corresponding group factor).
#' @param type character. Either \code{"EFAtools"} (default) or \code{"psych"}.
#' @param g_load numeric. A vector of general factor loadings from an SL solution.
#' This needs only be specified if \code{model} is left \code{NULL}.
#' @param s_load matrix. A matrix of group factor loadings from an SL solution. This needs
#' only be specified if \code{model} is left \code{NULL}
#' @param u2 numeric. A vector of uniquenesses from an SL solution. This needs
#' only be specified if \code{model} is left \code{NULL}.
#' @param cormat matrix. A correlation matrix to be used when \code{type = "psych"}.
#' If left \code{NULL}, the correlation matrix is found based on the
#' pattern matrix and Phi using \code{\link[psych:factor.model]{psych::factor.model}}.
#' If the correlation matrix is available, \code{cormat} should be specified instead
#' of \code{Phi} and \code{pattern}.
#' @param Phi matrix. Factor intercorrelations from an oblique factor solution.
#' @param pattern matrix. Pattern coefficients from an oblique factor solution.
#' @param variance character. If \code{"correlation"}, then total variances for the whole
#' scale as well as for the subscale composites are calculated based on the correlation
#' matrix. If \code{"sums_load"}, then total variances are calculated using the
#' squared sums of general factor loadings and group factor loadings and
#' the sum of uniquenesses (see details).
#'
#'
#' @section How to combine arguments:
#'
#' If \code{model} is specified and of class \code{\link{lavaan}},
#' no other arguments need to be specified.
#' If \code{model} is of class
#' \code{\link{SL}} or \code{\link{schmid}}, only the argument \code{factor_corres}
#' needs to be specified additionally. There is, however, the option to reproduce
#' Watkins' Omega or \code{\link[psych:omega]{psych::omega}} results by setting the
#' \code{type} argument to \code{"Watkins"} or \code{"psych"}.
#' If \code{model = NULL} and \code{type = "EFAtools"}(default), the arguments
#' \code{var_names}, \code{factor_corres}, \code{g_load}, \code{s_load}, and \code{u2}
#' need to be specified.
#' If \code{type = "psych"} or \code{type = "EFAtools"}, either \code{cormat}
#' (recommended) or \code{Phi} and \code{pattern} must be specified.
#' Additionally, the argument \code{factor_corres} should be left NULL to
#' replicate \code{\link[psych:omega]{psych::omega}}
#' results, where variable-to-factor correspondences are found by taking the highest
#' group factor loading for each variable as the relevant group factor loading.
#' If \code{type = "Watkins"}, the \code{u2} argument
#' should be left \code{NULL} to replicate results from Watkins' Omega program, where
#' uniquenesses are found based on the general factor loadings and relevant group factor
#' loadings only. If, however, the argument \code{u2} is specified, the specified
#' uniquenesses are taken, with a warning.
#'
#' @section Calculation of omega for different types:
#'
#' The main differences between the types concern the calculation of the total
#' variance (for the whole scale as well as the subscale composites) as well as
#' the finding of variable-to-factor correspondences. The former aspect
#' can also be controlled individually by specifying the variance argument, the
#' latter by specifying the factor_corres argument.
#' For \code{type = "EFAtools"}, total variances are found using the correlation
#' matrix, and variable-to-factor correspondences have to be specified manually.
#' The only difference for \code{type = "psych"} is that it takes the highest
#' group factor loading for each variable as the relevant group factor loading.
#' To mimik results from Watkins' Omega program, for \code{type = "Watkins"}
#' for each variable only the general factor loading and the relevant group-factor
#' loadings according to the specified variable-to-factor correspondences is
#' taken into account. The other loadings are set to zero. Uniquenesses are found
#' based on these two loadings per variable only and total variance is calculated
#' based on all using the squared sums of general loadings and group factor loadings
#' and the sum of these uniquenesses.
#'
#' @return A matrix with omegas for the whole scale and for the subscales.
#' \item{tot}{Omega total.}
#' \item{hier}{Omega hierarchical.}
#' \item{sub}{Omega subscale.}
#'
#' @export
#'
#' @source
#' @example
#' \dontrun{
#' ## Use with a lavaan output
#' library(lavaan)
#'
#' # Create and fit model in lavaan
#' mod <- 'AVR =~ GS + PL + CM + EP
#'         PS =~ b*TC + b*CB
#'         ASTM =~ c*NL + c*NLM
#'         VSSTM =~ d*GF + d*RGF
#'         SLTM =~ CA + OP + RS + DP
#'         g =~ GS + PL + TC + CB + NL + NLM + GF + RGF + CM + EP + CA + OP +
#'              RS + DP'
#' SD <- c(3.20, 3.18, 3.13, 3.15, 3.20, 3.11, 3.06, 3.07, 3.14, 3.18, 3.22,
#'         3.16, 3.16, 3.08) # Taken from Table A1 in Grieder & Grob (2019)
#' IDS2_cov <- lavaan::cor2cov(IDS2_R, sds = SD, names = colnames(IDS2_R))
#' fit <- cfa(mod, sample.cov = IDS2_cov, sample.nobs = 1991, estimator = "ml",
#'            orthogonal = TRUE)
#'
#' # Calculate omega
#' OMEGA(fit)
#'
#' # Use with an output from the SL function, with type = "EFAtools" (default)
#'
#' efa_mod <- EFA(IDS2_R,  N = 1991, n_factors = 5, type = "EFAtools",
#'                method = "PAF", rotation = "promax")
#' sl_mod <- SL(efa_mod, type = "EFAtools", method = "PAF")
#'
#' OMEGA(sl_mod, factor_corres = c(1, 1, 5, 5, 2, 2, 4, 4, 1, 1, 3, 3, 3, 4),
#'       type = "EFAtools")
#'
#' # Manually specify components (useful, if omegas should be computed for a SL
#' # or bifactor solution found with another program)
#' OMEGA()
#' }
OMEGA <- function(model = NULL, g_name = NULL, var_names = NULL, fac_names = NULL,
                  factor_corres = NULL, g_load = NULL,
                  s_load = NULL, u2 = NULL, Phi = NULL, pattern = NULL,
                  cormat = NULL, variance = NULL, type = c("EFAtools", "psych")){

  type <- match.arg(type)

  if(!is.null(model) & (!is.null(var_names) || !is.null(g_load) || !is.null(s_load)
                        || !is.null(u2))){

    warning("You entered a model and specified at least one of the arguments
            var_names, g_load, s_load, or u2. These arguments are ignored.
            To use specific values for these, leave model = NULL and specify
            all arguments separately.")

  }

  if(all(class(model) == "lavaan")){

    .OMEGA_LAVAAN(model = model)

  } else if(all(class(model) == c("psych", "schmid")) || all(class(model) == "SL")) {

    if(is.null(factor_corres & all(class(model) == c("psych", "schmid")))){

      stop("Please specify the argument 'factor_corres'")

    } else {

      .OMEGA_FLEX(model = model, var_names = var_names, fac_names = fac_names,
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

          .OMEGA_FLEX(model = model, var_names = var_names, fac_names = fac_names,
                 factor_corres = factor_corres, g_load = g_load, s_load = s_load,
                 u2 = u2, Phi = Phi, pattern = pattern, cormat = cormat,
                 variance = variance, type = type)

        }

      }
}
