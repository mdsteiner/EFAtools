#' McDonald's omega
#'
#' This function finds omega total, omega hierarchical, and omega subscale
#' from a Schmid-Leiman (SL) solution or lavaan single factor, second-order (see below),
#' or bifactor solution. The SL-based omegas can either be found from a
#' \code{\link[psych:schmid]{psych::schmid}}, \code{\link{SL}}, or,
#' in a more flexible way, by leaving
#' \code{model = NULL} and specifying additional arguments. By setting the
#' \code{type} argument, results from \code{\link[psych:omega]{psych::omega}}
#' can be reproduced.
#'
#' @param model class \code{\link{SL}}, class \code{\link{schmid}}, or class
#' \code{lavaan} object. That is, an output object from \code{\link{SL}} or
#' \code{\link[psych:schmid]{psych::schmid}}, or a \code{lavaan} fit object with a
#' single factor, second-order, or bifactor solution. If of class \code{lavaan},
#' only \code{g_name} needs to be specified additionally. If of class
#' \code{\link{SL}} or \code{\link{schmid}}, only the arguments \code{factor_corres}
#' and \code{cormat} need to be specified additionally.
#' @param type character. Either \code{"EFAtools"} (default) or \code{"psych"}
#' (see details)
#' @param g_name character. The name of the general factor from the lavaan solution.
#' This needs only be specified if \code{model} is a \code{lavaan} second-order
#' or bifactor solution. Default is "g".
#' @param group_names character. An optional vector of group names. The length
#' must correspond to the number of groups for which the \code{lavaan} model
#' was fitted.
#' @param factor_corres matrix. A logical matrix or a numeric matrix containing
#' 0's and 1's that indicates which variable corresponds to which group factor.
#' Must have the same dimensions as the matrix of group factor loadings from the
#' SL solution. Cross-loadings are allowed here. See examples for use.
#' @param var_names character. A vector with subtest names in the order
#' of the rows from the SL solution. This needs only be specified if \code{model}
#' is left \code{NULL}.
#' @param fac_names character. An optional vector of group factor names in the
#' order of the columns of the SL solution. If left \code{NULL}, names of the
#' group factors from the entered solution are taken.
#' @param g_load numeric. A vector of general factor loadings from an SL solution.
#' This needs only be specified if \code{model} is left \code{NULL}.
#' @param s_load matrix. A matrix of group factor loadings from an SL solution.
#' This needs only be specified if \code{model} is left \code{NULL}.
#' @param u2 numeric. A vector of uniquenesses from an SL solution. This needs
#' only be specified if \code{model} is left \code{NULL}.
#' @param cormat matrix. A correlation matrix to be used when
#' \code{variance = "correlation"}. If left \code{NULL} and an \code{\link{SL}}
#' output is entered in \code{model}, the correlation matrix is taken from the
#' output. If left \code{NULL} and a \code{\link[psych:schmid]{psych::schmid}}
#' output is entered, the correlation matrix will be found based on the pattern
#' matrix and Phi from the \code{\link[psych:schmid]{psych::schmid}} output
#' using \code{\link[psych:factor.model]{psych::factor.model}}.
#' If left \code{NULL} and model is also left \code{NULL}, the correlation matrix
#' is found based on the pattern matrix and Phi entered. However, if the
#' correlation matrix is available, \code{cormat} should be specified instead
#' of \code{Phi} and \code{pattern}.
#' @param pattern matrix. Pattern coefficients from an oblique factor solution.
#' This needs only be specified if \code{model} is left \code{NULL},
#' \code{variance = "correlation"} and \code{cormat} is also left \code{NULL}.
#' @param Phi matrix. Factor intercorrelations from an oblique factor solution.
#' This needs only be specified if \code{model} is left \code{NULL},
#' \code{variance = "correlation"} and \code{cormat} is also left \code{NULL}.
#' @param variance character. If \code{"correlation"} (default), then total
#' variances for the whole scale as well as for the subscale composites are
#' calculated based on the correlation
#' matrix. If \code{"sums_load"}, then total variances are calculated using the
#' squared sums of general factor loadings and group factor loadings and
#' the sum of uniquenesses (see details).
#'
#' @details If \code{model} is a \code{lavaan} second-order or bifactor solution,
#' only the name of the general factor from the lavaan model needs to be specified
#' additionally with the \code{g_name} argument. It is then determined whether this
#' general factor is a second-order factor (second-order model assumed) or a breadth
#' factor (bifactor model assumed). In case of a second-order solution, a
#' Schmid-Leiman transformation is performed on the first- and second-order loadings
#' and omega coefficents are obtained from the transformed (orthogonalized) solution
#' (see \code{\link{SL}} for more information on Schmid-Leiman transformation).
#' There is also the possibility to enter a \code{lavaan} single factor solution.
#' In this case, \code{g_name} is not needed. Finally, if a solution from a
#' \code{lavaan} multiple group analysis is entered, the omegas are computed for
#' each group.
#' The type argument is not evaluated if \code{model} is of class
#' \code{lavaan}.
#'
#' If \code{model} is of class \code{\link{SL}} or
#' \code{\link[psych:schmid]{psych::schmid}} only the
#' \code{type} and, depending on the type (see below), the \code{factor_corres}
#' arguments need to be specified additionally. If model is of class
#' \code{\link[psych:schmid]{psych::schmid}} and \code{variance = "correlation"}
#' (default), it is
#' recommended to also provide the original correlation matrix in \code{cormat}
#' to get more accurate results. Otherwise, the correlation matrix will be found
#' based on the pattern matrix and Phi from the
#' \code{\link[psych:schmid]{psych::schmid}} output
#' using the \code{\link[psych:factor.model]{psych::factor.model}} function.
#'
#' If \code{model = NULL}, the arguments \code{type}, \code{factor_corres}
#' (depending on the type, see below), \code{var_names}, \code{g_load}, \code{s_load},
#' and \code{u2} and either \code{cormat} (recommended) or \code{Phi} and
#' \code{pattern} need to be specified. If \code{Phi} and \code{pattern} are
#' specified instead of \code{cormat}, the correlation matrix is found using
#' the \code{\link[psych:factor.model]{psych::factor.model}} function.
#'
#' The only difference between \code{type = "EFAtools"} and \code{type = "psych"}
#' is the determination of variable-to-factor correspondences. \code{type = "psych"}
#' reproduces the \code{\link[psych:omega]{psych::omega}} results, where
#' variable-to-factor correspondences are found by taking the highest
#' group factor loading for each variable as the relevant group factor loading.
#' To do this, \code{factor_corres} must be left \code{NULL}.
#'
#' The calculation of the total variance (for the whole scale as well as the
#' subscale composites) can also be controlled in this function using the
#' \code{variance} argument. For both types---\code{"EFAtools"} and \code{"psych"}
#' ---\code{variance} is set to \code{"correlation"} by default, which means that
#' total variances are found using the correlation matrix. If
#' \code{variance = "sums_load"} the total variance is calculated using the
#' squared sums of general loadings and group factor loadings and the sum of the
#' uniquenesses. This will only get comparable results to
#' \code{variance = "correlation"} if no cross-loadings are present and simple
#' structure is well-achieved in general with the SL solution (i.e., the
#' uniquenesses should capture almost all of the variance not explained by the
#' general factor and the variable's allocated group factor).
#'
#' @return If found for an SL or \code{lavaan} second-order of bifactor solution
#' without multiple groups:
#' A matrix with omegas for the whole scale and for the subscales.
#' \item{tot}{Omega total.}
#' \item{hier}{Omega hierarchical.}
#' \item{sub}{Omega subscale.}
#'
#' If found for a \code{lavaan} single factor solution without multiple groups:
#' A vector with omega total for the single factor.
#'
#' If found for a \code{lavaan} output from a multiple group analysis: A list
#' containing the output described above for each group.
#'
#' @source McDonald, R. P. (1978). Generalizability in factorable domains: ‘‘Domain
#' validity and generalizability’’. Educational and Psychological Measurement,
#' 38, 75–79.
#' @source McDonald, R. P. (1985). Factor analysis and related methods. Hillsdale,
#' NJ: Erlbaum.
#' @source McDonald, R. P. (1999). Test theory: A unified treatment. Mahwah,
#' NJ: Erlbaum.
#' @source Gignac, G. E. (2014). On the Inappropriateness of Using Items to
#' Calculate Total Scale Score Reliability via Coefficient Alpha for Multidimensional
#' Scales. European Journal of Psychological Assessment, 30, 130-139.
#'
#' @export
#'
#' @examples
#' \donttest{
#' ## Use with lavaan outputs
#'
#' # Create and fit bifactor model in lavaan (assume all variables have SDs of 1)
#' mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
#'         F2 =~ V7 + V8 + V9 + V10 + V11 + V12
#'         F3 =~ V13 + V14 + V15 + V16 + V17 + V18
#'         g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
#'              V13 + V14 + V15 + V16 + V17 + V18'
#' fit_bi <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
#'                       sample.nobs = 500, estimator = "ml", orthogonal = TRUE)
#'
#' # Compute omega for bifactor solution
#' OMEGA(fit_bi, g_name = "g")
#'
#' # Create and fit second-order model in lavaan (assume all variables have SDs of 1)
#' mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
#'         F2 =~ V7 + V8 + V9 + V10 + V11 + V12
#'         F3 =~ V13 + V14 + V15 + V16 + V17 + V18
#'         g =~ F1 + F2 + F3'
#' fit_ho <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
#'                       sample.nobs = 500, estimator = "ml")
#'
#' # Compute omega for second-order solution
#' OMEGA(fit_ho, g_name = "g")
#' }
#'
#' ## Use with an output from the SL function, with type EFAtools
#' efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                type = "EFAtools", method = "PAF", rotation = "promax")
#' sl_mod <- SL(efa_mod, type = "EFAtools", method = "PAF")
#'
#' # Two examples how to specify the indicator-to-factor correspondences:
#'
#' # Based on a specific salience threshold for the loadings (here: .20):
#' factor_corres_1 <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2
#'
#' # Or more flexibly (could also be TRUE and FALSE instead of 0 and 1):
#' factor_corres_2 <- matrix(c(rep(0, 12), rep(1, 6), rep(0, 6), rep(1, 6),
#'                          rep(0, 6), rep(1, 6), rep(0, 12)), ncol = 3,
#'                          byrow = FALSE)
#'
#' OMEGA(sl_mod, type = "EFAtools", factor_corres = factor_corres_1)
#'
#' ## Use with an output from the psych::schmid function, with type psych for
#' ## OMEGA
#' schmid_mod <- psych::schmid(test_models$baseline$cormat, nfactors = 3,
#'                             n.obs = 500, fm = "pa", rotate = "Promax")
#' # Find correlation matrix from phi and pattern matrix from psych::schmid output
#' OMEGA(schmid_mod, type = "psych")
#' # Use specified correlation matrix
#' OMEGA(schmid_mod, type = "psych", cormat = test_models$baseline$cormat)
#'
#' ## Manually specify components (useful if omegas should be computed for a SL
#' ## or bifactor solution found with another program)
#' ## As an example, we extract the elements from an SL output here. This gives
#' ## the same results as in the second example above.
#'
#' efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                type = "EFAtools", method = "PAF", rotation = "promax")
#' sl_mod <- SL(efa_mod, type = "EFAtools", method = "PAF")
#'
#' factor_corres <- matrix(c(rep(0, 12), rep(1, 6), rep(0, 6), rep(1, 6),
#'                         rep(0, 6), rep(1, 6), rep(0, 12)), ncol = 3,
#'                         byrow = FALSE)
#'
#' OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
#'       g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
#'       u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
#'       factor_corres = factor_corres)
#'
OMEGA <- function(model = NULL, type = c("EFAtools", "psych"), g_name = "g",
                  group_names = NULL, factor_corres = NULL,
                  var_names = NULL, fac_names = NULL,
                  g_load = NULL, s_load = NULL, u2 = NULL, cormat = NULL,
                  pattern = NULL, Phi = NULL, variance = c("correlation",
                                                           "sums_load")){

  # Perform argument checks
  type <- match.arg(type)
  checkmate::assert_string(g_name)
  checkmate::assert_character(group_names, null.ok = TRUE)
  # Check for factor_corres in OMEGA_helper
  checkmate::assert_character(var_names, null.ok = TRUE)
  checkmate::assert_character(fac_names, null.ok = TRUE)
  checkmate::assert_numeric(g_load, null.ok = TRUE)
  if(!is.null(s_load) && !inherits(s_load, c("matrix", "SLLOADINGS"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Specification of 's_load' was invalid. Please either leave this 'NULL' if you enter a model input or specify a matrix of loadings from a Schmid-Leiman solution of class matrix or SLLOADINGS.\n"))

  }
  checkmate::assert_numeric(u2, null.ok = TRUE)
  checkmate::assert_matrix(cormat, null.ok = TRUE)
  if(!is.null(pattern) && !inherits(pattern, c("matrix", "loadings", "LOADINGS"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Specification of 'pattern' was invalid. Please either leave this NULL or specify a matrix of pattern coefficients form an oblique factor solution of class matrix, loadings, or LOADINGS.\n"))

  }
  checkmate::assert_matrix(Phi, null.ok = TRUE)
  variance <- match.arg(variance)

  # Determine which function to use
  if(!is.null(model) & (!is.null(var_names) || !is.null(g_load) || !is.null(s_load)
     || !is.null(u2))){

    warning(crayon::yellow$bold("!"), crayon::yellow(" You entered a model and specified at least one of the arguments 'var_names', 'g_load', 's_load', or 'u2'. These arguments are ignored. To use specific values for these, leave model = NULL and specify all arguments separately.\n"))

  }

  if(inherits(model, "lavaan")){

    .OMEGA_LAVAAN(model = model, g_name = g_name, group_names = group_names)

  } else if(inherits(model, c("schmid", "SL"))) {

     .OMEGA_FLEX(model = model, type = type, factor_corres = factor_corres,
                 var_names = var_names, fac_names = fac_names, g_load = g_load,
                 s_load = s_load, u2 = u2, cormat = cormat, pattern = pattern,
                 Phi = Phi, variance = variance)
  } else {

      if(!is.null(model)){

        stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Invalid input for model. Either enter a lavaan, psych::schmid or SL object or specify the arguments 'var_names', 'g_load', and 's_load'.\n"))

      } else if(is.null(var_names) || is.null(g_load) || is.null(s_load) || is.null(u2)){

      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Please specify all of the following arguments: 'var_names', 'g_load', 's_load', 'u2'\n"))

        } else {

          .OMEGA_FLEX(model = model, type = type, factor_corres = factor_corres,
                      var_names = var_names, fac_names = fac_names, g_load = g_load,
                      s_load = s_load, u2 = u2, cormat = cormat, pattern = pattern,
                      Phi = Phi, variance = variance)

        }

      }
}
