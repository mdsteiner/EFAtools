#' McDonald's omega
#'
#' This function finds omega total, hierarchical, and subscale, as well as additional
#' model-based indices of interpretive relevance (H index, ECV, PUC)
#' from a Schmid-Leiman (SL) solution or lavaan single factor, second-order (see below),
#' or bifactor solution. The SL-based omegas can either be found from a
#' [psych::schmid()], [SL()], or,
#' in a more flexible way, by leaving
#' `model = NULL` and specifying additional arguments. By setting the
#' `type` argument, results from [psych::omega()]
#' can be reproduced.
#'
#' @param model class [SL()], class `schmid`, or class
#' `lavaan` object. That is, an output object from [SL()] or
#' [psych::schmid()], or a `lavaan` fit object with a
#' single factor, second-order, or bifactor solution. If of class `lavaan`,
#' only `g_name` needs to be specified additionally. If of class
#' [SL()] or `schmid`, only the arguments `factor_corres`
#' and `cormat` need to be specified additionally.
#' @param type character. Either `"EFAtools"` (default) or `"psych"`
#' (see details)
#' @param g_name character. The name of the general factor from the lavaan solution.
#' This needs only be specified if `model` is a `lavaan` second-order
#' or bifactor solution. Default is "g".
#' @param group_names character. An optional vector of group names. The length
#' must correspond to the number of groups for which the `lavaan` model
#' was fitted.
#' @param add_ind logical. Whether additional indices (H index, ECV, PUC) should
#' be calculated or not (see details for these indices). If FALSE, only omegas
#' are returned. Default is `TRUE`.
#' @param factor_corres matrix. A logical matrix or a numeric matrix containing
#' 0's and 1's that indicates which variable corresponds to which group factor.
#' Must have the same dimensions as the matrix of group factor loadings from the
#' SL solution. Cross-loadings are allowed here. See examples for use.
#' @param var_names character. A vector with subtest names in the order
#' of the rows from the SL solution. This needs only be specified if `model`
#' is left `NULL`.
#' @param fac_names character. An optional vector of group factor names in the
#' order of the columns of the SL solution. If left `NULL`, names of the
#' group factors from the entered solution are taken.
#' @param g_load numeric. A vector of general factor loadings from an SL solution.
#' This needs only be specified if `model` is left `NULL`.
#' @param s_load matrix. A matrix of group factor loadings from an SL solution.
#' This needs only be specified if `model` is left `NULL`.
#' @param u2 numeric. A vector of uniquenesses from an SL solution. This needs
#' only be specified if `model` is left `NULL`.
#' @param cormat matrix. A correlation matrix to be used when
#' `variance = "correlation"`. If left `NULL` and an [SL()]
#' output is entered in `model`, the correlation matrix is taken from the
#' output. If left `NULL` and a [psych::schmid()]
#' output is entered, the correlation matrix will be found based on the pattern
#' matrix and Phi from the [psych::schmid()] output
#' using [psych::factor.model()].
#' If left `NULL` and model is also left `NULL`, the correlation matrix
#' is found based on the pattern matrix and Phi entered. However, if the
#' correlation matrix is available, `cormat` should be specified instead
#' of `Phi` and `pattern`.
#' @param pattern matrix. Pattern coefficients from an oblique factor solution.
#' This needs only be specified if `model` is left `NULL`,
#' `variance = "correlation"` and `cormat` is also left `NULL`.
#' @param Phi matrix. Factor intercorrelations from an oblique factor solution.
#' This needs only be specified if `model` is left `NULL`,
#' `variance = "correlation"` and `cormat` is also left `NULL`.
#' @param variance character. If `"correlation"` (default), then total
#' variances for the whole scale as well as for the subscale composites are
#' calculated based on the correlation
#' matrix. If `"sums_load"`, then total variances are calculated using the
#' squared sums of general factor loadings and group factor loadings and
#' the sum of uniquenesses (see details).
#'
#' @details ## What this function does
#'
#' This function calculates McDonald's omegas (McDonald, 1978, 1985, 1999),
#' the H index (Hancock & Mueller, 2001), the explained common variance (ECV;
#' Sijtsma, 2009), and the percent of uncontaminated correlations (PUC; Bonifay
#' et al., 2015; Reise et al., 2013).
#'
#' All types of omegas (total, hierarchical, and subscale) are calculated for
#' the general factor as well as for the subscales / group factors (see, e.g.,
#' Gignac, 2014; Rodriguez et al., 2016a, 2016b). Omegas refer to the correlation
#' between a factor and a unit-weighted composite score and thus the
#' true score variance in a unit-weighted composite based on the respective
#' indicators. Omega total is the total true score variance in a composite.
#' Omega hierarchical is the true score variance in a composite that is attributable
#' to the general factor, and omega subscale is the true score variance in a
#' composite attributable to all subscales / group factors (for the whole scale)
#' or to the specific subscale / group factor (for subscale composites).
#'
#' The H index (also construct reliability or replicability index) is the
#' correlation between an optimally-weighted composite score
#' and a factor (Hancock & Mueller, 2001; Rodriguez et al., 2016a, 2016b). It, too,
#' can be calculated for the whole scale / general factor as well as for the
#' subscales / grouup factors. Low values indicate that a latent variable is not well
#' defined by its indicators.
#'
#' The ECV (Sijtsma, 2009, Rodriguez et al., 2016a, 2016b) is the ratio of the
#' variance explained by the general factor and the variance explained by the
#' general factor and the group factors.
#'
#' The PUC (Bonifay et al., 2015; Reise et al., 2013, Rodriguez et al., 2016a,
#' 2016b) refers to the proportion
#' of correlations in the underlying correlation matrix that is not contaminated
#' by variance of both the general factor and the group factors (i.e., correlations
#' between indicators from different group factors, which reflect only general
#' factor variance). The higher the PUC, the more similar a general factor from
#' a multidimensional model will be to the single factor from a unidimensional
#' model.
#'
#' ## How to use this function
#'
#' If `model` is a `lavaan` second-order or bifactor solution,
#' only the name of the general factor from the lavaan model needs to be specified
#' additionally with the `g_name` argument. It is then determined whether this
#' general factor is a second-order factor (second-order model with one second-order
#' factor assumed) or a breadth factor (bifactor model assumed). Please note that
#' this function only works for second-order models if they contain no more than
#' one second-order factor. In case of a second-order solution, a
#' Schmid-Leiman transformation is performed on the first- and second-order loadings
#' and omega coefficents are obtained from the transformed (orthogonalized) solution
#' (see [SL()] for more information on Schmid-Leiman transformation).
#' There is also the possibility to enter a `lavaan` single factor solution.
#' In this case, `g_name` is not needed. Finally, if a solution from a
#' `lavaan` multiple group analysis is entered, the indices are computed for
#' each group.
#' The type argument is not evaluated if `model` is of class
#' `lavaan`.
#'
#' If `model` is of class [SL()] or
#' [psych::schmid()] only the
#' `type` and, depending on the type (see below), the `factor_corres`
#' arguments need to be specified additionally. If model is of class
#' [psych::schmid()] and `variance = "correlation"`
#' (default), it is
#' recommended to also provide the original correlation matrix in `cormat`
#' to get more accurate results. Otherwise, the correlation matrix will be found
#' based on the pattern matrix and Phi from the
#' [psych::schmid()] output
#' using the [psych::factor.model()] function.
#'
#' If `model = NULL`, the arguments `type`, `factor_corres`
#' (depending on the type, see below), `var_names`, `g_load`, `s_load`,
#' and `u2` and either `cormat` (recommended) or `Phi` and
#' `pattern` need to be specified. If `Phi` and `pattern` are
#' specified instead of `cormat`, the correlation matrix is found using
#' the [psych::factor.model()] function.
#'
#' The only difference between `type = "EFAtools"` and `type = "psych"`
#' is the determination of variable-to-factor correspondences. `type = "psych"`
#' reproduces the [psych::omega()] results, where
#' variable-to-factor correspondences are found by taking the highest
#' group factor loading for each variable as the relevant group factor loading.
#' To do this, `factor_corres` must be left `NULL`.
#'
#' The calculation of the total variance (for the whole scale as well as the
#' subscale composites) can also be controlled in this function using the
#' `variance` argument. For both types---`"EFAtools"` and `"psych"`
#' ---`variance` is set to `"correlation"` by default, which means that
#' total variances are found using the correlation matrix. If
#' `variance = "sums_load"` the total variance is calculated using the
#' squared sums of general loadings and group factor loadings and the sum of the
#' uniquenesses. This will only get comparable results to
#' `variance = "correlation"` if no cross-loadings are present and simple
#' structure is well-achieved in general with the SL solution (i.e., the
#' uniquenesses should capture almost all of the variance not explained by the
#' general factor and the variable's allocated group factor).
#'
#'
#' @return If found for an SL or `lavaan` second-order of bifactor solution
#' without multiple groups:
#' A matrix with omegas for the whole scale and for the subscales and (only if
#' `add_ind = TRUE`) with the H index, ECV, and PUC.
#' \item{tot}{Omega total.}
#' \item{hier}{Omega hierarchical.}
#' \item{sub}{Omega subscale.}
#' \item{H}{H index.}
#' \item{ECV}{Explained common variance.}
#' \item{PUC}{Percent of uncontaminated correlations.}
#'
#' If found for a `lavaan` single factor solution without multiple groups:
#' A (named) vector with omega total and (if `add_ind = TRUE`) the H index
#' for the single factor.
#'
#' If found for a `lavaan` output from a multiple group analysis: A list
#' containing the output described above for each group.
#'
#' @source McDonald, R. P. (1978). Generalizability in factorable domains: ‘‘Domain
#' validity and generalizability’’. Educational and Psychological Measurement,
#' 38, 75–79.
#' @source McDonald, R. P. (1985). Factor analysis and related methods. Hillsdale,
#' NJ: Erlbaum.
#' @source McDonald, R. P. (1999). Test theory: A unified treatment. Mahwah,
#' NJ: Erlbaum.
#' @source Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016a). Applying bifactor
#' statistical indices in the evaluation of psychological measures. Journal of
#' Personality Assessment, 98, 223-237.
#' @source Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016b). Evaluating
#' bifactor models: Calculating and interpreting statistical indices.
#' Psychological Methods, 21, 137-150.
#' @source Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct reliability
#' within latent variable systems. In R. Cudeck, S. du Toit, & D. Sörbom (Eds.),
#' Structural equation modeling: Present and future—A Festschrift in honor of Karl
#' Jöreskog (pp. 195–216). Lincolnwood, IL: Scientific Software International.
#' @source Sijtsma, K. (2009). On the use, the misuse, and the very limited usefulness
#' of Cronbach’s alpha. Psychometrika, 74, 107–120.
#' @source Reise, S. P., Scheines, R., Widaman, K. F., & Haviland, M. G. (2013).
#' Multidimensionality and structural coefficient bias in structural equation
#' modeling: A bifactor perspective. Educational and Psychological Measurement,
#' 73, 5–26.
#' @source Bonifay, W. E., Reise, S. P., Scheines, R., & Meijer, R. R. (2015).
#' When are multidimensional data unidimensional enough for structural equation
#' modeling?: An evaluation of the DETECT multidimensionality index. Structural
#' Equation Modeling, 22, 504—516.
#' @source Gignac, G. E. (2014). On the Inappropriateness of Using Items to
#' Calculate Total Scale Score Reliability via Coefficient Alpha for Multidimensional
#' Scales. European Journal of Psychological Assessment, 30, 130-139.
#'
#' @export
#'
#' @examples
#' \donttest{
#' ## Use with lavaan outputs
#' if (requireNamespace("lavaan", quietly = TRUE)) {
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
#' # Compute omegas and additional indices for bifactor solution
#' OMEGA(fit_bi, g_name = "g")
#'
#' # Compute only omegas
#' OMEGA(fit_bi, g_name = "g", add_ind = FALSE)
#'
#' # Create and fit second-order model in lavaan (assume all variables have SDs of 1)
#' mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
#'         F2 =~ V7 + V8 + V9 + V10 + V11 + V12
#'         F3 =~ V13 + V14 + V15 + V16 + V17 + V18
#'         g =~ F1 + F2 + F3'
#' fit_ho <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
#'                       sample.nobs = 500, estimator = "ml")
#'
#' # Compute omegas and additional indices for second-order solution
#' OMEGA(fit_ho, g_name = "g")
#' }
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
                  group_names = NULL, add_ind = TRUE, factor_corres = NULL,
                  var_names = NULL, fac_names = NULL,
                  g_load = NULL, s_load = NULL, u2 = NULL, cormat = NULL,
                  pattern = NULL, Phi = NULL, variance = c("correlation",
                                                           "sums_load")){

  # Perform argument checks
  type <- match.arg(type)
  checkmate::assert_string(g_name)
  checkmate::assert_logical(add_ind)
  checkmate::assert_character(group_names, null.ok = TRUE)
  # Check for factor_corres in OMEGA_helper
  checkmate::assert_character(var_names, null.ok = TRUE)
  checkmate::assert_character(fac_names, null.ok = TRUE)
  checkmate::assert_numeric(g_load, null.ok = TRUE)
  if(!is.null(s_load) && !inherits(s_load, c("matrix", "SLLOADINGS"))){

    cli::cli_abort(
      c("Invalid {.arg s_load}.",
        "i" = "Leave it {.code NULL} for a model input, or supply a Schmid-Leiman loading matrix of class {.cls matrix} or {.cls SLLOADINGS}."),
      class = "efa_omega_bad_s_load"
    )

  }
  checkmate::assert_numeric(u2, null.ok = TRUE)
  checkmate::assert_matrix(cormat, null.ok = TRUE)
  if(!is.null(pattern) && !inherits(pattern, c("matrix", "loadings", "LOADINGS"))){

    cli::cli_abort(
      c("Invalid {.arg pattern}.",
        "i" = "Leave it {.code NULL}, or supply a pattern matrix from an oblique solution of class {.cls matrix}, {.cls loadings}, or {.cls LOADINGS}."),
      class = "efa_omega_bad_pattern"
    )

  }
  checkmate::assert_matrix(Phi, null.ok = TRUE)
  variance <- match.arg(variance)

  # Determine which function to use
  if(!is.null(model) & (!is.null(var_names) || !is.null(g_load) || !is.null(s_load)
     || !is.null(u2))){

    cli::cli_warn(
      c("You entered a {.arg model} and also specified {.arg var_names}, {.arg g_load}, {.arg s_load}, or {.arg u2}; these are ignored.",
        "i" = "To use specific values, leave {.code model = NULL} and specify all arguments separately."),
      class = "efa_omega_model_args_ignored"
    )

  }

  if(inherits(model, "lavaan")){

    .require_lavaan()

    .OMEGA_LAVAAN(model = model, g_name = g_name, group_names = group_names,
                  add_ind = add_ind)

  } else if(inherits(model, c("schmid", "SL"))) {

     .OMEGA_FLEX(model = model, type = type, factor_corres = factor_corres,
                 var_names = var_names, fac_names = fac_names, g_load = g_load,
                 s_load = s_load, u2 = u2, cormat = cormat, pattern = pattern,
                 Phi = Phi, variance = variance, add_ind = add_ind)
  } else {

      if(!is.null(model)){

        cli::cli_abort(
          c("Invalid {.arg model}.",
            "i" = "Enter a lavaan, {.fn psych::schmid}, or SL object, or specify {.arg var_names}, {.arg g_load}, and {.arg s_load}."),
          class = "efa_omega_invalid_model"
        )

      } else if(is.null(var_names) || is.null(g_load) || is.null(s_load) || is.null(u2)){

      cli::cli_abort("Specify all of {.arg var_names}, {.arg g_load}, {.arg s_load}, and {.arg u2}.",
                     class = "efa_omega_missing_args")

        } else {

          .OMEGA_FLEX(model = model, type = type, factor_corres = factor_corres,
                      var_names = var_names, fac_names = fac_names, g_load = g_load,
                      s_load = s_load, u2 = u2, cormat = cormat, pattern = pattern,
                      Phi = Phi, variance = variance, add_ind = add_ind)

        }

      }
}
