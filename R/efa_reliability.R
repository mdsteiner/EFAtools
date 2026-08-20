#' Reliability and common-variance coefficients for a factor solution
#'
#' Computes model-based reliability coefficients for a factor solution:
#' McDonald's omega (total, hierarchical, and subscale), standardized Cronbach's
#' alpha, and the H index. For a bifactor solution, it also computes two
#' common-variance indices, ECV and PUC, for the general factor. The result is a
#' tidy, long-format table with one row per coefficient.
#'
#' The function reads many kinds of input: a Schmid-Leiman solution
#' ([efa_schmid_leiman()] or [psych::schmid()]), an oblique [efa_fit()]
#' (correlated-factors) solution, a `lavaan` fit (single-factor,
#' correlated-factors, second-order, or bifactor), a raw bifactor loading matrix,
#' an oblique pattern matrix given with its factor intercorrelations, or manually
#' supplied components.
#'
#' @details ## Coefficients
#'
#' The reliability coefficients are McDonald's omegas (McDonald, 1978, 1985,
#' 1999; Zinbarg et al., 2005, 2006, for omega hierarchical specifically),
#' standardized Cronbach's alpha (Cronbach, 1951), and the H index (construct
#' replicability; Hancock & Mueller, 2001).
#'
#' The omegas give the share of true score variance in a unit-weighted
#' composite. Omega total is the share due to all factors together. Omega
#' hierarchical is the share due to the general factor only. Omega subscale is
#' the share due to the group factors: for the whole scale, or for one specific
#' group factor in a subscale composite.
#'
#' Alpha is the standardized coefficient, computed from the correlation matrix
#' of the items. Where no such matrix is available -- for `lavaan` input, and
#' for components supplied without one -- alpha is computed from the
#' model-implied correlation matrix instead, so it then reflects the fitted
#' model rather than the raw data.
#'
#' The H index is the reliability of an optimally weighted composite. A low
#' value means the factor is not well defined by its indicators.
#'
#' All of these coefficients describe the raw sum of the variables as fitted,
#' without reverse-coding. If some items are keyed in the opposite direction --
#' for example, a reverse-worded item that was not reverse-scored -- they lower
#' that sum's true-score variance. The coefficients then look poor even though
#' the model fits well. The function warns when it detects this. Reverse-code
#' such items before fitting the solution (Flora, 2020).
#'
#' The sum these coefficients describe is not always the sum of the raw answers.
#' Polychoric and tetrachoric correlations describe the continuous latent
#' responses assumed to underlie ordinal answers; a `lavaan` fit that declares
#' its indicators `ordered` does the same. In that case the loadings, the
#' uniquenesses, the correlation matrix, and the coefficients are all on that
#' latent-response metric. They give the reliability of the unit-weighted sum of
#' the latent responses, which is not observed -- not the reliability of the
#' ordinal sum score the user actually computes from the answers. The two can
#' differ substantially, especially where the answers use few categories or are
#' strongly skewed. Green and Yang (2009) give an omega for the ordinal sum
#' score itself, computed from the fitted model and its thresholds; this
#' package does not compute it. Pearson correlations raise no such distinction,
#' because their metric is the answers as scored.
#'
#' Omega total is lower when a solution reproduces the observed correlations
#' poorly, because residual covariance does not count as true score.
#' [psych::omega()] computes the whole-scale omega total differently:
#' residually, from the observed total-score variance. The two agree when the
#' model reproduces the correlations exactly, and diverge otherwise.
#'
#' The three coefficients are not generally additive. Omega total need not
#' equal omega hierarchical plus omega subscale, except on the whole-scale row
#' under `variance = "sums_load"`.
#'
#' A single-factor solution is scored as such on every input route, and
#' returns omega total, alpha, and the H index. Alpha assumes essentially
#' tau-equivalent items, an assumption nested within a one-factor model -- so a
#' single factor is the case where reporting alpha is defensible, not merely
#' possible.
#'
#' The other coefficients are omitted because a single factor does not define
#' them: omega subscale is the variance due to the group factors, and there
#' are none; omega hierarchical would equal omega total, since the one factor
#' accounts for all common variance; and ECV and PUC would both be 1 by
#' construction, which reflects the number of factors in the model rather than
#' evidence of unidimensionality.
#'
#' Each coefficient answers a different question:
#'
#' - Omega hierarchical: can the total score be read as a measure of one
#'   construct?
#' - Omega subscale: does a subscale score add anything beyond the general
#'   factor?
#' - The H index: is a factor well defined by its indicators?
#' - ECV together with PUC: is a unidimensional model defensible?
#'
#' Alpha assumes essentially tau-equivalent items. Factor analysis instead
#' yields congeneric solutions, for which alpha is only a lower bound. For a
#' multidimensional scale, alpha is rarely the coefficient to report (Gignac,
#' 2014).
#'
#' Composite reliability and average variance extracted (AVE) are not among
#' the coefficients this function computes. Use `semTools::compRelSEM()` and
#' `semTools::AVE()` to compute them from a `lavaan` fit.
#'
#' The common-variance indices ECV and PUC (Bonifay et al., 2015; Reise et al.,
#' 2013; Rodriguez et al., 2016a, 2016b) describe the general factor, so they
#' are reported for the general factor only. ECV is the share of the common
#' variance explained by the general factor. PUC is the proportion of
#' correlations that reflect general-factor variance alone -- correlations
#' between indicators of different group factors. The higher the PUC, the more
#' the general factor resembles the single factor of a unidimensional model.
#'
#' ## Input
#'
#' `efa_reliability()` reads several kinds of input, illustrated in the
#' examples below.
#'
#' - **An oblique [efa_fit()] solution** is scored as the correlated-factors
#'   model it is. It has no general factor, so it omits the bifactor indices
#'   (omega hierarchical, ECV, and PUC). Its whole-scale row is labelled
#'   `"total"` rather than `"g"`.
#' - **An [efa_schmid_leiman()] solution** -- or a `schmid` object from
#'   [psych::schmid()], or a raw bifactor loading matrix with the general
#'   factor in its first column -- is scored as a bifactor solution, with the
#'   whole-scale row labelled `"g"`. Indicator-to-factor correspondences come
#'   from `factor_map` if supplied (see `factor_map` below). Without one, an
#'   [efa_schmid_leiman()] or [psych::schmid()] solution is mapped by each
#'   variable's strongest group loading. A raw bifactor matrix instead
#'   defaults to its exact zero pattern -- supply `factor_map` explicitly if
#'   that matrix has no exact zeros (e.g. an estimated rather than a target
#'   matrix).
#' - **A `lavaan` fit** -- single-factor, correlated-factors, second-order, or
#'   bifactor -- is scored per its structure. The structure is detected
#'   automatically, not taken from `g_name`. The general-factor coefficients
#'   need the latent factors to be uncorrelated: fit a bifactor model with
#'   `orthogonal = TRUE` (`lavaan`'s default leaves the factors correlated),
#'   and leave a second-order model's first-order factor covariances at zero.
#'   A fit whose factors correlate is rejected rather than scored as though
#'   they did not. `variance` is not used for `lavaan` input: its composite
#'   variances are always model-implied, computed separately per group for a
#'   multiple-group fit.
#' - **A pattern matrix with its `Phi`** -- for example, the loadings and
#'   `Phi` of an oblique [efa_fit()] solution -- is scored as the
#'   correlated-factors solution it represents. Uniquenesses are derived from
#'   the loadings under `Phi`, unless `u2` is given. Supply the pattern
#'   matrix, not the structure matrix -- the structure matrix would be read
#'   as a pattern too, giving the wrong result. Unlike a fitted solution, this
#'   pair carries no correlation matrix. `variance = "correlation"` then
#'   needs one, supplied in `cormat`, and errors without it. `variance =
#'   "sums_load"` needs none.
#'
#' A one-factor solution -- a one-column loading matrix, a single-factor
#' `lavaan` fit, a single-factor [efa_fit()] solution, or single-factor
#' components -- is scored the same way regardless of route (see Coefficients
#' above for what it returns). Its row is labelled with the factor's own name,
#' with `fac_names` if supplied, or with `"F1"` if neither names it; it is
#' never labelled `"g"`.
#'
#' Manually supplied components (`g_load`, `s_load`, `u2`, `var_names`, `Phi`)
#' follow the same reading rules as the matrix routes above. Do not pass a
#' correlation matrix as `s_load` (or as `model`, for the matrix routes) --
#' supply it as `cormat` instead.
#'
#' @param model an [efa_schmid_leiman()], `schmid` ([psych::schmid()]), [efa_fit()]
#'   (oblique), or `lavaan` object; a raw bifactor loading matrix (general factor
#'   first), the loading table of an [efa_schmid_leiman()] solution, or the pattern
#'   matrix of an oblique solution together with its `Phi`; or `NULL` to
#'   supply the components manually via `g_load`, `s_load`, `u2`, and `var_names`.
#' @param coefficients character. An optional subset of the coefficients to
#'   return, any of `"omega_total"`, `"omega_hierarchical"`, `"omega_subscale"`,
#'   `"alpha"`, `"H"`, `"ECV"`, and `"PUC"`. Default `NULL` returns all of them.
#' @param g_name character. The name of the general factor in the `lavaan`
#'   solution. Only needed for a `lavaan` second-order or bifactor fit. A fit in
#'   which every variable loads on a single factor has no general factor, and
#'   reads none. No other input is affected by it. Default is `"g"`.
#' @param group_names character. An optional vector of group names for a `lavaan`
#'   multiple-group fit. Its length must match the number of groups. Not used for
#'   any other input -- including a single-group `lavaan` fit, since every one of
#'   those is scored as a single ungrouped solution with no group label.
#' @param factor_map matrix. A logical or 0/1 matrix indicating which variable
#'   belongs to which group factor, with the same dimensions as the group
#'   loading matrix (cross-loadings are allowed). Match its columns to the
#'   group factors by position -- a map given in a different factor order than
#'   the solution still runs, but produces meaningless subscale coefficients.
#'   The function warns if a mapped item loads weakly on its assigned factor
#'   and more strongly on another one. If `NULL` (default), each variable is
#'   assigned to the group factor on which it loads most strongly. Not used
#'   for `lavaan` input.
#' @param variance character. The total-variance denominator for the
#'   coefficients. `"correlation"` (default) takes each composite's variance
#'   from the correlation matrix, giving the observed-score omega. `"sums_load"`
#'   uses the model-implied composite variance from the loadings and the
#'   uniquenesses instead; it needs no correlation matrix, so it is the way to
#'   score a bare loading matrix or manual components given without one.
#'   `lavaan` input fixes the convention: its composite variances are always
#'   model-implied, and include any freed residual covariance as well as the
#'   residual variances. Neither convention changes the metric the coefficients
#'   are on -- with polychoric or tetrachoric correlations, or an `ordered`
#'   `lavaan` fit, both describe the latent-response composite rather than the
#'   ordinal sum score (see Details).
#' @param var_names character. Subtest names in the row order of the loadings.
#'   Only needed when `model` is `NULL`.
#' @param fac_names character. An optional vector of group-factor names in the
#'   column order of the loadings. Taken from the input if `NULL`. A
#'   single-factor solution has no group factors; `fac_names` then labels its
#'   one factor instead. If neither `fac_names` nor the solution names that
#'   factor, it is labelled `"F1"`. Not used for `lavaan` input, whose factor
#'   labels come from the model syntax.
#' @param g_load numeric. General-factor loadings. Only needed when `model` is
#'   `NULL`.
#' @param s_load matrix. Group-factor loadings. Only needed when `model` is
#'   `NULL`.
#' @param u2 numeric. Uniquenesses. Only needed when `model` is `NULL`, or to
#'   override the communality-based default for a loading matrix. Under
#'   `variance = "correlation"`, the coefficients follow from the loadings and
#'   the correlation matrix alone, so `u2` only enters the check for an improper
#'   solution. Under `"sums_load"`, it is part of every composite's variance.
#' @param cormat matrix or data.frame. A correlation matrix used when
#'   `variance = "correlation"`. It must hold the same variables as the
#'   solution: named variables in a different order are reordered to match; a
#'   different set of variables, or a different number of them, is an error.
#'   The matrix must use the solution's own variable names -- for manually
#'   supplied components, that means the row names of `s_load`, not
#'   `var_names`, which only labels the output. Without row names, the
#'   variables are matched by position, so supply the matrix in the row order
#'   of the loadings. If `NULL`, it is taken from the input object, or
#'   reconstructed from `pattern` and `Phi` where possible. Supply the matrix
#'   on the same metric as the loadings: a polychoric or tetrachoric matrix
#'   keeps the coefficients on the latent-response metric, while a Pearson
#'   matrix given with loadings fitted to a polychoric one mixes the two (see
#'   Details).
#' @param pattern matrix. Pattern coefficients from a separate oblique solution, used
#'   with `Phi` to reconstruct a correlation matrix when `model` is `NULL`. Supply it for
#'   a Schmid-Leiman input, whose `s_load` holds the orthogonalized group loadings rather
#'   than the oblique ones. It is an alternative to `cormat`, not a supplement: giving
#'   both is an error.
#' @param Phi matrix. Factor intercorrelations. `NULL` (default) means
#'   uncorrelated group factors, as a Schmid-Leiman or bifactor solution has.
#'   It cannot be combined with a general factor: supplied together with a
#'   non-zero `g_load`, or with the loading table of an
#'   [efa_schmid_leiman()] solution, it is an error.
#'
#'   Supply it together with `s_load` (manually supplied components, with
#'   `g_load` zero throughout), or with a loading matrix of two or more
#'   factors in `model`, to score the input as a correlated-factors solution
#'   instead.
#'
#'   Without `Phi`, a loading matrix in `model` is read as a bifactor
#'   solution (general factor first). The one exception is a matrix that
#'   still carries the loading class [efa_fit()] gives its output -- a class
#'   that subsetting a matrix or reordering its rows drops. Such a matrix is
#'   rejected instead; supply `Phi` to score it as a correlated-factors
#'   solution.
#'
#'   With `Phi`, a matrix in `model` that carries no loading class is read as
#'   a correlated-factors solution, and a warning names this reading (drop
#'   `Phi` to read it as a bifactor matrix instead).
#'
#'   Given together with `pattern` instead, see `pattern` above. With a
#'   fitted solution already in `model`, `Phi` is ignored, with a warning,
#'   since that solution carries its own factor intercorrelations.
#'
#' @returns An object of class `efa_reliability`: a long-format data frame with
#'   one row per computed coefficient, with columns
#'   \item{coefficient}{the coefficient name (e.g. `"omega_total"`).}
#'   \item{level}{`"general"` for the general-factor row, `"total"` for the
#'     whole-scale row of a correlated-factors solution, and `"group"` for
#'     every other row. A single-factor solution has only the general-factor
#'     row.}
#'   \item{factor}{the factor label: `"g"` for the general factor of a
#'     solution with group factors. A single-factor solution's one factor
#'     takes its own name, or `"F1"` if neither the input nor `fac_names`
#'     names it. A correlated-factors solution has no general factor, so its
#'     first row is labelled `"total"` instead -- it describes the composite
#'     of every variable, not a factor of the model.}
#'   \item{group}{the group label, or `NA` for a single ungrouped solution.}
#'   \item{value}{the coefficient value.}
#'   Structurally undefined cells (for example, ECV and PUC on a group factor)
#'   are omitted. The object also carries a `settings` attribute (the
#'   total-variance convention used, and whether the solution has a general
#'   factor) and a `kind` attribute tagging each coefficient as a reliability
#'   coefficient or a common-variance index. It has a
#'   [print.efa_reliability()] method.
#'
#' @source McDonald, R. P. (1978). Generalizability in factorable domains: Domain
#'   validity and generalizability. Educational and Psychological Measurement, 38,
#'   75-79.
#' @source McDonald, R. P. (1985). Factor analysis and related methods. Hillsdale,
#'   NJ: Erlbaum.
#' @source McDonald, R. P. (1999). Test theory: A unified treatment. Mahwah,
#'   NJ: Erlbaum.
#' @source Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
#'   tests. Psychometrika, 16, 297-334.
#' @source Gignac, G. E. (2014). On the inappropriateness of using items to
#'   calculate total scale score reliability via coefficient alpha for
#'   multidimensional scales. European Journal of Psychological Assessment, 30,
#'   130-139.
#' @source Flora, D. B. (2020). Your coefficient alpha is probably wrong, but which
#'   coefficient omega is right? A tutorial on using R to obtain better reliability
#'   estimates. Advances in Methods and Practices in Psychological Science, 3,
#'   484-501.
#' @source Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
#'   structural equation modeling: An alternative to coefficient alpha. Psychometrika,
#'   74, 155-167.
#' @source Zinbarg, R. E., Revelle, W., Yovel, I., & Li, W. (2005). Cronbach's
#'   alpha, Revelle's beta, and McDonald's omega H: Their relations with each
#'   other and two alternative conceptualizations of reliability. Psychometrika,
#'   70, 123-133.
#' @source Zinbarg, R. E., Yovel, I., Revelle, W., & McDonald, R. P. (2006).
#'   Estimating generalizability to a latent variable common to all of a scale's
#'   indicators: A comparison of estimators for omega H. Applied Psychological
#'   Measurement, 30, 121-144.
#' @source Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
#'   reliability within latent variable systems. In R. Cudeck, S. du Toit, & D.
#'   Sörbom (Eds.), Structural equation modeling: Present and future - A
#'   Festschrift in honor of Karl Jöreskog (pp. 195-216). Lincolnwood, IL:
#'   Scientific Software International.
#' @source Bonifay, W. E., Reise, S. P., Scheines, R., & Meijer, R. R. (2015). When
#'   are multidimensional data unidimensional enough for structural equation
#'   modeling? An evaluation of the DETECT multidimensionality index. Structural
#'   Equation Modeling, 22, 504-516.
#' @source Reise, S. P., Scheines, R., Widaman, K. F., & Haviland, M. G. (2013).
#'   Multidimensionality and structural coefficient bias in structural equation
#'   modeling: A bifactor perspective. Educational and Psychological Measurement,
#'   73, 5-26.
#' @source Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016a). Applying bifactor
#'   statistical indices in the evaluation of psychological measures. Journal of
#'   Personality Assessment, 98, 223-237.
#' @source Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016b). Evaluating
#'   bifactor models: Calculating and interpreting statistical indices.
#'   Psychological Methods, 21, 137-150.
#'
#' @family reliability coefficients
#' @seealso [efa_fit()] for the solution these are computed from, and [OMEGA()], the
#'   superseded function that returns these same coefficients in a wide, per-factor layout.
#'
#' @export
#'
#' @examples
#' ## From an oblique EFA (correlated-factors) solution. With no factor_map, each
#' ## item is auto-assigned to its highest-loading factor.
#' efa_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                    estimator = "PAF", rotation = "promax")
#' efa_reliability(efa_mod)
#'
#' ## From a Schmid-Leiman solution, with an explicit indicator-to-factor map.
#' sl_mod <- efa_schmid_leiman(efa_mod, estimator = "PAF")
#' fc <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2
#' efa_reliability(sl_mod, factor_map = fc)
#'
#' ## Request a subset of the coefficients only.
#' efa_reliability(sl_mod, factor_map = fc,
#'                 coefficients = c("omega_total", "alpha"))
#'
#' ## From an oblique pattern matrix and its factor intercorrelations. This is
#' ## the same correlated-factors solution, and gives the same coefficients.
#' efa_reliability(efa_mod$rot_loadings, Phi = efa_mod$Phi,
#'                 cormat = test_models$baseline$cormat)
#'
#' ## From lavaan fits: a bifactor solution, and a correlated-factors one.
#' \donttest{
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#' mod_cf <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
#'            F2 =~ V7 + V8 + V9 + V10 + V11 + V12
#'            F3 =~ V13 + V14 + V15 + V16 + V17 + V18'
#' mod <- paste(mod_cf, 'g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
#'                            V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18',
#'              sep = "\n")
#' fit <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
#'                    sample.nobs = 500, estimator = "ml", orthogonal = TRUE)
#' print(efa_reliability(fit, g_name = "g"))
#'
#' ## No general factor: omega hierarchical, ECV, and PUC are omitted.
#' fit_cf <- lavaan::cfa(mod_cf, sample.cov = test_models$baseline$cormat,
#'                       sample.nobs = 500, estimator = "ml")
#' efa_reliability(fit_cf)
#' }
#' }
#'
efa_reliability <- function(model = NULL,
                            coefficients = NULL,
                            g_name = "g",
                            group_names = NULL,
                            factor_map = NULL,
                            variance = c("correlation", "sums_load"),
                            var_names = NULL,
                            fac_names = NULL,
                            g_load = NULL,
                            s_load = NULL,
                            u2 = NULL,
                            cormat = NULL,
                            pattern = NULL,
                            Phi = NULL) {

  variance <- .match_arg_ci(variance)
  checkmate::assert_string(g_name)
  checkmate::assert_character(group_names, null.ok = TRUE)

  # Validate the coefficient selector against the surfaced menu.
  menu <- .reliability_registry()$coefficient
  if (!is.null(coefficients)) {
    checkmate::assert_character(coefficients, any.missing = FALSE, min.len = 1L)
    bad <- setdiff(coefficients, menu)
    if (length(bad) > 0) {
      cli::cli_abort(
        c("Unknown {.arg coefficients}: {.val {bad}}.",
          "i" = "Choose from {.val {menu}}."),
        class = "efa_reliability_bad_coefficients"
      )
    }
  }

  # Read by both argument reports below, the second of which covers a `NULL` `model` too.
  lavaan_fit <- inherits(model, "lavaan")

  # These arguments describe manually supplied components; a solution passed in `model`
  # carries its own and the adapters never read them. Report them rather than drop them: a
  # user who expects a `Phi` to define the group factors' correlations, or a `u2` to replace
  # the solution's own, would otherwise get an orthogonal or unchanged solution's
  # coefficients with nothing to say why.
  if (!is.null(model)) {
    # Two of them are read by a matrix `model` and so are not reported for one: `u2`, as an
    # override of the uniquenesses the adapter would derive, and `Phi`, which a matrix never
    # merely ignores -- it makes the columns of an oblique pattern matrix correlated factors,
    # and with a matrix read as a hierarchy it is refused rather than dropped, where the
    # dispatch can say which of the two readings was taken. A `lavaan` fit reads neither the
    # correspondences nor a correlation matrix -- it is scored per group from its own -- so
    # those two are discarded there as well, and pointing at them would send the user
    # nowhere.
    unused <- c("g_load", "s_load", "var_names", "pattern")
    if (!is.matrix(model)) unused <- c(unused, "Phi", "u2")
    if (lavaan_fit) unused <- c(unused, "cormat", "factor_map")
    given <- unused[!vapply(mget(unused), is.null, logical(1))]
    if (length(given) > 0) {
      cli::cli_warn(
        c("{.arg {given}} {?is/are} ignored for a fitted solution.",
          "i" = "{cli::qty(given)}{?It/They} describe{?s/} a solution supplied through the components; {.arg model} carries its own.",
          "i" = if (lavaan_fit) {
            "A {.cls lavaan} fit is scored per group from its own loadings, uniquenesses, and model-implied composite variances."
          } else {
            "Supply {.arg cormat} to score the solution against an observed correlation matrix, and {.arg factor_map} to set the indicator-to-factor correspondences."
          }),
        class = "efa_reliability_args_ignored"
      )
    }
  }

  # `fac_names` and `group_names` label the result rather than describe a solution, so the
  # reason above -- that `model` carries its own components -- is not theirs: each is read on
  # one kind of input only. `fac_names` names the factors of every solution but a `lavaan`
  # fit, whose factor labels come from the model syntax; `group_names` names the blocks of a
  # `lavaan` multiple-group fit, the only input scored as more than one solution. Reported for
  # the same reason the components are: a caller who labelled their factors or their groups
  # would otherwise get the input's own labels back with nothing to say why -- and a
  # `fac_names` of the wrong length is not even refused for a `lavaan` fit, the length check
  # below belonging to the routes that read it. Both are checked outside the block above,
  # which a `NULL` `model` does not enter although it reads no `group_names` either.
  if (lavaan_fit && !is.null(fac_names)) {
    cli::cli_warn(
      c("{.arg fac_names} is not used for a {.cls lavaan} fit.",
        "i" = "Its factor labels come from the model syntax; rename the latent variables there to change them."),
      class = "efa_reliability_arg_not_applicable"
    )
  }
  if (!is.null(group_names)) {
    # A single-group `lavaan` fit reads the labels no more than a non-`lavaan` input does:
    # .rel_adapt_lavaan() checks their length and the result builder then drops them, one
    # block needing no key, so its `group` column stays NA like every ungrouped solution's.
    # Keyed on the group count rather than on the class for exactly that reason. lavaan is
    # asked for the count only when there are labels to place.
    if (lavaan_fit) .require_lavaan()
    multigroup <- lavaan_fit && lavaan::lavInspect(model, "ngroups") > 1L
    if (!multigroup) {
      cli::cli_warn(
        c("{.arg group_names} is only used for a {.cls lavaan} multiple-group fit.",
          "i" = "Every other input is scored as a single ungrouped solution, whose {.field group} column is {.code NA}."),
        class = "efa_reliability_arg_not_applicable"
      )
    }
  }

  # Dispatch to the right adapter, then score each factor solution through the
  # shared reliability core (add_rel = TRUE so standardized alpha is available;
  # the registry drops the core's CR/AVE columns). lavaan input is scored per
  # group with the model-implied composite variance.
  no_general <- FALSE

  if (inherits(model, "lavaan")) {

    .require_lavaan()
    adapt <- .rel_adapt_lavaan(model, g_name = g_name, group_names = group_names)
    used_variance <- adapt$variance
    no_general <- isTRUE(adapt$correlated)

    # Distinct group labels: the group column keys each block of the tidy result,
    # so duplicate labels would silently merge groups.
    if (length(adapt$groups) > 1L && anyDuplicated(adapt$group_names)) {
      cli::cli_abort(
        c("{.arg group_names} must be unique.",
          "i" = "Duplicate labels would merge groups in the result."),
        class = "efa_reliability_group_names"
      )
    }

    if (isTRUE(adapt$higher_order)) {
      cli::cli_inform(
        c("i" = "The specified general factor is a second-order factor; coefficients are computed on the Schmid-Leiman transformed second-order solution."),
        class = "efa_reliability_g_second_order"
      )
    }

    if (isTRUE(adapt$few_loadings)) {
      cli::cli_inform(
        c("i" = "Some variables have fewer than two loadings; did you enter a bifactor model? Provide a bifactor model, a second-order model, or a single-factor model."),
        class = "efa_reliability_few_loadings"
      )
    }

    # No variable loads on more than one factor, so the fit carries no general factor over
    # and above its factors. Said rather than left to be read off a missing column: a user
    # who fitted what they think of as a bifactor model, or who named `g_name`, would
    # otherwise only see a table without the general-factor coefficients.
    if (no_general) {
      cli::cli_inform(
        c("i" = "Each variable loads on one factor only, so the fit is scored as a correlated-factors solution; omega hierarchical, ECV, and PUC are not defined for it and are omitted."),
        class = "efa_reliability_no_general_factor"
      )
    }

    mats <- vector("list", length(adapt$groups))
    informed_single <- FALSE
    for (i in seq_along(adapt$groups)) {
      grp <- adapt$groups[[i]]
      if (isTRUE(grp$single)) {
        if (!informed_single) {
          .rel_inform_single_factor()
          informed_single <- TRUE
        }
        # A single factor is scored through the core on the spec the adapter normalized it
        # to, as every other group is; the coefficients it does not define are then dropped,
        # exactly as on the other input routes.
        mats[[i]] <- .rel_drop_single_factor(
          .reliability_core(grp, used_variance, add_ind = TRUE, add_rel = TRUE,
                            arg = "factor_map"),
          fac_names = grp$fac_label)
      } else {
        mats[[i]] <- .reliability_core(grp, used_variance, add_ind = TRUE,
                                       add_rel = TRUE, arg = "factor_map")
        # A correlated-factors fit defines only some of the coefficients the core
        # computes; the rest are dropped, per group, as on the non-lavaan routes.
        if (no_general) mats[[i]] <- .rel_drop_general(mats[[i]])
      }
    }

    # One matrix (single ungrouped fit) or a named list (multiple groups).
    x <- if (length(mats) == 1L) mats[[1]] else
      stats::setNames(mats, adapt$group_names)

  } else {

    spec <- .rel_dispatch_spec(model, factor_map = factor_map,
                               var_names = var_names, fac_names = fac_names,
                               g_load = g_load, s_load = s_load, u2 = u2,
                               cormat = cormat, pattern = pattern, Phi = Phi)

    # A solution with exactly one factor is scored as the single factor it is, whichever slot
    # that factor arrived in -- as the general factor of a one-column bifactor matrix, or as
    # the only group factor of a one-factor `efa_fit()` solution or of manual components with
    # a zero `g_load`. Normalizing them to one spec is what makes every route report the same
    # coefficients for the same solution (see .rel_single_factor_spec).
    single_spec <- .rel_single_factor_spec(spec)
    single <- !is.null(single_spec)

    # One label per factor, checked once the solution is in hand and its factors counted.
    # `fac_names` labels the group-factor rows, of which a single-factor solution has one --
    # its own factor, whichever slot it arrived in -- and any other solution has one per
    # group loading column; the general-factor row is labelled "g" and takes none. Nothing
    # downstream can say so: too few or too many reach `rownames()` on the coefficient
    # matrix and fail there in base R's terms, naming neither the argument nor the count it
    # needed, and a single-factor solution does not fail at all -- .rel_single_factor_label()
    # reads a vector that is not of length one as no name and falls back to "F1", so the
    # labels the caller asked for are dropped without a word.
    n_fac <- if (single) 1L else ncol(as.matrix(spec$s_load))
    if (!is.null(fac_names) && length(fac_names) != n_fac) {
      cli::cli_abort(
        c("{.arg fac_names} must have one name per factor.",
          "x" = "It has {length(fac_names)} name{?s}, but the solution has {n_fac} factor{?s}.",
          "i" = if (single) {
            "A single-factor solution has no group factors; {.arg fac_names} labels its one factor."
          } else {
            "Name the group factors in the column order of the loadings; the general factor is labelled {.val g}."
          }),
        class = "efa_reliability_fac_names_length"
      )
    }

    single_fac_names <- NULL
    if (single) {
      # Read off the spec as it arrived, before it is replaced: the normalized one has no
      # group factor left to carry the name the input gave this one.
      single_fac_names <- spec$fac_names
      spec <- single_spec
    }

    # An all-zero general factor marks a solution with no general factor -- an oblique EFA,
    # or a bifactor/manual input whose general column is zero. It says nothing about
    # whether the group factors are correlated: only a spec that carries the factor
    # intercorrelations describes correlated ones, and the coefficients account for them
    # in either variance convention. A single factor is not such a solution: it is itself
    # the general factor, and the normalization above has put it there.
    no_general <- !single && isTRUE(all(spec$g_load == 0))

    # Only correlation mode needs a correlation matrix, whether or not the solution has a
    # general factor. The model-implied "sums_load" variance needs none, and is correct for
    # every solution the adapters produce -- it reads the factor intercorrelations where the
    # spec carries them. Abort with an actionable message when correlation mode is asked for
    # and no matrix is available (the adapters supply one when they can; an SL object from a
    # bare loading matrix, or manual components without cormat/pattern/Phi, leave it NULL)
    # rather than divide the omega denominators by sum(NULL) = 0 and return non-finite values.
    if (variance == "correlation" && is.null(spec$cormat)) {
      cli::cli_abort(
        c("A correlation matrix is required for {.code variance = \"correlation\"} but none is available.",
          "i" = "Supply {.arg cormat} (or {.arg pattern} and {.arg Phi}), or use {.code variance = \"sums_load\"}."),
        class = "efa_reliability_need_cormat"
      )
    }

    if (single) .rel_inform_single_factor()

    used_variance <- variance
    x <- .reliability_core(spec, variance, add_ind = TRUE, add_rel = TRUE,
                           arg = "factor_map")

    # A correlated-factors solution (no general factor) does not identify the
    # general-factor decomposition, and a single-factor one leaves the general-factor
    # decomposition with nothing to decompose, so each drops the coefficients it does not
    # define (see .rel_drop_general() and .rel_drop_single_factor() for which, and why).
    # NA-ing the cells lets the result builder omit them as any other undefined coefficient.
    if (single) {
      x <- .rel_drop_single_factor(x, fac_names = single_fac_names)
    } else if (no_general) {
      x <- .rel_drop_general(x)
    }

  }

  res <- .reliability_result(x, settings = list(variance = used_variance,
                                                no_general = no_general))

  if (!is.null(coefficients)) {
    # Note any requested coefficient this solution does not define (e.g. omega
    # hierarchical for a correlated-factors EFA, or omega subscale for a single-factor
    # solution), so a reduced or empty selection is not silent.
    undefined <- setdiff(coefficients, unique(res$coefficient))
    if (length(undefined) > 0) {
      cli::cli_inform(
        c("i" = "{cli::qty(undefined)}The coefficient{?s} {.val {undefined}} {?is/are} not defined for this solution and {?is/are} omitted."),
        class = "efa_reliability_coef_undefined"
      )
    }
    res <- .reliability_select(res, coefficients)
  }

  res

}


# Read a loading matrix and its factor intercorrelations as one correlated-factors
# solution, returning the normalized reliability spec. Shared by the two matrix routes
# that reach this reading, so a pattern matrix gives the same coefficients whether or not
# it still carries the loading class.
.rel_read_pattern <- function(model, Phi, u2, factor_map, cormat, fac_names) {
  # The loadings and any supplied uniquenesses get the checks the bifactor-matrix route
  # applies, and for the same reasons.
  .assert_args({
    checkmate::assert_numeric(model, any.missing = FALSE, finite = TRUE)
    if (!is.null(u2)) {
      checkmate::assert_numeric(u2, any.missing = FALSE, finite = TRUE,
                                len = nrow(model))
    }
  })
  L <- as.matrix(unclass(model))
  # Checked here because the uniquenesses are derived from it just below, which needs
  # it conformable with the loadings; .rel_adapt_manual() checks it again on its own.
  .rel_check_phi(Phi, ncol(L))
  # 1 - diag(L Phi L'), the communalities the oblique pattern implies under these
  # factor correlations, unless the caller supplied uniquenesses of their own. The
  # same expression the oblique `efa_fit()` route uses, so both reach the same ones.
  if (is.null(u2)) u2 <- 1 - diag(L %*% Phi %*% t(L))
  var_names <- rownames(L)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(nrow(L)))
  if (is.null(fac_names)) {
    fac_names <- colnames(L)
    if (is.null(fac_names)) fac_names <- paste0("F", seq_len(ncol(L)))
  }
  # `pattern` is reported as ignored for any `model`, so it is not passed on: here
  # `Phi` describes the loading matrix that was given, not a separate solution.
  .rel_adapt_manual(g_load = rep(0, nrow(L)), s_load = L, u2 = u2,
                    var_names = var_names, factor_corres = factor_map,
                    type = "psych", cormat = cormat, Phi = Phi,
                    fac_names = fac_names)
}


# Dispatch a non-lavaan input to its front-end adapter, returning the normalized
# reliability spec. Mirrors OMEGA()'s model dispatch, adding the oblique EFA and
# raw bifactor-matrix paths that OMEGA does not expose. The internal adapters take a
# `type` that gates how the item-to-factor map is derived when `factor_map` is
# absent; efa_reliability always auto-maps to the highest-loading factor (adapter
# type "psych"), so a supplied map is used and an omitted one is derived.
.rel_dispatch_spec <- function(model, factor_map, var_names, fac_names,
                               g_load, s_load, u2, cormat, pattern, Phi) {

  if (inherits(model, "SL")) {
    .rel_adapt_SL(model, factor_corres = factor_map, type = "psych",
                  cormat = cormat, fac_names = fac_names)
  } else if (inherits(model, "schmid")) {
    .rel_adapt_schmid(model, factor_corres = factor_map, type = "psych",
                      cormat = cormat, fac_names = fac_names)
  } else if (inherits(model, "EFA")) {
    .rel_adapt_efa(model, factor_corres = factor_map, type = "psych",
                   cormat = cormat, fac_names = fac_names)
  } else if (is.matrix(model)) {
    if (!is.numeric(model)) {
      cli::cli_abort(
        c("A matrix {.arg model} must be a numeric bifactor loading matrix (general factor in the first column).",
          "i" = "Enter numeric loadings, or pass a fitted model object instead."),
        class = "efa_reliability_bad_matrix"
      )
    }
    # Recorded before the Schmid-Leiman branch below replaces `model` with its factor
    # columns, which drops the class that says so. A Schmid-Leiman table states that it is
    # a hierarchy; a matrix that carries no loading class merely fails to state anything,
    # which is not the same evidence, and the two are told apart where `Phi` arrives.
    declared_hierarchy <- inherits(model, c("efa_sl_loadings", "SLLOADINGS"))
    if (declared_hierarchy) {
      # A Schmid-Leiman loading table is a bifactor matrix plus the communality and
      # uniqueness columns it prints. Score it as such, dropping those two -- they would
      # otherwise be read as two more group factors -- and taking the uniquenesses from
      # the table rather than recomputing them from the loadings.
      sl <- unclass(model)
      nms <- colnames(sl)
      # The table is always built as the general factor, then the group factors, then h2
      # and u2, so an unlabelled one is read off that layout rather than left with its
      # two display columns scored as group factors.
      if (is.null(nms) && ncol(sl) > 3L) {
        display <- seq_len(ncol(sl)) > ncol(sl) - 2L
        u2_col <- ncol(sl)
      } else {
        # A table too short to carry that layout keeps every column, as before.
        if (is.null(nms)) nms <- character(ncol(sl))
        display <- nms %in% c("h2", "u2")
        u2_col <- match("u2", nms)
      }
      if (is.null(u2) && !is.na(u2_col)) u2 <- sl[, u2_col]
      model <- sl[, !display, drop = FALSE]
      # An estimated Schmid-Leiman table carries no structural zeros, so the raw-matrix
      # default -- every non-zero group loading -- would assign every variable to every
      # group factor. Derive the correspondences the same way the Schmid-Leiman object
      # does instead, so both routes to one solution report the same coefficients.
      if (is.null(factor_map)) {
        factor_map <- .rel_map(model[, -1, drop = FALSE], type = "psych")
      }
    } else if (inherits(model, c("efa_loadings", "LOADINGS"))) {
      # An ordinary loading table is a pattern or structure matrix: its first column is
      # simply the first factor, not a general factor. Read as a bifactor matrix it
      # yields plausible but meaningless coefficients.
      #
      # With `Phi` it is a complete correlated-factors solution -- the oblique pattern and
      # the correlations of the factors it belongs to -- and is scored as one, through the
      # same components the manual route takes for such a solution. Without it nothing says
      # the columns are correlated factors rather than a hierarchy, so a matrix of two or
      # more factors is refused and the hierarchy asked for instead. The class states that
      # the matrix is such a solution's pattern, so reading it as one overrides nothing and
      # is not reported, where the same reading of a matrix that states nothing is.
      if (!is.null(Phi)) {
        return(.rel_read_pattern(model, Phi = Phi, u2 = u2, factor_map = factor_map,
                                 cormat = cormat, fac_names = fac_names))
      }
      # A single column is the one case in which there is nothing to confuse: one factor is
      # no hierarchy, and a solution with one factor has no factor intercorrelations to
      # supply, so neither remedy the refusal offers is one its caller could take. It falls
      # through to the bare-matrix route below, which reads a one-column matrix as a general
      # factor with no group factors; the front end then rewrites that spec into the
      # single-factor one (.rel_single_factor_spec), as it does for every other route that
      # can carry one factor, so the classed matrix, the unclassed one, and the one-factor
      # `efa_fit()` object it came from all reach the same coefficients.
      if (ncol(model) > 1L) {
        cli::cli_abort(
          c("{.arg model} is a loading matrix of an ordinary factor solution, not a bifactor loading matrix.",
            "i" = "Pass the {.fn efa_fit} object itself for a correlated-factors solution, or, if this is its pattern matrix, supply the matching factor intercorrelations in {.arg Phi} to score it as one.",
            "i" = "To read it as a hierarchy instead, transform it with {.fn efa_schmid_leiman} first.",
            "i" = "If it really is a bifactor solution, pass {.code unclass(x)} to read its first column as the general factor."),
          class = "efa_reliability_pattern_matrix"
        )
      }
    } else if (!anyNA(model) && .is_cormat(model)) {
      # `.is_cormat()` aborts on a missing value in a matrix that is otherwise shaped like a
      # correlation matrix, and does so in terms of its own argument, hence the `anyNA()`
      # guard: an incomplete matrix is named by the loading check below rather than by a
      # misdirected error.
      #
      # A correlation matrix is square, symmetric, and has a unit diagonal at once; a
      # bifactor loading matrix is essentially never all three, since a unit diagonal would
      # mean every variable loads exactly 1 on one factor. Read as loadings it takes the
      # first variable's correlations for general-factor loadings and returns plausible but
      # meaningless coefficients, so it is refused where the mistake can still be named.
      cli::cli_abort(
        c("{.arg model} is a correlation matrix, not a bifactor loading matrix.",
          "i" = "Fit a solution first, e.g. with {.fn efa_fit} or {.fn efa_schmid_leiman}, and pass that.",
          "i" = "To score a solution against these correlations, pass the matrix as {.arg cormat}."),
        class = "efa_reliability_cormat_as_model"
      )
    }
    # Everything still here is read as a hierarchy unless `Phi` says otherwise: a bare
    # bifactor loading matrix, or the loading table of a Schmid-Leiman solution reduced to
    # its factor columns above. Read as one, its general factor is the first column and its
    # group factors are orthogonal by construction, so there are no correlated group factors
    # for `Phi` to be the correlation matrix of, and the general and group parts would no
    # longer partition the composite.
    # A Schmid-Leiman table states that it is a hierarchy, so the combination contradicts
    # itself and is refused, as it is through the components. A matrix that carries no
    # loading class states nothing: the class is not a reliable signal here, since `[`,
    # `rbind()` and a round-trip through a data frame all drop it while the argument stays,
    # so a caller who subsets an oblique pattern and keeps its `Phi` holds a complete
    # correlated-factors solution in two objects. `Phi` settles the reading for that matrix
    # -- the group factors of a hierarchy are orthogonal and have no factor correlations to
    # supply -- so it is read as the pattern it can only be. The reading is reported rather
    # than taken quietly, because the two give different coefficients: a hierarchy also
    # defines omega hierarchical, the ECV and the PUC, which a correlated-factors solution
    # does not, and a caller who meant the bifactor reading would otherwise see them absent
    # with nothing to say why.
    #
    # One column is the exception this route makes throughout: one factor is not a hierarchy,
    # and the only correlation matrix its factor can have is the 1 by 1 identity, which says
    # nothing either way. It is checked and then changes nothing, as on the classed route.
    if (!is.null(Phi)) {
      # Before the conformability check, so a Schmid-Leiman table is told that its own
      # reading refuses `Phi` rather than that this `Phi` is the wrong shape: no `Phi` of
      # any shape belongs with one, and naming the shape would send the caller to fix it.
      if (ncol(model) > 1L && declared_hierarchy) {
        cli::cli_abort(
          c("{.arg Phi} describes correlated group factors, which the loading table of a Schmid-Leiman solution does not have.",
            "i" = "To score the pattern matrix of an oblique solution, pass the {.fn efa_fit} object, or supply the components: the pattern in {.arg s_load}, a {.arg g_load} that is zero throughout, and {.arg u2} and {.arg var_names}, with this {.arg Phi}.",
            "i" = "To reconstruct a correlation matrix for a hierarchy instead, supply its components together with the {.arg pattern} of the separate oblique solution {.arg Phi} belongs to.",
            "i" = "Otherwise drop {.arg Phi}: the group factors of a hierarchy are uncorrelated."),
          class = "efa_reliability_phi_with_general"
        )
      }
      if (ncol(model) > 1L) {
        # Read first, report after. Everything the reading can refuse -- a `Phi` that
        # cannot be these columns' correlation matrix, a missing or infinite loading or
        # uniqueness -- is named on its own that way, rather than under a warning about a
        # reading that is then abandoned. The spec is built before the warning and returned
        # after it, so the two cannot come apart.
        spec <- .rel_read_pattern(model, Phi = Phi, u2 = u2, factor_map = factor_map,
                                  cormat = cormat, fac_names = fac_names)
        cli::cli_warn(
          c("{.arg model} is read as the pattern matrix of a correlated-factors solution, because {.arg Phi} gives its factor correlations.",
            "i" = "Every column is a group factor, so the coefficients a general factor defines -- omega hierarchical, the ECV, and the PUC -- are not reported.",
            "i" = "Drop {.arg Phi} to read the first column as a general factor and score the matrix as a bifactor one."),
          class = "efa_reliability_phi_as_pattern"
        )
        return(spec)
      }
      # One column is not a hierarchy and has no correlated group factors, so it falls
      # through to the bifactor route below; the check is all that applies to it.
      .rel_check_phi(Phi, ncol(model))
    }
    # Uniquenesses reach this route either from the caller, overriding the ones the adapter
    # would derive from the loadings, or off a Schmid-Leiman table's own `u2` column. Both
    # get the check the manual route's do, and for the same reason: a missing or infinite
    # one is scored rather than refused, and the table it produces looks ordinary -- an
    # infinite uniqueness drives its composite's omega total to a flat 0, a missing one
    # drops the coefficients of every composite the variable belongs to -- while nothing
    # downstream can name it.
    # The loadings get the same check, and not only because a missing or infinite one is no
    # more a loading than a missing uniqueness is a uniqueness: when `u2` is left out the
    # adapter derives it as 1 - rowSums(loadings^2), so an unchecked loading carries a
    # missing or infinite uniqueness straight past the check above and into the coefficients.
    .assert_args({
      checkmate::assert_numeric(model, any.missing = FALSE, finite = TRUE)
      if (!is.null(u2)) {
        checkmate::assert_numeric(u2, any.missing = FALSE, finite = TRUE,
                                  len = nrow(model))
      }
    })

    .rel_adapt_bifactor(model, factor_corres = factor_map, u2 = u2,
                        cormat = cormat, fac_names = fac_names)
  } else if (is.null(model)) {
    if (is.null(var_names) || is.null(g_load) || is.null(s_load) || is.null(u2)) {
      cli::cli_abort(
        "Specify all of {.arg var_names}, {.arg g_load}, {.arg s_load}, and {.arg u2}.",
        class = "efa_reliability_missing_args"
      )
    }
    # Validate the manually supplied components, as OMEGA() does: a non-matrix
    # s_load would otherwise be coerced (e.g. a vector to a one-column matrix) and
    # silently yield wrong coefficients instead of an error.
    if (!inherits(s_load, c("matrix", "SLLOADINGS"))) {
      cli::cli_abort(
        c("Invalid {.arg s_load}.",
          "i" = "Supply a Schmid-Leiman or bifactor group-loading matrix of class {.cls matrix} or {.cls efa_sl_loadings}."),
        class = "efa_reliability_bad_s_load"
      )
    }
    # A loading or a uniqueness that is missing or infinite is not one, and nothing
    # downstream can say so: every check that might is blind to exactly these values. The
    # Heywood report tests `u2 < eps`, which skips NA; the keying check and the
    # standardization check both compare totals that are themselves non-finite and so
    # return early. What the caller gets instead is not an empty table but a plausible one
    # -- an infinite uniqueness drives its composite's omega total to a flat 0, and an
    # infinite general loading leaves an ordinary-looking coefficient -- so these have to
    # be refused where they enter. Checking the group loadings as a numeric vector also
    # catches a matrix that is not numeric at all, which the class test above admits.
    .assert_args({
      checkmate::assert_numeric(s_load, any.missing = FALSE, finite = TRUE)
      checkmate::assert_numeric(g_load, any.missing = FALSE, finite = TRUE)
      checkmate::assert_numeric(u2, any.missing = FALSE, finite = TRUE)
    })
    # Square, symmetric, and a unit diagonal at once identifies a correlation matrix: a
    # group-loading matrix is essentially never all three, since a unit diagonal would mean
    # every variable loads exactly 1 on one group factor. Nothing else here would say so --
    # the class test above admits it, the numeric checks pass, and its p rows match p
    # variables, so it is scored as loadings and yields plausible but meaningless
    # coefficients. Named where the mistake is made, as for a `model` matrix. After the
    # numeric checks, so a `s_load` carrying missing values is reported as such rather than
    # by `.is_cormat()`, which aborts on one in terms of its own argument.
    if (.is_cormat(s_load)) {
      cli::cli_abort(
        c("{.arg s_load} is a correlation matrix, not a group-factor loading matrix.",
          "i" = "Supply the group-factor loadings of a Schmid-Leiman or bifactor solution.",
          "i" = "To score the components against these correlations, pass the matrix as {.arg cormat}."),
        class = "efa_reliability_cormat_as_s_load"
      )
    }
    # The whole printed Schmid-Leiman table, rather than its group-factor columns: its `h2`
    # and `u2` columns would be scored as two more group factors, and its general column as
    # a third, on top of the `g_load` given separately. The matrix route reads that table by
    # its layout because it is given as the entire solution; here the components are named
    # one by one, so the columns to pass are the ones the caller means.
    display_cols <- intersect(c("h2", "u2"), colnames(s_load))
    if (length(display_cols) > 0) {
      cli::cli_abort(
        c("{.arg s_load} carries the {.val {display_cols}} display column{?s} of a Schmid-Leiman loading table, which {?is/are} not group factors.",
          "i" = "Pass the group-factor columns only, e.g. {.code s_load = x$sl[, c(\"F1\", \"F2\")]}.",
          "i" = "To score the whole table as one solution, pass it as {.arg model} instead."),
        class = "efa_reliability_s_load_display_cols"
      )
    }
    # The components must agree in length so they combine into one item-by-factor spec;
    # otherwise they only fail deep in the algebra with an opaque base-R error.
    n_items <- nrow(as.matrix(s_load))
    if (length(g_load) != n_items || length(u2) != n_items ||
        length(var_names) != n_items) {
      cli::cli_abort(
        c("{.arg g_load}, {.arg u2}, and {.arg var_names} must each have one entry per variable.",
          "i" = "{.arg s_load} has {n_items} row{?s}, but {.arg g_load} has {length(g_load)}, {.arg u2} has {length(u2)}, and {.arg var_names} has {length(var_names)}."),
        class = "efa_reliability_length_mismatch"
      )
    }
    .rel_adapt_manual(g_load = g_load, s_load = s_load, u2 = u2,
                      var_names = var_names, factor_corres = factor_map,
                      type = "psych", cormat = cormat, pattern = pattern, Phi = Phi,
                      fac_names = fac_names)
  } else {
    cli::cli_abort(
      c("Invalid {.arg model}.",
        "i" = "Enter a {.cls lavaan}, {.fn psych::schmid}, {.cls SL}, or {.cls EFA} object, a bifactor loading matrix, or specify {.arg var_names}, {.arg g_load}, {.arg s_load}, and {.arg u2}.",
        "i" = "From an {.cls efa_average} object, use the averaged loading matrix in
               {.code $loadings$average}, with {.code $Phi$average} as {.arg Phi} and
               {.code $orig_R} as {.arg cormat}."),
      class = "efa_reliability_invalid_model"
    )
  }

}


# Keep only the requested coefficients in an efa_reliability result, preserving the
# object's schema, settings, and kind tagging (base data.frame subsetting drops the
# custom attributes).
.reliability_select <- function(res, coefficients) {
  settings <- attr(res, "settings")
  kind <- attr(res, "kind")
  out <- res[res$coefficient %in% coefficients, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "settings") <- settings
  attr(out, "kind") <- kind[names(kind) %in% coefficients]
  class(out) <- c("efa_reliability", "data.frame")
  out
}
