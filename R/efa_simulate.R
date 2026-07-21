#' Simulate data from a common-factor population model
#'
#' Draws data from a population correlation matrix, given either directly or
#' built from a factor model. The population correlation is either supplied in
#' `R`, or assembled from a loading matrix `Lambda`, the factor intercorrelations
#' `Phi`, and the unique variances `Psi` as
#' \eqn{R = Lambda\, Phi\, Lambda' + Psi}, standardized to a correlation matrix.
#' By default (`marginals = "normal"`) cases are drawn with normal marginals via
#' a matrix square root of the population correlation (a Cholesky factor, or a
#' symmetric eigen square root for a positive-semidefinite but singular
#' population). With `marginals = "empirical"` the cases instead reproduce the
#' population correlation while carrying the empirical marginal distributions of
#' a supplied data set. With `marginals = "VM"` or `"IG"` the cases reproduce the
#' population correlation while carrying non-normal marginals with a prescribed
#' `skewness` and `kurtosis`. Setting `categories` additionally discretizes the
#' drawn data into ordered categories, optionally so that the population polychoric
#' correlation of the categorized data equals the target correlation. Setting
#' `target_rmsea` or `target_cfi` perturbs the population with model error, so the
#' factor model fits it only approximately (a more realistic simulation target).
#'
#' @details
#' Provide the population either as a ready correlation matrix in `R`, or through
#' the model components `Lambda`, `Phi`, and `Psi`; the two ways are mutually
#' exclusive. When the model components are used, `Phi` defaults to the identity
#' matrix (orthogonal factors) and `Psi` defaults to the unique variances that
#' make the population a correlation matrix (\eqn{1 - \mathrm{diag}(Lambda\, Phi\,
#' Lambda')}); the assembled covariance is standardized with
#' [stats::cov2cor()] so a non-standardized `Psi` still yields a correlation
#' matrix. With the default `Psi`, a factor model whose implied communalities
#' exceed 1 (a Heywood case) leaves no unique variance and is rejected; a `Psi`
#' you supply is instead only required to give positive variances and a
#' positive-semidefinite population.
#'
#' With `marginals = "empirical"`, the iterative rank-matching algorithm of Ruscio
#' and Kaczetow (2008) reproduces the population correlation while each variable
#' takes the empirical marginal distribution of the matching column of
#' `marginal_data` (resampled with replacement). Only the marginals of
#' `marginal_data` are used; its own correlations are ignored, and the drawn columns
#' follow the population's variables, not those of `marginal_data`.
#'
#' With `marginals = "VM"` (Vale-Maurelli, 1983) or `"IG"` (the independent-generator
#' method; Foldnes & Olsson, 2016), the cases reproduce the population correlation
#' while carrying non-normal marginals with the target `skewness` and (excess)
#' `kurtosis`. The Vale-Maurelli family does not span every valid non-normal
#' distribution (Foldnes & Grønneberg, 2015); `"IG"` covers distributions `"VM"`
#' cannot. Both accept `skewness` and `kurtosis` as a single value (used for every
#' variable) or one value per variable, defaulting the unset one to 0. Not every
#' (`skewness`, `kurtosis`) pair is attainable: every distribution needs excess
#' kurtosis of at least `skewness^2 - 2`, and the method covers a smaller region
#' still, so an unreachable request is rejected. For `marginals = "VM"`, the
#' intermediate correlation matrix used for the draw can itself be
#' non-positive-definite; it is rejected unless `force_pd = TRUE`, which projects it
#' to the nearest correlation matrix (via [psych::cor.smooth()]) with a warning.
#'
#' With `categories`, the drawn data are discretized into ordered categories (an
#' integer code `1` to `K`). `categories` gives either the number of equally probable
#' categories (one count for every variable, or one per variable) or, as a list of
#' proportion vectors, the marginal category proportions per variable. The cut points are
#' the thresholds that reproduce the requested proportions (Olsson, 1979): the
#' standard-normal quantiles for `marginals = "normal"`, and for `marginals = "VM"` those
#' quantiles mapped through the same Fleishman cubic the draw uses, so the requested
#' proportions are reproduced on the non-normal scale too. Under `marginals = "IG"` the
#' thresholds stay on the standard-normal scale while the data do not, so the achieved
#' proportions depart from the request systematically rather than by sampling noise; the
#' departure grows with the non-normality, and only the *number* of categories is
#' guaranteed. The same holds for a `"VM"` variable whose Fleishman cubic is not increasing
#' over its own thresholds and the tails beyond them, which keeps the normal-scale ones and
#' is reported with a warning; this arises when a turning point of the cubic sits at or
#' near an outer threshold -- under a strongly platykurtic marginal, or under substantial
#' skewness or kurtosis combined with a small outer-category proportion.
#' Because categorization attenuates product-moment correlations, the
#' categorized data's Pearson correlation is smaller in magnitude than the population
#' correlation; under non-normal marginals its polychoric correlation departs from the
#' population as well. `match` changes none of this -- it asserts an intent rather than
#' selecting a computation, as described under that argument. Ordinal output is not
#' available with `marginals = "empirical"`. Empty categories left by a draw are reported
#' with a warning, as they destabilize the polychoric correlation and the factor analysis.
#'
#' With `missing`, missing values are introduced into the drawn data under a chosen
#' mechanism (Rubin, 1976), each variable holed at a target expected rate
#' `missing_prop`. `"MCAR"` draws an independent mask, so missingness is unrelated to
#' the data. `"MAR"` and `"MNAR"` set each case's missing probability by a logistic
#' model of a standardized predictor: another variable for `"MAR"` (chosen by
#' `missing_predictor`) or the variable's own value for `"MNAR"`, with slope
#' `missing_strength`. The mechanism acts on the drawn (latent) values, so when
#' `categories` also discretizes the data the missingness is keyed on the underlying
#' value, not the category code. For `"MAR"` the predictor is evaluated on the complete
#' drawn values, but every variable is holed at rate `missing_prop`, so a variable's MAR
#' predictor is itself missing for roughly a `missing_prop` fraction of the cases whose
#' missingness it drove. The mechanism is therefore MAR conditional on the *complete*
#' data, and **not** ignorable for an analyst who sees only the observed data: estimators
#' that are consistent under ignorable MAR, such as `cor_method = "fiml"` in [efa_fit()]
#' and the multiple imputation behind [efa_mi()], keep a residual bias here that grows
#' with `missing_prop` and `missing_strength`. The returned matrix carries the `NA`s,
#' which the correlation estimators handle downstream.
#'
#' With `model_error`, the population is perturbed away from the exact factor
#' structure so the `q`-factor model (`q = ncol(Lambda)`) fits it only approximately,
#' at a prescribed misfit; exact factor structures are unrealistic and overstate
#' recovery in simulation studies (MacCallum, 2003). The perturbation is applied once
#' to the population, and the achieved fit is computed with the same fit-index
#' formulas [efa_fit()] uses and returned in the `model_error` element. It is applied only
#' when a target is supplied (`target_rmsea` and/or `target_cfi`), needs a
#' factor-model population (`Lambda`) with residual degrees of freedom and an exact
#' factor structure (a diagonal `Psi`), and is orthogonal to the marginal, ordinal,
#' and missing-data options. Three methods are available. `"CB"` (Cudeck & Browne,
#' 1992) matches the target RMSEA to numerical precision and keeps the `q`-factor
#' model the exact minimizer (the CFI follows as a derived quantity). `"TKL"` (Tucker,
#' Koopman & Linn, 1969) adds minor common factors tuned so the achieved RMSEA -- and,
#' optionally, CFI -- match the target(s); with a single target the match is close,
#' with both it is a compromise. `"WB"` (Wu & Browne, 2015) draws the population from
#' an inverse-Wishart distribution around the model-implied correlation; its calibration
#' applies to the *best-fitting* model, so the reported misfit of the *generating* model
#' is systematically larger than the target -- by roughly \eqn{\sqrt{(p(p-1)/2)/df}},
#' about 1.4 times for 12 variables and 3 factors. Use `"CB"` when the reported RMSEA
#' must equal the target. `"CB"` and
#' `"WB"` target the RMSEA only; `"TKL"` can target the RMSEA and/or the CFI. The
#' reported RMSEA/CFI is the misfit of the specified generating model.
#'
#' Replicated draws (`n_datasets > 1`) are generated in parallel across
#' replicates with \pkg{future.apply}; a parallel plan can be selected with
#' [future::plan()] (the default plan runs sequentially). Each replicate is
#' assigned its own reproducible random-number stream, so with a fixed `seed` the
#' output is identical regardless of the number of workers.
#'
#' @param N numeric. Number of cases (rows) to draw per dataset. Required unless
#'   `return_pop = TRUE`.
#' @param Lambda matrix. A `p` by `m` matrix of factor loadings. Supply this
#'   (optionally with `Phi` and `Psi`) instead of `R` to build the population
#'   from a factor model.
#' @param Phi matrix. The `m` by `m` factor intercorrelation matrix. Only used
#'   with `Lambda`. Default is `NULL`, in which case the factors are orthogonal
#'   (an identity matrix).
#' @param Psi numeric vector or matrix. The unique variances: either a length-`p`
#'   vector or a `p` by `p` matrix (added as the residual covariance). Only used
#'   with `Lambda`. Default is `NULL`, in which case the unique variances that
#'   standardize the population to a correlation matrix are used.
#' @param R matrix. A `p` by `p` population correlation matrix to draw from
#'   directly. Supply this instead of `Lambda`/`Phi`/`Psi`.
#' @param model_error character. The method used to perturb the population so the factor model
#'   fits it imperfectly ("model error"): one of `"CB"` (Cudeck-Browne, the default), `"TKL"`
#'   (Tucker-Koopman-Linn), `"WB"` (Wu-Browne), or `"none"`. Model error is only applied when a
#'   target is supplied in `target_rmsea` or `target_cfi`; without one the population is exact,
#'   whatever `model_error`. Only used with a factor-model population (`Lambda`).
#' @param target_rmsea numeric. The population RMSEA the factor model should have relative to
#'   the perturbed population, a single number strictly in `(0, 1)`. Supplying it activates model
#'   error. Simulating from an exact population overstates recovery, so a realistic value (around
#'   `0.05`) is recommended for simulation studies (MacCallum, 2003). Default is `NULL` (an
#'   exact population; do not pass `0`). Required for `"CB"` and `"WB"`; optional for `"TKL"`.
#' @param target_cfi numeric. Only used with `model_error = "TKL"`: the population CFI to
#'   target, a single number strictly in `(0, 1)`, on its own or together with `target_rmsea`
#'   (TKL then trades the two off). Default is `NULL` (do not pass `1`). `"CB"` and `"WB"` target
#'   the RMSEA only.
#' @param marginals character. The marginal distribution of the drawn data: one
#'   of `"normal"` (the default), which draws normal marginals; `"empirical"`,
#'   which reproduces the population correlation while preserving the empirical
#'   marginals supplied in `marginal_data`; or `"VM"` (Vale-Maurelli) and `"IG"`
#'   (independent generator), which draw non-normal marginals with the target
#'   `skewness` and `kurtosis`.
#' @param marginal_data matrix or data frame. Only used with
#'   `marginals = "empirical"`, where it is required: a data set with one numeric
#'   column per variable (`p` columns), each with at least two distinct values,
#'   whose per-column distributions are resampled as the marginals of the drawn
#'   data. Its correlations are ignored. Default is `NULL`.
#' @param n_factors numeric. Only used with `marginals = "empirical"`: the number
#'   of factors the rank-matching reproduction fits. Default is `NULL`, in which
#'   case it is the number of columns of `Lambda` when the population is built
#'   from a factor model; it must be given when the population is supplied via `R`.
#' @param skewness numeric. Only used with `marginals = "VM"` or `"IG"`: the
#'   target marginal skewness, as a single value applied to every variable or a
#'   length-`p` vector. Default is `NULL` (0, a symmetric marginal). At least one
#'   of `skewness` or `kurtosis` must be given for these marginals.
#' @param kurtosis numeric. Only used with `marginals = "VM"` or `"IG"`: the
#'   target marginal excess kurtosis (0 for a normal marginal), as a single value
#'   applied to every variable or a length-`p` vector. Default is `NULL` (0).
#' @param force_pd logical. Used with `marginals = "VM"` and with Cudeck-Browne model error
#'   (`model_error = "CB"`). If the Vale-Maurelli intermediate correlation matrix (or, for CB,
#'   the perturbed population at a large `target_rmsea`) is not positive definite, `FALSE` (the
#'   default) rejects it with an error, while `TRUE` projects it to the nearest correlation
#'   matrix (with a warning). Has no effect for the `"TKL"` or `"WB"` methods.
#' @param categories numeric or list. Requests ordinal output by discretizing each
#'   variable into ordered categories. Either a count of equally probable categories
#'   (a single value applied to every variable or a length-`p` vector), or a
#'   length-`p` list of numeric vectors giving the marginal category proportions per
#'   variable (each strictly positive and summing to 1). Default is `NULL`, which
#'   returns the continuous data.
#' @param match character. Only used with `categories`: an assertion about how the
#'   categorization relates to the population correlation. With a normal latent, cutting
#'   at the normal-scale thresholds already leaves the population polychoric correlation
#'   of the categorized data equal to the target correlation, so both values compute the
#'   same thresholds and produce identical data whenever both are legal. `"thresholds"`
#'   (the default) also cuts the `"VM"` and `"IG"` draws, whose ordinal Pearson and
#'   polychoric correlations then both depart from the population; `"polychoric"` states
#'   that the polychoric match is required and therefore rejects non-normal marginals.
#'   Not available with `marginals = "empirical"`. The value is matched
#'   case-insensitively. Default is `NULL` (`"thresholds"` when `categories` is set).
#' @param missing character. An optional missing-data mechanism to impose on the drawn
#'   data: one of `"none"` (the default, complete data), `"MCAR"` (missing completely at
#'   random), `"MAR"` (missing at random, depending on another variable), or `"MNAR"`
#'   (missing not at random, depending on the variable's own value). Introduced values
#'   become `NA`. Every variable is holed, so under `"MAR"` each variable's predictor is
#'   itself subject to missingness: the mechanism is MAR given the *complete* data and is
#'   not ignorable for an analyst who sees only the observed data (see Details).
#' @param missing_prop numeric. Only used when `missing` is not `"none"`, where it is
#'   required: the target marginal proportion of missing values per variable, a single
#'   number strictly between 0 and 1. This is the expected rate; the realized rate of a
#'   given draw varies around it.
#' @param missing_strength numeric. Only used with `missing = "MAR"` or `"MNAR"`: the
#'   slope of the logistic missingness model, setting how strongly the missing
#'   probability depends on the predictor. Default is `NULL` (1, a moderate dependence);
#'   `0` removes the dependence (equivalent to MCAR at the same rate) and large
#'   magnitudes make missingness nearly deterministic.
#' @param missing_predictor integer or character. Only used with `missing = "MAR"`:
#'   which variable drives each variable's missingness, as one column index or variable
#'   name per variable. Each must reference another variable, not itself (so a single
#'   shared predictor is not allowed -- it would predict its own missingness). Default is
#'   `NULL`, in which case each variable's missingness is driven by the next variable
#'   cyclically (so variable order matters; supply this explicitly when the order is
#'   arbitrary).
#' @param n_datasets numeric. The number of datasets to draw. Default is 1. With
#'   more than one, a list of datasets is returned.
#' @param seed numeric. Optional seed for reproducible draws. When supplied, the
#'   caller's random-number stream is saved and restored, so the call leaves the
#'   global RNG state unchanged. Default is `NULL` (no seeding).
#' @param return_pop logical. If `TRUE`, return only the population correlation
#'   matrix and draw no data. Default is `FALSE`.
#'
#' @returns An object of class `efa_simulated`: a list with elements `data` (the simulated
#'   data -- an `N` by `p` numeric matrix, an integer matrix of category codes when `categories`
#'   is set, or a length-`n_datasets` list of these when `n_datasets > 1`; `NULL` when
#'   `return_pop = TRUE`), `population` (the `p` by `p` population correlation matrix drawn from,
#'   model-error-perturbed when requested; with `force_pd = TRUE` and `marginals = "VM"` it stays
#'   the target matrix, from which the realized correlations of the draw can drift), `model_error`
#'   (`NULL`, or a list of the method and
#'   the target and achieved RMSEA/CFI when model error was applied), and `settings`. Printing
#'   the object shows a compact summary.
#'
#' @references Cudeck, R., & Browne, M. W. (1992). Constructing a covariance matrix
#'   that yields a specified minimizer and a specified minimum discrepancy function
#'   value. Psychometrika, 57(3), 357-369. \doi{10.1007/BF02295424}
#'
#'   Fleishman, A. I. (1978). A method for simulating non-normal
#'   distributions. Psychometrika, 43(4), 521-532. \doi{10.1007/BF02293811}
#'
#'   Foldnes, N., & Grønneberg, S. (2015). How general is the Vale-Maurelli
#'   simulation approach? Psychometrika, 80(4), 1066-1083.
#'   \doi{10.1007/s11336-014-9414-0}
#'
#'   Foldnes, N., & Olsson, U. H. (2016). A simple simulation technique for
#'   nonnormal data with prespecified skewness, kurtosis, and covariance matrix.
#'   Multivariate Behavioral Research, 51(2-3), 207-219.
#'   \doi{10.1080/00273171.2015.1133274}
#'
#'   MacCallum, R. C. (2003). 2001 Presidential Address: Working with imperfect
#'   models. Multivariate Behavioral Research, 38(1), 113-139.
#'   \doi{10.1207/S15327906MBR3801_5}
#'
#'   Olsson, U. (1979). Maximum likelihood estimation of the polychoric correlation
#'   coefficient. Psychometrika, 44(4), 443-460. \doi{10.1007/BF02296207}
#'
#'   Olvera Astivia, O. L., & Zumbo, B. D. (2019). A note on the solution
#'   multiplicity of the Vale-Maurelli intermediate correlation equation. Journal
#'   of Educational and Behavioral Statistics, 44(2), 127-143.
#'   \doi{10.3102/1076998618803381}
#'
#'   Rubin, D. B. (1976). Inference and missing data. Biometrika, 63(3), 581-592.
#'   \doi{10.1093/biomet/63.3.581}
#'
#'   Ruscio, J., & Kaczetow, W. (2008). Simulating multivariate nonnormal data
#'   using an iterative algorithm. Multivariate Behavioral Research, 43(3),
#'   355-381. \doi{10.1080/00273170802285693}
#'
#'   Tucker, L. R., Koopman, R. F., & Linn, R. L. (1969). Evaluation of factor
#'   analytic research procedures by means of simulated correlation matrices.
#'   Psychometrika, 34(4), 421-459. \doi{10.1007/BF02290601}
#'
#'   Vale, C. D., & Maurelli, V. A. (1983). Simulating multivariate nonnormal
#'   distributions. Psychometrika, 48(3), 465-471. \doi{10.1007/BF02293687}
#'
#'   Wu, H., & Browne, M. W. (2015). Quantifying adventitious error in a covariance
#'   structure as a random effect. Psychometrika, 80(3), 571-600.
#'   \doi{10.1007/s11336-015-9451-3}
#'
#' @family data simulation
#'
#' @export
#'
#' @examples
#' # Build a population from a shipped loading pattern and factor correlations
#' Lambda <- population_models$loadings$baseline
#' Phi <- population_models$phis_3$moderate
#'
#' # Draw one normal dataset of 500 cases (the data live in $data)
#' sim <- efa_simulate(N = 500, Lambda = Lambda, Phi = Phi, seed = 42)
#' dim(sim$data)
#'
#' # Return only the population correlation matrix
#' R_pop <- efa_simulate(Lambda = Lambda, Phi = Phi, return_pop = TRUE)$population
#'
#' # Draw several datasets at once from a supplied correlation matrix
#' sims <- efa_simulate(N = 500, R = R_pop, n_datasets = 3, seed = 42)
#' length(sims$data)
#'
#' # Reproduce the population correlation but with skewed, empirical marginals
#' # (here from a chi-squared source with one column per variable)
#' src <- matrix(rchisq(200 * nrow(Lambda), df = 3), ncol = nrow(Lambda))
#' dat_emp <- efa_simulate(N = 500, Lambda = Lambda, Phi = Phi,
#'                         marginals = "empirical", marginal_data = src, seed = 42)
#'
#' # Draw skewed, leptokurtic data with the Vale-Maurelli method
#' dat_vm <- efa_simulate(N = 500, Lambda = Lambda, Phi = Phi, marginals = "VM",
#'                        skewness = 1.5, kurtosis = 4, seed = 42)
#'
#' # Draw five-category ordinal data whose polychoric correlation matches R
#' dat_ord <- efa_simulate(N = 500, Lambda = Lambda, Phi = Phi,
#'                         categories = 5, match = "polychoric", seed = 42)
#'
#' # Draw data with 15% missing at random, driven by a neighbouring item
#' dat_mar <- efa_simulate(N = 500, Lambda = Lambda, Phi = Phi, missing = "MAR",
#'                         missing_prop = 0.15, seed = 42)
#' colMeans(is.na(dat_mar$data))
#'
#' # Add realistic model error: a population the model fits with RMSEA of about .05
#' # (Cudeck-Browne, the default method; the achieved fit is reported)
#' sim_me <- efa_simulate(N = 500, Lambda = Lambda, Phi = Phi,
#'                        target_rmsea = 0.05, seed = 42)
#' sim_me$model_error$rmsea
efa_simulate <- function(N = NULL, Lambda = NULL, Phi = NULL, Psi = NULL,
                         R = NULL,
                         model_error = c("CB", "TKL", "WB", "none"),
                         target_rmsea = NULL, target_cfi = NULL,
                         marginals = c("normal", "empirical", "VM", "IG"),
                         marginal_data = NULL, n_factors = NULL,
                         skewness = NULL, kurtosis = NULL, force_pd = FALSE,
                         categories = NULL, match = NULL,
                         missing = c("none", "MCAR", "MAR", "MNAR"),
                         missing_prop = NULL, missing_strength = NULL,
                         missing_predictor = NULL,
                         n_datasets = 1L, seed = NULL, return_pop = FALSE) {

  # Eigenvalues below -tol mark the population as indefinite (not a valid
  # covariance); kept equal to the compiled core's default so both agree.
  tol <- 1e-8

  checkmate::assert_flag(return_pop)
  checkmate::assert_flag(force_pd)
  checkmate::assert_count(n_datasets, positive = TRUE)
  checkmate::assert_int(seed, null.ok = TRUE)

  marginals <- .match_arg_ci(marginals)

  # `model_error` selects the method used to perturb the population so a factor model fits it
  # imperfectly; `target_rmsea`/`target_cfi` set the amount and activate it. A target is what
  # turns model error on -- without one the population is exact, whatever `model_error`. CB and
  # WB target the RMSEA only; TKL can target the RMSEA and/or the CFI. Compatibility is checked
  # here; the perturbation itself runs once the population is built (below).
  model_error <- .match_arg_ci(model_error)
  # The targets are strictly inside (0, 1): a target RMSEA of 0 (or CFI of 1) means no misfit,
  # which is the exact population -- omit the target instead of asking a model-error method to
  # inject zero error (which would only degenerate the solvers).
  if (!is.null(target_rmsea) &&
      (!checkmate::test_number(target_rmsea) || target_rmsea <= 0 || target_rmsea >= 1)) {
    cli::cli_abort(
      c("{.arg target_rmsea} must be a single number in (0, 1).",
        "i" = "Omit {.arg target_rmsea} for an exact (zero-misfit) population."),
      class = "efa_simulate_input")
  }
  if (!is.null(target_cfi) &&
      (!checkmate::test_number(target_cfi) || target_cfi <= 0 || target_cfi >= 1)) {
    cli::cli_abort(
      c("{.arg target_cfi} must be a single number in (0, 1).",
        "i" = "Omit {.arg target_cfi} for an exact (perfect-fit) population."),
      class = "efa_simulate_input")
  }
  have_target <- !is.null(target_rmsea) || !is.null(target_cfi)
  do_model_error <- model_error != "none" && have_target
  if (model_error == "none" && have_target) {
    cli::cli_abort(
      c("{.arg target_rmsea} and {.arg target_cfi} require a model-error method.",
        "x" = "You supplied a target but left {.code model_error = \"none\"}.",
        "i" = "Choose {.code model_error = \"CB\"}, {.val TKL}, or {.val WB}."),
      class = "efa_simulate_input"
    )
  }
  # CB and WB target the RMSEA only.
  if (do_model_error && model_error %in% c("CB", "WB")) {
    if (is.null(target_rmsea)) {
      cli::cli_abort(
        c("{.code model_error = \"{model_error}\"} requires {.arg target_rmsea}.",
          "i" = "{.val CB} and {.val WB} target the RMSEA."),
        class = "efa_simulate_input"
      )
    }
    if (!is.null(target_cfi)) {
      cli::cli_abort(
        c("{.code model_error = \"{model_error}\"} targets the RMSEA only.",
          "x" = "You also supplied {.arg target_cfi}.",
          "i" = "Use {.code model_error = \"TKL\"} to target the CFI (with or without the RMSEA)."),
        class = "efa_simulate_input"
      )
    }
  }

  # marginal_data and n_factors configure the empirical rank-matching path only;
  # skewness and kurtosis configure the non-normal VM/IG paths, and force_pd the
  # Vale-Maurelli path only.
  if (marginals != "empirical" && (!is.null(marginal_data) || !is.null(n_factors))) {
    supplied <- c("marginal_data", "n_factors")[c(!is.null(marginal_data),
                                                  !is.null(n_factors))]
    cli::cli_abort(
      c("{.arg marginal_data} and {.arg n_factors} only apply when {.code marginals = \"empirical\"}.",
        "x" = "You set {.code marginals = \"{marginals}\"} but supplied {.arg {supplied}}."),
      class = "efa_simulate_input"
    )
  }
  if (!(marginals %in% c("VM", "IG")) &&
      (!is.null(skewness) || !is.null(kurtosis))) {
    supplied <- c("skewness", "kurtosis")[c(!is.null(skewness), !is.null(kurtosis))]
    cli::cli_abort(
      c("{.arg skewness} and {.arg kurtosis} only apply when {.code marginals} is {.val VM} or {.val IG}.",
        "x" = "You set {.code marginals = \"{marginals}\"} but supplied {.arg {supplied}}."),
      class = "efa_simulate_input"
    )
  }
  # force_pd rescues the Vale-Maurelli intermediate correlation matrix and, for Cudeck-Browne
  # model error, a population that cannot reach the target RMSEA while staying positive definite;
  # it has no role otherwise (the IG path needs a positive definite population outright).
  if (marginals != "VM" && !do_model_error && isTRUE(force_pd)) {
    cli::cli_abort(
      c("{.arg force_pd} only applies when {.code marginals = \"VM\"} or model error is requested.",
        "x" = "You set {.code force_pd = TRUE} but {.code marginals = \"{marginals}\"} with no model-error target."),
      class = "efa_simulate_input"
    )
  }

  # `categories` requests ordered-category (ordinal) output; `match` selects how
  # that categorization relates to the population correlation and is meaningless
  # without it.
  if (is.null(categories)) {
    if (!is.null(match)) {
      cli::cli_abort(
        c("{.arg match} only applies when {.arg categories} requests ordinal output.",
          "x" = "You set {.arg match} but left {.arg categories} unset, so continuous data are returned."),
        class = "efa_simulate_input"
      )
    }
  } else {
    match <- .match_arg_ci(match, c("thresholds", "polychoric"))
    # The polychoric correlation is defined by a bivariate-normal latent, so it can
    # only equal the target correlation when the latent is normal. Pairing it with
    # non-normal marginals is a contradiction rather than an approximation.
    if (match == "polychoric" && marginals != "normal") {
      cli::cli_abort(
        c("{.code match = \"polychoric\"} requires normal latents, but {.code marginals = \"{marginals}\"} draws non-normal marginals.",
          "x" = "The polychoric correlation assumes a bivariate-normal latent, so it cannot equal the target correlation under non-normal marginals.",
          "i" = "Use {.code match = \"thresholds\"} to cut the non-normal data (its ordinal Pearson correlation is then attenuated), or keep {.code marginals = \"normal\"}."),
        class = "efa_simulate_match_conflict"
      )
    }
    # Thresholds sit on the standard-normal scale, which categorizes the standardized
    # draws (normal, and the unit-variance VM/IG marginals) as intended. The empirical
    # draw instead carries the arbitrary location and scale of `marginal_data`, so
    # those thresholds would not reproduce the requested categories.
    if (marginals == "empirical") {
      cli::cli_abort(
        c("{.arg categories} is not supported with {.code marginals = \"empirical\"}.",
          "x" = "Ordinal thresholds sit on the standard-normal scale, but the empirical draw keeps the arbitrary scale of {.arg marginal_data}, so the requested categories are not reproduced.",
          "i" = "Use {.code marginals = \"normal\"}, {.val VM}, or {.val IG} to draw ordinal data with {.arg categories}."),
        class = "efa_simulate_input"
      )
    }
  }

  # `missing` selects an optional missing-data mechanism; `missing_prop` sets the
  # target marginal missing rate, `missing_strength` the strength of the value
  # dependence in MAR/MNAR, and `missing_predictor` the MAR predictor. Each companion
  # is meaningful only under the mechanisms that use it (compatibility is checked here;
  # value validation runs later, once the data dimension is known).
  missing <- .match_arg_ci(missing)
  if (missing == "none") {
    supplied <- c("missing_prop", "missing_strength", "missing_predictor")[
      c(!is.null(missing_prop), !is.null(missing_strength), !is.null(missing_predictor))]
    if (length(supplied)) {
      cli::cli_abort(
        c("Missing-data settings only apply when {.arg missing} requests a mechanism.",
          "x" = "You supplied {.arg {supplied}} but left {.code missing = \"none\"}."),
        class = "efa_simulate_input"
      )
    }
  } else {
    if (is.null(missing_prop)) {
      cli::cli_abort(
        c("{.arg missing_prop} is required when {.arg missing} is not {.val none}.",
          "i" = "Set the target marginal proportion of missing values, e.g. {.code missing_prop = 0.1}."),
        class = "efa_simulate_input"
      )
    }
    # MCAR removes values independently of the data, so it has no strength or predictor;
    # MNAR's predictor is fixed to each variable's own value, so none may be supplied.
    if (missing == "MCAR" &&
        (!is.null(missing_strength) || !is.null(missing_predictor))) {
      supplied <- c("missing_strength", "missing_predictor")[
        c(!is.null(missing_strength), !is.null(missing_predictor))]
      cli::cli_abort(
        c("{.code missing = \"MCAR\"} removes values completely at random.",
          "x" = "You supplied {.arg {supplied}}, but MCAR has no predictor to depend on."),
        class = "efa_simulate_input"
      )
    }
    if (missing == "MNAR" && !is.null(missing_predictor)) {
      cli::cli_abort(
        c("{.arg missing_predictor} does not apply when {.code missing = \"MNAR\"}.",
          "x" = "MNAR makes each variable's missingness depend on its own value, not another predictor."),
        class = "efa_simulate_input"
      )
    }
  }

  # The population is given either as a ready correlation matrix or as model
  # components, never both and never neither.
  have_lambda <- !is.null(Lambda)
  have_R <- !is.null(R)
  if (have_lambda == have_R) {
    cli::cli_abort(
      c("Specify the population in exactly one way.",
        "x" = "Provide either {.arg R}, or {.arg Lambda} (with optional {.arg Phi}/{.arg Psi})."),
      class = "efa_simulate_input"
    )
  }
  if (have_R && (!is.null(Phi) || !is.null(Psi))) {
    cli::cli_abort(
      c("{.arg Phi} and {.arg Psi} only apply when the population is built from {.arg Lambda}.",
        "x" = "You supplied {.arg R} together with {.arg Phi}/{.arg Psi}."),
      class = "efa_simulate_input"
    )
  }

  if (have_lambda) {

    Lambda <- as.matrix(Lambda)
    checkmate::assert_matrix(Lambda, mode = "numeric", any.missing = FALSE,
                             min.rows = 1, min.cols = 1)
    # assert_matrix's any.missing rejects NA/NaN but not Inf; a non-finite loading
    # would propagate NaN into the communalities and turn the Heywood guard below
    # into an uninformative `if (NA)` error, so reject it up front.
    if (!all(is.finite(Lambda))) {
      cli::cli_abort("{.arg Lambda} must contain only finite values.",
                     class = "efa_simulate_input")
    }
    p <- nrow(Lambda)
    m <- ncol(Lambda)

    if (is.null(Phi)) {
      Phi <- diag(m)
    } else {
      Phi <- as.matrix(Phi)
      checkmate::assert_matrix(Phi, mode = "numeric", nrows = m, ncols = m,
                               any.missing = FALSE)
      if (!isSymmetric(unname(Phi))) {
        cli::cli_abort(c("{.arg Phi} must be a symmetric {m} by {m} matrix."),
                       class = "efa_simulate_input")
      }
    }

    common <- Lambda %*% Phi %*% t(Lambda)

    if (is.null(Psi)) {
      # Unique variances that standardize the population to a correlation matrix.
      # An implied communality above 1 (a Heywood case) leaves no unique variance
      # to simulate from; a communality of exactly 1 gives a drawable, if singular,
      # population and is allowed.
      psi_vec <- 1 - diag(common)
      if (any(psi_vec < 0)) {
        cli::cli_abort(
          c("The factor model implies a communality above 1 (a Heywood case).",
            "x" = "There is no unique variance left to simulate from.",
            "i" = "Provide unique variances in {.arg Psi}, or reduce the loadings."),
          class = "efa_simulate_heywood"
        )
      }
      Psi_mat <- diag(psi_vec, p)
    } else if (is.matrix(Psi)) {
      checkmate::assert_matrix(Psi, mode = "numeric", nrows = p, ncols = p,
                               any.missing = FALSE)
      # A residual covariance must be finite and symmetric; a non-finite or
      # non-symmetric matrix would otherwise be silently averaged into shape by the
      # symmetrization below (as Phi and R are also checked to be symmetric).
      if (!all(is.finite(Psi)) || !isSymmetric(unname(Psi))) {
        cli::cli_abort("{.arg Psi} must be a finite, symmetric {p} by {p} matrix.",
                       class = "efa_simulate_input")
      }
      Psi_mat <- Psi
    } else {
      checkmate::assert_numeric(Psi, len = p, any.missing = FALSE, finite = TRUE)
      Psi_mat <- diag(Psi, p)
    }

    Sigma <- common + Psi_mat
    if (any(diag(Sigma) <= 0)) {
      cli::cli_abort(
        c("The assembled population has a non-positive variance and cannot be standardized.",
          "x" = "Check {.arg Lambda}, {.arg Phi}, and {.arg Psi}."),
        class = "efa_simulate_heywood"
      )
    }
    R_pop <- stats::cov2cor(Sigma)
    var_names <- rownames(Lambda)

  } else {

    R <- as.matrix(R)
    checkmate::assert_matrix(R, mode = "numeric", any.missing = FALSE,
                             min.rows = 1)
    # `any.missing` does not reject Inf, so check finiteness explicitly; otherwise a
    # non-finite R reaches eigen() below and aborts with an unclassed base error.
    if (nrow(R) != ncol(R) || !all(is.finite(R)) || !isSymmetric(unname(R))) {
      cli::cli_abort(
        c("{.arg R} must be a square, symmetric correlation matrix with finite values."),
        class = "efa_simulate_input")
    }
    p <- nrow(R)
    R_pop <- R
    # Variable names live on either dimension of a correlation matrix; fall back to
    # the column names when only those are present.
    var_names <- rownames(R)
    if (is.null(var_names)) var_names <- colnames(R)

  }

  # Remove any floating-point asymmetry so the matrix square root sees an exactly
  # symmetric matrix, and attach variable names.
  R_pop <- (R_pop + t(R_pop)) / 2
  if (is.null(var_names)) var_names <- paste0("V", seq_len(p))
  dimnames(R_pop) <- list(var_names, var_names)

  # Reject an indefinite population: a matrix with a markedly negative eigenvalue
  # is not a valid covariance to simulate from. This runs before the `return_pop`
  # exit too, so `return_pop = TRUE` only ever returns a valid correlation matrix.
  # A positive-semidefinite but singular population (smallest eigenvalue in
  # [-tol, 0]) is allowed and handled by the compiled core's eigen square-root
  # fallback.
  min_ev <- min(eigen(R_pop, symmetric = TRUE, only.values = TRUE)$values)
  if (min_ev < -tol) {
    cli::cli_abort(
      c("The population correlation matrix is not positive semi-definite.",
        "x" = "Its smallest eigenvalue is {.val {min_ev}}.",
        "i" = "A valid population model must have a positive semi-definite correlation matrix."),
      class = "efa_simulate_not_pd"
    )
  }

  # A fixed seed makes the whole call reproducible and independent of the number of parallel
  # workers. The caller's RNG stream is saved and restored afterwards -- or, if none existed yet,
  # the stream set.seed() creates is removed again -- so efa_simulate() leaves no side effect on
  # it (mirrors EFA()'s bootstrap seeding). Seeding here (before the optional model-error draw and
  # before the replicate loop) puts both under one umbrella; the model-error draw runs once in the
  # main process, and nothing between it and the replicate loop consumes random numbers, so a call
  # without model error draws the identical stream it did before this argument existed.
  .set_local_seed(seed)

  # Inject model error: perturb the population so the q-factor model fits it with the target
  # RMSEA (CB/WB) or RMSEA and/or CFI (TKL), and record the achieved fit. Requires the factor-model
  # population (there is no specified model to misfit against a bare `R`, and CB needs the model
  # Jacobian) and residual degrees of freedom. Runs once here so the perturbed population feeds the
  # `return_pop` exit, the VM/IG setup, and every replicate.
  model_error_info <- NULL
  if (do_model_error) {
    if (!have_lambda) {
      cli::cli_abort(
        c("Model error requires a factor-model population.",
          "x" = "You supplied {.arg R}, so there is no factor model to misfit against.",
          "i" = "Provide {.arg Lambda} (with optional {.arg Phi}/{.arg Psi}) to use {.arg model_error}."),
        class = "efa_simulate_input"
      )
    }
    q_me <- ncol(Lambda)
    if (.efa_df(p, q_me) <= 0) {
      cli::cli_abort(
        c("The {q_me}-factor model has no residual degrees of freedom, so no model error can be injected.",
          "i" = "Model error needs {.code ((p - q)^2 - (p + q)) / 2 > 0}."),
        class = "efa_simulate_input"
      )
    }
    # All three methods invert the population (solve()/chol()); a positive-semidefinite but
    # singular population (e.g. a communality of exactly 1) has no inverse, so require it to be
    # positive definite here (min_ev is the smallest eigenvalue computed above) rather than let a
    # base-R linear-algebra error escape unclassed.
    if (min_ev < tol) {
      cli::cli_abort(
        c("Model error requires a positive-definite population correlation matrix.",
          "x" = "Its smallest eigenvalue is {.val {min_ev}} (the population is singular or near-singular).",
          "i" = "Reduce the communalities (a communality of 1 makes the population singular), or omit the model-error target."),
        class = "efa_simulate_model_error"
      )
    }
    me <- .apply_model_error(R_pop, common, Sigma, q_me, model_error,
                            target_rmsea, target_cfi, force_pd, tol)
    R_pop <- (me$population + t(me$population)) / 2
    dimnames(R_pop) <- list(var_names, var_names)
    model_error_info <- me$model_error
    # The perturbed population must still be a valid (positive semi-definite) correlation matrix.
    min_ev <- min(eigen(R_pop, symmetric = TRUE, only.values = TRUE)$values)
    if (min_ev < -tol) {
      cli::cli_abort(
        c("The model-error population correlation matrix is not positive semi-definite.",
          "x" = "Its smallest eigenvalue is {.val {min_ev}}.",
          "i" = "Lower the target RMSEA/CFI, or use a different {.arg model_error} method."),
        class = "efa_simulate_not_pd"
      )
    }
  }

  if (return_pop) {
    return(.new_efa_simulated(
      data = NULL, population = R_pop, model_error = model_error_info,
      N = NA_integer_, n_datasets = NA_integer_, marginals = marginals,
      categories = categories, match = match, missing = missing, seed = seed))
  }

  checkmate::assert_count(N, positive = TRUE)

  # The empirical path needs a marginal source and a factor count for the
  # rank-matching reproduction; validate both against the population dimension.
  if (marginals == "empirical") {
    # The reproduction sorts and correlates the drawn cases, which needs at least
    # two of them: a single case has no correlation and collapses the working
    # matrix to a vector.
    if (N < 2L) {
      cli::cli_abort(
        c("{.arg N} must be at least 2 when {.code marginals = \"empirical\"}.",
          "x" = "You requested {N} case{?s}."),
        class = "efa_simulate_input"
      )
    }
    if (is.null(marginal_data)) {
      cli::cli_abort(
        c("{.arg marginal_data} is required when {.code marginals = \"empirical\"}.",
          "i" = "Supply a matrix or data frame with {p} column{?s} giving the empirical marginals."),
        class = "efa_simulate_input"
      )
    }
    marginal_data <- as.matrix(marginal_data)
    checkmate::assert_matrix(marginal_data, mode = "numeric", any.missing = FALSE,
                             min.rows = 1)
    # any.missing = FALSE rejects NA/NaN but not Inf, so check finiteness too.
    if (!all(is.finite(marginal_data))) {
      cli::cli_abort("{.arg marginal_data} must contain only finite values.",
                     class = "efa_simulate_input")
    }
    if (ncol(marginal_data) != p) {
      cli::cli_abort(
        c("{.arg marginal_data} must have one column per variable ({p}).",
          "x" = "It has {ncol(marginal_data)} column{?s}."),
        class = "efa_simulate_input"
      )
    }
    # Each column must offer at least two distinct values: a constant column is
    # not a distribution to resample and cannot reproduce a correlation (its
    # sample correlation is undefined), and a single-row source would resample a
    # length-one vector -- both of which yield silently wrong or degenerate draws.
    n_distinct <- apply(marginal_data, 2L, function(col) length(unique(col)))
    if (any(n_distinct < 2L)) {
      bad <- as.character(which(n_distinct < 2L))
      cli::cli_abort(
        c("Every column of {.arg marginal_data} must have at least two distinct values.",
          "x" = "{cli::qty(bad)}Constant column{?s}: {bad}."),
        class = "efa_simulate_input"
      )
    }

    # Default the reproduction's factor count to the model's when the population
    # was built from Lambda; a supplied R carries no such count, so require it.
    if (is.null(n_factors)) {
      if (have_lambda) {
        n_factors <- ncol(Lambda)
      } else {
        cli::cli_abort(
          c("{.arg n_factors} is required when the population is supplied via {.arg R}.",
            "i" = "It sets the factor count of the rank-matching reproduction."),
          class = "efa_simulate_input"
        )
      }
    }
    checkmate::assert_count(n_factors, positive = TRUE)
    # The iterative principal-axis step extracts n_factors factors from a p by p
    # matrix, which requires n_factors < p.
    if (n_factors >= p) {
      cli::cli_abort(
        c("{.arg n_factors} must be smaller than the number of variables ({p}).",
          "x" = "The reproduction uses {n_factors} factor{?s} for {p} variable{?s}.",
          "i" = "When not supplied, {.arg n_factors} defaults to the number of columns of {.arg Lambda}."),
        class = "efa_simulate_input"
      )
    }
  }

  # The VM/IG non-normal paths need target marginal moments; validate them,
  # recycle to length p, and precompute the Fleishman transform once. None of the
  # solves below draw random numbers, so they live outside the replicate loop.
  vm_ftable <- vm_Rinter <- ig_ftable <- ig_L <- ig_Lt <- NULL
  if (marginals %in% c("VM", "IG")) {
    if (is.null(skewness) && is.null(kurtosis)) {
      cli::cli_abort(
        c("{.arg skewness} and/or {.arg kurtosis} are required when {.code marginals} is {.val VM} or {.val IG}.",
          "i" = "The unset one defaults to 0 (the normal-marginal value)."),
        class = "efa_simulate_input"
      )
    }
    skewness <- .moment_vec(skewness, p, "skewness")
    kurtosis <- .moment_vec(kurtosis, p, "kurtosis")
    # Universal moment bound (Pearson, 1916): every distribution has excess
    # kurtosis at least skewness^2 - 2. Fleishman's cubic covers a smaller region
    # still, which the per-variable coefficient solve screens.
    infeasible <- kurtosis < skewness^2 - 2
    if (any(infeasible)) {
      bad <- as.character(which(infeasible))
      cli::cli_abort(
        c("The requested marginal moments are not attainable by any distribution.",
          "x" = "{cli::qty(bad)}Variable{?s} {bad}: excess {.arg kurtosis} is below {.arg skewness}^2 - 2.",
          "i" = "Every distribution must have excess kurtosis >= skewness^2 - 2 (Pearson, 1916)."),
        class = "efa_simulate_infeasible_moments"
      )
    }
  }

  if (marginals == "VM") {
    # Vale-Maurelli (1983): draw normal deviates at an intermediate correlation
    # chosen so the per-variable Fleishman cubic reproduces the target correlation.
    vm_ftable <- .fleishman_table(skewness, kurtosis)
    vm_Rinter <- .vm_intermediate_cor(R_pop, vm_ftable)
    # The pairwise cubic solves ignore joint positive-definiteness, so the
    # intermediate matrix can be indefinite even when R_pop is not. Reject it, or
    # project it to the nearest correlation matrix on request.
    min_ev_i <- min(eigen(vm_Rinter, symmetric = TRUE, only.values = TRUE)$values)
    if (min_ev_i < -tol) {
      if (!force_pd) {
        cli::cli_abort(
          c("The Vale-Maurelli intermediate correlation matrix is not positive definite.",
            "x" = "Its smallest eigenvalue is {.val {min_ev_i}}.",
            "i" = "Set {.code force_pd = TRUE} to project it to the nearest correlation matrix, or lower the target correlations or non-normality."),
          class = "efa_simulate_intermediate_not_pd"
        )
      }
      # Muffle cor.smooth's routine note; re-surface a single classed warning.
      vm_Rinter <- withCallingHandlers(
        psych::cor.smooth(vm_Rinter),
        warning = function(w) {
          if (grepl("smoothing was done", conditionMessage(w), fixed = TRUE)) {
            invokeRestart("muffleWarning")
          }
        }
      )
      vm_Rinter <- (vm_Rinter + t(vm_Rinter)) / 2
      cli::cli_warn(
        c("The Vale-Maurelli intermediate correlation matrix was not positive definite; it has been projected to the nearest correlation matrix.",
          "i" = "Projection was applied via {.fun psych::cor.smooth}; the achieved correlations may drift from the target."),
        class = "efa_simulate_pd_forced"
      )
    }

  } else if (marginals == "IG") {
    # Independent generator (Foldnes & Olsson, 2016): X = L Z with L the lower
    # Cholesky factor of R_pop and Z independent non-normal generators.
    ig_L <- tryCatch(t(chol(R_pop)), error = function(e) NULL)
    if (is.null(ig_L)) {
      cli::cli_abort(
        c("{.code marginals = \"IG\"} needs a positive definite population correlation matrix.",
          "x" = "The population is singular, so it has no Cholesky factor.",
          "i" = "Use {.code marginals = \"VM\"} or {.code marginals = \"empirical\"} for a singular population."),
        class = "efa_simulate_input"
      )
    }
    # Cumulant additivity turns the target marginal moments into the independent
    # generators' moments via lower-triangular solves; these depend on the
    # population, so a feasible request can still imply an infeasible generator.
    gen <- .ig_generator_moments(ig_L, skewness, kurtosis)
    # Report generator infeasibility in terms of the request the user made -- the
    # implied generator moments below are internal and not directly controllable.
    # Fleishman's region is tighter than the universal bound, so also translate a
    # solve failure (raised inside .fleishman_table) into this same message.
    ig_infeasible <- function() {
      cli::cli_abort(
        c("No independent generator can produce the requested marginals from this population.",
          "x" = "The moments the generators would need are not attainable.",
          "i" = "Lower the target correlations or the non-normality, or use {.code marginals = \"VM\"}."),
        class = "efa_simulate_infeasible_moments"
      )
    }
    if (any(gen$kurtosis < gen$skewness^2 - 2)) ig_infeasible()
    ig_ftable <- tryCatch(
      .fleishman_table(gen$skewness, gen$kurtosis),
      efa_simulate_infeasible_moments = function(e) ig_infeasible()
    )
    # Transpose of the Cholesky factor for the X = Z L' mixing; loop-invariant, so
    # precompute it once here rather than re-forming it in every replicate.
    ig_Lt <- t(ig_L)
  }

  # Resolve and validate the missing-data settings once (they draw no random numbers).
  # The mask itself is applied per replicate inside the draw so it shares each
  # replicate's reproducible stream. `miss_b` is the logistic slope (MAR/MNAR) and
  # `pred_idx` the per-variable MAR predictor column; both stay NULL when unused.
  miss_b <- NULL
  pred_idx <- NULL
  if (missing != "none") {
    # A proportion strictly inside (0, 1): 0 or 1 is not a mechanism to calibrate, and
    # the boundary would send qlogis() to +-Inf and break the intercept root-find.
    if (!checkmate::test_number(missing_prop) ||
        missing_prop <= 0 || missing_prop >= 1) {
      cli::cli_abort(
        c("{.arg missing_prop} must be a single number strictly between 0 and 1.",
          "x" = "It is the target marginal proportion of missing values."),
        class = "efa_simulate_input"
      )
    }
    if (missing %in% c("MAR", "MNAR")) {
      # The slope defaults to 1 (a moderate dependence); 0 removes the dependence, which
      # is MCAR at the same rate; large magnitudes make missingness near-deterministic.
      miss_b <- if (is.null(missing_strength)) 1 else missing_strength
      if (!checkmate::test_number(miss_b, finite = TRUE)) {
        cli::cli_abort(
          c("{.arg missing_strength} must be a single finite number.",
            "i" = "It is the slope of the logistic missingness model; larger magnitudes make missingness depend more strongly on the predictor."),
          class = "efa_simulate_input"
        )
      }
    }
    if (missing == "MAR") {
      pred_idx <- .resolve_missing_predictor(missing_predictor, p, var_names)
    }
  }

  # With `categories`, the drawn data are cut into ordered categories at thresholds
  # on the standard-normal scale. The thresholds draw no random numbers and are the
  # same for every replicate, so they are built once here; the cutting itself runs
  # after the draw (below), keeping the replicate RNG streams untouched.
  ord_thresholds <- if (!is.null(categories)) {
    .ordinal_thresholds(categories, p)
  }
  # The Vale-Maurelli draw is a cubic map of a standard normal, so the normal-scale
  # thresholds would cut it at the wrong quantiles; push them through the same cubic
  # so the requested proportions are reproduced on the non-normal scale. A variable
  # whose cubic does not increase across its own thresholds, or folds tail mass back
  # below an outermost cut, cannot be mapped faithfully and keeps the normal ones.
  if (!is.null(ord_thresholds) && marginals == "VM") {
    mapped <- .vm_thresholds(ord_thresholds, vm_ftable)
    ord_thresholds <- mapped$thresholds
    if (length(mapped$non_monotone)) {
      kept <- var_names[mapped$non_monotone]
      cli::cli_warn(
        c("{cli::qty(kept)}Variable{?s} {.val {kept}} kept standard-normal category thresholds.",
          "x" = "{cli::qty(kept)}{?Its/Their} Fleishman polynomial{?s} {?is/are} not increasing over the requested thresholds and the tails beyond them, so mapping the thresholds would misassign category proportions.",
          "i" = "The achieved category proportions of {cli::qty(kept)}{?this/these} variable{?s} depart from the request; less extreme category proportions or a milder {.arg skewness}/{.arg kurtosis} restore the mapping."),
        class = "efa_simulate_threshold_fallback"
      )
    }
  }

  # Draw each dataset in its own future; future.seed = TRUE assigns every
  # replicate a reproducible L'Ecuyer-CMRG stream, so the result is identical at
  # any number of workers (mirrors the np-boot fit in .boot_fun()).
  datasets <- future.apply::future_lapply(
    seq_len(n_datasets),
    function(i) {
      dat <- if (marginals == "empirical") {
        .simulate_cfm_empirical(R_pop, marginal_data, n_factors, N)
      } else if (marginals == "VM") {
        # normal deviates at the intermediate correlation, then the cubic map.
        .fleishman_transform(.simulate_cfm_mvn(vm_Rinter, N), vm_ftable)
      } else if (marginals == "IG") {
        # independent Fleishman generators mixed by the Cholesky factor: X = Z L'.
        z <- .fleishman_transform(matrix(stats::rnorm(N * p), N, p), ig_ftable)
        z %*% ig_Lt
      } else {
        .simulate_cfm_mvn(R_pop, N)
      }
      colnames(dat) <- var_names
      # Introduce missing values on the drawn (continuous) data; runs inside the future
      # so the mask draws from this replicate's stream (reproducible across worker
      # counts). When `categories` is set, the NAs propagate into the codes below.
      if (missing != "none") {
        dat <- .apply_missing(dat, missing, missing_prop, miss_b, pred_idx)
      }
      dat
    },
    future.seed = TRUE
  )

  # Discretize to ordered categories when requested. Cutting is deterministic given
  # the continuous draw, so it runs in the main process (leaving the draw and its
  # per-replicate streams unchanged) with the same thresholds for every replicate. A
  # draw can leave a requested category empty -- more likely at small `N` or, under
  # `match = "thresholds"`, with strongly non-normal marginals -- which destabilizes
  # the polychoric correlation; report any such variable once.
  if (!is.null(ord_thresholds)) {
    k_req <- vapply(ord_thresholds, function(tau) length(tau) + 1L, integer(1L))
    datasets <- lapply(datasets, .ordinalize, thresholds = ord_thresholds)
    empty <- character(0L)
    for (d in datasets) {
      # Count non-empty categories per variable with an O(N) tabulate over the 1..K
      # codes; tabulate ignores any NA codes (from masked values), so no missing-value
      # filtering is needed.
      filled <- vapply(seq_len(ncol(d)),
                       function(j) sum(tabulate(d[, j], nbins = k_req[j]) > 0L),
                       integer(1L))
      empty <- union(empty, colnames(d)[filled < k_req])
    }
    if (length(empty)) {
      cli::cli_warn(
        c("{cli::qty(empty)}Variable{?s} {.val {empty}} did not fill every requested category.",
          "i" = "Empty or near-empty categories destabilize the polychoric correlation and the factor analysis; increase {.arg N} or request fewer categories."),
        class = "efa_simulate_empty_category"
      )
    }
  }

  data <- if (n_datasets == 1L) datasets[[1L]] else datasets
  .new_efa_simulated(
    data = data, population = R_pop, model_error = model_error_info,
    N = N, n_datasets = n_datasets, marginals = marginals,
    categories = categories, match = match, missing = missing, seed = seed)

}


# Construct the `efa_simulated` result object: the drawn `data` (an N x p matrix, a length-
# n_datasets list of them, or NULL for return_pop = TRUE), the `population` correlation matrix
# actually drawn from (model-error-perturbed when requested), the `model_error` record (NULL when
# the population is exact), and a `settings` list for provenance and printing.
.new_efa_simulated <- function(data, population, model_error, N, n_datasets,
                               marginals, categories, match, missing, seed) {
  structure(
    list(
      data = data,
      population = population,
      model_error = model_error,
      settings = list(N = N, n_datasets = n_datasets, marginals = marginals,
                      categories = categories, match = match, missing = missing,
                      seed = seed)
    ),
    class = "efa_simulated"
  )
}


# Ruscio-Kaczetow iterative rank-matching: draw `N` cases that reproduce the
# target correlation `R` while carrying the empirical column marginals of `x`
# (each column is resampled with replacement from the matching column of `x`).
# The target correlation and the marginal source are separate arguments so an
# arbitrary population `R` can be paired with any supplied marginals; the shared
# rank-matching loop otherwise follows Ruscio & Kaczetow (2008). `stats::cor()`,
# `sort()`, `sort.list()`, and `.paf_iter()` draw no random numbers, so with
# `R = stats::cor(x, method = cor_method)` the stream, and hence the output, is
# identical to the observed-data special case that `CD()` uses to draw its
# comparison populations. `x` may be a data frame or a matrix (`CD()` passes its
# raw input, `efa_simulate()` a matrix); `x[, i]` and `sample()` treat both
# identically.
# Ruscio, J., & Kaczetow, W. (2008). Simulating multivariate nonnormal data
# using an iterative algorithm. Multivariate Behavioral Research, 43(3),
# 355-381. doi: 10.1080/00273170802285693
.simulate_cfm_empirical <- function(R, x, n_factors, N, cor_method = "pearson",
                                    max_trials = 5, initial_multiplier = 1,
                                    max_iter = 100) {

  k <- ncol(x)
  sim_dat <- matrix(0, nrow = N, ncol = k)         # Simulated data
  dists <- matrix(0, nrow = N, ncol = k)           # Each variable's marginal
  best_RMSE <- 1                                    # Lowest RMSE correlation
  t_no_impr <- 0                                    # Trial counter

  # Resample the supplied empirical marginals (sorted) for each variable.
  for (i in 1:k) {
    dists[, i] <- sort(sample(x[, i], size = N, replace = TRUE))
  }

  # The target correlation is supplied; keep a working and a best-so-far copy so
  # the no-improvement branch below is safe even if the first iteration does not
  # lower the initial RMSE (best_RMSE = 1).
  R_inter <- R
  R_best <- R_inter
  res_best <- matrix(0, nrow = k, ncol = k)

  # Random normal shared and unique components; factor loadings applied below.
  shared_comp <- matrix(stats::rnorm(N * n_factors, 0, 1), nrow = N,
                        ncol = n_factors)
  unique_comp <- matrix(stats::rnorm(N * k, 0, 1), nrow = N, ncol = k)
  shared_load <- matrix(0, nrow = k, ncol = n_factors)
  unique_load <- matrix(0, nrow = k, ncol = 1)

  # Iterate until a set number of trials pass without improving the RMSE between
  # the reproduced and target correlations.
  while (t_no_impr < max_trials) {

    # Factor loadings for the current intermediate correlation.
    L <- suppressWarnings(.paf_iter(rep(1, k), criterion = .001, R = R_inter,
                   n_fac = n_factors, abs_eig = TRUE, crit_type = 1,
                   max_iter = max_iter)$L)

    shared_load[, 1:n_factors] <- L

    # get rid of Heywood cases
    shared_load[shared_load > 1] <- 1
    shared_load[shared_load < -1] <- -1

    if (shared_load[1, 1] < 0) {
      shared_load <- shared_load * -1
    }

    for (i in seq_len(k)) {
      if (sum(shared_load[i, ] * shared_load[i, ]) < 1) {
        unique_load[i, 1] <-
          (1 - sum(shared_load[i, ] * shared_load[i, ]))
      } else {
        unique_load[i, 1] <- 0
      }
    }

    unique_load <- sqrt(unique_load)

    for (i in seq_len(k)) {
      sim_dat[, i] <- (shared_comp %*% t(shared_load))[, i] +
        unique_comp[, i] * unique_load[i, 1]
    }

    # Replace the normal scores with the empirical marginals via rank matching.
    for (i in seq_len(k)) {
      sim_dat <- sim_dat[sort.list(sim_dat[, i]), ]
      sim_dat[, i] <- dists[, i]
    }

    # RMSE between reproduced and target correlations; accept or shrink the step.
    R_rep <- stats::cor(sim_dat, method = cor_method)
    R_res <- R - R_rep
    RMSE <- sqrt(sum(R_res[lower.tri(R_res)] * R_res[lower.tri(R_res)]) /
                   (.5 * (k * k - k)))

    if (RMSE < best_RMSE) {
      best_RMSE <- RMSE
      R_best <- R_inter
      res_best <- R_res
      R_inter <- R_inter + initial_multiplier * R_res
      t_no_impr <- 0
    } else {
      t_no_impr <- t_no_impr + 1
      current_multiplier <- initial_multiplier * .5 ^ t_no_impr
      R_inter <- R_best + current_multiplier * res_best
    }
  }

  # Reconstruct the data set from the best intermediate correlation.
  L <- suppressWarnings(.paf_iter(rep(1, k), criterion = .001, R = R_best,
                n_fac = n_factors, abs_eig = TRUE, crit_type = 1,
                max_iter = max_iter)$L)
  shared_load[, seq_len(n_factors)] <- L

  shared_load[shared_load > 1] <- 1
  shared_load[shared_load < -1] <- -1
  if (shared_load[1, 1] < 0) {
    shared_load <- shared_load * -1
  }

  for (i in seq_len(k)) {
    if (sum(shared_load[i, ] * shared_load[i, ]) < 1) {
      unique_load[i, 1] <-
        (1 - sum(shared_load[i, ] * shared_load[i, ]))
    } else {
      unique_load[i, 1] <- 0
    }
  }

  unique_load <- sqrt(unique_load)

  for (i in seq_len(k)) {
    sim_dat[, i] <- (shared_comp %*% t(shared_load))[, i] +
      unique_comp[, i] * unique_load[i, 1]
  }

  for (i in seq_len(k)) {
    sim_dat <- sim_dat[sort.list(sim_dat[, i]), ]
    sim_dat[, i] <- dists[, i]
  }

  sim_dat
}


# Recycle a length-1 or length-p marginal-moment argument to length p (NULL gives
# zeros, the normal-marginal value), validating type and finiteness.
.moment_vec <- function(x, p, arg) {
  if (is.null(x)) return(rep(0, p))
  # A matrix passes is.numeric() and, if it has p entries, the length check below,
  # so reject any dimensioned object here rather than flatten it column-major.
  if (!is.numeric(x) || !is.null(dim(x)) || !all(is.finite(x))) {
    cli::cli_abort("{.arg {arg}} must be a finite numeric vector.",
                   class = "efa_simulate_input")
  }
  if (length(x) == 1L) return(rep(as.numeric(x), p))
  if (length(x) != p) {
    cli::cli_abort(
      c("{.arg {arg}} must have length 1 or {p} (one value per variable).",
        "x" = "It has length {length(x)}."),
      class = "efa_simulate_input"
    )
  }
  as.numeric(x)
}

# Fleishman (1978) cubic coefficients (a, b, c, d) with a = -c that transform a
# standard normal Z into a unit-variance variable with the target skewness and
# excess kurtosis via a + bZ + cZ^2 + dZ^3. The three power-method equations are
# minimised as a sum of squares (as in lavaan/semTools); a non-zero minimum means
# the target lies outside Fleishman's feasible region and is rejected.
# Fleishman, A. I. (1978). A method for simulating non-normal distributions.
# Psychometrika, 43(4), 521-532.
.fleishman_coef <- function(skewness, kurtosis) {
  sys <- function(x) {
    b. <- x[1L]; c. <- x[2L]; d. <- x[3L]
    eq1 <- b.^2 + 6 * b. * d. + 2 * c.^2 + 15 * d.^2 - 1
    eq2 <- 2 * c. * (b.^2 + 24 * b. * d. + 105 * d.^2 + 2) - skewness
    eq3 <- 24 * (b. * d. + c.^2 * (1 + b.^2 + 28 * b. * d.) +
                 d.^2 * (12 + 48 * b. * d. + 141 * c.^2 + 225 * d.^2)) - kurtosis
    eq1^2 + eq2^2 + eq3^2
  }
  out <- stats::nlminb(start = c(1, 0, 0), objective = sys, scale = 10,
                       control = list(trace = 0L))
  if (out$objective > 1e-6) {
    cli::cli_abort(
      c("The marginal moments (skewness {.val {skewness}}, excess kurtosis {.val {kurtosis}}) are not representable by a Fleishman polynomial.",
        "i" = "Fleishman's cubic covers a smaller region than excess kurtosis >= skewness^2 - 2; choose moments further inside it."),
      class = "efa_simulate_infeasible_moments"
    )
  }
  b. <- out$par[1L]; c. <- out$par[2L]; d. <- out$par[3L]
  c(-c., b., c., d.)
}

# Fleishman coefficients for every variable as a p x 4 matrix of (a, b, c, d).
.fleishman_table <- function(skewness, kurtosis) {
  t(vapply(seq_along(skewness),
           function(i) .fleishman_coef(skewness[i], kurtosis[i]),
           numeric(4L)))
}

# Apply each variable's Fleishman cubic to the columns of Z (a p x 4 coefficient
# table): X_i = a_i + b_i Z_i + c_i Z_i^2 + d_i Z_i^3.
.fleishman_transform <- function(Z, ftable) {
  for (i in seq_len(ncol(Z))) {
    z <- Z[, i]
    Z[, i] <- ftable[i, 1L] + ftable[i, 2L] * z + ftable[i, 3L] * z^2 +
      ftable[i, 4L] * z^3
  }
  Z
}

# Vale-Maurelli (1983) intermediate normal correlation for each variable pair:
# solve the cubic in rho relating the intermediate normal correlation to the
# target correlation R_ij given the two variables' Fleishman coefficients. The
# cubic can admit several real roots in [-1, 1]; take the one nearest the target
# (Olvera Astivia & Zumbo, 2019). Base polyroot() is used rather than a numeric
# minimiser so all roots are found explicitly.
# Vale, C. D., & Maurelli, V. A. (1983). Simulating multivariate nonnormal
# distributions. Psychometrika, 48(3), 465-471.
.vm_intermediate_cor <- function(R, ftable) {
  p <- ncol(R)
  b <- ftable[, 2L]; cc <- ftable[, 3L]; d <- ftable[, 4L]
  ICOR <- diag(p)
  for (j in seq_len(p - 1L)) {
    for (i in (j + 1L):p) {
      target <- R[i, j]
      if (target == 0) next                # an uncorrelated pair stays at rho = 0
      # coefficients of c0 + c1 rho + c2 rho^2 + c3 rho^3
      coefs <- c(-target,
                 b[i] * b[j] + 3 * b[i] * d[j] + 3 * d[i] * b[j] + 9 * d[i] * d[j],
                 2 * cc[i] * cc[j],
                 6 * d[i] * d[j])
      # drop trailing zero coefficients so polyroot sees the true polynomial degree
      while (length(coefs) > 1L && coefs[length(coefs)] == 0) {
        coefs <- coefs[-length(coefs)]
      }
      roots <- polyroot(coefs)
      real <- Re(roots)[abs(Im(roots)) < 1e-8]
      if (!length(real)) real <- Re(roots)[which.min(abs(Im(roots)))]
      in_range <- real[real >= -1 & real <= 1]
      cand <- if (length(in_range)) in_range else real
      # A degenerate cubic (all higher coefficients zero) leaves no root; fall back
      # to the target so the assignment below never sees a zero-length value.
      if (!length(cand)) cand <- target
      rho <- cand[which.min(abs(cand - target))]
      ICOR[i, j] <- ICOR[j, i] <- rho
    }
  }
  ICOR
}

# Independent-generator moments (Foldnes & Olsson, 2016). With X = L Z, L the
# lower Cholesky factor of the target correlation and Z independent with mean 0
# and variance 1, the additivity and homogeneity of cumulants give
# skewness(X) = (L^3) skewness(Z) and excess-kurtosis(X) = (L^4) kurtosis(Z);
# solve the lower-triangular systems for the generators' moments.
# Foldnes, N., & Olsson, U. H. (2016). A simple simulation technique for
# nonnormal data with prespecified skewness, kurtosis, and covariance matrix.
# Multivariate Behavioral Research, 51(2-3), 207-219.
.ig_generator_moments <- function(L, skewness, kurtosis) {
  list(skewness = forwardsolve(L^3, skewness),
       kurtosis = forwardsolve(L^4, kurtosis))
}


# Interior category thresholds on the standard-normal scale for each variable, from
# a `categories` specification: either a scalar or length-p count of equally probable
# categories, or a length-p list of per-variable category proportions. Cutting a
# standard normal at tau = qnorm(cumsum(prop)), dropping the trailing +Inf, yields
# categories with population proportions `prop` (Olsson, 1979). Returns a length-p
# list of numeric threshold vectors (a variable with K categories has K - 1 of them).
# Olsson, U. (1979). Maximum likelihood estimation of the polychoric correlation
# coefficient. Psychometrika, 44(4), 443-460.
.ordinal_thresholds <- function(categories, p) {
  if (is.list(categories)) {
    if (length(categories) != p) {
      cli::cli_abort(
        c("{.arg categories}, as a list of category proportions, must have one element per variable ({p}).",
          "x" = "It has {length(categories)} element{?s}."),
        class = "efa_simulate_input"
      )
    }
    lapply(seq_along(categories), function(j) {
      prop <- categories[[j]]
      # Reject a dimensioned object (matrix/array): it passes is.numeric() and would
      # be flattened column-major by cumsum() instead of read as a proportion vector
      # (mirrors the guard in .moment_vec()).
      if (!is.numeric(prop) || !is.null(dim(prop)) || length(prop) < 2L ||
          !all(is.finite(prop))) {
        cli::cli_abort(
          c("Each element of {.arg categories} must be a finite numeric vector of at least two category proportions.",
            "x" = "Element {j} is not."),
          class = "efa_simulate_input"
        )
      }
      if (any(prop <= 0)) {
        cli::cli_abort(
          c("Category proportions in {.arg categories} must be strictly positive.",
            "x" = "Element {j} has a non-positive proportion (a zero proportion is an empty category with no threshold)."),
          class = "efa_simulate_input"
        )
      }
      if (abs(sum(prop) - 1) > 1e-6) {
        cli::cli_abort(
          c("Category proportions in {.arg categories} must sum to 1.",
            "x" = "Element {j} sums to {.val {sum(prop)}}."),
          class = "efa_simulate_input"
        )
      }
      # Renormalize onto the exact simplex: floating-point slack in `prop` (within
      # the tolerance above) could otherwise push an interior cumulative proportion
      # to or past 1, which qnorm would turn into a NaN/Inf threshold. The outer gate
      # has already rejected any genuine mismatch.
      prop <- prop / sum(prop)
      stats::qnorm(cumsum(prop)[-length(prop)])
    })
  } else {
    # Reject a dimensioned object (matrix/array): it passes is.numeric() and, with p
    # entries, the length check below, so it would be flattened column-major into
    # per-variable counts (mirrors the guard in .moment_vec()).
    if (!is.numeric(categories) || !is.null(dim(categories)) ||
        !all(is.finite(categories)) || any(categories != round(categories))) {
      cli::cli_abort(
        "{.arg categories} must be a whole number of categories per variable, or a list of category proportions.",
        class = "efa_simulate_input"
      )
    }
    if (length(categories) != 1L && length(categories) != p) {
      cli::cli_abort(
        c("{.arg categories}, as a count of categories, must have length 1 or {p} (one per variable).",
          "x" = "It has length {length(categories)}."),
        class = "efa_simulate_input"
      )
    }
    if (any(categories < 2L)) {
      cli::cli_abort(
        c("{.arg categories} must request at least two categories per variable.",
          "x" = "The smallest requested count is {.val {min(categories)}}."),
        class = "efa_simulate_input"
      )
    }
    counts <- if (length(categories) == 1L) rep(categories, p) else categories
    lapply(counts, function(k) stats::qnorm(seq_len(k - 1L) / k))
  }
}


# Exact probability mass a Fleishman cubic f assigns at or below the cut f(t0): decompose
# the real line at the real roots of f(z) = f(t0) and add the standard-normal mass of every
# interval on which f - f(t0) <= 0 (sign tested at an interior point). Used to detect tail
# fold-back: for a non-monotone cubic, P(f(Z) <= f(t0)) can differ from pnorm(t0) even when
# f increases across the thresholds themselves, because a tail segment beyond t0 dips back
# below the cut. Returns NA for a degenerate (constant) map.
.vm_cut_mass <- function(t0, b, cc, d) {
  # f(z) - f(t0) in increasing degree; the constant Fleishman term cancels
  co <- c(-(b * t0 + cc * t0^2 + d * t0^3), b, cc, d)
  while (length(co) > 1L && co[length(co)] == 0) co <- co[-length(co)]
  if (length(co) <= 1L) return(NA_real_)
  r <- polyroot(co)
  re <- sort(unique(Re(r[abs(Im(r)) < 1e-8])))
  cuts <- c(-Inf, re, Inf)
  mass <- 0
  for (i in seq_len(length(cuts) - 1L)) {
    lo <- cuts[i]; hi <- cuts[i + 1L]
    mid <- if (is.infinite(lo) && is.infinite(hi)) 0
    else if (is.infinite(lo)) hi - 1
    else if (is.infinite(hi)) lo + 1
    else (lo + hi) / 2
    if (b * mid + cc * mid^2 + d * mid^3 -
        (b * t0 + cc * t0^2 + d * t0^3) <= 0) {
      mass <- mass + (stats::pnorm(hi) - stats::pnorm(lo))
    }
  }
  mass
}

# Map standard-normal ordinal thresholds onto the Vale-Maurelli scale. The VM draw is
# X = f(Z) with f the per-variable Fleishman cubic a + bZ + cZ^2 + dZ^3, so where f is
# order-preserving over (essentially) all the normal mass, P(X <= f(tau)) = P(Z <= tau) and
# cutting X at f(tau) reproduces exactly the proportions tau was built for. Two screens
# check this. First, f must increase across the requested thresholds themselves, or the cut
# points come out of order: the derivative b + 2cz + 3dz^2 is a parabola, so its minimum
# over the closed threshold range is attained at an endpoint or, when the parabola opens
# upward (d > 0), at the vertex -c/(3d) if that falls inside. Second, the tails beyond the
# outermost thresholds must not fold back below the cuts: a non-monotone cubic (every
# platykurtic marginal has d < 0; a skewed one can dip on one side) re-crosses an outer cut
# value at some more extreme z, and the mass beyond that re-crossing is silently assigned to
# the wrong side -- 75% of a 1.4% top category at excess kurtosis -1, 84% of a 5% bottom
# category at skewness 1.5 / excess kurtosis 3. The exact misassigned mass at the two
# outermost thresholds is compared against the requested outer-category proportions
# (checking the outermost cuts suffices: an interior cut's re-crossing lies further out and
# folds strictly less mass) and a variable fails when it is off by more than 5% of the
# proportion -- departures below that are negligible against the error the normal-scale
# fallback itself incurs. The indices of the failing variables are returned; those keep
# their normal-scale thresholds. Infinite thresholds are passed through: an increasing map
# sends +-Inf to itself, and evaluating the polynomial there could produce a NaN cut point.
# Fleishman, A. I. (1978). A method for simulating non-normal distributions.
# Psychometrika, 43(4), 521-532.
.vm_thresholds <- function(thresholds, ftable) {
  non_monotone <- integer(0L)
  for (j in seq_along(thresholds)) {
    tau <- thresholds[[j]]
    fin <- is.finite(tau)
    if (!any(fin)) next                 # nothing to map
    b <- ftable[j, 2L]; cc <- ftable[j, 3L]; d <- ftable[j, 4L]
    rng <- range(tau[fin])
    z <- rng
    if (d > 0) {
      zv <- -cc / (3 * d)
      if (zv > z[1L] && zv < z[2L]) z <- c(z, zv)
    }
    if (any(b + 2 * cc * z + 3 * d * z^2 <= 0)) {
      non_monotone <- c(non_monotone, j)
      next
    }
    p_lo <- stats::pnorm(rng[1L])
    p_hi <- 1 - stats::pnorm(rng[2L])
    m_lo <- .vm_cut_mass(rng[1L], b, cc, d)
    m_hi <- .vm_cut_mass(rng[2L], b, cc, d)
    if (is.na(m_lo) || is.na(m_hi) ||
        abs(m_lo - p_lo) > 0.05 * p_lo ||
        abs((1 - m_hi) - p_hi) > 0.05 * p_hi) {
      non_monotone <- c(non_monotone, j)
      next
    }
    tau[fin] <- .fleishman_transform(matrix(tau[fin], ncol = 1L),
                                     ftable[j, , drop = FALSE])[, 1L]
    thresholds[[j]] <- tau
  }
  list(thresholds = thresholds, non_monotone = non_monotone)
}


# Discretize the columns of a continuous draw into ordered category codes 1..K,
# cutting each column at its interior thresholds (a length-p list). findInterval
# returns the count of thresholds at or below each value (0..K-1); +1 shifts to
# 1-based codes.
.ordinalize <- function(dat, thresholds) {
  out <- matrix(NA_integer_, nrow(dat), ncol(dat), dimnames = dimnames(dat))
  for (j in seq_len(ncol(dat))) {
    out[, j] <- findInterval(dat[, j], thresholds[[j]]) + 1L
  }
  out
}


# Resolve the MAR predictor specification to a length-p vector of column indices: for
# each variable, which other variable's value drives its missingness. NULL defaults to
# the next variable cyclically (variable p is driven by variable 1), so every variable
# has a distinct, observed predictor. A supplied spec (column indices or names, one per
# variable) must reference existing columns other than the variable itself -- a variable
# predicting its own missingness would be MNAR, not MAR (which also rules out a single
# recycled predictor, as it would predict itself). The cyclic default is order-dependent,
# so an explicit spec is preferable when the column order is arbitrary.
.resolve_missing_predictor <- function(spec, p, var_names) {
  if (p < 2L) {
    cli::cli_abort(
      c("{.code missing = \"MAR\"} needs at least two variables.",
        "x" = "With a single variable there is no other predictor for the missingness.",
        "i" = "Use {.code missing = \"MNAR\"} or {.val MCAR} for a single variable."),
      class = "efa_simulate_input"
    )
  }
  # Cyclic next variable: 2, 3, ..., p, 1 (never a variable's own position).
  if (is.null(spec)) return(c(seq_len(p)[-1L], 1L))

  if (is.character(spec)) {
    idx <- match(spec, var_names)
    if (anyNA(idx)) {
      bad <- spec[is.na(idx)]
      cli::cli_abort(
        c("{.arg missing_predictor} names a variable not in the data.",
          "x" = "{cli::qty(bad)}Unknown name{?s}: {.val {bad}}."),
        class = "efa_simulate_input"
      )
    }
  } else {
    # Reject a dimensioned or non-integer spec rather than flatten/round it silently.
    if (!is.numeric(spec) || !is.null(dim(spec)) || !all(is.finite(spec)) ||
        any(spec != round(spec))) {
      cli::cli_abort(
        "{.arg missing_predictor} must be whole-number column indices or variable names.",
        class = "efa_simulate_input"
      )
    }
    idx <- as.integer(spec)
    if (any(idx < 1L | idx > p)) {
      cli::cli_abort(
        c("{.arg missing_predictor} indexes a column outside 1 to {p}.",
          "x" = "It has value{?s} {.val {idx[idx < 1L | idx > p]}}."),
        class = "efa_simulate_input"
      )
    }
  }

  if (length(idx) != p) {
    cli::cli_abort(
      c("{.arg missing_predictor} must give one predictor per variable ({p}).",
        "x" = "It has length {length(idx)}.",
        "i" = "A single predictor for every variable is not allowed, as that variable would predict its own missingness (MNAR)."),
      class = "efa_simulate_input"
    )
  }
  if (any(idx == seq_len(p))) {
    self <- as.character(which(idx == seq_len(p)))
    cli::cli_abort(
      c("{.arg missing_predictor} must reference another variable, not the variable itself.",
        "x" = "{cli::qty(self)}Variable{?s} {self} {?is/are} set to predict {?its/their} own missingness, which is MNAR rather than MAR.",
        "i" = "Use {.code missing = \"MNAR\"} for self-dependent missingness."),
      class = "efa_simulate_input"
    )
  }
  idx
}


# Solve for the logistic intercept a so the mean missing probability
# mean(plogis(a + b * z)) equals the target `prop`. The mean is strictly increasing in a
# from 0 to 1, so a unique root exists; the bracket qlogis(prop) +- (|b| max|z| + 10)
# forces a + b*z <= qlogis(prop) - 10 at the low end (mean below prop) and
# >= qlogis(prop) + 10 at the high end (mean above prop), so the ends are opposite-signed
# for any prop in (0, 1) and any finite b. uniroot draws no random numbers.
.calibrate_missing_intercept <- function(z, b, prop) {
  f <- function(a) mean(stats::plogis(a + b * z)) - prop
  # max|z| lives at an extremum, so range() avoids materializing a length-n abs(z).
  reach <- abs(b) * max(abs(range(z))) + 10
  centre <- stats::qlogis(prop)
  stats::uniroot(f, interval = c(centre - reach, centre + reach),
                 tol = .Machine$double.eps^0.5)$root
}


# Introduce missing values into a drawn data matrix under a chosen mechanism (Rubin,
# 1976). Each variable is holed at (in expectation) `prop` of its cases. MCAR draws a
# Bernoulli(prop) mask independent of the data. MAR and MNAR set each case's missing
# probability by a logistic model, logit(p_i) = a + b * z_i, where z standardizes a
# predictor: another variable (MAR, given by `pred_idx`) or the variable's own value
# (MNAR). The intercept a is calibrated per variable so the mean missing probability is
# `prop`; the slope b (`strength`) sets the dependence. Predictors are read from a single
# pre-loop copy of the complete draw so an earlier column's holes cannot corrupt a later
# column's predictor (R copy-on-modify leaves `complete` untouched once `dat` is written).
# A constant predictor (sd 0) contributes z = 0, degrading that variable to MCAR. Exactly
# one runif(n) is drawn per column in fixed order, never short-circuited, so the RNG
# stream -- and hence the reproducibility across worker counts -- does not depend on the
# mechanism or the achieved probabilities.
# Rubin, D. B. (1976). Inference and missing data. Biometrika, 63(3), 581-592.
.apply_missing <- function(dat, mechanism, prop, strength, pred_idx) {
  n <- nrow(dat)
  complete <- dat
  for (j in seq_len(ncol(dat))) {
    prob <- if (mechanism == "MCAR") {
      prop                                 # recycled against runif(n) below
    } else {
      xk <- complete[, if (mechanism == "MAR") pred_idx[j] else j]
      s <- stats::sd(xk)
      z <- if (s > 0) (xk - mean(xk)) / s else rep(0, n)
      a <- .calibrate_missing_intercept(z, strength, prop)
      stats::plogis(a + strength * z)
    }
    dat[stats::runif(n) < prob, j] <- NA
  }
  dat
}


# vech of a symmetric matrix: its lower triangle (including the diagonal) in column-major
# order (matching MBESS's vech and .duplication-weighted inner product used below).
.vech <- function(x) x[!upper.tri(x)]


# Perturb a population correlation matrix so a `q`-factor model fits it with a prescribed
# discrepancy, injecting "model error" (all models are approximations, so an exact-fit
# population overstates recovery; MacCallum, 2003). Recovers the model's unrotated orthogonal
# loadings from the common part (`common` = Lambda Phi Lambda', `Sigma_cov` the assembled
# covariance) -- the standardized common matrix `C_std` has the communalities on its diagonal,
# and its top-`q` eigenpairs give loadings `L` with `L L' + diag(psi) == R_pop` -- then hands
# them to the chosen method. Returns the perturbed correlation and the achieved fit -- the misfit
# of the specified q-factor (generating) model to the perturbed population, measured with the
# shared population-limit wrapper (.efa_population_fit) using EFA()'s fit-index formulas. For CB
# the generating model is the minimizer, so this equals the target; for TKL/WB it is the
# generating model's misfit (as in fungible/noisemaker), which a re-fitted EFA of the data may
# recover as somewhat smaller since it re-optimizes the loadings.
.apply_model_error <- function(R_pop, common, Sigma_cov, q, method,
                               target_rmsea, target_cfi, force_pd, tol = 1e-8) {
  p <- nrow(R_pop)
  h2 <- diag(common) / diag(Sigma_cov)     # standardized communalities
  u <- 1 - h2                              # uniquenesses
  C_std <- R_pop
  diag(C_std) <- h2
  ev <- eigen(C_std, symmetric = TRUE)
  L <- ev$vectors[, seq_len(q), drop = FALSE] %*%
    diag(sqrt(pmax(ev$values[seq_len(q)], 0)), q)
  psi <- 1 - rowSums(L^2)

  # The recovered L reproduces R_pop only when the population is exactly a rank-q common factor
  # structure with diagonal, non-negative uniquenesses. A correlated-residual (off-diagonal) Psi
  # or a communality above 1 breaks that: the eigen-recovered model would misrepresent the
  # population, so CB would centre on the wrong model (missing its exact-target guarantee), TKL
  # would take the square root of a negative uniqueness, and the reported achieved fit would be
  # measured against the wrong model. Reject such populations rather than silently mislead.
  if (any(psi < -tol) ||
      max(abs(tcrossprod(L) + diag(psi, p) - R_pop)) > 1e-8) {
    cli::cli_abort(
      c("Model error requires a population that is exactly a {q}-factor structure with non-negative uniquenesses.",
        "x" = "The population has correlated residuals (an off-diagonal {.arg Psi}) or a communality above 1, which the {q}-factor model cannot represent.",
        "i" = "Use a factor model with a diagonal {.arg Psi} and communalities below 1, or omit the model-error target."),
      class = "efa_simulate_model_error"
    )
  }

  me <- switch(
    method,
    CB  = .model_error_cb(R_pop, L, q, target_rmsea, force_pd),
    TKL = .model_error_tkl(R_pop, u, q, target_rmsea, target_cfi),
    WB  = .model_error_wb(R_pop, target_rmsea)
  )

  # Achieved fit of the major q-factor model on the perturbed population, on the same
  # (population-limit) scale a fit would report. For CB this returns the target to machine
  # precision (the model is the minimizer by construction); for TKL/WB it is the realized value.
  fit <- .efa_population_fit(L, me$population)
  model_error <- c(
    list(method = method, target_rmsea = target_rmsea, target_cfi = target_cfi,
         rmsea = fit$rmsea, cfi = fit$cfi, df = fit$df),
    me$params
  )
  list(population = me$population, model_error = model_error)
}


# Cudeck & Browne (1992) model error: construct Sigma* = Omega + kappa E such that the q-factor
# model (unrotated loadings `L`; its uniquenesses follow from `R_pop`) remains the minimizer of the
# ML discrepancy and the minimum equals delta = df * target_rmsea^2, so the achieved RMSEA equals
# the target.
# Follows MBESS::Sigma.2.SigmaStar at the correlation scale: differentiate the model-implied
# Sigma(theta) = Lt Lt' + diag(psi) with respect to theta = (vec(L), psi), form the ML-metric
# Jacobian columns B[,i] = -d * vech(W^-1 dSigma_i W^-1) (W = Omega, d the vech weights
# that make an ordinary residual orthogonal in the ML metric), and take E from the part of a
# random symmetric matrix orthogonal to B. G = W^-1 E has the same eigenvalues s as the symmetric
# W^{-1/2} E W^{-1/2}, so the discrepancy is kappa*sum(s) - sum(log(1 + kappa*s)); it is monotone
# in kappa up to the positive-definite boundary 1 + kappa*min(s) > 0. kappa is the exact root
# (uniroot on the monotone equation, unlike MBESS's nlm on |.|, which is less precise near the
# kink and would miss the target here). A direction that cannot reach delta while staying positive
# definite is redrawn; if none succeeds, force_pd = TRUE returns the closest positive-definite
# approximation (a smaller RMSEA, with a warning), otherwise the call aborts.
# Cudeck, R., & Browne, M. W. (1992). Constructing a covariance matrix that yields a specified
# minimizer and a specified minimum discrepancy function value. Psychometrika, 57(3), 357-369.
.model_error_cb <- function(R_pop, L, q, target_rmsea, force_pd,
                            max_tries = 50L) {
  p <- nrow(R_pop)
  df <- .efa_df(p, q)
  delta <- target_rmsea^2 * df
  pstar <- p * (p + 1) / 2

  W <- R_pop
  Wi <- solve(W)
  # vech inner-product weights (1 on the diagonal, 2 off it): with these, an ordinary residual
  # against B is orthogonal in the ML metric, which is what keeps theta the minimizer.
  wmat <- matrix(2, p, p)
  diag(wmat) <- 1
  dvec <- .vech(wmat)

  # The Jacobian columns in closed form. Elementwise Sigma_ij(theta) = sum_k L_ik L_jk +
  # [i == j] psi_i, so
  #   dSigma / dL[a, b] = e_a t(L[, b]) + L[, b] t(e_a)     -> column (b - 1) * p + a
  #   dSigma / dpsi[a]  = e_a t(e_a)                        -> column p * q + a,
  # ordered to match theta = (vec(L), psi), which R fills column-major. Each ML sandwich is then
  # a sum of outer products,
  #   Wi (e_a t(y) + y t(e_a)) Wi = Wi[, a] t(t(Wi) y) + (Wi y) t(Wi[a, ]),
  # whose vech picks x[ii] * y[jj] over the lower-triangle index pairs, so no p x p product is
  # formed per column. Wi[, a] and Wi[a, ] are kept apart because solve() need not return an
  # exactly symmetric inverse.
  #
  # Differentiating in closed form rather than by a finite difference is what keeps the draw
  # below reproducible: differencing divides a cancelling subtraction by the step size, which
  # lifts the q(q - 1)/2 rotational null directions of the model off zero to roughly eps/step.
  # Their orientation is then set purely by rounding, qr() counts them as real columns, and
  # qr.resid() projects them out -- so E, and with it the perturbed population, would swing by
  # tens of percent between BLAS implementations for one and the same seed.
  idx <- which(!upper.tri(matrix(0, p, p)), arr.ind = TRUE)   # the vech pairs, in .vech order
  ii <- idx[, 1L]
  jj <- idx[, 2L]
  WL <- Wi %*% L            # columns: Wi %*% L[, b]
  WtL <- crossprod(Wi, L)   # columns: t(Wi) %*% L[, b]
  npar <- p * q + p
  B <- matrix(0, pstar, npar)
  for (a in seq_len(p)) {
    wi_i <- Wi[ii, a]       # (Wi e_a)[ii]
    wi_j <- Wi[a, jj]       # (t(e_a) Wi)[jj]
    for (b in seq_len(q)) {
      B[, (b - 1L) * p + a] <- -dvec * (wi_i * WtL[jj, b] + WL[ii, b] * wi_j)
    }
    B[, p * q + a] <- -dvec * (wi_i * wi_j)
  }
  Bqr <- qr(B)

  # symmetric whitening W^{-1/2} for the eigenvalues s that drive both the discrepancy and the
  # positive-definite bound on kappa.
  Weig <- eigen(W, symmetric = TRUE)
  Wih <- Weig$vectors %*% diag(1 / sqrt(Weig$values), p) %*% t(Weig$vectors)

  best_pop <- NULL       # closest positive-definite approximation, for the force_pd escape
  best_gap <- -Inf       # its (negative) shortfall in discrepancy relative to delta
  for (tr in seq_len(max_tries)) {
    e <- qr.resid(Bqr, stats::rnorm(pstar))
    E <- matrix(0, p, p)
    E[lower.tri(E, diag = TRUE)] <- e
    E <- E + t(E) - diag(diag(E))
    s <- eigen(Wih %*% E %*% Wih, symmetric = TRUE)$values
    trG <- sum(s)
    fdisc <- function(k) k * trG - sum(log1p(k * s)) - delta
    hi <- if (min(s) < 0) -1 / min(s) * (1 - 1e-9) else 1e6
    gap <- fdisc(hi)                          # discrepancy reachable at the PD boundary, minus delta
    if (gap >= 0) {
      kappa <- stats::uniroot(fdisc, c(1e-12, hi),
                              tol = .Machine$double.eps^0.5)$root
      return(list(population = stats::cov2cor(R_pop + kappa * E),
                  params = list(kappa = kappa)))
    }
    if (gap > best_gap) {                      # track the direction that comes closest to delta
      best_gap <- gap
      best_pop <- stats::cov2cor(R_pop + hi * E)
    }
  }

  # No direction reached delta while staying positive definite: the target is too large for this
  # population. Either accept the closest positive-definite approximation or abort.
  if (isTRUE(force_pd) && !is.null(best_pop)) {
    cli::cli_warn(
      c("The requested {.arg target_rmsea} could not be reached while keeping the population positive definite; the closest positive-definite approximation was used.",
        "i" = "The achieved RMSEA is below the target; lower {.arg target_rmsea} or use {.code model_error = \"TKL\"} for a larger misfit."),
      class = "efa_simulate_pd_forced"
    )
    return(list(population = best_pop, params = list(kappa = NA_real_)))
  }
  cli::cli_abort(
    c("The Cudeck-Browne model error cannot reach {.arg target_rmsea} = {target_rmsea} while keeping the population positive definite.",
      "i" = "Lower {.arg target_rmsea}, use {.code model_error = \"TKL\"} or {.val WB}, or set {.code force_pd = TRUE} to accept the closest positive-definite approximation."),
    class = "efa_simulate_model_error"
  )
}


# Tucker, Koopman & Linn (1969) model error: add a set of `n_minor` minor common factors to the
# population so the q-factor major model misfits it by a target RMSEA (and, optionally, CFI). The
# minor loadings `W0` are standard normal; a geometric column decay (1 - eps)^(0:(n_minor-1))
# concentrates variance in the first minor factors, and each row is rescaled so the minor common
# variance of variable i equals v * u_i (v is the proportion of its uniqueness u_i reallocated to
# minor factors). The perturbed population is R_pop + W W' (diagonal reset to 1). The two knobs
# (v, eps) are chosen by bounded optimization (L-BFGS-B with random restarts) to minimise the
# weighted relative squared deviation from the target(s); with a single target the fit is near
# exact, with both targets it is a compromise. Ft (major-model discrepancy) and Fb (independence
# baseline) use the same ML discrepancy as .gof()/.efa_population_fit().
# Tucker, L. R., Koopman, R. F., & Linn, R. L. (1969). Evaluation of factor analytic research
# procedures by means of simulated correlation matrices. Psychometrika, 34(4), 421-459.
.model_error_tkl <- function(R_pop, u, q, target_rmsea, target_cfi,
                             n_minor = 50L, max_tries = 25L) {
  p <- nrow(R_pop)
  df <- .efa_df(p, q)
  W0 <- matrix(stats::rnorm(p * n_minor), p, n_minor)
  Rinv <- solve(R_pop)
  ldRpop <- as.numeric(determinant(R_pop, logarithm = TRUE)$modulus)
  has_r <- !is.null(target_rmsea)
  has_c <- !is.null(target_cfi)

  obj <- function(par, return_vals = FALSE) {
    v <- par[1L]
    eps <- par[2L]
    # Column j of the minor loadings decays by (1 - eps)^(j-1); scale columns directly rather
    # than multiplying by a dense n_minor x n_minor diagonal on every objective evaluation.
    Wc <- sweep(W0, 2L, (1 - eps)^(0:(n_minor - 1L)), "*")
    wsq <- rowSums(Wc^2)
    Wc <- sqrt(v * u / wsq) * Wc          # scale row i so rowSums(Wc^2) = v * u_i
    RpopME <- R_pop + tcrossprod(Wc)
    diag(RpopME) <- 1
    ldME <- as.numeric(determinant(RpopME, logarithm = TRUE)$modulus)
    Ft <- ldRpop - ldME + sum(Rinv * RpopME) - p     # ML discrepancy of the major model
    Fb <- -ldME                                      # independence baseline
    rmsea <- sqrt(max(0, Ft) / df)
    cfi <- 1 - Ft / Fb
    if (return_vals) {
      return(list(population = RpopME, v = v, eps = eps))
    }
    val <- 0
    if (has_r) val <- val + (rmsea - target_rmsea)^2 / target_rmsea^2
    if (has_c) val <- val + ((1 - cfi) - (1 - target_cfi))^2 / (1 - target_cfi)^2
    if (!is.finite(val)) 1e6 else val
  }

  best_par <- NULL
  best_val <- Inf
  for (tr in seq_len(max_tries)) {
    start <- c(stats::runif(1, 0.02, 0.9), stats::runif(1, 0, 0.8))
    opt <- tryCatch(
      stats::optim(start, obj, method = "L-BFGS-B",
                   lower = c(0.001, 0), upper = c(1, 1)),
      error = function(e) NULL)
    if (is.null(opt)) next
    if (opt$value < best_val) {
      best_val <- opt$value
      best_par <- opt$par
    }
    if (opt$convergence == 0 && opt$value < 1e-8) break   # target(s) essentially hit
  }
  if (is.null(best_par)) {
    cli::cli_abort(
      c("The Tucker-Koopman-Linn optimization failed to find minor-factor weights for the requested target(s).",
        "i" = "Try a different {.arg target_rmsea}/{.arg target_cfi}, or {.code model_error = \"CB\"}."),
      class = "efa_simulate_model_error")
  }
  res <- obj(best_par, return_vals = TRUE)
  list(population = res$population, params = list(v = res$v, eps = res$eps))
}


# Wu & Browne (2015) model error: draw the population from an inverse-Wishart distribution
# centered on the model-implied correlation Omega, Sigma ~ IW(m, m Omega), where the precision
# m = 1 / target_rmsea^2 sets the expected RMSEA to approximately the target (via the Wu-Browne
# small-error relation v ~= epsilon^2; the realized RMSEA is biased slightly upward, more so for
# larger targets). The drawn covariance is
# standardized to a correlation matrix; the achieved RMSEA is random around the target (this is a
# single random draw, not an exact construction), so the realized value is reported. Requires
# m >= p for a proper draw (a very large target RMSEA is rejected).
# Wu, H., & Browne, M. W. (2015). Quantifying adventitious error in a covariance structure as a
# random effect. Psychometrika, 80(3), 571-600.
.model_error_wb <- function(R_pop, target_rmsea) {
  p <- nrow(R_pop)
  m <- 1 / target_rmsea^2
  if (m < p) {
    cli::cli_abort(
      c("{.arg target_rmsea} = {target_rmsea} is too large for the Wu-Browne method with {p} variable{?s}.",
        "i" = "The precision m = 1 / target_rmsea^2 = {round(m, 1)} must be at least the number of variables ({p}); use a smaller {.arg target_rmsea}."),
      class = "efa_simulate_model_error")
  }
  Sigma <- .riwish(m, m * R_pop)
  list(population = stats::cov2cor(Sigma), params = list(m = m))
}


# One inverse-Wishart draw Sigma ~ IW(nu, S) via the Bartlett decomposition (no external
# dependency): S^-1 = C C', and A A' ~ Wishart(nu, S^-1) with A = C T, T lower-triangular whose
# diagonal is sqrt(chi-square(nu - i + 1)) and whose strictly-lower entries are standard normal;
# Sigma = (A A')^-1. Draws its randomness (rchisq, rnorm) from the caller's stream, so it is
# reproducible under a fixed seed.
.riwish <- function(nu, S) {
  p <- nrow(S)
  C <- t(chol(solve(S)))
  A <- matrix(0, p, p)
  diag(A) <- sqrt(stats::rchisq(p, df = nu - seq_len(p) + 1))
  A[lower.tri(A)] <- stats::rnorm(p * (p - 1) / 2)
  CA <- C %*% A
  Sig <- solve(tcrossprod(CA))
  (Sig + t(Sig)) / 2
}
