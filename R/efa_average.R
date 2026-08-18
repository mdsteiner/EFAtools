#' Model averaging across different EFA estimators and types
#'
#' Not all EFA procedures always arrive at the same solution. This function allows
#' you perform a number of EFAs from different estimators (e.g., Maximum Likelihood
#' and Principal Axis Factoring), with different implementations (e.g., the SPSS
#' and psych implementations of Principal Axis Factoring), and across different
#' rotations of the same type (e.g., multiple oblique rotations, like promax and
#' oblimin). `efa_average()` will then run all these EFAs (using the [efa_fit()]
#' function) and provide a summary across the different solutions.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations. If raw data is entered, the correlation matrix is found from the
#' data.
#' @param n_factors numeric. Number of factors to extract.
#' @param N numeric. The number of observations. Needs only be specified if a
#' correlation matrix is used. If input is a correlation matrix and `N` = NA
#' (default), not all fit indices can be computed.
#' @param estimator character vector. Any combination of  "PAF", "ML", and "ULS",
#' to use principal axis factoring, maximum likelihood, or unweighted least
#' squares, respectively, to fit the EFAs. "MINRES" is accepted as a synonym for
#' "ULS" (the same estimator). Default is "PAF".
#' "DWLS", which [efa_fit()] does accept, is deliberately not
#' offered here: it weights each residual correlation by the inverse of its
#' asymptotic variance, which is only available from raw ordinal data analysed
#' with `cor_method = "poly"` or `"tetra"`, whereas every EFA in the grid is fitted
#' to the single correlation matrix computed once from `x`. Fit a DWLS solution
#' with [efa_fit()] directly.
#' @param rotation character vector. Either perform no rotation ("none"),
#' any combination of orthogonal rotations ("varimax", "equamax", "quartimax", "geominT",
#' "bentlerT", and "bifactorT"; using "orthogonal" runs all of these), or of
#' oblique rotations ("promax", "oblimin", "quartimin", "simplimax", "bentlerQ",
#' "geominQ", and "bifactorQ"; using "oblique" runs all of these). Rotation types
#' (no rotation, orthogonal rotations, and oblique rotations) cannot be mixed.
#' Default is "promax".
#' @param type character vector. Any combination of "none" (default), "EFAtools",
#' "psych", and "SPSS" can be entered. "none" allows the specification of various
#' combinations of the arguments controlling both factor extraction methods and
#' the rotations. The others ("EFAtools", "psych", and "SPSS") take the extraction
#' and rotation tuning of the respective implementation: this package's default
#' procedure, the psych package's, and SPSS's. A specific psych implementation
#' exists for PAF, ML, varimax, and promax. The SPSS implementation exists for
#' PAF, varimax, and promax. For details, see [efa_fit()]. The factor ordering is
#' the one setting a named `type` does not bring here: every solution in the grid
#' is fitted with the eigenvalue-based ordering, so that the solutions can be
#' realigned to a common target before averaging.
#' @param averaging character. One of "mean" (default), and "median". Controls
#' whether the different results should be averaged using the (trimmed) mean,
#' or the median.
#' @param trim numeric. If averaging is set to "mean", this argument controls
#' the trimming of extremes (for details see [base::mean()]).
#' By default no trimming is done (i.e., trim = 0).
#' @param salience_threshold numeric. The threshold to use to classify a pattern
#' coefficient or loading as salient (i.e., substantial enough to assign it to
#' a factor). Default is 0.3. Indicator-to-factor correspondences will be inferred
#' based on this threshold. Note that this may not be meaningful if rotation = "none"
#' and n_factors > 1 are used, as no simple structure is present there.
#' @param max_iter numeric. The maximum number of iterations to perform after which
#' the iterative PAF procedure is halted with a warning. Default is 10,000. It is
#' only evaluated for the "PAF" solutions run under `type` "none": a named `type`
#' brings the iteration cap that defines it ("SPSS" 25, "psych" 50, and
#' "EFAtools" 300), and "ML" and "ULS" do not iterate this way. Note
#' that non-converged procedures are excluded from the averaging procedure.
#' @param init_comm character vector. Any combination of "smc", "mac", and "unity".
#' Controls the methods to estimate the initial communalities in `PAF` if
#' "none" is among the specified types. "smc" will use squared multiple
#' correlations, "mac" will use maximum absolute correlations, "unity" will use
#' 1s (for details see [efa_fit()]). Default is `c("smc", "mac", "unity")`.
#' @param criterion numeric vector. The convergence criterion used for PAF if
#' "none" is among the specified types.
#' If the change in communalities from one iteration to the next is smaller than
#' this criterion the solution is accepted and the procedure ends.
#' Default is `0.001`.
#' @param criterion_type character vector. Any combination of "max_individual" and
#' "sum". Type of convergence criterion used for PAF if "none" is among the
#' specified types. "max_individual" selects the maximum change in any of the
#' communalities from one iteration to the next and tests it against the
#' specified criterion. "sum" takes the difference of
#' the sum of all communalities in one iteration and the sum of all communalities
#' in the next iteration and tests this against the criterion
#' (for details see [efa_fit()]). Default is `c("sum", "max_individual")`.
#' @param abs_eigen logical vector. Any combination of TRUE and FALSE.
#' Which algorithm to use in the PAF iterations if "none" is among the specified
#' types. If FALSE, the loadings are computed from the eigenvalues. This is also
#' used by the [psych::fa()] function. If TRUE the
#' loadings are computed with the absolute eigenvalues as done by SPSS
#' (for details see [efa_fit()]). Default is `TRUE`.
#' @param varimax_type character vector. Any combination of "svd" and "kaiser".
#' The type of the varimax rotation performed if "none" is among the specified
#' types and "varimax", "promax", "orthogonal", or "oblique" is among the specified
#' rotations. "svd" uses singular value decomposition, as
#' [stats::varimax()] does, and "kaiser" uses the varimax
#' procedure performed in SPSS. This is the original procedure from Kaiser (1958),
#' but with slight alterations in the varimax criterion (for details, see
#' [efa_fit()] and Grieder & Steiner, 2022).
#' Default is `c("svd", "kaiser")`.
#' @param normalize logical vector. Any combination of TRUE and FALSE.
#' `TRUE` performs a kaiser normalization before the specified rotation(s).
#' Default is `TRUE`.
#' @param k_promax numeric vector. The power used for computing the target matrix
#' P in the promax rotation if "none" is among the specified types and "promax"
#' or "oblique" is among the specified rotations. Default is `2:4`.
#' @param k_simplimax numeric. The number of 'close to zero loadings' for the
#' simplimax rotation if "simplimax" or "oblique" is among the specified rotations. Default
#' is `ncol(x)`, where x is the entered data. It counts loadings, so each value must be a
#' whole number no larger than the number of loadings in the solution; a simplimax fit
#' given anything else fails and is reported as an errored solution in the grid.
#' @param p_type character vector. Any combination of "norm" and "unnorm".
#' This specifies how the target matrix P is computed in promax rotation if
#' "none" is among the specified types and "promax" or "oblique" is among the
#' specified rotations. "unnorm" will use the unnormalized target matrix as
#' originally done in Hendrickson and White (1964). "norm" will use a
#' normalized target matrix (for details see [efa_fit()]).
#' Default is `c("norm", "unnorm")`.
#' @param precision numeric vector. The tolerance for stopping in the rotation
#' procedure(s). Default is 10^-5.
#' @param start_method character vector. Any combination of "psych" and "factanal".
#' How to specify the starting values for the optimization procedure for ML.
#' "psych" takes the starting values specified in [psych::fa()].
#' "factanal" takes the starting values specified in the
#' [stats::factanal()] function. Default is
#' `c("psych", "factanal")`.
#' @param use character. Passed to [stats::cor()] if raw data
#' is given as input. Default is "pairwise.complete.obs". It is ignored when
#' `cor_method = "fiml"`, which handles the missingness itself, so every case
#' contributes.
#' @param cor_method character. Correlation computed from raw data: `"pearson"`,
#'   `"spearman"`, or `"kendall"` (passed to [stats::cor()]), `"poly"` /
#'   `"tetra"` for polychoric / tetrachoric correlations of ordinal / binary data
#'   (a two-step estimator), or `"fiml"`
#'   for a two-stage full-information maximum-likelihood correlation from raw data
#'   with missing values. With `"fiml"` the saturated multivariate-normal mean and
#'   covariance are estimated by an EM algorithm assuming the data are missing at
#'   random and the standardized covariance is analysed, reproducing
#'   `psych::corFiml()` followed by `psych::fa()` and
#'   `lavaan(missing = "two.stage")`, *not* `lavaan::efa(missing = "ml")` (see
#'   [efa_fit()] and the details). Default is "pearson".
#' @param show_progress logical. Whether a progress bar should be shown in the
#' console. Default is TRUE.
#' @param seed numeric or `NULL`. An optional seed making the averaging run
#' reproducible and independent of the number of parallel workers (the grid runs with
#' [future_lapply()][future.apply::future_lapply], for which a parallel plan can be set
#' via [future::plan()]). It matters whenever a criterion-based rotation is included,
#' since those draw random starts. When supplied, the caller's random-number stream is
#' restored afterwards, leaving no side effect. Default is `NULL`.
#' @param P_type `r lifecycle::badge("superseded")` Former name of `p_type`. Still
#' accepted (silently) for backwards compatibility; please use `p_type`.
#'
#' @details
#'
#' As a first step in this function, a grid is produced containing the setting
#' combinations for the to-be-performed EFAs. These settings are then entered as
#' arguments to the [efa_fit()] function and the EFAs are run in a second
#' step. After all EFAs are run, the factor solutions are averaged and their
#' variability determined in a third step.
#'
#' When raw data are supplied, the correlation matrix is computed once before the
#' grid is run and reused for every EFA in it. Under `cor_method = "fiml"` this
#' means the saturated multivariate-normal moments are EM-estimated a single time
#' (from the raw data with missing values, assuming the data are missing at
#' random) and the resulting two-stage correlation is analysed by every solution
#' in the grid; the EM is not re-run per solution. Under `cor_method = "fiml"`,
#' `use` does not select cases (every case contributes to the EM). The averaged
#' loadings and communalities are then the two-stage FIML estimates, but the
#' averaged Chi-Square and the indices derived from it (CFI, TLI, RMSEA, AIC, BIC,
#' ECVI) are the ordinary ML/ULS discrepancy statistics on the EM correlation, not
#' the corrected two-stage (Satorra-Bentler) statistics that a standalone
#' [efa_fit()] with `cor_method = "fiml"` reports; in particular the averaged AIC and
#' BIC are finite here rather than `NA`.
#'
#' The grid containing the setting combinations is produced based on the entries
#' to the respective arguments. To this end, all possible combinations resulting
#' in unique EFA models are considered: combinations that resolve to the same
#' model are run only once. Two combinations are the same model only if every
#' setting the fit consumes agrees, and for "PAF" that includes the iteration cap
#' `max_iter`. Since a named `type` brings its own cap, a `type` of
#' `c("none", "SPSS")` whose specific settings match the SPSS combination in
#' every other respect still gives two "PAF" models, unless `max_iter` is also
#' set to the cap of the SPSS implementation. We include here a list
#' of arguments that are only evaluated under specific conditions:
#'
#' The arguments `init_comm`, `criterion`, `criterion_type`,
#' `abs_eigen`, and `max_iter` are only evaluated if "PAF" is included in
#' `estimator` and "none" is included in `type`.
#'
#' The argument `varimax_type` is only evaluated if "varimax", "promax",
#' "oblique", or "orthogonal" is included in `rotation` and "none" is
#' included in `type`.
#'
#' The argument `normalize` is only evaluated if `rotation` is not
#' set to "none" and "none" is included in `type`.
#'
#' The argument `k_simplimax` is only evaluated if "simplimax" or "oblique"
#' is included in `rotation`.
#'
#' The arguments `k_promax` and `p_type` are only evaluated if
#' "promax" or "oblique" is included in `rotation` and "none" is included
#' in `type`.
#'
#' The argument `start_method` is only evaluated if "ML" is included in
#' `estimator`.
#'
#' Every solution in the grid is fitted with the eigenvalue-based factor ordering,
#' including under a named `type`: the solutions are realigned to a common target
#' before averaging, so a per-fit ordering (SPSS orders by the sum of squared
#' loadings) would not survive into the averaged result. That target is the first
#' solution the grid retains, in grid order; solutions dropped for an error,
#' non-convergence, or a Heywood case cannot become it. The averaged loadings are
#' therefore in the factor order and sign of that solution. This is the one setting a
#' named `type` does not carry, and it is visible in the two places the individual
#' fits are: the solutions returned in `efa_list` are eigenvalue-ordered, and so is
#' the single [efa_fit()] object returned when the grid collapses to one row. Their
#' loadings can therefore appear in a different factor order than the same fit run
#' through [efa_fit()] under that `type`, even though the solutions are the same.
#'
#' To avoid a bias in the averaged factor solutions from problematic solutions,
#' these are excluded prior to averaging. A solution is deemed problematic if
#' at least one of the following is true: an error occurred, the model did not
#' converge, or there is at least one Heywood (improper) case (a communality at
#' or above 1, or, for `ML`/`ULS`, a uniqueness pinned at the estimator's lower
#' bound).
#' Information on errors, convergence, and Heywood cases are returned in the
#' implementations_grid and a summary of these is given when printing the output.
#' In addition to these, information on the admissibility of the factor solutions
#' is also included. A solution was deemed admissible if (1) no error occurred,
#' (2) the model converged, (3) no Heywood cases are present, and (4) there are
#' at least two salient loadings (i.e., loadings exceeding the specified
#' `salience_threshold`) for each factor. So, solutions failing one of the
#' first three of these criteria of admissibility are also deemed problematic and
#' therefore excluded from averaging. However, solutions failing only
#' the fourth criterion of admissibility are still included for averaging.
#' Finally, if all solutions are problematic (e.g., all solutions contain
#' Heywood cases), no averaging is performed and the respective outputs are NA.
#' In this case, the implementations_grid should be inspected to see if there
#' are any error messages, and the separate EFA solutions that are also included
#' in the output can be inspected as well, for example, to see where Heywood
#' cases occurred.
#'
#' A core output of this function includes the average, minimum, and maximum
#' loadings derived from all non-problematic (see above) factor solutions. Please
#' note that these are not entire solutions, but the matrices include the average,
#' minimum, or maximum value for each cell (i.e., each loading separately). This
#' means that, for example, the matrix with the minimum loadings will contain
#' the minimum value in any of the factor solutions for each specific loading,
#' and therefore most likely contains loadings from different factor solutions.
#' The matrices containing the minimum and maximum factor solutions can
#' therefore not be interpreted as whole factor solutions.
#'
#' The averaged loading matrix is likewise a cell-wise summary rather than a fitted
#' solution: it is not itself the solution of any EFA, does not in general reproduce
#' the correlation matrix, and need not reproduce the averaged communalities. The
#' fit indices described below are correspondingly the mean (or, under
#' `averaging = "median"`, the median) of the per-solution fit indices, not the fit
#' of the averaged loadings, so the averaged loadings and the reported fit do not
#' describe one and the same model.
#'
#' The output also includes information on the average, minimum, maximum, and
#' variability of the fit indices across the non-problematic factor solutions.
#' It is important to note that not all fit indices are computed for all fit
#' methods: For ML and ULS, all fit indices can be computed, while for PAF the
#' chi-square-based indices (the chi-square statistic and its significance, CFI,
#' TLI, RMSEA, AIC, BIC, and ECVI) are NA. The common part accounted for (CAF)
#' index (Lorenzo-Seva, Timmerman, & Kiers, 2011) and the residual-based SRMR and
#' RMSR are still computed for PAF. As a consequence, if only "PAF" is included in
#' the `estimator` argument, averaging is performed for the CAF, SRMR, and RMSR, while
#' the chi-square-based indices are NA. If a combination of "PAF" and "ML" and/or
#' "ULS" are included in the `estimator` argument, the CAF, SRMR, and RMSR are
#' averaged across all non-problematic factor solutions, while the chi-square-based
#' indices are only averaged across the ML and ULS solutions. The user should
#' therefore keep in mind that the number of EFAs across which the fit indices are
#' averaged can diverge for the CAF, SRMR, and RMSR compared to the chi-square-based
#' indices.
#'
#' Each reported fit index is summarised across the (non-problematic) solutions in
#' the same descriptive way: the average, standard deviation, minimum, and maximum
#' of the per-solution values. This includes the chi-square significance level
#' (`p_chi`), which is therefore the mean (or median) of the per-solution p-values
#' and is purely descriptive; it is not the p-value of any pooled chi-square test.
#'
#' @return A list of class `c("efa_average", "EFA_AVERAGE")` containing the
#' components below. Throughout, `range` is the width `maximum - minimum` of each
#' cell across the factor solutions, not the interval `[minimum, maximum]`, and
#' `average` is the (trimmed) mean or the median, following `averaging`.
#' \item{orig_R}{Original correlation matrix.}
#' \item{h2}{A list with the average, standard deviation, minimum, maximum, and
#' range of the final communality estimates across the factor solutions.}
#' \item{loadings}{A list with the average, standard deviation, minimum, maximum,
#' and range of the final loadings across the factor solutions. If rotation was
#' "none", the unrotated loadings, otherwise the rotated loadings (pattern
#' coefficients).}
#' \item{Phi}{A list with the average, standard deviation, minimum, maximum, and
#' range of the factor intercorrelations across factor solutions obtained with
#' oblique rotations.}
#' \item{ind_fac_corres}{A matrix with each cell containing the proportion of
#' the factor solutions in which the respective indicator-to-factor correspondence
#' occurred, i.e., in which the loading exceeded the specified salience threshold.
#' Note: Rowsums can exceed 1 due to cross-loadings.}
#' \item{vars_accounted}{A list with the average, standard deviation, minimum,
#' maximum, and range of explained variances and sums of squared loadings across
#' the factor solutions. Based on the unrotated loadings if rotation was "none"
#' or only one factor was extracted, otherwise on the rotated loadings. Each entry
#' is a matrix with rows `"SS loadings"`, `"Prop Tot Var"`, and `"Prop Comm Var"`;
#' the last is omitted for a single-factor solution, where it is identically 1, so
#' those matrices have two rows there and three otherwise.}
#' \item{fit_indices}{A matrix containing the average, standard deviation,
#' minimum, maximum, and range for all applicable fit indices across the respective
#' factor solutions, and the degrees of freedom (df). If the estimator argument
#' contains ML or ULS: Fit indices derived
#' from the unrotated factor loadings: Chi Square (chisq), including significance
#' level, Comparative Fit Index (CFI), Tucker-Lewis Index (TLI), Root Mean Square
#' Error of Approximation (RMSEA), Akaike Information Criterion (AIC), Bayesian
#' Information Criterion (BIC), Expected Cross-Validation Index (ECVI), and the
#' common part accounted for (CAF) index as proposed by Lorenzo-Seva, Timmerman,
#' & Kiers (2011). The residual-based Standardized Root Mean Square Residual
#' (SRMR) and Root Mean Square Residual (RMSR) and the CAF are also computed for
#' PAF; for PAF the remaining (chi-square-based) indices are not available (see
#' details).}
#' \item{implementations_grid}{A matrix containing, for each performed EFA,
#' the setting combination, if an error occurred (logical), the error message
#' (character), an integer convergence code (0 = converged; for ML and ULS the
#' same codes as [`stats::optim()`][stats::optim]'s "L-BFGS-B", for PAF 1 if the
#' maximum number of iterations was reached without meeting the convergence
#' criterion and 0 otherwise),
#' if heywood cases occurred (logical, see details for definition), if the
#' solution was admissible (logical, see details for definition), and the fit
#' indices.}
#' \item{efa_list}{A list containing the outputs of all performed EFAs. The names
#' correspond to the rownames from the implementations_grid.}
#' \item{settings}{A list of the settings used, including `seed` (`NULL` when none was
#'   supplied).}
#'
#' If the supplied arguments admit only a single EFA, there is nothing to average
#' across: that one [efa_fit()] object is returned instead, with a warning. Its
#' `settings` are that fit's, with `seed` recording the seed the run was governed by.
#' They therefore describe the concrete arguments the fit ran under rather than the
#' `type` that supplied them: a row taken from a named preset records `type = "none"`
#' together with the preset's resolved values (for example `max_iter = 25` for
#' "SPSS"), and `order_type = "eigen"` as for every other row in the grid.
#'
#' @source Grieder, S., & Steiner, M. D. (2022). Algorithmic jingle jungle: A comparison
#' of implementations of principal axis factoring and promax rotation in R and SPSS.
#' Behavior Research Methods, 54, 54–74. doi: 10.3758/s13428-021-01581-x
#' @source Hendrickson, A. E., & White, P. O. (1964). Promax: A quick method for
#' rotation to oblique simple structure. British Journal of Statistical Psychology,
#' 17 , 65–70. doi: 10.1111/j.2044-8317.1964.tb00244.x
#' @source Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. L. (2011). The
#' Hull Method for Selecting the Number of Common Factors, Multivariate Behavioral
#' Research, 46, 340-364, doi: 10.1080/00273171.2011.564527
#' @source Kaiser, H. F. (1958). The varimax criterion for analytic rotation in
#' factor analysis. Psychometrika, 23, 187–200. doi: 10.1007/BF02289233
#'
#' @family factor analysis
#'
#' @export
#'
#' @examples
#' # Averaging across one implementation each of PAF (EFAtools type), ULS, and
#' # ML with one implementation of promax (EFAtools type) (3 EFAs)
#' Aver_meth <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                          estimator = c("PAF", "ULS", "ML"), type = "EFAtools",
#'                          start_method = "psych")
#'
#' \donttest{
#' # Averaging across different implementations of PAF and promax rotation (72 EFAs)
#' Aver_PAF <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500)
#'
#' # Use median instead of mean for averaging (72 EFAs)
#' Aver_PAF_md <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                            averaging = "median")
#'
#' # Averaging across different implementations of PAF and promax rotation,
#' # and across ULS and different versions of ML (108 EFAs)
#' Aver_meth_ext <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                              estimator = c("PAF", "ULS", "ML"))
#'
#' # Averaging across different oblique rotation methods, using one implementation
#' # of ML and one implementation of promax (EFAtools type) (7 EFAs)
#' Aver_rot <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                          estimator = "ML", rotation = "oblique", type = "EFAtools",
#'                          start_method = "psych")
#'}
#'
#' \donttest{
#' # Two-stage FIML correlations from raw data with missing values: the EM
#' # saturated moments are estimated once and the resulting correlation is
#' # averaged across the grid of EFAs.
#' x_miss <- GRiPS_raw
#' x_miss[cbind(1:20, 1)] <- NA
#' Aver_fiml <- efa_average(x_miss, n_factors = 1, estimator = c("PAF", "ML"),
#'                          cor_method = "fiml")
#' }
#'
efa_average <- function(x, n_factors, N = NA, estimator = "PAF", rotation = "promax",
                        type = "none", averaging = c("mean", "median"), trim = 0,
                        salience_threshold = .3,
                        max_iter = 1e4,
                        init_comm = c("smc", "mac", "unity"),
                        criterion = c(1e-3),
                        criterion_type = c("sum", "max_individual"),
                        abs_eigen = c(TRUE),
                        varimax_type = c("svd", "kaiser"),
                        normalize = TRUE,
                        k_promax = 2:4, k_simplimax = ncol(x),
                        p_type = c("norm", "unnorm"), precision = 1e-5,
                        start_method = c("psych", "factanal"),
                        use = c("pairwise.complete.obs", "all.obs",
                                "complete.obs", "everything", "na.or.complete"),
                        cor_method = c("pearson", "spearman", "kendall", "poly",
                                       "tetra", "fiml"),
                        show_progress = TRUE, seed = NULL,
                        P_type = lifecycle::deprecated()) {

  # Perform argument checks
  .assert_cor_input(x)

  # `P_type` was renamed to the snake_case `p_type`; the old spelling is still
  # accepted (silently) for backwards compatibility and mapped onto `p_type` here.
  if (lifecycle::is_present(P_type)) p_type <- P_type

  checkmate::assert_count(n_factors)
  # Matches efa_fit()'s contract: NA means "sample size unknown", 0 is an impossible one.
  # Checked here as well as there, because the grid runs each fit inside try(), which would
  # otherwise turn the argument error into an unexplained all-solutions-failed result.
  # Raised through .assert_args() so the rejection carries the same condition class here as
  # it does from efa_fit().
  .assert_args(checkmate::assert_count(N, positive = TRUE, na.ok = TRUE))
  # The vector-valued choice arguments are case-insensitive: map them onto the
  # canonical spellings first, so the subset checks and the type grid below only
  # ever see the canonical values.
  estimator <- .map_subset_ci(estimator, c("PAF", "ML", "ULS", "MINRES"))
  checkmate::assert_subset(estimator, c("PAF", "ML", "ULS", "MINRES"),
                           empty.ok = FALSE)
  # "MINRES" is a synonym for "ULS" (same estimator); resolve to the canonical
  # name so the type grid below is keyed once per requested estimator.
  estimator[estimator == "MINRES"] <- "ULS"
  estimator <- unique(estimator)
  rotation_choices <- c("none", "orthogonal", "oblique", "varimax",
                        "equamax", "quartimax", "geominT", "bentlerT",
                        "bifactorT", "promax", "oblimin", "quartimin",
                        "simplimax", "bentlerQ", "geominQ", "bifactorQ")
  rotation <- .map_subset_ci(rotation, rotation_choices)
  checkmate::assert_subset(rotation, rotation_choices, empty.ok = FALSE)
  # Drop duplicates the case-folding may have collapsed (e.g. "promax" and "Promax"), as for
  # the estimators above, so the averaging grid runs each requested rotation and type once.
  rotation <- unique(rotation)
  type <- .map_subset_ci(type, c("none", "EFAtools", "psych", "SPSS"))
  checkmate::assert_subset(type, c("none", "EFAtools", "psych", "SPSS"),
                           empty.ok = FALSE)
  type <- unique(type)
  averaging <- .match_arg_ci(averaging)
  checkmate::assert_number(trim, lower = 0, upper = 0.5)
  # `trim` only reaches mean(); the median has no trimming to do. Say so rather than
  # silently recording a trim in the settings that never affected a single value.
  if (averaging == "median" && trim > 0) {
    cli::cli_inform(
      c("i" = "{.arg trim} is only used when {.code averaging = \"mean\"}; it is ignored for the median."),
      class = "efa_avg_trim_ignored"
    )
  }
  checkmate::assert_number(salience_threshold, lower = 0, upper = 1)
  checkmate::assert_count(max_iter)
  init_comm <- .map_subset_ci(init_comm, c("smc", "mac", "unity"))
  checkmate::assert_subset(init_comm, c("smc", "mac", "unity"),
                           empty.ok = FALSE)
  init_comm <- unique(init_comm)
  checkmate::assert_vector(criterion, strict = TRUE, any.missing = FALSE,
                           min.len = 1)
  checkmate::assert_true(all(criterion > 0 & criterion < 1))
  criterion_type <- .map_subset_ci(criterion_type, c("max_individual", "sum"))
  checkmate::assert_subset(criterion_type, c("max_individual", "sum"),
                           empty.ok = FALSE)
  criterion_type <- unique(criterion_type)
  checkmate::assert_subset(abs_eigen, c(TRUE, FALSE),
                           empty.ok = FALSE)
  varimax_type <- .map_subset_ci(varimax_type, c("svd", "kaiser"))
  checkmate::assert_subset(varimax_type, c("svd", "kaiser"),
                           empty.ok = FALSE)
  varimax_type <- unique(varimax_type)
  checkmate::assert_subset(normalize, c(TRUE, FALSE),
                           empty.ok = FALSE)
  checkmate::assert_vector(k_promax, strict = TRUE, any.missing = FALSE,
                           min.len = 1)
  checkmate::assert_vector(k_simplimax, strict = TRUE, any.missing = FALSE,
                           min.len = 1)
  p_type <- .map_subset_ci(p_type, c("unnorm", "norm"))
  checkmate::assert_subset(p_type, c("unnorm", "norm"),
                           empty.ok = FALSE)
  p_type <- unique(p_type)
  checkmate::assert_vector(precision, strict = TRUE, any.missing = FALSE,
                           min.len = 1)
  checkmate::assert_true(all(precision > 0 & precision < 1))
  start_method <- .map_subset_ci(start_method, c("psych", "factanal"))
  checkmate::assert_subset(start_method, c("psych", "factanal"),
                           empty.ok = FALSE)
  start_method <- unique(start_method)
  use <- .match_arg_ci(use)
  cor_method <- .match_arg_ci(cor_method)
  checkmate::assert_flag(show_progress)



  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, N = N, use = use, cor_method = cor_method,
                             N_policy = "optional")
  R <- prep$R
  N <- prep$N

  # Check if model is identified

  # calculate degrees of freedom
  m <- ncol(R)
  df <- .efa_df(m, n_factors)

  if(df < 0){

    cli::cli_abort(
      c("The model is underidentified.",
        "i" = "Use fewer factors or more indicators, then try again."),
      class = "efa_underidentified"
    )

  } else if (df == 0){

    cli::cli_warn(
      c("The model is just identified ({.code df = 0}).",
        "i" = "Consider fewer factors or more indicators."),
      class = "efa_just_identified"
    )

  }

  if (n_factors == 1 && !all(rotation == "none")) {
    cli::cli_inform(
      c("i" = "{.arg n_factors} is 1 but {.arg rotation} is not {.val none}; setting {.arg rotation} to {.val none}, as single-factor solutions cannot be rotated."),
      class = "efa_avg_single_factor_rotation"
    )
    rotation <- "none"
  }

  ### create the grid with all combinations of the input arguments

  arg_grid <- .build_avg_grid(
    estimator = estimator, type = type, rotation = rotation, init_comm = init_comm,
    criterion = criterion, criterion_type = criterion_type, abs_eigen = abs_eigen,
    max_iter = max_iter, start_method = start_method, k_promax = k_promax,
    normalize = normalize, P_type = p_type, precision = precision,
    varimax_type = varimax_type, k_simplimax = k_simplimax)

  # The grid leaves precision NA for unrotated rows, where it is not applicable;
  # those rows still need a valid tolerance to pass to EFA(), so fall back to its
  # default. Kept as one constant so both the single- and multi-EFA paths agree.
  default_precision <- 1e-5

  # Iteration budget handed to the rotation engine, well above its own default so a
  # criterion-based rotation is not stopped short of the requested precision. Kept as one
  # constant for the same reason as `default_precision`: a rotation that converges inside
  # the grid must not fail to converge when the same settings collapse to a single row.
  rotation_maxit <- 5e4

  ### Run all efas

  # A supplied `seed` makes the averaging run reproducible and leaves the caller's RNG
  # stream untouched: it is saved and restored on exit (or, if none existed, the state
  # set.seed() creates is removed again). The grid's fits are stochastic whenever a
  # criterion-based rotation draws random starts, and future_lapply(future.seed = TRUE)
  # derives its per-fit L'Ecuyer streams from the state set here, which also makes the
  # grid independent of the number of parallel workers. Set before the single-combination
  # branch below so that shortcut is covered too. Mirrors the seed handling in efa_fit().
  .set_local_seed(seed)

  if (nrow(arg_grid) == 1) {

    cli::cli_warn("There was only one combination of arguments; returning a normal EFA output.",
                  class = "efa_avg_single_combination")

    single <- efa_fit(R, n_factors, N = N, estimator = arg_grid$estimator,
        rotation = arg_grid$rotation,
        estimate_control = estimate_control(type = "none",
            init_comm = arg_grid$init_comm, criterion = arg_grid$criterion,
            # The grid carries the row's iteration cap: the type's own on a named-type
            # PAF row (SPSS 25, psych 50, EFAtools 300), the max_iter argument on a PAF
            # row of type "none", and NA on an ML/ULS row that never uses it.
            criterion_type = arg_grid$criterion_type, max_iter = arg_grid$max_iter,
            abs_eigen = arg_grid$abs_eigen,
            # The grid leaves start_method NA for the non-ML rows that ignore it; the
            # estimate control requires a concrete value, so fall back to the default.
            start_method = if (is.na(arg_grid$start_method)) "psych" else arg_grid$start_method),
        rotate_control = rotate_control(type = "none",
            # normalize and precision are NA on the unrotated rows that ignore them; the
            # rotate control requires concrete values, so fall back to the defaults there.
            normalize = if (is.na(arg_grid$normalize)) TRUE else arg_grid$normalize,
            precision = if (is.na(arg_grid$precision)) default_precision else arg_grid$precision,
            # order_type is fixed to "eigen" (the ordering only affects the factor
            # order of this single returned fit, not any averaged quantity).
            order_type = "eigen", varimax_type = arg_grid$varimax_type,
            p_type = arg_grid$P_type,
            k = ifelse(arg_grid$rotation == "promax", arg_grid$k_promax, arg_grid$k_simplimax),
            # engine extra carried in the control, not efa_fit()'s dots: the
            # unrotated rows never consume it, and bare dots are rejected there
            maxit = rotation_maxit))
    # The fit drew from the stream .set_local_seed() just seeded, so this shortcut
    # result was produced under `seed` even though efa_fit() was not handed it; record
    # it, or the one documented return path that is not an efa_average object would be
    # the only one that cannot say how it was seeded.
    # (single-bracket assignment: `$<-` with a NULL value would DELETE the field
    # rather than set it, dropping it from the schema whenever no seed was supplied)
    single$settings["seed"] <- list(seed)
    return(single)

  }

  # Set the progress handler only for the duration of this call, restoring the
  # user's original setting on exit (capturing the raw option restores an unset
  # state to NULL rather than to an explicit empty handler list).
  old_handlers <- getOption("progressr.handlers")
  progressr::handlers(if (isTRUE(show_progress)) "cli" else "void")
  on.exit(options(progressr.handlers = old_handlers), add = TRUE)

  progressr::with_progress({
    n_efa <- nrow(arg_grid)
    if (n_efa <= 10) {
      stepsize <- 1
    } else {
      stepsize <- round(n_efa / 100 * 10)
    }

    # One declared step per update actually sent: the loop below reports on every
    # stepsize-th EFA and on the last one, which is `ceiling(n_efa / stepsize)` updates,
    # plus the closing message. Declaring fewer completes the bar while EFAs are still
    # being fitted, and an update sent to an already-completed progressor is dropped
    # (the closing message would never be shown) or makes with_progress() warn that it
    # is no longer listening to this progressor.
    efa_progress_bar <- progressr::progressor(steps = ceiling(n_efa / stepsize) + 1)
    efa_progress_bar("Running EFAs:", class = "sticky", amount = 0)
    efa_list <- future.apply::future_lapply(1:n_efa,
                                            function(i, estimators, rotations,
                                                     init_comms, criteria,
                                                     criterion_types, abs_eigens,
                                                     max_iters, varimax_types, k_ps,
                                                     k_ss, normalizes, P_types,
                                                     precisions, start_methods) {
      # Report every stepsize-th EFA, and the last one whatever its position, so the
      # bar tracks the grid to its end instead of completing on the last whole group.
      if (i %% stepsize == 0 || i == n_efa){
        efa_progress_bar(message = sprintf("Running EFA %g of %g", i, n_efa))
      }

      # precision only affects rotated solutions; the grid leaves it NA for the
      # unrotated rows, so fall back to the default tolerance there.
      precision_i <- if (is.na(precisions[i])) default_precision else precisions[i]

      # Per-EFA warnings (e.g. Heywood cases) are summarised once after the grid
      # is run, so they are suppressed here to avoid one warning per model.
      try(suppressWarnings(
          efa_fit(R, n_factors, N = N, estimator = estimators[i], rotation = rotations[i],
              estimate_control = estimate_control(type = "none", init_comm = init_comms[i],
                  criterion = criteria[i], criterion_type = criterion_types[i],
                  # Per-row iteration cap: the type's own on a named-type PAF row
                  # (SPSS 25, psych 50, EFAtools 300), the max_iter argument on a PAF
                  # row of type "none", and NA on the ML/ULS rows that never use it.
                  max_iter = max_iters[i], abs_eigen = abs_eigens[i],
                  # The grid leaves start_method NA for the non-ML rows that ignore it; the
                  # estimate control requires a concrete value, so fall back to the default.
                  start_method = if (is.na(start_methods[i])) "psych" else start_methods[i]),
              rotate_control = rotate_control(type = "none",
                  # normalize is NA on the unrotated rows that ignore it; the rotate control
                  # requires a concrete value, so fall back to the default there.
                  normalize = if (is.na(normalizes[i])) TRUE else normalizes[i],
                  # order_type is fixed to "eigen": the averaging re-aligns every solution
                  # to the first by congruence, so the per-fit factor order is irrelevant.
                  precision = precision_i, order_type = "eigen",
                  varimax_type = varimax_types[i], p_type = P_types[i],
                  k = ifelse(rotations[i] == "promax", k_ps[i], k_ss[i]),
                  # engine extra carried in the control, not efa_fit()'s dots: the
                  # unrotated rows never consume it, and bare dots are rejected there
                  maxit = rotation_maxit))),
          silent = TRUE)
    }, estimators = arg_grid$estimator, rotations = arg_grid$rotation,
    init_comms = arg_grid$init_comm, criteria = arg_grid$criterion,
    criterion_types = arg_grid$criterion_type, abs_eigens = arg_grid$abs_eigen,
    max_iters = arg_grid$max_iter, varimax_types = arg_grid$varimax_type,
    k_ps = arg_grid$k_promax, k_ss = arg_grid$k_simplimax,
    normalizes = arg_grid$normalize, P_types = arg_grid$P_type,
    precisions = arg_grid$precision, start_methods = arg_grid$start_method,
    future.seed = TRUE)

    # The last of the declared steps, so the bar completes on this message rather than
    # on the final EFA. It has to consume a step: progressr drops an update sent to a
    # progressor that has already reached its last one.
    efa_progress_bar(message = "Done Running EFAs")
  })

  names(efa_list) <- rownames(arg_grid)

  ### Extract relevant information from EFA outputs
  if (isTRUE(show_progress)) {
    # The post-grid stages run sequentially in this process, so they are reported with
    # cli's own progress steps (each call ticks off the previous one). The grid above
    # keeps progressr, which is what relays progress out of the future workers.
    cli::cli_progress_step("Extracting data")
  }
  ext_list <- .extract_data(efa_list, R, n_factors, n_efa, rotation, salience_threshold)

  bad <- ext_list$for_grid$converged != 0 | ext_list$for_grid$errors |
         ext_list$for_grid$heywood

  # Summarise the per-model problems once (the individual EFA warnings were
  # suppressed during the grid run to avoid one warning per model).
  n_excluded <- sum(bad %in% TRUE)
  if (n_excluded > 0L) {
    cli::cli_warn(
      c(paste("Excluded {n_excluded} of {length(bad)} EFA{?s} from averaging due to",
              "errors, non-convergence, or Heywood cases."),
        "i" = "Inspect {.field implementations_grid} in the returned object for the per-model details."),
      class = "efa_avg_excluded_solutions"
    )
  }

  if (all(bad %in% TRUE)) {

    av_list <- list(
      h2 = NA,
      loadings = NA,
      phi = NA,
      ind_fac_corres = NA,
      vars_accounted = NA,
      fit_indices = NA)

  } else {

    if (n_factors > 1) {
      if (isTRUE(show_progress)) {
        cli::cli_progress_step("Reordering factors")
      }

      re_list <- .array_reorder(ext_list$vars_accounted, ext_list$L, ext_list$L_corres,
                                ext_list$phi, ext_list$extract_phi, n_factors)
    } else {
      re_list <- ext_list
    }


    if (isTRUE(show_progress)) {
      cli::cli_progress_step("Averaging data")
    }
    av_list <- .average_values(re_list$vars_accounted, re_list$L, re_list$L_corres,
                               ext_list$h2, re_list$phi, ext_list$extract_phi,
                               averaging, trim,
                               ext_list$for_grid[, c("chisq", "p_chi", "caf", "cfi",
                                                     "rmsea", "aic", "bic", "srmr",
                                                     "tli", "ecvi", "rmsr")], df,
                               colnames(R))

  }


  arg_grid <- cbind(arg_grid, ext_list$for_grid)


  settings <- list(
    estimator = estimator,
    # back-compat alias, as for the frozen P_type key
    method = estimator,
    rotation = rotation,
    type = type,
    n_factors = n_factors,
    N = N,
    init_comm = init_comm,
    criterion = criterion,
    criterion_type = criterion_type,
    abs_eigen = abs_eigen,
    varimax_type = varimax_type,
    normalize = normalize,
    k_promax = k_promax,
    k_simplimax = k_simplimax,
    P_type = p_type,
    precision = precision,
    start_method = start_method,
    use = use,
    cor_method = cor_method,
    max_iter = max_iter,
    averaging = averaging,
    trim = trim,
    salience_threshold = salience_threshold,
    seed = seed
  )

  # Create output
  output <- list(
    orig_R = R,
    h2 = av_list$h2,
    loadings = av_list$loadings,
    Phi = av_list$phi,
    ind_fac_corres = av_list$ind_fac_corres,
    vars_accounted = av_list$vars_accounted,
    fit_indices = av_list$fit_indices,
    implementations_grid = arg_grid,
    efa_list = efa_list,
    settings = settings
  )

  class(output) <- c("efa_average", "EFA_AVERAGE")

  if (isTRUE(show_progress)) {
    # Ticks off the last step that was started; which one that is depends on how far
    # the branches above got (no solution to average leaves it at "Extracting data").
    cli::cli_progress_done()
  }

  return(output)

}

