#' Exploratory factor analysis (EFA)
#'
#' This function does an EFA with either `PAF`, `ML`,
#' or `ULS` with or without subsequent rotation.
#' All arguments with default value `NA` can be left to default if `type`
#' is set to one of "EFAtools", "SPSS", or "psych". The respective specifications are
#' then handled according to the specified type (see details). For all rotations
#' except varimax and promax, the `GPArotation` package is needed.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations. If raw data is entered, the correlation matrix is found from the
#' data.
#' @param n_factors numeric. Number of factors to extract.
#' @param N numeric. The number of observations. Needs only be specified if a
#' correlation matrix is used. If input is a correlation matrix and `N` = NA
#' (default), not all fit indices can be computed. When raw data with missing
#' values are entered and `use` is `"complete.obs"` or `"na.or.complete"`, rows
#' are deleted listwise, so `N` is taken as the number of complete cases.
#' @param method character. One of "PAF", "ML", or "ULS" to use principal axis
#' factoring, maximum likelihood, or unweighted least squares, respectively, to fit
#' the EFA. "MINRES" is accepted as a synonym for "ULS": minimum residual and
#' unweighted least squares are two names for the same estimator and return identical
#' results.
#' @param rotation character. Either perform no rotation ("none"; default),
#' an orthogonal rotation ("varimax", "equamax", "quartimax", "geominT",
#' "bentlerT", or "bifactorT"), or an oblique rotation ("promax", "oblimin",
#' "quartimin", "simplimax", "bentlerQ", "geominQ", or "bifactorQ").
#' @param se character. Whether and how standard errors should be computed.
#'  Currently, only "np-boot" for non-parametric bootstrap is available
#'  and needs raw data to work. Default is "none".
#' @param type character. If one of "EFAtools" (default), "psych", or "SPSS" is
#'  used, and the following arguments with default NA are left with
#'  NA, these implementations are executed according to the respective program
#'  ("psych" and "SPSS") or according to the best solution found in Grieder &
#'  Steiner (2020; "EFAtools"). Individual properties can be adapted using one of
#'  the three types and specifying some of the following arguments. If set to
#'  "none" additional arguments must be specified depending on the `method`
#'  and `rotation` used (see details).
#' @param max_iter numeric. The maximum number of iterations to perform after which
#' the iterative PAF procedure is halted with a warning. If `type` is one of
#' "EFAtools", "SPSS", or "psych", this is automatically specified if `max_iter` is
#' left to be `NA`, but can be overridden by entering a number. Default is
#' `NA`.
#' @param init_comm character. The method to estimate the initial communalities
#' in `PAF`. "smc" will use squared multiple correlations, "mac" will use
#' maximum absolute correlations, "unity" will use 1s (see details).
#' Default is `NA`.
#' @param criterion numeric. The convergence criterion used for PAF.
#' If the change in communalities from one iteration to the next is smaller than
#' this criterion the solution is accepted and the procedure ends.
#' Default is `NA`.
#' @param criterion_type character. Type of convergence criterion used for
#' PAF. "max_individual" selects the maximum change in any of the
#' communalities from one iteration to the next and tests it against the
#' specified criterion. This is also used by SPSS. "sum" takes the difference of
#' the sum of all communalities in one iteration and the sum of all communalities
#' in the next iteration and tests this against the criterion. This procedure is
#' used by the [psych::fa()] function. Default is `NA`.
#' @param abs_eigen logical. Which algorithm to use in the PAF
#' iterations. If FALSE, the loadings are computed from the eigenvalues. This is
#' also used by the [psych::fa()] function. If TRUE the
#' loadings are computed with the absolute eigenvalues as done by SPSS.
#' Default is `NA`.
#' @param use character. Passed to [stats::cor()] if raw data
#' is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. Passed to [stats::cor()].
#' Default is "pearson".
#' @param k numeric. Either the power used for computing the target matrix P in
#' the promax rotation or the number of 'close to zero loadings' for the simplimax
#' rotation (see \code{\link[GPArotation:GPA]{GPArotation::GPFoblq}}). If left to
#' `NA` (default), the value for promax depends on the specified type.
#' For simplimax, `nrow(L)`, where L is the matrix of unrotated loadings,
#' is used by default.
#' @param normalize logical. If `TRUE`, a kaiser normalization is
#' performed before the specified rotation. Default is `TRUE`.
#' @param P_type character. This specifies how the target
#' matrix P is computed in promax rotation. If "unnorm" it will use the
#' unnormalized target matrix as originally done in Hendrickson and White (1964).
#' This is also used in the psych and stats packages. If "norm" it will use the
#' normalized target matrix as used in SPSS. Default is `NA`.
#' @param precision numeric. The tolerance for stopping in the rotation
#' procedure. Default is 10^-5 for all rotation methods.
#' @param varimax_type character. The type of the varimax rotation performed.
#' If "svd", singular value decomposition is used, as [stats::varimax()] does. If "kaiser", the varimax procedure performed in SPSS is used.
#' This is the original procedure from Kaiser (1958), but with slight alterations
#' in the varimax criterion (see details, and Grieder & Steiner, 2020). Default is `NA`.
#' @param order_type character. How to order the factors. "eigen" will reorder
#' the factors according to the largest to lowest eigenvalues of the matrix of
#' rotated loadings. "ss_factors" will reorder the factors according to descending
#' sum of squared factor loadings per factor. Default is `NA`.
#' @param start_method character. How to specify the starting values for the
#' optimization procedure for ML. Default is "psych" which takes the
#' starting values specified in [psych::fa()]. "factanal" takes the
#' starting values specified in the [stats::factanal()] function.
#' Solutions are very similar.
#' @param b_boot numeric. The number of bootstrap samples to draw. Default is 1000.
#' @param ci numeric. The confidence interval to create from the bootstrap samples.
#'  Must be between 0 and 1. Default ist .95 for 95% CIs.
#' @param randomStarts numeric. The number of random starts to use in rotations
#'  that use the `GPArotation` package. Some rotations are prone to produce
#'  local minima and sometimes many random starts are needed (see the GPArotation
#'  package documentation for details). Default is 10.
#' @param seed numeric. An optional seed for the random-number generator used by the
#'  non-parametric bootstrap (`se = "np-boot"`), i.e. for the case resampling, the
#'  rotation random starts, and the Procrustes random starts. Setting it makes the
#'  bootstrap reproducible and, importantly, independent of the number of parallel
#'  workers (see Details); the caller's random-number stream is restored afterwards,
#'  so supplying a seed leaves no lasting effect on it. Default is `NULL`, which uses
#'  (and advances) the current state of the generator.
#' @param ... Additional arguments passed to rotation functions from the `GPArotation` package (e.g., `maxit` for maximum number of iterations).
#'
#' @details There are two main ways to use this function. The easiest way is to
#' use it with a specified `type` (see above), which sets most of the other
#' arguments accordingly. Another way is to use it more flexibly by explicitly
#' specifying all arguments used and set `type` to "none" (see examples).
#' A mix of the two can also be done by specifying a `type` as well as
#' additional arguments. However, this will throw warnings to avoid unintentional
#' deviations from the implementations according to the specified `type`.
#'
#' The `type` argument is evaluated for PAF and for all rotations (mainly
#' important for the varimax and promax rotations). The type-specific settings
#' for these functions are detailed below.
#'
#' For PAF, the values of `init_comm`, `criterion`, `criterion_type`,
#' and `abs_eigen` depend on the `type` argument.
#'
#' `type = "EFAtools"` will use the following argument specification:
#' `init_comm = "smc", criterion = .001, criterion_type = "sum",
#' abs_eigen = TRUE`.
#'
#' `type = "psych"` will use the following argument specification:
#' `init_comm = "smc", criterion = .001, criterion_type = "sum",
#' abs_eigen = FALSE`.
#'
#' `type = "SPSS"` will use the following argument specification:
#' `init_comm = "smc", criterion = .001, criterion_type = "max_individual",
#' abs_eigen = TRUE`.
#'
#' If SMCs fail, SPSS takes "mac". However, as SPSS takes absolute eigenvalues,
#' this is hardly ever the case. Psych, on the other hand, takes "unity" if SMCs
#' fail, but uses the Moore-Penrose Psudo Inverse of a matrix, thus, taking "unity"
#' is only necessary if negative eigenvalues occur afterwards in the iterative
#' PAF procedure. The EFAtools type setting combination was the best in terms of accuracy
#' and number of Heywood cases compared to all the
#' other setting combinations tested in simulation studies in Grieder & Steiner
#' (2020), which is why this type is used as a default here.
#'
#' For varimax, the values of `varimax_type` and `order_type` depend on
#' the `type` argument.
#'
#' `type = "EFAtools"` will use the following argument specification:
#' `varimax_type = "kaiser", order_type = "eigen"`.
#'
#' `type = "psych"` will use the following argument specification:
#' `varimax_type = "svd", order_type = "eigen"`.
#'
#' `type = "SPSS"` will use the following argument specification:
#' `varimax_type = "kaiser", order_type = "ss_factors"`.
#'
#' For promax, the values of `P_type`,
#' `order_type`, and `k` depend on the `type` argument.
#'
#' `type = "EFAtools"` will use the following argument specification:
#' `P_type = "norm", order_type = "eigen", k = 4`.
#'
#' `type = "psych"` will use the following argument specification:
#' `P_type = "unnorm", order_type = "eigen", k = 4`.
#'
#' `type = "SPSS"` will use the following argument specification:
#' `P_type = "norm", order_type = "ss_factors", k = 4`.
#'
#' The `P_type` argument can take two values, "unnorm" and "norm". It controls
#' which formula is used to compute the target matrix P in the promax rotation.
#' "unnorm" uses the formula from Hendrickson and White (1964), specifically:
#' `P = abs(A^(k + 1)) / A`,
#' where A is the unnormalized matrix containing varimax rotated loadings.
#' "norm" uses the normalized varimax rotated loadings. Specifically it used the
#' following formula, which can be found in the SPSS 23 and SPSS 27 Algorithms manuals:
#' `P = abs(A / sqrt(rowSums(A^2))) ^(k + 1) * (sqrt(rowSums(A^2)) / A)`.
#' As for PAF, the EFAtools type setting combination for promax was the best
#' compared to the other setting combinations tested in simulation studies in
#' Grieder & Steiner (2020).
#'
#' The `varimax_type` argument can take two values, "svd", and "kaiser". "svd" uses
#' singular value decomposition, by calling [stats::varimax()]. "kaiser"
#' performs the varimax procedure as described in the SPSS 23 Algorithms manual and as described
#' by Kaiser (1958). However, there is a slight alteration in computing the varimax criterion, which
#' we found to better align with the results obtain from SPSS. Specifically, the original varimax
#' criterion as described in the SPSS 23 Algorithms manual is
#' `sum(n*colSums(lambda ^ 4) - colSums(lambda ^ 2) ^ 2) / n ^ 2`, where n is the
#' number of indicators, and lambda is the rotated loadings matrix. However, we found the following
#' to produce results more similar to those of SPSS:
#' `sum(n*colSums(abs(lambda)) - colSums(lambda ^ 4) ^ 2) / n^2`.
#'
#' For all other rotations except varimax and promax, the `type` argument
#' only controls the `order_type` argument with the same values as stated
#' above for the varimax and promax rotations. For these other rotations, the
#' `GPArotation` package is needed. Additional arguments can also be
#' specified and will be passed to the respective `GPArotation` function
#' (e.g., maxit to change the maximum number of iterations for the rotation procedure).
#'
#' The `type` argument has no effect on ULS and ML. For ULS, no additional
#' arguments are needed. For ML, an additional argument
#' `start_method` is needed to determine the starting values for the
#' optimization procedure. Default for this argument is "psych" which takes
#' the starting values specified in [psych::fa()].
#'
#' When `se = "np-boot"`, the bootstrap replicate fits are run in parallel across
#' replicates with the `future` framework. By default they run sequentially; to run
#' them in parallel, register a plan with [future::plan()] (e.g.
#' `future::plan(future::multisession, workers = 2)`; see examples). With a fixed
#' `seed` the bootstrap is reproducible and yields the same result regardless of the
#' number of workers. Each worker runs its own (Armadillo) linear algebra, so if your
#' `BLAS` is multi-threaded, limit the number of workers (or the BLAS threads) to
#' avoid over-subscribing the available cores.
#'
#'
#' @return A list of class EFA containing (a subset of) the following:
#'
#' \item{orig_R}{Original correlation matrix.}
#' \item{h2_init}{Initial communality estimates from PAF.}
#' \item{h2}{Final communality estimates from the unrotated solution.}
#' \item{orig_eigen}{Eigen values of the original correlation matrix.}
#' \item{init_eigen}{Initial eigenvalues, obtained from the correlation matrix
#'  with the initial communality estimates as diagonal in PAF.}
#' \item{final_eigen}{Eigenvalues obtained from the correlation matrix
#'  with the final communality estimates as diagonal.}
#' \item{iter}{The number of iterations needed for convergence.}
#' \item{convergence}{Integer convergence code (0 = converged). For ML and ULS
#'  this is the code returned by [`stats::optim()`][stats::optim]; for PAF it is
#'  1 if the maximum number of iterations was reached without meeting the
#'  convergence criterion and 0 otherwise.}
#' \item{heywood}{A named integer vector indicating which variables have a
#'  Heywood (improper) case in the unrotated solution; empty if there are none.}
#' \item{unrot_loadings}{Loading matrix containing the final unrotated loadings.}
#' \item{vars_accounted}{Matrix of explained variances and sums of squared loadings. Based on the unrotated loadings.}
#' \item{fit_indices}{For ML and ULS: Fit indices derived from the unrotated
#' factor loadings: Chi Square, including significance level, degrees of freedom
#' (df), Comparative Fit Index (CFI), Tucker-Lewis Index (TLI, also called the
#' non-normed fit index; Tucker & Lewis, 1973), Root Mean Square Error of
#' Approximation (RMSEA), including its 90% confidence interval, Akaike
#' Information Criterion (AIC), Bayesian Information Criterion (BIC), Expected
#' Cross-Validation Index (ECVI; Browne & Cudeck, 1989), Root Mean Squared
#' Residual (RMSR), Standardized Root Mean Squared Residual (SRMR; Bentler,
#' 1995), and the common part accounted for (CAF) index as proposed by
#' Lorenzo-Seva, Timmerman, & Kiers (2011). The Chi Square, CFI, TLI, RMSEA,
#' AIC, BIC, and ECVI are based on the Bartlett-corrected Chi Square (matching
#' [stats::factanal()] for ML). For PAF, only the CAF, RMSR, SRMR, and dfs are
#' returned. Note that while in Lorenzo-Seva, Timmerman, & Kiers (2011) the CAF is introduced as ranging between 0 and 1, with values close to 1 indicating close fit, this does not match the formula they introduce for calculating CAF: `1 - KMO(residuals)`, which only works if the diagonal of the residual matrix is set to 1s and will then approximate 0.5 with close fit.}
#' \item{model_implied_R}{The model implied correlation
#' matrix.}
#' \item{residuals}{Residual correlations, i.e., orig_R - model_implied_R}
#' \item{standardized_residuals}{Residual correlations standardized by their
#'  bootstrap standard errors. Only returned, if `se = "np-boot"`.}
#' \item{rot_loadings}{Loading matrix containing the final rotated loadings
#' (pattern matrix).}
#' \item{Phi}{The factor intercorrelations (only for oblique rotations).}
#' \item{Structure}{The structure matrix (only for oblique rotations).}
#' \item{rotmat}{The rotation matrix.}
#' \item{vars_accounted_rot}{Matrix of explained variances and sums of squared
#' loadings. Based on rotated loadings and, for oblique rotations, the factor
#' intercorrelations.}
#' \item{settings}{A list of the settings used.}
#' \item{boot.SE}{A list bootstrap standard errors for loadings (rotated and unrotated), structure coefficients (if rotated obliquely), factor correlations (Phi, only if rotated), and fit indices. Only returned, if `se = "np-boot"`.}
#' \item{boot.CI}{A list bootstrap confidence intervals of width `ci` for loadings (rotated and unrotated), structure coefficients (if rotated obliquely), factor correlations (Phi, if obliquely rotated), and fit indices. Only returned, if `se = "np-boot"`.}
#' \item{boot.arrays}{A list of arrays with the bootstrapped loadings (aligned rotated and unrotated), aligned structure coefficients (if rotated obliquely), aligned factor correlations (Phi, if obliquely rotated), and fit indices. Only returned, if `se = "np-boot"`.}
#'
#' @source Grieder, S., & Steiner, M.D. (2020). Algorithmic Jingle Jungle:
#' A Comparison of Implementations of Principal Axis Factoring and Promax Rotation
#'  in R and SPSS. Manuscript in Preparation.
#' @source Hendrickson, A. E., & White, P. O. (1964). Promax: A quick method for
#' rotation to oblique simple structure. British Journal of Statistical Psychology,
#' 17 , 65–70. doi: 10.1111/j.2044-8317.1964.tb00244.x
#' @source Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. L. (2011). The
#' Hull Method for Selecting the Number of Common Factors, Multivariate Behavioral
#' Research, 46, 340-364, doi: 10.1080/00273171.2011.564527
#' @source Kaiser, H. F. (1958). The varimax criterion for analytic rotation in
#' factor analysis. Psychometrika, 23, 187–200. doi: 10.1007/BF02289233
#'
#' @export
#'
#' @examples
#' # A type EFAtools (as presented in Steiner and Grieder, 2020) EFA
#' EFAtools_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                     type = "EFAtools", method = "PAF", rotation = "none")
#'
#' # A type SPSS EFA to mimick the SPSS implementation (this will throw a warning,
#' # see below)
#' SPSS_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                 type = "SPSS", method = "PAF", rotation = "none")
#'
#' # A type psych EFA to mimick the psych::fa() implementation
#' psych_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                  type = "psych", method = "PAF", rotation = "none")
#'
#' # Use ML instead of PAF with type EFAtools
#' EFAtools_ML <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                    type = "EFAtools", method = "ML", rotation = "none")
#'
#' # Use oblimin rotation instead of no rotation with type EFAtools
#' EFAtools_oblim <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                       type = "EFAtools", method = "PAF", rotation = "oblimin")
#'
#' # Do a PAF without rotation without specifying a type, so the arguments
#' # can be flexibly specified (this is only recommended if you know what you're
#' # doing)
#' PAF_none <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                 type = "none", method = "PAF", rotation = "none",
#'                 max_iter = 500, init_comm = "mac", criterion = 1e-4,
#'                 criterion_type = "sum", abs_eigen = FALSE)
#'
#' # Add a promax rotation
#' PAF_pro <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                type = "none", method = "PAF", rotation = "promax",
#'                max_iter = 500, init_comm = "mac", criterion = 1e-4,
#'                criterion_type = "sum", abs_eigen = FALSE, k = 3,
#'                P_type = "unnorm", precision= 1e-5, order_type = "eigen",
#'                varimax_type = "svd")
#'
#' \dontrun{
#' # Bootstrap standard errors from raw data, reproducible via a fixed seed and run
#' # in parallel across replicates. Keep the number of workers small (here 2) so the
#' # workers do not over-subscribe the cores against a multi-threaded BLAS.
#' future::plan(future::multisession, workers = 2)
#' EFA_boot <- EFA(GRiPS_raw, n_factors = 1, method = "PAF", rotation = "none",
#'                 se = "np-boot", b_boot = 1000, seed = 42)
#' future::plan(future::sequential)
#' }
#'
EFA <- function(x, n_factors, N = NA, method = c("PAF", "ML", "ULS", "MINRES"),
                rotation = c("none", "varimax", "equamax", "quartimax", "geominT",
                             "bentlerT", "bifactorT", "promax", "oblimin",
                             "quartimin", "simplimax", "bentlerQ", "geominQ",
                             "bifactorQ"),
                se = c("none", "np-boot"),
                type = c("EFAtools", "psych", "SPSS", "none"), max_iter = NA,
                init_comm = NA, criterion = NA, criterion_type = NA,
                abs_eigen = NA, use = c("pairwise.complete.obs", "all.obs",
                                          "complete.obs", "everything",
                                          "na.or.complete"),
                varimax_type = NA,
                k = NA, normalize = TRUE, P_type = NA, precision = 1e-5,
                order_type = NA, start_method = "psych",
                cor_method = c("pearson", "spearman", "kendall"),
                b_boot = 1000, ci = .95,
                randomStarts = 10, seed = NULL,
                ...) {

  # Perform argument checks
  .assert_cor_input(x)

  method <- match.arg(method)
  # "MINRES" is a synonym for "ULS" (same estimator); resolve it once here so the
  # rest of EFA() and the reported settings use the single canonical name.
  if (method == "MINRES") method <- "ULS"
  rotation <- match.arg(rotation)
  se <- match.arg(se)
  np_boot <- se == "np-boot"
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  type <- match.arg(type)
  start_method <- checkmate::matchArg(start_method, c("psych", "factanal", NA))

  if (is.na(start_method) && method == "ML") {
    cli::cli_abort(
      c("{.arg start_method} must be set when {.code method = \"ML\"}.",
        "i" = "Set {.arg start_method} to {.val psych} or {.val factanal}."),
      class = "efa_ml_start_missing"
    )
  }

  checkmate::assert_count(n_factors)
  checkmate::assert_count(N, na.ok = TRUE)
  checkmate::assert_count(max_iter, na.ok = TRUE)
  checkmate::assert_choice(init_comm, c("smc", "mac", "unity", NA))
  checkmate::assert_number(criterion, lower = 0, upper = 1, na.ok = TRUE)
  if (!is.na(criterion) && criterion >= 1) {
    cli::cli_abort(
      c("{.arg criterion} must be smaller than 1.",
        "x" = "You supplied {.arg criterion} = {criterion}.",
        "i" = "Use a small positive convergence tolerance such as {.val 0.001}."),
      class = "efa_criterion_too_large"
    )
  }
  checkmate::assert_choice(criterion_type, c("max_individual", "sum", NA))
  checkmate::assert_flag(abs_eigen, na.ok = TRUE)
  checkmate::assert_number(k, na.ok = TRUE)
  checkmate::assert_choice(varimax_type, c("svd", "kaiser", NA))
  checkmate::assert_flag(normalize, na.ok = TRUE)
  checkmate::assert_choice(P_type, c("unnorm", "norm", NA))
  checkmate::assert_number(precision, lower = 0, upper = 1)
  checkmate::assert_choice(order_type, c("eigen", "ss_factors", NA))
  checkmate::assert_integerish(b_boot, len = 1, any.missing = FALSE)
  checkmate::assert_number(ci, lower = 0, upper = 1)
  checkmate::assert_int(seed, null.ok = TRUE)

  # The common-factor model requires fewer factors than variables; with
  # n_factors >= n_variables it is not identified and the eigenvalue-based
  # extraction in the ML, ULS, and PAF fitters reads past the available
  # eigenvalues (undefined behaviour in an unchecked build).
  n_vars <- ncol(x)
  if (n_factors >= n_vars) {
    cli::cli_abort(
      c("{.arg n_factors} must be smaller than the number of variables.",
        "x" = "You requested {n_factors} factor{?s} for {n_vars} variable{?s}.",
        "i" = "Extract fewer factors."),
      class = "efa_too_many_factors"
    )
  }

  # A correlation matrix cannot yield bootstrap SEs; resolve this before
  # validating/smoothing R so the warning is not pre-empted by a singular or
  # non-positive-definite matrix.
  is_cormat <- .is_cormat(x)

  if (is_cormat) {

    if (isTRUE(np_boot)) {
      cli::cli_warn(
        c("Cannot compute bootstrap standard errors from correlation matrix.",
        "x" = "You've supplied {.var se} = {.val {se}}, but {.var x} is a correlation matrix.",
        "i" = "Setting {.var se} to {.val none}. Rerun with raw data to calculate bootstrap SEs.")
      )
    }

    np_boot <- FALSE
    se <- "none"

  }

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, N = N, use = use, cor_method = cor_method,
                             N_policy = "optional",
                             check_singular = type != "psych",
                             posdef_abort = type == "SPSS")
  R <- prep$R
  N <- prep$N

  if (!is_cormat && isTRUE(np_boot)) {

    # A fixed seed makes the whole bootstrap reproducible and independent of the
    # number of parallel workers: the case resampling, the rotation random starts,
    # and the Procrustes random starts all run from this state. The case resampling
    # advances the global RNG by a b_boot-dependent amount; the parallel replicate
    # fit then adds a fixed, worker-count-independent step when future.seed = TRUE
    # derives a per-replicate L'Ecuyer stream. Both advances are deterministic given
    # the seed, so the downstream draws -- and the result -- are identical at any
    # number of workers. The caller's RNG stream is saved and restored afterwards --
    # or, if none existed yet, the seed set.seed() creates is removed again -- so
    # EFA() leaves no side effect on it.
    if (!is.null(seed)) {
      seed_existed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
      saved_seed <- if (seed_existed) {
        get(".Random.seed", envir = globalenv(), inherits = FALSE)
      }
      on.exit({
        if (seed_existed) {
          assign(".Random.seed", saved_seed, envir = globalenv())
        } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
          rm(".Random.seed", envir = globalenv())
        }
      }, add = TRUE)
      set.seed(seed)
    }

    m <- ncol(R)
    # Resample the cases the correlation matrix was actually built from. Under
    # listwise deletion that is the complete cases (N of them), not the first N
    # row positions, so the case bootstrap stays a faithful resample of the
    # estimator that produced R (Efron & Tibshirani, 1993).
    rows <- if (use %in% c("complete.obs", "na.or.complete")) {
      which(stats::complete.cases(x))
    } else {
      seq_len(nrow(x))
    }

    # create bootstrap samples and from these, correlation matrices
    R_boot_array <- array(NA_real_, c(m, m, b_boot), dimnames = list(colnames(x),
                                                               colnames(x),
                                                               NULL))

    for (boot_i in seq_len(b_boot)) {
      ind <- sample(rows, size = N, replace = TRUE)
      R_boot_array[,, boot_i] <- stats::cor(x[ind,], use = use, method = cor_method)
    }

  }

  # Check if model is identified

  # calculate degrees of freedom
  m <- ncol(R)
  df <- ((m - n_factors)**2 - (m + n_factors)) / 2

  if(df < 0){

    cli::cli_warn(
      c("The model is underidentified; no fit indices were computed.",
        "i" = "Use fewer factors or more indicators."),
      class = "efa_underidentified"
    )

  } else if (df == 0){

    cli::cli_warn(
      c("The model is just identified ({.code df = 0}).",
        "i" = "Consider fewer factors or more indicators."),
      class = "efa_just_identified"
    )

  }

  # run factor analysis with respective fit method

  if (method %in% c("ML", "ULS")) {

    if (type == "SPSS") {

      cli::cli_warn(
        c("Only {.val PAF} is validated against the SPSS implementation.",
          "i" = "{.val {method}} results may differ from those returned by SPSS."),
        class = "efa_spss_method_untested"
      )

    }

    if (is.na(N)) {

      cli::cli_warn(
        c("{.arg N} is {.val NA}; not all fit indices could be computed.",
          "i" = "Provide {.arg N} or raw data to compute all fit indices."),
        class = "efa_fit_na_n"
      )

    }

  }

  fit_out <- .estimate_model(R, method = method, n_factors = n_factors, N = N,
                             type = type, max_iter = max_iter,
                             init_comm = init_comm, criterion = criterion,
                             criterion_type = criterion_type,
                             abs_eigen = abs_eigen, start_method = start_method)

  # Surface Heywood cases from the point-estimate solution (the detector runs in
  # .finalize_fit for every fit). This fires once per EFA() call; EFA_AVERAGE,
  # which fits one EFA per grid cell, suppresses these per-model warnings and
  # reports a single summary instead.
  if (length(fit_out$heywood) > 0) {
    heywood_vars <- names(fit_out$heywood)
    cli::cli_warn(
      c(paste("{cli::qty(heywood_vars)}Heywood case{?s} detected for {.val {heywood_vars}}:",
              "the solution is improper (a communality at or above 1, or a uniqueness",
              "fixed at the estimation boundary)."),
        "i" = "Interpret the affected loadings and uniquenesses with caution."),
      class = "efa_heywood"
    )
  }

  if (isTRUE(np_boot)) {

    boot_fits <- .boot_fun(R_boot_array, b_boot, .estimate_model,
                           # .estimate_model arguments:
                           method = method, n_factors = n_factors, N = N,
                           type = type, max_iter = max_iter,
                           init_comm = init_comm, criterion = criterion,
                           criterion_type = criterion_type,
                           abs_eigen = abs_eigen, start_method = start_method,
                           # Each replicate fits only the quantities the bootstrap
                           # aggregation consumes (see .finalize_fit()).
                           lean = TRUE)

  }

  # rotate factor analysis results
  if (rotation == "none") {

    output <- fit_out
    boot_rot <- "none"

  } else {

    # Only promax and varimax are validated against the SPSS implementation.
    if (type == "SPSS" && !rotation %in% c("promax", "varimax")) {

      cli::cli_warn(
        c("Only the {.val promax} and {.val varimax} rotations are validated against the SPSS implementation.",
          "i" = "{.val {rotation}} results may differ from those returned by SPSS."),
        class = "efa_spss_rotation_untested"
      )

    }

    rot_out <- .rotate_model(fit_out, rotation = rotation, type = type,
                             normalize = normalize, precision = precision,
                             order_type = order_type, varimax_type = varimax_type,
                             P_type = P_type, k = k, randomStarts = randomStarts,
                             ...)

    boot_rot <- .rotation_family(rotation)

  }

  if (rotation != "none"){

    if(method == "ULS"){

      settings <- rot_out$settings
      output <- c(fit_out, within(rot_out, rm(settings)),
                  settings = list(settings))

    } else {

      settings <- c(fit_out$settings, rot_out$settings)
      output <- c(within(fit_out, rm(settings)), within(rot_out, rm(settings)),
                  settings = list(settings))

    }

  }

  # Add settings used to output
  settings_EFA <- list(
    method = method,
    rotation = rotation,
    type = type,
    n_factors = n_factors,
    N = N,
    use = use,
    cor_method = cor_method,
    se = se,
    b_boot = b_boot,
    ci = ci
  )

  if(method == "ULS" & rotation == "none"){

    output <- c(output, settings = list(settings_EFA))

  } else {

    settings <- c(settings_EFA, output$settings)

    output <- c(within(output, rm(settings)),
                settings = list(settings))

  }

  if (isTRUE(np_boot)) {
    if (rotation == "none") {
      L_rot <- NULL
    } else {
      L_rot <- rot_out$rot_loadings
    }
    boot_out <- .boot_se_ci(fit_out, L_rot,
                            boot_fits, boot_rot, ci, b_boot)
    output <- c(output, boot = boot_out)
    output$standardized_residuals <- output$residuals / boot_out$SE$residuals
  }

  class(output) <- "EFA"

  return(output)

}



.boot_fun <- function(x, b, call_fun, ...) {

  # The per-replicate fits are independent and estimation is RNG-free, so they are
  # run in parallel across replicates at the R/process level with future.apply. A
  # parallel processing plan can be selected with future::plan(); the default plan
  # runs sequentially. future.seed = TRUE assigns each replicate its own reproducible
  # L'Ecuyer-CMRG stream, so the bootstrap is reproducible and independent of the
  # number of workers (any RNG a fitter might draw is bound to the replicate index,
  # never to the worker).
  #
  # A replicate whose (possibly degenerate) resampled correlation matrix cannot be
  # fit returns NULL and is skipped later, rather than aborting the whole call.
  #
  # Per-replicate warnings are suppressed here: they repeat identical information b
  # times. The type/preset-override notice depends only on the (type, pinned
  # arguments) combination, which is the same for every replicate and already
  # surfaced once by the point-estimate fit; the iterative fitter's max-iteration
  # warning would otherwise fire once per non-converged replicate. Non-convergence
  # is instead tallied and reported once, after all replicates have been fitted.
  boot_list <- future.apply::future_lapply(seq_len(b), function(boot_i) {
    tryCatch(suppressWarnings(call_fun(x[,, boot_i], ...)),
             error = function(e) NULL)
  }, future.seed = TRUE)

  n_nonconverged <- sum(vapply(boot_list,
                               function(fit_i) isTRUE(fit_i$convergence != 0),
                               logical(1)))

  if (n_nonconverged > 0L) {
    cli::cli_warn(
      c("{n_nonconverged} of {b} bootstrap replicate{?s} did not converge.",
        "i" = "Their bootstrap standard errors and confidence intervals may be unreliable."),
      class = "efa_boot_nonconvergence"
    )
  }

  boot_list

}

.boot_se_ci <- function(fit_target, L_rot, boot_fit, boot_rot, ci, b) {

  l_ci <- (1 - ci) / 2
  ps <- c(l_ci, ci + l_ci)

  ### calculate stats for unrot loadings and gof measures
  L_unrot <- fit_target$unrot_loadings

  ncol_L <- ncol(L_unrot)
  nrow_L <- nrow(L_unrot)
  colnam_L <- colnames(L_unrot)
  rownam_L <- rownames(L_unrot)

  L_unrot_boot <- array(NA_real_, c(nrow_L, ncol_L, b),
                        dimnames = list(rownam_L, colnam_L,
                                        NULL))
  gof_boot <- matrix(NA_real_, ncol = length(fit_target$fit_indices), nrow = b)
  residuals_boot <- array(NA_real_, c(nrow_L, nrow_L, b),
                          dimnames = list(rownam_L, rownam_L,
                                          NULL))

  # Track replicates that could not be fit (NULL) or aligned, so they are
  # excluded from the bootstrap statistics rather than aborting the whole call.
  failed <- vapply(boot_fit, is.null, logical(1))

  for (boot_i in seq_len(b)) {

    if (failed[boot_i]) next

    # save aligned loading matrix
    aligned <- tryCatch(
      .align_solution(L_unrot, boot_fit[[boot_i]]$unrot_loadings),
      error = function(e) NULL
    )
    if (is.null(aligned)) {
      failed[boot_i] <- TRUE
      next
    }

    L_unrot_boot[,, boot_i] <- aligned$loadings
    gof_boot[boot_i, ] <- unlist(boot_fit[[boot_i]]$fit_indices)
    residuals_boot[,, boot_i] <- boot_fit[[boot_i]]$residuals

  }

  n_failed <- sum(failed)
  if (n_failed == b) {
    cli::cli_abort(
      c("All {b} bootstrap replicates failed; no bootstrap standard errors could be computed.",
        "i" = "The resampled correlation matrices may be degenerate; try more observations or fewer factors."),
      class = "efa_boot_all_failed"
    )
  }
  if (n_failed > 0) {
    cli::cli_warn(
      c("{n_failed} bootstrap replicate{?s} failed and {?was/were} excluded.",
        "i" = "Bootstrap standard errors and confidence intervals are based on {b - n_failed} replicate{?s}."),
      class = "efa_boot_replicate_failed"
    )
  }

  # se = sd of bootstrap replications (Zhang, 2014, Estimating Standard Errors
  # in Exploratory Factor Analysis)
  L_unrot_se_ci <- .array_se_ci(L_unrot_boot, ps)
  gof_se_ci <- .array_se_ci(gof_boot, ps, M = 2)
  residuals_se_ci <- .array_se_ci(residuals_boot, ps)
  names(gof_se_ci$se) <- names(fit_target$fit_indices)
  names(gof_se_ci$ci$lower) <- names(fit_target$fit_indices)
  names(gof_se_ci$ci$upper) <- names(fit_target$fit_indices)

  if(boot_rot == "oblique") {

    colnam_L <- colnames(L_rot)

    L_rot_boot <- array(NA_real_, c(nrow_L, ncol_L, b),
                        dimnames = list(rownam_L, colnam_L,
                                        NULL))
    Phi_rot_boot <- array(NA_real_, c(ncol_L, ncol_L, b),
                          dimnames = list(colnam_L, colnam_L,
                                          NULL))
    Structure_boot <- array(NA_real_, c(nrow_L, ncol_L, b),
                        dimnames = list(rownam_L, colnam_L,
                                        NULL))

    failed_rot <- 0

    # Align every successfully fit replicate to the target in a single compiled
    # call over the loading cube, rather than one PROCRUSTES() round trip per
    # replicate. Replicates that failed to fit or could not be aligned to the
    # point estimate are excluded from the cube and stay NA in the output arrays.
    keep <- which(!failed)
    if (length(keep) > 0) {
      A_cube <- array(NA_real_, c(nrow_L, ncol_L, length(keep)))
      for (j in seq_along(keep)) {
        A_cube[, , j] <- boot_fit[[keep[j]]]$unrot_loadings
      }

      aligned <- .oblique_procrustes_batch(A_cube, L_rot, random_starts = 5)

      for (j in seq_along(keep)) {
        boot_i <- keep[j]
        # Exclude a replicate only when no valid alignment could be produced. The
        # best multi-start fit is kept even if it did not formally converge: it is
        # the lowest-objective alignment available and its loadings are well-defined.
        if (!isTRUE(aligned$valid[j])) {
          failed_rot <- failed_rot + 1
          next
        }
        L_j <- matrix(aligned$loadings[, , j], nrow_L, ncol_L)
        Phi_j <- matrix(aligned$Phi[, , j], ncol_L, ncol_L)
        L_rot_boot[, , boot_i] <- L_j
        Phi_rot_boot[, , boot_i] <- Phi_j
        Structure_boot[, , boot_i] <- L_j %*% Phi_j
      }
    }

    valid_rot <- b - n_failed - failed_rot
    if (failed_rot > 0) {
      cli::cli_warn(c("{failed_rot} target rotation{?s} in the bootstrap procedure could not be aligned.",
                    "i" = "Bootstrap SE and CI of rotated loadings, factor correlations and structure coefficients are based on {valid_rot} bootstrap sample{?s}."))
    }

    L_rot_se_ci <- .array_se_ci(L_rot_boot, ps)
    Phi_rot_se_ci <- .array_se_ci(Phi_rot_boot, ps)
    Structure_se_ci <- .array_se_ci(Structure_boot, ps)

    out <- list(
      SE = list(
        unrot_loadings = L_unrot_se_ci$se,
        rot_loadings = L_rot_se_ci$se,
        Phi = Phi_rot_se_ci$se,
        Structure = Structure_se_ci$se,
        fit_indices = gof_se_ci$se,
        residuals = residuals_se_ci$se,
        valid_target_rotations = valid_rot
      ),
      CI = list(
        unrot_loadings = L_unrot_se_ci$ci,
        rot_loadings = L_rot_se_ci$ci,
        Phi = Phi_rot_se_ci$ci,
        Structure = Structure_se_ci$ci,
        fit_indices = gof_se_ci$ci,
        residuals = residuals_se_ci$ci
      ),
      arrays = list(
        unrot_loadings = L_unrot_boot,
        rot_loadings = L_rot_boot,
        Phi = Phi_rot_boot,
        Structure = Structure_boot,
        fit_indices = gof_boot,
        residuals = residuals_boot
      )
    )

  } else if (boot_rot == "orthogonal") {
    colnam_L <- colnames(L_rot)

    L_rot_boot <- array(NA_real_, c(nrow_L, ncol_L, b),
                        dimnames = list(rownam_L, colnam_L,
                                        NULL))

    failed_rot <- 0

    for (boot_i in seq_len(b)) {

      if (failed[boot_i]) next

      # save target-rotated loading matrix
      aligned_i <- tryCatch(
        PROCRUSTES(boot_fit[[boot_i]]$unrot_loadings,
                   Target = L_rot, rotation = "orthogonal"),
        error = function(e) NULL
      )
      if (is.null(aligned_i)) {
        failed_rot <- failed_rot + 1
        next
      }
      L_rot_boot[,, boot_i] <- aligned_i$loadings
    }

    valid_rot <- b - n_failed - failed_rot
    if (failed_rot > 0) {
      cli::cli_warn(c("{failed_rot} target rotation{?s} in the bootstrap procedure could not be aligned.",
                    "i" = "Bootstrap SE and CI of rotated loadings are based on {valid_rot} bootstrap sample{?s}."))
    }

    L_rot_se_ci <- .array_se_ci(L_rot_boot, ps)

    out <- list(
      SE = list(
        unrot_loadings = L_unrot_se_ci$se,
        rot_loadings = L_rot_se_ci$se,
        fit_indices = gof_se_ci$se,
        residuals = residuals_se_ci$se,
        valid_target_rotations = valid_rot
      ),
      CI = list(
        unrot_loadings = L_unrot_se_ci$ci,
        rot_loadings = L_rot_se_ci$ci,
        fit_indices = gof_se_ci$ci,
        residuals = residuals_se_ci$ci
      ),
      arrays = list(
        unrot_loadings = L_unrot_boot,
        rot_loadings = L_rot_boot,
        fit_indices = gof_boot,
        residuals = residuals_boot
      )
    )


  } else {

    out <- list(
      SE = list(
        unrot_loadings = L_unrot_se_ci$se,
        fit_indices = gof_se_ci$se,
        residuals = residuals_se_ci$se
      ),
      CI = list(
        unrot_loadings = L_unrot_se_ci$ci,
        fit_indices = gof_se_ci$ci,
        residuals = residuals_se_ci$ci
      ),
      arrays = list(
        unrot_loadings = L_unrot_boot,
        fit_indices = gof_boot,
        residuals = residuals_boot
      )
    )

  }


  out

}

.array_se_ci <- function(x, probs, M = c(1, 2)) {
  se <- apply(x, M, stats::sd, na.rm = TRUE)
  ci <- lapply(probs, function(p) {
    apply(x, M, stats::quantile, probs = p, na.rm = TRUE)
  })
  names(ci) <- c("lower", "upper")

  list(
    se = se,
    ci = ci
  )
}
