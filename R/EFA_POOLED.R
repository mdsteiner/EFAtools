#' Exploratory factor analysis on multiple data imputations
#'
#' @author Andreas Soteriades, Markus Steiner
#'
#' @description
#' Fits [EFA()] to each of several imputed datasets, aligns the
#' factor solutions to a common factor space, and pools the resulting estimates
#' and selected fit quantities across imputations.
#'
#' @details
#' `EFA_POOLED()` is the multiple-imputation route to handling missing data:
#' several imputed datasets are each fitted and the solutions pooled. A
#' single-fit alternative is full-information maximum likelihood, available
#' directly in [EFA()] as `cor_method = "fiml"`, which EM-estimates a two-stage
#' correlation from one raw dataset with missing values and analyses it in a
#' single fit. Both feed the same correlation-scale EFA core; they differ only in
#' how the missingness is handled (pooling across imputations here, a single EM
#' correlation there). FIML is intentionally not routed through `EFA_POOLED()`,
#' which is a multi-fit pooler by construction.
#'
#' ## Standard-error pooling routes
#'
#' The standard-error pooling pathway is selected automatically from the `se`
#' method recorded on the component [EFA()] fits, which must be identical across
#' imputations:
#'
#' - `se = "none"`: no standard errors are pooled.
#' - `se = "information"`: the per-imputation expected-information standard
#'   errors are pooled with Rubin's rules (Wald intervals, plain Rubin (1987)
#'   degrees of freedom).
#' - `se = "sandwich"`: the two-stage pooled-inputs (MI2S) approach fits a
#'   single model on the Rubin-pooled correlation matrix and asymptotic
#'   covariance.
#' - `se = "np-boot"`: the non-parametric bootstrap replicates are re-aligned to
#'   the multiple-imputation target and Rubin-pooled.
#'
#' Component fits that mix `se` methods abort with `efa_pooled_mixed_se`. On the
#' information and np-boot routes, a non-trivial failure to produce pooled
#' standard errors (for example an unreliable analytic covariance or too few
#' bootstrap replicates) falls back to point-estimate-only pooling, downgrades
#' `settings$se` to `"none"`, and signals the classed warning
#' `efa_pooled_se_unavailable` (see the Conditions section). The MI2S route is
#' the exception: its single fit fuses the point estimates and standard errors
#' through the pooled asymptotic covariance, so a structural failure such as a
#' non-positive-definite pooled covariance aborts directly rather than falling
#' back.
#'
#' ## Aligning solutions across imputations
#'
#' The function first fits the same [EFA()] model to each imputed
#' dataset. Unrotated loading matrices can optionally be put into a common
#' signed/permuted factor order before averaging. Rotated loading matrices are
#' aligned with either a consensus Procrustes target or with the first
#' imputation's rotated solution as a fixed target. For oblique solutions,
#' factor intercorrelations are transformed/aligned together with the loading
#' matrices so that the factor model remains internally consistent.
#'
#' `target_method` controls how rotated solutions are aligned. `"first_target"` (the
#' default) aligns all imputations to the first imputation's rotated solution via one
#' Procrustes rotation per imputation. `"consensus"` instead refines a centroid target
#' via Generalized Procrustes Analysis (Gower 1975; van Ginkel & Kroonenberg 2014;
#' Lorenzo-Seva & Van Ginkel 2016), iteratively re-Procrustes-rotating all solutions
#' toward it and re-averaging until the target stabilises. For orthogonal rotations the
#' two methods converge to the same pooled estimate (consensus is just a more expensive
#' route to the same answer); for oblique rotations `"consensus"` is not supported and
#' aborts with a classed condition (`efa_consensus_oblique_unsupported`), because the
#' centroid iteration is degenerate for oblique transforms with more than one factor.
#' Anchoring on the first imputation can understate the genuine between-imputation
#' variability when the imputations disagree substantially; `"consensus"` (orthogonal
#' rotations) is more robust to an atypical first imputation (van Ginkel & Kroonenberg
#' 2014).
#'
#' `align_unrotated` controls how unrotated loadings are aligned before pooling.
#' `"signed_tucker_congruence"` (the default) preserves the unrotated axes up to factor
#' reordering and sign changes using Tucker congruence. `"procrustes"` aligns the
#' unrotated matrices to the first imputation by orthogonal Procrustes rotation. `"none"`
#' averages unrotated loadings as returned by `EFA`.
#'
#' ## Pooling point estimates
#'
#' Point estimates are pooled by arithmetic averaging after alignment. For
#' oblique rotations, the returned structure matrix is computed from the pooled
#' aligned pattern matrix and the pooled factor correlation matrix,
#' \eqn{Structure = \Lambda \Phi}. Communalities are always computed as the
#' diagonal of the common-factor reproduced correlation matrix,
#' \eqn{diag(\Lambda \Phi \Lambda')} for oblique rotations and
#' \eqn{diag(\Lambda \Lambda')} otherwise.
#'
#' Residuals are not averaged from the per-imputation residual matrices. Instead,
#' the observed correlation matrices are averaged across imputations and residuals
#' are calculated from the pooled observed correlation matrix minus the
#' model-implied correlation matrix of the pooled solution. Consequently,
#' residual-based fit indices such as RMSR/SRMR are based on pooled residuals.
#' Both are returned for programmatic use, but the print and summary methods show
#' SRMR only because it is the more common reported residual fit index and differs
#' from RMSR only by a fixed scaling for a fixed number of variables.
#'
#' ## Pooling the model chi-square and fit indices
#'
#' The model chi-square and the indices derived from it (RMSEA, ECVI, and the
#' descriptive AIC/BIC) are not arithmetic means of the per-imputation values: the
#' complete-data chi-squares are pooled with the D2 rule and the asymptotic
#' chi-square approximation to D2 is used downstream. The reported pooled
#' chi-square is this asymptotic statistic \eqn{df \cdot F_{D2}}{df * F_D2}, but
#' its p-value is the D2 reference F-test tail
#' \eqn{\Pr(F_{df_1, df_2} > F_{D2})}{P(F(df1, df2) > F_D2)}, not
#' \eqn{\Pr(\chi^2_{df} > df \cdot F_{D2})}{P(chi^2_df > df * F_D2)}; the two
#' coincide only as the D2
#' denominator degrees of freedom grow (no between-imputation variance), and the
#' printed fit block flags the chi-square line accordingly. Because D2 deflates the
#' pooled chi-square by its between-imputation average relative increase in
#' variance, the pooled RMSEA inherits that shrinkage and can fall below the mean
#' of the per-imputation RMSEAs (as it does in `lavaan.mi`); read it together with
#' the per-imputation fit. The incremental indices CFI
#' (Bentler, 1990) and TLI (Tucker & Lewis, 1973) are reported instead as the
#' average of the per-imputation indices, which keeps them consistent with the
#' component fits (CFI stays in \[0, 1\]; TLI, like a single `EFA()` fit, is left
#' unclamped). Pooling the model and baseline chi-squares separately, as
#' `lavaan.mi`/`semTools` do, shrinks the two noncentralities unequally and can
#' drive these non-linear indices outside \[0, 1\]. The separately pooled
#' noncentralities used by that reference remain available in `mi_diagnostics`.
#' AIC and BIC, if returned, are chi-square-derived descriptive quantities based on
#' the D2 approximation and should not be interpreted as likelihood-based MI
#' information criteria. The chi-square pooling rule is D2 (Li, Meng,
#' Raghunathan & Rubin, 1991) on the information and bootstrap routes and the
#' single fit's scaled chi-square on the sandwich/MI2S route.
#'
#' ## Bootstrap pooling (np-boot)
#'
#' If each component `EFA` call was run with `se = "np-boot"` and
#' returned `replicates`, pooled bootstrap SEs and Wald-type MI confidence
#' intervals are computed for loadings, communalities, residuals, and, when
#' applicable, factor correlations and structure coefficients. The
#' rotated bootstrap loading matrices stored by the component `EFA` calls
#' are not reused directly, because they were aligned to imputation-specific
#' targets. Instead, the unrotated bootstrap loading matrices are re-aligned to
#' the final MI target before estimating within-imputation bootstrap covariance
#' matrices. Rubin-type MI pooling is then applied with
#' \eqn{T = Ubar + (1 + 1 / m) B}. Confidence intervals for loadings, Phi,
#' communalities, residuals, and structure coefficients are Wald-type MI
#' intervals. The confidence level of these pooled intervals is controlled by
#' `p`; the `ci` argument passed through `...` to the component
#' `EFA` calls is not used for the pooled intervals.
#' `CI$fit_indices_descriptive`, when available, holds Rubin-Wald MI
#' summaries of the per-imputation bootstrap fit indices. The chi-square-based
#' fit indices themselves are pooled with the D2 rule, whose reference
#' distribution (e.g. the RMSEA confidence interval) supplies their uncertainty.
#'
#' ## Analytic pooling (information)
#'
#' When the component `EFA` calls were run with
#' `se = "information"`, the unrotated-loading and
#' uniqueness SE matrices that those calls return analytically are pooled
#' element-wise with Rubin's rules. The within-imputation variance \eqn{U_d} is
#' taken directly from \eqn{SE_d^2}, the between-imputation variance \eqn{B} is
#' the sample variance of the per-imputation point estimates (after column
#' alignment), and \eqn{T = Ubar + (1 + 1/m) B}. Wald confidence intervals use
#' the plain Rubin (1987) degrees of freedom \eqn{\nu = (m - 1)(1 + 1/r)^2} (with
#' \eqn{r} the relative increase in variance): the analytic loadings are
#' asymptotically normal, so their complete-data reference is normal
#' (\eqn{\nu_{com} = \infty}) and the Barnard-Rubin (1999) small-sample adjustment
#' reduces to this form, matching `lavaan.mi`. NA propagation is fail-closed: if
#' any imputation
#' carries NA at an element, all pooled outputs (SE, CI bounds, RIV, FMI, df)
#' for that element are NA. When a rotation was requested, the rotated loadings,
#' communalities, and (for oblique rotations) factor correlations and structure
#' coefficients are pooled as well (see below); residual SE pooling remains
#' available only on the bootstrap path. Under
#' `align_unrotated = "procrustes"`, the per-imputation orthogonal transform
#' \eqn{Q_d} mixes loading columns, so marginal SEs alone no longer determine
#' the aligned marginal SE matrix. The full unrotated covariance persisted on
#' each fit (`vcov_unrot_loadings`, populated by `se = "information"`) is
#' propagated via the column-major Kronecker identity
#' \deqn{Var(vec(L_d Q_d)) = (Q_d' \otimes I_p) V_d (Q_d \otimes I_p),}
#' from which the aligned row-marginal variances are read off and pooled.
#' \eqn{Q_d} is treated as fixed (Schoenemann 1966; mirrors the lavaan / mice
#' convention of not differentiating the per-imputation Procrustes solution).
#' If any component fit lacks `vcov_unrot_loadings`, the pool aborts with the
#' classed condition `efa_pooled_no_vcov`; if any component fit carries an
#' NA-filled `vcov_unrot_loadings` (the upstream signal of a Heywood case or
#' singular bordered information matrix), it aborts with the classed condition
#' `efa_pooled_unreliable_vcov`.
#'
#' ### Rotated-loading standard errors
#'
#' Rotated quantities are pooled in a common multiple-imputation rotational
#' gauge. A rotated-loading standard error is conditional on the rotation
#' criterion (Jennrich 1973, 1974; Browne 2001; Zhang & Preacher 2015), so for
#' orthogonal and oblique rotations alike the within-imputation variance is each
#' component fit's own criterion-aware delta-method rotated standard error (the
#' quantity `EFA()` returns), reused after a signed-permutation alignment to the
#' multiple-imputation target; the between-imputation variance is the sample
#' variance of the aligned per-imputation rotated loadings. This shares the
#' estimand of the bootstrap pool, which likewise aligns each replicate to the
#' target, and is flagged by `MI$<param>$method = "signed_permutation_approx"`. It
#' is a deliberate approximation: the per-fit standard error is conditional on
#' that fit's rotation-criterion optimum rather than on the common gauge. A
#' fixed-rotation linear map of the unrotated covariance is deliberately not used
#' for rotated-loading standard errors -- an orthogonal rotation preserves each
#' row's total variance, so it would merely redistribute the rotational
#' indeterminacy variance the criterion otherwise pins down, inflating the
#' reported standard errors. The same signed-permutation alignment carries the
#' factor-correlation and structure-coefficient standard errors for oblique
#' rotations. Communalities are rotation-invariant and pool element-wise with no
#' alignment. Researchers who need a fully gauge-consistent uncertainty for
#' rotated solutions should cross-check with `se = "np-boot"`, which re-aligns the
#' unrotated replicates to the target before estimating the within-imputation
#' covariance; the multiple-imputation rotation literature otherwise summarises
#' rotated uncertainty with generalised Procrustes centroids rather than
#' per-element standard errors (van Ginkel & Kroonenberg 2014; Lorenzo-Seva & Van
#' Ginkel 2016).
#'
#' ## Two-stage pooling (sandwich / MI2S)
#'
#' When the component `EFA` calls were run with `se = "sandwich"` (robust
#' standard errors from a polychoric/tetrachoric or continuous-Pearson
#' asymptotic covariance), pooling instead follows the two-stage,
#' pooled-inputs approach (Chung & Cai 2019; Sriutaisuk, Liu, Chung, Kim & Gu
#' 2025). Rather than fitting per imputation and pooling the estimates, the
#' correlation matrix and the asymptotic covariance of its off-diagonal entries
#' are pooled across imputations with Rubin's (1987) rules,
#' \deqn{\bar r = \frac{1}{m}\sum_d r_d, \qquad
#'       \tilde\Gamma = \Gamma_W + \left(1 + \frac{1}{m}\right)\Gamma_B,}
#' with \eqn{\Gamma_W = \frac{1}{m}\sum_d \Gamma_d} the within-imputation
#' covariance (the average of the per-imputation \eqn{\Gamma_d}) and
#' \eqn{\Gamma_B} the between-imputation covariance of the off-diagonal
#' correlations (sample covariance, \eqn{m - 1} divisor). A single `EFA` model
#' is then fitted to \eqn{\bar r} with \eqn{\tilde\Gamma} as the robust meat
#' (and, for `method = "DWLS"`, its diagonal as the weights). Because there is
#' only one fit there is a single rotational gauge, so this route bypasses the
#' per-imputation Procrustes/target alignment entirely (`target_method` and
#' `align_unrotated` do not apply; supplying a non-default value raises the classed
#' warning `efa_pooled_mi2s_alignment_ignored`, while the defaults are ignored
#' silently).
#' The fitted object carries native scaled-shifted chi-square test statistics
#' and sandwich standard errors that already reflect the multiple-imputation
#' uncertainty; the chi-square is not D2-pooled, and the
#' likelihood-ratio-based AIC/BIC/ECVI are left `NA` (as for any sandwich fit).
#' The single fit is returned in the `mi_fit` slot, and the per-imputation
#' `fits` are retained for diagnostics. The pooled asymptotic covariance is not
#' guaranteed positive semidefinite at small \eqn{m}; rather than projecting it
#' to the nearest positive-definite matrix (which would silently distort the
#' robust meat), an indefinite \eqn{\tilde\Gamma} aborts with the classed
#' condition `efa_pooled_mi2s_acov_not_psd`, which signals that more
#' imputations are needed. Sriutaisuk et al. (2025) recommend at least 20
#' imputations for the scaled-shifted statistic, and more (around 100) at
#' higher rates of missingness; fewer than 20 raises a classed warning. The
#' polychoric/tetrachoric (ordinal) case is the primary, best-evaluated target;
#' the continuous-Pearson case uses the same recipe but is less benchmarked
#' under multiple imputation.
#'
#' @param data_list A list of length \eqn{m}, where \eqn{m} is the number of
#' imputations. Each list element is a data frame or matrix of raw data, or a
#' correlation matrix. See argument `x` in [EFA()].
#' @param p Numeric in \eqn{(0, 1)}. One minus the confidence level used for
#' pooled Wald-type bootstrap/MI confidence intervals when bootstrap replicates
#' are available. For example, `p = .05` gives 95% intervals.
#' @param target_method Character. How rotated solutions are aligned across imputations
#' before pooling: `"first_target"` (the default) aligns every imputation to the first
#' imputation's rotated solution, while `"consensus"` refines a centroid target by
#' Generalized Procrustes Analysis (orthogonal rotations only). See *Aligning solutions
#' across imputations* in Details.
#' @param align_unrotated Character. How unrotated loadings are aligned before pooling:
#' `"signed_tucker_congruence"` (the default; sign/permutation via Tucker congruence),
#' `"procrustes"` (orthogonal Procrustes to the first imputation), or `"none"`. See
#' *Aligning solutions across imputations* in Details.
#' @param fit_pool_method Character. Currently only `"D2"` is implemented
#' for chi-square-type fit. If no chi-square is available, only residual-based
#' fit and descriptive quantities are returned. See *Pooling the model chi-square and
#' fit indices* in Details.
#' @param consensus_args List of additional arguments controlling the
#' GPA-consensus iteration when `target_method = "consensus"`. See the
#' source of `.gpa_consensus_target` for the available tuning parameters
#' (convergence tolerance, maximum iterations, multi-start options).
#' @param procrustes_args List of additional arguments passed to `PROCRUSTES`
#' for fixed-target alignment.
#' @param rmsea_ci_level Numeric. Confidence level for the RMSEA CI.
#' @param rmsr_upper Logical. If `TRUE`, compute RMSR from the unique
#' off-diagonal residual correlations. If `FALSE`, use the full off-diagonal
#' matrix.
#' @param ... Additional arguments passed to [EFA()] (e.g. `method`, `rotation`, `se`,
#' `n_factors`, `N`). These select the estimator, rotation, standard-error method, and
#' fit indices used for every imputation; see [EFA()] for the available options, their
#' properties, and which combinations are valid.
#'
#' @return A list of class `"EFA_POOLED"` containing pooled estimates,
#' residuals, fit indices, the individual fits, and MI diagnostics. In
#' addition to the slots inherited from [EFA()] (including `SE`, `CI`, and,
#' on the bootstrap path, `replicates`), the object carries:
#' \describe{
#' \item{MI}{Multiple-imputation diagnostics for each pooled parameter family.
#' On the bootstrap path: `unrot_loadings`, `h2`, `residuals`, optionally
#' `rot_loadings`, `Phi`, `Structure`, and `fit_indices_descriptive`, plus
#' integer vectors `bootstrap_source_failures` (replicates the component `EFA`
#' could not fit), `bootstrap_rotation_failures` (replicates whose Procrustes
#' alignment to the target was invalid), and `bootstrap_rotation_valid` (those
#' that entered the pool, `B - source - rotation` failures). Both paths use the
#' plain Rubin (1987) df. On the analytic path (`se = "information"`):
#' `unrot_loadings` and `uniquenesses`, plus, when a rotation was requested,
#' `rot_loadings`, `h2`, and (oblique) `Phi` and `Structure`. Each per-family
#' entry is a list with `RIV` (relative increase in variance), `FMI` (the
#' fraction of missing information, reported as Rubin's asymptotic
#' \eqn{\lambda = RIV / (1 + RIV)}, equal to `lavaan.mi`'s `fmi`), and `df`; the
#' rotated families on the analytic path additionally carry a `method` string
#' recording the gauge alignment used (`"gauge_invariant"` for communalities and
#' `"signed_permutation_approx"` for rotated loadings and, for oblique rotations,
#' factor correlations and structure coefficients). `fit_indices_descriptive`, on
#' the bootstrap path, pools every per-imputation fit index, so the structural
#' constants among them (`df`, `df_null`) appear with a standard error of 0.}
#' \item{mi_fit}{On the `se = "sandwich"` (MI2S) path only: the single [EFA()]
#' fit on the pooled correlation matrix \eqn{\bar r} and pooled asymptotic
#' covariance \eqn{\tilde\Gamma}. Its `orig_R` is \eqn{\bar r} and its `Gamma`
#' is \eqn{\tilde\Gamma}; the pooled `SE`, `CI`, and `fit_indices` are taken
#' from it. `MI` is `NULL` on this path because the imputation uncertainty is
#' carried by \eqn{\tilde\Gamma} rather than by per-parameter Rubin pooling.}
#' }
#'
#'
#' @references
#' Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
#' multiple imputation. *Biometrika*, 86(4), 948-955.
#'
#' Bentler, P. M. (1990). Comparative fit indexes in structural models.
#' *Psychological Bulletin*, 107(2), 238-246.
#'
#' Browne, M. W. (2001). An overview of analytic rotation in exploratory factor
#' analysis. *Multivariate Behavioral Research*, 36(1), 111-150.
#'
#' Chan, K. W., & Meng, X.-L. (2022). Multiple improvements of multiple
#' imputation likelihood ratio tests. *Statistica Sinica*, 32, 1489-1514.
#'
#' Chung, S., & Cai, L. (2019). Alternative multiple imputation inference for
#' categorical structural equation modeling. *Multivariate Behavioral
#' Research*, 54(3), 323-337.
#'
#' Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika*,
#' 40(1), 33-51.
#'
#' Li, K. H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
#' Significance levels from repeated p-values with multiply-imputed data.
#' *Statistica Sinica*, 1(1), 65-92.
#'
#' Jennrich, R. I. (1973). Standard errors for obliquely rotated factor
#' loadings. *Psychometrika*, 38(4), 593-604.
#'
#' Jennrich, R. I. (1974). Simplified formulae for standard errors in
#' maximum-likelihood factor analysis. *British Journal of Mathematical and
#' Statistical Psychology*, 27(1), 122-131.
#'
#' Lorenzo-Seva, U., & Van Ginkel, J. R. (2016). Multiple imputation of missing
#' values in exploratory factor analysis of multidimensional scales. *Anales de
#' Psicologia*, 32(2), 596-608.
#'
#' Rubin, D. B. (1987). *Multiple imputation for nonresponse in surveys*. Wiley.
#'
#' Schoenemann, P. H. (1966). A generalized solution of the orthogonal
#' Procrustes problem. *Psychometrika*, 31(1), 1-10.
#'
#' Sriutaisuk, S., Liu, Y., Chung, S., Kim, H., & Gu, F. (2025). Evaluating
#' imputation-based fit statistics in structural equation modeling with ordinal
#' data: The MI2S approach. *Educational and Psychological Measurement*, 85(1),
#' 82-113.
#'
#' Tucker, L. R., & Lewis, C. (1973). A reliability coefficient for maximum
#' likelihood factor analysis. *Psychometrika*, 38(1), 1-10.
#'
#' van Ginkel, J. R., & Kroonenberg, P. M. (2014). Using generalized Procrustes
#' analysis for multiple imputation in principal component analysis. *Journal of
#' Classification*, 31(2), 242-269.
#'
#' Zhang, G., & Preacher, K. J. (2015). Factor rotation and standard errors in
#' exploratory factor analysis. *Journal of Educational and Behavioral
#' Statistics*, 40(6), 579-603.
#'
#' @export
#'
#' @examples
#'
#' # create a list of three datasets, mimicking a list you would obtain from
#' # e.g. mice.
#' dat_list <- lapply(1:3, function(x) GRiPS_raw[sample(1:nrow(GRiPS_raw), replace = TRUE),])
#' mod <- EFA_POOLED(dat_list, n_factors = 1, method = "ML")
#' mod
#'
#' \donttest{
#' # add computation of standard errors and CIs
#' mod <- EFA_POOLED(dat_list, n_factors = 1, method = "ML", se = "np-boot")
#' mod
#' }
EFA_POOLED <- function(data_list,
                       p = 0.05,
                       target_method = c("first_target", "consensus"),
                       align_unrotated = c("signed_tucker_congruence", "none", "procrustes"),
                       fit_pool_method = c("D2"),
                       consensus_args = list(),
                       procrustes_args = list(),
                       rmsea_ci_level = .90,
                       rmsr_upper = TRUE,
                       ...) {

  efa_args <- list(...)

  checkmate::assert_list(data_list, min.len = 2, null.ok = FALSE)
  lapply(data_list, checkmate::assert_multi_class, c("matrix", "data.frame"))
  checkmate::assert_number(p, na.ok = FALSE, lower = 0, upper = 1)
  checkmate::assert_list(consensus_args, null.ok = FALSE)
  checkmate::assert_list(procrustes_args, null.ok = FALSE)
  checkmate::assert_number(rmsea_ci_level, na.ok = FALSE, lower = 0, upper = 1)
  checkmate::assert_flag(rmsr_upper)

  if (p <= 0 || p >= 1) {
    cli::cli_abort("{.arg p} must be strictly between 0 and 1.", class = "efa_pooled_bad_p")
  }
  if (rmsea_ci_level <= 0 || rmsea_ci_level >= 1) {
    cli::cli_abort("{.arg rmsea_ci_level} must be strictly between 0 and 1.", class = "efa_pooled_bad_ci_level")
  }
  if (!is.null(efa_args$ci) && length(efa_args$ci) == 1L &&
      is.finite(efa_args$ci) &&
      abs(as.numeric(efa_args$ci) - (1 - p)) > sqrt(.Machine$double.eps)) {
    cli::cli_warn("{.fn EFA_POOLED} uses {.arg p}, not the component {.fn EFA} argument {.arg ci}, to set pooled bootstrap/MI confidence intervals.",
                  class = "efa_pooled_ci_ignored")
  }

  target_method <- match.arg(target_method)
  align_unrotated <- match.arg(align_unrotated)
  fit_pool_method <- match.arg(fit_pool_method)

  m_imp <- length(data_list)

  ## -------------------------------------------------------------------------
  ## Fit EFA to each imputed dataset
  ## -------------------------------------------------------------------------

  fits <- lapply(data_list, function(data_list_subset) {
    do.call(EFA, c(list(x = data_list_subset), efa_args))
  })

  .efa_pooled_check_fits(fits)

  # Select the SE-pooling pathway from the component fits' shared se method. A
  # heterogeneous se cannot be pooled into a single MI estimate, so abort rather
  # than silently producing an uninterpretable mixture.
  route <- .efa_pooled_route(fits)
  if (identical(route, "mixed")) {
    cli::cli_abort(
      c("The component {.fn EFA} fits use different {.arg se} methods, so their standard errors cannot be pooled.",
        "i" = "Re-fit every imputation with the same {.arg se} (all {.val none}, {.val information}, {.val sandwich}, or {.val np-boot})."),
      class = "efa_pooled_mixed_se"
    )
  }

  settings <- fits[[1]]$settings
  method <- settings$method
  rotation <- settings$rotation
  var_names <- rownames(fits[[1]]$orig_R)
  if (is.null(var_names)) {
    var_names <- colnames(data_list[[1]])
  }

  rotation_type <- if (is.null(rotation)) "none" else .rotation_family(rotation)

  ## -------------------------------------------------------------------------
  ## MI2S route (Chung & Cai 2019; Sriutaisuk et al. 2025): when the component
  ## fits carry sandwich (robust) standard errors, pool the correlation matrix
  ## and the asymptotic covariance of its off-diagonal entries across imputations
  ## and fit a single model on the pooled inputs. The single fit yields native
  ## scaled-shifted chi-square and sandwich SEs that reflect the multiple-
  ## imputation uncertainty, so this path bypasses the per-imputation rotation
  ## alignment and Rubin pooling used by the bootstrap/information routes below.
  ## -------------------------------------------------------------------------
  if (identical(route, "sandwich")) {
    return(.efa_pooled_mi2s(
      fits = fits, data_list = data_list, efa_args = efa_args,
      settings = settings, method = method, rotation = rotation,
      rotation_type = rotation_type, target_method = target_method,
      align_unrotated = align_unrotated, fit_pool_method = fit_pool_method,
      p = p, rmsea_ci_level = rmsea_ci_level, rmsr_upper = rmsr_upper
    ))
  }

  ## -------------------------------------------------------------------------
  ## Extract and align unrotated loadings
  ## -------------------------------------------------------------------------

  unrot_loadings <- lapply(.extract_list_object(fits, "unrot_loadings"),
                           .change_class, "matrix")

  aligned_unrot <- .efa_pooled_align_unrotated_list(
    unrot_loadings = unrot_loadings,
    align_unrotated = align_unrotated,
    return_meta = TRUE
  )
  unrot_loadings_aligned <- aligned_unrot$loadings
  unrot_align_meta       <- aligned_unrot$meta

  mean_unrot_loadings <- .average_matrices(unrot_loadings_aligned)

  ## -------------------------------------------------------------------------
  ## Align rotated loadings and Phi, if a rotation was requested
  ## -------------------------------------------------------------------------

  # A single-factor solution cannot be rotated and has no factor correlations, so
  # an oblique request on one factor is treated as the no-Phi case throughout,
  # matching a single-fit EFA() (which returns Phi = NULL for one factor).
  oblique <- rotation_type == "oblique" && ncol(unrot_loadings_aligned[[1]]) >= 2L

  phis <- NULL
  structure_loadings <- NULL

  if (rotation_type != "none") {
    rot_loadings_initial <- lapply(.extract_list_object(fits, "rot_loadings"),
                                   .change_class, "matrix")

    if (target_method == "first_target") {
      target_rotations <- vector("list", m_imp)
      rot_loadings <- vector("list", m_imp)
      phis <- vector("list", m_imp)
      # One factor is aligned by sign only; an oblique Procrustes is neither
      # needed nor well defined for a single column.
      proc_rotation <- if (oblique) "oblique" else "orthogonal"

      rot_loadings[[1]] <- rot_loadings_initial[[1]]
      if (oblique) {
        phis[[1]] <- fits[[1]]$Phi
      }

      point_rotation_failures <- logical(m_imp)
      for (d in 2:m_imp) {
        target_rotations[[d]] <- do.call(
          PROCRUSTES,
          c(list(A = unrot_loadings[[d]],
                 Target = rot_loadings_initial[[1]],
                 rotation = proc_rotation),
            procrustes_args)
        )
        point_rotation_failures[d] <- isFALSE(target_rotations[[d]]$valid)
        rot_loadings[[d]] <- target_rotations[[d]]$loadings
        if (oblique) {
          phis[[d]] <- target_rotations[[d]]$Phi
        }
      }

      if (any(point_rotation_failures)) {
        cli::cli_warn(
          c("At least one fixed-target Procrustes alignment could not be aligned to a valid rotation.",
            "i" = "The pooled point estimates still use the best available {.fn PROCRUSTES} alignment; inspect {.code alignment$point_rotation_failures}."),
          class = "efa_pooled_align_failed"
        )
      }

      final_target <- rot_loadings_initial[[1]]
      inner_converged <- vapply(target_rotations[-1L], function(x) {
        isTRUE(x$valid) && (is.null(x$convergence) || isTRUE(x$convergence))
      }, logical(1L))
      alignment <- list(method = "first_target",
                        target = final_target,
                        target_rotations = target_rotations,
                        point_rotation_failures = which(point_rotation_failures),
                        converged = all(inner_converged))

    } else if (target_method == "consensus") {
      consensus <- do.call(
        .gpa_consensus_target,
        c(list(unrotated_list = unrot_loadings,
               init_targets = rot_loadings_initial,
               rotation = rotation_type),
          consensus_args)
      )

      if (!isTRUE(consensus$converged)) {
        cli::cli_warn(
          c("Consensus Procrustes alignment did not meet its convergence criterion.",
            "i" = "Inspect {.code alignment$history} and consider stricter or multi-start {.arg consensus_args}."),
          class = "efa_pooled_align_failed"
        )
      }

      rot_loadings <- consensus$aligned_loadings
      phis <- consensus$aligned_phi
      final_target <- consensus$target
      alignment <- consensus
    }

    mean_rot_loadings <- .average_matrices(rot_loadings)
    rownames(mean_rot_loadings) <- var_names

    if (oblique) {
      mean_phis <- .average_matrices(phis)
      # Ensure mean_phis is symmetric (needed, as due to floating point imprecision,
      # the averaged matrix can be asymmetric, leading to issues in downstream checks
      # of matrix symmetry)
      mean_phis <- (mean_phis + t(mean_phis)) / 2
      structure_loadings <- Map(function(L, Phi) L %*% Phi, rot_loadings, phis)
      # Keep Structure parallel to the returned pooled pattern matrix and Phi:
      # it is the plug-in structure of the pooled solution, not the arithmetic
      # mean of the imputation-specific structure matrices.
      mean_structure_loadings <- mean_rot_loadings %*% mean_phis
      dimnames(mean_structure_loadings) <- dimnames(mean_rot_loadings)
    } else {
      mean_phis <- NULL
      mean_structure_loadings <- NULL
    }

  } else {
    rot_loadings <- NULL
    mean_rot_loadings <- NULL
    mean_phis <- NULL
    mean_structure_loadings <- NULL
    final_target <- NULL
    alignment <- NULL
  }

  ## -------------------------------------------------------------------------
  ## Pooled observed and model-implied correlation matrices; residuals and RMSR
  ## -------------------------------------------------------------------------

  orig_R_list <- .extract_list_object(fits, "orig_R")

  pooled_orig_R <- .average_matrices(orig_R_list)

  if (oblique) {
    common_R <- mean_rot_loadings %*% mean_phis %*% t(mean_rot_loadings)
  } else if (rotation_type != "none") {
    common_R <- mean_rot_loadings %*% t(mean_rot_loadings)
  } else {
    common_R <- mean_unrot_loadings %*% t(mean_unrot_loadings)
  }

  h2 <- diag(common_R)
  names(h2) <- rownames(common_R)

  model_implied_R <- common_R
  diag(model_implied_R) <- 1
  dimnames(model_implied_R) <- list(var_names, var_names)

  residuals <- pooled_orig_R - model_implied_R
  diag(residuals) <- 0
  dimnames(residuals) <- list(var_names, var_names)

  RMSR <- .rmsr(residuals, upper = rmsr_upper)

  ## -------------------------------------------------------------------------
  ## Variance-accounted tables
  ## -------------------------------------------------------------------------

  mean_vars_accounted <- .compute_vars(
    L_unrot = mean_unrot_loadings,
    L_rot = mean_unrot_loadings,
    Phi = NULL
  )

  if (oblique) {
    mean_vars_accounted_rot <- .compute_vars(
      L_unrot = mean_unrot_loadings,
      L_rot = mean_rot_loadings,
      Phi = mean_phis
    )
  } else if (rotation_type != "none") {
    mean_vars_accounted_rot <- .compute_vars(
      L_unrot = mean_unrot_loadings,
      L_rot = mean_rot_loadings,
      Phi = NULL
    )
  } else {
    mean_vars_accounted_rot <- NULL
  }

  ## -------------------------------------------------------------------------
  ## Fit indices: residual-based from pooled residuals, chi-square-based via D2
  ## -------------------------------------------------------------------------

  Ns <- .efa_pooled_get_Ns(data_list, fits, efa_args)
  Ns_ok <- Ns[is.finite(Ns)]
  if (length(Ns_ok) == 0L) {
    N_pool <- NA_real_
    if (!identical(method, "PAF")) {
      cli::cli_warn(
        c("N could not be recovered for any imputation.",
          "i" = "Chi-square-based fit indices cannot be computed."),
        class = "efa_pooled_no_n"
      )
    }
  } else {
    if (length(Ns_ok) < length(Ns)) {
      cli::cli_warn(
        c("N could not be recovered for every imputation.",
          "i" = "Fit indices use the mean of the available Ns."),
        class = "efa_pooled_partial_n"
      )
    }
    if (length(unique(Ns_ok)) > 1L) {
      cli::cli_warn(
        c("The imputed datasets appear to have different N.",
          "i" = "Fit indices use the mean N across imputations."),
        class = "efa_pooled_unequal_n"
      )
    }
    N_pool <- mean(Ns_ok)
  }

  fit_indices <- .efa_pooled_fit_indices(
    fits = fits,
    pooled_R = pooled_orig_R,
    residuals = residuals,
    RMSR = RMSR,
    N = N_pool,
    Ns = Ns,
    method = method,
    pool_method = fit_pool_method,
    rmsea_ci_level = rmsea_ci_level
  )

  ## -------------------------------------------------------------------------
  ## SEs and CIs for pooled MI estimates. Two routes: a bootstrap pool that
  ## requires component-level replicate cubes, and an analytic pool that
  ## consumes the per-imputation closed-form SE matrices populated by
  ## se = "information" via Rubin's (1987) rules.
  ## -------------------------------------------------------------------------

  # The route fixes which pool runs: "information" consumes the per-imputation
  # closed-form SE matrices, "np-boot" the replicate cubes. A pooling failure on
  # either route is non-fatal and falls back to point-estimate-only pooling,
  # signalling the umbrella `efa_pooled_se_unavailable` condition. The "sandwich"
  # route returned via .efa_pooled_mi2s() above, and "none" requests no SEs at all.
  boot_pooled <- NULL
  se_failed <- FALSE
  # Shared fallback: a classed analytic-pool abort becomes a point-estimate-only
  # solution flagged with the umbrella condition. `route` (captured here) resolves
  # to the active route inside each branch, so it labels the warning correctly.
  se_unavailable_handler <- function(e) {
    se_failed <<- TRUE
    .efa_pooled_se_unavailable(route, parent = e)
    NULL
  }
  if (identical(route, "information")) {
    # Catch only the analytic pool's own classed aborts (an unreliable/absent
    # unrotated covariance, or missing alignment metadata) and fall back to
    # point-estimate-only pooling; an unexpected error is a real bug and must
    # propagate rather than be silently downgraded to "se unavailable".
    boot_pooled <- tryCatch(
      .efa_pooled_analytic_pool(
        fits = fits,
        unrot_loadings_aligned = unrot_loadings_aligned,
        align_meta = unrot_align_meta,
        ci = 1 - p,
        align_unrotated = align_unrotated,
        rotation_type = rotation_type,
        rot_loadings = rot_loadings,
        phis = phis,
        structure_loadings = structure_loadings,
        mean_structure_loadings = mean_structure_loadings,
        mean_phis = mean_phis,
        h2 = h2
      ),
      efa_pooled_unreliable_vcov = se_unavailable_handler,
      efa_pooled_no_vcov = se_unavailable_handler,
      efa_pooled_analytic_align_meta_missing = se_unavailable_handler
    )
  } else if (identical(route, "np-boot")) {
    # The bootstrap pool soft-returns NULL (with its own classed warning) when it
    # cannot pool; that NULL is handled by the fallback check below. It has no
    # classed aborts of its own, so call it directly and let any unexpected error
    # propagate rather than masking it as "se unavailable".
    boot_pooled <- .efa_pooled_bootstrap_pool(
      fits = fits,
      orig_R_list = orig_R_list,
      unrot_loadings_aligned = unrot_loadings_aligned,
      mean_unrot_loadings = mean_unrot_loadings,
      rot_loadings = rot_loadings,
      phis = phis,
      structure_loadings = structure_loadings,
      mean_structure_loadings = mean_structure_loadings,
      final_target = final_target,
      rotation_type = rotation_type,
      align_unrotated = align_unrotated,
      procrustes_args = procrustes_args,
      h2 = h2,
      residuals = residuals,
      pooled_orig_R = pooled_orig_R,
      N = N_pool,
      method = method,
      pool_method = fit_pool_method,
      rmsea_ci_level = rmsea_ci_level,
      alpha = p,
      rmsr_upper = rmsr_upper
    )
  }

  # A soft failure (the helper returned NULL after emitting its own specific
  # condition, e.g. too few bootstrap replicates) still downgrades to no SE
  # pooling; surface the umbrella condition unless the tryCatch handler already
  # did so for a classed abort.
  if ((identical(route, "information") || identical(route, "np-boot")) &&
      is.null(boot_pooled) && !se_failed) {
    .efa_pooled_se_unavailable(route)
  }

  ## -------------------------------------------------------------------------
  ## Return object
  ## -------------------------------------------------------------------------

  mean_unrot_loadings <- .change_class(mean_unrot_loadings, "LOADINGS")
  if (!is.null(mean_rot_loadings)) {
    mean_rot_loadings <- .change_class(mean_rot_loadings, "LOADINGS")
  }
  if (!is.null(mean_structure_loadings)) {
    mean_structure_loadings <- .change_class(mean_structure_loadings, "LOADINGS")
  }

  settings_pooled <- settings
  settings_pooled$N <- N_pool
  settings_pooled$pooled_N <- N_pool
  settings_pooled$pooled <- TRUE
  settings_pooled$component_se <- settings$se
  # Downgrade se only on the fallback path: a requested se = "information" /
  # "np-boot" whose pooled SEs could not be produced (an `efa_pooled_se_unavailable`
  # warning was emitted above). The se = "none" route never requested SEs, so it
  # keeps "none" without a warning or downgrade. `component_se` retains the
  # requested method so print.EFA_POOLED() can still describe what the component
  # fits used.
  if (is.null(boot_pooled) && !identical(route, "none")) {
    settings_pooled$se <- "none"
  }
  settings_pooled$n_imputations <- m_imp
  settings_pooled$target_method <- target_method
  settings_pooled$align_unrotated <- align_unrotated
  settings_pooled$fit_pool_method <- fit_pool_method
  settings_pooled$p <- p
  settings_pooled$ci <- 1 - p
  settings_pooled$rmsea_ci_level <- rmsea_ci_level
  settings_pooled$rmsr_upper <- rmsr_upper
  if (!is.null(boot_pooled) && !is.null(boot_pooled$n_boot)) {
    settings_pooled$b_boot <- boot_pooled$n_boot
  }

  results <- list(
    h2 = h2,
    unrot_loadings = mean_unrot_loadings,
    vars_accounted = mean_vars_accounted,
    fit_indices = fit_indices,
    model_implied_R = model_implied_R,
    residuals = residuals,
    orig_R = pooled_orig_R,
    settings = settings_pooled,
    fits = fits,
    alignment = alignment,
    mi_diagnostics = fit_indices$mi_diagnostics
  )

  if (rotation_type != "none") {
    results$rot_loadings <- mean_rot_loadings
    results$vars_accounted_rot <- mean_vars_accounted_rot
  }

  if (oblique) {
    results$Phi <- mean_phis
    results$Structure <- mean_structure_loadings
  }

  if (!is.null(boot_pooled)) {
    results$SE <- boot_pooled$SE
    results$CI <- boot_pooled$CI
    results$MI <- boot_pooled$MI
    # Single-bracket list assignment preserves a present-but-NULL `replicates`
    # slot on the analytic path, matching the EFA() schema contract pinned in
    # test-EFA-fields.R; `results$replicates <- NULL` would remove the slot.
    results["replicates"] <- list(boot_pooled$replicates)
    if (!is.null(boot_pooled$SE$residuals)) {
      # The residual diagonal is fixed at 0 with SE 0, so 0/0 would yield NaN on
      # the diagonal; the off-diagonal standardised residuals are the meaningful
      # quantities, so set the diagonal to 0 (matches single-fit EFA()).
      std_resid <- results$residuals / boot_pooled$SE$residuals
      diag(std_resid) <- 0
      results$standardized_residuals <- std_resid
    }
  }

  class(results) <- c("EFA_POOLED", "EFA")
  results
}

## =============================================================================
## Internal helpers for EFA_POOLED
## =============================================================================

.efa_pooled_setting_chr <- function(f, name) {
  # A component fit's settings value as a length-1 character (NA when absent), for
  # comparing settings across imputations.
  val <- f$settings[[name]]
  if (is.null(val)) NA_character_ else as.character(val)
}

.efa_pooled_route <- function(fits) {
  # Decide which multiple-imputation SE-pooling pathway a set of component fits
  # selects, from the standard-error method recorded on each fit. Pooling is only
  # defined when every imputation used the same se, so a heterogeneous (or
  # absent) se returns "mixed" and the caller aborts. Returns one of "none",
  # "np-boot", "information", "sandwich", or "mixed".
  ses <- vapply(fits, .efa_pooled_setting_chr, character(1), name = "se")
  if (anyNA(ses) || length(unique(ses)) != 1L) {
    return("mixed")
  }
  ses[[1L]]
}

.efa_pooled_se_unavailable <- function(route, parent = NULL) {
  # Umbrella condition for "the requested standard-error method could not be
  # pooled". Emitted whenever EFA_POOLED() falls back to point-estimate-only
  # pooling on the information / np-boot routes, layered over any specific
  # condition the pooling helper raised (carried as `parent` for a classed
  # abort, NULL when the helper soft-returned NULL after its own warning).
  cli::cli_warn(
    c("Pooled standard errors could not be produced for the {.val {route}} standard-error method; the pooled point estimates are returned without standard errors.",
      "i" = "Inspect the component {.fn EFA} fits and any preceding condition for the cause."),
    parent = parent, class = "efa_pooled_se_unavailable"
  )
}

## -----------------------------------------------------------------------------
## MI2S (multiple imputation, two-stage) pooled-inputs route for sandwich SEs
## -----------------------------------------------------------------------------

.efa_pooled_mi2s_inputs <- function(fits) {
  # Pool the per-imputation correlation matrices and their asymptotic covariances
  # into a single correlation matrix r_bar and a single ACOV Gamma_tilde via the
  # two-stage pooled-input ACOV of Chung & Cai (2019) and Sriutaisuk, Liu, Chung,
  # Kim & Gu (2025, Educ. Psychol. Meas.). This is the ACOV-pooling step of MI2S,
  # distinct from the parameter-pooling Rubin's-rules routes (which carry FMI/df):
  #   r_bar       = (1/m) sum_d R_d                         (element-wise mean)
  #   Gamma_W     = (1/m) sum_d Gamma_d                     (within-imputation ACOV)
  #   Gamma_B     = (1/(m-1)) sum_d (v_d - v_bar)(v_d - v_bar)'   (between, = cov())
  #   Gamma_tilde = Gamma_W + (1 + 1/m) Gamma_B             (total pooled ACOV)
  # All quantities stay on the per-fit variance scale (Var(rho-hat)) and in the
  # utils::combn(p, 2) off-diagonal order that EFA()'s $Gamma uses, so the pooled
  # ACOV feeds the sandwich core unchanged.
  m <- length(fits)
  R_list <- lapply(fits, function(f) as.matrix(f$orig_R))
  G_list <- lapply(fits, function(f) as.matrix(f$Gamma))
  p <- ncol(R_list[[1]])

  idx <- utils::combn(p, 2L)
  pair_ij <- cbind(idx[1L, ], idx[2L, ])
  q <- ncol(idx)

  r_bar <- Reduce(`+`, R_list) / m
  r_bar <- (r_bar + t(r_bar)) / 2          # symmetrise away round-off
  diag(r_bar) <- 1

  # Stack the off-diagonal correlations (combn order) for the between-imputation
  # covariance; stats::cov() uses the m - 1 divisor required by Rubin's rule.
  V <- t(vapply(R_list, function(R) R[pair_ij], numeric(q)))

  Gamma_W <- Reduce(`+`, G_list) / m
  Gamma_B <- stats::cov(V)
  Gamma_tilde <- Gamma_W + (1 + 1 / m) * Gamma_B
  Gamma_tilde <- (Gamma_tilde + t(Gamma_tilde)) / 2

  # Gamma_tilde = Gamma_W + (1 + 1/m) Gamma_B is in exact arithmetic a sum of two
  # positive-semidefinite matrices (Gamma_W averages PSD per-fit ACOVs; Gamma_B is
  # a sample covariance), so it cannot become indefinite from a small m alone. A
  # negative eigenvalue therefore signals floating-point round-off or a degenerate
  # / corrupted per-fit Gamma rather than imputation noise. Abort rather than
  # silently project to the nearest PD matrix (which would distort the sandwich
  # meat and the test statistic without warning).
  ev <- eigen(Gamma_tilde, symmetric = TRUE, only.values = TRUE)$values
  tol <- sqrt(.Machine$double.eps) * max(abs(ev))
  if (min(ev) < -tol) {
    cli::cli_abort(
      c("The pooled asymptotic covariance is not positive semidefinite.",
        "x" = "Its smallest eigenvalue is {.val {min(ev)}}, beyond round-off; a per-imputation asymptotic covariance is likely degenerate or corrupted.",
        "i" = "Check the per-imputation {.code $Gamma} matrices (e.g. a near-singular polychoric correlation); re-fitting the affected imputations usually resolves it."),
      class = "efa_pooled_mi2s_acov_not_psd"
    )
  }

  # Preserve the off-diagonal labelling/order from the per-fit Gamma so the
  # sandwich core and any downstream consumer see the same schema.
  dimnames(Gamma_tilde) <- dimnames(G_list[[1]])
  dimnames(r_bar) <- dimnames(R_list[[1]])

  list(r_bar = r_bar, Gamma_tilde = Gamma_tilde, m = m)
}

.efa_pooled_mi2s <- function(fits, data_list, efa_args, settings, method,
                             rotation, rotation_type, target_method,
                             align_unrotated, fit_pool_method, p,
                             rmsea_ci_level, rmsr_upper) {
  # Two-stage (MI2S) pooled-inputs route for component fits carrying sandwich
  # (robust) SEs. Pools the correlation matrix and its asymptotic covariance
  # across imputations (.efa_pooled_mi2s_inputs), then fits the model once on the
  # pooled inputs via the shared .efa_core(), so the resulting object carries
  # native scaled-shifted chi-square and sandwich SEs that already reflect the
  # multiple-imputation uncertainty. No per-imputation rotation alignment or
  # Rubin pooling of estimates is performed on this path.
  m_imp <- length(fits)

  ## ---- Fail closed on inconsistent inputs --------------------------------
  ses <- vapply(fits, .efa_pooled_setting_chr, character(1), name = "se")
  cor_methods <- vapply(fits, .efa_pooled_setting_chr, character(1),
                        name = "cor_method")
  gamma_ok <- vapply(fits, function(f) {
    !is.null(f$Gamma) && !anyNA(f$Gamma)
  }, logical(1))
  if (anyNA(ses) || !all(ses == "sandwich") ||
      length(unique(cor_methods)) != 1L ||
      !cor_methods[[1L]] %in% c("poly", "tetra", "pearson") ||
      !all(gamma_ok)) {
    cli::cli_abort(
      c("MI2S pooling requires every imputation to be fitted with {.code se = \"sandwich\"} and the same {.arg cor_method} ({.val poly}, {.val tetra}, or {.val pearson}), each carrying a valid asymptotic covariance.",
        "i" = "Re-fit every imputation with the same {.code se = \"sandwich\"} setting and {.arg cor_method}."),
      class = "efa_pooled_mi2s_inputs_inconsistent"
    )
  }
  cor_method <- cor_methods[[1L]]
  n_factors <- settings$n_factors
  use_setting <- settings$use
  if (is.null(use_setting)) use_setting <- "pairwise.complete.obs"
  type_setting <- settings$type
  if (is.null(type_setting)) type_setting <- "EFAtools"

  ## ---- Sample size -------------------------------------------------------
  Ns <- .efa_pooled_get_Ns(data_list, fits, efa_args)
  Ns_ok <- Ns[is.finite(Ns)]
  if (length(Ns_ok) == 0L) {
    cli::cli_abort(
      c("MI2S sandwich pooling requires the sample size {.arg N}.",
        "i" = "Supply raw data or {.arg N} for the imputations."),
      class = "efa_pooled_mi2s_no_n"
    )
  }
  if (length(unique(Ns_ok)) > 1L) {
    cli::cli_warn(
      c("The imputed datasets appear to have different N.",
        "i" = "MI2S pooling uses the mean N across imputations."),
      class = "efa_pooled_unequal_n"
    )
  }
  N_pool <- mean(Ns_ok)

  ## ---- Imputation-count guidance -----------------------------------------
  if (m_imp < 20L) {
    cli::cli_warn(
      c("MI2S pooling was run with only {m_imp} imputation{?s}.",
        "i" = "The scaled-shifted statistic is calibrated for 20 or more imputations (more at higher rates of missingness; Sriutaisuk et al. 2025); interpret the pooled fit and SEs with caution."),
      class = "efa_pooled_mi2s_n_too_small"
    )
  }

  ## ---- Alignment settings do not apply to a single pooled fit ------------
  if (!identical(target_method, "first_target") ||
      !identical(align_unrotated, "signed_tucker_congruence")) {
    cli::cli_warn(
      c("{.arg target_method} and {.arg align_unrotated} are ignored on the MI2S sandwich path.",
        "i" = "MI2S pools the inputs and fits a single solution, so no per-imputation rotation alignment is performed."),
      class = "efa_pooled_mi2s_alignment_ignored"
    )
  }

  ## ---- Pool inputs and fit once ------------------------------------------
  pooled <- .efa_pooled_mi2s_inputs(fits)
  r_bar <- pooled$r_bar
  Gamma_tilde <- pooled$Gamma_tilde

  weights <- if (method == "DWLS") {
    # A non-positive pooled asymptotic variance makes .poly_weight_matrix() abort
    # with the low-level efa_dwls_degenerate_weight; relabel it to the documented
    # MI2S condition with the "increase imputations" remediation.
    tryCatch(
      .poly_weight_matrix(diag(Gamma_tilde), ncol(r_bar)),
      efa_dwls_degenerate_weight = function(e) {
        cli::cli_abort(
          c("The pooled asymptotic covariance has a non-positive diagonal entry, so DWLS weights cannot be formed.",
            "i" = "Increase the number of imputations (Sriutaisuk et al. 2025 recommend 20 or more, and more at higher rates of missingness)."),
          class = "efa_pooled_mi2s_acov_not_psd", parent = e
        )
      }
    )
  } else {
    NULL
  }

  # Reuse the same estimate -> rotate -> SE pipeline EFA() runs, defaulting the
  # estimation/rotation tuning to whatever the component fits used (the extra
  # arguments forwarded through `...`); .efa_core()'s argument defaults supply the
  # rest. The arguments resolved explicitly here are dropped to avoid passing them
  # twice through do.call().
  drop_args <- c("x", "data", "n_factors", "N", "method", "rotation", "se",
                 "type", "use", "cor_method", "ci", "seed")
  extra_args <- efa_args[setdiff(names(efa_args), drop_args)]
  mi_fit <- do.call(.efa_core, c(
    list(R = r_bar, N = N_pool, weights = weights, Gamma = Gamma_tilde,
         np_boot = FALSE, method = method, rotation = rotation,
         type = type_setting, n_factors = n_factors, se = "sandwich",
         ci = 1 - p, use = use_setting, cor_method = cor_method),
    extra_args
  ))

  # .efa_core()/EFA() compute RMSEA confidence bounds at a fixed 90% level; honor
  # the requested rmsea_ci_level on the pooled object by recomputing the bounds
  # from the single fit's pooled test statistic, matching how the information and
  # bootstrap routes apply rmsea_ci_level (a no-op when rmsea_ci_level == 0.90).
  fi <- mi_fit$fit_indices
  if (!is.null(fi) && is.finite(fi$chi) && is.finite(fi$df) && fi$df > 0) {
    rmsea_ci <- .efa_pooled_rmsea_ci(fi$chi, fi$df, N_pool, level = rmsea_ci_level)
    fi$RMSEA_LB <- unname(rmsea_ci[["lower"]])
    fi$RMSEA_UB <- unname(rmsea_ci[["upper"]])
    mi_fit$fit_indices <- fi
  }

  ## ---- Assemble the pooled object ----------------------------------------
  settings_pooled <- settings
  settings_pooled$N <- N_pool
  settings_pooled$pooled_N <- N_pool
  settings_pooled$pooled <- TRUE
  settings_pooled$component_se <- "sandwich"
  settings_pooled$se <- "sandwich"
  settings_pooled$n_imputations <- m_imp
  settings_pooled$target_method <- target_method
  settings_pooled$align_unrotated <- align_unrotated
  settings_pooled$fit_pool_method <- fit_pool_method
  settings_pooled$p <- p
  settings_pooled$ci <- 1 - p
  settings_pooled$rmsea_ci_level <- rmsea_ci_level
  settings_pooled$rmsr_upper <- rmsr_upper

  results <- list(
    h2 = mi_fit$h2,
    unrot_loadings = mi_fit$unrot_loadings,
    vars_accounted = mi_fit$vars_accounted,
    fit_indices = mi_fit$fit_indices,
    model_implied_R = mi_fit$model_implied_R,
    residuals = mi_fit$residuals,
    orig_R = mi_fit$orig_R,
    settings = settings_pooled,
    fits = fits,
    alignment = NULL,
    mi_diagnostics = NULL,
    mi_fit = mi_fit
  )

  if (rotation_type != "none") {
    results$rot_loadings <- mi_fit$rot_loadings
    results$vars_accounted_rot <- mi_fit$vars_accounted_rot
  }
  if (rotation_type == "oblique") {
    results$Phi <- mi_fit$Phi
    results$Structure <- mi_fit$Structure
  }

  results$SE <- mi_fit$SE
  results$CI <- mi_fit$CI
  # No per-parameter Rubin pooling on the MI2S path: the imputation uncertainty
  # is carried by Gamma_tilde, not by element-wise between-imputation variance.
  # The slots stay present-but-NULL so the schema matches the other pooled paths.
  results["MI"] <- list(NULL)
  results["replicates"] <- list(NULL)

  class(results) <- c("EFA_POOLED", "EFA")
  results
}

.efa_pooled_get_Ns <- function(data_list, fits, efa_args) {
  # Recover the N used in each EFA fit. Correlation-matrix input may not carry N,
  # so return NA unless N was supplied to EFA() or stored in settings.
  vapply(seq_along(data_list), function(d) {
    if (!is.null(fits[[d]]$settings$N) && !is.na(fits[[d]]$settings$N)) {
      return(as.numeric(fits[[d]]$settings$N))
    }
    if (!is.null(efa_args$N) && !is.na(efa_args$N)) {
      return(as.numeric(efa_args$N))
    }
    if (!.is_cormat(data_list[[d]])) {
      return(nrow(data_list[[d]]))
    }
    NA_real_
  }, numeric(1))
}

.efa_pooled_check_fits <- function(fits) {
  # Fail early if the fitted EFA objects are not conformable. Pooling only makes
  # sense when all imputations estimate the same model on the same variables.
  if (length(fits) < 2L) {
    cli::cli_abort("At least two EFA fits are required for MI pooling.", class = "efa_pooled_min_fits")
  }

  dims <- vapply(fits, function(x) {
    paste(dim(as.matrix(x$unrot_loadings)), collapse = "x")
  }, character(1))
  if (length(unique(dims)) != 1L) {
    cli::cli_abort("All unrotated loading matrices must have the same dimensions.", class = "efa_pooled_dim_mismatch")
  }

  var_names <- lapply(fits, function(x) rownames(as.matrix(x$orig_R)))
  if (!all(vapply(var_names[-1], identical, logical(1), var_names[[1]]))) {
    cli::cli_abort("All imputations must contain the same variables in the same order.", class = "efa_pooled_var_mismatch")
  }

  for (nm in c("method", "rotation", "n_factors")) {
    vals <- vapply(fits, .efa_pooled_setting_chr, character(1), name = nm)
    vals <- vals[!is.na(vals)]
    if (length(unique(vals)) > 1L) {
      cli::cli_abort("All imputations must use the same {.arg {nm}}.", class = "efa_pooled_setting_mismatch")
    }
  }

  invisible(TRUE)
}

.efa_pooled_align_unrotated_list <- function(unrot_loadings,
                                             align_unrotated = c("signed_tucker_congruence", "none", "procrustes"),
                                             return_meta = FALSE) {
  # Unrotated factor axes are arbitrary up to ordering and signs. This helper
  # puts them into a common orientation before simple arithmetic averaging.
  # Unlike the rotated solution below, this step should not seek simple
  # structure; it only removes indeterminacy in the unrotated axes.
  #
  # When return_meta = TRUE, a parallel list of per-imputation alignment metadata
  # is returned alongside the aligned loadings so that downstream consumers can
  # apply the same column permutation (and, when signed gauges matter, sign
  # vector) to other per-imputation matrices that share the loading gauge, e.g.
  # marginal-SE matrices used by Rubin's-rules pooling of analytic SEs.
  align_unrotated <- match.arg(align_unrotated)

  D <- length(unrot_loadings)
  if (align_unrotated == "none") {
    if (!return_meta) return(unrot_loadings)
    k_first <- if (D > 0L) ncol(as.matrix(unrot_loadings[[1]])) else 0L
    identity_meta <- list(type         = "none",
                          factor_order = seq_len(k_first),
                          factor_sign  = rep(1, k_first))
    meta <- rep(list(identity_meta), D)
    return(list(loadings = unrot_loadings, meta = meta))
  }

  target <- as.matrix(unrot_loadings[[1]])
  k <- ncol(target)

  out  <- vector("list", D)
  meta <- vector("list", D)
  out[[1]] <- target
  meta[[1]] <- if (align_unrotated == "signed_tucker_congruence") {
    list(type         = "signed_tucker_congruence",
         factor_order = seq_len(k),
         factor_sign  = rep(1, k))
  } else {
    # Procrustes anchor gauge: Q_1 = I_k, so the aligned analytic-SE propagation
    # reduces to the identity on the first imputation.
    list(type = "procrustes", Q = diag(k))
  }

  for (d in seq_along(unrot_loadings)[-1]) {
    Ld <- as.matrix(unrot_loadings[[d]])

    if (align_unrotated == "signed_tucker_congruence") {
      # Preserve the axes up to a one-to-one permutation/sign change based on
      # Tucker congruence. .align_solution() implements the assignment and sign
      # correction without applying a continuous rotation.
      al <- .align_solution(L_target = target, L = Ld)
      out[[d]]  <- al$loadings
      meta[[d]] <- list(type         = "signed_tucker_congruence",
                        factor_order = as.integer(al$factor_order),
                        factor_sign  = as.numeric(al$factor_sign))
    } else if (align_unrotated == "procrustes") {
      # Continuous orthogonal alignment: Q_d minimises ||L_d Q_d - L_1||_F over
      # orthogonal Q. The transform is retained alongside the aligned loadings
      # so the analytic-SE pool can propagate the per-imputation full unrotated
      # covariance via the Kronecker identity Var(vec(L_d Q_d)) = (Q_d' (x) I_p)
      # V_d (Q_d (x) I_p).
      pr <- PROCRUSTES(A = Ld, Target = target, rotation = "orthogonal")
      out[[d]]  <- .change_class(pr$loadings, "matrix")
      meta[[d]] <- list(type = "procrustes", Q = pr$T)
    }
  }

  if (!return_meta) return(out)
  list(loadings = out, meta = meta)
}

.efa_pooled_rmsea_ci <- function(chi, df, N, level = .90) {
  # RMSEA CI from the noncentral chi-square inversion used by EFAtools.
  # The chi-square supplied here should already be the pooled test statistic.
  if (is.na(chi) || is.na(df) || is.na(N) || df <= 0 || N <= 1) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  alpha <- 1 - level
  lower_goal <- 1 - alpha / 2
  upper_goal <- alpha / 2

  pfun <- function(lambda, goal) {
    goal - stats::pchisq(chi, df = df, ncp = lambda)
  }

  lambda_l <- 0
  lambda_u <- 0

  if (stats::pchisq(chi, df = df, ncp = 0) >= lower_goal) {
    lambda_l <- tryCatch(
      stats::uniroot(pfun, interval = c(1e-10, 10000), goal = lower_goal,
                     extendInt = "upX", maxiter = 100L)$root,
      error = function(e) NA_real_
    )
  }

  if (stats::pchisq(chi, df = df, ncp = 0) >= upper_goal) {
    lambda_u <- tryCatch(
      stats::uniroot(pfun, interval = c(1e-10, 10000), goal = upper_goal,
                     extendInt = "upX", maxiter = 100L)$root,
      error = function(e) NA_real_
    )
  }

  denom <- df * (N - 1)
  if (!is.finite(denom) || denom <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  # Cap the bounds at 1 like the non-pooled .chi_fit_indices() RMSEA CI, so all
  # pooled routes (information, bootstrap, MI2S) report RMSEA in [0, 1].
  c(lower = min(sqrt(lambda_l / denom), 1),
    upper = min(sqrt(lambda_u / denom), 1))
}

.efa_pooled_D2 <- function(chis, df) {
  # D2 pools complete-data chi-square statistics across imputations. The
  # returned chi value is the asymptotic chi-square approximation df * F_D2,
  # which is used downstream for RMSEA/CFI.
  chis <- chis[is.finite(chis)]
  M <- length(chis)

  if (M < 2L || !is.finite(df) || df <= 0) {
    return(NULL)
  }

  Tbar <- mean(chis)

  if (Tbar <= 0) {
    return(list(
      F = 0,
      chi = 0,
      df1 = df,
      df2 = Inf,
      p = 1,
      ARIV = 0,
      FMI = 0,
      Tbar = Tbar,
      M = M
    ))
  }

  # Average relative increase in variance (Li, Meng, Raghunathan & Rubin, 1991):
  # the between-imputation variance of the sqrt-transformed statistics is centred
  # on mean(sqrt(chi^2)) = stats::var(sqrt(chis)), not on sqrt(mean(chi^2)).
  ARIV <- (1 + 1 / M) * stats::var(sqrt(chis))

  if (ARIV <= .Machine$double.eps) {
    F_D2 <- Tbar / df
    df2 <- Inf
    p_val <- stats::pchisq(Tbar, df = df, lower.tail = FALSE)
  } else {
    F_D2 <- (Tbar / df - ((M + 1) / (M - 1)) * ARIV) / (1 + ARIV)
    df2 <- df^(-3 / M) * (M - 1) * (1 + 1 / ARIV)^2
    p_val <- if (F_D2 < 0) 1 else stats::pf(F_D2, df1 = df, df2 = df2, lower.tail = FALSE)
  }

  chi_D2 <- max(0, df * F_D2)
  FMI <- ARIV / (1 + ARIV)

  list(
    F = F_D2,
    chi = chi_D2,
    df1 = df,
    df2 = df2,
    p = p_val,
    ARIV = ARIV,
    FMI = FMI,
    Tbar = Tbar,
    M = M
  )
}

.efa_pooled_fit_indices <- function(fits,
                                    pooled_R,
                                    residuals,
                                    RMSR,
                                    N,
                                    Ns = NULL,
                                    method,
                                    pool_method = "D2",
                                    rmsea_ci_level = .90) {

  fit_list <- .extract_list_object(fits, "fit_indices")
  keep <- !vapply(fit_list, is.null, logical(1))
  fit_list <- fit_list[keep]

  # Per-imputation N aligned to the kept fits. Each imputation's chi-square was
  # computed in .gof() with that imputation's own N, so the rescaling to the
  # common (N - 1) scale used by the CFI/TLI diagnostics below must use the same
  # per-imputation N.
  Ns_kept <- if (is.null(Ns)) rep(NA_real_, length(fit_list)) else as.numeric(Ns)[keep]

  p_vars <- nrow(pooled_R)
  df <- NA_real_
  chis <- numeric(0)
  chis_null <- numeric(0)

  if (length(fit_list) > 0L) {
    dfs <- vapply(fit_list, function(x) {
      if (!is.null(x$df)) x$df else NA_real_
    }, numeric(1))
    finite_dfs <- unique(stats::na.omit(dfs))
    if (length(finite_dfs) > 1L) {
      cli::cli_abort("Cannot D2-pool chi-square fit because the imputation-specific dfs differ.",
                     class = "efa_pooled_chisq_df")
    }
    if (length(finite_dfs) == 1L) {
      df <- finite_dfs[[1L]]
    }

    chis <- vapply(fit_list, function(x) {
      if (!is.null(x$chi)) x$chi else NA_real_
    }, numeric(1))
    chis_null <- vapply(fit_list, function(x) {
      if (!is.null(x$chi_null)) x$chi_null else NA_real_
    }, numeric(1))
  }

  # D2-pool a vector of per-imputation chi-squares at `dfree` degrees of freedom,
  # or NULL when pooling does not apply (non-D2 request, no statistics, or a
  # degenerate df). Shared by the reported model/baseline statistics and their
  # (N - 1)-scaled diagnostic counterparts (chi_cfi / chi_null_cfi).
  .pool_chi <- function(chi_stats, dfree) {
    if (identical(pool_method, "D2") && length(chi_stats) > 0L &&
        is.finite(dfree) && dfree > 0) {
      .efa_pooled_D2(chi_stats, dfree)
    } else {
      NULL
    }
  }

  D2 <- .pool_chi(chis, df)
  if (!is.null(D2)) {
    chi <- D2$chi
    p_chi <- D2$p
  } else {
    chi <- NA_real_
    p_chi <- NA_real_
  }

  df_null <- p_vars * (p_vars - 1) / 2
  D2_null <- .pool_chi(chis_null, df_null)

  if (!is.null(D2_null)) {
    chi_null <- D2_null$chi
    p_null <- D2_null$p
  } else if (is.finite(N) && N > 1) {
    chi_null <- .null_chisq(pooled_R, N)
    p_null <- stats::pchisq(chi_null, df_null, lower.tail = FALSE)
  } else {
    chi_null <- NA_real_
    p_null <- NA_real_
  }

  # The incremental indices CFI and TLI are reported as the average of the
  # per-imputation indices rather than formed from separately pooled model/baseline
  # chi-squares; see the rationale in the function documentation (@details). Average
  # via .efa_pooled_col_means(), the same finite-mean primitive that centres the
  # bootstrap fit_indices_descriptive CI, so the reported point estimate and that CI
  # share one definition.
  cfis <- vapply(fit_list, function(x) {
    if (!is.null(x$CFI)) x$CFI else NA_real_
  }, numeric(1))
  tlis <- vapply(fit_list, function(x) {
    if (!is.null(x$TLI)) x$TLI else NA_real_
  }, numeric(1))
  # Columns are positional (col_means drops dimnames): column 1 is CFI, column 2 TLI.
  inc_means <- .efa_pooled_col_means(cbind(cfis, tlis))
  CFI <- inc_means[[1L]]
  TLI <- inc_means[[2L]]

  # The model and baseline chi-squares D2-pooled on the common (N - 1) noncentrality
  # scale, the basis on which lavaan.mi/semTools form pooled incremental indices.
  # chi_cfi is the (N - 1)-scaled pooled model statistic that drives the reported RMSEA
  # below; chi_null_cfi is exposed via mi_diagnostics so the pooled fit can be reconciled
  # against that reference (the reported CFI/TLI are averaged per imputation, not formed
  # from these). Each imputation's Bartlett-corrected chi is mapped back to the
  # (N - 1) scale -- chi = F * mult implies F * (N - 1) = chi * (N - 1) / mult, with
  # .bartlett_mult() shared with .gof()/.null_chisq() so the multiplier cannot drift
  # -- before pooling, since D2 is non-linear in a constant rescaling. The model and
  # baseline use different multipliers (the factor-count term (2q)/3 enters only the
  # model), so a small-N imputation can pass one finiteness mask (ok_m / ok_n) and
  # not the other; the model mask ok_m gates the RMSEA statistic, ok_n the chi_null_cfi
  # diagnostic. q is read from the first fit, valid because .efa_pooled_check_fits()
  # pins an equal factor count across imputations.
  chi_cfi <- NA_real_
  chi_null_cfi <- NA_real_
  if (identical(pool_method, "D2") && length(chis) > 0L && is.finite(df) && df > 0) {
    q <- ncol(as.matrix(fits[[1L]]$unrot_loadings))
    mult_model <- .bartlett_mult(Ns_kept, p_vars, q)
    mult_null <- .bartlett_mult(Ns_kept, p_vars)
    chis_cfi <- rep(NA_real_, length(chis))
    ok_m <- is.finite(Ns_kept) & is.finite(mult_model) & mult_model > 0
    chis_cfi[ok_m] <- chis[ok_m] * (Ns_kept[ok_m] - 1) / mult_model[ok_m]
    chis_null_cfi <- rep(NA_real_, length(chis_null))
    ok_n <- is.finite(Ns_kept) & is.finite(mult_null) & mult_null > 0
    chis_null_cfi[ok_n] <- chis_null[ok_n] * (Ns_kept[ok_n] - 1) / mult_null[ok_n]

    D2_cfi <- .pool_chi(chis_cfi, df)
    if (!is.null(D2_cfi)) chi_cfi <- D2_cfi$chi
    D2_null_cfi <- .pool_chi(chis_null_cfi, df_null)
    if (!is.null(D2_null_cfi)) chi_null_cfi <- D2_null_cfi$chi
  }

  if (is.finite(chi) && is.finite(df) && is.finite(N) && N > 1) {
    n_params <- p_vars * (p_vars + 1) / 2 - df
    ECVI <- (chi + 2 * n_params) / (N - 1)
  } else {
    ECVI <- NA_real_
  }

  # RMSEA is built on the uncorrected (N - 1) discrepancy scale (as in
  # .chi_fit_indices()), so it uses the (N - 1)-scaled pooled model statistic chi_cfi,
  # not the Bartlett-corrected chi (which drives the reported model chi-square test).
  rmsea_denom <- df * (N - 1)
  if (is.finite(chi_cfi) && is.finite(df) && is.finite(N) &&
      df > 0 && N > 1 && is.finite(rmsea_denom) && rmsea_denom > 0) {
    RMSEA <- .rmsea_point(chi_cfi, df, N)
    rmsea_ci <- .efa_pooled_rmsea_ci(chi_cfi, df, N, level = rmsea_ci_level)
  } else {
    RMSEA <- NA_real_
    rmsea_ci <- c(lower = NA_real_, upper = NA_real_)
  }

  # These mirror the chi-square-derived quantities used elsewhere in EFAtools.
  # Under MI/D2 pooling they are descriptive only, not likelihood-based MI AIC/BIC.
  AIC <- if (is.finite(chi) && is.finite(df)) chi - 2 * df else NA_real_
  BIC <- if (is.finite(chi) && is.finite(df) && is.finite(N)) chi - log(N) * df else NA_real_

  ## CAF in EFAtools is 1 - KMO(delta_hat), with diagonal set to 1.
  # Ensure symmetry of residuals by (residuals + t(residuals)) / 2. Asymmetry
  # can arise due to floating point imprecision from averaging the matrices.
  delta_hat <- (residuals + t(residuals)) / 2
  diag(delta_hat) <- 1
  CAF <- .compute_caf(delta_hat)

  ## SRMR (standardized RMR; Bentler, 1995) over the pooled residuals: the
  ## model-implied diagonal is 1, so only the off-diagonal residuals contribute,
  ## divided by the count of non-redundant elements p(p + 1)/2. Distinct from RMSR
  ## (off-diagonal mean), which uses the p(p - 1)/2 denominator.
  SRMR <- sqrt(sum(residuals[upper.tri(residuals)]^2) /
                 (nrow(residuals) * (nrow(residuals) + 1) / 2))

  out <- list(
    chi = chi,
    df = df,
    p_chi = p_chi,
    CAF = CAF,
    CFI = CFI,
    TLI = TLI,
    RMSEA = RMSEA,
    RMSEA_LB = rmsea_ci[["lower"]],
    RMSEA_UB = rmsea_ci[["upper"]],
    AIC = AIC,
    BIC = BIC,
    ECVI = ECVI,
    RMSR = RMSR,
    SRMR = SRMR,
    chi_null = chi_null,
    df_null = df_null,
    p_null = p_null,
    pool_method = pool_method,
    mi_diagnostics = list(
      D2_F = if (!is.null(D2)) D2$F else NA_real_,
      D2_df1 = if (!is.null(D2)) D2$df1 else NA_real_,
      D2_df2 = if (!is.null(D2)) D2$df2 else NA_real_,
      D2_chi_asymptotic = if (!is.null(D2)) D2$chi else NA_real_,
      chi_bar_naive = if (!is.null(D2)) D2$Tbar else if (length(chis) > 0L) mean(chis, na.rm = TRUE) else NA_real_,
      D2_null_F = if (!is.null(D2_null)) D2_null$F else NA_real_,
      D2_null_df1 = if (!is.null(D2_null)) D2_null$df1 else NA_real_,
      D2_null_df2 = if (!is.null(D2_null)) D2_null$df2 else NA_real_,
      D2_null_chi_asymptotic = if (!is.null(D2_null)) D2_null$chi else NA_real_,
      chi_null_bar_naive = if (!is.null(D2_null)) D2_null$Tbar else if (length(chis_null) > 0L) mean(chis_null, na.rm = TRUE) else NA_real_,
      ARIV = if (!is.null(D2)) D2$ARIV else NA_real_,
      FMI = if (!is.null(D2)) D2$FMI else NA_real_,
      ARIV_null = if (!is.null(D2_null)) D2_null$ARIV else NA_real_,
      FMI_null = if (!is.null(D2_null)) D2_null$FMI else NA_real_,
      # Pooled model and baseline chi-squares on the common (N - 1) noncentrality
      # scale -- the lavaan.mi/semTools basis for pooled incremental indices, kept
      # for reconciliation against that reference. The reported CFI/TLI average the
      # per-imputation indices instead; these are diagnostic and distinct from the
      # reported chi / chi_null, which keep their Bartlett multipliers.
      chi_cfi = chi_cfi,
      chi_null_cfi = chi_null_cfi,
      m = if (!is.null(D2)) D2$M else length(fits)
    )
  )

  out
}

## -----------------------------------------------------------------------------
## Bootstrap pooling helpers
## -----------------------------------------------------------------------------

.efa_pooled_has_replicates <- function(fits) {
  # Bootstrap MI pooling requires the raw unrotated bootstrap loading replicates.
  # Scalar bootstrap SEs/CIs are not enough because the replicates must be
  # re-expressed in the final MI target space.
  vapply(fits, function(x) {
    !is.null(x$replicates) &&
      !is.null(x$replicates$unrot_loadings)
  }, logical(1))
}

.efa_pooled_vec <- function(M) {
  # Column-major vectorization so matrix estimates can be Rubin-pooled.
  as.vector(as.matrix(M))
}

.efa_pooled_vech <- function(M) {
  # Vectorize the lower triangle of a symmetric Phi matrix without duplicates.
  M <- as.matrix(M)
  M[lower.tri(M, diag = TRUE)]
}

.efa_pooled_unvech <- function(v, k) {
  # Reconstruct a symmetric matrix from vech() output.
  M <- matrix(NA_real_, k, k)
  M[lower.tri(M, diag = TRUE)] <- v
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  M
}

.efa_pooled_communalities <- function(Lambda, Phi = NULL) {
  # Communalities are the diagonal of the common-factor reproduced matrix; read it
  # directly instead of forming the full p x p product:
  # diag(L L') = rowSums(L^2) and diag(L Phi L') = rowSums((L Phi) * L).
  Lambda <- as.matrix(Lambda)
  if (is.null(Phi)) {
    unname(rowSums(Lambda^2))
  } else {
    unname(rowSums((Lambda %*% as.matrix(Phi)) * Lambda))
  }
}

.efa_pooled_model_implied <- function(Lambda, Phi = NULL) {
  # Correlation-model implied matrix: common-factor part plus uniquenesses,
  # implemented by setting the diagonal of Lambda Phi Lambda' to one.
  Lambda <- as.matrix(Lambda)
  common_R <- if (is.null(Phi)) {
    Lambda %*% t(Lambda)
  } else {
    Lambda %*% as.matrix(Phi) %*% t(Lambda)
  }
  implied <- common_R
  diag(implied) <- 1
  implied
}

.efa_pooled_residual_from_solution <- function(R, Lambda, Phi = NULL) {
  # Residuals for one imputation in the same target space as its aligned solution.
  E <- as.matrix(R) - .efa_pooled_model_implied(Lambda, Phi)
  diag(E) <- 0
  E
}

.efa_pooled_make_ci <- function(est, se, df, alpha) {
  # Wald-type confidence intervals using Rubin degrees of freedom when finite.
  # `est` may be a plug-in estimate (communalities, Structure, symmetrised Phi)
  # supplied via the callers' `est_override`, in which case the interval is
  # centred on the plug-in while its half-width (se, df) comes from the pooled
  # column-mean's between/within variance; the two coincide when the plug-in
  # equals that column mean and differ only by the (typically small) plug-in gap.
  crit <- ifelse(
    is.finite(df),
    stats::qt(1 - alpha / 2, df = df),
    stats::qnorm(1 - alpha / 2)
  )
  list(
    lower = est - crit * se,
    upper = est + crit * se
  )
}


.efa_pooled_col_means <- function(x) {
  # Column means that return NA for columns without finite values.
  x <- as.matrix(x)
  if (ncol(x) < 1L) {
    return(numeric(0))
  }
  vapply(seq_len(ncol(x)), function(j) {
    vals <- x[, j]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 1L) {
      NA_real_
    } else {
      mean(vals)
    }
  }, numeric(1))
}

.efa_pooled_col_vars <- function(x) {
  # Column variances that return NA for columns with fewer than two finite values.
  # Rubin pooling only needs diagonal within- and between-imputation variances,
  # so computing full covariance matrices is unnecessary and can be prohibitively
  # memory-intensive for residual matrices.
  x <- as.matrix(x)
  if (ncol(x) < 1L) {
    return(numeric(0))
  }
  if (nrow(x) < 2L) {
    return(rep(NA_real_, ncol(x)))
  }
  vapply(seq_len(ncol(x)), function(j) {
    vals <- x[, j]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2L) {
      NA_real_
    } else {
      stats::var(vals)
    }
  }, numeric(1))
}

.efa_pooled_rubin_core <- function(Qbar, Ubar, B, m, N_pool = Inf,
                                   alpha = 0.05) {
  # Element-wise Rubin-style pooling core, shared by the bootstrap (analytic SEs
  # absent) and analytic (per-imputation marginal SEs available) paths.
  #
  # Qbar  pooled point estimate per parameter (vector or matrix-as-vector).
  # Ubar  mean within-imputation variance (vector).
  # B     between-imputation variance (sample, m-1 denominator) (vector).
  # m     number of imputations.
  # N_pool complete-data sample size for Barnard-Rubin small-sample df. When
  #        non-finite, the plain Rubin (1987) df ν_old = (m-1)(1 + 1/r)^2 is
  #        used (lavaan.mi convention). When finite, the Barnard & Rubin
  #        (1999, Biometrika 86) adjustment is applied:
  #           ν_obs = ((ν_com + 1)/(ν_com + 3)) ν_com (1 - λ),  ν_com = N_pool - 1
  #           ν_BR  = 1 / (1/ν_old + 1/ν_obs)
  # alpha 1 - confidence level for the Wald interval.
  #
  # Returns a list(se, ci, RIV, FMI, df). Bootstrap callers receive the same
  # df as before (N_pool defaults to Inf for byte-identical behavior).
  Tot <- Ubar + (1 + 1 / m) * B
  Tot[is.finite(Tot) & Tot < 0 & abs(Tot) < sqrt(.Machine$double.eps)] <- 0
  se <- sqrt(Tot)

  r <- rep(0, length(se))
  has_U <- is.finite(Ubar) & Ubar > 0
  r[has_U] <- ((1 + 1 / m) * B[has_U]) / Ubar[has_U]
  r[!has_U & is.finite(B) & B > 0] <- Inf
  r[!is.finite(B) & !is.infinite(r)] <- NA_real_

  df_old <- rep(Inf, length(se))
  finite_pos_r <- is.finite(r) & r > 0
  df_old[finite_pos_r] <- (m - 1) * (1 + 1 / r[finite_pos_r])^2
  df_old[is.infinite(r)] <- m - 1
  df_old[is.na(r)] <- NA_real_

  # Reported FMI is Rubin's asymptotic lambda = r / (1 + r) (the limiting fraction
  # of missing information; equals lavaan.mi's `fmi`), not mice's finite-m gamma
  # = (r + 2/(df_old + 3)) / (r + 1), which adds a small-m correction. The two
  # agree as m grows; lambda is the conventional EFA/SEM-MI reporting choice.
  fmi <- r / (1 + r)
  fmi[is.infinite(r)] <- 1
  fmi[is.nan(fmi)] <- 0

  if (is.finite(N_pool) && N_pool > 1) {
    nu_com <- N_pool - 1
    nu_obs <- ((nu_com + 1) / (nu_com + 3)) * nu_com * (1 - fmi)
    # ν_BR = 1/(1/ν_old + 1/ν_obs). When ν_old is Inf (no between-variance),
    # 1/ν_old = 0 and ν_BR collapses to ν_obs (the complete-data t reference).
    df <- 1 / (1 / df_old + 1 / nu_obs)
    df[is.infinite(df_old)] <- nu_obs[is.infinite(df_old)]
    # When fmi -> 1 (Ubar = 0, B > 0), ν_obs -> 0 and the BR formula collapses
    # to 0, sending qt(., 0) to NaN. Fall back to the plain-Rubin ν_old in that
    # regime so the Wald interval is well-defined (wide, but finite).
    no_obs <- is.finite(nu_obs) & nu_obs == 0
    df[no_obs] <- df_old[no_obs]
    df[is.na(df_old)] <- NA_real_
  } else {
    df <- df_old
  }

  ci <- .efa_pooled_make_ci(Qbar, se, df, alpha)

  list(
    se  = se,
    Tot = Tot,
    ci  = ci,
    RIV = r,
    FMI = fmi,
    df  = df
  )
}

.efa_pooled_rubin_pool <- function(q_list, boot_mat_list, alpha = 0.05,
                                   est_override = NULL) {
  # Rubin-style pooling for one parameter family.
  #
  # q_list contains one aligned point-estimate vector per imputation. boot_mat_list
  # contains aligned bootstrap replicates for estimating the within-imputation
  # variances U_d. Only diagonal variances are required for the returned SEs, CIs,
  # RIV, FMI, and df, so the implementation intentionally avoids building full
  # p x p covariance matrices for each parameter vector.
  m <- length(q_list)
  q_mat <- do.call(rbind, q_list)
  q_bar <- .efa_pooled_col_means(q_mat)

  U_diag_mat <- do.call(rbind, lapply(boot_mat_list, .efa_pooled_col_vars))
  U_bar <- colMeans(U_diag_mat)
  B_mi <- if (m > 1L) {
    .efa_pooled_col_vars(q_mat)
  } else {
    rep(0, ncol(q_mat))
  }

  est <- if (is.null(est_override)) q_bar else est_override
  core <- .efa_pooled_rubin_core(est, U_bar, B_mi, m, N_pool = Inf, alpha = alpha)

  list(
    est = est,
    se  = core$se,
    ci  = core$ci,
    RIV = core$RIV,
    FMI = core$FMI,
    df  = core$df
  )
}

.efa_pooled_propagate_procrustes_vcov <- function(V, Q, p, k) {
  # Marginal SEs of L_d Q_d implied by the per-imputation full unrotated
  # covariance V_d (p k x p k over column-major vec(L_d)) and a fixed orthogonal
  # Q_d (k x k). Under the column-major identity vec(L_d Q_d) = (Q_d' (x) I_p)
  # vec(L_d), V_d^aligned = (Q_d' (x) I_p) V_d (Q_d (x) I_p). For element (i, j)
  # only the variance of L_d Q_d at row i is required, so the implementation
  # extracts the k x k row sub-block V_d^(i) at indices (0:(k-1)) p + i and
  # returns sqrt(diag(Q_d' V_d^(i) Q_d)) per row. Cost is O(p k^3).
  #
  # Q_d is treated as fixed (Schoenemann 1966; mirrors the lavaan / mice
  # convention of not differentiating the per-imputation orthogonal Procrustes
  # solution back through the SVD on vec(L_d)).
  #
  # NA propagation is fail-closed at the row level: if any element of the row
  # block of V or of Q is non-finite, the entire row of the returned p x k
  # marginal-SE matrix is NA so the per-element NA mask in the analytic pool
  # can blank out the affected pooled outputs.
  q_finite <- all(is.finite(Q))
  out <- matrix(NA_real_, p, k)
  if (!q_finite) return(out)
  for (i in seq_len(p)) {
    idx <- (seq_len(k) - 1L) * p + i
    block <- V[idx, idx, drop = FALSE]
    if (!all(is.finite(block))) next
    # diag(Q' (block Q)) = colSums(Q * (block Q)); avoids forming the full k x k
    # product only to keep its diagonal.
    BQ <- block %*% Q
    out[i, ] <- sqrt(pmax(colSums(Q * BQ), 0))
  }
  out
}

.efa_pooled_analytic_marginal <- function(q_mat, se_mat, m, alpha,
                                          est_override = NULL) {
  # Rubin pool of one parameter family from aligned per-imputation point estimates
  # (rows of `q_mat`) and aligned per-imputation marginal SEs (rows of `se_mat`, so
  # the within-imputation variance U_d is `se_mat^2`). Only diagonal variances are
  # needed. Fail-closed: any per-imputation NA in q or se for an element blanks the
  # pooled outputs (SE, CI, RIV, FMI, df) at that element. `est_override` supplies
  # the pooled point estimate for plug-in quantities (communalities, structure)
  # whose returned estimate is computed from the pooled solution, not the mean of
  # the per-imputation values.
  Qbar <- .efa_pooled_col_means(q_mat)
  Ubar <- .efa_pooled_col_means(se_mat^2)
  B    <- if (m > 1L) .efa_pooled_col_vars(q_mat) else rep(0, ncol(q_mat))

  na <- apply(is.na(q_mat) | is.na(se_mat), 2, any)
  Qbar[na] <- NA_real_; Ubar[na] <- NA_real_; B[na] <- NA_real_

  est <- if (is.null(est_override)) Qbar else est_override
  # Plain Rubin (1987) df (N_pool = Inf): analytic information-matrix loadings are
  # asymptotically normal (Wald-z), so the complete-data reference is normal
  # (nu_com = Inf) and the Barnard-Rubin small-sample adjustment does not apply -- it
  # reduces to plain Rubin. This matches lavaan.mi, which reports df -> Inf for
  # asymptotically-normal SEM parameters.
  .efa_pooled_rubin_core(est, Ubar, B, m, N_pool = Inf, alpha = alpha)
}

.efa_pooled_assemble_family <- function(pool, reshape, method = NULL) {
  # Reshape one Rubin-pooled family's SE, Wald CI bounds, and MI diagnostics with a
  # single reshape function (loading matrix, uniqueness/communality vector, or
  # vech-to-symmetric). `method` tags the gauge-alignment provenance when supplied.
  mi <- list(RIV = reshape(pool$RIV), FMI = reshape(pool$FMI), df = reshape(pool$df))
  if (!is.null(method)) mi$method <- method
  list(
    SE = reshape(pool$se),
    CI = list(lower = reshape(pool$ci$lower), upper = reshape(pool$ci$upper)),
    MI = mi
  )
}

.efa_pooled_analytic_pool <- function(fits,
                                      unrot_loadings_aligned,
                                      align_meta,
                                      ci = 0.95,
                                      align_unrotated = c("signed_tucker_congruence",
                                                          "none", "procrustes"),
                                      rotation_type = c("none", "orthogonal",
                                                        "oblique"),
                                      rot_loadings = NULL,
                                      phis = NULL,
                                      structure_loadings = NULL,
                                      mean_structure_loadings = NULL,
                                      mean_phis = NULL,
                                      h2 = NULL) {
  # Rubin's-rules pool of analytic per-imputation marginal SEs for the unrotated
  # loadings and uniquenesses. Serves component fits run with se = "information":
  # each fit carries $SE$unrot_loadings (p x k) and $SE$uniquenesses (length p),
  # and there are no bootstrap replicate cubes to estimate within-imputation
  # variance from. (The sandwich route is handled separately by
  # .efa_pooled_mi2s() and never reaches this helper.)
  #
  # Under `align_unrotated = "signed_tucker_congruence"` or `"none"`, the
  # per-imputation gauge transform is a signed permutation that preserves
  # marginal SEs up to a column permutation, so $SE$unrot_loadings is reused
  # directly after applying the alignment's `factor_order`.
  #
  # Under `align_unrotated = "procrustes"`, the per-imputation orthogonal Q_d
  # MIXES loading columns, so marginal SEs alone no longer determine the aligned
  # marginal SE matrix. The full unrotated covariance persisted on each fit
  # (`vcov_unrot_loadings`, populated by `se = "information"` and
  # `se = "sandwich"`) is propagated by
  # `.efa_pooled_propagate_procrustes_vcov()`, which implements
  # diag(Q' V_d^(i) Q) per row i. Q_d is treated as fixed.
  #
  # Per-element NA propagation is fail-closed: if any imputation has NA for an
  # element (in either Q_d or SE_d), every pooled output at that position (SE,
  # CI bounds, RIV, FMI, df) is NA.
  #
  # The returned shape mirrors .efa_pooled_bootstrap_pool() so the EFA_POOLED()
  # assembly block can route either pool through the same downstream code:
  # list(SE, CI, replicates = NULL, MI, n_boot = NULL). Empty inputs or fits
  # missing the analytic SE slot return NULL so the caller falls through to the
  # existing silent-downgrade branch. Fits missing the full vcov slot under the
  # procrustes branch abort with a classed condition (fail-closed; the marginal
  # SE alone cannot represent the aligned gauge).
  align_unrotated <- match.arg(align_unrotated)
  rotation_type <- match.arg(rotation_type)
  if (length(fits) == 0L) return(NULL)
  has_se <- vapply(fits, function(x) {
    !is.null(x$SE) &&
      !is.null(x$SE$unrot_loadings) &&
      !is.null(x$SE$uniquenesses)
  }, logical(1L))
  if (!all(has_se)) {
    if (any(has_se)) {
      cli::cli_warn(
        c("Analytic SE slots were found for only some imputations; pooled analytic SEs/CIs were not computed.",
          "i" = "Inspect each fit's {.code $SE} slot to identify the imputations with missing analytic SEs."),
        class = "efa_pooled_analytic_partial_se"
      )
    }
    return(NULL)
  }

  # Whole-fit reliability gate (every alignment method). `EFA()` NA-fills the
  # entire `vcov_unrot_loadings` when the analytic covariance is unreliable (a
  # Heywood case or a singular bordered information matrix), NA-filling the
  # marginal SE matrices in lockstep and firing `efa_se_unreliable` at fit time.
  # Pooling those would yield an all-NA SE bundle with no diagnostic, so abort
  # here; `EFA_POOLED()` converts this into a classed `efa_pooled_se_unavailable`
  # warning and falls back to point-estimate-only pooling. The gate keys on the
  # whole-matrix vcov NA fill, not on individual `$SE` elements, so a genuine
  # per-element NA (one imputation NA at a single cell) still fail-closes per
  # element via the mask below rather than dropping the entire pool.
  vcov_reliable <- vapply(fits, function(x) {
    is.null(x$vcov_unrot_loadings) || !anyNA(x$vcov_unrot_loadings)
  }, logical(1L))
  if (!all(vcov_reliable)) {
    bad <- which(!vcov_reliable)
    n_bad <- length(bad)
    cli::cli_abort(
      c("Pooled analytic standard errors require a reliable unrotated loading covariance on every imputation.",
        "i" = "{cli::qty(n_bad)}Imputation{?s} {.val {bad}} {cli::qty(n_bad)}carr{?ies/y} an NA-filled {.code $vcov_unrot_loadings} (Heywood case or singular information). Drop or re-fit {cli::qty(n_bad)}{?it/them}."),
      class = "efa_pooled_unreliable_vcov"
    )
  }

  if (align_unrotated == "procrustes") {
    # The Procrustes branch propagates the full covariance through Q_d, so the
    # slot must be present (not merely reliable) on every imputation.
    has_vcov <- vapply(fits, function(x) {
      !is.null(x$vcov_unrot_loadings)
    }, logical(1L))
    if (!all(has_vcov)) {
      bad <- which(!has_vcov)
      n_bad <- length(bad)
      cli::cli_abort(
        c("{.code align_unrotated = \"procrustes\"} with {.code se = \"information\"} requires the full unrotated loading covariance on every imputation.",
          "i" = "{cli::qty(n_bad)}Imputation{?s} {.val {bad}} {cli::qty(n_bad)}{?is/are} missing {.code $vcov_unrot_loadings}; re-fit those imputations or use {.code align_unrotated = \"signed_tucker_congruence\"}."),
        class = "efa_pooled_no_vcov"
      )
    }
  }

  m <- length(fits)
  L_target <- as.matrix(unrot_loadings_aligned[[1]])
  p_vars   <- nrow(L_target)
  k_fac    <- ncol(L_target)
  loading_dn <- dimnames(L_target)
  var_names  <- rownames(L_target)

  # Per-imputation aligned loadings and aligned SEs.
  # - signed_tucker_congruence / none: column-permute the marginal SE matrix by
  #   the alignment's factor_order.
  # - procrustes: propagate the full vcov through Q_d row-by-row to obtain the
  #   aligned marginal variances, then take sqrt to land back on the same
  #   variance-of-vec(L_d Q_d) scale as Ubar_L (which squares SE_mat).
  Q_mat   <- matrix(NA_real_, nrow = m, ncol = p_vars * k_fac)
  SE_mat  <- matrix(NA_real_, nrow = m, ncol = p_vars * k_fac)
  Q_psi   <- matrix(NA_real_, nrow = m, ncol = p_vars)
  SE_psi_mat <- matrix(NA_real_, nrow = m, ncol = p_vars)
  for (d in seq_len(m)) {
    Q_mat[d, ] <- as.vector(as.matrix(unrot_loadings_aligned[[d]]))

    if (identical(align_meta[[d]]$type, "procrustes")) {
      V_d  <- fits[[d]]$vcov_unrot_loadings
      Q_d  <- align_meta[[d]]$Q
      se_d <- .efa_pooled_propagate_procrustes_vcov(V_d, Q_d, p_vars, k_fac)
      SE_mat[d, ] <- as.vector(se_d)
    } else {
      se_d <- as.matrix(fits[[d]]$SE$unrot_loadings)
      fo <- align_meta[[d]]$factor_order
      if (is.null(fo)) {
        # All currently-supported producers of `align_meta` populate
        # `factor_order` (the `none` and `signed_tucker_congruence` branches of
        # `.efa_pooled_align_unrotated_list`); procrustes meta is routed into
        # the V_d / Q_d branch above by the `type` test and never reaches
        # here. Catch any future bypass at the helper boundary rather than
        # silently falling back to an identity permutation that would pool
        # SEs against unaligned loadings.
        cli::cli_abort(
          "Per-imputation alignment metadata is missing for imputation {d}; analytic-SE pooling requires a column-permutation/sign alignment.",
          class = "efa_pooled_analytic_align_meta_missing"
        )
      }
      SE_mat[d, ] <- as.vector(se_d[, fo, drop = FALSE])
    }

    psi_d <- .efa_pooled_communalities(unrot_loadings_aligned[[d]], NULL)
    Q_psi[d, ] <- 1 - as.numeric(psi_d)
    SE_psi_mat[d, ] <- as.numeric(fits[[d]]$SE$uniquenesses)
  }

  # Fail-closed NA mask and Rubin pooling are shared with the rotated families
  # below via .efa_pooled_analytic_marginal().
  alpha <- 1 - ci
  pool_L   <- .efa_pooled_analytic_marginal(Q_mat, SE_mat,     m, alpha)
  pool_psi <- .efa_pooled_analytic_marginal(Q_psi, SE_psi_mat, m, alpha)

  reshape_loading <- function(v) {
    M <- matrix(v, p_vars, k_fac)
    dimnames(M) <- loading_dn
    M
  }
  reshape_psi <- function(v) {
    out <- as.numeric(v)
    if (!is.null(var_names)) names(out) <- var_names
    out
  }

  unrot_asm <- .efa_pooled_assemble_family(pool_L,   reshape_loading)
  psi_asm   <- .efa_pooled_assemble_family(pool_psi, reshape_psi)

  SE <- list(unrot_loadings = unrot_asm$SE, uniquenesses = psi_asm$SE)
  CI <- list(unrot_loadings = unrot_asm$CI, uniquenesses = psi_asm$CI)
  MI <- list(unrot_loadings = unrot_asm$MI, uniquenesses = psi_asm$MI)

  ## Rotated quantities -------------------------------------------------------
  ## Pool the rotated loadings, communalities, and (oblique) factor correlations
  ## and structure coefficients, each in a common MI rotational gauge. A
  ## rotated-loading standard error is conditional on the rotation criterion
  ## (Jennrich 1973, 1974; Browne 2001; Zhang & Preacher 2015), so for orthogonal
  ## and oblique rotations alike the within-imputation variance is each fit's own
  ## criterion-aware delta-method rotated SE (the quantity EFA() returns), reused
  ## after a signed-permutation alignment to the MI target. The between-imputation
  ## spread of the aligned point estimates supplies B, so the pooled SE shares the
  ## estimand of the np-boot pool, which likewise aligns each replicate to the
  ## target. This is a deliberate approximation -- the per-fit SE is conditional on
  ## that fit's rotation optimum rather than the common gauge -- and is flagged in
  ## MI$<param>$method = "signed_permutation_approx". Communalities are
  ## rotation-invariant and pool without any alignment.
  if (rotation_type != "none") {
    # Phi and Structure exist only for a genuine oblique rotation (a single-factor
    # solution carries no factor correlations); gate on k_fac so the oblique branch
    # never dereferences the NULL Phi/Structure slots of a one-factor fit.
    oblique <- rotation_type == "oblique" && k_fac >= 2L

    q_rot  <- matrix(NA_real_, m, p_vars * k_fac)
    se_rot <- matrix(NA_real_, m, p_vars * k_fac)
    q_h2   <- matrix(NA_real_, m, p_vars)
    se_h2  <- matrix(NA_real_, m, p_vars)
    if (oblique) {
      vlen   <- k_fac * (k_fac + 1L) / 2L
      q_phi  <- matrix(NA_real_, m, vlen)
      se_phi <- matrix(NA_real_, m, vlen)
      q_str  <- matrix(NA_real_, m, p_vars * k_fac)
      se_str <- matrix(NA_real_, m, p_vars * k_fac)
    }

    for (d in seq_len(m)) {
      Lr_d  <- as.matrix(rot_loadings[[d]])
      Phi_d <- if (oblique) as.matrix(phis[[d]]) else NULL
      q_rot[d, ] <- .efa_pooled_vec(Lr_d)
      q_h2[d, ]  <- .efa_pooled_communalities(Lr_d, Phi_d)
      se_h2[d, ] <- as.numeric(fits[[d]]$SE$communalities)

      # Within-imputation variance for the rotated loadings = this fit's
      # criterion-aware delta-method rotated SE, reused after a signed-permutation
      # alignment to THIS imputation's pooled rotated loading Lr_d. Aligning to
      # Lr_d (rather than the global target) keeps the SE column order matched to
      # the order q_rot/q_phi/q_str use; aligning to the global target can pick a
      # different order when factors are weakly separated, pairing a cell's pooled
      # estimate with another factor's within-imputation variance. Only the column
      # order is used (SE magnitudes are sign-invariant).
      fo <- .align_solution(
        L_target = Lr_d,
        L        = as.matrix(fits[[d]]$rot_loadings)
      )$factor_order
      se_rot[d, ] <- .efa_pooled_vec(
        as.matrix(fits[[d]]$SE$rot_loadings)[, fo, drop = FALSE]
      )

      if (oblique) {
        # Apply the same column alignment to the per-fit Phi and Structure SEs.
        se_phi[d, ] <- .efa_pooled_vech(
          as.matrix(fits[[d]]$SE$Phi)[fo, fo, drop = FALSE]
        )
        q_phi[d, ]  <- .efa_pooled_vech(Phi_d)
        se_str[d, ] <- .efa_pooled_vec(
          as.matrix(fits[[d]]$SE$Structure)[, fo, drop = FALSE]
        )
        q_str[d, ]  <- .efa_pooled_vec(as.matrix(structure_loadings[[d]]))
      }
    }

    # A component fit whose rotated-SE Jacobian could not be reproduced carries an
    # all-NA $SE$rot_loadings (the upstream efa_se_unreliable signal) even when its
    # unrotated covariance is fine, so the whole-fit reliability gate above does not
    # catch it. Such a fit leaves an all-NA row in se_rot, which the fail-closed mask
    # in .efa_pooled_analytic_marginal() then propagates to every cell of the pooled
    # rotated-loading (and, for oblique fits, Phi/Structure) SE family. Name the
    # affected imputations so the resulting all-NA rotated SE block is diagnosable
    # rather than silent.
    rot_se_unreliable <- which(apply(se_rot, 1L, function(r) all(is.na(r))))
    if (length(rot_se_unreliable) > 0L) {
      n_bad <- length(rot_se_unreliable)
      affected <- if (oblique) {
        "rotated-loading, factor-correlation, and structure-coefficient"
      } else {
        "rotated-loading"
      }
      cli::cli_warn(
        c("Pooled rotated standard errors could not be produced for {n_bad} imputation{?s}.",
          "i" = "Affected: {.val {rot_se_unreliable}} (each carried an NA-filled rotated standard error, the upstream {.code efa_se_unreliable} signal).",
          "i" = "The pooled {affected} standard errors are returned as {.code NA}."),
        class = "efa_pooled_rotated_se_unreliable"
      )
    }

    pool_rot <- .efa_pooled_analytic_marginal(q_rot, se_rot, m, alpha)
    pool_h2  <- .efa_pooled_analytic_marginal(q_h2, se_h2, m, alpha,
                                              est_override = as.numeric(h2))

    # Rotated loadings and structure are p x k like the unrotated loadings, so they
    # share the same factor-column labels (reshape_loading carries loading_dn).
    rot_asm <- .efa_pooled_assemble_family(
      pool_rot, reshape_loading,
      method = "signed_permutation_approx"
    )
    SE$rot_loadings <- rot_asm$SE; CI$rot_loadings <- rot_asm$CI
    MI$rot_loadings <- rot_asm$MI

    h2_asm <- .efa_pooled_assemble_family(pool_h2, reshape_psi,
                                          method = "gauge_invariant")
    SE$h2 <- h2_asm$SE; CI$h2 <- h2_asm$CI; MI$h2 <- h2_asm$MI

    if (oblique) {
      phi_dn <- dimnames(as.matrix(phis[[1]]))
      unvech_phi <- function(v) {
        M <- .efa_pooled_unvech(v, k_fac)
        dimnames(M) <- phi_dn
        M
      }
      # Center the Phi CI on the same symmetrised pooled estimate reported as $Phi
      # (mean_phis, always present on this oblique branch), matching the est_override
      # the h2 and Structure families use; without it the off-diagonal CI would
      # center on the unsymmetrised column mean of vech(Phi), which can differ from
      # the reported (symmetrised) estimate.
      pool_phi <- .efa_pooled_analytic_marginal(
        q_phi, se_phi, m, alpha,
        est_override = .efa_pooled_vech(mean_phis)
      )
      phi_asm <- .efa_pooled_assemble_family(pool_phi, unvech_phi,
                                             method = "signed_permutation_approx")
      # The fixed unit diagonal of Phi is not an estimated parameter: zero its SE
      # and NA its MI diagnostics so the fixed cells stay out of the printed
      # FMI/RIV summaries.
      diag(phi_asm$SE) <- 0
      diag(phi_asm$MI$RIV) <- NA_real_
      diag(phi_asm$MI$FMI) <- NA_real_
      diag(phi_asm$MI$df)  <- NA_real_
      SE$Phi <- phi_asm$SE; CI$Phi <- phi_asm$CI; MI$Phi <- phi_asm$MI

      pool_str <- .efa_pooled_analytic_marginal(
        q_str, se_str, m, alpha,
        est_override = .efa_pooled_vec(mean_structure_loadings)
      )
      str_asm <- .efa_pooled_assemble_family(pool_str, reshape_loading,
                                             method = "signed_permutation_approx")
      SE$Structure <- str_asm$SE; CI$Structure <- str_asm$CI
      MI$Structure <- str_asm$MI
    }
  }

  list(SE = SE, CI = CI, replicates = NULL, MI = MI, n_boot = NULL)
}


.efa_pooled_align_unrot_boot <- function(L_boot, target, align_unrotated) {
  # Align one bootstrap unrotated loading matrix to the pooled unrotated target.
  if (align_unrotated == "none") {
    return(as.matrix(L_boot))
  }

  if (align_unrotated == "procrustes") {
    return(.change_class(PROCRUSTES(
      A = as.matrix(L_boot),
      Target = as.matrix(target),
      rotation = "orthogonal"
    )$loadings, "matrix"))
  }

  .align_solution(
    L_target = as.matrix(target),
    L = as.matrix(L_boot)
  )$loadings
}

.efa_pooled_align_rot_boot <- function(L_boot_unrot, final_target, rotation_type,
                                       procrustes_args = list()) {
  # Re-rotate/re-align the bootstrap unrotated loadings to the final MI target.
  # This is the key step that makes within-imputation bootstrap variances refer
  # to the same estimand as the pooled rotated solution.
  if (rotation_type == "none") {
    return(NULL)
  }

  pr <- do.call(
    PROCRUSTES,
    c(list(A = as.matrix(L_boot_unrot),
           Target = as.matrix(final_target),
           rotation = rotation_type),
      procrustes_args)
  )

  # Drop the replicate only when no valid alignment was produced. A best
  # multi-start fit that did not formally converge is still the lowest-objective
  # alignment available, with well-defined loadings, and is retained.
  if (isFALSE(pr$valid)) {
    return(NULL)
  }

  list(
    Lambda = .change_class(pr$loadings, "matrix"),
    Phi = if (!is.null(pr$Phi)) as.matrix(pr$Phi) else NULL
  )
}

.efa_pooled_rubin_matrix_result <- function(pool, nrow, ncol, dimnames = NULL) {
  # Convert vectorized Rubin output back into matrix-shaped estimates.
  est <- matrix(pool$est, nrow = nrow, ncol = ncol)
  se <- matrix(pool$se, nrow = nrow, ncol = ncol)
  lower <- matrix(pool$ci$lower, nrow = nrow, ncol = ncol)
  upper <- matrix(pool$ci$upper, nrow = nrow, ncol = ncol)

  if (!is.null(dimnames)) {
    dimnames(est) <- dimnames
    dimnames(se) <- dimnames
    dimnames(lower) <- dimnames
    dimnames(upper) <- dimnames
  }

  list(est = est, se = se, ci = list(lower = lower, upper = upper))
}

.efa_pooled_rubin_symmetric_result <- function(pool, k, dimnames = NULL) {
  # Convert vech-based Rubin output back into symmetric matrix form.
  est <- .efa_pooled_unvech(pool$est, k)
  se <- .efa_pooled_unvech(pool$se, k)
  lower <- .efa_pooled_unvech(pool$ci$lower, k)
  upper <- .efa_pooled_unvech(pool$ci$upper, k)

  if (!is.null(dimnames)) {
    dimnames(est) <- dimnames
    dimnames(se) <- dimnames
    dimnames(lower) <- dimnames
    dimnames(upper) <- dimnames
  }

  list(est = est, se = se, ci = list(lower = lower, upper = upper))
}

.efa_pooled_bootstrap_pool <- function(fits,
                                       orig_R_list,
                                       unrot_loadings_aligned,
                                       mean_unrot_loadings,
                                       rot_loadings = NULL,
                                       phis = NULL,
                                       structure_loadings = NULL,
                                       mean_structure_loadings = NULL,
                                       final_target = NULL,
                                       rotation_type = c("none", "orthogonal", "oblique"),
                                       align_unrotated = c("signed_tucker_congruence", "none", "procrustes"),
                                       procrustes_args = list(),
                                       h2,
                                       residuals,
                                       pooled_orig_R,
                                       N,
                                       method,
                                       pool_method = "D2",
                                       rmsea_ci_level = .90,
                                       alpha = 0.05,
                                       rmsr_upper = TRUE) {

  rotation_type <- match.arg(rotation_type)
  align_unrotated <- match.arg(align_unrotated)

  has_boot <- .efa_pooled_has_replicates(fits)
  if (!all(has_boot)) {
    if (any(has_boot)) {
      cli::cli_warn("Bootstrap replicates were found for only some imputations; pooled bootstrap SEs/CIs were not computed.",
                    class = "efa_pooled_partial_boot")
    }
    return(NULL)
  }

  m <- length(fits)
  # Coerce the loop-invariant alignment targets once; mean_unrot_loadings is a
  # LOADINGS object whose as.matrix() coercion is otherwise repeated on every
  # bootstrap replicate inside the d/b loops.
  mean_unrot_target <- as.matrix(mean_unrot_loadings)
  final_target_m <- if (is.null(final_target)) NULL else as.matrix(final_target)
  p_vars <- nrow(mean_unrot_target)
  k <- ncol(mean_unrot_target)
  loading_dimnames <- dimnames(mean_unrot_target)

  # A single factor cannot be rotated and has no factor correlations, so an
  # oblique request on one factor degenerates to the no-Phi case (matching a
  # single-fit EFA()); align such a solution by sign only.
  oblique <- rotation_type == "oblique" && k >= 2L
  boot_proc_rotation <- if (oblique) "oblique" else "orthogonal"

  B_vec <- vapply(fits, function(x) dim(x$replicates$unrot_loadings)[3], integer(1))
  if (length(unique(B_vec)) > 1L) {
    B_use <- min(B_vec)
    cli::cli_warn("The number of bootstrap replicates differs across imputations; using the minimum number available in all imputations.",
                  class = "efa_pooled_unequal_boot")
  } else {
    B_use <- B_vec[[1]]
  }

  if (!is.finite(B_use) || B_use < 2L) {
    cli::cli_warn("At least two bootstrap replicates per imputation are required for pooled bootstrap SEs/CIs.",
                  class = "efa_pooled_min_boot")
    return(NULL)
  }

  q_unrot <- lapply(unrot_loadings_aligned, .efa_pooled_vec)
  boot_unrot <- vector("list", m)

  q_rot <- boot_rot <- NULL
  q_phi <- boot_phi <- NULL
  q_structure <- boot_structure <- NULL
  q_h2 <- vector("list", m)
  boot_h2 <- vector("list", m)
  q_residuals <- vector("list", m)
  boot_residuals <- vector("list", m)

  if (rotation_type != "none") {
    q_rot <- lapply(rot_loadings, .efa_pooled_vec)
    boot_rot <- vector("list", m)
  }
  if (oblique) {
    q_phi <- lapply(phis, .efa_pooled_vech)
    boot_phi <- vector("list", m)
    q_structure <- lapply(structure_loadings, .efa_pooled_vec)
    boot_structure <- vector("list", m)
  }

  nonconv_procrustes <- integer(m)
  boot_failures <- integer(m)

  for (d in seq_len(m)) {
    arr <- fits[[d]]$replicates
    unrot_arr <- arr$unrot_loadings

    boot_unrot[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars * k)

    if (rotation_type != "none") {
      boot_rot[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars * k)
    }
    if (oblique) {
      boot_phi[[d]] <- matrix(NA_real_, nrow = B_use, ncol = k * (k + 1) / 2)
      boot_structure[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars * k)
    }

    boot_h2[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars)
    boot_residuals[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars * p_vars)

    if (oblique) {
      q_h2[[d]] <- .efa_pooled_communalities(rot_loadings[[d]], phis[[d]])
      q_residuals[[d]] <- .efa_pooled_vec(.efa_pooled_residual_from_solution(
        orig_R_list[[d]], rot_loadings[[d]], phis[[d]]
      ))
    } else if (rotation_type != "none") {
      q_h2[[d]] <- .efa_pooled_communalities(rot_loadings[[d]], NULL)
      q_residuals[[d]] <- .efa_pooled_vec(.efa_pooled_residual_from_solution(
        orig_R_list[[d]], rot_loadings[[d]], NULL
      ))
    } else {
      q_h2[[d]] <- .efa_pooled_communalities(unrot_loadings_aligned[[d]], NULL)
      q_residuals[[d]] <- .efa_pooled_vec(.efa_pooled_residual_from_solution(
        orig_R_list[[d]], unrot_loadings_aligned[[d]], NULL
      ))
    }

    for (b in seq_len(B_use)) {
      Lb_unrot0 <- unrot_arr[, , b]

      # A component EFA NA-fills a bootstrap replicate it could not fit; skip it
      # (its pre-initialized rows stay NA and are dropped by the NA-tolerant
      # pooling) and tally it, rather than letting the NA reach the alignment /
      # Procrustes step and abort the whole pooled bootstrap.
      if (!all(is.finite(Lb_unrot0))) {
        boot_failures[d] <- boot_failures[d] + 1L
        next
      }

      Lb_unrot <- .efa_pooled_align_unrot_boot(
        L_boot = Lb_unrot0,
        target = mean_unrot_target,
        align_unrotated = align_unrotated
      )
      boot_unrot[[d]][b, ] <- .efa_pooled_vec(Lb_unrot)

      if (rotation_type != "none") {
        # Do not use arr$rot_loadings here. Those matrices were aligned to the
        # imputation-specific EFA target. MI pooling requires every bootstrap
        # replicate to be represented in the final MI target space.
        rot_b <- .efa_pooled_align_rot_boot(
          L_boot_unrot = Lb_unrot0,
          final_target = final_target_m,
          rotation_type = boot_proc_rotation,
          procrustes_args = procrustes_args
        )

        if (!is.null(rot_b)) {
          boot_rot[[d]][b, ] <- .efa_pooled_vec(rot_b$Lambda)

          if (oblique) {
            boot_phi[[d]][b, ] <- .efa_pooled_vech(rot_b$Phi)
            boot_structure[[d]][b, ] <- .efa_pooled_vec(rot_b$Lambda %*% rot_b$Phi)
            boot_h2[[d]][b, ] <- .efa_pooled_communalities(rot_b$Lambda, rot_b$Phi)
          } else {
            boot_h2[[d]][b, ] <- .efa_pooled_communalities(rot_b$Lambda, NULL)
          }
        } else {
          nonconv_procrustes[d] <- nonconv_procrustes[d] + 1L
        }
      } else {
        boot_h2[[d]][b, ] <- .efa_pooled_communalities(Lb_unrot, NULL)
      }

      if (!is.null(arr$residuals)) {
        boot_residuals[[d]][b, ] <- .efa_pooled_vec(arr$residuals[, , b])
      }
    }
  }

  # Insufficient implies some failed, so test it first and emit a single warning.
  if (any(B_use - boot_failures < 2L)) {
    cli::cli_warn("Too few valid bootstrap replicates remain after dropping failures; pooled bootstrap SEs/CIs could not be computed.",
                  class = "efa_pooled_boot_insufficient")
    return(NULL)
  }

  if (any(boot_failures > 0L)) {
    cli::cli_warn("Some bootstrap replicates failed in the component EFAs; pooled bootstrap SEs/CIs are based on the valid replicates.",
                  class = "efa_pooled_boot_failed")
  }

  if (any(nonconv_procrustes > 0L)) {
    cli::cli_warn("Some bootstrap Procrustes rotations did not converge; pooled bootstrap SEs/CIs are based on the valid aligned replicates.",
                  class = "efa_pooled_boot_nonconv")
  }

  pool_unrot <- .efa_pooled_rubin_pool(q_unrot, boot_unrot, alpha = alpha)
  pool_h2 <- .efa_pooled_rubin_pool(
    q_h2, boot_h2, alpha = alpha,
    est_override = as.vector(h2)
  )
  pool_residuals <- .efa_pooled_rubin_pool(
    q_residuals, boot_residuals, alpha = alpha,
    est_override = .efa_pooled_vec(residuals)
  )

  unrot_res <- .efa_pooled_rubin_matrix_result(
    pool_unrot, p_vars, k, loading_dimnames
  )
  residual_res <- .efa_pooled_rubin_matrix_result(
    pool_residuals, p_vars, p_vars, dimnames(as.matrix(residuals))
  )

  SE <- list(
    h2 = pool_h2$se,
    unrot_loadings = unrot_res$se,
    residuals = residual_res$se
  )

  CI <- list(
    h2 = list(lower = pool_h2$ci$lower, upper = pool_h2$ci$upper),
    unrot_loadings = unrot_res$ci,
    residuals = residual_res$ci
  )

  arrays <- list(
    unrot_loadings = boot_unrot,
    h2 = boot_h2,
    residuals = boot_residuals
  )

  MI <- list(
    unrot_loadings = list(RIV = pool_unrot$RIV, FMI = pool_unrot$FMI, df = pool_unrot$df),
    h2 = list(RIV = pool_h2$RIV, FMI = pool_h2$FMI, df = pool_h2$df),
    residuals = list(RIV = pool_residuals$RIV, FMI = pool_residuals$FMI, df = pool_residuals$df),
    # Per-imputation alignment outcomes. `boot_failures` are replicates the
    # component EFA could not fit (NA-filled, dropped before rotation);
    # `bootstrap_rotation_failures` are replicates whose Procrustes alignment to
    # the MI target did not produce a valid rotation. The valid count subtracts
    # both, so it never overstates the replicates that actually entered the pool.
    bootstrap_source_failures = boot_failures,
    bootstrap_rotation_failures = nonconv_procrustes,
    bootstrap_rotation_valid = B_use - boot_failures - nonconv_procrustes
  )

  if (rotation_type != "none") {
    pool_rot <- .efa_pooled_rubin_pool(q_rot, boot_rot, alpha = alpha)
    rot_res <- .efa_pooled_rubin_matrix_result(pool_rot, p_vars, k, loading_dimnames)
    SE$rot_loadings <- rot_res$se
    CI$rot_loadings <- rot_res$ci
    arrays$rot_loadings <- boot_rot
    MI$rot_loadings <- list(RIV = pool_rot$RIV, FMI = pool_rot$FMI, df = pool_rot$df)
  }

  if (oblique) {
    pool_phi <- .efa_pooled_rubin_pool(q_phi, boot_phi, alpha = alpha)
    phi_res <- .efa_pooled_rubin_symmetric_result(
      pool_phi, k, dimnames(as.matrix(phis[[1]]))
    )
    SE$Phi <- phi_res$se
    CI$Phi <- phi_res$ci
    arrays$Phi <- boot_phi
    MI$Phi <- list(RIV = pool_phi$RIV, FMI = pool_phi$FMI, df = pool_phi$df)

    pool_structure <- .efa_pooled_rubin_pool(
      q_structure, boot_structure, alpha = alpha,
      est_override = .efa_pooled_vec(mean_structure_loadings)
    )
    structure_res <- .efa_pooled_rubin_matrix_result(
      pool_structure, p_vars, k, loading_dimnames
    )
    SE$Structure <- structure_res$se
    CI$Structure <- structure_res$ci
    arrays$Structure <- boot_structure
    MI$Structure <- list(
      RIV = pool_structure$RIV,
      FMI = pool_structure$FMI,
      df = pool_structure$df
    )
  }

  ## Fit-index summaries ------------------------------------------------------
  ## fit_indices_descriptive: Rubin-Wald MI summaries of the per-imputation
  ## bootstrap fit indices. Chi-square-type fit is pooled with the D2 rule and its
  ## reference distribution supplies the point-estimate uncertainty, so no
  ## bootstrap-percentile interval is constructed for the pooled fit statistic.

  fit_q_names <- Reduce(intersect, lapply(fits, function(x) names(x$fit_indices)))
  fit_q_names <- fit_q_names[vapply(fits[[1]]$fit_indices[fit_q_names], function(x) {
    is.numeric(x) && length(x) == 1L
  }, logical(1))]

  can_describe_fit_boot <- length(fit_q_names) > 0L &&
    all(vapply(fits, function(x) {
      fit_arr <- x$replicates$fit_indices
      !is.null(fit_arr) &&
        nrow(as.matrix(fit_arr)) >= B_use &&
        ncol(as.matrix(fit_arr)) >= length(fit_q_names)
    }, logical(1)))

  if (can_describe_fit_boot) {
    q_fit <- lapply(fits, function(x) unlist(x$fit_indices[fit_q_names]))
    boot_fit <- vector("list", m)

    for (d in seq_len(m)) {
      fit_arr <- as.matrix(fits[[d]]$replicates$fit_indices)
      if (is.null(colnames(fit_arr))) {
        colnames(fit_arr) <- names(fits[[d]]$fit_indices)
      }
      boot_fit[[d]] <- fit_arr[seq_len(B_use), fit_q_names, drop = FALSE]
    }

    pool_fit_desc <- .efa_pooled_rubin_pool(q_fit, boot_fit, alpha = alpha)
    SE$fit_indices_descriptive <- stats::setNames(pool_fit_desc$se, fit_q_names)
    CI$fit_indices_descriptive <- list(
      lower = stats::setNames(pool_fit_desc$ci$lower, fit_q_names),
      upper = stats::setNames(pool_fit_desc$ci$upper, fit_q_names)
    )
    arrays$fit_indices_descriptive <- boot_fit
    MI$fit_indices_descriptive <- list(
      RIV = stats::setNames(pool_fit_desc$RIV, fit_q_names),
      FMI = stats::setNames(pool_fit_desc$FMI, fit_q_names),
      df = stats::setNames(pool_fit_desc$df, fit_q_names)
    )
  }

  list(
    SE = SE,
    CI = CI,
    replicates = arrays,
    MI = MI,
    n_boot = B_use
  )
}
