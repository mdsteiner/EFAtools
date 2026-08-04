#' Exploratory factor analysis on multiple data imputations
#'
#' @author Andreas Soteriades, Markus Steiner
#'
#' @description
#' Fits [efa_fit()] to each of several imputed datasets, aligns the
#' factor solutions to a common factor space, and pools the resulting estimates
#' and selected fit quantities across imputations.
#'
#' @details
#' `efa_mi()` is the multiple-imputation route to handling missing data:
#' several imputed datasets are each fitted with [efa_fit()] and the solutions pooled.
#' A single-fit alternative is full-information maximum likelihood, available
#' directly in [efa_fit()] as `cor_method = "fiml"`, which EM-estimates a two-stage
#' correlation from one raw dataset with missing values. Both feed the same
#' correlation-scale EFA core and differ only in how the missingness is handled;
#' FIML is intentionally not routed through `efa_mi()`, which is a multi-fit
#' pooler by construction.
#'
#' Both routes assume the values are missing at random (MAR). Which one to prefer
#' is largely practical: FIML is a single, efficient fit and is the simpler
#' default when the analysis model is the whole story, whereas multiple imputation
#' is more flexible when the imputation model should draw on auxiliary variables
#' not in the factor model, or when the same imputations feed several downstream
#' analyses.
#'
#' ## Standard-error pooling routes
#'
#' The pooling pathway is selected automatically from the `se` method recorded on
#' the component [efa_fit()] fits, which must be identical across imputations:
#'
#' - `se = "none"`: no standard errors are pooled.
#' - `se = "information"`: the per-imputation expected-information standard errors
#'   are pooled with Rubin's (1987) rules (Wald intervals).
#' - `se = "sandwich"`: the two-stage pooled-inputs (MI2S) approach fits a single
#'   model on the Rubin-pooled correlation matrix and asymptotic covariance.
#' - `se = "np-boot"`: the non-parametric bootstrap replicates are re-aligned to
#'   the multiple-imputation target and Rubin-pooled.
#'
#' On the information and np-boot routes, if pooled standard errors cannot be
#' produced (for example an unreliable analytic covariance or too few bootstrap
#' replicates) the pool falls back to point-estimate-only pooling and downgrades
#' `settings$se` to `"none"`. The MI2S route is the exception: its single fit
#' fuses the point estimates and standard errors through the pooled asymptotic
#' covariance, so a structural failure aborts directly rather than falling back.
#'
#' ## Aligning solutions across imputations
#'
#' The same [efa_fit()] model is fitted to each imputed dataset and the solutions are
#' put into a common factor space before averaging. For oblique solutions the
#' factor intercorrelations are aligned together with the loadings so the model
#' stays internally consistent.
#'
#' `target_method` controls how rotated solutions are aligned. `"first_target"`
#' (the default) aligns every imputation to the first imputation's rotated
#' solution by one Procrustes rotation each. `"consensus"` instead refines a
#' centroid target by Generalized Procrustes Analysis (Gower 1975; van Ginkel &
#' Kroonenberg 2014; Lorenzo-Seva & Van Ginkel 2016). The two give the same pooled
#' estimate for orthogonal rotations (consensus is just more expensive), and
#' `"consensus"` is only supported there. Anchoring on the first imputation can
#' understate the between-imputation variability when the imputations disagree
#' substantially, whereas `"consensus"` is more robust to an atypical first
#' imputation (van Ginkel & Kroonenberg 2014).
#'
#' `align_unrotated` controls how unrotated loadings are aligned before pooling:
#' `"signed_tucker_congruence"` (the default) matches them up to factor reordering
#' and sign changes, `"procrustes"` aligns them to the first imputation by
#' orthogonal Procrustes rotation, and `"none"` averages them as returned by
#' [efa_fit()].
#'
#' The default anchors that matching on the *medoid* imputation -- the one closest
#' in aligned squared distance to all the others -- rather than on whichever
#' imputation happens to come first, so the pooled *unrotated* solution does not
#' depend on the order of `data_list`. The rotated solution is aligned separately,
#' against a reference chosen by `target_method`, and still depends on that
#' reference. The pooled matrix is then returned in the same gauge every
#' component fit uses, by restoring the constraint that identifies the unrotated
#' solution: diagonal \eqn{L'L} for a principal-axis extraction, diagonal
#' \eqn{L' \Psi^{-1} L} for maximum likelihood (Anderson & Rubin 1956; Lawley &
#' Maxwell 1971). Which one
#' applies is read off the component fits themselves, and a solution in neither
#' gauge -- an improper one, say -- is left as aligned. Without this step the
#' average of several aligned solutions sits in a gauge no single fit uses and
#' cannot be compared element-by-element with an [efa_fit()] solution. The
#' rotation is orthogonal and
#' common to all imputations, so communalities, the total variance accounted for,
#' the model-implied correlation matrix, the residuals and RMSR are unchanged by
#' it; only the split of variance across factors moves. `"procrustes"` and
#' `"none"` keep their first-imputation anchor and are returned as aligned.
#'
#' ## Pooling point estimates
#'
#' Point estimates are pooled by arithmetic averaging after alignment. For oblique
#' rotations the structure matrix is recomputed from the pooled pattern matrix and
#' pooled factor correlations, \eqn{Structure = \Lambda \Phi}, and communalities
#' are the diagonal of the reproduced correlation matrix,
#' \eqn{diag(\Lambda \Phi \Lambda')} for oblique rotations and
#' \eqn{diag(\Lambda \Lambda')} otherwise. Residuals are not averaged across
#' imputations; they are the pooled observed correlation matrix minus the
#' model-implied correlation of the pooled solution, so RMSR/SRMR are based on
#' these pooled residuals. Both are returned, though the print and summary methods
#' show SRMR only.
#'
#' ## Pooling the model chi-square and fit indices
#'
#' The model chi-square and the indices derived from it (ECVI and the descriptive
#' AIC/BIC) are pooled with the D2 rule (Li, Meng, Raghunathan &
#' Rubin, 1991), not arithmetically averaged. RMSEA is pooled by the same rule but
#' from a second D2 pool of the per-imputation discrepancies taken on the
#' uncorrected \eqn{N - 1} scale, because RMSEA places the model noncentrality on
#' the scale on which it is defined, so the Bartlett small-sample correction enters
#' only the chi-square test and not this approximate-fit index (as in
#' [efa_fit()]). The printed RMSEA therefore does not reconcile by hand with the
#' printed chi-square; the statistic it is formed from is `chi_cfi` in
#' `mi_diagnostics`. Because D2 shrinks the pooled
#' chi-square in proportion to the between-imputation variability, the pooled RMSEA
#' can fall below the mean of the per-imputation RMSEAs (as it does in
#' `lavaan.mi`); read it together with the per-imputation fit. The incremental
#' indices CFI (Bentler, 1990) and TLI (Tucker & Lewis, 1973) are instead the
#' average of the per-imputation indices, which keeps them consistent with the
#' component fits and avoids the out-of-range values that separately pooling the
#' model and baseline noncentralities (as `lavaan.mi`/`semTools` do) can produce;
#' the separately pooled model and baseline chi-squares those indices would be
#' formed from remain available in `mi_diagnostics`.
#' AIC and BIC, if returned, are chi-square-derived descriptive quantities and are
#' not likelihood-based MI information criteria. They are reported only where the
#' component fits report them: whenever a component withholds them -- any
#' `cor_method = "fiml"` fit, and any fit whose chi-square is a scaled statistic,
#' such as `se = "sandwich"` -- the pooled AIC, BIC, and ECVI are `NA` too, matching
#' what [efa_fit()] returns for a single such fit. On the sandwich/MI2S route the
#' chi-square is the single fit's scaled statistic rather than a D2 pool.
#'
#' ## Bootstrap pooling (np-boot)
#'
#' If each component [efa_fit()] call was run with `se = "np-boot"`, pooled bootstrap SEs
#' and Wald-type MI confidence intervals are computed for loadings, communalities,
#' residuals, and, when applicable, factor correlations and structure
#' coefficients. The unrotated bootstrap replicates are re-aligned to the final MI
#' target before the within-imputation covariance is estimated, and Rubin pooling
#' is applied with \eqn{T = Ubar + (1 + 1/m) B}. The confidence level of the pooled
#' intervals is set by `p`, not by the component [efa_fit()] calls' `ci`.
#'
#' ## Analytic pooling (information)
#'
#' With `se = "information"`, the analytic unrotated-loading and uniqueness SEs
#' returned by each fit are pooled element-wise with Rubin's rules
#' (\eqn{T = Ubar + (1 + 1/m) B}), with Wald intervals on the plain Rubin (1987)
#' degrees of freedom (the analytic loadings are asymptotically normal, so the
#' Barnard-Rubin (1999) adjustment reduces to this form, matching `lavaan.mi`). NA
#' propagation is fail-closed: if any imputation is NA at an element, all pooled
#' outputs for that element are NA. When a rotation was requested, the rotated
#' loadings, communalities, and (for oblique rotations) factor correlations and
#' structure coefficients are pooled as well; residual SE pooling is available only
#' on the bootstrap path. Under `align_unrotated = "procrustes"` the full unrotated
#' covariance `vcov_unrot_loadings` (populated by `se = "information"`) is
#' propagated through the alignment, so it must be present and reliable on every
#' fit. The default alignment also mixes loading columns, through the common
#' canonical-gauge rotation, and so propagates the same covariance; where a fit
#' does not carry it, the unrotated standard errors are returned as `NA` rather
#' than aborting, and the remaining families still pool.
#'
#' A rotated-loading standard error is conditional on the rotation criterion
#' (Archer & Jennrich 1973; Jennrich 1973, 1974; Zhang, Preacher, & Jennrich
#' 2012; Zhang & Preacher 2015). For both orthogonal
#' and oblique rotations the within-imputation variance is therefore each fit's own
#' criterion-aware delta-method rotated SE (the quantity `efa_fit()` returns), reused
#' after a signed-permutation alignment to the MI target, and the
#' between-imputation variance is the sample variance of the aligned rotated
#' loadings. This is a deliberate approximation -- each SE is conditional on its
#' own fit's rotation optimum rather than on a common gauge -- and is flagged by
#' `MI$<param>$method = "signed_permutation_approx"`. Communalities are
#' rotation-invariant and pool element-wise. For a fully gauge-consistent rotated
#' uncertainty, cross-check with `se = "np-boot"`.
#'
#' ## Two-stage pooling (sandwich / MI2S)
#'
#' With `se = "sandwich"` (robust SEs from a polychoric/tetrachoric or
#' continuous-Pearson asymptotic covariance), pooling follows the two-stage,
#' pooled-inputs approach (Chung & Cai 2019; Sriutaisuk, Liu, Chung, Kim & Gu
#' 2025): the correlation matrix and the asymptotic covariance of its off-diagonal
#' entries are Rubin-pooled across imputations,
#' \deqn{\tilde\Gamma = \Gamma_W + \left(1 + \frac{1}{m}\right)\Gamma_B,}
#' and a single `EFA` model is fitted to the pooled correlation with
#' \eqn{\tilde\Gamma} as the robust meat (its diagonal as the weights for
#' `estimator = "DWLS"`). Because there is only one fit and one rotational gauge, this
#' route bypasses the per-imputation alignment: `target_method` and
#' `align_unrotated` do not apply. The fitted object carries native scaled-shifted
#' chi-square statistics and sandwich SEs that already reflect the
#' multiple-imputation uncertainty, so the chi-square is not D2-pooled and the
#' likelihood-ratio-based AIC/BIC/ECVI are `NA`; it is returned in the `mi_fit`
#' slot, with the per-imputation `fits` retained for diagnostics. The pooled fit
#' uses the same [estimate_control()] and [rotate_control()] tuning (including
#' any rotation-engine extras) as the per-imputation fits. At least 20
#' imputations are recommended for the scaled-shifted statistic, and more (around
#' 100) at higher rates of missingness (Sriutaisuk et al. 2025). The
#' polychoric/tetrachoric (ordinal) case is the primary, best-evaluated target; the
#' continuous-Pearson case uses the same recipe but is less benchmarked.
#'
#' @section Conditions:
#' All conditions are classed (prefix `efa_pooled_`, or `efa_consensus_` for the
#' consensus target; the dots validation shared with [efa_fit()] signals
#' `efa_flat_knob_in_dots` and `efa_renamed_arg`) so they can be caught by
#' class. The ones most likely to be encountered:
#'
#' - **Inputs.** `efa_pooled_min_fits` (at least two fits are required);
#'   `efa_pooled_mixed_se` (every imputation must use the same `se`).
#' - **Alignment.** `efa_consensus_oblique_unsupported` (`target_method =
#'   "consensus"` is orthogonal-only).
#' - **Standard errors.** `efa_pooled_se_unavailable` (a warning: pooled SEs could
#'   not be produced, so only point estimates are returned); `efa_pooled_no_vcov`
#'   and `efa_pooled_unreliable_vcov` (the analytic `"procrustes"` path needs a
#'   reliable `vcov_unrot_loadings` on every fit).
#' - **Two-stage (`se = "sandwich"`).** `efa_pooled_mi2s_inputs_inconsistent`
#'   (every imputation must use `se = "sandwich"` with the same `cor_method`);
#'   `efa_pooled_mi2s_n_too_small` (a warning below 20 imputations);
#'   `efa_pooled_mi2s_acov_not_psd` (the pooled covariance is indefinite -- use
#'   more imputations); `efa_pooled_mi2s_alignment_ignored` (a warning that the
#'   alignment arguments do not apply here).
#'
#' The remaining conditions concern partial or insufficient bootstrap replicates
#' and unequal sample sizes across imputations.
#'
#' @param data_list A list of length \eqn{m}, where \eqn{m} is the number of
#' imputations. Each list element is a data frame or matrix of raw data, or a
#' correlation matrix. See argument `x` in [efa_fit()].
#' @param p Numeric in \eqn{(0, 1)}. One minus the confidence level used for
#' pooled Wald-type bootstrap/MI confidence intervals when bootstrap replicates
#' are available. For example, `p = .05` gives 95% intervals.
#' @param target_method Character. How rotated solutions are aligned across imputations
#' before pooling: `"first_target"` (the default) aligns every imputation to the first
#' imputation's rotated solution, while `"consensus"` refines a centroid target by
#' Generalized Procrustes Analysis (orthogonal rotations only). See *Aligning solutions
#' across imputations* in Details.
#' @param align_unrotated Character. How unrotated loadings are aligned before pooling:
#' `"signed_tucker_congruence"` (the default; sign/permutation via Tucker congruence,
#' anchored on the medoid imputation and returned in the extraction's canonical
#' gauge), `"procrustes"` (orthogonal Procrustes to the first imputation), or
#' `"none"`. See *Aligning solutions across imputations* in Details.
#' @param fit_pool_method Character. Currently only `"D2"` is implemented
#' for chi-square-type fit. If no chi-square is available, only residual-based
#' fit and descriptive quantities are returned. See *Pooling the model chi-square and
#' fit indices* in Details.
#' @param consensus_args List of additional arguments controlling the
#' GPA-consensus iteration when `target_method = "consensus"`. Recognised tuning
#' parameters include the convergence tolerances `tol` and `loss_tol`, the
#' iteration bounds `min_iter` and `max_iter`, the target-update damping `alpha`,
#' and the multi-start controls `multi_start` and `starts`.
#' @param procrustes_args List of additional arguments passed to [efa_procrustes()]
#' for fixed-target alignment.
#' @param rmsea_ci_level Numeric. Confidence level for the RMSEA CI.
#' @param rmsr_upper `r lifecycle::badge("deprecated")` Accepted and ignored. It
#' selected between computing RMSR from the unique off-diagonal residual correlations
#' and from the full off-diagonal matrix. The two element sets hold each residual pair
#' once and twice respectively, so their sums and counts double together and the mean
#' square is the same number whenever the residual matrix is symmetric, which the pooled
#' residuals always are. RMSR is therefore always the root mean square of the unique
#' off-diagonal residuals; SRMR, which divides the same sum by the number of
#' non-redundant elements, is reported alongside it. Supplying it to [efa_mi()]
#' signals a deprecation warning; the superseded [EFA_POOLED()] accepts it silently.
#' @param ... Additional arguments passed to [efa_fit()] (e.g. `estimator`, `rotation`, `se`,
#' `n_factors`, `N`). These select the estimator, rotation, standard-error method, and
#' fit indices used for every imputation; see [efa_fit()] for the available options, their
#' properties, and which combinations are valid.
#'
#' @return A list of class `c("efa_mi", "EFA_POOLED", "efa", "EFA")` containing
#' pooled estimates, residuals, fit indices, the individual fits, and MI
#' diagnostics. The trailing legacy classes keep `inherits(x, "EFA_POOLED")` and
#' the single-fit EFA accessors and S3 dispatch working. In
#' addition to the slots inherited from [efa_fit()] (including `SE`, `CI`, and,
#' on the bootstrap path, `replicates`), the object carries:
#' \describe{
#' \item{MI}{Multiple-imputation diagnostics for each pooled parameter family.
#' On the bootstrap path: `unrot_loadings`, `h2`, `residuals`, optionally
#' `rot_loadings`, `Phi`, `Structure`, and `fit_indices_descriptive`, plus
#' integer vectors `bootstrap_source_failures` (replicates the component [efa_fit()]
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
#' \item{mi_fit}{On the `se = "sandwich"` (MI2S) path only: the single [efa_fit()]
#' fit on the pooled correlation matrix \eqn{\bar r} and pooled asymptotic
#' covariance \eqn{\tilde\Gamma}. Its `orig_R` is \eqn{\bar r} and its `Gamma`
#' is \eqn{\tilde\Gamma}; the pooled `SE`, `CI`, and `fit_indices` are taken
#' from it. `MI` is `NULL` on this path because the imputation uncertainty is
#' carried by \eqn{\tilde\Gamma} rather than by per-parameter Rubin pooling.}
#' \item{mi_diagnostics}{Diagnostics for the pooled model fit, `NULL` on the
#' `se = "sandwich"` (MI2S) path, where there is one fit and no D2 pool. `m` is the
#' number of imputations that entered the pool. `D2_F`, `D2_df1`, `D2_df2`,
#' `D2_chi_asymptotic`, `ARIV` and `FMI` describe the D2 pool of the model
#' chi-square (the average relative increase in variance and the fraction of
#' missing information it implies), and `chi_bar_naive` is the plain mean of the
#' per-imputation statistics for comparison; the `*_null` entries are the same
#' quantities for the independence baseline. `chi_cfi` and `chi_null_cfi` are the
#' pooled model and baseline **chi-squares** on the common \eqn{N - 1}
#' noncentrality scale. `chi_cfi` is the statistic the reported RMSEA is formed
#' from; the pair is also the basis on which `lavaan.mi`/`semTools` form pooled
#' incremental indices, so the reference CFI is
#' `1 - (chi_cfi - df) / (chi_null_cfi - df_null)` (and analogously for TLI) --
#' a different quantity from the reported CFI/TLI, which average the
#' per-imputation indices.}
#' \item{fits}{The list of \eqn{m} component [efa_fit()] fits, in the order of
#' `data_list`, kept for per-imputation diagnostics. On the MI2S path these are the
#' per-imputation fits whose inputs were pooled, not the pooled fit itself (which
#' is `mi_fit`).}
#' \item{alignment}{Metadata from aligning the rotated solutions, `NULL` when no
#' rotation was requested or on the MI2S path (one fit, one gauge). Under
#' `target_method = "first_target"`: the `method` used, the `target` it aligned to,
#' the per-imputation `target_rotations`, the indices of any
#' `point_rotation_failures`, and whether every inner alignment `converged`. Under
#' `target_method = "consensus"` it is the full [efa_procrustes()]-based GPA record:
#' the converged `target`, the `aligned_loadings` and `aligned_phi`, the iteration
#' `history`, convergence flags, and the multi-start summary.}
#' \item{settings}{The component fits' [efa_fit()] settings with the pooling
#' settings added: `pooled` (always `TRUE`), `pooled_N` and `N` (the mean N across
#' imputations), `n_imputations`, `component_se` (the `se` the component fits
#' used), `target_method`, `align_unrotated`, `fit_pool_method`, `p`, `ci` and
#' `rmsea_ci_level`. `se` records what was actually pooled, so it
#' is `"none"` when pooled standard errors could not be produced although the
#' component fits computed them (`component_se` keeps the request).}
#' }
#'
#'
#' @references
#' Anderson, T. W., & Rubin, H. (1956). Statistical inference in factor analysis.
#' In *Proceedings of the Third Berkeley Symposium on Mathematical Statistics and
#' Probability* (Vol. 5, pp. 111-150). University of California Press.
#'
#' Archer, C. O., & Jennrich, R. I. (1973). Standard errors for rotated factor
#' loadings. *Psychometrika*, 38(4), 581-592.
#'
#' Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
#' multiple imputation. *Biometrika*, 86(4), 948-955.
#'
#' Bentler, P. M. (1990). Comparative fit indexes in structural models.
#' *Psychological Bulletin*, 107(2), 238-246.
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
#' Zhang, G., Preacher, K. J., & Jennrich, R. I. (2012). The infinitesimal
#' jackknife with exploratory factor analysis. *Psychometrika*, 77(4), 634-648.
#'
#' @family factor analysis
#'
#' @export
#'
#' @examples
#'
#' # create a list of three datasets, mimicking a list you would obtain from
#' # e.g. mice.
#' dat_list <- lapply(1:3, function(x) GRiPS_raw[sample(1:nrow(GRiPS_raw), replace = TRUE),])
#' mod <- efa_mi(dat_list, n_factors = 1, estimator = "ML")
#' mod
#'
#' \donttest{
#' # add computation of standard errors and CIs
#' mod <- efa_mi(dat_list, n_factors = 1, estimator = "ML", se = "np-boot")
#' mod
#' }
efa_mi <- function(data_list,
                   p = 0.05,
                   target_method = c("first_target", "consensus"),
                   align_unrotated = c("signed_tucker_congruence", "none", "procrustes"),
                   fit_pool_method = c("D2"),
                   consensus_args = list(),
                   procrustes_args = list(),
                   rmsea_ci_level = .90,
                   rmsr_upper = lifecycle::deprecated(),
                   ...) {

  efa_args <- list(...)

  # A flat tuning knob or the former `method` spelling in the dots would only surface from
  # inside the first per-imputation fit, with the error attributed to efa_fit(); reject it
  # here so the message names the function the user called.
  .reject_flat_knobs(...names(), fn = "efa_mi")

  # The per-imputation do.call(efa_fit, ...) binds the dots by R's argument matching,
  # which accepts unique abbreviations of efa_fit()'s formals (e.g. `rotate_c =` for
  # `rotate_control`), but the pooling routes below read `efa_args` by exact name.
  # Canonicalize abbreviated names up front with the same exact-first, unique-prefix
  # matching do.call() applies; ambiguous or unknown names are left untouched for
  # efa_fit() to reject.
  arg_names <- names(efa_args)
  if (length(arg_names)) {
    efa_fit_formals <- setdiff(names(formals(efa_fit)), "...")
    canonical <- pmatch(arg_names, efa_fit_formals)
    names(efa_args)[!is.na(canonical)] <- efa_fit_formals[canonical[!is.na(canonical)]]
  }

  checkmate::assert_list(data_list, null.ok = FALSE)
  # Arity is asserted here rather than through checkmate's `min.len`, so that too few
  # imputations signal the documented `efa_pooled_min_fits` condition instead of a
  # plain checkmate error.
  if (length(data_list) < 2L) {
    cli::cli_abort(
      c("At least two imputations are required for MI pooling.",
        "i" = "{.arg data_list} has {length(data_list)} element{?s}."),
      class = "efa_pooled_min_fits"
    )
  }
  lapply(data_list, checkmate::assert_multi_class, c("matrix", "data.frame"))
  checkmate::assert_number(p, na.ok = FALSE, lower = 0, upper = 1)
  checkmate::assert_list(consensus_args, null.ok = FALSE)
  checkmate::assert_list(procrustes_args, null.ok = FALSE)
  checkmate::assert_number(rmsea_ci_level, na.ok = FALSE, lower = 0, upper = 1)

  # `rmsr_upper` never changed a returned value; the argument's documentation gives the
  # argument, and the warning below repeats it for the caller.
  if (lifecycle::is_present(rmsr_upper)) {
    lifecycle::deprecate_warn(
      when = "1.1.0",
      what = "efa_mi(rmsr_upper)",
      details = c(
        "i" = paste("It selected between two RMSR conventions that coincide for a symmetric",
                    "residual matrix, so it never affected the reported RMSR."),
        "i" = paste("RMSR is the root mean square of the unique off-diagonal residuals; SRMR,",
                    "which divides the same sum by the non-redundant elements, is reported",
                    "alongside it.")
      )
    )
  }

  if (p <= 0 || p >= 1) {
    cli::cli_abort("{.arg p} must be strictly between 0 and 1.", class = "efa_pooled_bad_p")
  }
  if (rmsea_ci_level <= 0 || rmsea_ci_level >= 1) {
    cli::cli_abort("{.arg rmsea_ci_level} must be strictly between 0 and 1.", class = "efa_pooled_bad_ci_level")
  }
  if (!is.null(efa_args$ci) && length(efa_args$ci) == 1L &&
      is.finite(efa_args$ci) &&
      abs(as.numeric(efa_args$ci) - (1 - p)) > sqrt(.Machine$double.eps)) {
    cli::cli_warn("{.fn efa_mi} uses {.arg p}, not the component {.fn efa_fit} argument {.arg ci}, to set pooled bootstrap/MI confidence intervals.",
                  class = "efa_pooled_ci_ignored")
  }

  target_method <- .match_arg_ci(target_method)
  align_unrotated <- .match_arg_ci(align_unrotated)
  fit_pool_method <- .match_arg_ci(fit_pool_method)

  m_imp <- length(data_list)

  ## -------------------------------------------------------------------------
  ## Fit EFA to each imputed dataset
  ## -------------------------------------------------------------------------

  fits <- lapply(data_list, function(data_list_subset) {
    do.call(efa_fit, c(list(x = data_list_subset), efa_args))
  })

  .efa_pooled_check_fits(fits)

  # Select the SE-pooling pathway from the component fits' shared se method. A
  # heterogeneous se cannot be pooled into a single MI estimate, so abort rather
  # than silently producing an uninterpretable mixture.
  route <- .efa_pooled_route(fits)
  if (identical(route, "mixed")) {
    cli::cli_abort(
      c("The component {.fn efa_fit} fits use different {.arg se} methods, so their standard errors cannot be pooled.",
        "i" = "Re-fit every imputation with the same {.arg se} (all {.val none}, {.val information}, {.val sandwich}, or {.val np-boot})."),
      class = "efa_pooled_mixed_se"
    )
  }

  settings <- fits[[1]]$settings
  estimator <- settings$estimator
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
      settings = settings, estimator = estimator, rotation = rotation,
      rotation_type = rotation_type, target_method = target_method,
      align_unrotated = align_unrotated, fit_pool_method = fit_pool_method,
      p = p, rmsea_ci_level = rmsea_ci_level
    ))
  }

  ## -------------------------------------------------------------------------
  ## Extract and align unrotated loadings
  ## -------------------------------------------------------------------------

  unrot_loadings <- lapply(.extract_list_object(fits, "unrot_loadings"),
                           .change_class, "matrix")

  aligned_unrot <- .efa_pooled_align_unrotated_list(
    unrot_loadings = unrot_loadings,
    align_unrotated = align_unrotated
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
          efa_procrustes,
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
            "i" = "The pooled point estimates still use the best available {.fn efa_procrustes} alignment; inspect {.code alignment$point_rotation_failures}."),
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

  RMSR <- .rmsr(residuals)

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
    if (!identical(estimator, "PAF")) {
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

  pooled_fit <- .efa_pooled_fit_indices(
    fits = fits,
    pooled_R = pooled_orig_R,
    residuals = residuals,
    RMSR = RMSR,
    N = N_pool,
    Ns = Ns,
    pool_method = fit_pool_method,
    rmsea_ci_level = rmsea_ci_level
  )
  fit_indices <- pooled_fit$fit_indices
  mi_diagnostics <- pooled_fit$mi_diagnostics

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
      alpha = p
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

  mean_unrot_loadings <- .change_class(mean_unrot_loadings,
                                       c("efa_loadings", "LOADINGS"))
  if (!is.null(mean_rot_loadings)) {
    mean_rot_loadings <- .change_class(mean_rot_loadings,
                                       c("efa_loadings", "LOADINGS"))
  }
  if (!is.null(mean_structure_loadings)) {
    mean_structure_loadings <- .change_class(mean_structure_loadings,
                                             c("efa_loadings", "LOADINGS"))
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
  # requested method so print.efa_mi() can still describe what the component
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
    mi_diagnostics = mi_diagnostics
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

  # The trailing legacy classes are load-bearing: "EFA_POOLED" keeps
  # `inherits(x, "EFA_POOLED")` resolving, and "efa"/"EFA" keep the shared print,
  # format, summary, and residuals methods (and the single-fit accessors) working.
  class(results) <- c("efa_mi", "EFA_POOLED", "efa", "EFA")
  results
}

## =============================================================================
## Internal helpers for efa_mi
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
  # pooled". Emitted whenever efa_mi() falls back to point-estimate-only
  # pooling on the information / np-boot routes, layered over any specific
  # condition the pooling helper raised (carried as `parent` for a classed
  # abort, NULL when the helper soft-returned NULL after its own warning).
  cli::cli_warn(
    c("Pooled standard errors could not be produced for the {.val {route}} standard-error method; the pooled point estimates are returned without standard errors.",
      "i" = "Inspect the component {.fn efa_fit} fits and any preceding condition for the cause."),
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

.efa_pooled_mi2s <- function(fits, data_list, efa_args, settings, estimator,
                             rotation, rotation_type, target_method,
                             align_unrotated, fit_pool_method, p,
                             rmsea_ci_level) {
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
    # The settings can be perfectly consistent and still unsupported here (e.g.
    # cor_method = "fiml", whose sandwich carries no poolable correlation ACOV);
    # advising the user to "re-fit with the same settings" would then be wrong.
    unsupported_cor <- !anyNA(ses) && all(ses == "sandwich") &&
      length(unique(cor_methods)) == 1L &&
      !cor_methods[[1L]] %in% c("poly", "tetra", "pearson")
    hint <- if (unsupported_cor) {
      c("i" = "The imputations were fitted with {.code cor_method = \"{cor_methods[[1L]]}\"}, which carries no poolable asymptotic correlation covariance; use {.code se = \"np-boot\"} to pool such fits.")
    } else {
      c("i" = "Re-fit every imputation with the same {.code se = \"sandwich\"} setting and {.arg cor_method}.")
    }
    cli::cli_abort(
      c("MI2S pooling requires every imputation to be fitted with {.code se = \"sandwich\"} and the same {.arg cor_method} ({.val poly}, {.val tetra}, or {.val pearson}), each carrying a valid asymptotic covariance.",
        hint),
      class = "efa_pooled_mi2s_inputs_inconsistent"
    )
  }
  cor_method <- cor_methods[[1L]]
  n_factors <- settings$n_factors
  use_setting <- settings$use
  if (is.null(use_setting)) use_setting <- "pairwise.complete.obs"

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

  weights <- if (estimator == "DWLS") {
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

  # Reuse the same estimate -> rotate -> SE pipeline efa_fit() runs, with the same
  # estimation/rotation tuning the component fits used: unbundle the two control
  # objects from the component call (defaulting to the constructors' presets,
  # exactly as efa_fit() does) and pass the knobs to .efa_core() under its
  # historical formal names. The arguments resolved explicitly here are dropped
  # from the dots to avoid passing them twice through do.call(); the remaining
  # dots are merged over the rotation-engine extras stored in rotate_control(),
  # with the dots winning on a name clash (mirroring efa_fit()).
  ec <- efa_args$estimate_control
  if (is.null(ec)) ec <- estimate_control()
  rc <- efa_args$rotate_control
  if (is.null(rc)) rc <- rotate_control()
  drop_args <- c("n_factors", "N", "estimator", "rotation", "se", "use",
                 "cor_method", "ci", "seed",
                 "estimate_control", "rotate_control")
  extra_args <- modifyList(rc$extra_args,
                           efa_args[setdiff(names(efa_args), drop_args)])
  core_args <- list(R = r_bar, N = N_pool, weights = weights, Gamma = Gamma_tilde,
                    np_boot = FALSE, estimator = estimator, rotation = rotation,
                    type = ec$type, rot_type = rc$type, n_factors = n_factors,
                    se = "sandwich", ci = 1 - p, use = use_setting,
                    cor_method = cor_method, max_iter = ec$max_iter,
                    init_comm = ec$init_comm, criterion = ec$criterion,
                    criterion_type = ec$criterion_type, abs_eigen = ec$abs_eigen,
                    start_method = ec$start_method, normalize = rc$normalize,
                    precision = rc$precision, order_type = rc$order_type,
                    varimax_type = rc$varimax_type, P_type = rc$p_type, k = rc$k,
                    randomStarts = rc$random_starts)
  # Mirror efa_fit()'s extras splice: an extra whose name is already an explicit
  # core argument would abort do.call() with "matched by multiple actual
  # arguments", so drop it (the resolved value wins).
  extra_args <- extra_args[setdiff(names(extra_args), names(core_args))]
  mi_fit <- do.call(.efa_core, c(core_args, extra_args))

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

  # Same load-bearing legacy classes as the Rubin/bootstrap return above.
  class(results) <- c("efa_mi", "EFA_POOLED", "efa", "EFA")
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
  # sense when all imputations estimate the same model on the same variables. The
  # arity of the list is already asserted on `data_list` by the caller.
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

  for (nm in c("estimator", "rotation", "n_factors")) {
    vals <- vapply(fits, .efa_pooled_setting_chr, character(1), name = nm)
    vals <- vals[!is.na(vals)]
    if (length(unique(vals)) > 1L) {
      cli::cli_abort("All imputations must use the same {.arg {nm}}.", class = "efa_pooled_setting_mismatch")
    }
  }

  invisible(TRUE)
}
