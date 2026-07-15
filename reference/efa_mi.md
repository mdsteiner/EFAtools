# Exploratory factor analysis on multiple data imputations

Fits
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
to each of several imputed datasets, aligns the factor solutions to a
common factor space, and pools the resulting estimates and selected fit
quantities across imputations.

## Usage

``` r
efa_mi(
  data_list,
  p = 0.05,
  target_method = c("first_target", "consensus"),
  align_unrotated = c("signed_tucker_congruence", "none", "procrustes"),
  fit_pool_method = c("D2"),
  consensus_args = list(),
  procrustes_args = list(),
  rmsea_ci_level = 0.9,
  rmsr_upper = TRUE,
  ...
)
```

## Arguments

- data_list:

  A list of length \\m\\, where \\m\\ is the number of imputations. Each
  list element is a data frame or matrix of raw data, or a correlation
  matrix. See argument `x` in
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).

- p:

  Numeric in \\(0, 1)\\. One minus the confidence level used for pooled
  Wald-type bootstrap/MI confidence intervals when bootstrap replicates
  are available. For example, `p = .05` gives 95% intervals.

- target_method:

  Character. How rotated solutions are aligned across imputations before
  pooling: `"first_target"` (the default) aligns every imputation to the
  first imputation's rotated solution, while `"consensus"` refines a
  centroid target by Generalized Procrustes Analysis (orthogonal
  rotations only). See *Aligning solutions across imputations* in
  Details.

- align_unrotated:

  Character. How unrotated loadings are aligned before pooling:
  `"signed_tucker_congruence"` (the default; sign/permutation via Tucker
  congruence), `"procrustes"` (orthogonal Procrustes to the first
  imputation), or `"none"`. See *Aligning solutions across imputations*
  in Details.

- fit_pool_method:

  Character. Currently only `"D2"` is implemented for chi-square-type
  fit. If no chi-square is available, only residual-based fit and
  descriptive quantities are returned. See *Pooling the model chi-square
  and fit indices* in Details.

- consensus_args:

  List of additional arguments controlling the GPA-consensus iteration
  when `target_method = "consensus"`. Recognised tuning parameters
  include the convergence tolerances `tol` and `loss_tol`, the iteration
  bounds `min_iter` and `max_iter`, the target-update damping `alpha`,
  and the multi-start controls `multi_start` and `starts`.

- procrustes_args:

  List of additional arguments passed to
  [`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)
  for fixed-target alignment.

- rmsea_ci_level:

  Numeric. Confidence level for the RMSEA CI.

- rmsr_upper:

  Logical. If `TRUE`, compute RMSR from the unique off-diagonal residual
  correlations. If `FALSE`, use the full off-diagonal matrix.

- ...:

  Additional arguments passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  (e.g. `method`, `rotation`, `se`, `n_factors`, `N`). These select the
  estimator, rotation, standard-error method, and fit indices used for
  every imputation; see
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  for the available options, their properties, and which combinations
  are valid.

## Value

A list of class `c("efa_mi", "EFA_POOLED", "efa", "EFA")` containing
pooled estimates, residuals, fit indices, the individual fits, and MI
diagnostics. The trailing legacy classes keep
`inherits(x, "EFA_POOLED")` and the single-fit EFA accessors and S3
dispatch working. In addition to the slots inherited from
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
(including `SE`, `CI`, and, on the bootstrap path, `replicates`), the
object carries:

- MI:

  Multiple-imputation diagnostics for each pooled parameter family. On
  the bootstrap path: `unrot_loadings`, `h2`, `residuals`, optionally
  `rot_loadings`, `Phi`, `Structure`, and `fit_indices_descriptive`,
  plus integer vectors `bootstrap_source_failures` (replicates the
  component
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  could not fit), `bootstrap_rotation_failures` (replicates whose
  Procrustes alignment to the target was invalid), and
  `bootstrap_rotation_valid` (those that entered the pool,
  `B - source - rotation` failures). Both paths use the plain
  Rubin (1987) df. On the analytic path (`se = "information"`):
  `unrot_loadings` and `uniquenesses`, plus, when a rotation was
  requested, `rot_loadings`, `h2`, and (oblique) `Phi` and `Structure`.
  Each per-family entry is a list with `RIV` (relative increase in
  variance), `FMI` (the fraction of missing information, reported as
  Rubin's asymptotic \\\lambda = RIV / (1 + RIV)\\, equal to
  `lavaan.mi`'s `fmi`), and `df`; the rotated families on the analytic
  path additionally carry a `method` string recording the gauge
  alignment used (`"gauge_invariant"` for communalities and
  `"signed_permutation_approx"` for rotated loadings and, for oblique
  rotations, factor correlations and structure coefficients).
  `fit_indices_descriptive`, on the bootstrap path, pools every
  per-imputation fit index, so the structural constants among them
  (`df`, `df_null`) appear with a standard error of 0.

- mi_fit:

  On the `se = "sandwich"` (MI2S) path only: the single
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fit on the pooled correlation matrix \\\bar r\\ and pooled asymptotic
  covariance \\\tilde\Gamma\\. Its `orig_R` is \\\bar r\\ and its
  `Gamma` is \\\tilde\Gamma\\; the pooled `SE`, `CI`, and `fit_indices`
  are taken from it. `MI` is `NULL` on this path because the imputation
  uncertainty is carried by \\\tilde\Gamma\\ rather than by
  per-parameter Rubin pooling.

## Details

`efa_mi()` is the multiple-imputation route to handling missing data:
several imputed datasets are each fitted with
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
and the solutions pooled. A single-fit alternative is full-information
maximum likelihood, available directly in
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
as `cor_method = "fiml"`, which EM-estimates a two-stage correlation
from one raw dataset with missing values. Both feed the same
correlation-scale EFA core and differ only in how the missingness is
handled; FIML is intentionally not routed through `efa_mi()`, which is a
multi-fit pooler by construction.

### Standard-error pooling routes

The pooling pathway is selected automatically from the `se` method
recorded on the component
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
fits, which must be identical across imputations:

- `se = "none"`: no standard errors are pooled.

- `se = "information"`: the per-imputation expected-information standard
  errors are pooled with Rubin's (1987) rules (Wald intervals).

- `se = "sandwich"`: the two-stage pooled-inputs (MI2S) approach fits a
  single model on the Rubin-pooled correlation matrix and asymptotic
  covariance.

- `se = "np-boot"`: the non-parametric bootstrap replicates are
  re-aligned to the multiple-imputation target and Rubin-pooled.

On the information and np-boot routes, if pooled standard errors cannot
be produced (for example an unreliable analytic covariance or too few
bootstrap replicates) the pool falls back to point-estimate-only pooling
and downgrades `settings$se` to `"none"`. The MI2S route is the
exception: its single fit fuses the point estimates and standard errors
through the pooled asymptotic covariance, so a structural failure aborts
directly rather than falling back.

### Aligning solutions across imputations

The same
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
model is fitted to each imputed dataset and the solutions are put into a
common factor space before averaging. For oblique solutions the factor
intercorrelations are aligned together with the loadings so the model
stays internally consistent.

`target_method` controls how rotated solutions are aligned.
`"first_target"` (the default) aligns every imputation to the first
imputation's rotated solution by one Procrustes rotation each.
`"consensus"` instead refines a centroid target by Generalized
Procrustes Analysis (Gower 1975; van Ginkel & Kroonenberg 2014;
Lorenzo-Seva & Van Ginkel 2016). The two give the same pooled estimate
for orthogonal rotations (consensus is just more expensive), and
`"consensus"` is only supported there. Anchoring on the first imputation
can understate the between-imputation variability when the imputations
disagree substantially, whereas `"consensus"` is more robust to an
atypical first imputation (van Ginkel & Kroonenberg 2014).

`align_unrotated` controls how unrotated loadings are aligned before
pooling: `"signed_tucker_congruence"` (the default) matches them up to
factor reordering and sign changes, `"procrustes"` aligns them to the
first imputation by orthogonal Procrustes rotation, and `"none"`
averages them as returned by
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).

### Pooling point estimates

Point estimates are pooled by arithmetic averaging after alignment. For
oblique rotations the structure matrix is recomputed from the pooled
pattern matrix and pooled factor correlations, \\Structure = \Lambda
\Phi\\, and communalities are the diagonal of the reproduced correlation
matrix, \\diag(\Lambda \Phi \Lambda')\\ for oblique rotations and
\\diag(\Lambda \Lambda')\\ otherwise. Residuals are not averaged across
imputations; they are the pooled observed correlation matrix minus the
model-implied correlation of the pooled solution, so RMSR/SRMR are based
on these pooled residuals. Both are returned, though the print and
summary methods show SRMR only.

### Pooling the model chi-square and fit indices

The model chi-square and the indices derived from it (RMSEA, ECVI, and
the descriptive AIC/BIC) are pooled with the D2 rule (Li, Meng,
Raghunathan & Rubin, 1991), not arithmetically averaged. Because D2
shrinks the pooled chi-square in proportion to the between-imputation
variability, the pooled RMSEA can fall below the mean of the
per-imputation RMSEAs (as it does in `lavaan.mi`); read it together with
the per-imputation fit. The incremental indices CFI (Bentler, 1990) and
TLI (Tucker & Lewis, 1973) are instead the average of the per-imputation
indices, which keeps them consistent with the component fits and avoids
the out-of-range values that separately pooling the model and baseline
noncentralities (as `lavaan.mi`/`semTools` do) can produce; those
separately pooled noncentralities remain available in `mi_diagnostics`.
AIC and BIC, if returned, are chi-square-derived descriptive quantities
and are not likelihood-based MI information criteria. On the
sandwich/MI2S route the chi-square is the single fit's scaled statistic
rather than a D2 pool.

### Bootstrap pooling (np-boot)

If each component
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
call was run with `se = "np-boot"`, pooled bootstrap SEs and Wald-type
MI confidence intervals are computed for loadings, communalities,
residuals, and, when applicable, factor correlations and structure
coefficients. The unrotated bootstrap replicates are re-aligned to the
final MI target before the within-imputation covariance is estimated,
and Rubin pooling is applied with \\T = Ubar + (1 + 1/m) B\\. The
confidence level of the pooled intervals is set by `p`, not by the
component
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
calls' `ci`.

### Analytic pooling (information)

With `se = "information"`, the analytic unrotated-loading and uniqueness
SEs returned by each fit are pooled element-wise with Rubin's rules (\\T
= Ubar + (1 + 1/m) B\\), with Wald intervals on the plain Rubin (1987)
degrees of freedom (the analytic loadings are asymptotically normal, so
the Barnard-Rubin (1999) adjustment reduces to this form, matching
`lavaan.mi`). NA propagation is fail-closed: if any imputation is NA at
an element, all pooled outputs for that element are NA. When a rotation
was requested, the rotated loadings, communalities, and (for oblique
rotations) factor correlations and structure coefficients are pooled as
well; residual SE pooling is available only on the bootstrap path. Under
`align_unrotated = "procrustes"` the full unrotated covariance
`vcov_unrot_loadings` (populated by `se = "information"`) is propagated
through the alignment, so it must be present and reliable on every fit.

A rotated-loading standard error is conditional on the rotation
criterion (Jennrich 1973, 1974; Browne 2001; Zhang & Preacher 2015). For
both orthogonal and oblique rotations the within-imputation variance is
therefore each fit's own criterion-aware delta-method rotated SE (the
quantity
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
returns), reused after a signed-permutation alignment to the MI target,
and the between-imputation variance is the sample variance of the
aligned rotated loadings. This is a deliberate approximation – each SE
is conditional on its own fit's rotation optimum rather than on a common
gauge – and is flagged by
`MI$<param>$method = "signed_permutation_approx"`. Communalities are
rotation-invariant and pool element-wise. For a fully gauge-consistent
rotated uncertainty, cross-check with `se = "np-boot"`.

### Two-stage pooling (sandwich / MI2S)

With `se = "sandwich"` (robust SEs from a polychoric/tetrachoric or
continuous-Pearson asymptotic covariance), pooling follows the
two-stage, pooled-inputs approach (Chung & Cai 2019; Sriutaisuk, Liu,
Chung, Kim & Gu 2025): the correlation matrix and the asymptotic
covariance of its off-diagonal entries are Rubin-pooled across
imputations, \$\$\tilde\Gamma = \Gamma_W + \left(1 +
\frac{1}{m}\right)\Gamma_B,\$\$ and a single `EFA` model is fitted to
the pooled correlation with \\\tilde\Gamma\\ as the robust meat (its
diagonal as the weights for `method = "DWLS"`). Because there is only
one fit and one rotational gauge, this route bypasses the per-imputation
alignment: `target_method` and `align_unrotated` do not apply. The
fitted object carries native scaled-shifted chi-square statistics and
sandwich SEs that already reflect the multiple-imputation uncertainty,
so the chi-square is not D2-pooled and the likelihood-ratio-based
AIC/BIC/ECVI are `NA`; it is returned in the `mi_fit` slot, with the
per-imputation `fits` retained for diagnostics. At least 20 imputations
are recommended for the scaled-shifted statistic, and more (around 100)
at higher rates of missingness (Sriutaisuk et al. 2025). The
polychoric/tetrachoric (ordinal) case is the primary, best-evaluated
target; the continuous-Pearson case uses the same recipe but is less
benchmarked.

## Conditions

All conditions are classed (prefix `efa_pooled_`, or `efa_consensus_`
for the consensus target) so they can be caught by class. The ones most
likely to be encountered:

- **Inputs.** `efa_pooled_min_fits` (at least two fits are required);
  `efa_pooled_mixed_se` (every imputation must use the same `se`).

- **Alignment.** `efa_consensus_oblique_unsupported`
  (`target_method = "consensus"` is orthogonal-only).

- **Standard errors.** `efa_pooled_se_unavailable` (a warning: pooled
  SEs could not be produced, so only point estimates are returned);
  `efa_pooled_no_vcov` and `efa_pooled_unreliable_vcov` (the analytic
  `"procrustes"` path needs a reliable `vcov_unrot_loadings` on every
  fit).

- **Two-stage (`se = "sandwich"`).**
  `efa_pooled_mi2s_inputs_inconsistent` (every imputation must use
  `se = "sandwich"` with the same `cor_method`);
  `efa_pooled_mi2s_n_too_small` (a warning below 20 imputations);
  `efa_pooled_mi2s_acov_not_psd` (the pooled covariance is indefinite –
  use more imputations); `efa_pooled_mi2s_alignment_ignored` (a warning
  that the alignment arguments do not apply here).

The remaining conditions concern partial or insufficient bootstrap
replicates and unequal sample sizes across imputations.

## References

Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
multiple imputation. *Biometrika*, 86(4), 948-955.

Bentler, P. M. (1990). Comparative fit indexes in structural models.
*Psychological Bulletin*, 107(2), 238-246.

Browne, M. W. (2001). An overview of analytic rotation in exploratory
factor analysis. *Multivariate Behavioral Research*, 36(1), 111-150.

Chan, K. W., & Meng, X.-L. (2022). Multiple improvements of multiple
imputation likelihood ratio tests. *Statistica Sinica*, 32, 1489-1514.

Chung, S., & Cai, L. (2019). Alternative multiple imputation inference
for categorical structural equation modeling. *Multivariate Behavioral
Research*, 54(3), 323-337.

Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika*,
40(1), 33-51.

Li, K. H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
Significance levels from repeated p-values with multiply-imputed data.
*Statistica Sinica*, 1(1), 65-92.

Jennrich, R. I. (1973). Standard errors for obliquely rotated factor
loadings. *Psychometrika*, 38(4), 593-604.

Jennrich, R. I. (1974). Simplified formulae for standard errors in
maximum-likelihood factor analysis. *British Journal of Mathematical and
Statistical Psychology*, 27(1), 122-131.

Lorenzo-Seva, U., & Van Ginkel, J. R. (2016). Multiple imputation of
missing values in exploratory factor analysis of multidimensional
scales. *Anales de Psicologia*, 32(2), 596-608.

Rubin, D. B. (1987). *Multiple imputation for nonresponse in surveys*.
Wiley.

Schoenemann, P. H. (1966). A generalized solution of the orthogonal
Procrustes problem. *Psychometrika*, 31(1), 1-10.

Sriutaisuk, S., Liu, Y., Chung, S., Kim, H., & Gu, F. (2025). Evaluating
imputation-based fit statistics in structural equation modeling with
ordinal data: The MI2S approach. *Educational and Psychological
Measurement*, 85(1), 82-113.

Tucker, L. R., & Lewis, C. (1973). A reliability coefficient for maximum
likelihood factor analysis. *Psychometrika*, 38(1), 1-10.

van Ginkel, J. R., & Kroonenberg, P. M. (2014). Using generalized
Procrustes analysis for multiple imputation in principal component
analysis. *Journal of Classification*, 31(2), 242-269.

Zhang, G., & Preacher, K. J. (2015). Factor rotation and standard errors
in exploratory factor analysis. *Journal of Educational and Behavioral
Statistics*, 40(6), 579-603.

## See also

Other factor analysis:
[`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md),
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
[`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md),
[`plot.efa_group()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_group.md),
[`print.efa_group()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_group.md)

## Author

Andreas Soteriades, Markus Steiner

## Examples

``` r

# create a list of three datasets, mimicking a list you would obtain from
# e.g. mice.
dat_list <- lapply(1:3, function(x) GRiPS_raw[sample(1:nrow(GRiPS_raw), replace = TRUE),])
mod <- efa_mi(dat_list, n_factors = 1, method = "ML")
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
mod
#> 
#> Pooled EFA across 3 imputations performed with type = 'EFAtools', method = 'ML', and rotation = 'none'.
#> Pooling settings: align_unrotated = 'signed_tucker_congruence', fit_pool_method = 'D2'.
#> 
#> ── Unrotated Loadings ──────────────────────────────────────────────────────────
#> 
#>             F1    h2    u2
#> fun        .799  .638  .362
#> friends    .843  .710  .290
#> enjoy      .879  .773  .227
#> hurt       .764  .584  .416
#> part       .821  .674  .326
#> commonly   .836  .699  .301
#> chances    .796  .634  .366
#> attracted  .844  .712  .288
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                 F1
#> SS loadings   5.423
#> Prop Tot Var   .678
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> D2-pooled χ²(20) = 0.00, p = 1.000
#> CFI: .98
#> TLI: .98
#> RMSEA [90% CI]: .00 [.00; .00]
#> AIC: -40.00
#> BIC: -133.94
#> ECVI: 0.04
#> CAF: .50
#> SRMR: .02
#> Note: the pooled χ² is the D2 statistic; its p uses the D2 reference F(20, .2),
#> not the χ²(20) tail.

# \donttest{
# add computation of standard errors and CIs
mod <- efa_mi(dat_list, n_factors = 1, method = "ML", se = "np-boot")
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
mod
#> 
#> Pooled EFA across 3 imputations performed with type = 'EFAtools', method = 'ML', and rotation = 'none'.
#> Pooling settings: align_unrotated = 'signed_tucker_congruence', fit_pool_method = 'D2'.
#> 
#> ── Unrotated Loadings ──────────────────────────────────────────────────────────
#> 
#>             F1    h2    u2
#> fun        .799  .638  .362
#> friends    .843  .710  .290
#> enjoy      .879  .773  .227
#> hurt       .764  .584  .416
#> part       .821  .674  .326
#> commonly   .836  .699  .301
#> chances    .796  .634  .366
#> attracted  .844  .712  .288
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                 F1
#> SS loadings   5.423
#> Prop Tot Var   .678
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> D2-pooled χ²(20) = 0.00, p = 1.000
#> CFI [95% bootstrap/MI-CI]: .98 [.96, 1.01]
#> TLI [95% bootstrap/MI-CI]: .98 [.94, 1.01]
#> RMSEA [90% CI] [95% bootstrap/MI-CI]: .00 [.00; .00] [.02, .12]
#> AIC [95% bootstrap/MI-CI]: -40.00 [-49.60, 186.84]
#> BIC [95% bootstrap/MI-CI]: -133.94 [-143.54, 92.90]
#> ECVI [95% bootstrap/MI-CI]: 0.04 [.03, .32]
#> CAF [95% bootstrap/MI-CI]: .50 [.46, .53]
#> SRMR [95% bootstrap/MI-CI]: .02 [.01, .03]
#> Note: the pooled χ² is the D2 statistic; its p uses the D2 reference F(20, .2),
#> not the χ²(20) tail.
#> 
#> Note: Bootstrap/MI CIs based on 1000 bootstrap samples per imputation.
# }
```
