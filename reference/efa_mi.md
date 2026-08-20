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
  rmsr_upper = lifecycle::deprecated(),
  ...
)
```

## Arguments

- data_list:

  A list of length \\m\\, where \\m\\ is the number of imputations. Each
  list element is a data frame or matrix of raw data, or a correlation
  matrix. See argument `x` in
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).
  A `mids` object from mice must be converted first, with
  `mice::complete(x, "all")`.

- p:

  Numeric in \\(0, 1)\\. One minus the confidence level for the pooled
  confidence intervals, whichever `se` method produced them
  (`"information"`, `"np-boot"`, or `"sandwich"`). For example,
  `p = .05` gives 95% intervals.

- target_method:

  Character. How rotated solutions are aligned across imputations before
  pooling: `"first_target"` (the default) aligns every imputation to the
  first imputation's rotated solution, while `"consensus"` refines a
  centroid target by Generalized Procrustes Analysis, started from the
  medoid imputation so that the pooled rotated solution does not depend
  on the order of `data_list` (orthogonal rotations only). See *Aligning
  solutions across imputations* in Details.

- align_unrotated:

  Character. How unrotated loadings are aligned before pooling:
  `"signed_tucker_congruence"` (the default; sign/permutation via Tucker
  congruence, anchored on the medoid imputation and returned in the
  extraction's canonical gauge), `"procrustes"` (orthogonal Procrustes
  to the first imputation), or `"none"`. See *Aligning solutions across
  imputations* in Details.

- fit_pool_method:

  Character. Only `"D2"` is implemented for pooling chi-square-type fit.
  If no chi-square is available, only residual-based fit and descriptive
  quantities are returned. See *Pooling the model chi-square and fit
  indices* in Details.

- consensus_args:

  List of additional arguments controlling the GPA-consensus iteration
  when `target_method = "consensus"`. Recognised tuning parameters
  include the convergence tolerances `tol` and `loss_tol`, the iteration
  bounds `min_iter` and `max_iter`, the target-update damping `alpha`,
  the multi-start controls `multi_start` and `starts`, and `start`,
  which overrides the medoid imputation the iteration is otherwise
  started from.

- procrustes_args:

  List of
  [`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)
  algorithm controls for fixed-target alignment, for example
  `oblique_maxit` or `oblique_random_starts`. The loadings `A`, the
  alignment `Target`, the `rotation` family, and the cross-product `S`
  are derived from the imputations and cannot be set here.

- rmsea_ci_level:

  Numeric. Confidence level for the RMSEA CI.

- rmsr_upper:

  **\[deprecated\]** Deprecated and ignored. `efa_mi()` now always
  computes RMSR the same way, from the unique off-diagonal residuals;
  SRMR is reported alongside it. Supplying it to `efa_mi()` signals a
  deprecation warning; the superseded
  [`EFA_POOLED()`](https://mdsteiner.github.io/EFAtools/reference/EFA_POOLED.md)
  accepts it silently.

- ...:

  Additional arguments passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  (e.g. `estimator`, `rotation`, `se`, `n_factors`, `N`). These select
  the estimator, rotation, standard-error method, and fit indices used
  for every imputation; see
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  for the available options, their properties, and which combinations
  are valid. Two of them shape the pooled object rather than a single
  fit: `seed` sets the random state once for the whole `efa_mi()` call –
  every component bootstrap and every random-start rotation draws from
  it, so a seeded call is reproducible as a whole, and the caller's
  random stream is restored afterwards – and `b_boot` sets the number of
  bootstrap replicates drawn per imputation under `se = "np-boot"`,
  which is what the pooled within-imputation variances are estimated
  from and is recorded in `settings$b_boot`. The
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  and
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  objects are accepted through `...` as well, although they are not
  declared formals: pass them as `estimate_control =` /
  `rotate_control =` exactly as you would to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).

## Value

A list of class `c("efa_mi", "EFA_POOLED", "efa", "EFA")` containing
pooled estimates, residuals, fit indices, the individual fits, and MI
diagnostics. The trailing legacy classes keep
`inherits(x, "EFA_POOLED")` and the single-fit EFA accessors and S3
dispatch working. In addition to the slots inherited from
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
(including `SE`, `CI`, and, on the bootstrap path, `replicates`), the
object carries:

- SE, CI:

  Pooled standard errors and confidence intervals, named as
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  names them: where a pooled communality standard error and interval are
  produced, they are `SE$communalities` and `CI$communalities` on every
  route. The Rubin routes (`se = "information"`, `se = "np-boot"`)
  additionally return them under the compatibility alias `h2`, which
  holds the same values. The analytic route builds the communality
  family only when a rotation was requested; an unrotated analytic pool
  reports `uniquenesses` instead.

- fit_indices:

  The pooled fit indices. Every route reports `chi`, `df`, `p_chi`,
  `CAF`, `CFI`, `TLI`, `RMSEA`, `RMSEA_LB`, `RMSEA_UB`, `AIC`, `BIC`,
  `ECVI`, `RMSR`, `SRMR`, `chi_null`, `df_null`, `p_null`, and
  `pool_method` under those names and in that order. `pool_method`
  records the rule the model chi-square was pooled with (`"D2"`); it is
  `NA` on the `se = "sandwich"` (MI2S) path, which fits once on the
  pooled inputs and reports that fit's own scaled statistic rather than
  pooling several. On that path a few extra scaled-statistic fields are
  appended after the common block, for advanced diagnostics.

- standardized_residuals:

  The pooled residuals divided by their pooled bootstrap standard
  errors, with a zero diagonal. Returned on the `se = "np-boot"` path
  only, the one route that pools a residual standard error.

- MI:

  Multiple-imputation diagnostics for each pooled parameter family. On
  the bootstrap path: `unrot_loadings`, `communalities`, `residuals`,
  optionally `rot_loadings`, `Phi`, `Structure`, and
  `fit_indices_descriptive`, plus integer vectors
  `bootstrap_source_failures` (replicates the component
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  could not fit), `bootstrap_rotation_failures` (replicates whose
  Procrustes alignment to the target was invalid), and
  `bootstrap_rotation_valid` (those that entered the pool,
  `B - source - rotation` failures). Both paths use the plain
  Rubin (1987) df. On the analytic path (`se = "information"`):
  `unrot_loadings` and `uniquenesses`, plus, when a rotation was
  requested, `rot_loadings`, `communalities`, and (oblique) `Phi` and
  `Structure`. The communality family is keyed by its canonical name
  here, without the `SE`/`CI` alias, so each family is counted once in
  the printed FMI/RIV summary. Each per-family entry is a list with
  `RIV` (relative increase in variance), `FMI` (the fraction of missing
  information, reported as Rubin's asymptotic \\\lambda = RIV / (1 +
  RIV)\\), and `df`; the rotated families on the analytic path
  additionally carry a `method` string recording the gauge alignment
  used (`"gauge_invariant"` for communalities and
  `"signed_permutation_approx"` for rotated loadings and, for oblique
  rotations, factor correlations and structure coefficients).
  `fit_indices_descriptive`, on the bootstrap path, pools every fit
  index the bootstrap replicates carry, so the structural constants
  among them (`df`, `df_null`) appear with a standard error of 0. The
  RMSEA confidence bounds are not among them: the replicate fits run
  without confidence intervals, so no per-replicate value exists to
  pool.

- mi_fit:

  On the `se = "sandwich"` (MI2S) path only: the single
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fit on the pooled correlation matrix \\\bar r\\ and pooled asymptotic
  covariance \\\tilde\Gamma\\. Its `orig_R` is \\\bar r\\ and its
  `Gamma` is \\\tilde\Gamma\\; the pooled `SE` and `CI` are taken from
  it, as are the pooled `fit_indices` (put into the common order above
  and extended with `pool_method`, while `mi_fit` keeps
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)'s
  own layout). `MI` is `NULL` on this path because the imputation
  uncertainty is carried by \\\tilde\Gamma\\ rather than by
  per-parameter Rubin pooling.

- mi_diagnostics:

  Diagnostics for the pooled model fit, `NULL` on the `se = "sandwich"`
  (MI2S) path, where there is one fit and no D2 pool. `m` is the number
  of imputations that entered the pool. `D2_F`, `D2_df1`, `D2_df2`,
  `D2_chi_asymptotic`, `ARIV` and `FMI` describe the D2 pool of the
  model chi-square (the average relative increase in variance and the
  fraction of missing information it implies), and `chi_bar_naive` is
  the plain mean of the per-imputation statistics for comparison; the
  `*_null` entries are the same quantities for the independence
  baseline. `D2_F` is the rule's raw statistic and is reported
  unfloored, so it is negative whenever the between-imputation
  variability of the component statistics exceeds the pooled discrepancy
  – a diagnostic of the pool rather than a fit statistic. The reported
  fit is not affected: the pooled chi-square is floored at zero and its
  p-value is 1 in that case. `chi_cfi` and `chi_null_cfi` are the pooled
  model and baseline **chi-squares** on the common \\N - 1\\
  noncentrality scale. `chi_cfi` is the statistic the reported RMSEA is
  formed from; the pair also gives a reference CFI formed the
  conventional way, `1 - (chi_cfi - df) / (chi_null_cfi - df_null)` (and
  analogously for TLI) – a different quantity from the reported CFI/TLI,
  which average the per-imputation indices.

- mi_admissibility:

  Admissibility and convergence of the component fits, kept on the
  pooled object so a saved solution carries the record independently of
  `fits`: `m` (the number of imputations that entered the pool),
  `heywood_imputations` (the indices of the fits with at least one
  Heywood case – an improper solution where a variable's communality is
  at or above 1, or its uniqueness is fixed at the estimation boundary),
  `n_heywood_items` (the number of flagged variables per imputation),
  `nonconverged` (the indices whose extraction reported a non-zero
  convergence code), and `iter` (the iterations each extraction used).
  Averaging aligned solutions pulls boundary communalities back inside
  the admissible range, so a pooled matrix with no Heywood case can
  still rest on component fits that had them;
  [`summary()`](https://rdrr.io/r/base/summary.html) reports the pooled
  count together with these.

- fits:

  The list of \\m\\ component
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fits, in the order of `data_list`, kept for per-imputation
  diagnostics. On the MI2S path these are the per-imputation fits whose
  inputs were pooled, not the pooled fit itself (which is `mi_fit`).

- alignment:

  Metadata from aligning the rotated solutions, `NULL` when no rotation
  was requested or on the MI2S path (one fit, one gauge). Under
  `target_method = "first_target"`: the `method` used, the `target` it
  aligned to, the per-imputation `target_rotations`, the indices of any
  `point_rotation_failures`, and whether every inner alignment
  `converged`. Under `target_method = "consensus"` it is the full
  [`efa_procrustes()`](https://mdsteiner.github.io/EFAtools/reference/efa_procrustes.md)-based
  GPA record: the converged `target`, the `aligned_loadings` and
  `aligned_phi`, the iteration `history`, convergence flags, and the
  multi-start summary.

- settings:

  The component fits'
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  settings with the pooling settings added: `pooled` (always `TRUE`),
  `pooled_N` and `N` (the mean N across imputations), `n_imputations`,
  `component_se` (the `se` the component fits used), `target_method`,
  `align_unrotated`, `fit_pool_method`, `p`, `ci` and `rmsea_ci_level`.
  `se` records what was actually pooled, so it is `"none"` when pooled
  standard errors could not be produced although the component fits
  computed them (`component_se` keeps the request).

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

Both routes assume the values are missing at random (MAR). Which one to
prefer is largely practical: FIML is a single, efficient fit and is the
simpler default when the analysis model is the whole story, whereas
multiple imputation is more flexible when the imputation model should
draw on auxiliary variables not in the factor model, or when the same
imputations feed several downstream analyses.

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
Lorenzo-Seva & Van Ginkel 2016), starting from the *medoid* imputation's
rotated solution – the one closest in aligned squared distance to all
the others. `"consensus"` is supported for orthogonal rotations only.

Factor loadings are only unique up to rotation; "gauge" here means which
particular rotation, or orientation, a solution is expressed in. The two
methods differ in the rotational gauge the pooled solution ends up in,
and so in how it responds to the order of `data_list`. The GPA iteration
moves its target toward the centroid but keeps the gauge of the solution
it started from, so starting it at the medoid – a property of the set,
not of the list order – makes the pooled rotated solution invariant to
that order, as the pooled unrotated solution already is.
`"first_target"` anchors on the first imputation by construction, so an
atypical first imputation fixes the orientation for every other one.
Permuting `data_list` therefore moves the pooled pattern – by a few
hundredths of a loading unit when imputations are similar, more when
they disagree. For oblique rotations, the factor correlations move with
it. Where the two anchors coincide – and, more generally, where the
imputations agree – the two methods give effectively the same pooled
estimate and `"consensus"` is simply the more expensive; where they do
not, the pooled patterns differ by the rotation between the two gauges.
Passing `start` through `consensus_args` overrides the medoid anchor and
makes the consensus order dependent again.

`align_unrotated` controls how unrotated loadings are aligned before
pooling: `"signed_tucker_congruence"` (the default) matches them up to
factor reordering and sign changes, `"procrustes"` aligns them to the
first imputation by orthogonal Procrustes rotation, and `"none"`
averages them as returned by
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).

The default anchors the matching on the *medoid* imputation (defined
above) rather than on whichever imputation happens to come first, so the
pooled *unrotated* solution does not depend on the order of `data_list`.
The rotated solution is aligned separately, against a reference chosen
by `target_method`, and still depends on that reference.

The pooled unrotated matrix is then returned in the same gauge as a
single
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
fit. Its identifying constraint differs by extraction method – a
principal-axis extraction and a maximum-likelihood extraction fix the
rotation differently (Anderson & Rubin 1956; Lawley & Maxwell 1971) –
and is read off each component fit automatically, so the pooled matrix
can be compared element-by-element with an
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
solution. A solution that meets neither constraint – an improper one,
say – is left as aligned. The correction is a common orthogonal
rotation, so communalities, the total variance accounted for, the
model-implied correlation matrix, the residuals, and RMSR are unchanged;
only the split of variance across factors moves. `"procrustes"` and
`"none"` keep their first-imputation anchor and are returned as aligned.

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

The model chi-square and the indices derived from it (ECVI and the
descriptive AIC/BIC) are pooled with the D2 rule (Li, Meng, Raghunathan
& Rubin, 1991), not arithmetically averaged. RMSEA is pooled by the same
rule but from a second D2 pool of the per-imputation discrepancies taken
on the uncorrected \\N - 1\\ scale. The printed RMSEA therefore does not
reconcile by hand with the printed chi-square; the statistic it is
formed from is `chi_cfi` in `mi_diagnostics`. Because D2 shrinks the
pooled chi-square in proportion to the between-imputation variability,
the pooled RMSEA can fall below the mean of the per-imputation RMSEAs;
read it together with the per-imputation fit. The incremental indices
CFI (Bentler, 1990) and TLI (Tucker & Lewis, 1973) are instead the
average of the per-imputation indices, which keeps them consistent with
the component fits and avoids out-of-range values; the separately pooled
model and baseline chi-squares those indices would be formed from remain
available in `mi_diagnostics`. AIC and BIC, if returned, are
chi-square-derived descriptive quantities and are not likelihood-based
MI information criteria. They are reported only where the component fits
report them: whenever a component withholds them – any
`cor_method = "fiml"` fit, and any fit whose chi-square is a scaled
statistic, such as `se = "sandwich"` – the pooled AIC, BIC, and ECVI are
`NA` too, matching what
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
returns for a single such fit. On the sandwich/MI2S route the chi-square
is the single fit's scaled statistic rather than a D2 pool.

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
the Barnard-Rubin (1999) adjustment reduces to this form). NA
propagation is fail-closed: if any imputation is NA at an element, all
pooled outputs for that element are NA. When a rotation was requested,
the rotated loadings, communalities, and (for oblique rotations) factor
correlations and structure coefficients are pooled as well; residual SE
pooling is available only on the bootstrap path. Under
`align_unrotated = "procrustes"` the full unrotated covariance
`vcov_unrot_loadings` (populated by `se = "information"`) is propagated
through the alignment, so it must be present and reliable on every fit.
The default alignment also mixes loading columns, through the common
canonical-gauge rotation, and so propagates the same covariance; where a
fit does not carry it, the unrotated standard errors are returned as
`NA` rather than aborting, and the remaining families still pool.

A rotated-loading standard error is conditional on the rotation
criterion used (see References). For both orthogonal and oblique
rotations the within-imputation variance is therefore each fit's own
criterion-aware delta-method rotated SE (the quantity
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
imputations, and a single `EFA` model is fitted to the pooled
correlation with the pooled covariance, \\\tilde\Gamma\\, as the robust
meat (its diagonal as the weights for `estimator = "DWLS"`). Because
there is only one fit and one rotational gauge, this route bypasses the
per-imputation alignment: `target_method` and `align_unrotated` do not
apply. The fitted object carries native scaled-shifted chi-square
statistics and sandwich SEs that already reflect the multiple-imputation
uncertainty, so the chi-square is not D2-pooled and the
likelihood-ratio-based AIC/BIC/ECVI are `NA`; it is returned in the
`mi_fit` slot, with the per-imputation `fits` retained for diagnostics.
The pooled fit uses the same
[`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
and
[`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
tuning (including any rotation-engine extras) as the per-imputation
fits. At least 20 imputations are recommended for the scaled-shifted
statistic, and more (around 100) at higher rates of missingness
(Sriutaisuk et al. 2025). The polychoric/tetrachoric (ordinal) case is
the primary, best-evaluated target; the continuous-Pearson case uses the
same recipe but is less benchmarked.

## Conditions

Errors and warnings raised by `efa_mi()` are classed, with an
`efa_pooled_` prefix (`efa_consensus_` for the consensus target) –
except the dots validation shared with
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
which signals `efa_flat_knob_in_dots` or `efa_renamed_arg` – so they can
be caught programmatically. The message shown explains what went wrong
and, where relevant, how to fix it.

## References

Anderson, T. W., & Rubin, H. (1956). Statistical inference in factor
analysis. In *Proceedings of the Third Berkeley Symposium on
Mathematical Statistics and Probability* (Vol. 5, pp. 111-150).
University of California Press.

Archer, C. O., & Jennrich, R. I. (1973). Standard errors for rotated
factor loadings. *Psychometrika*, 38(4), 581-592.

Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
multiple imputation. *Biometrika*, 86(4), 948-955.

Bentler, P. M. (1990). Comparative fit indexes in structural models.
*Psychological Bulletin*, 107(2), 238-246.

Chung, S., & Cai, L. (2019). Alternative multiple imputation inference
for categorical structural equation modeling. *Multivariate Behavioral
Research*, 54(3), 323-337.

Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika*,
40(1), 33-51.

Jennrich, R. I. (1973). Standard errors for obliquely rotated factor
loadings. *Psychometrika*, 38(4), 593-604.

Jennrich, R. I. (1974). Simplified formulae for standard errors in
maximum-likelihood factor analysis. *British Journal of Mathematical and
Statistical Psychology*, 27(1), 122-131.

Lawley, D. N., & Maxwell, A. E. (1971). *Factor analysis as a
statistical method* (2nd ed.). Butterworths.

Li, K. H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
Significance levels from repeated p-values with multiply-imputed data.
*Statistica Sinica*, 1(1), 65-92.

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

Zhang, G., Preacher, K. J., & Jennrich, R. I. (2012). The infinitesimal
jackknife with exploratory factor analysis. *Psychometrika*, 77(4),
634-648.

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
mod <- efa_mi(dat_list, n_factors = 1, estimator = "ML")
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
mod
#> 
#> Pooled EFA across 3 imputations performed with estimator = 'ML' and
#>   rotation = 'none'.
#> Pooling settings: align_unrotated = 'signed_tucker_congruence',
#>   fit_pool_method = 'D2'.
#> 
#> ── Unrotated Loadings ──────────────────────────────────────────────────────────
#> 
#>             F1    h2    u2
#> fun        .808  .652  .348
#> friends    .837  .701  .299
#> enjoy      .871  .759  .241
#> hurt       .757  .573  .427
#> part       .814  .662  .338
#> commonly   .818  .669  .331
#> chances    .783  .613  .387
#> attracted  .844  .713  .287
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                 F1
#> SS loadings   5.343
#> Prop Tot Var   .668
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> D2-pooled χ²(20) = 0.00, p = 1.000
#> CFI (avg. over imputations): .99
#> TLI (avg. over imputations): .98
#> RMSEA [90% CI]: .00 [.00; .00]
#> AIC: -40.00
#> BIC: -133.94
#> ECVI: 0.04
#> CAF: .50
#> SRMR: .02
#> Note: the pooled χ² is the D2 statistic; its p uses the D2 reference F(20, .1),
#> not the χ²(20) tail.
#> Note: CFI and TLI are averaged over the imputations, not formed from the
#> separately pooled model and baseline statistics in `mi_diagnostics`.

# \donttest{
# add computation of standard errors and CIs
mod <- efa_mi(dat_list, n_factors = 1, estimator = "ML", se = "np-boot")
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
mod
#> 
#> Pooled EFA across 3 imputations performed with estimator = 'ML' and
#>   rotation = 'none'.
#> Pooling settings: align_unrotated = 'signed_tucker_congruence',
#>   fit_pool_method = 'D2'.
#> 
#> ── Unrotated Loadings ──────────────────────────────────────────────────────────
#> 
#>             F1    h2    u2
#> fun        .808  .652  .348
#> friends    .837  .701  .299
#> enjoy      .871  .759  .241
#> hurt       .757  .573  .427
#> part       .814  .662  .338
#> commonly   .818  .669  .331
#> chances    .783  .613  .387
#> attracted  .844  .713  .287
#> 
#> Legend:
#>   bold = |loading| >= .300
#>   grey = below cutoff
#>   red h2/u2 = Heywood-relevant value
#> 
#> ── Variances Accounted for ─────────────────────────────────────────────────────
#> 
#>                 F1
#> SS loadings   5.343
#> Prop Tot Var   .668
#> 
#> ── Model Fit ───────────────────────────────────────────────────────────────────
#> 
#> D2-pooled χ²(20) = 0.00, p = 1.000
#> CFI (avg. over imputations) [95% bootstrap/MI-CI]: .99 [.96, 1.02]
#> TLI (avg. over imputations) [95% bootstrap/MI-CI]: .98 [.94, 1.02]
#> RMSEA [90% CI] [95% bootstrap/MI-CI]: .00 [.00; .00] [-.01, .14]
#> AIC [95% bootstrap/MI-CI]: -40.00 [-93.87, 197.52]
#> BIC [95% bootstrap/MI-CI]: -133.94 [-187.82, 103.58]
#> ECVI [95% bootstrap/MI-CI]: 0.04 [-.03, .33]
#> CAF [95% bootstrap/MI-CI]: .50 [.48, .52]
#> SRMR [95% bootstrap/MI-CI]: .02 [.00, .03]
#> Note: the pooled χ² is the D2 statistic; its p uses the D2 reference F(20, .1),
#> not the χ²(20) tail.
#> Note: CFI and TLI are averaged over the imputations, not formed from the
#> separately pooled model and baseline statistics in `mi_diagnostics`.
#> 
#> Note: Bootstrap/MI CIs based on 1000 bootstrap samples per imputation.
# }
```
