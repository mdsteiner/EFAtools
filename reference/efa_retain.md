# Various factor retention criteria

Choosing the number of factors to retain is one of the most important
decisions in an exploratory factor analysis (EFA). Many criteria exist
to help with this choice. This function runs several of them together,
and can also check whether the data are suitable for factor analysis.

## Usage

``` r
efa_retain(
  x,
  criteria = c("CD", "EKC", "HULL", "MAP", "NEST", "PARALLEL"),
  suitability = TRUE,
  N = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  n_factors_max = NA,
  N_pop = 10000,
  N_samples = 500,
  alpha = lifecycle::deprecated(),
  ...,
  alpha_CD = 0.3,
  max_iter_CD = 50,
  n_fac_theor = NA,
  estimator = c("ML", "PAF", "ULS"),
  gof = "CAF",
  eigen_type_HULL = c("PCA", "SMC", "EFA"),
  decision_rule_HULL = c("percentile", "means", "crawford"),
  eigen_type_other = "PCA",
  n_factors = 1,
  n_datasets = 1000,
  percent = 95,
  decision_rule = c("means", "percentile", "crawford"),
  ekc_type = lifecycle::deprecated(),
  n_datasets_nest = 1000,
  alpha_nest = 0.05,
  show_progress = FALSE,
  estimate_control = NULL
)
```

## Source

Auerswald, M., & Moshagen, M. (2019). How to determine the number of
factors to retain in exploratory factor analysis: A comparison of
extraction methods under realistic conditions. Psychological Methods,
24(4), 468–491. https://doi.org/10.1037/met0000200

Lim, S., & Jahng, S. (2019). Determining the number of factors using
parallel analysis and its recent variants. Psychological Methods, 24(4),
452–467. https://doi.org/10.1037/met0000230

Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011). The Hull
method for selecting the number of common factors. Multivariate
Behavioral Research, 46(2), 340-364.
https://doi.org/10.1080/00273171.2011.564527

Ruscio, J., & Roche, B. (2012). Determining the number of factors to
retain in an exploratory factor analysis using comparison data of known
factorial structure. Psychological Assessment, 24, 282–292.
https://doi.org/10.1037/a0025697

## Arguments

- x:

  data.frame or matrix. Raw data, or a correlation matrix. `"CD"` needs
  raw data.

- criteria:

  character. The factor retention criteria to run: one or more of
  `"CD"`, `"EKC"`, `"HULL"`, `"KGC"`, `"MAP"`, `"NEST"`, `"PARALLEL"`,
  `"SCREE"`, and `"SMT"`. The default runs a subset of commonly used,
  well-performing criteria (see details).

- suitability:

  logical. Whether the data should be checked for suitability for factor
  analysis using Bartlett's test of sphericity and the
  Kaiser-Meyer-Olkin criterion (see details). Default is `TRUE`.

- N:

  numeric. The number of observations. Only needed if x is a correlation
  matrix.

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is `"pairwise.complete.obs"`. `"CD"` always
  uses only the complete cases: it removes the incomplete rows with a
  warning.

- cor_method:

  character. Correlation computed from raw data: `"pearson"`,
  `"spearman"`, or `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)), or `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations (a two-step
  estimator). With `"poly"` / `"tetra"`, `"HULL"`, `"NEST"`,
  `"PARALLEL"`, and `"SMT"` are skipped, and `"CD"` and Bartlett's test
  use the Pearson correlations of the raw data (see details). Default is
  `"pearson"`.

- n_factors_max:

  numeric. Passed to
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md).
  The maximum number of factors to test against. Larger numbers will
  increase the duration the procedure takes, but test more possible
  solutions. If left NA (default), the maximum number of factors for
  which the model is still over-identified (df \> 0) is used.

- N_pop:

  numeric. Passed to
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md).
  Size of finite populations of comparison data. Default is 10000.

- N_samples:

  numeric. Passed to
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md).
  Number of samples drawn from each population. Default is 500.

- alpha:

  **\[deprecated\]** Use `alpha_CD` instead.

- ...:

  Further arguments passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  in the criteria that fit a model, except
  [`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md).
  Unknown arguments, rotation settings, estimation settings (use
  `estimate_control`), and the standard-error arguments (`se`, `b_boot`,
  `ci`, `seed`) cause an error. Arguments after `...` must be given by
  their full name.

- alpha_CD:

  numeric. Passed to
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md)
  as `alpha`. The alpha level used to test the significance of the
  improvement added by an additional factor. Default is .30.

- max_iter_CD:

  numeric. Passed to
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md).
  The maximum number of iterations to perform after which the iterative
  PAF procedure is halted. Default is 50.

- n_fac_theor:

  numeric. Passed to
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md).
  Theoretical number of factors to retain. The Hull method uses one plus
  the larger of this number and the number of factors suggested by
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  as its upper bound.

- estimator:

  character. The estimator for the criteria that fit a model: `"ML"`
  (default), `"PAF"`, or `"ULS"`. Some criteria have a different default
  when called on their own (for example,
  [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md)
  uses `"PAF"`).
  [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
  [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
  and
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  use it only for `"EFA"` eigenvalues.
  [`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)
  always uses maximum likelihood.

- gof:

  character. Passed to
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md).
  The goodness of fit index to use: `"CAF"` (default), `"CFI"`, or
  `"RMSEA"`, or any combination of them. With the `"PAF"` estimator, the
  Hull method uses only the CAF. In Lorenzo-Seva, Timmerman, and Kiers
  (2011), with factor correlations up to .40, the CFI was the most
  accurate index. With strongly correlated factors (e.g., .70), CFI and
  RMSEA can make the Hull method retain too few factors. For details on
  the CAF, see Lorenzo-Seva et al. (2011).

- eigen_type_HULL:

  character. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  in
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md).
  The eigenvalues of the parallel analysis: `"PCA"` (default) uses the
  correlation matrix, `"SMC"` puts the squared multiple correlations on
  its diagonal, and `"EFA"` puts the communalities of an EFA solution on
  its diagonal (see details).

- decision_rule_HULL:

  character. Passed to
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)
  as `decision_rule`. The rule of the parallel analysis that sets the
  upper bound of the Hull method, with the same choices as
  `decision_rule`. Default is `"percentile"` (see details).

- eigen_type_other:

  character. Passed to
  [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
  [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
  and
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md).
  The same types as `eigen_type_HULL`, but one or more of them. Default
  is `"PCA"`.

- n_factors:

  numeric. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  (also within
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)),
  [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
  and
  [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md).
  Number of factors to extract if `"EFA"` is included in
  `eigen_type_HULL` or `eigen_type_other`. Default is 1.

- n_datasets:

  numeric. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  (also within
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)).
  The number of reference correlation matrices to simulate. Default is
  1000.

- percent:

  numeric. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  (also within
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)).
  The percentile to take from the simulated eigenvalues. Default is 95.

- decision_rule:

  character. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  for the `"PARALLEL"` criterion only (the Hull method uses
  `decision_rule_HULL`). Which rule to use to determine the number of
  factors to retain. Default is `"means"`, which uses the average
  simulated eigenvalues. `"percentile"` uses the percentiles specified
  in `percent`. `"crawford"` uses the 95th percentile for the first
  factor and the mean afterwards (based on Crawford et al., 2010).

- ekc_type:

  **\[deprecated\]** Accepted and ignored. It used to select between two
  ways to compute the
  [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md)
  reference values. The `"AM2019"` reference values do not depend on the
  observed eigenvalues. They therefore skip the empirical correction
  that defines the criterion, so they are no longer computed.

- n_datasets_nest:

  numeric. Passed to
  [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md).
  The number of reference correlation matrices to simulate. Default is
  1000.

- alpha_nest:

  numeric. Passed to
  [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md).
  The alpha level to use. The reference values are the eigenvalues at
  the (1 - alpha_nest) percentile. Default is .05.

- show_progress:

  logical. Whether a progress bar should be shown in the console.
  Default is FALSE.

- estimate_control:

  an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object with the estimation settings for the unrotated
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fits of the criteria that fit a model (all except
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md),
  [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md),
  and
  [`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md)).
  `NULL` (default) uses the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  defaults.
  [`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)
  always uses maximum likelihood, so only `start_method` applies there.

## Value

A list of class `c("efa_retain", "N_FACTORS")`, the trailing class
keeping `inherits(x, "N_FACTORS")` working for code written against the
superseded name. It contains

- suitability:

  A list with the results from
  [`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)
  and
  [`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md)
  (`bartlett` and `kmo`), or `NULL` if `suitability = FALSE`. `bartlett`
  is `NULL` if `N` is missing.

- outputs:

  A named list with one `efa_retention` object per factor retention
  criterion that was run, in alphabetical order (see, e.g.,
  [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md)).

- n_factors:

  A named numeric vector with the suggested number of factors per
  criterion and variant (e.g. `EKC_BvA2017` or `PARALLEL_PCA`), without
  the scree plot.

- not_run:

  A named character vector with the criteria that were skipped or failed
  and the reason, or `NULL` if all requested criteria ran.

- settings:

  A list of the settings used. `gof_used` records the Hull goodness of
  fit indices that were computed (`NA` if the Hull method did not run).

## Details

By default, `efa_retain()` checks if the data are suitable for factor
analysis. It uses Bartlett's test of sphericity
([`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)),
which needs `N`, and the Kaiser-Meyer-Olkin criterion
([`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md)).
A significant Bartlett's test only shows that the correlation matrix is
not an identity matrix. This is almost always true for realistic sample
sizes. With `cor_method = "poly"` or `"tetra"` and raw data, Bartlett's
test uses the Pearson correlations of the raw data. Its chi-square
distribution is not valid for polychoric correlations, because these
have a larger sampling variance. Zero polychoric correlations give zero
Pearson correlations, so the test has the same null hypothesis. The KMO
uses the analysed correlations.

The available factor retention criteria are the following (see the
respective documentation for details):

- Comparison data (see
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md))

- Empirical Kaiser criterion (see
  [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md))

- Hull method (see
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md))

- Kaiser-Guttman criterion (see
  [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md))

- Velicer's minimum average partial, MAP (see
  [`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md))

- Next Eigenvalue Sufficiency Test, NEST (see
  [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md))

- Parallel analysis (see
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md))

- Scree plot (see
  [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md))

- Sequential chi-square model tests, RMSEA lower bound, and AIC (see
  [`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md))

The default `criteria` are comparison data, the empirical Kaiser
criterion, the Hull method, MAP, NEST, and parallel analysis. No single
criterion is the most accurate in all conditions. `efa_retain()`
therefore runs several criteria together, and the printed summary gives
the range of their suggestions and the most common one (one vote per
criterion). Auerswald and Moshagen (2019) compare the criteria and give
guidance on the selection.

Parallel analysis, also in the Hull method, uses PCA eigenvalues by
default. With the eigenvalues of the reduced correlation matrix
(`"SMC"`), parallel analysis often retains too many factors,
particularly when the model does not fit the data exactly. In Lim and
Jahng (2019), parallel analysis with PCA eigenvalues was the most
accurate in most conditions. The `"PARALLEL"` criterion compares the
observed eigenvalues with the mean simulated eigenvalues by default. To
set the upper bound of the Hull method, parallel analysis uses the 95th
percentile instead, as Lorenzo-Seva et al. (2011) suggest. This helps to
prevent overfactoring (see
[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)).
The Hull method uses the CAF by default (see `gof`).

A criterion that cannot use the input is skipped: `"CD"` with a
correlation matrix, and `"HULL"`, `"NEST"`, `"PARALLEL"`, and `"SMT"`
with `cor_method = "poly"` or `"tetra"`. If you do not set `criteria`,
this occurs without a warning. With a correlation matrix and no `N`, the
criteria that need `N` (`"EKC"`, `"HULL"`, `"NEST"`, `"PARALLEL"`, and
`"SMT"`) are skipped with a warning. The printed result shows all
criteria that did not run. With `cor_method = "poly"` or `"tetra"`,
comparison data use the Pearson correlations of the raw data (see
[`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md)).

Comparison data, parallel analysis, NEST, and the Hull method use
simulated reference data, so their suggestions can change slightly from
run to run. The Hull method uses simulation only for the parallel
analysis that sets its upper bound. If parallel analysis also runs and
`eigen_type_other` includes `eigen_type_HULL`, the Hull method uses its
simulated eigenvalues with `decision_rule_HULL`, and does not simulate
again. To make the results reproducible, call
[`base::set.seed()`](https://rdrr.io/r/base/Random.html) first.

## See also

[`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md)
for data screening before retention, and
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
to extract the chosen number of factors.

Other factor retention criteria:
[`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md),
[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md),
[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
[`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
[`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md),
[`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
[`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)

## Examples

``` r
# \donttest{
# The simulation-based criteria are seeded to make the runs reproducible
set.seed(42)

# Correlation matrix: the default criteria except "CD", which needs raw data
efa_retain(test_models$baseline$cormat, N = 500, n_datasets = 100,
           n_datasets_nest = 100)
#> ── Sampling adequacy and sphericity ────────────────────────────────────────────
#> 
#> ✔ Overall KMO = 0.916: marvellous (Kaiser & Rice, 1974).
#> ✔ Bartlett's test of sphericity: χ²(153) = 2173.28, p < .001 (significant).
#> 
#> ── Suggested number of factors ─────────────────────────────────────────────────
#> 
#> 6 suggestions from 5 criteria, ranging from 1 to 3 factors (most common: 3).
#> 
#> Empirical Kaiser Criterion
#> • Braeken & van Assen (2017): 3
#> 
#> Hull method
#> Estimator: ML
#> • CAF: 3
#> 
#> Minimum average partial
#> • Original implementation (TR2): 1
#> • Revised implementation (TR4): 3
#> 
#> Next Eigenvalue Sufficiency Test
#> • Suggested number of factors: 3
#> 
#> Parallel analysis
#> Eigenvalues found using PCA; 100 simulated datasets.
#> • PCA: 3
#> 
#> ── Criteria not run ────────────────────────────────────────────────────────────
#> 
#> • CD: needs raw data, but a correlation matrix was supplied

# Raw data: "CD" runs too
efa_retain(GRiPS_raw, N_pop = 500, N_samples = 20, n_datasets = 100,
           n_datasets_nest = 100)
#> ── Sampling adequacy and sphericity ────────────────────────────────────────────
#> 
#> ✔ Overall KMO = 0.955: marvellous (Kaiser & Rice, 1974).
#> ✔ Bartlett's test of sphericity: χ²(28) = 5054.06, p < .001 (significant).
#> 
#> ── Suggested number of factors ─────────────────────────────────────────────────
#> 
#> 7 suggestions from 6 criteria, all suggesting 1 factor.
#> 
#> Comparison data
#> • Suggested number of factors: 1
#> 
#> Empirical Kaiser Criterion
#> • Braeken & van Assen (2017): 1
#> 
#> Hull method
#> Estimator: ML
#> • CAF: 1
#> 
#> Minimum average partial
#> • Original implementation (TR2): 1
#> • Revised implementation (TR4): 1
#> 
#> Next Eigenvalue Sufficiency Test
#> • Suggested number of factors: 1
#> 
#> Parallel analysis
#> Eigenvalues found using PCA; 100 simulated datasets.
#> • PCA: 1
# }
```
