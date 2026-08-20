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
  alpha = 0.3,
  ...,
  max_iter_CD = 50,
  n_fac_theor = NA,
  estimator = c("ML", "PAF", "ULS"),
  gof = c("CAF", "CFI", "RMSEA"),
  eigen_type_HULL = c("SMC", "PCA", "EFA"),
  eigen_type_other = c("SMC"),
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

## Arguments

- x:

  data.frame or matrix. Raw data, or a correlation matrix. If `"CD"` is
  included as a criterion, x must be raw data.

- criteria:

  character. Which factor retention methods to run: one or more of
  `"CD"`, `"EKC"`, `"HULL"`, `"KGC"`, `"MAP"`, `"NEST"`, `"PARALLEL"`,
  `"SCREE"`, and `"SMT"` (see details). The default runs a subset of
  commonly used, well-performing methods, listed in the details.

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
  given as input. Default is `"pairwise.complete.obs"`.

- cor_method:

  character. Correlation computed from raw data: `"pearson"`,
  `"spearman"`, or `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)), or `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations (a two-step
  estimator). `CD`, `PARALLEL`, `NEST`, `HULL`, and `SMT` do not support
  `"poly"` / `"tetra"` and are skipped automatically if you request them
  together. Default is `"pearson"`.

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

  numeric. Passed to
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md).
  The alpha level used to test the significance of the improvement added
  by an additional factor. Default is .30.

- ...:

  Further arguments passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  in
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
  [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
  [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
  [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
  and (through its parallel analysis and its own candidate fits)
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md).
  An argument that
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  does not recognize causes an error. The estimation tuning knobs are
  not passed here; they live in `estimate_control`. The standard-error
  arguments (`se`, `b_boot`, `ci`, `seed`) are not accepted, because the
  criterion fits are internal steps whose standard errors are not
  reported. Arguments listed after `...` must be given by their full
  name (R matches an abbreviated name only against the arguments before
  `...`), so a tuning knob such as `max_iter` cannot be mistaken for
  `max_iter_CD`.

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

  character. Passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  in
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
  [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
  [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
  and
  [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md).
  The estimator to use. One of `"PAF"`, `"ULS"`, or `"ML"`, for
  principal axis factoring, unweighted least squares, and maximum
  likelihood, respectively. The default here is `"ML"`. Some criteria
  default to something else when called on their own (for example,
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)
  defaults to `"PAF"`), so results from `efa_retain()` can differ from
  calling that criterion directly unless you set `estimator` to match.
  In
  [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
  [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
  and
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  it only takes effect when the respective `eigen_type` includes
  `"EFA"`.

- gof:

  character. Passed to
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md).
  The goodness of fit index to use. Either `"CAF"`, `"CFI"`, or
  `"RMSEA"`, or any combination of them. With the `"PAF"` estimator,
  only the CAF can be used as goodness of fit index. For details on the
  CAF, see Lorenzo-Seva, Timmerman, and Kiers (2011).

- eigen_type_HULL:

  character. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  in
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md).
  What the eigenvalues in the parallel analysis are based on. One of
  `"SMC"`, `"PCA"`, or `"EFA"` – different ways of estimating how much
  variance each indicator shares with the others before the eigenvalues
  are computed. `"SMC"` (default) uses each indicator's squared multiple
  correlation with the others (its diagonal value in the correlation
  matrix). `"PCA"` leaves the diagonal at 1, so each indicator's total
  variance – not just the shared part – feeds into the eigenvalues.
  `"EFA"` uses the communalities from a fitted EFA solution instead.

- eigen_type_other:

  character. Passed to
  [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
  [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
  and
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md).
  The same as eigen_type_HULL, but multiple inputs are possible here
  (any combination of `"PCA"`, `"SMC"`, and `"EFA"`). Default is
  `"SMC"`.

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
  The number of datasets to simulate. Default is 1000.

- percent:

  numeric. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  (also within
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)).
  The percentile to take from the simulated eigenvalues. Default is 95.

- decision_rule:

  character. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  (also within
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)).
  Which rule to use to determine the number of factors to retain.
  Default is `"means"`, which uses the average simulated eigenvalues.
  `"percentile"` uses the percentiles specified in `percent`.
  `"crawford"` uses the 95th percentile for the first factor and the
  mean afterwards (based on Crawford et al., 2010).

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
  The number of datasets to simulate. Default is 1000.

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
  object with the estimation settings for the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fits run by the criteria that fit a model
  ([`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
  [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
  [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
  [`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
  and
  [`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)).
  `NULL` (default) uses the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  defaults. It only applies to criteria that fit a model, and only to
  the parts of that fit each criterion actually runs.
  [`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md),
  [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md),
  and
  [`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md)
  fit no model.
  [`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
  [`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
  and
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  only fit one when their `eigen_type` includes `"EFA"`.
  [`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)
  fits with maximum likelihood by definition, so only `start_method`
  takes effect there. All fits are unrotated, so no rotation settings
  apply.

## Value

A list of class `c("efa_retain", "N_FACTORS")`, the trailing class
keeping `inherits(x, "N_FACTORS")` working for code written against the
superseded name. It contains

- suitability:

  A list with the results from
  [`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)
  and
  [`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md)
  (`bartlett` and `kmo`), or `NULL` if `suitability = FALSE`.

- outputs:

  A named list with one `efa_retention` object per factor retention
  criterion that was run (see, e.g.,
  [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md)).

- n_factors:

  A named numeric vector with the suggested number of factors per
  criterion and, where a criterion has several variants, per variant
  (e.g. `EKC_BvA2017` or `PARALLEL_SMC`). Criteria without a numeric
  suggestion (the scree plot) are not included. The printed summary's
  "most common" value is based on each criterion's own most frequent
  (modal) suggestion among its variants, not a plain tally of this
  vector, so counting entries here by hand can give a different answer.

- not_run:

  A named character vector with the criteria that were skipped or failed
  and the reason, or `NULL` if all requested criteria ran.

- settings:

  A list of the settings used. Its `criteria` element records the
  requested criteria, in the order they were given, while `outputs` and
  `n_factors` are in the order in which the criteria were run. `gof`
  records the requested Hull goodness-of-fit indices, and `gof_used`
  records the ones the Hull method actually computed (it reduces them to
  `"CAF"` for the PAF estimator). `gof_used` is `NA` when HULL was not
  requested, was skipped, or failed.

## Details

By default, the entered data are checked for suitability for factor
analysis using the following methods (see the respective documentation
for details):

- Bartlett's test of sphericity (see
  [`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md))

- Kaiser-Meyer-Olkin criterion (see
  [`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md))

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
the range of their suggestions and the most common one. Auerswald and
Moshagen (2019) compare the criteria and give guidance on the selection.

The comparison data, parallel analysis, and NEST criteria compare the
data against simulated reference data, so their suggested numbers of
factors vary slightly from run to run. The Hull method also varies,
because it calls
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
to set its upper bound. Call
[`base::set.seed()`](https://rdrr.io/r/base/Random.html) before
`efa_retain()` to make the results reproducible.

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
# Default criteria, with correlation matrix and estimator "ML" (where needed)
# This will throw a warning for CD, as no raw data were specified
# The simulation-based criteria are seeded to make the run reproducible
set.seed(42)
nfac_all <- efa_retain(test_models$baseline$cormat, N = 500, estimator = "ML",
                       n_datasets = 100, n_datasets_nest = 100)
#> Warning: `x` is a correlation matrix, but "CD" needs raw data.
#> ℹ Skipping "CD".

# The same as above, but without "CD"
nfac_wo_CD <- efa_retain(test_models$baseline$cormat, criteria = c("EKC",
                         "HULL", "PARALLEL", "NEST"), N = 500,
                         estimator = "ML", n_datasets = 100,
                         n_datasets_nest = 100)

# Use PAF instead of ML (this will take longer). PAF only supports "CAF" as
# gof for the Hull method, so set it explicitly to avoid the automatic message.
nfac_PAF <- efa_retain(test_models$baseline$cormat, criteria = c("EKC",
                       "HULL", "PARALLEL", "NEST"), N = 500,
                       estimator = "PAF", gof = "CAF", n_datasets = 100,
                       n_datasets_nest = 100)

# Back to the default ML estimator (unlike above), with only "PCA" type eigenvalues
nfac_PCA <- efa_retain(test_models$baseline$cormat, criteria = c("EKC",
                       "HULL", "PARALLEL", "NEST"), N = 500,
                       estimator = "ML", eigen_type_other = "PCA",
                       n_datasets = 100, n_datasets_nest = 100)

# Use raw data, such that CD can also be performed
nfac_raw <- efa_retain(GRiPS_raw, estimator = "ML", N_pop = 500,
                        N_samples = 20, n_datasets = 100,
                        n_datasets_nest = 100)
#> Warning: The suggested maximum number of factors was 2, but the Hull method needs at
#> least 3.
#> ℹ Setting it to 3.
# }
```
