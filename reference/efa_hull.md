# Hull method for determining the number of factors to retain

Implementation of the Hull method suggested by Lorenzo-Seva, Timmerman,
and Kiers (2011), with an extension to principal axis factoring. See
details for parallelization.

## Usage

``` r
efa_hull(
  x,
  N = NA,
  n_fac_theor = NA,
  estimator = c("PAF", "ULS", "ML"),
  gof = c("CAF", "CFI", "RMSEA"),
  eigen_type = c("SMC", "PCA", "EFA"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  n_datasets = 1000,
  percent = 95,
  decision_rule = c("means", "percentile", "crawford"),
  n_factors = 1,
  estimate_control = NULL,
  ...
)
```

## Source

Auerswald, M., & Moshagen, M. (2019). How to determine the number of
factors to retain in exploratory factor analysis: A comparison of
extraction methods under realistic conditions. Psychological Methods,
24(4), 468–491. https://doi.org/10.1037/met0000200

Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011). The Hull
method for selecting the number of common factors. Multivariate
Behavioral Research, 46(2), 340-364.

## Arguments

- x:

  matrix or data.frame. Dataframe or matrix of raw data or matrix with
  correlations.

- N:

  numeric. Number of cases in the data. This is passed to
  [efa_parallel](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md).
  Only has to be specified if x is a correlation matrix, otherwise it is
  determined based on the dimensions of x.

- n_fac_theor:

  numeric. Theoretical number of factors to retain. One plus the larger
  of this number and the number of factors suggested by
  [efa_parallel](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  is used as the upper bound *J* of factors to extract in the Hull
  method.

- estimator:

  character. The estimator to use. One of `"PAF"`, `"ULS"`, or `"ML"`,
  for principal axis factoring, unweighted least squares, and maximum
  likelihood, respectively.

- gof:

  character. The goodness of fit index to use. Either `"CAF"`, `"CFI"`,
  or `"RMSEA"`, or any combination of them. With the `"PAF"` estimator,
  only the CAF can be used as goodness of fit index. For details on the
  CAF, see Lorenzo-Seva, Timmerman, and Kiers (2011).

- eigen_type:

  character. On what the eigenvalues should be found in the parallel
  analysis. Can be one of `"SMC"`, `"PCA"`, or `"EFA"`. If using `"SMC"`
  (default), the diagonal of the correlation matrices is replaced by the
  squared multiple correlations (SMCs) of the indicators. If using
  `"PCA"`, the diagonal values of the correlation matrices are left to
  be 1. If using `"EFA"`, eigenvalues are found on the correlation
  matrices with the final communalities of an EFA solution as diagonal.
  This is passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md).

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is `"pairwise.complete.obs"`.

- cor_method:

  character. One of `"pearson"`, `"spearman"`, or `"kendall"`, passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). `"poly"` and
  `"tetra"` are not supported because `HULL` derives its factor-search
  bound from an internal parallel analysis against continuous reference
  data. Default is `"pearson"`.

- n_datasets:

  numeric. The number of datasets to simulate. Default is 1000. This is
  passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md).

- percent:

  numeric. The percentile to take from the simulated eigenvalues.
  Default is 95. This is passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md).

- decision_rule:

  character. Which rule to use to determine the number of factors to
  retain. Default is `"means"`, which will use the average simulated
  eigenvalues. `"percentile"`, uses the percentiles specified in
  percent. `"crawford"` uses the 95th percentile for the first factor
  and the mean afterwards (based on Crawford et al, 2010). This is
  passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md).

- n_factors:

  numeric. Number of factors to extract if `"EFA"` is included in
  `eigen_type`. Default is 1. This is passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md).

- estimate_control:

  an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object with the estimation settings for the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fits of the 0 to *J* factor solutions, and for the fit inside the
  internal
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  call. `NULL` (default) uses the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  defaults. This object carries estimation settings only; the fits are
  always unrotated, which the hull statistics (CFI, RMSEA, CAF) do not
  depend on.

- ...:

  Further arguments passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
  also in
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md).
  The estimation tuning knobs are not passed here; they live in
  `estimate_control`, a rotation setting is not accepted because the
  fits are unrotated, and neither are the standard-error arguments
  (`se`, `b_boot`, `ci`, `seed`), because the fits are internal steps
  whose standard errors are not reported.

## Value

An object of class `efa_retention` (see
[`print.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_retention.md)
and
[`plot.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_retention.md)
for the print and plot methods). Its main fields are:

- n_factors:

  A named numeric vector with the suggested number of factors for each
  requested goodness-of-fit index (`"CAF"`, `"CFI"`, and/or `"RMSEA"`).

- results:

  A list with one record per goodness-of-fit index, each holding the
  goodness-of-fit values, the degrees of freedom, the hull membership,
  and the retained solution used for printing and plotting.

- settings:

  A list of the settings used, including `n_fac_max`, the upper bound
  *J* of the number of factors to extract (see details).

## Details

The Hull method aims to find a model with an optimal balance between
model fit and number of parameters, retaining only major factors
(Lorenzo-Seva, Timmerman, & Kiers, 2011). It fits 0 to *J* factors –
where *J* is the number of factors suggested by parallel analysis (or
`n_fac_theor`, if that is larger), plus one – keeps the solutions on the
upper boundary of the convex hull of goodness-of-fit against degrees of
freedom, and selects the one at the sharpest elbow, i.e. with the
highest *st* value.

Because it trades fit against parsimony instead of testing against a
null model of uncorrelated variables, the Hull method does not lose
accuracy for the correlated-factor structures where parallel analysis
([`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md))
tends to under-extract; the CAF variant in particular was among the more
accurate criteria in Auerswald and Moshagen (2019). It needs at least
six indicators and fits a model at every candidate factor count, so it
is comparatively slow and is not an option for very short scales.

The
[efa_parallel](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
function and the principal axis factoring of the different number of
factors can be parallelized using the future framework, by calling the
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
function. The examples provide example code on how to enable parallel
processing.

The upper bound *J* comes from
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
which compares against simulated data, so the suggested number of
factors varies slightly from run to run; a criterion-based rotation
passed through `...` adds its own random starts. Call
[`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand to make a
run reproducible; the result is then also independent of the parallel
plan.

Note that if `gof = "RMSEA"` is used, 1 - RMSEA is actually used to
compare the different solutions. This is necessary due to how the
heuristic to locate the elbow of the hull works.

The solutions are fitted without inequality constraints, so a solution
can be inadmissible (a Heywood case, or a fit that did not converge).
Only the selected solution is checked for this; if it is inadmissible a
warning is raised and the retained number of factors should be
interpreted with caution.

The ML estimation method uses the
[`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) starting values.
See also the
[efa_fit](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
documentation.

## See also

[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
as a wrapper function for this and the other factor retention criteria.

Other factor retention criteria:
[`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md),
[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md),
[`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
[`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md),
[`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md),
[`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)

## Examples

``` r
# \donttest{
# using PAF (this will throw a warning if gof is not specified manually
# and CAF will be used automatically)
efa_hull(test_models$baseline$cormat, N = 500, gof = "CAF")
#> ── Hull method ─────────────────────────────────────────────────────────────────
#> Estimation method: PAF
#> 
#> • CAF: 3

# using ML with all available fit indices (CAF, CFI, and RMSEA)
efa_hull(test_models$baseline$cormat, N = 500, estimator = "ML")
#> ── Hull method ─────────────────────────────────────────────────────────────────
#> Estimation method: ML
#> 
#> • CAF: 3
#> • CFI: 1
#> • RMSEA: 1

# using ULS with only RMSEA
efa_hull(test_models$baseline$cormat, N = 500, estimator = "ULS", gof = "RMSEA")
#> ── Hull method ─────────────────────────────────────────────────────────────────
#> Estimation method: ULS
#> 
#> • RMSEA: 1
# }

if (FALSE) { # \dontrun{
# using parallel processing (Note: plans can be adapted, see the future
# package for details)
future::plan(future::multisession)
efa_hull(test_models$baseline$cormat, N = 500, gof = "CAF")
} # }
```
