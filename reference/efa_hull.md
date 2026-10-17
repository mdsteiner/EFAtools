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
  estimator = c("ML", "PAF", "ULS"),
  gof = "CAF",
  eigen_type = c("PCA", "SMC", "EFA"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  n_datasets = 1000,
  percent = 95,
  decision_rule = c("percentile", "means", "crawford"),
  n_factors = 1,
  estimate_control = NULL,
  ...
)
```

## Source

Horn, J. L. (1965). A rationale and test for the number of factors in
factor analysis. Psychometrika, 30(2), 179–185.
https://doi.org/10.1007/BF02289447

Lim, S., & Jahng, S. (2019). Determining the number of factors using
parallel analysis and its recent variants. Psychological Methods, 24(4),
452–467. https://doi.org/10.1037/met0000230

Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011). The Hull
method for selecting the number of common factors. Multivariate
Behavioral Research, 46(2), 340-364.
https://doi.org/10.1080/00273171.2011.564527

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

  character. The estimator to use: `"ML"` (default), `"PAF"`, or
  `"ULS"`, for maximum likelihood, principal axis factoring, and
  unweighted least squares, respectively.

- gof:

  character. The goodness of fit index to use: `"CAF"` (default),
  `"CFI"`, or `"RMSEA"`, or any combination of them. With the `"PAF"`
  estimator, the Hull method uses only the CAF. In Lorenzo-Seva,
  Timmerman, and Kiers (2011), with factor correlations up to .40, the
  CFI was the most accurate index. With strongly correlated factors
  (e.g., .70), CFI and RMSEA can make the Hull method retain too few
  factors. For details on the CAF, see Lorenzo-Seva et al. (2011).

- eigen_type:

  character. The eigenvalues of the parallel analysis that sets the
  upper bound *J*: `"PCA"` (default) uses the correlation matrix,
  `"SMC"` puts the squared multiple correlations on its diagonal, and
  `"EFA"` puts the communalities of an EFA solution on its diagonal (see
  details). This is passed to
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

  numeric. The number of reference correlation matrices to simulate.
  Must be at least 1. Default is 1000. This is passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md).

- percent:

  numeric. The percentile to take from the simulated eigenvalues.
  Default is 95. This is passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md).

- decision_rule:

  character. The rule of the parallel analysis that sets the upper bound
  *J*. Default is `"percentile"`, which uses the percentile of the
  simulated eigenvalues that `percent` sets (see details). `"means"`
  uses the mean simulated eigenvalues. `"crawford"` uses the 95th
  percentile for the first factor and the mean for the other factors
  (based on Crawford et al., 2010). This is passed to
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
  and the retained solution used for printing and plotting. Each record
  also carries `st`, the elbow sharpness of every solution (see
  details): the retained solution has the largest value, and the
  runner-up shows how close the selection was. `st` is `NA` for the
  solutions where it is undefined, that is for those not on the hull and
  for the two hull endpoints, which have no neighbouring hull solution
  on one side. When fewer than three solutions remain on the hull, `st`
  is undefined throughout and the whole vector is `NA`; the retained
  solution is then the one with the highest goodness of fit, and a
  warning says so.

- settings:

  A list of the settings used, including `n_fac_max`, the upper bound
  *J* of the number of factors to extract, and the `decision_rule` and
  `percent` of the parallel analysis that set it (see details). For
  backwards compatibility the estimator is also repeated in
  `settings$method`.

For backwards compatibility the per-index suggestions are additionally
available as the top-level fields `n_fac_CAF`, `n_fac_CFI` and
`n_fac_RMSEA`, each `NA` if that index was not requested in `gof`. New
code should read them from `n_factors` instead.

## Details

The Hull method aims to find a model with an optimal balance between
model fit and number of parameters, retaining only major factors
(Lorenzo-Seva, Timmerman, & Kiers, 2011). It fits 0 to *J* factors –
where *J* is the number of factors suggested by parallel analysis (or
`n_fac_theor`, if that is larger), plus one – keeps the solutions on the
upper boundary of the convex hull of goodness-of-fit against degrees of
freedom, and selects the one at the sharpest elbow, i.e. with the
highest *st* value. *J* is at least 3 (or, with few indicators, the
largest number of factors that gives an over-identified model), so that
the hull has sufficient solutions.

The upper bound *J* helps to prevent overfactoring. The selected elbow
is never the first or the last solution of the range, so the Hull method
retains at most *J* - 1 factors. If fewer than three solutions are on
the hull, there is no elbow; the solution with the highest goodness of
fit is then selected, and a warning says so.

As Lorenzo-Seva et al. (2011) suggest, this parallel analysis (Horn,
1965) uses PCA eigenvalues and the 95th percentile of the simulated
eigenvalues by default. With the eigenvalues of the reduced correlation
matrix (`"SMC"`), parallel analysis often retains too many factors,
particularly when the model does not fit the data exactly. In Lim and
Jahng (2019), parallel analysis with PCA eigenvalues was the most
accurate in most conditions.

In the simulation of Lorenzo-Seva et al. (2011), with factor
correlations up to .40, all Hull variants were more accurate overall
than parallel analysis. The Hull method needs at least six indicators
and fits a model at every candidate factor count, so it is comparatively
slow and is not an option for very short scales.

The
[efa_parallel](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
function and the fits of the different numbers of factors can be
parallelized using the future framework, by calling the
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
function. The examples provide example code on how to enable parallel
processing.

The upper bound *J* comes from
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
which compares against simulated data, so the suggested number of
factors varies slightly from run to run; a criterion-based rotation
passed through `...` adds its own random starts. Call
[`base::set.seed()`](https://rdrr.io/r/base/Random.html) beforehand to
make a run reproducible; the result is then also independent of the
parallel plan.

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
# using the defaults: ML, the CAF, and PCA eigenvalues with the 95th percentile in
# the parallel analysis
efa_hull(test_models$baseline$cormat, N = 500, n_datasets = 100)
#> ── Hull method ─────────────────────────────────────────────────────────────────
#> Estimator: ML
#> 
#> • CAF: 3

# using ML with all available fit indices (CAF, CFI, and RMSEA)
efa_hull(test_models$baseline$cormat, N = 500, gof = c("CAF", "CFI", "RMSEA"),
         n_datasets = 100)
#> ── Hull method ─────────────────────────────────────────────────────────────────
#> Estimator: ML
#> 
#> • CAF: 3
#> • CFI: 1
#> • RMSEA: 1

# using PAF, with which the Hull method uses only the CAF
efa_hull(test_models$baseline$cormat, N = 500, estimator = "PAF", n_datasets = 100)
#> ── Hull method ─────────────────────────────────────────────────────────────────
#> Estimator: PAF
#> 
#> • CAF: 3

# using ULS with only RMSEA
efa_hull(test_models$baseline$cormat, N = 500, estimator = "ULS", gof = "RMSEA",
         n_datasets = 100)
#> ── Hull method ─────────────────────────────────────────────────────────────────
#> Estimator: ULS
#> 
#> • RMSEA: 1
# }

if (FALSE) { # \dontrun{
# using parallel processing (Note: plans can be adapted, see the future
# package for details). future::plan() returns the plan it replaces, so
# on.exit() puts the session back as it was -- also if the call fails.
local({
  old_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)
  efa_hull(test_models$baseline$cormat, N = 500)
})
} # }
```
