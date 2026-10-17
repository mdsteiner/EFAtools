# Hull method

**\[superseded\]**

`HULL()` has been superseded by
[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
HULL(
  x,
  N = NA,
  n_fac_theor = NA,
  method = c("PAF", "ULS", "ML"),
  gof = c("CAF", "CFI", "RMSEA"),
  eigen_type = c("SMC", "PCA", "EFA"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  n_datasets = 1000,
  percent = 95,
  decision_rule = c("means", "percentile", "crawford"),
  n_factors = 1,
  ...
)
```

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

- method:

  character. The estimator to use; passed to
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)
  as its `estimator` argument. One of `"PAF"` (default), `"ULS"`, or
  `"ML"`.

- gof:

  character. Passed to
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md).
  The goodness of fit index to use: `"CAF"`, `"CFI"`, or `"RMSEA"`, or
  any combination of them. Default is all three. With the `"PAF"`
  method, the Hull method uses only the CAF.

- eigen_type:

  character. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  in
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md).
  One of `"SMC"` (default), `"PCA"`, or `"EFA"`.

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

  character. Passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
  in
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md).
  Which rule to use to determine the number of factors to retain.
  Default is `"means"`, which uses the average simulated eigenvalues.
  `"percentile"` uses the percentiles specified in `percent`.
  `"crawford"` uses the 95th percentile for the first factor and the
  mean afterwards (based on Crawford et al., 2010).

- n_factors:

  numeric. Number of factors to extract if `"EFA"` is included in
  `eigen_type`. Default is 1. This is passed to
  [`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md).

- ...:

  Further arguments passed on to the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fits, including the estimation tuning knobs (`type`, `init_comm`,
  `criterion`, `criterion_type`, `max_iter`, `abs_eigen`,
  `start_method`), which are repacked into an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object so that they tune the fits exactly as they always did. The
  estimator is selected with `method`.

## Value

An object of class `efa_retention`, identical to the value of
[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md);
see there for the components.

## See also

[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)
