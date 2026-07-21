# Parallel analysis

**\[superseded\]**

`PARALLEL()` has been superseded by
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
PARALLEL(
  x = NULL,
  N = NA,
  n_vars = NA,
  n_datasets = 1000,
  percent = 95,
  eigen_type = c("PCA", "SMC", "EFA"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  decision_rule = c("means", "percentile", "crawford"),
  n_factors = 1,
  ...
)
```

## Arguments

- x:

  matrix or data.frame. The real data to compare the simulated
  eigenvalues against. Must not contain variables of classes other than
  numeric. Can be a correlation matrix or raw data.

- N:

  numeric. The number of cases / observations to simulate. Only has to
  be specified if `x` is either a correlation matrix or `NULL`. If x
  contains raw data, `N` is found from the dimensions of `x`.

- n_vars:

  numeric. The number of variables / indicators to simulate. Only has to
  be specified if `x` is left as `NULL` as otherwise the dimensions are
  taken from `x`.

- n_datasets:

  numeric. The number of datasets to simulate. Default is 1000.

- percent:

  numeric. The percentile to take from the simulated eigenvalues.
  Default is 95.

- eigen_type:

  character. On what the eigenvalues should be found. Can be either
  "SMC", "PCA", or "EFA". If using "SMC", the diagonal of the
  correlation matrix is replaced by the squared multiple correlations
  (SMCs) of the indicators. If using "PCA", the diagonal values of the
  correlation matrices are left to be 1. If using "EFA", eigenvalues are
  found on the correlation matrices with the final communalities of an
  EFA solution as diagonal.

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is "pairwise.complete.obs".

- cor_method:

  character. One of `"pearson"`, `"spearman"`, or `"kendall"`, passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). `"poly"` and
  `"tetra"` are not supported because `PARALLEL` compares the data
  against simulated continuous reference data. Default is "pearson".

- decision_rule:

  character. Which rule to use to determine the number of factors to
  retain. Default is `"means"`, which will use the average simulated
  eigenvalues. `"percentile"`, uses the percentiles specified in
  percent. `"crawford"` uses the 95th percentile for the first factor
  and the mean afterwards (based on Crawford et al, 2010). The `"means"`
  rule retains a factor whenever its real eigenvalue exceeds the average
  simulated one and thus tends to retain more factors than the more
  conservative `"percentile"` rule (Glorfeld, 1995).

- n_factors:

  numeric. Number of factors to extract if "EFA" is included in
  `eigen_type`. Default is 1.

- ...:

  Further arguments passed on to the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fits. For example, `estimator`, to change the estimator (default is
  "PAF"; PAF is more robust, but it will take longer compared to "ML"
  and "ULS"), or one of the estimation tuning knobs (`type`,
  `init_comm`, `criterion`, `criterion_type`, `max_iter`, `abs_eigen`,
  `start_method`), which are repacked into an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object so that they tune the fits exactly as they always did.

## Value

An object of class `efa_retention`, identical to the value of
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md);
see there for the components.

## See also

[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)
