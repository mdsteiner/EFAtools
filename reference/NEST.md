# Next eigenvalue sufficiency test

**\[superseded\]**

`NEST()` has been superseded by
[`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
NEST(
  x,
  N = NA,
  alpha = 0.05,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  n_datasets = 1000,
  ...
)
```

## Arguments

- x:

  data.frame or matrix. data.frame or matrix of raw data or matrix with
  correlations.

- N:

  numeric. The number of observations. Only needed if x is a correlation
  matrix.

- alpha:

  numeric. The alpha level to use (i.e., 1-alpha percentile of
  eigenvalues is used for reference values).

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is `"pairwise.complete.obs"`.

- cor_method:

  character. One of `"pearson"`, `"spearman"`, or `"kendall"`, passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). `"poly"` and
  `"tetra"` are not supported because `NEST` compares the data against
  simulated continuous reference data. Default is `"pearson"`.

- n_datasets:

  numeric. The number of datasets to simulate. Default is 1000.

- ...:

  Additional arguments passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).
  For example, the extraction method can be changed here (default is
  "PAF"). PAF is more robust, but it will take longer compared to the
  other estimation methods available ("ML" and "ULS"). The estimation
  tuning knobs are not passed here; they live in `estimate_control`.

## Value

An object of class `efa_retention`, identical to the value of
[`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md);
see there for the components.

## See also

[`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md)
