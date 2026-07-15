# Minimum average partial

**\[superseded\]**

`MAP()` has been superseded by
[`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
MAP(
  x,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")
)
```

## Arguments

- x:

  A numeric `matrix` or `data.frame`. Can be either (a) a correlation
  matrix, or (b) raw data (rows = observations, columns = variables)
  from which correlations are computed.

- use:

  Character string specifying the treatment of missing values when
  computing correlations. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). Defaults to
  `"pairwise.complete.obs"`.

- cor_method:

  Character string specifying the correlation coefficient to be computed
  if raw data are supplied. One of `"pearson"`, `"spearman"`, or
  `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)), or `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations of ordinal /
  binary data (a two-step estimator with no empty-cell continuity
  correction). Defaults to `"pearson"`.

## Value

An object of class `efa_retention`, identical to the value of
[`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md);
see there for the components.

## See also

[`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md)
