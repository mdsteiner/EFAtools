# Kaiser-Meyer-Olkin criterion

**\[superseded\]**

`KMO()` has been superseded by
[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
KMO(
  x,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")
)
```

## Arguments

- x:

  data.frame or matrix. Dataframe or matrix of raw data or matrix with
  correlations.

- use:

  character. The missing-data policy for raw data. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) for `"pearson"`,
  `"spearman"`, and `"kendall"`; for `"poly"` / `"tetra"` the same
  policies are applied to the raw data before the polychoric estimation,
  where `"all.obs"` and `"everything"` abort on a missing value instead
  of returning `NA` correlations. Default is "pairwise.complete.obs".

- cor_method:

  character. Correlation computed from raw data: `"pearson"`,
  `"spearman"`, or `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)), or `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations of ordinal /
  binary data (a two-step estimator). Default is "pearson".

## Value

A list of class `c("efa_kmo", "KMO")`, identical to the value of
[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md);
see there for the components.

## See also

[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md)
