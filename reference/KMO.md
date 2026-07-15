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

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is "pairwise.complete.obs".

- cor_method:

  character. Correlation computed from raw data: `"pearson"`,
  `"spearman"`, or `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)), or `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations of ordinal /
  binary data (a two-step estimator with no empty-cell continuity
  correction). Default is "pearson".

## Value

A list of class `c("efa_kmo", "KMO")`, identical to the value of
[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md);
see there for the components.

## See also

[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md)
