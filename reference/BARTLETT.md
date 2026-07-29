# Bartlett's test of sphericity

**\[superseded\]**

`BARTLETT()` has been superseded by
[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
BARTLETT(
  x,
  N = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")
)
```

## Arguments

- x:

  data.frame or matrix. Dataframe or matrix of raw data or matrix with
  correlations.

- N:

  numeric. The number of observations. Needs only be specified if a
  correlation matrix is used.

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is "pairwise.complete.obs".

- cor_method:

  character. Correlation computed from raw data: `"pearson"`,
  `"spearman"`, or `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)), or `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations of ordinal /
  binary data (a two-step estimator). Default is "pearson".

## Value

A list of class `c("efa_bartlett", "BARTLETT")`, identical to the value
of
[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md);
see there for the components.

## See also

[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)
