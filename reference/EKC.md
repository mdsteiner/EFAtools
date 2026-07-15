# Empirical Kaiser criterion

**\[superseded\]**

`EKC()` has been superseded by
[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
EKC(
  x,
  N = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  type = "BvA2017"
)
```

## Arguments

- x:

  data.frame or matrix. data.frame or matrix of raw data or matrix with
  correlations.

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
  `"tetra"` for polychoric / tetrachoric correlations of ordinal /
  binary data (a two-step estimator with no empty-cell continuity
  correction). Default is `"pearson"`.

- type:

  character. The calculation of EKC. type `"BvA2017"` is the original
  implementation; type `"AM2019"` differs from the original
  implementation but was used in simulation studies (Auerswald &
  Moshagen, 2019; Caron, 2025). See details. Use
  `type = c("BvA2017", "AM2019")` for both implementations. Make sure to
  report which version you used.

## Value

An object of class `efa_retention`, identical to the value of
[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md);
see there for the components.

## See also

[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md)
