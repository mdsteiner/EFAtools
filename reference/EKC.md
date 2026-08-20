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
  type = lifecycle::deprecated()
)
```

## Arguments

- x:

  data.frame or matrix. Dataframe or matrix of raw data or matrix with
  correlations.

- N:

  numeric. The number of observations. Only needed if x is a correlation
  matrix. Must be larger than the number of variables.

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is "pairwise.complete.obs".

- cor_method:

  character. Correlation computed from raw data: `"pearson"`,
  `"spearman"`, or `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)), or `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations of ordinal /
  binary data (a two-step estimator). Default is `"pearson"`. Note that
  the EKC reference values rest on the Marchenko-Pastur law for the
  eigenvalues of a sample correlation matrix of independent variables,
  which assumes the sampling behaviour of product-moment correlations;
  with `"poly"` / `"tetra"` (and, to a lesser degree, the rank-based
  methods) the reference series is therefore an approximation.

- type:

  **\[deprecated\]** Accepted and ignored. It selected between two ways
  to compute the reference values. The `"AM2019"` reference values do
  not depend on the observed eigenvalues, so they do not apply the
  empirical correction that defines the criterion, and they are no
  longer computed.

## Value

An object of class `efa_retention`, identical to the value of
[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md);
see there for the components.

## See also

[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md)
