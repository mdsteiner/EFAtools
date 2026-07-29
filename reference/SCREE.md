# Scree plot

**\[superseded\]**

`SCREE()` has been superseded by
[`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
SCREE(
  x,
  eigen_type = c("PCA", "SMC", "EFA"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  n_factors = 1,
  ...
)
```

## Arguments

- x:

  data.frame or matrix. Dataframe or matrix of raw data or matrix with
  correlations.

- eigen_type:

  character. On what the eigenvalues should be found. Can be either
  "PCA", "SMC", or "EFA", or some combination of them. If using "PCA",
  the diagonal values of the correlation matrices are left to be 1. If
  using "SMC", the diagonal of the correlation matrices is replaced by
  the squared multiple correlations (SMCs) of the indicators. If using
  "EFA", eigenvalues are found on the correlation matrices with the
  final communalities of an exploratory factor analysis solution
  (default is principal axis factoring extracting 1 factor) as diagonal.

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

- n_factors:

  numeric. Number of factors to extract if "EFA" is included in
  `eigen_type`. Default is 1.

- ...:

  Further arguments passed on to the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fit. For example, `estimator`, to change the estimator (PAF is
  default), or one of the estimation tuning knobs (`type`, `init_comm`,
  `criterion`, `criterion_type`, `max_iter`, `abs_eigen`,
  `start_method`), which are repacked into an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object so that they tune the fit exactly as they always did.

## Value

An object of class `efa_retention`, identical to the value of
[`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md);
see there for the components.

## See also

[`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md)
