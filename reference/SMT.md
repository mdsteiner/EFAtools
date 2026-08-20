# Sequential model tests

**\[superseded\]**

`SMT()` has been superseded by
[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
SMT(
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
  correlation matrix is used. Must be larger than the number of
  variables.

- use:

  character. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) if raw data is
  given as input. Default is "pairwise.complete.obs".

- cor_method:

  character. One of `"pearson"`, `"spearman"`, or `"kendall"`, passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). `"poly"` and
  `"tetra"` are not supported because `SMT` rests on a normal-theory
  chi-square test that is not valid for polychoric / tetrachoric
  correlations. Default is `"pearson"`.

## Value

An object of class `efa_retention`, identical to the value of
[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md);
see there for the components.

## See also

[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)
