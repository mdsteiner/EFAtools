# Comparison data

**\[superseded\]**

`CD()` has been superseded by
[`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md),
which is the recommended interface going forward. It remains available
and unchanged so existing code keeps working.

## Usage

``` r
CD(
  x,
  n_factors_max = NA,
  N_pop = 10000,
  N_samples = 500,
  alpha = 0.3,
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  max_iter = 50
)
```

## Arguments

- x:

  data.frame or matrix. Dataframe or matrix of raw data.

- n_factors_max:

  numeric. The maximum number of factors to test against. Larger numbers
  will increase the duration the procedure takes, but test more possible
  solutions. If left NA (default) the maximum number of factors for
  which the model is still over-identified (df \> 0) is used.

- N_pop:

  numeric. Size of finite populations of comparison data. Default is
  10000.

- N_samples:

  numeric. Number of samples drawn from each population. Default is 500.

- alpha:

  numeric. The alpha level used to test the significance of the
  improvement added by an additional factor. Default is .30.

- cor_method:

  character. One of `"pearson"`, `"spearman"`, or `"kendall"`, passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). `"poly"` and
  `"tetra"` are not supported because `CD` compares the data against
  simulated continuous reference data. Default is "pearson".

- max_iter:

  numeric. The maximum number of iterations after which the iterative
  PAF procedure inside the comparison-data generation is halted; it does
  not cap an EFA of `x`. Default is 50.

## Value

An object of class `efa_retention`, identical to the value of
[`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md);
see there for the components.

## See also

[`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md)
