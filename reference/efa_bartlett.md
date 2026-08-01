# Bartlett's test of sphericity

This function tests whether a correlation matrix is significantly
different from an identity matrix (Bartlett, 1951). If the Bartlett's
test is not significant, the correlation matrix is not suitable for
factor analysis because the variables show too little covariance.

## Usage

``` r
efa_bartlett(
  x,
  N = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")
)
```

## Source

Bartlett, M. S. (1951). The effect of standardization on a Chi-square
approximation in factor analysis. Biometrika, 38, 337-344.

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

A list containing

- chisq:

  The chi square statistic, or `NA` if `N` is too small for the Bartlett
  correction (i.e. \\N - 1 - (2p + 5)/6 \le 0\\).

- p_value:

  The p value of the chi square statistic, or `NA` when `chisq` is `NA`.

- df:

  The degrees of freedom for the chi square statistic.

- settings:

  A list of the settings used.

## Details

Bartlett (1951) proposed this statistic to determine a correlation
matrix' suitability for factor analysis. The statistic is approximately
chi square distributed with \\df = \frac{p(p - 1)}{2}\\ and is given by

\$\$chi^2 = -log(det(R)) (N - 1 - (2 \* p + 5)/6)\$\$

where \\det(R)\\ is the determinant of the correlation matrix, \\N\\ is
the sample size, and \\p\\ is the number of variables.

This tests requires multivariate normality. If this condition is not
met, the Kaiser-Meyer-Olkin criterion
([`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md))
can still be used.

This function was heavily influenced by the
[`psych::cortest.bartlett()`](https://rdrr.io/pkg/psych/man/cortest.bartlett.html)
function from the psych package.

The `efa_bartlett` function can also be called together with the
([`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md))
function and with factor retention criteria in the
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
function.

## See also

[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md)
for another measure to determine suitability for factor analysis.

[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
as a wrapper function for this function,
[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md)
and several factor retention criteria.

Other factor analysis suitability:
[`efa_kmo()`](https://mdsteiner.github.io/EFAtools/reference/efa_kmo.md),
[`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md),
[`print.efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_screen.md)

## Examples

``` r
efa_bartlett(test_models$baseline$cormat, N = 500)
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05.
#> These data are probably suitable for factor analysis.
#> 
#> χ²(153) = 2173.28, p < .001
```
