# Kaiser-Meyer-Olkin criterion

This function computes the Kaiser-Meyer-Olkin (KMO) criterion overall
and for each variable in a correlation matrix. The KMO represents the
degree to which each observed variable is predicted by the other
variables in the dataset and with this indicates the suitability for
factor analysis.

## Usage

``` r
efa_kmo(
  x,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")
)
```

## Source

Kaiser, H. F. (1970). A second generation little jiffy. Psychometrika,
35, 401-415.

Kaiser, H. F. & Rice, J. (1974). Little jiffy, mark IV. Educational and
Psychological Measurement, 34, 111-117.

Cureton, E. E. & D'Agostino, R. B. (1983). Factor analysis: An applied
approach. Hillsdale, N.J.: Lawrence Erlbaum Associates, Inc.

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

A list containing

- KMO:

  Overall KMO.

- KMO_i:

  KMO for each variable.

- settings:

  A list of the settings used.

## Details

Kaiser (1970) proposed this index, originally called measure of sampling
adequacy (MSA), that indicates how near the inverted correlation matrix
\\R^{-1}\\ is to a diagonal matrix to determine a given correlation
matrix's (\\R\\) suitability for factor analysis. The index is \$\$KMO =
\frac{\sum\_{i \neq j} r\_{ij}^2}{\sum\_{i \neq j} r\_{ij}^2 + \sum\_{i
\neq j} q\_{ij}^2}\$\$ with \\Q = SR^{-1}S\\ and S = \\(diag
R^{-1})^{-1/2}\\ where \\\sum\_{i \neq j} r\_{ij}^2\\ is the sum of
squares of the off-diagonal elements of \\R\\ and \\\sum\_{i \neq j}
q\_{ij}^2\\ is the sum of squares of the off-diagonal elements of \\Q\\
(see also Cureton & D'Agostino, 1983).

So KMO varies between 0 and 1, with larger values indicating higher
suitability for factor analysis. Kaiser and Rice (1974) suggest that KMO
should at least exceed .50 for a correlation matrix to be suitable for
factor analysis.

This function was heavily influenced by the
[`psych::KMO()`](https://rdrr.io/pkg/psych/man/KMO.html) function.

See also
[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)
for another test of suitability for factor analysis.

The `efa_kmo` function can also be called together with the
[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)
function and with factor retention criteria in the
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
function.

## See also

[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)
for another measure to determine suitability for factor analysis.

[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
as a wrapper function for this function,
[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)
and several factor retention criteria.

Other factor analysis suitability:
[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md),
[`efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/efa_screen.md),
[`print.efa_screen()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_screen.md)

## Examples

``` r
efa_kmo(test_models$baseline$cormat)
#> 
#> ── Kaiser-Meyer-Olkin criterion (KMO) ──────────────────────────────────────────
#> 
#> ✔ The overall KMO value for your data is marvellous.
#> These data are probably suitable for factor analysis.
#> 
#> Overall: 0.916
#> 
#> For each variable:
#>    V1    V2    V3    V4    V5    V6    V7    V8    V9   V10   V11   V12   V13 
#> 0.900 0.914 0.924 0.932 0.923 0.891 0.928 0.919 0.916 0.892 0.928 0.908 0.922 
#>   V14   V15   V16   V17   V18 
#> 0.905 0.924 0.934 0.907 0.923 
```
