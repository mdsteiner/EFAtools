# Kaiser-Guttman criterion

Probably the most popular factor retention criterion. Kaiser and Guttman
suggested to retain as many factors as there are sample eigenvalues
greater than 1. This is why the criterion is also known as
eigenvalues-greater-than-one rule.

## Usage

``` r
efa_kgc(
  x,
  eigen_type = c("PCA", "SMC", "EFA"),
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  n_factors = 1,
  estimate_control = NULL,
  ...
)
```

## Source

Auerswald, M., & Moshagen, M. (2019). How to determine the number of
factors to retain in exploratory factor analysis: A comparison of
extraction methods under realistic conditions. Psychological Methods,
24(4), 468–491. https://doi.org/10.1037/met0000200

Guttman, L. (1954). Some necessary conditions for common-factor
analysis. Psychometrika, 19, 149 –161.
http://dx.doi.org/10.1007/BF02289162

Kaiser, H. F. (1960). The application of electronic computers to factor
analysis. Educational and Psychological Measurement, 20, 141–151.
http://dx.doi.org/10.1177/001316446002000116

Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
determining the number of components to retain. Psychological Bulletin,
99, 432–442. http://dx.doi.org/10.1037/0033-2909.99.3.432

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
  binary data (a two-step estimator with no empty-cell continuity
  correction). Default is "pearson".

- n_factors:

  numeric. Number of factors to extract if "EFA" is included in
  `eigen_type`. Default is 1.

- estimate_control:

  an
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  object with the estimation settings for the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  fit that provides the communalities when `"EFA"` is included in
  `eigen_type`. `NULL` (default) uses the
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md)
  defaults. The fit is unrotated, so no rotation settings apply.

- ...:

  Additional arguments passed to
  [`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md).
  For example, `estimator`, to change the estimator (PAF is default).
  The estimation tuning knobs are not passed here; they live in
  `estimate_control`.

## Value

An object of class `efa_retention` (see
[`print.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_retention.md)
and
[`plot.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_retention.md)
for the print and plot methods). Its main fields are:

- n_factors:

  A named numeric vector with the suggested number of factors for each
  requested eigenvalue type (`"PCA"`, `"SMC"`, and/or `"EFA"`).

- results:

  A list with one record per eigenvalue type, each holding the
  eigenvalues and the retained solution used for printing and plotting.

- settings:

  A list of the settings used.

## Details

Originally, the Kaiser-Guttman criterion was intended for the use with
principal components, hence with eigenvalues derived from the original
correlation matrix. This can be done here by setting `eigen_type` to
"PCA". However, it is well-known that this criterion is often inaccurate
and that it tends to overestimate the number of factors, especially for
unidimensional or orthogonal factor structures (e.g., Zwick & Velicer,
1986).

The criterion's inaccuracy in these cases is somewhat addressed if it is
applied on the correlation matrix with communalities in the diagonal,
either initial communalities estimated from SMCs (done setting
`eigen_type` to "SMC") or final communality estimates from an EFA (done
setting `eigen_type` to "EFA"; see Auerswald & Moshagen, 2019). However,
although this variant of the KGC is more accurate in some cases compared
to the traditional KGC, it is at the same time less accurate than the
PCA-variant in other cases, and it is still often less accurate than
other factor retention methods, for example parallel analysis
([`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)),
the Hull method
[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
or sequential \\chi^2\\ model tests
([`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md);
see Auerswald & Moshagen, 2019).

The `efa_kgc` function can also be called together with other factor
retention criteria in the
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
function.

## See also

[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
as a wrapper function for this and the other factor retention criteria.

Other factor retention criteria:
[`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md),
[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md),
[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
[`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md),
[`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md),
[`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)

## Examples

``` r
efa_kgc(test_models$baseline$cormat, eigen_type = c("PCA", "SMC"))
#> ── Kaiser-Guttman criterion ────────────────────────────────────────────────────
#> 
#> • PCA eigenvalues: 3
#> • SMC eigenvalues: 1
```
