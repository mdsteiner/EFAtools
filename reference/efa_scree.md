# Scree plot

The scree plot was originally introduced by Cattell (1966) to perform
the scree test. In a scree plot, the eigenvalues of the factors /
components are plotted against the index of the factors / components,
ordered from 1 to N factors components, hence from largest to smallest
eigenvalue. According to the scree test, the number of factors /
components to retain is the number of factors / components to the left
of the "elbow" (where the curve starts to level off) in the scree plot.

## Usage

``` r
efa_scree(
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

Cattell, R. B. (1966). The scree test for the number of factors.
Multivariate Behavioral Research, 1(2), 245–276.
https://doi.org/10.1207/s15327906mbr0102_10

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
  binary data (a two-step estimator). Default is "pearson".

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
for the print and plot methods). The scree plot is a visual criterion,
so it returns no numeric suggestion. Its main fields are:

- results:

  A list with one record per requested eigenvalue type, each holding the
  eigenvalues used for the scree plot.

- settings:

  A list of the settings used.

## Details

As the scree test requires visual examination, the test has been
especially criticized for its subjectivity and with this low inter-rater
reliability. Moreover, a scree plot can be ambiguous if there are either
no clear "elbow" or multiple "elbows", making it difficult to judge just
where the eigenvalues do level off. Finally, the scree test has also
been found to be less accurate than other factor retention criteria. For
all these reasons, the scree test has been recommended against, at least
for exclusive use as a factor retention criterion (Zwick & Velicer,
1986)

The `efa_scree` function can also be called together with other factor
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
[`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
[`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md),
[`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md),
[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)

## Examples

``` r
efa_scree(test_models$baseline$cormat, eigen_type = c("PCA", "SMC"))
#> ── Scree plot ──────────────────────────────────────────────────────────────────
#> Eigenvalues found using PCA and SMC.
#> 
#> ℹ Scree plot is a visual criterion; inspect the plot to identify the elbow.
```
