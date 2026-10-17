# Kaiser-Guttman criterion

Retains as many factors as there are eigenvalues above a cutoff. For the
correlation matrix (`"PCA"`), the cutoff is 1: the classic
eigenvalues-greater-than-one rule (Kaiser, 1960). For the reduced
correlation matrix, with communality estimates in the diagonal (`"SMC"`
and `"EFA"`), the cutoff is 0 (Guttman, 1954).

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
analysis. Psychometrika, 19, 149–161. https://doi.org/10.1007/BF02289162

Kaiser, H. F. (1960). The application of electronic computers to factor
analysis. Educational and Psychological Measurement, 20, 141–151.
https://doi.org/10.1177/001316446002000116

Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
determining the number of components to retain. Psychological Bulletin,
99, 432–442. https://doi.org/10.1037/0033-2909.99.3.432

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
  Default is `c("PCA", "SMC", "EFA")`, i.e. all three; `"EFA"` is the
  only one that fits a model.

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
  `estimate_control`, and the standard-error arguments (`se`, `b_boot`,
  `ci`, `seed`) are not accepted because the fit is an internal step
  that keeps only its communalities.

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

The two cutoffs apply to different matrices. The eigenvalues of the
correlation matrix sum to the number of variables, and an eigenvalue
greater than 1 identifies a component that explains more variance than
one standardized variable. The eigenvalues of the reduced correlation
matrix sum to the sum of the communalities, so the cutoff 1 does not
apply to them. With SMCs in the diagonal, the number of positive
eigenvalues is a lower bound to the number of common factors in the
population. This bound is the strongest of the three lower bounds of
Guttman (1954). `"EFA"` applies the same cutoff to the final
communalities of an EFA solution.

Both rules are often inaccurate in samples. Because of sampling error,
the reduced correlation matrix of a sample usually has many small
positive eigenvalues, so the cutoff 0 tends to suggest too many factors.
The classic rule also tends to suggest too many factors (e.g., Zwick &
Velicer, 1986; see Auerswald & Moshagen, 2019, for a comparison with
other criteria). Which criteria are informative depends on the data at
hand, so rather than substituting one for another, run several of them
together and compare their suggestions.

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
#> Eigenvalues found using PCA and SMC; cutoff 1 (PCA), 0 (SMC).
#> 
#> • PCA: 3
#> • SMC: 8
```
