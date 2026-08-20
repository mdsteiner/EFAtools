# Empirical Kaiser criterion

The empirical Kaiser criterion incorporates random sampling variations
of the eigenvalues from the Kaiser-Guttman criterion
([`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md);
see Auerswald & Moshagen, 2019; Braeken & van Assen, 2017). The
implementation follows Braeken and van Assen (2017).

## Usage

``` r
efa_ekc(
  x,
  N = NA,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
  type = lifecycle::deprecated()
)
```

## Source

Auerswald, M., & Moshagen, M. (2019). How to determine the number of
factors to retain in exploratory factor analysis: A comparison of
extraction methods under realistic conditions. Psychological Methods,
24(4), 468–491. https://doi.org/10.1037/met0000200

Braeken, J., & van Assen, M. A. (2017). An empirical Kaiser criterion.
Psychological Methods, 22, 450 – 466. https://doi.org/10.1037/met0000074

Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
determining the number of components to retain. Psychological Bulletin,
99, 432–442. https://doi.org/10.1037/0033-2909.99.3.432

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

An object of class `efa_retention` (see
[`print.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_retention.md)
and
[`plot.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_retention.md)
for the print and plot methods). Its main fields are:

- n_factors:

  A numeric vector of length one, named `"BvA2017"`, with the suggested
  number of factors. The factors up to the first observed eigenvalue
  that fails to exceed its reference value are retained. The
  "all-exceed" convention of parallel analysis
  ([`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md)),
  which retains all *J* factors when no such crossing is found, cannot
  be reached here: the reference values are never below 1, while the
  eigenvalues of a correlation matrix sum to *J* and are sorted
  downwards, so the last of them is never above 1.

- results:

  A list with one record, holding the eigenvalues, the reference
  eigenvalues, and the retained solution used for printing and plotting.

- settings:

  A list with the settings used.

## Details

The Kaiser-Guttman criterion was defined with the intend that a factor
should only be extracted if it explains at least as much variance as a
single factor (see
[`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md)).
However, this only applies to population-level correlation matrices. Due
to sampling variation, the KGC strongly overestimates the number of
factors to retrieve (e.g., Zwick & Velicer, 1986). To account for this
and to introduce a factor retention method that performs well with small
number of indicators and correlated factors (cases where the performance
of parallel analysis, see
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
is known to deteriorate) Braeken and van Assen (2017) introduced the
empirical Kaiser criterion in which a series of reference eigenvalues is
created as a function of the variables-to-sample-size ratio and the
observed eigenvalues.

Braeken and van Assen (2017) showed that "(a) EKC performs about as well
as parallel analysis for data arising from the null, 1-factor, or
orthogonal factors model; and (b) clearly outperforms parallel analysis
for the specific case of oblique factors, particularly whenever factor
intercorrelation is moderate to high and the number of variables per
factor is small, which is characteristic of many applications these
days" (p.463-464).

## See also

[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
as a wrapper function for this and the other factor retention criteria.

Other factor retention criteria:
[`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md),
[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
[`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
[`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md),
[`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md),
[`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)

## Examples

``` r
efa_ekc(test_models$baseline$cormat, N = 500)
#> ── Empirical Kaiser Criterion ──────────────────────────────────────────────────
#> 
#> • Braeken & van Assen (2017): 3
```
