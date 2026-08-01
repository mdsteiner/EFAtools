# Velicer's minimum average partial (MAP) criterion

Computes Velicer's Minimum Average Partial (MAP) criterion for
determining the number of factors/components to retain. The function
implements the original MAP criterion (Velicer, 1976), expressed via the
\\\mathrm{TR2}\\ representation, and the revised \\\mathrm{TR4}\\
variant proposed by Velicer, Eaton, and Fava (2000).

## Usage

``` r
efa_map(
  x,
  use = c("pairwise.complete.obs", "all.obs", "complete.obs", "everything",
    "na.or.complete"),
  cor_method = c("pearson", "spearman", "kendall", "poly", "tetra")
)
```

## Source

Auerswald, M., & Moshagen, M. (2019). How to determine the number of
factors to retain in exploratory factor analysis: A comparison of
extraction methods under realistic conditions. *Psychological Methods,
24*(4), 468–491. https://doi.org/10.1037/met0000200

Velicer, W. F. (1976). Determining the number of components from the
matrix of partial correlations. *Psychometrika, 41*, 321–327.

Velicer, W. F., Eaton, C. A., & Fava, J. L. (2000). Construct
explication through factor or component analysis: A review and
evaluation of alternative procedures for determining the number of
factors or components. In Goffin, R. D. & Helmes, E. (Eds.), *Problems
and Solutions in Human Assessment: Honoring Douglas N. Jackson at
Seventy* (pp. 41–71). Boston: Kluwer.

Zwick, W. R., & Velicer, W. F. (1986). Comparison of five rules for
determining the number of components to retain. *Psychological Bulletin,
99*, 432–442. http://dx.doi.org/10.1037/0033-2909.99.3.432

## Arguments

- x:

  A numeric `matrix` or `data.frame`. Can be either (a) a correlation
  matrix, or (b) raw data (rows = observations, columns = variables)
  from which correlations are computed.

- use:

  Character string specifying the treatment of missing values when
  computing correlations. Passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). Defaults to
  `"pairwise.complete.obs"`.

- cor_method:

  Character string specifying the correlation coefficient to be computed
  if raw data are supplied. One of `"pearson"`, `"spearman"`, or
  `"kendall"` (passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)), or `"poly"` /
  `"tetra"` for polychoric / tetrachoric correlations of ordinal /
  binary data (a two-step estimator). Defaults to `"pearson"`.

## Value

An object of class `efa_retention` (see
[`print.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/print.efa_retention.md)
for the print method). MAP has no plot;
[`plot.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_retention.md)
returns `NULL` with a message for it. Its main elements are:

- `n_factors`: A named numeric vector (`"TR2"`, `"TR4"`) with the index
  \\m\\ that minimizes the original (TR2) and revised (TR4) MAP
  criterion.

- `results`: A list with one record per criterion, each holding the
  criterion values over \\m\\.

- `settings`: A list containing `use` and `cor_method`.

## Details

MAP partials successive principal components out of the correlation
matrix and, after removing \\m\\ components, summarizes the off-diagonal
partial correlations \\r^\*\_{ij}\\ that remain in the \\m\\-th partial
correlation matrix \\M\\ (which has a unit diagonal); the suggested
number of factors is the \\m\\ that minimizes the criterion. Two
criteria are returned, each rescaling the trace of a matrix power of
\\M\\ by the number of off-diagonal cells \\p(p-1)\\:

- **TR2 (original MAP; Velicer, 1976):** the average squared
  off-diagonal partial correlation, \$\$\mathrm{TR2}\_m =
  \frac{\mathrm{tr}(M^2) - p}{p(p-1)} = \frac{\sum\_{i \neq j}
  (r^\*\_{ij})^2}{p(p-1)},\$\$ where subtracting \\p\\ removes the \\p\\
  unit diagonal entries.

- **TR4 (revised MAP; Velicer, Eaton, & Fava, 2000):** the analogous
  fourth-power summary, formed from the trace of the fourth matrix
  power, \$\$\mathrm{TR4}\_m = \frac{\mathrm{tr}(M^4) - p}{p(p-1)}.\$\$
  Moving from the squared to the fourth power downweights the small
  partial correlations relative to the large ones, which can sharpen the
  minimum. Unlike TR2, \\\mathrm{tr}(M^4)\\ is *not* the sum of the
  fourth powers of the individual partial correlations; the matrix power
  is intended and is what Velicer, Eaton, and Fava (2000) describe.

MAP is most dependable when the components are well determined, that is
with many indicators per factor and substantial loadings. It has a
well-documented tendency to under-extract, particularly with few
indicators per factor or weak loadings (Zwick & Velicer, 1986; Auerswald
& Moshagen, 2019), so it is best read as a lower bound and paired with a
criterion that errs in the other direction, such as the Kaiser-Guttman
criterion
([`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md)).

A non-positive-definite input correlation matrix (e.g. from sampling
error) is smoothed with
[`psych::cor.smooth()`](https://rdrr.io/pkg/psych/man/cor.smooth.html).

## See also

[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
as a wrapper function for this and the other factor retention criteria.

Other factor retention criteria:
[`efa_cd()`](https://mdsteiner.github.io/EFAtools/reference/efa_cd.md),
[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md),
[`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md),
[`efa_kgc()`](https://mdsteiner.github.io/EFAtools/reference/efa_kgc.md),
[`efa_nest()`](https://mdsteiner.github.io/EFAtools/reference/efa_nest.md),
[`efa_parallel()`](https://mdsteiner.github.io/EFAtools/reference/efa_parallel.md),
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md),
[`efa_scree()`](https://mdsteiner.github.io/EFAtools/reference/efa_scree.md),
[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md)

## Examples

``` r
## Example with raw data
res <- efa_map(GRiPS_raw)
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
res
#> ── Minimum average partial ─────────────────────────────────────────────────────
#> 
#> • Original implementation (TR2): 1
#> • Revised implementation (TR4): 1

## Example with a correlation matrix
res2 <- efa_map(test_models$baseline$cormat)
res2
#> ── Minimum average partial ─────────────────────────────────────────────────────
#> 
#> • Original implementation (TR2): 1
#> • Revised implementation (TR4): 3
```
