# Print and format a multigroup factor analysis

[`print()`](https://rdrr.io/r/base/print.html) turns an
[`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md)
result into a sectioned report: a header recapping the groups, the
common number of factors, the estimator, the rotation, and the
alignment; a group-pair table of the matched Tucker congruences between
the aligned loadings; a per-pair summary of the cross-group loading
differences (with the salient and, when a bootstrap was run, the
confidence-interval flags); and, when `invariance = TRUE`, a group-pair
by factor grid of the approximate-invariance verdicts.
[`format()`](https://rdrr.io/r/base/format.html) assembles the same
report and returns it as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled (for example when captured
into a file or stripped with
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html)).
[`print()`](https://rdrr.io/r/base/print.html) does not draw a plot; use
[`plot.efa_group()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_group.md).

## Usage

``` r
# S3 method for class 'efa_group'
print(x, digits = 3, ...)

# S3 method for class 'efa_group'
format(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `efa_group` (output from
  [`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md)).

- digits:

  Integer. The number of decimal places the reported values are rounded
  to. Default is 3.

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector with the report lines.

## See also

Other factor analysis:
[`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md),
[`efa_fit()`](https://mdsteiner.github.io/EFAtools/reference/efa_fit.md),
[`efa_group()`](https://mdsteiner.github.io/EFAtools/reference/efa_group.md),
[`efa_mi()`](https://mdsteiner.github.io/EFAtools/reference/efa_mi.md),
[`plot.efa_group()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_group.md)

## Examples

``` r
g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
mg <- efa_group(GRiPS_raw, groups = g, n_factors = 1)
#> ℹ `x` is not a correlation matrix; computing correlations from the raw data.
mg
#> ── Multigroup exploratory factor analysis ──────────────────────────────────────
#> 
#> 2 groups (g1, g2) · 1 factor · PAF extraction · unrotated
#> Aligned to a symmetric consensus target.
#> N: g1 = 405, g2 = 405
#> 
#> ── Factor congruence ───────────────────────────────────────────────────────────
#> 
#> Tucker's congruence between the aligned group loadings (matched factors):
#> 
#>            F1
#> g1 vs g2  .999
#> 
#> ── Loading differences ─────────────────────────────────────────────────────────
#> 
#> Absolute loading differences by group pair (salience threshold |Δ| ≥ .1):
#> 
#>           mean  median   min   max  rmse  flagged
#> g1 vs g2  .025    .017  .002  .075  .033      0/8

# format() returns the same lines as a character vector:
writeLines(format(mg))
#> ── Multigroup exploratory factor analysis ──────────────────────────────────────
#> 
#> 2 groups (g1, g2) · 1 factor · PAF extraction · unrotated
#> Aligned to a symmetric consensus target.
#> N: g1 = 405, g2 = 405
#> 
#> ── Factor congruence ───────────────────────────────────────────────────────────
#> 
#> Tucker's congruence between the aligned group loadings (matched factors):
#> 
#>            F1
#> g1 vs g2  .999
#> 
#> ── Loading differences ─────────────────────────────────────────────────────────
#> 
#> Absolute loading differences by group pair (salience threshold |Δ| ≥ .1):
#> 
#>           mean  median   min   max  rmse  flagged
#> g1 vs g2  .025    .017  .002  .075  .033      0/8
```
