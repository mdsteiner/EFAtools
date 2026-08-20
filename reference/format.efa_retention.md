# Format method for efa_retention objects

Format method for efa_retention objects

## Usage

``` r
# S3 method for class 'efa_retention'
format(x, ...)
```

## Arguments

- x:

  an object of class efa_retention, returned by a factor-retention
  criterion (e.g.
  [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md)
  or
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)).

- ...:

  not used.

## Value

A character vector with the report lines (styled to the active console
theme; plain when colours are disabled).

## Examples

``` r
writeLines(format(efa_ekc(test_models$baseline$cormat, N = 500)))
#> ── Empirical Kaiser Criterion ──────────────────────────────────────────────────
#> 
#> • Braeken & van Assen (2017): 3
```
