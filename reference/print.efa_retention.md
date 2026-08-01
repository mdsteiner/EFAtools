# Print method for efa_retention objects

Print method for efa_retention objects

## Usage

``` r
# S3 method for class 'efa_retention'
print(x, ...)
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

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly; it is `cat(format(x), sep = "\n")`.

## Examples

``` r
efa_ekc(test_models$baseline$cormat, N = 500)
#> ── Empirical Kaiser Criterion ──────────────────────────────────────────────────
#> 
#> • Original implementation (Braeken & van Assen, 2017): 3
#> 
#> ℹ Multiple implementations of EKC exist; make sure to report which one you used
#> (see the efa_ekc help page for details).
```
