# Print method for efa_retain objects

Print method for efa_retain objects

## Usage

``` r
# S3 method for class 'efa_retain'
print(x, ...)
```

## Arguments

- x:

  an object of class efa_retain, returned by
  [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md).

- ...:

  not used.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly; it is `cat(format(x), sep = "\n")`.

## Examples

``` r
# \donttest{
efa_retain(test_models$baseline$cormat, criteria = c("EKC", "SMT"), N = 500)
#> ── Tests for the suitability of the data for factor analysis ───────────────────
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05:
#>   χ²(153) = 2173.28, p < .001. These data are probably suitable for factor
#>   analysis.
#> ✔ The Kaiser-Meyer-Olkin criterion is marvellous (KMO = 0.916). These data are
#>   probably suitable for factor analysis.
#> 
#> ── Number of factors suggested by the factor retention criteria ────────────────
#> 
#> Empirical Kaiser Criterion
#> • Original implementation (Braeken & van Assen, 2017): 3
#> 
#> Sequential model tests
#> • Sequential chi-square model tests: 3
#> • Lower bound of RMSEA 90% CI: 2
#> • Akaike Information Criterion: 3
# }
```
