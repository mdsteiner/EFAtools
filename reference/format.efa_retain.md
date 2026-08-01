# Format method for efa_retain objects

Format method for efa_retain objects

## Usage

``` r
# S3 method for class 'efa_retain'
format(x, ...)
```

## Arguments

- x:

  an object of class efa_retain, returned by
  [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md).

- ...:

  not used.

## Value

A character vector with the report lines (styled to the active console
theme; plain when colours are disabled).

## Examples

``` r
# \donttest{
nf <- efa_retain(test_models$baseline$cormat, criteria = c("EKC", "SMT"),
                 N = 500)
writeLines(format(nf))
#> ── Tests for the suitability of the data for factor analysis ───────────────────
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05:
#>   χ²(153) = 2173.28, p < .001. These data are probably suitable for factor
#>   analysis.
#> ✔ The Kaiser-Meyer-Olkin criterion is marvellous (KMO = 0.916). These data are
#>   probably suitable for factor analysis.
#> 
#> ── Suggested number of factors ─────────────────────────────────────────────────
#> 
#> 4 suggestions from 2 criteria, ranging from 2 to 3 factors (most common: 3).
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
