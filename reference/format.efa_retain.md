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
#> ── Sampling adequacy and sphericity ────────────────────────────────────────────
#> 
#> ✔ Overall KMO = 0.916: marvellous (Kaiser & Rice, 1974).
#> ✔ Bartlett's test of sphericity: χ²(153) = 2173.28, p < .001 (significant).
#> 
#> ── Suggested number of factors ─────────────────────────────────────────────────
#> 
#> 4 suggestions from 2 criteria, ranging from 2 to 3 factors (most common: 3).
#> 
#> Empirical Kaiser Criterion
#> • Braeken & van Assen (2017): 3
#> 
#> Sequential model tests
#> • Sequential chi-square model tests: 3
#> • Lower bound of RMSEA 90% CI: 2
#> • Akaike Information Criterion: 3
# }
```
