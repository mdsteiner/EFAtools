# Print and format an efa_bartlett object

[`print()`](https://rdrr.io/r/base/print.html) reports the outcome of
[`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)'s
test of sphericity: a verdict on whether the test was significant (and
what that implies for the suitability of the data for factor analysis),
followed by the chi-square statistic, its degrees of freedom, and the
p-value. [`format()`](https://rdrr.io/r/base/format.html) assembles the
same report and returns it as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled (for example when captured
into a file or stripped with
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html)).

## Usage

``` r
# S3 method for class 'efa_bartlett'
print(x, ...)

# S3 method for class 'efa_bartlett'
format(x, ...)
```

## Arguments

- x:

  An object of class `efa_bartlett` (output from
  [`efa_bartlett()`](https://mdsteiner.github.io/EFAtools/reference/efa_bartlett.md)).

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector with the report lines.

## Examples

``` r
bart <- efa_bartlett(test_models$baseline$cormat, N = 500)
bart
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05.
#> These data are probably suitable for factor analysis.
#> 
#> χ²(153) = 2173.28, p < .001

# format() returns the same lines as a character vector:
writeLines(format(bart))
#> 
#> ✔ The Bartlett's test of sphericity was significant at an alpha level of .05.
#> These data are probably suitable for factor analysis.
#> 
#> χ²(153) = 2173.28, p < .001
```
