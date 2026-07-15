# Print and format an efa_power object

[`print()`](https://rdrr.io/r/base/print.html) turns an
[`efa_power()`](https://mdsteiner.github.io/EFAtools/reference/efa_power.md)
result into a short report. For an RMSEA-mode object this is a header
naming the test, the null and alternative hypotheses with the
significance level and degrees of freedom, the headline result (the
power at the sample size, or the required sample size for the target
power), and the critical value and noncentrality parameters. For a
simulation-mode object it is instead the population and design, the
retention hit-rate per criterion, the structure-recovery rate, and the
convergence and Heywood-case rate.
[`format()`](https://rdrr.io/r/base/format.html) assembles the same
report and returns it as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled (for example when captured
into a file or stripped with
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html)).

## Usage

``` r
# S3 method for class 'efa_power'
print(x, digits = 3, ...)

# S3 method for class 'efa_power'
format(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `efa_power` (output from
  [`efa_power()`](https://mdsteiner.github.io/EFAtools/reference/efa_power.md)).

- digits:

  Integer. The number of decimal places the reported values are rounded
  to. Default is 3.

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector with the report lines (styled to the active console
theme; plain when colours are disabled).

## See also

Other power analysis:
[`efa_power()`](https://mdsteiner.github.io/EFAtools/reference/efa_power.md),
[`plot.efa_power()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_power.md)

## Examples

``` r
pw <- efa_power(df = 100, N = 200)
pw
#> 
#> ── RMSEA power analysis ────────────────────────────────────────────────────────
#> 
#> Test of close fit: H0 RMSEA ≤ .050 vs. H1 RMSEA = .080.
#> alpha = .050 · df = 100
#> 
#> Power = .955 at N = 200.
#> Critical value χ²(100) = 183.967 · noncentrality H0 = 49.750, H1 = 127.360.

# format() returns the same lines as plain text:
writeLines(format(pw))
#> 
#> ── RMSEA power analysis ────────────────────────────────────────────────────────
#> 
#> Test of close fit: H0 RMSEA ≤ .050 vs. H1 RMSEA = .080.
#> alpha = .050 · df = 100
#> 
#> Power = .955 at N = 200.
#> Critical value χ²(100) = 183.967 · noncentrality H0 = 49.750, H1 = 127.360.
```
