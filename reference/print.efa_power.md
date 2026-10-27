# Print and format an efa_power object

[`print()`](https://rdrr.io/r/base/print.html) shows an
[`efa_power()`](https://mdsteiner.github.io/EFAtools/reference/efa_power.md)
result as a short report, and
[`format()`](https://rdrr.io/r/base/format.html) gives the same report
as a character vector ([`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`).

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

  Integer. The number of decimal places. Default is 3.

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.
[`format()`](https://rdrr.io/r/base/format.html) returns a character
vector with the report lines.

## Details

For an RMSEA-mode object, the report shows the test, the hypotheses, the
significance level, the degrees of freedom, and the power at the sample
size or the required sample size.

For a simulation-mode object, the report shows the population and
design, the hit-rate per criterion, the structure recovery, the
convergence, and an upper bound for the Monte Carlo standard error of
the rates.

The lines follow the active console theme, so they are plain text when
colours are off, for example in a file, or after
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html).

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

# format() returns the same lines as a character vector:
writeLines(format(pw))
#> 
#> ── RMSEA power analysis ────────────────────────────────────────────────────────
#> 
#> Test of close fit: H0 RMSEA ≤ .050 vs. H1 RMSEA = .080.
#> alpha = .050 · df = 100
#> 
#> Power = .955 at N = 200.
```
