# Print and format a reliability object

[`print()`](https://rdrr.io/r/base/print.html) shows the reliability
coefficients for the general factor and the group factors, for a single
group or for each group: the reliability coefficients (omega total,
omega hierarchical, and omega subscale, standardized Cronbach's alpha,
and the H index) and the common-variance indices (the explained common
variance, ECV, and the percent of uncontaminated correlations, PUC).
[`format()`](https://rdrr.io/r/base/format.html) assembles the same
report and returns it as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled (for example when captured
into a file or stripped with
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html)).

## Usage

``` r
# S3 method for class 'efa_reliability'
print(x, digits = 3, ...)

# S3 method for class 'efa_reliability'
format(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `efa_reliability`.

- digits:

  Integer. The number of decimal places the coefficients are rounded to.
  Default is 3.

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector with the report lines (styled to the active console
theme; plain when colours are disabled).

## See also

Other reliability coefficients:
[`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md)

## Examples

``` r
efa_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                   estimator = "PAF", rotation = "promax")
rel <- efa_reliability(efa_mod)
rel
#> 
#> Total variance from the correlation matrix.
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>      tot   sub  alpha    H
#> g   .883         .868
#> F1  .734  .734   .768  .760
#> F2  .680  .680   .763  .753
#> F3  .667  .667   .743  .738

# format() returns the same lines as plain text:
writeLines(format(rel))
#> 
#> Total variance from the correlation matrix.
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>      tot   sub  alpha    H
#> g   .883         .868
#> F1  .734  .734   .768  .760
#> F2  .680  .680   .763  .753
#> F3  .667  .667   .743  .738
```
