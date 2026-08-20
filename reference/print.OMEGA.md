# Print and format an OMEGA object

[`print()`](https://rdrr.io/r/base/print.html) shows the omega
coefficients computed by
[`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md):
omega total (and, for multi-factor solutions, omega hierarchical, omega
subscale, the H index, the explained common variance, and the percent of
uncontaminated correlations) for the general factor and the group
factors, for a single group or for each group.
[`format()`](https://rdrr.io/r/base/format.html) assembles the same
report and returns it as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled (for example when captured
into a file or stripped with
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html)).

## Usage

``` r
# S3 method for class 'OMEGA'
print(x, digits = 3, ...)

# S3 method for class 'OMEGA'
format(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `OMEGA` (output from
  [`OMEGA()`](https://mdsteiner.github.io/EFAtools/reference/OMEGA.md)).

- digits:

  Integer. The number of decimal places the coefficients are rounded to
  (passed to [`base::round()`](https://rdrr.io/r/base/Round.html)).
  Default is 3.

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector with the report lines.

## Examples

``` r
efa_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                   estimator = "PAF", rotation = "promax")
sl_mod <- efa_schmid_leiman(efa_mod, estimator = "PAF")

om <- OMEGA(sl_mod, type = "EFAtools",
            factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2)
om
#> Omega total, omega hierarchical, omega subscale, H index, explained common
#> variance (ECV), and percent of uncontaminated correlations (PUC) for the
#> general factor (top row) and omegas and H index for the group factors:
#> 
#>      tot  hier   sub     H   ECV   PUC
#> g  0.883 0.740 0.125 0.842 0.652 0.706
#> F1 0.769 0.500 0.269 0.463            
#> F2 0.765 0.494 0.270 0.472            
#> F3 0.745 0.519 0.225 0.408            

# format() returns the same lines as a character vector:
writeLines(format(om))
#> Omega total, omega hierarchical, omega subscale, H index, explained common
#> variance (ECV), and percent of uncontaminated correlations (PUC) for the
#> general factor (top row) and omegas and H index for the group factors:
#> 
#>      tot  hier   sub     H   ECV   PUC
#> g  0.883 0.740 0.125 0.842 0.652 0.706
#> F1 0.769 0.500 0.269 0.463            
#> F2 0.765 0.494 0.270 0.472            
#> F3 0.745 0.519 0.225 0.408            
```
