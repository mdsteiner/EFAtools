# Print and format a reliability object

[`print()`](https://rdrr.io/r/base/print.html) shows the reliability
coefficients (omega total, hierarchical, and subscale, standardized
alpha, and the H index) of each factor and the common-variance indices
(ECV and PUC) of the general factor, for one solution or for each group.
[`format()`](https://rdrr.io/r/base/format.html) returns the same lines
as a character vector. The lines follow the active console theme, and
are plain when colours are disabled.

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
character vector with the report lines.

## See also

Other reliability coefficients:
[`efa_reliability()`](https://mdsteiner.github.io/EFAtools/reference/efa_reliability.md),
[`efa_schmid_leiman()`](https://mdsteiner.github.io/EFAtools/reference/efa_schmid_leiman.md)

## Examples

``` r
efa_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                   estimator = "PAF", rotation = "promax")
rel <- efa_reliability(efa_mod)
rel
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>         tot  alpha    H
#> total  .883   .868
#> F1     .769   .768  .760
#> F2     .765   .763  .753
#> F3     .745   .743  .738

# format() returns the same lines as a character vector:
writeLines(format(rel))
#> 
#> ── Reliability coefficients ────────────────────────────────────────────────────
#> 
#>         tot  alpha    H
#> total  .883   .868
#> F1     .769   .768  .760
#> F2     .765   .763  .753
#> F3     .745   .743  .738
```
