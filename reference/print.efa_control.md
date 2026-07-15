# Print and format a control object

[`print()`](https://rdrr.io/r/base/print.html) shows the chosen `type`
and each tuning knob, with an unset (`NA`) preset-driven knob marked as
resolved from the `type` preset.
[`format()`](https://rdrr.io/r/base/format.html) assembles the same
report and returns it as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled.

## Usage

``` r
# S3 method for class 'efa_estimate_control'
print(x, ...)

# S3 method for class 'efa_estimate_control'
format(x, ...)

# S3 method for class 'efa_rotate_control'
print(x, ...)

# S3 method for class 'efa_rotate_control'
format(x, ...)
```

## Arguments

- x:

  A control object from
  [`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)
  or
  [`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md).

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector with the report lines.

## See also

[`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md),
[`rotate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)

Other Control functions:
[`estimate_control()`](https://mdsteiner.github.io/EFAtools/reference/estimate_control.md)

## Examples

``` r
est <- estimate_control(type = "SPSS")
est
#> Estimation control (type: "SPSS")
#> 
#> init_comm: <from type preset>
#> criterion: <from type preset>
#> criterion_type: <from type preset>
#> max_iter: <from type preset>
#> abs_eigen: <from type preset>
#> start_method: psych
writeLines(format(est))
#> Estimation control (type: "SPSS")
#> 
#> init_comm: <from type preset>
#> criterion: <from type preset>
#> criterion_type: <from type preset>
#> max_iter: <from type preset>
#> abs_eigen: <from type preset>
#> start_method: psych
```
