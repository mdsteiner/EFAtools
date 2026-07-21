# Print and format an efa_compare object

[`print()`](https://rdrr.io/r/base/print.html) shows a summarised output
of the
[`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)
function: the mean (with its range), median, and root mean squared
distance (RMSE) of the differences, the number of decimals to which all
numbers agree, the minimum number of decimals provided, and (for
matrices) the number of differing indicator-to-factor correspondences,
followed (optionally) by the table of elementwise differences.
[`format()`](https://rdrr.io/r/base/format.html) assembles the same
report and returns it as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled (for example when captured
into a file or stripped with
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html)).

## Usage

``` r
# S3 method for class 'efa_compare'
print(x, ...)

# S3 method for class 'efa_compare'
format(x, ...)
```

## Arguments

- x:

  An object of class `efa_compare` (output from
  [`efa_compare()`](https://mdsteiner.github.io/EFAtools/reference/efa_compare.md)).

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector with the report lines (styled to the active console
theme; plain when colours are disabled).

## Examples

``` r
# A type SPSS EFA to mimick the SPSS implementation
EFA_SPSS_5 <- efa_fit(IDS2_R, n_factors = 5,
                      estimate_control = estimate_control(type = "SPSS"),
                      rotate_control = rotate_control(type = "SPSS"))
#> Warning: Reached the maximum number of iterations without convergence; results may not
#> be interpretable.

# A type psych EFA to mimick the psych::fa() implementation
EFA_psych_5 <- efa_fit(IDS2_R, n_factors = 5,
                       estimate_control = estimate_control(type = "psych"),
                       rotate_control = rotate_control(type = "psych"))

# compare the two
comp <- efa_compare(EFA_SPSS_5$unrot_loadings, EFA_psych_5$unrot_loadings,
                    x_labels = c("SPSS", "psych"))
comp
#> Mean [min, max] absolute difference:  0.0017 [ 0.0000,  0.0090]
#> Median absolute difference:  0.0009
#> Root mean squared distance (RMSE):  0.0025
#> Max decimals where all numbers agree in absolute value: 1
#> Minimum number of decimals provided: 17
#> Differing indicator-to-factor correspondences: 0 (highest loading), 0 (all |loadings| >= 0.3)
#> 
#>        F1      F2      F3      F4      F5
#> GS    .0004   .0002  -.0002   .0002  -.0034
#> PL    .0000   .0001   .0007  -.0003   .0001
#> TC    .0007   .0021  -.0058   .0001   .0024
#> CB   -.0007  -.0039  -.0005   .0043  -.0024
#> NL    .0020  -.0090  -.0056   .0032   .0016
#> NLM  -.0018   .0069   .0082   .0034  -.0011
#> GF   -.0002   .0003   .0036  -.0021  -.0003
#> RGF   .0004   .0011   .0049   .0020   .0011
#> CM   -.0001  -.0001   .0004  -.0009   .0011
#> EP   -.0001  -.0001   .0003  -.0010   .0009
#> CA   -.0001   .0006  -.0020  -.0027   .0009
#> OP   -.0001   .0007  -.0016  -.0025   .0008
#> RS   -.0001   .0012  -.0013  -.0032  -.0010
#> DP   -.0001   .0009   .0004  -.0024  -.0009

# format() returns the same lines as plain text:
writeLines(format(comp))
#> Mean [min, max] absolute difference:  0.0017 [ 0.0000,  0.0090]
#> Median absolute difference:  0.0009
#> Root mean squared distance (RMSE):  0.0025
#> Max decimals where all numbers agree in absolute value: 1
#> Minimum number of decimals provided: 17
#> Differing indicator-to-factor correspondences: 0 (highest loading), 0 (all |loadings| >= 0.3)
#> 
#>        F1      F2      F3      F4      F5
#> GS    .0004   .0002  -.0002   .0002  -.0034
#> PL    .0000   .0001   .0007  -.0003   .0001
#> TC    .0007   .0021  -.0058   .0001   .0024
#> CB   -.0007  -.0039  -.0005   .0043  -.0024
#> NL    .0020  -.0090  -.0056   .0032   .0016
#> NLM  -.0018   .0069   .0082   .0034  -.0011
#> GF   -.0002   .0003   .0036  -.0021  -.0003
#> RGF   .0004   .0011   .0049   .0020   .0011
#> CM   -.0001  -.0001   .0004  -.0009   .0011
#> EP   -.0001  -.0001   .0003  -.0010   .0009
#> CA   -.0001   .0006  -.0020  -.0027   .0009
#> OP   -.0001   .0007  -.0016  -.0025   .0008
#> RS   -.0001   .0012  -.0013  -.0032  -.0010
#> DP   -.0001   .0009   .0004  -.0024  -.0009
```
