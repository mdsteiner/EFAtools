# Print and format an efa_simulated object

[`print()`](https://rdrr.io/r/base/print.html) shows a compact summary
of the data simulated by
[`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md):
how many datasets were drawn and their dimensions, the marginal
distribution, whether the data were discretized into ordered categories
or given missing values, and – when model error was injected – the
method with the target and achieved RMSEA and CFI. The simulated data
themselves live in the `data` element and the population correlation
matrix in `population`. [`format()`](https://rdrr.io/r/base/format.html)
returns the same summary as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled.

## Usage

``` r
# S3 method for class 'efa_simulated'
print(x, digits = 3, ...)

# S3 method for class 'efa_simulated'
format(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `efa_simulated` (output from
  [`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)).

- digits:

  Integer. The number of decimal places the reported fit values are
  rounded to. Default is 3.

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector with the summary lines (styled to the active console
theme; plain when colours are disabled).

## See also

Other data simulation:
[`efa_simulate()`](https://mdsteiner.github.io/EFAtools/reference/efa_simulate.md)

## Examples

``` r
Lambda <- population_models$loadings$baseline
Phi <- population_models$phis_3$moderate
efa_simulate(N = 500, Lambda = Lambda, Phi = Phi, target_rmsea = 0.05, seed = 42)
#> Simulated data (<efa_simulated>)
#> 
#> 500 cases by 18 variables (continuous).
#> Marginals: normal.
#> 
#> Model error: CB
#> RMSEA: target .050, achieved .050
#> CFI: target --, achieved .936
```
