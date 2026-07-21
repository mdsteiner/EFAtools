# Plot efa_average object

Plot method showing a summarized output of the
[`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
function

## Usage

``` r
# S3 method for class 'efa_average'
plot(x, ...)
```

## Arguments

- x:

  list. An output from the
  [`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
  function.

- ...:

  not used.

## Value

A ggplot object showing, for each indicator and factor, the minimum,
maximum, and average (mean or median) loading across the averaged
solutions.

## Examples

``` r
# \donttest{
EFA_aver <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500)
#>                                                                                                                                                                  🏃 Extracting data...                                                                                                                                                                 🚶 Reordering factors...                                                                                                                                                                 🏃 Averaging data...                                                                                                                                                                 Done!
plot(EFA_aver)

# }
```
