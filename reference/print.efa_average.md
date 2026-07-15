# Print and format an efa_average object

[`print()`](https://rdrr.io/r/base/print.html) shows a summarised output
of the
[`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)
function: the averaging settings, the
error/convergence/Heywood/admissibility rates, the indicator-to-factor
correspondences, the averaged loadings (and, for oblique solutions, the
factor intercorrelations), the variances accounted for, and the model
fit. [`format()`](https://rdrr.io/r/base/format.html) assembles the same
report and returns it as a character vector;
[`print()`](https://rdrr.io/r/base/print.html) is
`cat(format(x), sep = "\n")`. The lines follow the active console theme,
so they are plain when colours are disabled (for example when captured
into a file or stripped with
[`cli::ansi_strip()`](https://cli.r-lib.org/reference/ansi_strip.html)).

## Usage

``` r
# S3 method for class 'efa_average'
print(x, stat = c("average", "range"), plot = FALSE, ...)

# S3 method for class 'efa_average'
format(x, stat = c("average", "range"), ...)
```

## Arguments

- x:

  An object of class `efa_average` (output from
  [`efa_average()`](https://mdsteiner.github.io/EFAtools/reference/efa_average.md)).

- stat:

  character. A vector with the statistics to print. Possible inputs are
  "average", "sd", "range", "min", and "max". Default is "average" and
  "range".

- plot:

  logical. Whether a plot of the average and min- max loadings should be
  created. Default is FALSE. If more than 10 factors are extracted, no
  plot is created. Only used by
  [`print()`](https://rdrr.io/r/base/print.html).

- ...:

  Not used; for consistency with the generic.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument `x`
invisibly. [`format()`](https://rdrr.io/r/base/format.html) returns a
character vector with the report lines (styled to the active console
theme; plain when colours are disabled).

## Examples

``` r
if (FALSE) { # \dontrun{
EFA_aver <- efa_average(test_models$baseline$cormat, n_factors = 3, N = 500)
EFA_aver

# format() returns the same lines as plain text:
writeLines(format(EFA_aver))
} # }
```
