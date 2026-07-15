# Plot method for efa_retain objects

Plots every factor-retention criterion in the
[`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md)
result that has a plottable outcome (see
[`plot.efa_retention()`](https://mdsteiner.github.io/EFAtools/reference/plot.efa_retention.md));
criteria without a plot (e.g.
[`efa_map()`](https://mdsteiner.github.io/EFAtools/reference/efa_map.md)
or
[`efa_smt()`](https://mdsteiner.github.io/EFAtools/reference/efa_smt.md))
are skipped.

## Usage

``` r
# S3 method for class 'efa_retain'
plot(x, ...)
```

## Arguments

- x:

  an object of class efa_retain, returned by
  [`efa_retain()`](https://mdsteiner.github.io/EFAtools/reference/efa_retain.md).

- ...:

  not used.

## Value

A named list of
[ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
objects, one per criterion with a plottable result, or invisibly `NULL`
if there is none.

## Examples

``` r
# \donttest{
nf <- efa_retain(test_models$baseline$cormat, criteria = c("EKC", "SMT"),
                 N = 500)
plot(nf)
#> $EKC

#> 
# }
```
