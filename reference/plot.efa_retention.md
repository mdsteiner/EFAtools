# Plot method for efa_retention objects

Plots the result of a factor-retention criterion. Eigenvalue-based
criteria (e.g.
[`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md))
are shown as an eigenvalue plot, the Hull method
([`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md))
as a convex-hull plot. Criteria with more than one sub-variant are
faceted.

## Usage

``` r
# S3 method for class 'efa_retention'
plot(x, ...)
```

## Arguments

- x:

  an object of class efa_retention, returned by a factor-retention
  criterion (e.g.
  [`efa_ekc()`](https://mdsteiner.github.io/EFAtools/reference/efa_ekc.md)
  or
  [`efa_hull()`](https://mdsteiner.github.io/EFAtools/reference/efa_hull.md)).

- ...:

  not used.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object, or invisibly `NULL` if the criterion has no plottable result.

## Examples

``` r
plot(efa_ekc(test_models$baseline$cormat, N = 500))
```
