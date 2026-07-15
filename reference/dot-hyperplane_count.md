# Count near-zero loadings

Hyperplane count is the number of loadings with absolute value smaller
than a user-specified cutoff.

## Usage

``` r
.hyperplane_count(L, cutoff = 0.15)
```

## Arguments

- L:

  Numeric loading matrix.

- cutoff:

  Numeric scalar. Loadings with `abs(L) < cutoff` are counted as being
  in the hyperplane.

## Value

A list with the total hyperplane count and counts by factor and item.
