# Closed-form orthogonal Procrustes rotation

Rotate `A` to the orthogonal target `B` by minimizing
`||A %*% T - B||_F^2` subject to `t(T) %*% T = I`.

## Usage

``` r
.orthogonal_procrustes(A, B)
```

## Arguments

- A:

  Numeric matrix to be rotated.

- B:

  Numeric target matrix with the same dimensions as `A`.

## Value

A list with the rotated loadings, orthogonal transformation matrix,
target criterion value, and basic diagnostics.

## References

Schoenemann, P. H. (1966). A generalized solution of the orthogonal
Procrustes problem. *Psychometrika*, 31, 1-10.
