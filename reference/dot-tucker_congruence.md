# Tucker congruence between factors

Compute the Tucker congruence matrix between the columns of two loading
matrices.

## Usage

``` r
.tucker_congruence(L1, L2)
```

## Arguments

- L1:

  Numeric matrix.

- L2:

  Numeric matrix with the same dimensions as `L1`.

## Value

A square matrix whose `(i, j)` entry is the Tucker congruence between
column `i` of `L1` and column `j` of `L2`.

## References

Lorenzo-Seva, U., and ten Berge, J. M. F. (2006). Tucker's congruence
coefficient as a meaningful index of factor similarity. *Methodology*,
2, 57-64.
