# Compute explained variances from loadings

From unrotated loadings compute the communalities and uniquenesses for
total variance. Compute explained variances per factor from rotated
loadings (and factor intercorrelations Phi if oblique rotation was
used).

## Usage

``` r
.compute_vars(L_unrot, L_rot, Phi = NULL)
```

## Arguments

- L_unrot:

  matrix. Unrotated factor loadings.

- L_rot:

  matrix. Rotated factor loadings.

- Phi:

  matrix. Factor intercorrelations. Provide only if oblique rotation is
  used.

## Value

A matrix with sum of squared loadings, proportion explained variance
from total variance per factor, same as previous but cumulative,
Proportion of explained variance from total explained variance, and same
as previous but cumulative. The three cumulative and common-variance
rows are omitted when `L_rot` has a single column, where they would only
repeat the two above them, so the result has two rows there and five
otherwise.
