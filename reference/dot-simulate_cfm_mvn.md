# Draw multivariate-normal data from a population correlation matrix.

Internal helper called from efa_simulate() and from the rank-correlation
reference simulations of efa_parallel(). Draws `N` cases from a
`p`-variate normal with correlation (or covariance) `R` by
post-multiplying a matrix of standard normal deviates by a matrix square
root `M` of `R` (with M' M = R, so the rows of Z \* M are N(0, R)). A
positive-definite `R` is factored by Cholesky; a positive-semidefinite
but singular `R` (which makes the Cholesky fail although it is still a
valid covariance, e.g. a no-factor block or a smoothed factor
intercorrelation matrix) falls back to a symmetric eigen square root.

## Usage

``` r
.simulate_cfm_mvn(R, N, tol = 1e-08)
```

## Arguments

- R:

  numeric matrix. Population correlation/covariance matrix.

- N:

  integer. Number of cases to draw.

- tol:

  numeric. Eigenvalues below `-tol` mark `R` as indefinite.
