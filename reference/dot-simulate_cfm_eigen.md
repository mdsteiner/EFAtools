# Reference eigenvalues for the efa_nest() simulation.

Internal helper called from efa_nest(). Draws `nreps` correlation
matrices with the sampling distribution of the Pearson correlation
matrix of `N` cases from the `(nf - 1)`-factor reference model, given
that model's loadings `Lambda` and uniquenesses `Psi`, and returns the
`nf`-th largest eigenvalue of each. The matrices are drawn from the
Wishart distribution with covariance `Lambda Lambda' + diag(Psi)` and
`N - 1` degrees of freedom and rescaled to correlations, which is exact
for normal data and does not simulate the data.

## Usage

``` r
.simulate_cfm_eigen(nf, N, Lambda, Psi, nreps = 1000L)
```

## Arguments

- nf:

  integer. Position of the empirical eigenvalue being tested (1-based);
  the `nf`-th largest simulated eigenvalue is returned per replicate.

- N:

  integer. Number of cases / observations per simulated dataset. Must be
  larger than the number of variables.

- Lambda:

  numeric matrix. Loadings of the `(nf - 1)`-factor reference model
  (`p x (nf - 1)`); pass a `p x 0` matrix for the `nf == 1` null
  (identity) model.

- Psi:

  numeric vector. Uniquenesses (`1 - h2`) of the reference model.

- nreps:

  integer. Number of correlation matrices to draw.
