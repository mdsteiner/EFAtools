# Reference eigenvalues for the efa_nest() simulation via the shared kernel.

Internal helper called from efa_nest(). Simulates `nreps` datasets from
an `(nf - 1)`-factor reference model, given that model's loadings
`Lambda` and uniquenesses `Psi`, and returns the `nf`-th largest
eigenvalue of each simulated correlation matrix. The data are drawn with
the shared Z \* M rule (see .simulate_cfm_mvn()) using the factor-score
square root `M = t([Lambda | diag(sqrt(Psi))])`, so a row
`randn(1, nf - 1 + p) * M` is N(0, Lambda Lambda' + diag(Psi)). Drawing
`nf - 1 + p` standard normals and post-multiplying by the factor-score
matrix is faster than forming the model- implied correlation matrix and
drawing from it, and matches the position at which efa_nest() reads the
reference eigenvalue.

## Usage

``` r
.simulate_cfm_eigen(nf, N, Lambda, Psi, nreps = 1000L)
```

## Arguments

- nf:

  integer. Position of the empirical eigenvalue being tested (1-based);
  the `nf`-th largest simulated eigenvalue is returned per replicate.

- N:

  integer. Number of cases / observations per simulated dataset.

- Lambda:

  numeric matrix. Loadings of the `(nf - 1)`-factor reference model
  (`p x (nf - 1)`); pass a `p x 0` matrix for the `nf == 1` null
  (identity) model.

- Psi:

  numeric vector. Uniquenesses (`1 - h2`) of the reference model.

- nreps:

  integer. Number of datasets to simulate.
