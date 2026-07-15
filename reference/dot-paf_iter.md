# Perform the iterative PAF procedure

Function called from within PAF so usually no call to this is needed by
the user. Provides a C++ implementation of the PAF procedure

## Usage

``` r
.paf_iter(h2, criterion, R, n_fac, abs_eig, crit_type, max_iter)
```

## Arguments

- h2:

  numeric. The initial communality estimates.

- criterion:

  double. The convergence criterion to use.

- R:

  matrix. The correlation matrix with the initial communality estimates
  in the diagonal.

- n_fac:

  numeric. The number of factors to extract.

- abs_eig:

  logical. Whether absolute eigenvalues should be used to compute the
  loadings.

- crit_type:

  numeric. Whether maximum absolute differences (crit_type = 1), or sum
  of differences (crit_type = 2) should be used

- max_iter:

  numeric. The number of iterations after which to end the procedure if
  no convergence has been reached by then.
