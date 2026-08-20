# Parallel analysis on simulated data.

Function called from within efa_parallel() so usually no call to this is
needed by the user. Provides a C++ implementation of the efa_parallel()
simulation procedure

## Usage

``` r
.parallel_sim(n_datasets, n_vars, N, eigen_type, maxit = 10000L)
```

## Arguments

- n_datasets:

  numeric. Number of datasets with dimensions (N, n_vars) to simulate.

- n_vars:

  numeric. Number of variables / indicators in dataset.

- N:

  numeric. Number of cases / observations in dataset.

- eigen_type:

  numeric. Whether PCA (eigen_type = 1; i.e., leaving diagonal of
  correlation matrix at 1), SMC (eigen_type = 2; i.e., setting diagonal
  of correlation matrix to SMCs), or both from the same simulated
  datasets (eigen_type = 3), in which case the returned matrix holds the
  PCA eigenvalues in the first n_vars columns and the SMC eigenvalues in
  the next n_vars.

- maxit:

  numeric. Maximum iterations to perform after which to abort.
