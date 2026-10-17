# Parallel analysis on simulated correlation matrices.

Function called from within efa_parallel() so usually no call to this is
needed by the user. Provides a C++ implementation of the efa_parallel()
simulation procedure for Pearson correlations. Each simulated
correlation matrix is drawn directly from the sampling distribution of
the Pearson correlation matrix of `N` cases of `n_vars` uncorrelated
normal variables: a Wishart matrix with identity covariance and `N - 1`
degrees of freedom, rescaled to correlations. This is exact and does not
simulate the data.

## Usage

``` r
.parallel_sim(n_datasets, n_vars, N, eigen_type, maxit = 10000L)
```

## Arguments

- n_datasets:

  numeric. Number of correlation matrices to draw.

- n_vars:

  numeric. Number of variables / indicators in dataset.

- N:

  numeric. Number of cases / observations in dataset. Must be larger
  than `n_vars`.

- eigen_type:

  numeric. Whether PCA (eigen_type = 1; i.e., leaving diagonal of
  correlation matrix at 1), SMC (eigen_type = 2; i.e., setting diagonal
  of correlation matrix to SMCs), or both from the same draws
  (eigen_type = 3), in which case the returned matrix holds the PCA
  eigenvalues in the first n_vars columns and the SMC eigenvalues in the
  next n_vars.

- maxit:

  numeric. Maximum iterations to perform after which to abort.
