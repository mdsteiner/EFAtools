// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "eig_utils.h"

using namespace Rcpp;
using namespace arma;

//' Parallel analysis on simulated data.
//'
//' Function called from within efa_parallel() so usually no call to this is needed by the user.
//' Provides a C++ implementation of the efa_parallel() simulation procedure
//'
//' @param n_datasets numeric. Number of datasets with dimensions (N, n_vars) to simulate.
//' @param n_vars numeric. Number of variables / indicators in dataset.
//' @param N numeric. Number of cases / observations in dataset.
//' @param eigen_type numeric. Whether PCA (eigen_type = 1; i.e., leaving diagonal of correlation matrix at 1), SMC (eigen_type = 2; i.e., setting diagonal of correlation matrix to SMCs), or both from the same simulated datasets (eigen_type = 3), in which case the returned matrix holds the PCA eigenvalues in the first n_vars columns and the SMC eigenvalues in the next n_vars.
//' @param maxit numeric. Maximum iterations to perform after which to abort.
// [[Rcpp::export(.parallel_sim)]]
arma::mat parallel_sim(const int n_datasets, const int n_vars, const int N,
                          const int eigen_type, const int maxit = 10000) {
  if (n_datasets < 1) {
    stop("n_datasets must be at least 1.");
  }
  if (n_vars < 1) {
    stop("n_vars must be at least 1.");
  }
  if (N < 2) {
    stop("N must be at least 2.");
  }
  if (eigen_type != 1 && eigen_type != 2 && eigen_type != 3) {
    stop("eigen_type must be 1 (PCA), 2 (SMC), or 3 (both).");
  }
  if (maxit < 1) {
    stop("maxit must be at least 1.");
  }

  const bool want_pca = (eigen_type != 2);
  const bool want_smc = (eigen_type != 1);

  if (want_smc && N <= n_vars) {
    stop("N must be larger than n_vars for SMC simulation.");
  }

  // Initialize needed objects. The simulated data `x` and their correlation matrix
  // `R` keep their dimensions across replicates, so they are allocated once here and
  // refilled in place instead of once per simulated dataset -- at the default
  // n_datasets = 1000 that is 1000 fewer N x n_vars allocations in the busiest
  // simulation loop of the package. x.randn() draws from the same random-number
  // stream in the same (column-major) order as arma::randn(N, n_vars), so the
  // reference eigenvalues are unchanged.
  arma::vec eigval(n_vars);
  arma::mat eig_vals(n_datasets, (want_pca + want_smc) * n_vars);
  arma::mat x(N, n_vars);
  arma::mat R(n_vars, n_vars);
  // Only the SMC series inverts the simulated matrix, so its workspace is left empty
  // otherwise.
  arma::mat temp;
  if (want_smc) {
    temp.set_size(n_vars, n_vars);
  }

  // One loop over the draws for both series: PCA and SMC differ only in the diagonal
  // substituted into the same simulated correlation matrix, so a single draw serves
  // both and the two series stay paired dataset by dataset. A draw the SMC series
  // cannot use -- a simulated matrix with no inverse, hence no squared multiple
  // correlations -- is therefore discarded for the PCA series as well, which keeps
  // that pairing exact. Only the SMC series can reject, so only it needs the maxit
  // bound; the PCA series runs all n_datasets draws whatever maxit says.
  const int smc_col = want_pca ? n_vars : 0;
  int success = 0;
  int iter = 0;

  while ((success < n_datasets) && (!want_smc || (iter < maxit))) {
    iter++;
    x.randn();
    R = cor(x);

    if (want_smc && !inv_sympd(temp, R)) {
      continue;
    }

    if (want_pca) {
      eig_sym_values_checked(eigval, R, "the parallel-analysis simulation");
      eig_vals.submat(success, 0, success, n_vars - 1) = flipud(eigval).t();
    }

    if (want_smc) {
      R.diag() = 1 - (1 / temp.diag());
      eig_sym_values_checked(eigval, R, "the parallel-analysis simulation");
      eig_vals.submat(success, smc_col, success, smc_col + n_vars - 1) =
        flipud(eigval).t();
    }

    success++;
  }

  if (success < n_datasets) {
    stop("Could not generate enough non-singular matrices.");
  }

  return eig_vals;

}
