// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "eig_utils.h"

using namespace Rcpp;
using namespace arma;

// Wishart draw of a simulated correlation matrix; defined and explained in sim.cpp.
void wishart_factor(arma::mat& U, const int n);
void cov2cor_inplace(arma::mat& W);

//' Parallel analysis on simulated correlation matrices.
//'
//' Function called from within efa_parallel() so usually no call to this is needed by the user.
//' Provides a C++ implementation of the efa_parallel() simulation procedure for Pearson
//' correlations. Each simulated correlation matrix is drawn directly from the sampling
//' distribution of the Pearson correlation matrix of `N` cases of `n_vars` uncorrelated
//' normal variables: a Wishart matrix with identity covariance and `N - 1` degrees of
//' freedom, rescaled to correlations. This is exact and does not simulate the data.
//'
//' @param n_datasets numeric. Number of correlation matrices to draw.
//' @param n_vars numeric. Number of variables / indicators in dataset.
//' @param N numeric. Number of cases / observations in dataset. Must be larger than `n_vars`.
//' @param eigen_type numeric. Whether PCA (eigen_type = 1; i.e., leaving diagonal of correlation matrix at 1), SMC (eigen_type = 2; i.e., setting diagonal of correlation matrix to SMCs), or both from the same draws (eigen_type = 3), in which case the returned matrix holds the PCA eigenvalues in the first n_vars columns and the SMC eigenvalues in the next n_vars.
//' @param maxit numeric. Maximum iterations to perform after which to abort.
//'
//' @keywords internal
// [[Rcpp::export(.parallel_sim)]]
arma::mat parallel_sim(const int n_datasets, const int n_vars, const int N,
                          const int eigen_type, const int maxit = 10000) {
  if (n_datasets < 1) {
    stop("n_datasets must be at least 1.");
  }
  if (n_vars < 1) {
    stop("n_vars must be at least 1.");
  }
  if (N <= n_vars) {
    stop("N must be larger than n_vars.");
  }
  if (eigen_type != 1 && eigen_type != 2 && eigen_type != 3) {
    stop("eigen_type must be 1 (PCA), 2 (SMC), or 3 (both).");
  }
  if (maxit < 1) {
    stop("maxit must be at least 1.");
  }

  const bool want_pca = (eigen_type != 2);
  const bool want_smc = (eigen_type != 1);

  // Initialize needed objects. The Wishart factor `U` and the simulated correlation
  // matrix `R` keep their dimensions across replicates, so they are allocated once here
  // and refilled in place instead of once per draw. Only the upper triangle of `U` is
  // redrawn; its strict lower triangle stays zero.
  arma::vec eigval(n_vars);
  arma::mat eig_vals(n_datasets, (want_pca + want_smc) * n_vars);
  arma::mat U(n_vars, n_vars, arma::fill::zeros);
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
    wishart_factor(U, N - 1);
    R = U.t() * U;
    cov2cor_inplace(R);

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
