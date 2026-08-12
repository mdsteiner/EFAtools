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
//' @param eigen_type numeric. Whether PCA (eigen_type = 1; i.e., leaving diagonal of correlation matrix at 1) or SMC (eigen_type = 2; i.e., setting diagonal of correlation matrix to SMCs).
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
  if (eigen_type != 1 && eigen_type != 2) {
    stop("eigen_type must be 1 (PCA) or 2 (SMC).");
  }
  if (maxit < 1) {
    stop("maxit must be at least 1.");
  }
  if (eigen_type == 2 && N <= n_vars) {
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
  arma::mat eig_vals(n_datasets, n_vars);
  arma::mat x(N, n_vars);
  arma::mat R(n_vars, n_vars);

  if (eigen_type == 1) { // PCA

    // perform simulations for n_datasets time
    for (int i = 0; i < n_datasets; i++) {
      x.randn();
      R = cor(x);
      eig_sym_values_checked(eigval, R, "the parallel-analysis simulation");
      eig_vals.row(i) = flipud(eigval).t();
    }

  } else if (eigen_type == 2) { // SMC
    arma::vec smc(n_vars);
    arma::mat temp(n_vars, n_vars);
    int success = 0;
    int iter = 0;

    while((success < n_datasets) && (iter < maxit)) {
      x.randn();
      R = cor(x);
      bool flag = inv_sympd(temp, R);
      if (!flag){
        iter++;
        continue;
      }
      R.diag() = 1 - (1 / temp.diag());
      eig_sym_values_checked(eigval, R, "the parallel-analysis simulation");
      eig_vals.row(success) = flipud(eigval).t();
      iter++;
      success++;
    }

    if ((iter == maxit) && (success < n_datasets)) {
      stop("Could not generate enough non-singular matrices.");
    }

  }

  return eig_vals;

}
