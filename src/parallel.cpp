// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Symmetric eigendecomposition that fails loudly instead of leaving the output
// empty (which the caller would then index out of bounds). A degenerate simulated
// matrix can make arma::eig_sym return false, so turn that into a catchable R error.
static void eig_sym_checked(arma::vec& eigval, const arma::mat& X) {
  if (!arma::eig_sym(eigval, X)) {
    Rcpp::stop("Eigendecomposition failed during the parallel-analysis simulation; "
               "the simulated correlation matrix is not finite or not symmetric.");
  }
}

//' Parallel analysis on simulated data.
//'
//' Function called from within PARALLEL so usually no call to this is needed by the user.
//' Provides a C++ implementation of the PARALLEL simulation procedure
//'
//' @param n_datasets numeric. Number of datasets with dimensions (N, n_vars) to simulate.
//' @param n_vars numeric. Number of variables / indicators in dataset.
//' @param N numeric. Number of cases / observations in dataset.
//' @param eigen_type numeric. Whether PCA (eigen_type = 1; i.e., leaving diagonal of correlation matrix at 1) or SMC (eigen_type = 2; i.e., setting diagonal of correlation matrix to SMCs).
//' @param maxit numeric. Maximum iterations to perform after which to abort.
// [[Rcpp::export(.parallel_sim)]]
arma::mat parallel_sim(const int n_datasets, const int n_vars, const int N,
                         const int eigen_type, const int maxit = 10000) {
  // initialize needed objects
  arma::vec Lambda(n_vars);
  arma::vec eigval(n_vars);
  arma::mat eig_vals(n_datasets, n_vars);

  if (eigen_type == 1) { // PCA

    // perform simulations for n_datasets time
    for (int i = 0; i < n_datasets; i++) {
      arma::mat x = randn(N, n_vars);
      arma::mat R = cor(x);
      eig_sym_checked(eigval, R);
      Lambda = flipud(eigval);
      eig_vals.row(i) = Lambda.t();
    }

  } else if (eigen_type == 2) { // SMC
    arma::vec smc(n_vars);
    arma::mat temp(n_vars, n_vars);
    int success = 0;
    int iter = 0;

    while((success < n_datasets) && (iter < maxit)) {
      arma::mat x = randn(N, n_vars);
      arma::mat R = cor(x);
      bool flag = inv_sympd(temp, R);
      if (!flag){
        iter++;
        continue;
      }
      R.diag() = 1 - (1 / temp.diag());
      eig_sym_checked(eigval, R);
      Lambda = flipud(eigval);
      eig_vals.row(success) = Lambda.t();
      iter++;
      success++;
    }

    if ((iter == maxit) && (success < n_datasets)) {
      stop("Could not generate enough non-singular matrices.");
    }

  }

  return eig_vals;

}
