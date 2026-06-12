// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Symmetric eigendecomposition that fails loudly instead of leaving the output
// empty (which the caller would then index out of bounds). A degenerate simulated
// matrix can make arma::eig_sym return false, so turn that into a catchable R error.
static void eig_sym_checked(arma::vec& eigval, const arma::mat& X) {
  if (!arma::eig_sym(eigval, X)) {
    Rcpp::stop("Eigendecomposition failed during the NEST simulation; the "
               "simulated correlation matrix is not finite or not symmetric.");
  }
}


//' Get reference values for nest.
//'
//' Function called from within NEST so usually no call to this is needed by the user.
//' Provides a C++ implementation of the NEST simulation procedure
//'
//' @param nf numeric. Current number of factors.
//' @param N numeric. Number of cases / observations in dataset.
//' @param M numeric matrix. A transposed matrix of cbind(loadings, diag(sqrt(uniquenesses))).
//' @param nreps numeric. Number of datasets to create.
// [[Rcpp::export(.nest_sym)]]
arma::vec nest_sym(const int nf, const int N, arma::mat M,
                    const int nreps = 1000) {
   // nf indexes into the eigenvalues via ind = ncol - nf; a too-large nf would
   // make ind negative and read out of bounds, so reject it up front.
   if (nf < 1 || static_cast<arma::uword>(nf) > M.n_cols) {
     Rcpp::stop("nf must be between 1 and the number of variables.");
   }

   // initialize needed objects
   const int ncm = M.n_rows;
   const int ind = M.n_cols - nf;
   arma::vec eigvals(M.n_cols);
   //arma::vec Lambda(M.n_cols);
   arma::vec ref_values(nreps);
   arma::mat dat(N, M.n_cols);
   arma::mat R(M.n_cols, M.n_cols);

   for (uword i = 0; i < nreps; i++) {

     // generate random data for nfactors + nvariables colums (ncm)
     dat = randn(N, ncm);
     R = cor(dat * M);
     eig_sym_checked(eigvals, R);
     //Lambda = flipud(eigvals);
     ref_values(i) = eigvals(ind);

   }

   return ref_values;

 }
