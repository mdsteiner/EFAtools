// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


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
     eig_sym(eigvals, R);
     //Lambda = flipud(eigvals);
     ref_values(i) = eigvals(ind);

   }

   return ref_values;

 }
