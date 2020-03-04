// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Perform the iterative PAF procedure
//'
//' Function called from within PAF so usually no call to this is needed by the user.
//' Provides a C++ implementation of the PAF procedure
//'
//' @param n_fac numeric. The number of factors to extract.
//' @param R matrix. The correlation matrix.
//' @param criterion double. The convergence criterion to use.
//' @param max_iter numeric. The number of iterations after which to end the procedure if no convergence has been reached by then.
//' @export
// [[Rcpp::export]]
arma::mat hull_paf(const int n_fac, arma::mat R, double criterion, int max_iter) {

  int iter = 1;
  double delta = 1.0;
  arma::vec tv(R.n_cols);
  arma::vec Lambda;
  arma::mat V;
  arma::mat L;
  arma::vec new_h2;
  arma::vec eigval;
  arma::mat eigvec;
  arma::mat Lt;
  arma::vec h2;
  arma::uvec idx(n_fac);
  idx.fill(true);

  // compute smcs
  arma::mat temp(R.n_cols, R.n_cols);
  temp = inv_sympd(R);
  R.diag() = 1 - (1 / temp.diag());
  h2 = temp.diag();

    if (n_fac > 1) {

        while (delta > criterion & iter <= max_iter) {
          //  compute the eigenvalues and eigenvectors
          eig_sym(eigval, eigvec, R);
          Lambda = arma::abs(flipud(eigval));
          Lambda = Lambda.elem(arma::find(idx));
          V = fliplr(eigvec);
          V = V.cols(arma::find(idx));

          // compute the loadings from the eigenvector matrix and diagonal
          // eigenvalue matrix
          L = V * arma::diagmat(sqrt(arma::abs(Lambda)));

          // get the new communality estimates from the loadings
          Lt = L * L.t();
          new_h2 = Lt.diag();

          // convergence criterion according to the psych package
          delta = std::abs(arma::accu(h2) - arma::accu(new_h2));

          // update diagonal of R with the new communality estimates
          R.diag() = new_h2;

          // update old communality estimates with new ones
          h2 = new_h2;

          // incerase iterator
          iter += 1;

        }

    } else {

        while (delta > criterion & iter <= max_iter) {
          //  compute the eigenvalues and eigenvectors
          eig_sym(eigval, eigvec, R);
          Lambda = arma::abs(flipud(eigval));
          Lambda = Lambda.elem(arma::find(idx));
          V = fliplr(eigvec);
          V = V.cols(arma::find(idx));

          // compute the loadings from the eigenvector matrix and diagonal
          // eigenvalue matrix
          Lambda = arma::sqrt(arma::abs(Lambda));
          tv.fill(Lambda[0]);
          L = V % tv;

          // get the new communality estimates from the loadings
          Lt = L * L.t();
          new_h2 = Lt.diag();

          // convergence criterion according to the psych package
          delta = std::abs(arma::accu(h2) - arma::accu(new_h2));

          // update diagonal of R with the new communality estimates
          R.diag() = new_h2;

          // update old communality estimates with new ones
          h2 = new_h2;

          // incerase iterator
          iter += 1;

        }

    }

  // break if after maximum iterations there was no convergence
  if (iter >= max_iter){
    warning("Reached maximum number of iterations without convergence. Results may not be interpretable.");
  }

  return L;
}

