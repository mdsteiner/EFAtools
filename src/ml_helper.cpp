// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export(.grad_ml)]]
arma::vec grad_ml(arma::vec psi, arma::mat R, const int n_fac) {
  // gradient function for maximum likelihood estimation, adapted from stats::factanal()

  arma::vec eigval;
  arma::mat eigvec;
  arma::uvec idx(n_fac);
  idx.fill(true);

  arma::mat sc = arma::diagmat(1 / sqrt(psi));
  arma::mat Rs = sc * R * sc;
  eig_sym(eigval, eigvec, Rs);
  arma::vec Lambda = flipud(eigval);
  Lambda = Lambda.elem(arma::find(idx));
  Lambda -= 1;
  // replace values smaller than 0
  arma::uvec idx2 = find(Lambda < 0);
  Lambda.elem(idx2).fill(0);
  // extract eigenvectors
  arma::mat V = fliplr(eigvec);
  V = V.cols(arma::find(idx));

  arma::mat load = V * arma::diagmat(sqrt(Lambda));
  load = arma::diagmat(sqrt(psi)) * load;
  arma::mat g = load * load.t() + diagmat(psi) - R;
  arma::vec out = g.diag();
  out = out / pow(psi, 2);

  return(out);
}


// [[Rcpp::export(.error_ml)]]
double error_ml(arma::vec psi, arma::mat R, const int n_fac) {
  // loss function for maximum likelihood fitting; adapted from stats::factanal()

  arma::vec eigval;
  arma::vec Lambda;

  arma::mat sc = arma::diagmat(1 / sqrt(psi));
  arma::mat Rs = sc * R * sc;
  eig_sym(eigval, Rs);
  Lambda = flipud(eigval);
  int nth = Lambda.n_elem - 1;
  arma::vec e = Lambda.rows(n_fac, nth);
  double out = accu(log(e) - e) - n_fac + R.n_rows;
  return(-out);
}
