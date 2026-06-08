// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Symmetric eigendecompositions that fail loudly instead of silently leaving
// the outputs empty (which the callers would then index out of bounds). During
// ML optimization the decomposition is taken on a rescaled matrix that can go
// non-finite or non-symmetric for some psi, so turn a failed arma::eig_sym into
// a catchable R error rather than undefined behaviour.
static void eig_sym_checked(arma::vec& eigval, arma::mat& eigvec,
                            const arma::mat& X) {
  if (!arma::eig_sym(eigval, eigvec, X)) {
    Rcpp::stop("Eigendecomposition failed during factor extraction; the "
               "correlation matrix is not finite or not symmetric.");
  }
}

static void eig_sym_checked(arma::vec& eigval, const arma::mat& X) {
  if (!arma::eig_sym(eigval, X)) {
    Rcpp::stop("Eigendecomposition failed during factor extraction; the "
               "correlation matrix is not finite or not symmetric.");
  }
}

// [[Rcpp::export(.grad_ml)]]
arma::vec grad_ml(arma::vec psi, arma::mat R, const int n_fac) {
  // gradient function for maximum likelihood estimation, adapted from stats::factanal()

  arma::vec eigval;
  arma::mat eigvec;
  arma::uvec idx(n_fac);
  idx.fill(true);

  arma::mat sc = arma::diagmat(1 / sqrt(psi));
  arma::mat Rs = sc * R * sc;
  eig_sym_checked(eigval, eigvec, Rs);
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
  eig_sym_checked(eigval, Rs);
  Lambda = flipud(eigval);
  int nth = Lambda.n_elem - 1;
  arma::vec e = Lambda.rows(n_fac, nth);
  double out = accu(log(e) - e) - n_fac + R.n_rows;
  return(-out);
}
