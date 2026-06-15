// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Symmetric eigendecomposition that fails loudly instead of silently leaving
// eigval/eigvec empty (which the callers would then index out of bounds). During
// ULS optimization the decomposition is taken on R - diag(psi), which can go
// non-finite or non-symmetric for some psi, so turn a failed arma::eig_sym into
// a catchable R error rather than undefined behaviour.
static void eig_sym_checked(arma::vec& eigval, arma::mat& eigvec,
                            const arma::mat& X) {
  if (!arma::eig_sym(eigval, eigvec, X)) {
    Rcpp::stop("Eigendecomposition failed during factor extraction; the "
               "correlation matrix is not finite or not symmetric.");
  }
}

// [[Rcpp::export(.grad_uls)]]
arma::vec grad_uls(arma::vec psi, const arma::mat& R, const int n_fac) {

  // Guard the eigen-extraction below: it reads the largest n_fac eigenpairs and
  // would index past the available eigenvalues if n_fac >= ncol(R) (undefined
  // behaviour / crash in an unchecked build). Convert that into a catchable error.
  if (n_fac < 1 || static_cast<arma::uword>(n_fac) >= R.n_cols) {
    Rcpp::stop("n_fac must be at least 1 and smaller than the number of "
               "variables (got n_fac = %d for %d variables).",
               n_fac, static_cast<int>(R.n_cols));
  }

  arma::vec eigval;
  arma::mat eigvec;
  arma::mat load;
  arma::mat g;
  arma::vec Lambda;
  arma::mat V;

  arma::mat Rs = R - arma::diagmat(psi);
  eig_sym_checked(eigval, eigvec, Rs);
  Lambda = flipud(eigval);
  Lambda = Lambda.head(n_fac);
  // replace values smaller than 0
  arma::uvec idx2 = find(Lambda < 0);
  Lambda.elem(idx2).fill(0);
  // extract eigenvectors
  V = fliplr(eigvec);
  V = V.head_cols(n_fac);

  load = V * arma::diagmat(sqrt(Lambda));
  g = load * load.t() + diagmat(psi) - R;
  arma::vec out = g.diag();

  return(out);

}


// [[Rcpp::export(.uls_residuals)]]
double uls_residuals(arma::vec psi, arma::mat R, const int n_fac) {

  // Guard the eigen-extraction below: it reads the largest n_fac eigenpairs and
  // would index past the available eigenvalues if n_fac >= ncol(R) (undefined
  // behaviour / crash in an unchecked build). Convert that into a catchable error.
  if (n_fac < 1 || static_cast<arma::uword>(n_fac) >= R.n_cols) {
    Rcpp::stop("n_fac must be at least 1 and smaller than the number of "
               "variables (got n_fac = %d for %d variables).",
               n_fac, static_cast<int>(R.n_cols));
  }

  arma::vec eigval;
  arma::mat eigvec;
  arma::mat loadings;
  arma::mat model;
  arma::mat residual;
  arma::vec Lambda;
  arma::mat V;
  arma::vec tv(R.n_cols);
  double error;


  R.diag() = 1 - psi;
  eig_sym_checked(eigval, eigvec, R);
  Lambda = flipud(eigval);
  Lambda = Lambda.head(n_fac);
  // replace values smaller than 0
  arma::uvec idx2 = find(Lambda <   datum::eps);
  Lambda.elem(idx2).fill(datum::eps * 100);
  // extract eigenvectors
  V = fliplr(eigvec);
  V = V.head_cols(n_fac);

  if (n_fac > 1) {
    loadings = V * arma::diagmat(sqrt(Lambda));
  } else {
    Lambda = arma::sqrt(Lambda);
    tv.fill(Lambda[0]);
    loadings = V % tv;
  }

  model = loadings * loadings.t();
  residual = R - model;
  // Sum the squared off-diagonal residuals only (strictly lower triangle). The
  // diagonal residual is the communality gap, which the free uniquenesses absorb
  // and which the analytic gradient .grad_uls() and the reported Fm both exclude;
  // keeping the objective off-diagonal makes all three consistent (minres).
  residual = arma::trimatl(residual, -1);
  residual = pow(residual, 2);
  error = accu(residual);
  return(error);

}
