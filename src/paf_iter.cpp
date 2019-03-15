// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Perform the iterative PAF procedure
//'
//' Function called from within PAF so usually no call to this is needed by the user.
//' Provides a C++ implementation of the PAF procedure
//'
//' @param h2 numeric. The initial communality estimates.
//' @param criterion double. The convergence criterion to use.
//' @param R matrix. The correlation matrix with the initial communality estimates in the diagonal.
//' @param n_fac numeric. The number of factors to extract.
//' @param abs_eig logical. Whether absolute eigenvalues should be used to compute the loadings.
//' @param crit_type numeric. Whether maximum absolute differences (crit_type = 1), or sum of differences (crit_type = 2) should be used
//' @param max_iter numeric. The number of iterations after which to end the procedure if no convergence has been reached by then.
//' @export
// [[Rcpp::export]]
Rcpp::List paf_iter(arma::vec h2, double criterion, arma::mat R,
                    const int n_fac, bool abs_eig, int crit_type,
                    int max_iter) {

  int iter = 1;
  double delta = 1.0;
  arma::vec Lambda;
  arma::mat V;
  arma::mat L;
  arma::vec new_h2;
  arma::vec eigval;
  arma::mat eigvec;
  arma::mat Lt;

  while (delta > criterion) {
    //  compute the eigenvalues and eigenvectors using the R eigen function
    //eigen_list = eigs_sym_c(R, n_fac);

    // for clarity, store the m eigenvalues and
    // n_fac eigenvectors separately
    arma::mat xx(R.begin(),R.n_rows,R.n_cols, false);
    eigs_sym(eigval, eigvec, conv_to<sp_mat>::from(xx), n_fac);
    Lambda = flipud(eigval);
    V = fliplr(eigvec);

    if (abs_eig == false) {
      if (any(Lambda < 0)) {
        stop("Negative Eigenvalues detected; cannot compute communality estimates. Try again with init_comm = 'unity' or 'mac'");
      }


      // compute the loadings from the eigenvector matrix and diagonal
      // eigenvalue matrix
      L = V * arma::diagmat(sqrt(Lambda));
      // if (n_fac > 1) {
      //   L = V * arma::diagmat(sqrt(Lambda));
      // } else {
      //   L = V % sqrt(Lambda);
      // }

      // get the new communality estimates from the loadings
      Lt = L * L.t();
      new_h2 = Lt.diag();

    } else if (abs_eig == true) {

      // this is the implementation according to SPSS, which uses
      // absolute eigenvalues to avoid problems when having to compute
      // the square root

      // compute the loadings from the eigenvector matrix and diagonal
      // eigenvalue matrix

      L = V * arma::diagmat(sqrt(arma::abs(Lambda)));
      // if (n_fac > 1) {
      //   L = V * arma::diagmat(sqrt(arma::abs(Lambda)));
      // } else {
      //   L = V % sqrt(arma::abs(Lambda));
      // }

      // get the new communality estimates from the loadings
      // in SPSS implemented as rowSums(V^2 %*% diag(abs(Lambda)))
      // which is equivalent to

      //new_h2 = MM(L).diag();
      Lt = L * L.t();
      new_h2 = Lt.diag();
    }

    if (crit_type == 1) { // critertion type = "max_individual"

      // save the maximum change in the communality estimates
      delta = arma::abs(h2 - new_h2).max();

    } else if (crit_type == 2){ // critertion type = "sums"

      // convergence criterion according to the psych package
      delta = std::abs(arma::accu(h2) - arma::accu(new_h2));

    }

    // update diagonal of R with the new communality estimates
    R.diag() = new_h2;

    // update old communality estimates with new ones
    h2 = new_h2;

    // break if after maximum iterations there was no convergence
    if (iter > max_iter){
      warning("Reached maximum number of iterations without convergence. Results may not be interpretable.");
      break;
    }

    // incerase iterator
    iter += 1;

  }


  return Rcpp::List::create(Rcpp::Named("h2") = h2,
                            Rcpp::Named("R") = R,
                            Rcpp::Named("iter") = iter,
                            Rcpp::Named("L") = L);
}

