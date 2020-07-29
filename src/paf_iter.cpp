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
// [[Rcpp::export(.paf_iter)]]
Rcpp::List paf_iter(arma::vec h2, double criterion, arma::mat R,
                    const int n_fac, bool abs_eig, int crit_type,
                    int max_iter) {


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
  arma::uvec idx(n_fac);
  idx.fill(true);

  if (abs_eig == false) {

    if (n_fac > 1) {

      if (crit_type == 1) {

        while ((delta > criterion) & (iter <= max_iter)) {
          //  compute the eigenvalues and eigenvectors
          eig_sym(eigval, eigvec, R);
          Lambda = flipud(eigval);
          Lambda = Lambda.elem(arma::find(idx));
          V = fliplr(eigvec);
          V = V.cols(arma::find(idx));

          if (any(Lambda < 0)) {
            stop("Negative Eigenvalues detected; cannot compute communality estimates. Try again with init_comm = 'unity' or 'mac'");
          }

          // compute the loadings from the eigenvector matrix and diagonal
          // eigenvalue matrix
          L = V * arma::diagmat(sqrt(Lambda));

          // get the new communality estimates from the loadings
          Lt = L * L.t();
          new_h2 = Lt.diag();

          // save the maximum change in the communality estimates
          delta = arma::abs(h2 - new_h2).max();

          // update diagonal of R with the new communality estimates
          R.diag() = new_h2;

          // update old communality estimates with new ones
          h2 = new_h2;

          // incerase iterator
          iter += 1;

        }

      } else if (crit_type == 2) {

        while ((delta > criterion) & (iter <= max_iter)) {
          //  compute the eigenvalues and eigenvectors
          eig_sym(eigval, eigvec, R);
          Lambda = flipud(eigval);
          Lambda = Lambda.elem(arma::find(idx));
          V = fliplr(eigvec);
          V = V.cols(arma::find(idx));

          if (any(Lambda < 0)) {
            stop("Negative Eigenvalues detected; cannot compute communality estimates. Try again with init_comm = 'unity' or 'mac'");
          }

          // compute the loadings from the eigenvector matrix and diagonal
          // eigenvalue matrix
          L = V * arma::diagmat(sqrt(Lambda));

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


    } else {


      if (crit_type == 1) {

        while ((delta > criterion) & (iter <= max_iter)) {
          //  compute the eigenvalues and eigenvectors
          eig_sym(eigval, eigvec, R);
          Lambda = flipud(eigval);
          Lambda = Lambda.elem(arma::find(idx));
          V = fliplr(eigvec);
          V = V.cols(arma::find(idx));

          if (any(Lambda < 0)) {
            stop("Negative Eigenvalues detected; cannot compute communality estimates. Try again with init_comm = 'unity' or 'mac'");
          }

          // compute the loadings from the eigenvector matrix and diagonal
          // eigenvalue matrix
          Lambda = arma::sqrt(Lambda);
          tv.fill(Lambda[0]);
          L = V % tv;

          // get the new communality estimates from the loadings
          Lt = L * L.t();
          new_h2 = Lt.diag();

          // save the maximum change in the communality estimates
          delta = arma::abs(h2 - new_h2).max();

          // update diagonal of R with the new communality estimates
          R.diag() = new_h2;

          // update old communality estimates with new ones
          h2 = new_h2;

          // incerase iterator
          iter += 1;

        }

      } else if (crit_type == 2) {

        while ((delta > criterion) & (iter <= max_iter)) {
          //  compute the eigenvalues and eigenvectors
          eig_sym(eigval, eigvec, R);
          Lambda = flipud(eigval);
          Lambda = Lambda.elem(arma::find(idx));
          V = fliplr(eigvec);
          V = V.cols(arma::find(idx));

          if (any(Lambda < 0)) {
            stop("Negative Eigenvalues detected; cannot compute communality estimates. Try again with init_comm = 'unity' or 'mac'");
          }

          // compute the loadings from the eigenvector matrix and diagonal
          // eigenvalue matrix
          Lambda = arma::sqrt(Lambda);
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

    }

  } else if (abs_eig == true) {

    if (n_fac > 1) {

      if (crit_type == 1) {

        while ((delta > criterion) & (iter <= max_iter)) {
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

          // save the maximum change in the communality estimates
          delta = arma::abs(h2 - new_h2).max();

          // update diagonal of R with the new communality estimates
          R.diag() = new_h2;

          // update old communality estimates with new ones
          h2 = new_h2;

          // incerase iterator
          iter += 1;

        }

      } else if (crit_type == 2) {

        while ((delta > criterion) & (iter <= max_iter)) {
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

      }


    } else {


      if (crit_type == 1) {

        while ((delta > criterion) & (iter <= max_iter)) {
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

          // save the maximum change in the communality estimates
          delta = arma::abs(h2 - new_h2).max();

          // update diagonal of R with the new communality estimates
          R.diag() = new_h2;

          // update old communality estimates with new ones
          h2 = new_h2;

          // incerase iterator
          iter += 1;

        }

      } else if (crit_type == 2) {

        while ((delta > criterion) & (iter <= max_iter)) {
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

    }

  }

  // break if after maximum iterations there was no convergence
  if (iter >= max_iter){
    warning("Reached maximum number of iterations without convergence. Results may not be interpretable.");
  }

  return Rcpp::List::create(Rcpp::Named("h2") = h2,
                            Rcpp::Named("R") = R,
                            Rcpp::Named("iter") = iter,
                            Rcpp::Named("L") = L);
}

