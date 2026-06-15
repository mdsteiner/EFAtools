// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Symmetric eigendecomposition that fails loudly instead of silently leaving
// eigval/eigvec empty (which the callers would then index out of bounds). A
// non-finite or non-symmetric matrix - e.g. a degenerate bootstrap correlation
// matrix - makes arma::eig_sym return false, so turn that into a catchable R
// error rather than undefined behaviour.
static void eig_sym_checked(arma::vec& eigval, arma::mat& eigvec,
                            const arma::mat& X) {
  if (!arma::eig_sym(eigval, eigvec, X)) {
    Rcpp::stop("Eigendecomposition failed during factor extraction; the "
               "correlation matrix is not finite or not symmetric.");
  }
}

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

  // Guard the eigen-extraction below: it reads the largest n_fac eigenpairs and
  // would index past the available eigenvalues if n_fac >= ncol(R) (undefined
  // behaviour / crash in an unchecked build). Convert that into a catchable error.
  if (n_fac < 1 || static_cast<arma::uword>(n_fac) >= R.n_cols) {
    Rcpp::stop("n_fac must be at least 1 and smaller than the number of "
               "variables (got n_fac = %d for %d variables).",
               n_fac, static_cast<int>(R.n_cols));
  }

  // Only crit_type 1 and 2 have an iteration loop below; any other value would
  // skip every branch and return an uninitialised (empty) loading matrix.
  if (crit_type != 1 && crit_type != 2) {
    Rcpp::stop("crit_type must be 1 (maximum individual change) or 2 (sum of "
               "changes).");
  }

  int iter = 1;
  // Start above the criterion so the iteration runs at least once even when
  // criterion >= 1 (EFA() rejects that, but keep the kernel self-contained).
  double delta = criterion + 1.0;
  // Flagged (instead of aborting here) when an eigenvalue turns negative in a
  // non-absolute branch, so PAF() can raise a classed R condition.
  bool neg_eigen = false;
  arma::vec tv(R.n_cols);
  arma::vec Lambda;
  arma::mat V;
  arma::mat L;
  arma::vec new_h2;
  arma::vec eigval;
  arma::mat eigvec;
  arma::mat Lt;

  if (abs_eig == false) {

    if (n_fac > 1) {

      if (crit_type == 1) {

        while ((delta > criterion) & (iter <= max_iter)) {
          //  compute the eigenvalues and eigenvectors
          eig_sym_checked(eigval, eigvec, R);
          Lambda = flipud(eigval);
          Lambda = Lambda.head(n_fac);
          V = fliplr(eigvec);
          V = V.head_cols(n_fac);

          if (any(Lambda < 0)) {
            neg_eigen = true;
            break;
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
          eig_sym_checked(eigval, eigvec, R);
          Lambda = flipud(eigval);
          Lambda = Lambda.head(n_fac);
          V = fliplr(eigvec);
          V = V.head_cols(n_fac);

          if (any(Lambda < 0)) {
            neg_eigen = true;
            break;
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
          eig_sym_checked(eigval, eigvec, R);
          Lambda = flipud(eigval);
          Lambda = Lambda.head(n_fac);
          V = fliplr(eigvec);
          V = V.head_cols(n_fac);

          if (any(Lambda < 0)) {
            neg_eigen = true;
            break;
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
          eig_sym_checked(eigval, eigvec, R);
          Lambda = flipud(eigval);
          Lambda = Lambda.head(n_fac);
          V = fliplr(eigvec);
          V = V.head_cols(n_fac);

          if (any(Lambda < 0)) {
            neg_eigen = true;
            break;
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
          eig_sym_checked(eigval, eigvec, R);
          Lambda = arma::abs(flipud(eigval));
          Lambda = Lambda.head(n_fac);
          V = fliplr(eigvec);
          V = V.head_cols(n_fac);

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
          eig_sym_checked(eigval, eigvec, R);
          Lambda = arma::abs(flipud(eigval));
          Lambda = Lambda.head(n_fac);
          V = fliplr(eigvec);
          V = V.head_cols(n_fac);

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
          eig_sym_checked(eigval, eigvec, R);
          Lambda = arma::abs(flipud(eigval));
          Lambda = Lambda.head(n_fac);
          V = fliplr(eigvec);
          V = V.head_cols(n_fac);

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
          eig_sym_checked(eigval, eigvec, R);
          Lambda = arma::abs(flipud(eigval));
          Lambda = Lambda.head(n_fac);
          V = fliplr(eigvec);
          V = V.head_cols(n_fac);

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

  // Return the final change in the communality estimates and the
  // negative-eigenvalue flag so PAF() can raise the convergence and
  // negative-eigenvalue conditions as classed R conditions.
  return Rcpp::List::create(Rcpp::Named("h2") = h2,
                            Rcpp::Named("R") = R,
                            Rcpp::Named("iter") = iter,
                            Rcpp::Named("L") = L,
                            Rcpp::Named("delta") = delta,
                            Rcpp::Named("neg_eigen") = neg_eigen);
}

