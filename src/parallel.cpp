// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Parallel analysis on simulated data.
//'
//' Function called from within PARALLEL so usually no call to this is needed by the user.
//' Provides a C++ implementation of the PARALLEL simulation procedure
//'
//' @param n_datasets numeric. Number of datasets with dimensions (N, n_vars) to simulate.
//' @param n_vars numeric. Number of variables / indicators in dataset.
//' @param N numeric. Number of cases / observations in dataset.
//' @param eigen_type numeric. Whether PCA (eigen_type = 1; i.e., leaving diagonal of correlation matrix at 1) or PAF (eigen_type = 2; i.e., setting diagonal of correlation matrix to SMCs).
//' @param maxit numeric. Maximum iterations to perform after which to abort.
// [[Rcpp::export(.parallel_sim)]]
arma::mat parallel_sim(const int n_datasets, const int n_vars, const int N,
                         const int eigen_type, const int maxit = 10000) {
  // initialize needed objects
  arma::vec Lambda(n_vars);
  arma::vec eigval(n_vars);
  arma::mat eig_vals(n_datasets, n_vars);
  arma::mat x(N, n_vars);
  arma::mat R(n_vars, n_vars);

  if (eigen_type == 1) { // PCA

    // perform simulations for n_datasets time
    for (uword i = 0; i < n_datasets; i++) {
      x = randn(N, n_vars);
      R = cor(x);
      eig_sym(eigval, R);
      Lambda = flipud(eigval);
      eig_vals.row(i) = Lambda.t();
    }

  } else if (eigen_type == 2) { // SMC
    arma::vec smc(n_vars);
    arma::mat temp(n_vars, n_vars);
    int success = 0;
    int iter = 0;

    while((success < n_datasets) && (iter < maxit)) {
      x = randn(N, n_vars);
      R = cor(x);
      bool flag = inv_sympd(temp, R);
      if (!flag){
        iter++;
        continue;
      }
      R.diag() = 1 - (1 / temp.diag());
      eig_sym(eigval, R);
      Lambda = flipud(eigval);
      eig_vals.row(success) = Lambda.t();
      iter++;
      success++;
    }

    if ((iter == maxit) && (success < n_datasets)) {
      stop("Could not generate enough non-singular matrices.");
    }

    // for (uword i = 0; i < n_datasets; i++) {
    //   x = randn(N, n_vars);
    //   R = cor(x);
    //   temp = inv(R); // previously inv_sympd()
    //   R.diag() = 1 - (1 / temp.diag());
    //   eig_sym(eigval, R);
    //   Lambda = flipud(eigval);
    //   eig_vals.row(i) = Lambda.t();
    // }

  }

  return eig_vals;

}


// =============================================================================


// //' Principal Axis Factoring to extract eigenvalues from a 1 factor solution
// //'
// //' Function called from within PARALLEL so usually no call to this is needed by the user.
// //' Provides a C++ implementation of 1 factor PAF. Returns the eigenvalues obtained
// //' from the correlation matrix with the final communality estimates as diagonal of
// //' R.
// //'
// //' @param R numeric matrix. Correlation matrix to perform PAF with 1 factor solution on.
// //' @param criterion double. Convergence criterion to use.
// //' @param crit_type integer. Whether max_individual (1) or sums (2).
// //' @param max_iter integer. The maximum number of iterations after which to stop the iterative procedure if no convergence is reached by then.
// //' @export
// // [[Rcpp::export]]
// arma::vec parallel_paf(arma::mat R, double criterion, int crit_type,
//                        int max_iter) {
//
//
//   int iter = 1;
//   double delta = 1.0;
//   arma::vec tv(R.n_cols);
//   arma::vec Lambda(R.n_cols);
//   arma::vec Lambda_o(1);
//   arma::mat V(R.n_cols, R.n_cols);
//   arma::mat V_o(R.n_cols, 1);
//   arma::mat L;
//   arma::vec new_h2;
//   arma::vec eigval;
//   arma::mat eigvec;
//   arma::mat Lt;
//   arma::vec h2;
//
//   // compute smcs
//   arma::mat temp(R.n_cols, R.n_cols);
//   temp = inv_sympd(R);
//   R.diag() = 1 - (1 / temp.diag());
//   h2 = temp.diag();
//
//   if (crit_type == 1) { // "max_individual"
//
//     while (delta > criterion & iter <= max_iter) {
//       //  compute the eigenvalues and eigenvectors
//       eig_sym(eigval, eigvec, R);
//       Lambda = flipud(eigval);
//       Lambda_o = Lambda(0);
//       V = fliplr(eigvec);
//       V_o = V.col(0);
//
//       if (any(Lambda_o < 0)) {
//         stop("Negative Eigenvalues detected; cannot compute communality estimates. Try again with init_comm = 'unity' or 'mac'");
//       }
//
//       // compute the loadings from the eigenvector matrix and diagonal
//       // eigenvalue matrix
//       Lambda_o = arma::sqrt(Lambda_o);
//       tv.fill(Lambda_o[0]);
//       L = V_o % tv;
//
//       // get the new communality estimates from the loadings
//       Lt = L * L.t();
//       new_h2 = Lt.diag();
//
//       // save the maximum change in the communality estimates
//       delta = arma::abs(h2 - new_h2).max();
//
//       // update diagonal of R with the new communality estimates
//       R.diag() = new_h2;
//
//       // update old communality estimates with new ones
//       h2 = new_h2;
//
//       // incerase iterator
//       iter += 1;
//
//     }
//
//   } else if (crit_type == 2) { // "sums"
//
//     while (delta > criterion & iter <= max_iter) {
//       //  compute the eigenvalues and eigenvectors
//       eig_sym(eigval, eigvec, R);
//       Lambda = flipud(eigval);
//       Lambda_o = Lambda(0);
//       V = fliplr(eigvec);
//       V_o = V.col(0);
//
//       if (any(Lambda_o < 0)) {
//         stop("Negative Eigenvalues detected; cannot compute communality estimates. Try again with init_comm = 'unity' or 'mac'");
//       }
//
//       // compute the loadings from the eigenvector matrix and diagonal
//       // eigenvalue matrix
//       Lambda_o = arma::sqrt(Lambda_o);
//       tv.fill(Lambda_o[0]);
//       L = V_o % tv;
//
//       // get the new communality estimates from the loadings
//       Lt = L * L.t();
//       new_h2 = Lt.diag();
//
//       // convergence criterion according to the psych package
//       delta = std::abs(arma::accu(h2) - arma::accu(new_h2));
//
//       // update diagonal of R with the new communality estimates
//       R.diag() = new_h2;
//
//       // update old communality estimates with new ones
//       h2 = new_h2;
//
//       // incerase iterator
//       iter += 1;
//
//     }
//
//   }
//
//   // break if after maximum iterations there was no convergence
//   if (iter >= max_iter){
//     warning("Reached maximum number of iterations without convergence. Results may not be interpretable.");
//   }
//
//   eig_sym(eigval, R);
//   Lambda = flipud(eigval);
//
//   return Lambda;
//
// }

// //' Parallel analysis on simulated data.
// //'
// //' Function called from within PARALLEL so usually no call to this is needed by the user.
// //' Provides a C++ implementation of the PARALLEL simulation procedure where eigenvalues
// //' are found using the parallel_paf function.
// //'
// //' @param n_datasets numeric. Number of datasets with dimensions (N, n_vars) to simulate.
// //' @param n_vars numeric. Number of variables / indicators in dataset.
// //' @param N numeric. Number of cases / observations in dataset.
// //' @param criterion double. Convergence criterion to use.
// //' @param crit_type integer. Whether max_individual (1) or sums (2).
// //' @param max_iter integer. The maximum number of iterations after which to stop the iterative procedure if no convergence is reached by then.
// //' @export
// // [[Rcpp::export]]
// arma::mat parallel_paf_sim(const int n_datasets, const int n_vars, const int N,
//                        double criterion, int crit_type, int max_iter) {
//   // initialize needed objects
//   arma::vec Lambda(n_vars);
//   arma::vec eigval(n_vars);
//   arma::mat eig_vals(n_datasets, n_vars);
//   arma::mat x(N, n_vars);
//   arma::mat R(n_vars, n_vars);
//   arma::vec smc(n_vars);
//   arma::mat temp(n_vars, n_vars);
//
//   for (uword i = 0; i < n_datasets; i++) {
//     x = randn(N, n_vars);
//     R = cor(x);
//     Lambda = parallel_paf(R,criterion, crit_type, max_iter);
//     eig_vals.row(i) = Lambda.t();
//   }
//
//   return eig_vals;
//
// }

// //' Parallel analysis on resampled real data.
// //'
// //' Function called from within PARALLEL so usually no call to this is needed by the user.
// //' Provides a C++ implementation of the PARALLEL resampling procedure where eigenvalues
// //' are found using the parallel_paf function.
// //'
// //' @param n_datasets numeric. Number of datasets to simulate.
// //' @param data numeric matrix. The real data matrix to perform resampling on.
// //' @param replace logical. Should resampling be done with replacement (TRUE) or without (FALSE).
// //' @param criterion double. Convergence criterion to use.
// //' @param crit_type integer. Whether max_individual (1) or sums (2).
// //' @param max_iter integer. The maximum number of iterations after which to stop the iterative procedure if no convergence is reached by then.
// //' @export
// // [[Rcpp::export]]
// arma::mat parallel_paf_resample(const int n_datasets, arma::mat data,
//                             const bool replace, double criterion, int crit_type,
//                             int max_iter) {
//   const int n_vars = data.n_cols;
//   const int N = data.n_rows;
//   arma::vec Lambda(n_vars);
//   arma::vec eigval(n_vars);
//   arma::mat eig_vals(n_datasets, n_vars);
//   arma::mat x(N, n_vars);
//   arma::mat R(n_vars, n_vars);
//
//   if (replace == false) { // shuffle within columns
//
//     arma::vec smc(n_vars);
//     arma::mat temp(n_vars, n_vars);
//
//     for (uword i = 0; i < n_datasets; i++) {
//       for (uword cc = 0; cc < n_vars; i++) {
//         x.col(i) = shuffle(data.col(i));
//       }
//       R = cor(x);
//       temp = inv_sympd(R);
//       R.diag() = 1 - (1 / temp.diag());
//       Lambda = parallel_paf(R,criterion, crit_type, max_iter);
//       eig_vals.row(i) = Lambda.t();
//     }
//
//
//   } else if (replace == true) { // sample with replacement
//
//     arma::vec smc(n_vars);
//     arma::mat temp(n_vars, n_vars);
//
//     for (uword i = 0; i < n_datasets; i++) {
//       for (uword rr = 0; rr < N; i++) {
//         for (uword cc = 0; cc < n_vars; i++) {
//           x(rr, cc) = data(randi(distr_param(0, N - 1)), cc);
//         }
//       }
//       R = cor(x);
//       temp = inv_sympd(R);
//       R.diag() = 1 - (1 / temp.diag());
//       Lambda = parallel_paf(R,criterion, crit_type, max_iter);
//       eig_vals.row(i) = Lambda.t();
//     }
//
//   }
//
//   return eig_vals;
//
// }
//
//
//
// // =============================================================================
//
//
//
// //' Parallel analysis on resampled real data.
// //'
// //' Function called from within PARALLEL so usually no call to this is needed by the user.
// //' Provides a C++ implementation of the PARALLEL resampling procedure
// //'
// //' @param n_datasets numeric. Number of datasets to simulate.
// //' @param data numeric matrix. The real data matrix to perform resampling on.
// //' @param eigen_type numeric. Whether PCA (eigen_type = 1; i.e., leaving diagonal of correlation matrix at 1) or PAF (eigen_type = 2; i.e., setting diagonal of correlation matrix to SMCs).
// //' @param replace logical. Should resampling be done with replacement (TRUE) or without (FALSE).
// //' @export
// // [[Rcpp::export]]
// arma::mat parallel_resample(const int n_datasets, arma::mat data, const int eigen_type,
//                             const bool replace) {
//   const int n_vars = data.n_cols;
//   const int N = data.n_rows;
//   arma::vec Lambda(n_vars);
//   arma::vec eigval(n_vars);
//   arma::mat eig_vals(n_datasets, n_vars);
//   arma::mat x(N, n_vars);
//   arma::mat R(n_vars, n_vars);
//
//   if (replace == false) { // shuffle within columns
//
//     if (eigen_type == 1) { // PCA
//
//       for (uword i = 0; i < n_datasets; i++) {
//
//         for (uword cc = 0; cc < n_vars; i++) {
//           x.col(i) = shuffle(data.col(i));
//         }
//         R = cor(x);
//         eig_sym(eigval, R);
//         Lambda = flipud(eigval);
//         eig_vals.row(i) = Lambda.t();
//       }
//
//     } else if (eigen_type == 2) { // PAF
//
//       arma::vec smc(n_vars);
//       arma::mat temp(n_vars, n_vars);
//
//       for (uword i = 0; i < n_datasets; i++) {
//         for (uword cc = 0; cc < n_vars; i++) {
//           x.col(i) = shuffle(data.col(i));
//         }
//         R = cor(x);
//         temp = inv_sympd(R);
//         R.diag() = 1 - (1 / temp.diag());
//         eig_sym(eigval, R);
//         Lambda = flipud(eigval);
//         eig_vals.row(i) = Lambda.t();
//       }
//
//     }
//
//   } else if (replace == true) { // sample with replacement
//
//     if (eigen_type == 1) { // PCA
//
//       for (uword i = 0; i < n_datasets; i++) {
//
//         for (uword rr = 0; rr < N; i++) {
//           for (uword cc = 0; cc < n_vars; i++) {
//             x(rr, cc) = data(randi(distr_param(0, N - 1)), cc);
//           }
//         }
//         R = cor(x);
//         eig_sym(eigval, R);
//         Lambda = flipud(eigval);
//         eig_vals.row(i) = Lambda.t();
//       }
//
//     } else if (eigen_type == 2) { // PAF
//
//       arma::vec smc(n_vars);
//       arma::mat temp(n_vars, n_vars);
//
//       for (uword i = 0; i < n_datasets; i++) {
//         for (uword rr = 0; rr < N; i++) {
//           for (uword cc = 0; cc < n_vars; i++) {
//             x(rr, cc) = data(randi(distr_param(0, N - 1)), cc);
//           }
//         }
//         R = cor(x);
//         temp = inv_sympd(R);
//         R.diag() = 1 - (1 / temp.diag());
//         eig_sym(eigval, R);
//         Lambda = flipud(eigval);
//         eig_vals.row(i) = Lambda.t();
//       }
//     }
//
//   }
//
//   return eig_vals;
//
// }

// //' Summarise the raw data from the .parallel_sim
// //'
// //' Function called from within PARALLEL so usually no call to this is needed by the user.
// //' Provides a C++ implementation to aggregate the eigenvalues from the simulations
// //' performed using .parallel_sim.
// //'
// //' @param eig_vals matrix. A matrix as returned by .parallel_sim.
// //' @param percent numeric. A vector of percentiles for which the eigenvalues should be returned.
// //' @param n_datasets integer. The number of datasets simulated in .parallel_sim.
// //' @param n_vars numeric. The number of variables / indicators per dataset.
// //' @export
// // [[Rcpp::export(.parallel_summarise)]]
// NumericMatrix parallel_summarise(NumericMatrix eig_vals, NumericVector percent,
//                                  const int n_datasets, const int n_vars) {
//   NumericMatrix results(n_vars, 1 + percent.length());
//   int ind;
//   NumericVector temp(eig_vals.nrow());
//   results(_, 0) = colMeans(eig_vals);
//
//   for (int root = 0; root < n_vars; root++) {
//     for (int perc_i = 0; perc_i < percent.length(); perc_i++) {
//       ind = round((percent(perc_i) * n_datasets) / 100) - 1;
//       temp = eig_vals.column(root);
//       results(root, 1 + perc_i) = temp.sort()[ind];
//     }
//   }
//
//   return results;
//
// }
