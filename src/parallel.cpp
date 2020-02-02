// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Parallel analysis on simulated data.
//'
//' Function called from within PARALLEL so usually no call to this is needed by the user.
//' Provides a C++ implementation of the PARALLEL simulation procedure
//'
//' @param ndatasets numeric. Number of datasets with dimensions (ncases, nvars) to simulate.
//' @param nvars numeric. Number of variables / indicators in dataset.
//' @param ncases numeric. Number of cases / observations in dataset.
//' @param kind numeric. Whether PCA (kind = 1; i.e., leaving diagonal of correlation matrix at 1) or PAF (kind = 2; i.e., setting diagonal of correlation matrix to SMCs).
//' @export
// [[Rcpp::export]]
arma::mat parallel_sim(const int ndatasets, const int nvars, const int ncases,
                         const int kind) {
  // initialize needed objects
  arma::vec Lambda(nvars);
  arma::vec eigval(nvars);
  arma::mat eig_vals(ndatasets, nvars);
  arma::mat x(ncases, nvars);
  arma::mat R(nvars, nvars);

  if (kind == 1) { // PCA

    // perform simulations for ndatasets time
    for (uword i = 0; i < ndatasets; i++) {
      x = randn(ncases, nvars);
      R = cor(x);
      eig_sym(eigval, R);
      Lambda = flipud(eigval);
      eig_vals.row(i) = Lambda.t();
    }

  } else if (kind == 2) { // PAF
    arma::vec smc(nvars);
    arma::mat temp(nvars, nvars);

    for (uword i = 0; i < ndatasets; i++) {
      x = randn(ncases, nvars);
      R = cor(x);
      temp = inv_sympd(R);
      R.diag() = 1 - (1 / temp.diag());
      eig_sym(eigval, R);
      Lambda = flipud(eigval);
      eig_vals.row(i) = Lambda.t();
    }

  }

  return eig_vals;

}


//' Parallel analysis on resampled real data.
//'
//' Function called from within PARALLEL so usually no call to this is needed by the user.
//' Provides a C++ implementation of the PARALLEL resampling procedure
//'
//' @param ndatasets numeric. Number of datasets to simulate.
//' @param data numeric matrix. The real data matrix to perform resampling on.
//' @param kind numeric. Whether PCA (kind = 1; i.e., leaving diagonal of correlation matrix at 1) or PAF (kind = 2; i.e., setting diagonal of correlation matrix to SMCs).
//' @param replace logical. Should resampling be done with replacement (TRUE) or without (FALSE).
//' @export
// [[Rcpp::export]]
arma::mat parallel_resample(const int ndatasets, arma::mat data, const int kind,
                            const bool replace) {
  const int nvars = data.n_cols;
  const int ncases = data.n_rows;
  arma::vec Lambda(nvars);
  arma::vec eigval(nvars);
  arma::mat eig_vals(ndatasets, nvars);
  arma::mat x(ncases, nvars);
  arma::mat R(nvars, nvars);

  if (replace == false) { // shuffle within columns

    if (kind == 1) { // PCA

      for (uword i = 0; i < ndatasets; i++) {

        for (uword cc = 0; cc < nvars; i++) {
          x.col(i) = shuffle(data.col(i));
        }
        R = cor(x);
        eig_sym(eigval, R);
        Lambda = flipud(eigval);
        eig_vals.row(i) = Lambda.t();
      }

    } else if (kind == 2) { // PAF

      arma::vec smc(nvars);
      arma::mat temp(nvars, nvars);

      for (uword i = 0; i < ndatasets; i++) {
        for (uword cc = 0; cc < nvars; i++) {
          x.col(i) = shuffle(data.col(i));
        }
        R = cor(x);
        temp = inv_sympd(R);
        R.diag() = 1 - (1 / temp.diag());
        eig_sym(eigval, R);
        Lambda = flipud(eigval);
        eig_vals.row(i) = Lambda.t();
      }

    }

  } else if (replace == true) { // sample with replacement

    if (kind == 1) { // PCA

      for (uword i = 0; i < ndatasets; i++) {

        for (uword rr = 0; rr < ncases; i++) {
          for (uword cc = 0; cc < nvars; i++) {
            x(rr, cc) = data(randi(distr_param(0, ncases - 1)), cc);
          }
        }
        R = cor(x);
        eig_sym(eigval, R);
        Lambda = flipud(eigval);
        eig_vals.row(i) = Lambda.t();
      }

    } else if (kind == 2) { // PAF

      arma::vec smc(nvars);
      arma::mat temp(nvars, nvars);

      for (uword i = 0; i < ndatasets; i++) {
        for (uword rr = 0; rr < ncases; i++) {
          for (uword cc = 0; cc < nvars; i++) {
            x(rr, cc) = data(randi(distr_param(0, ncases - 1)), cc);
          }
        }
        R = cor(x);
        temp = inv_sympd(R);
        R.diag() = 1 - (1 / temp.diag());
        eig_sym(eigval, R);
        Lambda = flipud(eigval);
        eig_vals.row(i) = Lambda.t();
      }
    }

  }

  return eig_vals;

}

//' Summarise the raw data from the \link{parallel_sim} and \link{parallel_resample}
//'
//' Function called from within PARALLEL so usually no call to this is needed by the user.
//' Provides a C++ implementation to aggregate the eigenvalues from the simulations
//' performed using \link{parallel_sim} and \link{parallel_resample}.
//'
//' @param eig_vals matrix. A matrix as returned by \link{parallel_sim} or \link{parallel_resample}.
//' @param percent numeric. A vector of percentiles for which the eigenvalues should be returned.
//' @param ndatasets integer. The number of datasets simulated in \link{parallel_sim} or \link{parallel_resample}.
//' @param nvars numeric. The number of variables / indicators per dataset.
//' @export
// [[Rcpp::export]]
NumericMatrix parallel_summarise(NumericMatrix eig_vals, NumericVector percent,
                                 const int ndatasets, const int nvars) {
  NumericMatrix results(nvars, 1 + percent.length());
  int ind;
  NumericVector temp;
  results(_, 0) = colMeans(eig_vals);

  for (int root = 0; root < nvars; root++) {
    for (int perc_i = 0; perc_i < percent.length(); perc_i++) {
      ind = round((percent(perc_i) * ndatasets) / 100) - 1;
      temp = eig_vals.column(root);
      results(root, 1 + perc_i) = temp.sort().operator[](ind);
    }
  }

  return results;

}
