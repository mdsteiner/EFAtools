#ifndef EFATOOLS_GPF_COMMON_H
#define EFATOOLS_GPF_COMMON_H

#include <RcppArmadillo.h>
#include <cmath>

// Leaf helpers shared by the gradient-projection rotation solvers (the oblique
// Procrustes path in oblique_procrustes.cpp and the orthogonal GPForth path in
// rotate.cpp). They carry no coupling to a particular rotation objective, so they
// live here and are included by both translation units.

// Whether every entry of a matrix is finite (no NaN/Inf).
inline bool all_finite_cpp(const arma::mat& X) {
  for (arma::uword i = 0; i < X.n_elem; ++i) {
    if (!std::isfinite(X[i])) {
      return false;
    }
  }
  return true;
}

// Whether a scalar control value is finite.
inline bool is_valid_scalar_cpp(const double x) {
  return std::isfinite(x);
}

// A uniformly distributed random k x k orthogonal matrix, drawn the same way as
// GPArotation::Random.Start(): a standard-normal matrix run through a QR
// decomposition, with the orthogonal factor's columns sign-fixed by the diagonal
// of R. Random draws use R::rnorm so they advance the calling process's RNG. Falls
// back to the identity if the decomposition fails.
inline arma::mat random_orthogonal_start_cpp(const unsigned int k) {
  arma::mat Z(k, k);
  for (unsigned int i = 0; i < k; ++i) {
    for (unsigned int j = 0; j < k; ++j) {
      Z(i, j) = R::rnorm(0.0, 1.0);
    }
  }

  arma::mat Q, Rm;
  bool ok = arma::qr_econ(Q, Rm, Z);
  if (!ok || Q.n_rows != k || Q.n_cols != k) {
    arma::mat Eye(k, k, arma::fill::eye);
    return Eye;
  }

  arma::vec d = Rm.diag();
  for (unsigned int j = 0; j < d.n_elem; ++j) {
    double sgn = (d(j) >= 0.0) ? 1.0 : -1.0;
    Q.col(j) *= sgn;
  }

  return Q;
}

#endif  // EFATOOLS_GPF_COMMON_H
