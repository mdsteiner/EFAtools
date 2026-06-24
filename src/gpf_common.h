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

// Scale each column of X to unit Euclidean norm, the projection onto the oblique
// constraint set diag(t(T) %*% T) = 1. A zero/non-finite column norm is replaced by a
// unit basis vector so the candidate stays on the constraint set; a singular candidate
// is subsequently rejected by the inverse check in the oblique manifold's compute().
inline arma::mat normalize_cols_cpp(const arma::mat& X, double eps = 1e-15) {
  arma::mat Y = X;

  for (arma::uword j = 0; j < Y.n_cols; ++j) {
    double nrm = arma::norm(Y.col(j), 2);

    if (!std::isfinite(nrm) || nrm < eps) {
      Y.col(j).zeros();
      Y(std::min(j, Y.n_rows - 1), j) = 1.0;
    } else {
      Y.col(j) /= nrm;
    }
  }

  return Y;
}

// Solve X * invX = I for a square, finite X, returning false (leaving invX untouched)
// if X is non-square, non-finite, or the solve fails or yields a non-finite result, so
// callers can reject a singular transformation rather than evaluate it through a
// pseudo-inverse.
inline bool inverse_checked_cpp(const arma::mat& X, arma::mat& invX) {
  if (X.n_rows != X.n_cols || !all_finite_cpp(X)) {
    return false;
  }

  try {
    bool ok = arma::solve(
      invX,
      X,
      arma::eye<arma::mat>(X.n_rows, X.n_cols),
      arma::solve_opts::fast
    );
    return ok && all_finite_cpp(invX);
  } catch (...) {
    return false;
  }
}

// Project the raw gradient G onto the tangent space of the oblique column-normalized
// manifold diag(t(T) %*% T) = 1: Gp = G - T diag(colSums(T .* G)); also report its
// Frobenius norm. Shared by the oblique solvers (the Procrustes objective and the
// criterion objective); only the upstream gradient G differs between them.
inline void oblique_tangent_project(const arma::mat& Tmat, const arma::mat& G,
                                    arma::mat& Gp, double& s) {
  arma::rowvec proj = arma::sum(Tmat % G, 0);
  Gp = G - Tmat * arma::diagmat(proj.t());
  s = arma::norm(arma::vectorise(Gp), 2);
}

// A uniformly distributed random k x k orthogonal matrix, constructed as in
// GPArotation::Random.Start(): a standard-normal matrix run through a QR
// decomposition, with the orthogonal factor's columns sign-fixed by the diagonal
// of R. The normal entries are filled row-major here rather than column-major as
// R's matrix() does, so for a given seed this yields a valid uniform orthogonal
// start but not the bit-identical matrix GPArotation would draw. Random draws use
// R::rnorm so they advance the calling process's RNG. Falls back to the identity
// if the decomposition fails.
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
