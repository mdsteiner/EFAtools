#ifndef EFATOOLS_EIG_UTILS_H
#define EFATOOLS_EIG_UTILS_H

#include <RcppArmadillo.h>

// Symmetric eigendecompositions that fail loudly instead of silently leaving their
// outputs empty (which the callers would then index out of bounds). A non-finite or
// non-symmetric matrix -- a rescaled matrix that went pathological mid-optimization, a
// degenerate bootstrap correlation matrix, or a degenerate simulated one -- makes
// arma::eig_sym return false, so turn that into a catchable R error rather than
// undefined behaviour. Shared by the factor extractors (estimate.cpp, paf_iter.cpp) and
// the simulation kernels (sim.cpp, parallel.cpp).

// Eigenvalues and eigenvectors, for the eigen-extraction estimators.
inline void eig_sym_checked(arma::vec& eigval, arma::mat& eigvec,
                            const arma::mat& X) {
  if (!arma::eig_sym(eigval, eigvec, X)) {
    Rcpp::stop("Eigendecomposition failed during factor extraction; the "
               "correlation matrix is not finite or not symmetric.");
  }
}

// Eigenvalues only, for the simulation kernels. `during` names the simulation in the
// error message. Deliberately a separate name rather than an overload of the above:
// arma::mat converts implicitly from const char*, so a three-argument overload set
// would make every call with a string literal ambiguous.
inline void eig_sym_values_checked(arma::vec& eigval, const arma::mat& X,
                                   const char* during) {
  if (!arma::eig_sym(eigval, X)) {
    Rcpp::stop("Eigendecomposition failed during %s; the simulated correlation "
               "matrix is not finite or not symmetric.", during);
  }
}

#endif
