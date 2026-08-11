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

// The leading n_fac eigenvectors, in descending eigenvalue order. arma::eig_sym returns
// ascending eigenvalues, so the trailing n_fac columns are the leading eigenpairs and
// reversing just those gives the same columns in the same order as reversing the whole
// matrix would -- without allocating and copying the full p x p reversal on every call,
// which matters because the callers sit in optimiser inner loops. The bound is checked
// here, beside the indexing it protects: an out-of-range n_fac would read past the
// available columns (undefined behaviour in an unchecked build), and n_fac is signed
// while the subview index is not. Callers validate their own n_fac up front and report
// it in their own terms; this is the backstop.
inline arma::mat top_eigvec(const arma::mat& eigvec, int n_fac) {
  if (n_fac < 1 || static_cast<arma::uword>(n_fac) > eigvec.n_cols) {
    Rcpp::stop("n_fac must be at least 1 and at most the number of variables "
               "(got n_fac = %d for %d variables).",
               n_fac, static_cast<int>(eigvec.n_cols));
  }
  return arma::fliplr(eigvec.tail_cols(n_fac));
}

// Loadings from the leading n_fac eigenpairs with negative eigenvalues clipped to zero:
// the least-squares extraction shared by the ULS objective and the DWLS warm start, and
// the extraction stats::eigen() would give at the optimum.
inline arma::mat top_eigen_loadings(const arma::vec& eigval, const arma::mat& eigvec,
                                    int n_fac) {
  arma::mat V = top_eigvec(eigvec, n_fac);   // checks the bound for the tail() below
  arma::vec lambda = arma::reverse(eigval.tail(n_fac));
  lambda.elem(arma::find(lambda < 0)).zeros();
  return V * arma::diagmat(arma::sqrt(lambda));
}

#endif
