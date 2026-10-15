// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "eig_utils.h"

using namespace Rcpp;
using namespace arma;

// Exact test for the null / no-factor model: whether the matrix square root M is the
// identity (R = I in simulate_cfm_mvn, or the nf == 1 identity reference in
// simulate_cfm_eigen). It allocates a p x p identity and a p x p comparison, and M is
// fixed for a whole run, so each caller evaluates it once rather than per replicate.
static inline bool cfm_is_identity(const arma::mat& M) {
  return M.is_square() &&
         arma::all(arma::vectorise(M == arma::eye<arma::mat>(M.n_rows, M.n_cols)));
}

// Matrix square root M with M' M = R: the upper Cholesky factor for a positive-definite R
// (arma::chol returns M with M.t() * M == R). A positive-semidefinite but singular R
// (which makes the Cholesky fail although it is still a valid covariance, e.g. a
// no-factor block, a smoothed factor intercorrelation matrix, or a reference model with
// a zero uniqueness) falls back to the symmetric eigen square root V diag(sqrt(d)) V',
// which handles the zero-eigenvalue case and keeps the columns aligned with those of R.
static arma::mat cfm_sqrt(const arma::mat& R, const double tol) {
  arma::mat M;
  if (arma::chol(M, R)) {
    return M;
  }
  arma::vec eigval;
  arma::mat eigvec;
  if (!arma::eig_sym(eigval, eigvec, R)) {
    Rcpp::stop("Eigendecomposition of the population matrix failed; it is not "
               "finite or not symmetric.");
  }
  // A markedly negative eigenvalue means R is indefinite and cannot be a covariance to
  // simulate from. The R callers screen for this, so this branch is only a defensive
  // backstop.
  if (eigval.min() < -tol) {
    Rcpp::stop("The population matrix is not positive semi-definite; it cannot "
               "be a covariance to simulate from.");
  }
  eigval.transform([](double v) { return v < 0.0 ? 0.0 : std::sqrt(v); });
  return eigvec * arma::diagmat(eigval) * eigvec.t();
}

// Simulated correlation matrices without simulated data. For N cases drawn iid from
// N_p(0, Sigma), the centred cross-product matrix W = (N - 1) S has the Wishart
// distribution W_p(n, Sigma) with n = N - 1 degrees of freedom, and the Pearson
// correlation matrix of the cases is exactly cov2cor(W). By the Bartlett decomposition
// (Bartlett, 1933; Anderson, 2003, An Introduction to Multivariate Statistical Analysis,
// 3rd ed., Sec. 7.2), U' U has the W_p(n, I) distribution for an upper-triangular U with
// independent entries U_jj = sqrt(chi^2 with n - j + 1 df), j = 1..p, and U_ij ~ N(0, 1)
// for i < j; M' U' U M then has the W_p(n, Sigma) distribution for any M with
// M' M = Sigma. One draw therefore costs O(p^3) whatever N is, while correlating
// simulated data costs O(N p^2). n >= p is needed for all df to be positive; the callers
// ensure N > p.
//
// wishart_factor() refills the upper triangle of U in place; the caller zeroes the strict
// lower triangle once. The draw order (column by column, the diagonal before the entries
// above it) is that of stats::rWishart(), so a draw uses R's random-number stream exactly
// as stats::rWishart(1, n, diag(p)) does: it respects set.seed() and the per-chunk
// future.seed streams, and the tests compare it with stats::rWishart() under one seed.
// Also called by parallel_sim() in parallel.cpp.
void wishart_factor(arma::mat& U, const int n) {
  const arma::uword p = U.n_cols;
  for (arma::uword j = 0; j < p; j++) {
    U(j, j) = std::sqrt(R::rchisq(n - static_cast<double>(j)));
    for (arma::uword i = 0; i < j; i++) {
      U(i, j) = R::norm_rand();
    }
  }
}

// Turns a cross-product matrix W into its correlation matrix in place (as
// stats::cov2cor()). Each off-diagonal value is computed once and mirrored, so the result
// is exactly symmetric with an exactly unit diagonal. The diagonal holds 1 / sqrt(W_jj)
// while the off-diagonal values are scaled, so no workspace is needed. Also called by
// parallel_sim() in parallel.cpp.
void cov2cor_inplace(arma::mat& W) {
  const arma::uword p = W.n_cols;
  for (arma::uword j = 0; j < p; j++) {
    W(j, j) = 1.0 / std::sqrt(W(j, j));
  }
  for (arma::uword j = 0; j < p; j++) {
    for (arma::uword i = j + 1; i < p; i++) {
      W(j, i) = W(i, j) *= W(i, i) * W(j, j);
    }
  }
  W.diag().ones();
}

//' Draw multivariate-normal data from a population correlation matrix.
//'
//' Internal helper called from efa_simulate() and from the rank-correlation
//' reference simulations of efa_parallel(). Draws `N` cases from a `p`-variate
//' normal with correlation (or covariance) `R` by post-multiplying a matrix of
//' standard normal deviates by a matrix square root `M` of `R` (with M' M = R, so
//' the rows of Z * M are N(0, R)). A positive-definite `R` is factored by Cholesky;
//' a positive-semidefinite but singular `R` (which makes the Cholesky fail although
//' it is still a valid covariance, e.g. a no-factor block or a smoothed factor
//' intercorrelation matrix) falls back to a symmetric eigen square root.
//'
//' @param R numeric matrix. Population correlation/covariance matrix.
//' @param N integer. Number of cases to draw.
//' @param tol numeric. Eigenvalues below `-tol` mark `R` as indefinite.
//'
//' @keywords internal
// [[Rcpp::export(.simulate_cfm_mvn)]]
arma::mat simulate_cfm_mvn(const arma::mat& R, const int N,
                           const double tol = 1e-8) {

  const arma::mat M = cfm_sqrt(R, tol);

  // arma::randn draws column-major through R's RNG (RcppArmadillo's ARMA_RNG_ALT), so the
  // draw respects set.seed(). For an identity M the product Z * M is a no-op and only the
  // multiply is skipped: the same deviates are drawn in the same order.
  if (cfm_is_identity(M)) {
    return arma::randn(N, M.n_rows);
  }
  return arma::randn(N, M.n_rows) * M;
}

//' Reference eigenvalues for the efa_nest() simulation.
//'
//' Internal helper called from efa_nest(). Draws `nreps` correlation matrices with
//' the sampling distribution of the Pearson correlation matrix of `N` cases from
//' the `(nf - 1)`-factor reference model, given that model's loadings `Lambda` and
//' uniquenesses `Psi`, and returns the `nf`-th largest eigenvalue of each. The
//' matrices are drawn from the Wishart distribution with covariance
//' `Lambda Lambda' + diag(Psi)` and `N - 1` degrees of freedom and rescaled to
//' correlations, which is exact for normal data and does not simulate the data.
//'
//' @param nf integer. Position of the empirical eigenvalue being tested (1-based);
//'   the `nf`-th largest simulated eigenvalue is returned per replicate.
//' @param N integer. Number of cases / observations per simulated dataset. Must be
//'   larger than the number of variables.
//' @param Lambda numeric matrix. Loadings of the `(nf - 1)`-factor reference model
//'   (`p x (nf - 1)`); pass a `p x 0` matrix for the `nf == 1` null (identity) model.
//' @param Psi numeric vector. Uniquenesses (`1 - h2`) of the reference model.
//' @param nreps integer. Number of correlation matrices to draw.
//'
//' @keywords internal
// [[Rcpp::export(.simulate_cfm_eigen)]]
arma::vec simulate_cfm_eigen(const int nf, const int N, const arma::mat& Lambda,
                             const arma::vec& Psi, const int nreps = 1000) {

  // nf indexes into the ascending eigenvalues via ind = p - nf; a too-large nf would
  // make ind negative and read out of bounds, so reject it up front.
  const arma::uword p = Psi.n_elem;
  if (nf < 1 || static_cast<arma::uword>(nf) > p) {
    Rcpp::stop("nf must be between 1 and the number of variables.");
  }
  if (N <= static_cast<int>(p)) {
    Rcpp::stop("N must be larger than the number of variables.");
  }

  // Square root M (M' M = Sigma) of the reference covariance Sigma = Lambda Lambda' +
  // diag(Psi). For the nf == 1 null model Lambda has no columns, Sigma is the identity
  // and so is M; the identity test is made once, outside the replicate loop.
  arma::mat Sigma = Lambda * Lambda.t();
  Sigma.diag() += Psi;
  const arma::mat M = cfm_sqrt(Sigma, 1e-8);
  const bool m_is_identity = cfm_is_identity(M);

  // eig_sym returns eigenvalues in ascending order, so ind = p - nf is the nf-th
  // largest. The workspaces keep their size across replicates and are refilled in place.
  const arma::uword ind = p - nf;
  const int n = N - 1;
  arma::vec eigvals(p);
  arma::vec ref_values(nreps);
  arma::mat U(p, p, arma::fill::zeros);
  arma::mat UM(p, p);
  arma::mat W(p, p);

  for (arma::uword i = 0; i < static_cast<arma::uword>(nreps); i++) {
    wishart_factor(U, n);
    if (m_is_identity) {
      W = U.t() * U;
    } else {
      UM = U * M;
      W = UM.t() * UM;
    }
    cov2cor_inplace(W);
    eig_sym_values_checked(eigvals, W, "the reference simulation");
    ref_values(i) = eigvals(ind);
  }

  return ref_values;
}
