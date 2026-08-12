// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "eig_utils.h"

using namespace Rcpp;
using namespace arma;

// Exact test for the null / no-factor model: whether the matrix square root M is the
// identity (R = I in simulate_cfm_mvn, or the nf == 1 identity reference in
// simulate_cfm_eigen). It allocates a p x p identity and a p x p comparison, and M is
// fixed for a whole run, so each caller evaluates it once and passes the answer into
// the draw rather than repeating it per replicate.
static inline bool cfm_is_identity(const arma::mat& M) {
  return M.is_square() &&
         arma::all(arma::vectorise(M == arma::eye<arma::mat>(M.n_rows, M.n_cols)));
}

// Z * M draw shared by the population-MVN kernel (simulate_cfm_mvn) and the NEST
// reference simulation (simulate_cfm_eigen): N rows of standard normal deviates
// post-multiplied by a matrix square root M (with M' M equal to the target covariance),
// so each row of the result is N(0, M' M). arma::randn draws column-major through R's
// RNG (RcppArmadillo's ARMA_RNG_ALT), so the draw respects set.seed() and the
// per-replicate future.seed streams. `m_is_identity` comes from cfm_is_identity(M): the
// Z * M product is then a no-op and only the multiply is skipped, so the same
// N * M.n_rows deviates are drawn in the same order and the result stays byte-identical.
static inline arma::mat cfm_draw(const int N, const arma::mat& M,
                                 const bool m_is_identity) {
  if (m_is_identity) {
    return arma::randn(N, M.n_rows);
  }
  return arma::randn(N, M.n_rows) * M;
}

//' Draw multivariate-normal data from a population correlation matrix.
//'
//' Internal helper called from efa_simulate(). Draws `N` cases from a
//' `p`-variate normal with correlation (or covariance) `R` by post-multiplying a
//' matrix of standard normal deviates by a matrix square root `M` of `R` (with
//' M' M = R, so the rows of Z * M are N(0, R)). This is the same Z * M rule used by
//' the NEST reference simulation (.simulate_cfm_eigen): there `M` is the transposed
//' factor-score matrix, here it is a Cholesky or eigen square root.
//' A positive-definite `R` is factored by Cholesky; a positive-semidefinite but
//' singular `R` (which makes the Cholesky fail although it is still a valid
//' covariance, e.g. a no-factor block or a smoothed factor intercorrelation
//' matrix) falls back to a symmetric eigen square root.
//'
//' @param R numeric matrix. Population correlation/covariance matrix.
//' @param N integer. Number of cases to draw.
//' @param tol numeric. Eigenvalues below `-tol` mark `R` as indefinite.
// [[Rcpp::export(.simulate_cfm_mvn)]]
arma::mat simulate_cfm_mvn(const arma::mat& R, const int N,
                           const double tol = 1e-8) {

  // Matrix square root M with M' M = R. For a positive-definite R this is the
  // upper Cholesky factor (arma::chol returns M with M.t() * M == R); both the
  // chol output form and the eigen fallback size M themselves.
  arma::mat M;
  if (!arma::chol(M, R)) {

    // The Cholesky fails for a singular (but still positive-semidefinite) R.
    // Fall back to a symmetric eigen square root V diag(sqrt(d)) V', which
    // handles the zero-eigenvalue case; the symmetric root (rather than
    // diag(sqrt(d)) V') keeps the columns aligned with those of R.
    arma::vec eigval;
    arma::mat eigvec;
    if (!arma::eig_sym(eigval, eigvec, R)) {
      Rcpp::stop("Eigendecomposition of the population matrix failed; it is not "
                 "finite or not symmetric.");
    }
    // A markedly negative eigenvalue means R is indefinite and cannot be a
    // covariance to simulate from. efa_simulate() screens for this and raises a
    // classed error, so this branch is only a defensive backstop.
    if (eigval.min() < -tol) {
      Rcpp::stop("The population matrix is not positive semi-definite; it cannot "
                 "be a covariance to simulate from.");
    }
    eigval.transform([](double v) { return v < 0.0 ? 0.0 : std::sqrt(v); });
    M = eigvec * arma::diagmat(eigval) * eigvec.t();
  }

  return cfm_draw(N, M, cfm_is_identity(M));
}

//' Reference eigenvalues for the efa_nest() simulation via the shared kernel.
//'
//' Internal helper called from efa_nest(). Simulates `nreps` datasets from an
//' `(nf - 1)`-factor reference model, given that model's loadings `Lambda` and
//' uniquenesses `Psi`, and returns the `nf`-th largest eigenvalue of each simulated
//' correlation matrix. The data are drawn with the shared Z * M rule (see
//' .simulate_cfm_mvn()) using the factor-score square root
//' `M = t([Lambda | diag(sqrt(Psi))])`, so a row `randn(1, nf - 1 + p) * M` is
//' N(0, Lambda Lambda' + diag(Psi)). Drawing `nf - 1 + p` standard normals and
//' post-multiplying by the factor-score matrix is faster than forming the model-
//' implied correlation matrix and drawing from it, and matches the position at which
//' efa_nest() reads the reference eigenvalue.
//'
//' @param nf integer. Position of the empirical eigenvalue being tested (1-based);
//'   the `nf`-th largest simulated eigenvalue is returned per replicate.
//' @param N integer. Number of cases / observations per simulated dataset.
//' @param Lambda numeric matrix. Loadings of the `(nf - 1)`-factor reference model
//'   (`p x (nf - 1)`); pass a `p x 0` matrix for the `nf == 1` null (identity) model.
//' @param Psi numeric vector. Uniquenesses (`1 - h2`) of the reference model.
//' @param nreps integer. Number of datasets to simulate.
// [[Rcpp::export(.simulate_cfm_eigen)]]
arma::vec simulate_cfm_eigen(const int nf, const int N, const arma::mat& Lambda,
                             const arma::vec& Psi, const int nreps = 1000) {

  // nf indexes into the ascending eigenvalues via ind = p - nf; a too-large nf would
  // make ind negative and read out of bounds, so reject it up front.
  const arma::uword p = Psi.n_elem;
  if (nf < 1 || static_cast<arma::uword>(nf) > p) {
    Rcpp::stop("nf must be between 1 and the number of variables.");
  }

  // Factor-score square root M = t([Lambda | diag(sqrt(Psi))]). For the nf == 1 null
  // model Lambda has no columns, so M is just diag(sqrt(Psi)) -- guard the zero-column
  // case rather than relying on join_rows() accepting an empty first argument.
  arma::mat D = arma::diagmat(arma::sqrt(Psi));
  arma::mat M = (Lambda.n_cols == 0 ? D : arma::join_rows(Lambda, D)).t();
  // M does not change inside the replicate loop, so the identity test is made once
  const bool m_is_identity = cfm_is_identity(M);

  // eig_sym returns eigenvalues in ascending order, so ind = p - nf is the nf-th
  // largest.
  const arma::uword ind = p - nf;
  arma::vec eigvals(p);
  arma::vec ref_values(nreps);

  for (arma::uword i = 0; i < static_cast<arma::uword>(nreps); i++) {
    eig_sym_values_checked(eigvals, arma::cor(cfm_draw(N, M, m_is_identity)),
                           "the reference simulation");
    ref_values(i) = eigvals(ind);
  }

  return ref_values;
}
