#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include "gpf_engine.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

using namespace Rcpp;
using namespace arma;

// Gradient-projection factor rotation along the orthogonal (Stiefel) manifold,
// generalized over an arbitrary rotation criterion (Bernaards & Jennrich, 2005). A
// criterion exposes the objective value f and its gradient Gq = dQ/dL at the rotated
// loadings L; the engine maps Gq to the gradient with respect to the orthogonal
// transformation T, projects it onto the tangent space, line-searches with a
// sufficient-decrease test, and retracts back onto the orthogonal group via a polar
// (SVD) projection. This mirrors GPArotation::GPForth(); only the criterion functor
// changes between rotation methods.

// A rotation criterion: value f and gradient Gq = dQ/dL at the rotated loadings L.
struct RotationCriterion {
  virtual ~RotationCriterion() {}
  virtual void eval(const arma::mat& L, double& f, arma::mat& Gq) const = 0;
};

// Crawford-Ferguson family, parameterized by kappa (GPArotation::vgQ.cf). With
// N = 1_k 1_k' - I_k and M = 1_p 1_p' - I_p the criterion balances row complexity
// (between-factor cross-products of squared loadings) against column complexity
// (between-variable cross-products). kappa = 0 is quartimax; kappa = k / (2 p) is
// equamax.
struct CfCriterion : public RotationCriterion {
  double kappa_;
  explicit CfCriterion(double kappa) : kappa_(kappa) {}

  void eval(const arma::mat& L, double& f, arma::mat& Gq) const override {
    // With N = 1_k 1_k' - I_k and M = 1_p 1_p' - I_p, the products L2 %*% N and
    // M %*% L2 are sum-broadcasts of L2: (L2 N)(i,j) = rowSum_i(L2) - L2(i,j) and
    // (M L2)(i,j) = colSum_j(L2) - L2(i,j). Forming them from the row/column sums
    // avoids building the p x p / k x k constant matrices and the dense products on
    // every evaluation (this is the innermost kernel of the multi-start solver).
    arma::mat L2 = square(L);
    arma::vec row_sums = sum(L2, 1);
    arma::rowvec col_sums = sum(L2, 0);
    arma::mat L2N = -L2;
    L2N.each_col() += row_sums;
    arma::mat ML2 = -L2;
    ML2.each_row() += col_sums;
    double f1 = (1.0 - kappa_) * accu(L2 % L2N) / 4.0;
    double f2 = kappa_ * accu(L2 % ML2) / 4.0;
    f = f1 + f2;
    Gq = (1.0 - kappa_) * (L % L2N) + kappa_ * (L % ML2);
  }
};

// Oblimin family, parameterized by gamma (GPArotation::vgQ.oblimin). With
// N = 1_k 1_k' - I_k the criterion penalizes the between-factor cross-products of squared
// loadings; gamma centers those cross-products across variables. gamma = 0 is the
// quartimin criterion. As in CfCriterion the dense N / constant matrices are avoided: the
// product X = L2 %*% N is the sum-broadcast X(i,j) = rowSum_i(L2) - L2(i,j), and when
// gamma != 0 the premultiplication by (I_p - (gamma / p) 1_p 1_p') subtracts the
// (gamma / p)-scaled column sums of X.
struct ObliminCriterion : public RotationCriterion {
  double gamma_;
  explicit ObliminCriterion(double gamma) : gamma_(gamma) {}

  void eval(const arma::mat& L, double& f, arma::mat& Gq) const override {
    const double p = static_cast<double>(L.n_rows);
    arma::mat L2 = square(L);
    arma::vec row_sums = sum(L2, 1);
    arma::mat X = -L2;
    X.each_col() += row_sums;
    if (gamma_ != 0.0) {
      arma::rowvec col_sums = sum(X, 0);
      X.each_row() -= (gamma_ / p) * col_sums;
    }
    f = accu(L2 % X) / 4.0;
    Gq = L % X;
  }
};

// Geomin family, parameterized by delta (GPArotation::vgQ.geomin). The criterion is the sum
// over variables of the geometric mean of their squared loadings, each offset by delta:
// f = sum_i exp(mean_j log(L2_ij)) with L2 = L^2 + delta. The offset keeps the row geometric
// means (and their logarithms) finite when a loading is zero; smaller delta rewards a sparser
// pattern but sharpens the local minima. delta = 0.01 is the GPArotation default. With the
// per-row geometric mean pro_i = exp(mean_j log(L2_ij)), the gradient is (2 / k)(L / L2) with
// each row scaled by pro_i.
struct GeominCriterion : public RotationCriterion {
  double delta_;
  explicit GeominCriterion(double delta) : delta_(delta) {}

  void eval(const arma::mat& L, double& f, arma::mat& Gq) const override {
    const double k = static_cast<double>(L.n_cols);
    // delta > 0 (validated by the entry points) keeps L2 strictly positive, so log(L2) and
    // the geometric means below are always finite.
    arma::mat L2 = square(L);
    L2 += delta_;
    arma::vec pro = exp(mean(log(L2), 1));  // dim = 1: per-row mean of log(L2) (rowMeans)
    f = accu(pro);
    Gq = (2.0 / k) * (L / L2);
    Gq.each_col() %= pro;  // scale row i by pro_i (pro broadcast across columns)
  }
};

// Bentler's invariant pattern simplicity criterion (GPArotation::vgQ.bentler). With the k x k
// cross-product of squared loadings M = L2' L2 and its diagonal D, the criterion is
// f = -(log|M| - log|D|) / 4: it is minimized when M is as close to diagonal as possible, i.e.
// the columns of squared loadings are mutually orthogonal (an invariant-pattern notion of simple
// structure). The gradient is Gq = -L .* (L2 (M^-1 - D^-1)). A degenerate M -- a rank-deficient
// cross-product or a zero column of squared loadings (a zero diagonal entry) -- leaves the
// criterion undefined; this is signalled by an infinite objective and a zero gradient so the
// manifold's validity check rejects the candidate, mirroring the invalid-state handling of the
// other criteria rather than throwing.
struct BentlerCriterion : public RotationCriterion {
  void eval(const arma::mat& L, double& f, arma::mat& Gq) const override {
    arma::mat L2 = square(L);
    arma::mat M = L2.t() * L2;  // k x k cross-product of squared loadings
    arma::vec dM = M.diag();

    double logDetM, sgn;
    bool det_ok = arma::log_det(logDetM, sgn, M);  // logDetM = log|det(M)|, sgn = sign(det(M))
    arma::mat Minv;
    bool inv_ok = inverse_checked_cpp(M, Minv);

    if (!det_ok || sgn <= 0.0 || !std::isfinite(logDetM) || !inv_ok || dM.min() <= 0.0) {
      f = std::numeric_limits<double>::infinity();
      Gq.zeros(L.n_rows, L.n_cols);
      return;
    }

    double logDetD = accu(log(dM));
    f = -(logDetM - logDetD) / 4.0;
    Gq = -(L % (L2 * (Minv - diagmat(1.0 / dM))));
  }
};


// Retract X back onto the orthogonal group via its polar factor U V' from the SVD,
// the Stiefel retraction used by GPForth. Returns false (leaving Tout untouched) if
// the decomposition fails or yields a non-finite matrix, so the candidate is rejected.
static bool retract_orth(const arma::mat& X, arma::mat& Tout) {
  arma::mat U, V;
  arma::vec s;
  if (!arma::svd(U, s, V, X)) {
    return false;
  }
  Tout = U * V.t();
  return all_finite_cpp(Tout);
}

// Orthogonal rotation manifold for the shared gradient-projection engine. The
// transformation T lives on the Stiefel manifold (T'T = I) and the rotated loadings
// are L = A %*% T. The criterion supplies the objective f and its gradient Gq = dQ/dL;
// this maps Gq to the gradient with respect to T (G = A' Gq), projects it onto the
// tangent space (Gp = G - T (T'G + G'T)/2), and reports its Frobenius norm. Candidates
// are retracted back onto the orthogonal group by the polar (SVD) projection. A
// non-finite criterion or gradient yields an infinite (invalid) objective.
struct OrthCriterionManifold {
  const arma::mat& A;
  const RotationCriterion& crit;

  GpfState compute(const arma::mat& Tmat) const {
    arma::mat L = A * Tmat;
    if (!all_finite_cpp(L)) {
      return gpf_invalid_state();
    }

    double f;
    arma::mat Gq;
    crit.eval(L, f, Gq);

    arma::mat G = A.t() * Gq;
    arma::mat M = Tmat.t() * G;
    arma::mat S = (M + M.t()) / 2.0;
    arma::mat Gp = G - Tmat * S;
    double s = norm(vectorise(Gp), 2);

    GpfState out;
    out.f = f;
    out.s = s;
    out.Gp = Gp;
    out.valid = std::isfinite(f) && std::isfinite(s) && all_finite_cpp(Gq) &&
      all_finite_cpp(G) && all_finite_cpp(Gp);

    if (!out.valid) {
      return gpf_invalid_state();
    }
    return out;
  }

  bool retract(const arma::mat& X, arma::mat& T_out) const {
    return retract_orth(X, T_out);
  }

  // Orthogonal starts (identity and random_orthogonal_start_cpp) are already on the
  // manifold, so the start needs no further conditioning.
  arma::mat normalize_start(const arma::mat& Tmat) const {
    return Tmat;
  }
};

// Oblique rotation manifold for the shared gradient-projection engine. The transformation
// T lives on the column-normalized manifold diag(T'T) = 1 and the rotated loadings are
// L = A %*% solve(t(T)). The criterion supplies the objective f and its gradient
// Gq = dQ/dL; this maps Gq to the gradient with respect to T
// (G = -t(t(L) %*% Gq %*% solve(T))), projects it onto the tangent space
// (Gp = G - T diag(colSums(T .* G))), and reports its Frobenius norm. Candidates are
// retracted back onto the manifold by column normalization. This mirrors GPArotation's
// GPFoblq(); only the criterion functor changes between rotation methods. A non-invertible
// T or a non-finite criterion/gradient yields an infinite (invalid) objective.
struct OblqCriterionManifold {
  const arma::mat& A;
  const RotationCriterion& crit;

  GpfState compute(const arma::mat& Tmat) const {
    arma::mat invT;
    if (!inverse_checked_cpp(Tmat, invT)) {
      return gpf_invalid_state();
    }
    arma::mat L = A * invT.t();
    if (!all_finite_cpp(L)) {
      return gpf_invalid_state();
    }

    double f;
    arma::mat Gq;
    crit.eval(L, f, Gq);

    arma::mat G = -(L.t() * Gq * invT).t();
    arma::mat Gp;
    double s;
    oblique_tangent_project(Tmat, G, Gp, s);

    GpfState out;
    out.f = f;
    out.s = s;
    out.Gp = Gp;
    out.valid = std::isfinite(f) && std::isfinite(s) && all_finite_cpp(Gq) &&
      all_finite_cpp(G) && all_finite_cpp(Gp);

    if (!out.valid) {
      return gpf_invalid_state();
    }
    return out;
  }

  // The column-normalization projection always returns a candidate on the manifold; a
  // singular candidate is rejected later by the inverse check in compute().
  bool retract(const arma::mat& X, arma::mat& T_out) const {
    T_out = normalize_cols_cpp(X);
    return true;
  }

  arma::mat normalize_start(const arma::mat& Tmat) const {
    return normalize_cols_cpp(Tmat);
  }
};

// Kaiser row weights for a single loading matrix, matching GPArotation's
// NormalizingWeight(): each row weight is max(sqrt(rowSumSq), machine eps). Dividing
// the rows by these weights before rotation and restoring afterwards is Kaiser
// normalization.
static arma::vec kaiser_weights(const arma::mat& A) {
  arma::vec W = sqrt(sum(square(A), 1));
  const double floor = std::numeric_limits<double>::epsilon();
  for (arma::uword i = 0; i < W.n_elem; ++i) {
    if (!std::isfinite(W(i)) || W(i) < floor) {
      W(i) = floor;
    }
  }
  return W;
}

//' Orthogonal Crawford-Ferguson factor rotation
//'
//' Rotate a loading matrix orthogonally under the Crawford-Ferguson criterion using a
//' gradient-projection optimizer along the orthogonal (Stiefel) manifold.
//'
//' The criterion value `f` and its gradient `dQ/dL` at the rotated loadings
//' `L = A %*% T` define the search; the engine maps the gradient to the orthogonal
//' transformation `T`, projects it onto the tangent space, performs a
//' sufficient-decrease line search, and retracts back onto the orthogonal group via a
//' polar (singular value) projection. `kappa = 0` is the quartimax criterion and
//' `kappa = ncol(A) / (2 * nrow(A))` is the equamax criterion.
//'
//' Additional random orthogonal starts may be requested. To bound runtime the solver
//' screens each random start by its objective, runs a short triage optimization on the
//' best-screened starts, and fully optimizes only those that improve on the current
//' incumbent by at least `triage_improve_tol`.
//'
//' @param L Numeric matrix. The unrotated loading matrix (variables by factors).
//' @param kappa Numeric scalar in `[0, 1]`. The Crawford-Ferguson parameter.
//' @param eps Numeric scalar. Convergence tolerance for the projected-gradient norm.
//' @param normalize Logical scalar. If `TRUE`, apply Kaiser normalization before
//'   rotation and reverse it afterwards.
//' @param random_starts Integer scalar. Number of additional random orthogonal starts.
//' @param maxit Integer scalar. Maximum number of projected-gradient updates.
//' @param max_line_search Integer scalar. Maximum number of step-halving attempts after
//'   the initial trial step in each line-search phase.
//' @param step0 Numeric scalar. Initial step size used in the projected-gradient update.
//' @param screen_keep Integer scalar. Number of screened random starts retained for
//'   triage optimization.
//' @param triage_maxit Integer scalar. Number of short optimization iterations used in
//'   the triage stage.
//' @param triage_improve_tol Numeric scalar. Relative improvement required for a triaged
//'   start to be promoted to full optimization.
//'
//' @returns A named list with the rotated loadings, the orthogonal rotation matrix `Th`
//'   (with `L %*% Th` reproducing the rotated loadings), the attained criterion value, and
//'   the convergence and validity flags.
//'
//' @references
//' Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms and
//' software for arbitrary rotation criteria in factor analysis. *Educational and
//' Psychological Measurement*, 65, 676-696.
//'
//' Crawford, C. B., & Ferguson, G. A. (1970). A general rotation criterion and its use
//' in orthogonal rotation. *Psychometrika*, 35, 321-332.
//'
// [[Rcpp::export(.rotate_cf_orth)]]
Rcpp::List rotate_cf_orth(const arma::mat& L,
                          double kappa,
                          double eps = 1e-5,
                          bool normalize = true,
                          int random_starts = 0,
                          int maxit = 1000,
                          int max_line_search = 10,
                          double step0 = 1.0,
                          int screen_keep = 2,
                          int triage_maxit = 25,
                          double triage_improve_tol = 0.0) {
  if (L.n_rows == 0 || L.n_cols == 0) {
    Rcpp::stop("L must be a non-empty matrix.");
  }
  if (L.n_cols < 2) {
    Rcpp::stop("Orthogonal rotation requires at least two factors; use the R wrapper "
               "for single-factor solutions.");
  }
  if (!all_finite_cpp(L)) {
    Rcpp::stop("L must contain only finite values.");
  }
  if (!is_valid_scalar_cpp(kappa) || kappa < 0.0 || kappa > 1.0) {
    Rcpp::stop("kappa must be a finite scalar in [0, 1].");
  }
  validate_gpf_scalars(eps, maxit, max_line_search, step0, random_starts,
                       screen_keep, triage_maxit, triage_improve_tol);

  const arma::uword k = L.n_cols;

  // Kaiser normalization scales rows before rotation; because the weights act on rows
  // and the rotation right-multiplies, they cancel in the final loadings, so the
  // rotated loadings are reported as L %*% Th (the documented reproduction identity).
  arma::mat A_work = L;
  arma::vec W;
  if (normalize) {
    W = kaiser_weights(L);
    A_work.each_col() /= W;
  }

  CfCriterion crit(kappa);
  OrthCriterionManifold manifold{A_work, crit};
  arma::mat T_primary(k, k, fill::eye);

  GpfSummary summary = run_gpf_multistart(
    manifold, T_primary, eps, maxit, max_line_search, step0,
    random_starts, screen_keep, triage_maxit, triage_improve_tol
  );
  const GpfFit& best_fit = summary.best_fit;

  arma::mat loadings = L * best_fit.T;

  return Rcpp::List::create(
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("Th") = best_fit.T,
    Rcpp::Named("value") = best_fit.value,
    Rcpp::Named("convergence") = best_fit.convergence,
    Rcpp::Named("valid") = best_fit.valid
  );
}

//' Oblique oblimin factor rotation
//'
//' Rotate a loading matrix obliquely under the oblimin criterion using a
//' gradient-projection optimizer along the oblique (column-normalized) manifold.
//'
//' The criterion value `f` and its gradient `dQ/dL` at the rotated loadings
//' `L = A %*% solve(t(T))` define the search; the engine maps the gradient to the
//' transformation `T` on the manifold `diag(t(T) %*% T) = 1`, projects it onto the tangent
//' space, performs a sufficient-decrease line search, and retracts back onto the manifold
//' by column normalization. `gam = 0` is the quartimin criterion.
//'
//' Additional random starts may be requested. To bound runtime the solver screens each
//' random start by its objective, runs a short triage optimization on the best-screened
//' starts, and fully optimizes only those that improve on the current incumbent by at
//' least `triage_improve_tol`.
//'
//' @param L Numeric matrix. The unrotated loading matrix (variables by factors).
//' @param gam Numeric scalar. The oblimin parameter; `gam = 0` is the quartimin criterion.
//' @param eps Numeric scalar. Convergence tolerance for the projected-gradient norm.
//' @param normalize Logical scalar. If `TRUE`, apply Kaiser normalization before rotation
//'   and reverse it afterwards.
//' @param random_starts Integer scalar. Number of additional random starts.
//' @param maxit Integer scalar. Maximum number of projected-gradient updates.
//' @param max_line_search Integer scalar. Maximum number of step-halving attempts after
//'   the initial trial step in each line-search phase.
//' @param step0 Numeric scalar. Initial step size used in the projected-gradient update.
//' @param screen_keep Integer scalar. Number of screened random starts retained for triage
//'   optimization.
//' @param triage_maxit Integer scalar. Number of short optimization iterations used in the
//'   triage stage.
//' @param triage_improve_tol Numeric scalar. Relative improvement required for a triaged
//'   start to be promoted to full optimization.
//'
//' @returns A named list with the rotated loadings, the transformation matrix `Th`
//'   (with `L %*% t(solve(Th))` reproducing the rotated loadings), the factor correlation
//'   matrix `Phi` (`t(Th) %*% Th`), the attained criterion value, and the convergence and
//'   validity flags.
//'
//' @references
//' Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms and
//' software for arbitrary rotation criteria in factor analysis. *Educational and
//' Psychological Measurement*, 65, 676-696.
//'
//' Jennrich, R. I., & Sampson, P. F. (1966). Rotation for simple loadings.
//' *Psychometrika*, 31, 313-323.
//'
// [[Rcpp::export(.rotate_oblimin)]]
Rcpp::List rotate_oblimin(const arma::mat& L,
                          double gam = 0.0,
                          double eps = 1e-5,
                          bool normalize = true,
                          int random_starts = 0,
                          int maxit = 1000,
                          int max_line_search = 10,
                          double step0 = 1.0,
                          int screen_keep = 2,
                          int triage_maxit = 25,
                          double triage_improve_tol = 0.0) {
  if (L.n_rows == 0 || L.n_cols == 0) {
    Rcpp::stop("L must be a non-empty matrix.");
  }
  if (L.n_cols < 2) {
    Rcpp::stop("Oblique rotation requires at least two factors; use the R wrapper "
               "for single-factor solutions.");
  }
  if (!all_finite_cpp(L)) {
    Rcpp::stop("L must contain only finite values.");
  }
  if (!is_valid_scalar_cpp(gam)) {
    Rcpp::stop("gam must be a finite scalar.");
  }
  validate_gpf_scalars(eps, maxit, max_line_search, step0, random_starts,
                       screen_keep, triage_maxit, triage_improve_tol);

  const arma::uword k = L.n_cols;

  // Kaiser normalization scales rows before rotation; because the weights act on rows and
  // the rotation (right-multiplying solve(t(T))) acts on columns, they cancel in the final
  // loadings, so the rotated loadings reproduce L %*% t(solve(Th)) (the documented
  // reproduction identity).
  arma::mat A_work = L;
  arma::vec W;
  if (normalize) {
    W = kaiser_weights(L);
    A_work.each_col() /= W;
  }

  ObliminCriterion crit(gam);
  OblqCriterionManifold manifold{A_work, crit};
  arma::mat T_primary(k, k, fill::eye);

  GpfSummary summary = run_gpf_multistart(
    manifold, T_primary, eps, maxit, max_line_search, step0,
    random_starts, screen_keep, triage_maxit, triage_improve_tol
  );
  const GpfFit& best_fit = summary.best_fit;

  // Reconstruct the rotated loadings (L = A %*% solve(t(T)), Kaiser weights restored) and
  // the factor correlations (Phi = t(T) %*% T) from the winning transformation, shared with
  // the Procrustes oblique entry point.
  arma::mat loadings;
  arma::mat Phi;
  finalize_oblique(best_fit, A_work, W, normalize, loadings, Phi);

  return Rcpp::List::create(
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("Th") = best_fit.T,
    Rcpp::Named("Phi") = Phi,
    Rcpp::Named("value") = best_fit.value,
    Rcpp::Named("convergence") = best_fit.convergence,
    Rcpp::Named("valid") = best_fit.valid
  );
}

//' Orthogonal geomin factor rotation
//'
//' Rotate a loading matrix orthogonally under the geomin criterion using a
//' gradient-projection optimizer along the orthogonal (Stiefel) manifold.
//'
//' The criterion value `f` and its gradient `dQ/dL` at the rotated loadings
//' `L = A %*% T` define the search; the engine maps the gradient to the orthogonal
//' transformation `T`, projects it onto the tangent space, performs a sufficient-decrease
//' line search, and retracts back onto the orthogonal group via a polar (singular value)
//' projection. The geomin criterion sums the per-variable geometric mean of the squared
//' loadings offset by `delta`; it is prone to local minima, so additional random starts are
//' recommended.
//'
//' Additional random orthogonal starts may be requested. To bound runtime the solver screens
//' each random start by its objective, runs a short triage optimization on the best-screened
//' starts, and fully optimizes only those that improve on the current incumbent by at least
//' `triage_improve_tol`.
//'
//' @param L Numeric matrix. The unrotated loading matrix (variables by factors).
//' @param delta Numeric scalar. The geomin offset added to the squared loadings; must be a
//'   positive finite scalar. `delta = 0.01` is the usual default.
//' @param eps Numeric scalar. Convergence tolerance for the projected-gradient norm.
//' @param normalize Logical scalar. If `TRUE`, apply Kaiser normalization before rotation and
//'   reverse it afterwards.
//' @param random_starts Integer scalar. Number of additional random orthogonal starts.
//' @param maxit Integer scalar. Maximum number of projected-gradient updates.
//' @param max_line_search Integer scalar. Maximum number of step-halving attempts after the
//'   initial trial step in each line-search phase.
//' @param step0 Numeric scalar. Initial step size used in the projected-gradient update.
//' @param screen_keep Integer scalar. Number of screened random starts retained for triage
//'   optimization.
//' @param triage_maxit Integer scalar. Number of short optimization iterations used in the
//'   triage stage.
//' @param triage_improve_tol Numeric scalar. Relative improvement required for a triaged start
//'   to be promoted to full optimization.
//'
//' @returns A named list with the rotated loadings, the orthogonal rotation matrix `Th`
//'   (with `L %*% Th` reproducing the rotated loadings), the attained criterion value, and the
//'   convergence and validity flags.
//'
//' @references
//' Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms and
//' software for arbitrary rotation criteria in factor analysis. *Educational and
//' Psychological Measurement*, 65, 676-696.
//'
//' Browne, M. W. (2001). An overview of analytic rotation in exploratory factor analysis.
//' *Multivariate Behavioral Research*, 36, 111-150.
//'
// [[Rcpp::export(.rotate_geomin_orth)]]
Rcpp::List rotate_geomin_orth(const arma::mat& L,
                              double delta = 0.01,
                              double eps = 1e-5,
                              bool normalize = true,
                              int random_starts = 0,
                              int maxit = 1000,
                              int max_line_search = 10,
                              double step0 = 1.0,
                              int screen_keep = 2,
                              int triage_maxit = 25,
                              double triage_improve_tol = 0.0) {
  if (L.n_rows == 0 || L.n_cols == 0) {
    Rcpp::stop("L must be a non-empty matrix.");
  }
  if (L.n_cols < 2) {
    Rcpp::stop("Orthogonal rotation requires at least two factors; use the R wrapper "
               "for single-factor solutions.");
  }
  if (!all_finite_cpp(L)) {
    Rcpp::stop("L must contain only finite values.");
  }
  if (!is_valid_scalar_cpp(delta) || delta <= 0.0) {
    Rcpp::stop("delta must be a positive finite scalar.");
  }
  validate_gpf_scalars(eps, maxit, max_line_search, step0, random_starts,
                       screen_keep, triage_maxit, triage_improve_tol);

  const arma::uword k = L.n_cols;

  // Kaiser normalization scales rows before rotation; because the weights act on rows and the
  // rotation right-multiplies, they cancel in the final loadings, so the rotated loadings are
  // reported as L %*% Th (the documented reproduction identity).
  arma::mat A_work = L;
  arma::vec W;
  if (normalize) {
    W = kaiser_weights(L);
    A_work.each_col() /= W;
  }

  GeominCriterion crit(delta);
  OrthCriterionManifold manifold{A_work, crit};
  arma::mat T_primary(k, k, fill::eye);

  GpfSummary summary = run_gpf_multistart(
    manifold, T_primary, eps, maxit, max_line_search, step0,
    random_starts, screen_keep, triage_maxit, triage_improve_tol
  );
  const GpfFit& best_fit = summary.best_fit;

  arma::mat loadings = L * best_fit.T;

  return Rcpp::List::create(
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("Th") = best_fit.T,
    Rcpp::Named("value") = best_fit.value,
    Rcpp::Named("convergence") = best_fit.convergence,
    Rcpp::Named("valid") = best_fit.valid
  );
}

//' Oblique geomin factor rotation
//'
//' Rotate a loading matrix obliquely under the geomin criterion using a gradient-projection
//' optimizer along the oblique (column-normalized) manifold.
//'
//' The criterion value `f` and its gradient `dQ/dL` at the rotated loadings
//' `L = A %*% solve(t(T))` define the search; the engine maps the gradient to the
//' transformation `T` on the manifold `diag(t(T) %*% T) = 1`, projects it onto the tangent
//' space, performs a sufficient-decrease line search, and retracts back onto the manifold by
//' column normalization. The geomin criterion sums the per-variable geometric mean of the
//' squared loadings offset by `delta`; it is prone to local minima, so additional random
//' starts are recommended.
//'
//' Additional random starts may be requested. To bound runtime the solver screens each random
//' start by its objective, runs a short triage optimization on the best-screened starts, and
//' fully optimizes only those that improve on the current incumbent by at least
//' `triage_improve_tol`.
//'
//' @param L Numeric matrix. The unrotated loading matrix (variables by factors).
//' @param delta Numeric scalar. The geomin offset added to the squared loadings; must be a
//'   positive finite scalar. `delta = 0.01` is the usual default.
//' @param eps Numeric scalar. Convergence tolerance for the projected-gradient norm.
//' @param normalize Logical scalar. If `TRUE`, apply Kaiser normalization before rotation and
//'   reverse it afterwards.
//' @param random_starts Integer scalar. Number of additional random starts.
//' @param maxit Integer scalar. Maximum number of projected-gradient updates.
//' @param max_line_search Integer scalar. Maximum number of step-halving attempts after the
//'   initial trial step in each line-search phase.
//' @param step0 Numeric scalar. Initial step size used in the projected-gradient update.
//' @param screen_keep Integer scalar. Number of screened random starts retained for triage
//'   optimization.
//' @param triage_maxit Integer scalar. Number of short optimization iterations used in the
//'   triage stage.
//' @param triage_improve_tol Numeric scalar. Relative improvement required for a triaged start
//'   to be promoted to full optimization.
//'
//' @returns A named list with the rotated loadings, the transformation matrix `Th`
//'   (with `L %*% t(solve(Th))` reproducing the rotated loadings), the factor correlation
//'   matrix `Phi` (`t(Th) %*% Th`), the attained criterion value, and the convergence and
//'   validity flags.
//'
//' @references
//' Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms and
//' software for arbitrary rotation criteria in factor analysis. *Educational and
//' Psychological Measurement*, 65, 676-696.
//'
//' Browne, M. W. (2001). An overview of analytic rotation in exploratory factor analysis.
//' *Multivariate Behavioral Research*, 36, 111-150.
//'
// [[Rcpp::export(.rotate_geomin_oblq)]]
Rcpp::List rotate_geomin_oblq(const arma::mat& L,
                              double delta = 0.01,
                              double eps = 1e-5,
                              bool normalize = true,
                              int random_starts = 0,
                              int maxit = 1000,
                              int max_line_search = 10,
                              double step0 = 1.0,
                              int screen_keep = 2,
                              int triage_maxit = 25,
                              double triage_improve_tol = 0.0) {
  if (L.n_rows == 0 || L.n_cols == 0) {
    Rcpp::stop("L must be a non-empty matrix.");
  }
  if (L.n_cols < 2) {
    Rcpp::stop("Oblique rotation requires at least two factors; use the R wrapper "
               "for single-factor solutions.");
  }
  if (!all_finite_cpp(L)) {
    Rcpp::stop("L must contain only finite values.");
  }
  if (!is_valid_scalar_cpp(delta) || delta <= 0.0) {
    Rcpp::stop("delta must be a positive finite scalar.");
  }
  validate_gpf_scalars(eps, maxit, max_line_search, step0, random_starts,
                       screen_keep, triage_maxit, triage_improve_tol);

  const arma::uword k = L.n_cols;

  // Kaiser normalization scales rows before rotation; because the weights act on rows and the
  // rotation (right-multiplying solve(t(T))) acts on columns, they cancel in the final
  // loadings, so the rotated loadings reproduce L %*% t(solve(Th)) (the documented
  // reproduction identity).
  arma::mat A_work = L;
  arma::vec W;
  if (normalize) {
    W = kaiser_weights(L);
    A_work.each_col() /= W;
  }

  GeominCriterion crit(delta);
  OblqCriterionManifold manifold{A_work, crit};
  arma::mat T_primary(k, k, fill::eye);

  GpfSummary summary = run_gpf_multistart(
    manifold, T_primary, eps, maxit, max_line_search, step0,
    random_starts, screen_keep, triage_maxit, triage_improve_tol
  );
  const GpfFit& best_fit = summary.best_fit;

  // Reconstruct the rotated loadings (L = A %*% solve(t(T)), Kaiser weights restored) and the
  // factor correlations (Phi = t(T) %*% T) from the winning transformation, shared with the
  // Procrustes oblique entry point.
  arma::mat loadings;
  arma::mat Phi;
  finalize_oblique(best_fit, A_work, W, normalize, loadings, Phi);

  return Rcpp::List::create(
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("Th") = best_fit.T,
    Rcpp::Named("Phi") = Phi,
    Rcpp::Named("value") = best_fit.value,
    Rcpp::Named("convergence") = best_fit.convergence,
    Rcpp::Named("valid") = best_fit.valid
  );
}

//' Orthogonal Bentler factor rotation
//'
//' Rotate a loading matrix orthogonally under Bentler's invariant pattern simplicity criterion
//' using a gradient-projection optimizer along the orthogonal (Stiefel) manifold.
//'
//' The criterion value `f` and its gradient `dQ/dL` at the rotated loadings
//' `L = A %*% T` define the search; the engine maps the gradient to the orthogonal
//' transformation `T`, projects it onto the tangent space, performs a sufficient-decrease line
//' search, and retracts back onto the orthogonal group via a polar (singular value) projection.
//' The Bentler criterion measures the departure of the cross-products of squared loadings from a
//' diagonal pattern; it is prone to local minima, so additional random starts are recommended.
//'
//' Additional random orthogonal starts may be requested. To bound runtime the solver screens
//' each random start by its objective, runs a short triage optimization on the best-screened
//' starts, and fully optimizes only those that improve on the current incumbent by at least
//' `triage_improve_tol`.
//'
//' @param L Numeric matrix. The unrotated loading matrix (variables by factors).
//' @param eps Numeric scalar. Convergence tolerance for the projected-gradient norm.
//' @param normalize Logical scalar. If `TRUE`, apply Kaiser normalization before rotation and
//'   reverse it afterwards.
//' @param random_starts Integer scalar. Number of additional random orthogonal starts.
//' @param maxit Integer scalar. Maximum number of projected-gradient updates.
//' @param max_line_search Integer scalar. Maximum number of step-halving attempts after the
//'   initial trial step in each line-search phase.
//' @param step0 Numeric scalar. Initial step size used in the projected-gradient update.
//' @param screen_keep Integer scalar. Number of screened random starts retained for triage
//'   optimization.
//' @param triage_maxit Integer scalar. Number of short optimization iterations used in the
//'   triage stage.
//' @param triage_improve_tol Numeric scalar. Relative improvement required for a triaged start
//'   to be promoted to full optimization.
//'
//' @returns A named list with the rotated loadings, the orthogonal rotation matrix `Th`
//'   (with `L %*% Th` reproducing the rotated loadings), the attained criterion value, and the
//'   convergence and validity flags.
//'
//' @references
//' Bentler, P. M. (1977). Factor simplicity index and transformations. *Psychometrika*, 42,
//' 277-295.
//'
//' Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms and
//' software for arbitrary rotation criteria in factor analysis. *Educational and
//' Psychological Measurement*, 65, 676-696.
//'
// [[Rcpp::export(.rotate_bentler_orth)]]
Rcpp::List rotate_bentler_orth(const arma::mat& L,
                               double eps = 1e-5,
                               bool normalize = true,
                               int random_starts = 0,
                               int maxit = 1000,
                               int max_line_search = 10,
                               double step0 = 1.0,
                               int screen_keep = 2,
                               int triage_maxit = 25,
                               double triage_improve_tol = 0.0) {
  if (L.n_rows == 0 || L.n_cols == 0) {
    Rcpp::stop("L must be a non-empty matrix.");
  }
  if (L.n_cols < 2) {
    Rcpp::stop("Orthogonal rotation requires at least two factors; use the R wrapper "
               "for single-factor solutions.");
  }
  if (!all_finite_cpp(L)) {
    Rcpp::stop("L must contain only finite values.");
  }
  validate_gpf_scalars(eps, maxit, max_line_search, step0, random_starts,
                       screen_keep, triage_maxit, triage_improve_tol);

  const arma::uword k = L.n_cols;

  // Kaiser normalization scales rows before rotation; because the weights act on rows and the
  // rotation right-multiplies, they cancel in the final loadings, so the rotated loadings are
  // reported as L %*% Th (the documented reproduction identity).
  arma::mat A_work = L;
  arma::vec W;
  if (normalize) {
    W = kaiser_weights(L);
    A_work.each_col() /= W;
  }

  BentlerCriterion crit;
  OrthCriterionManifold manifold{A_work, crit};
  arma::mat T_primary(k, k, fill::eye);

  GpfSummary summary = run_gpf_multistart(
    manifold, T_primary, eps, maxit, max_line_search, step0,
    random_starts, screen_keep, triage_maxit, triage_improve_tol
  );
  const GpfFit& best_fit = summary.best_fit;

  arma::mat loadings = L * best_fit.T;

  return Rcpp::List::create(
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("Th") = best_fit.T,
    Rcpp::Named("value") = best_fit.value,
    Rcpp::Named("convergence") = best_fit.convergence,
    Rcpp::Named("valid") = best_fit.valid
  );
}

//' Oblique Bentler factor rotation
//'
//' Rotate a loading matrix obliquely under Bentler's invariant pattern simplicity criterion
//' using a gradient-projection optimizer along the oblique (column-normalized) manifold.
//'
//' The criterion value `f` and its gradient `dQ/dL` at the rotated loadings
//' `L = A %*% solve(t(T))` define the search; the engine maps the gradient to the
//' transformation `T` on the manifold `diag(t(T) %*% T) = 1`, projects it onto the tangent
//' space, performs a sufficient-decrease line search, and retracts back onto the manifold by
//' column normalization. The Bentler criterion measures the departure of the cross-products of
//' squared loadings from a diagonal pattern; it is prone to local minima, so additional random
//' starts are recommended.
//'
//' Additional random starts may be requested. To bound runtime the solver screens each random
//' start by its objective, runs a short triage optimization on the best-screened starts, and
//' fully optimizes only those that improve on the current incumbent by at least
//' `triage_improve_tol`.
//'
//' @param L Numeric matrix. The unrotated loading matrix (variables by factors).
//' @param eps Numeric scalar. Convergence tolerance for the projected-gradient norm.
//' @param normalize Logical scalar. If `TRUE`, apply Kaiser normalization before rotation and
//'   reverse it afterwards.
//' @param random_starts Integer scalar. Number of additional random starts.
//' @param maxit Integer scalar. Maximum number of projected-gradient updates.
//' @param max_line_search Integer scalar. Maximum number of step-halving attempts after the
//'   initial trial step in each line-search phase.
//' @param step0 Numeric scalar. Initial step size used in the projected-gradient update.
//' @param screen_keep Integer scalar. Number of screened random starts retained for triage
//'   optimization.
//' @param triage_maxit Integer scalar. Number of short optimization iterations used in the
//'   triage stage.
//' @param triage_improve_tol Numeric scalar. Relative improvement required for a triaged start
//'   to be promoted to full optimization.
//'
//' @returns A named list with the rotated loadings, the transformation matrix `Th`
//'   (with `L %*% t(solve(Th))` reproducing the rotated loadings), the factor correlation
//'   matrix `Phi` (`t(Th) %*% Th`), the attained criterion value, and the convergence and
//'   validity flags.
//'
//' @references
//' Bentler, P. M. (1977). Factor simplicity index and transformations. *Psychometrika*, 42,
//' 277-295.
//'
//' Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms and
//' software for arbitrary rotation criteria in factor analysis. *Educational and
//' Psychological Measurement*, 65, 676-696.
//'
// [[Rcpp::export(.rotate_bentler_oblq)]]
Rcpp::List rotate_bentler_oblq(const arma::mat& L,
                               double eps = 1e-5,
                               bool normalize = true,
                               int random_starts = 0,
                               int maxit = 1000,
                               int max_line_search = 10,
                               double step0 = 1.0,
                               int screen_keep = 2,
                               int triage_maxit = 25,
                               double triage_improve_tol = 0.0) {
  if (L.n_rows == 0 || L.n_cols == 0) {
    Rcpp::stop("L must be a non-empty matrix.");
  }
  if (L.n_cols < 2) {
    Rcpp::stop("Oblique rotation requires at least two factors; use the R wrapper "
               "for single-factor solutions.");
  }
  if (!all_finite_cpp(L)) {
    Rcpp::stop("L must contain only finite values.");
  }
  validate_gpf_scalars(eps, maxit, max_line_search, step0, random_starts,
                       screen_keep, triage_maxit, triage_improve_tol);

  const arma::uword k = L.n_cols;

  // Kaiser normalization scales rows before rotation; because the weights act on rows and the
  // rotation (right-multiplying solve(t(T))) acts on columns, they cancel in the final
  // loadings, so the rotated loadings reproduce L %*% t(solve(Th)) (the documented
  // reproduction identity).
  arma::mat A_work = L;
  arma::vec W;
  if (normalize) {
    W = kaiser_weights(L);
    A_work.each_col() /= W;
  }

  BentlerCriterion crit;
  OblqCriterionManifold manifold{A_work, crit};
  arma::mat T_primary(k, k, fill::eye);

  GpfSummary summary = run_gpf_multistart(
    manifold, T_primary, eps, maxit, max_line_search, step0,
    random_starts, screen_keep, triage_maxit, triage_improve_tol
  );
  const GpfFit& best_fit = summary.best_fit;

  // Reconstruct the rotated loadings (L = A %*% solve(t(T)), Kaiser weights restored) and the
  // factor correlations (Phi = t(T) %*% T) from the winning transformation, shared with the
  // Procrustes oblique entry point.
  arma::mat loadings;
  arma::mat Phi;
  finalize_oblique(best_fit, A_work, W, normalize, loadings, Phi);

  return Rcpp::List::create(
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("Th") = best_fit.T,
    Rcpp::Named("Phi") = Phi,
    Rcpp::Named("value") = best_fit.value,
    Rcpp::Named("convergence") = best_fit.convergence,
    Rcpp::Named("valid") = best_fit.valid
  );
}
