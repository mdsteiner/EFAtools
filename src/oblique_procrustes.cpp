#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

using namespace Rcpp;
using namespace arma;

static bool all_finite_cpp(const arma::mat& X) {
  for (arma::uword i = 0; i < X.n_elem; ++i) {
    if (!std::isfinite(X[i])) {
      return false;
    }
  }
  return true;
}

static bool is_valid_scalar_cpp(const double x) {
  return std::isfinite(x);
}

static arma::mat normalize_cols_cpp(const arma::mat& X, double eps = 1e-15) {
  arma::mat Y = X;

  for (arma::uword j = 0; j < Y.n_cols; ++j) {
    double nrm = arma::norm(Y.col(j), 2);

    if (!std::isfinite(nrm) || nrm < eps) {
      // A zero column violates diag(t(T) %*% T) = 1. This fallback keeps the
      // candidate on the constraint set; singular candidates are subsequently
      // rejected by the inverse check in compute_state_kxk().
      Y.col(j).zeros();
      Y(std::min(j, Y.n_rows - 1), j) = 1.0;
    } else {
      Y.col(j) /= nrm;
    }
  }

  return Y;
}

static bool inverse_checked_cpp(const arma::mat& X, arma::mat& invX) {
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

struct State {
  double f;
  double s;
  arma::mat U;
  arma::mat G;
  arma::mat Gp;
  bool valid;
};

struct FitResult {
  arma::mat T;
  arma::mat U;
  arma::mat Table;
  double value;
  bool convergence;
  bool valid;
  bool line_search_failed;
  int iterations;
  double kappa_T;
};

struct Candidate {
  int start_index;
  arma::mat Tstart;
  double screen_value;
};

// Outcome of a full multi-start oblique fit: the winning fit plus the per-start
// diagnostics surfaced by the single-matrix entry point. Shared by the single
// and batched entry points so the multi-start orchestration is written once.
struct ObliqueFitSummary {
  FitResult best_fit;
  int best_start_index;
  int n_triaged;
  int n_fully_optimized;
  std::vector<double> all_values;
  std::vector<int> all_converged;
  std::vector<int> all_iterations;
  std::vector<int> all_start_indices;
  std::vector<double> screen_values;
  std::vector<int> screen_start_indices;
};

static State invalid_state_cpp(const unsigned int k) {
  State out;
  out.f = std::numeric_limits<double>::infinity();
  out.s = std::numeric_limits<double>::infinity();
  out.U.zeros(k, k);
  out.G.zeros(k, k);
  out.Gp.zeros(k, k);
  out.valid = false;
  return out;
}

static State compute_state_kxk(const arma::mat& Tmat,
                               const arma::mat& S,
                               const arma::mat& C,
                               double BtB) {
  const unsigned int k = Tmat.n_cols;
  arma::mat invT;
  if (!inverse_checked_cpp(Tmat, invT)) {
    return invalid_state_cpp(k);
  }

  arma::mat U = invT.t();
  arma::mat SU = S * U;
  arma::mat GU = SU - C;

  double f = 0.5 * (accu(U % SU) - 2.0 * accu(U % C) + BtB);
  arma::mat G = -U * GU.t() * U;

  arma::rowvec proj = sum(Tmat % G, 0);
  arma::mat Gp = G - Tmat * diagmat(proj.t());
  double s = norm(vectorise(Gp), 2);

  State out;
  out.f = f;
  out.s = s;
  out.U = U;
  out.G = G;
  out.Gp = Gp;
  out.valid = std::isfinite(f) && std::isfinite(s) && all_finite_cpp(U) &&
    all_finite_cpp(G) && all_finite_cpp(Gp);

  if (!out.valid) {
    return invalid_state_cpp(k);
  }

  return out;
}

static double cond_checked_cpp(const arma::mat& X) {
  if (!all_finite_cpp(X)) {
    return NA_REAL;
  }
  try {
    double out = arma::cond(X);
    return std::isfinite(out) ? out : NA_REAL;
  } catch (...) {
    return NA_REAL;
  }
}

static arma::mat random_orthogonal_start_cpp(const unsigned int k) {
  arma::mat Z(k, k);
  for (unsigned int i = 0; i < k; ++i) {
    for (unsigned int j = 0; j < k; ++j) {
      Z(i, j) = R::rnorm(0.0, 1.0);
    }
  }

  arma::mat Q, Rm;
  bool ok = arma::qr_econ(Q, Rm, Z);
  if (!ok || Q.n_rows != k || Q.n_cols != k) {
    arma::mat Eye(k, k, fill::eye);
    return Eye;
  }

  arma::vec d = Rm.diag();
  for (unsigned int j = 0; j < d.n_elem; ++j) {
    double sgn = (d(j) >= 0.0) ? 1.0 : -1.0;
    Q.col(j) *= sgn;
  }

  return Q;
}

// Closed-form orthogonal Procrustes transform T = U V' from the SVD of A'B, the
// orthogonal matrix minimizing ||A T - B||_F. Mirrors the R-level warm start used
// by PROCRUSTES() (.procrustes_orthogonal_T); the product U V' is invariant to the
// LAPACK sign/ordering of the singular vectors, so it matches the R warm start.
// Returns false (leaving `T` untouched) if the decomposition fails, so callers can
// treat that slice as unalignable rather than propagating a hard error.
static bool orthogonal_procrustes_T_cpp(const arma::mat& A, const arma::mat& B,
                                        arma::mat& T) {
  arma::mat U, V;
  arma::vec s;
  if (!arma::svd(U, s, V, A.t() * B)) {
    return false;
  }
  T = U * V.t();
  return all_finite_cpp(T);
}

static FitResult run_single_oblique_fit(const arma::mat& S,
                                        const arma::mat& C,
                                        double BtB,
                                        arma::mat Tmat,
                                        double eps,
                                        int maxit,
                                        int max_line_search,
                                        double step0) {
  const unsigned int k = Tmat.n_cols;
  Tmat = normalize_cols_cpp(Tmat);

  arma::mat table(maxit + 1, 4);
  table.fill(NA_REAL);
  int filled = 0;
  bool convergence = false;
  bool line_search_failed = false;
  double al = step0;

  State state = compute_state_kxk(Tmat, S, C, BtB);

  for (int iter = 0; iter <= maxit; ++iter) {
    table(filled, 0) = static_cast<double>(iter);
    table(filled, 1) = state.f;
    table(filled, 2) = std::log10(std::max(state.s, std::numeric_limits<double>::epsilon()));
    table(filled, 3) = al;
    ++filled;

    if (state.valid && state.s < eps) {
      convergence = true;
      break;
    }

    if (!state.valid) {
      line_search_failed = true;
      break;
    }

    if (iter == maxit) {
      break;
    }

    al = 2.0 * al;
    if (!std::isfinite(al) || al <= 0.0) {
      al = step0;
    }

    bool armijo_improved = false;
    bool any_decrease = false;
    double best_decrease_value = state.f;
    double best_decrease_step = al;
    arma::mat best_decrease_T = Tmat;
    State best_decrease_state = state;

    double trial_step = al;
    for (int i = 0; i <= max_line_search; ++i) {
      arma::mat X = Tmat - trial_step * state.Gp;
      arma::mat T_candidate = normalize_cols_cpp(X);
      State cand = compute_state_kxk(T_candidate, S, C, BtB);

      if (cand.valid && cand.f < best_decrease_value) {
        any_decrease = true;
        best_decrease_value = cand.f;
        best_decrease_step = trial_step;
        best_decrease_T = T_candidate;
        best_decrease_state = cand;
      }

      if (cand.valid && (state.f - cand.f) > 0.5 * state.s * state.s * trial_step) {
        Tmat = T_candidate;
        state = cand;
        al = trial_step;
        armijo_improved = true;
        break;
      }

      trial_step = trial_step / 2.0;
    }

    if (!armijo_improved) {
      if (any_decrease) {
        // Preserve monotone descent even when the sufficient-decrease test is
        // too strict because of numerical noise near convergence.
        Tmat = best_decrease_T;
        state = best_decrease_state;
        al = best_decrease_step;
      } else {
        line_search_failed = true;
        break;
      }
    }
  }

  FitResult out;
  out.T = Tmat;
  out.U = state.U;
  out.Table = table.rows(0, filled - 1);
  out.value = state.f;
  out.convergence = convergence;
  out.valid = state.valid;
  out.line_search_failed = line_search_failed;
  out.iterations = filled - 1;
  out.kappa_T = cond_checked_cpp(Tmat);
  return out;
}

static bool fit_is_better_cpp(const FitResult& candidate,
                              const FitResult& incumbent) {
  // Selection across multi-start fits is driven by the rotation criterion value:
  // the attained objective measures closeness to the target, whereas convergence
  // is only a numerical termination flag (the projected-gradient norm fell below
  // eps within maxit). Each start's reported value is the objective of an actually
  // attained transformation, so a lower value is a strictly better solution and is
  // never discarded merely because that start has not formally converged. Exact
  // ties are broken toward the converged fit. An invalid (non-invertible)
  // candidate has an infinite objective and can never win.
  if (!candidate.valid) {
    return false;
  }
  if (!incumbent.valid) {
    return true;
  }
  if (candidate.value < incumbent.value) {
    return true;
  }
  if (candidate.value == incumbent.value) {
    return candidate.convergence && !incumbent.convergence;
  }
  return false;
}

// Run the primary fit from `T_primary`, then optionally screen, triage, and fully
// optimize the requested random starts, keeping the best fit by objective value.
// This is the multi-start orchestration shared by the single-matrix and batched
// oblique entry points. Random starts are drawn serially with R::rnorm.
static ObliqueFitSummary run_oblique_multistart(const arma::mat& S,
                                                const arma::mat& C,
                                                double BtB,
                                                const arma::mat& T_primary,
                                                double eps,
                                                int maxit,
                                                int max_line_search,
                                                double step0,
                                                int random_starts,
                                                int screen_keep,
                                                int triage_maxit,
                                                double triage_improve_tol) {
  const unsigned int k = T_primary.n_cols;

  ObliqueFitSummary summary;
  FitResult best_fit = run_single_oblique_fit(S, C, BtB, T_primary, eps, maxit, max_line_search, step0);
  double incumbent_value = best_fit.value;

  summary.all_values.reserve(static_cast<std::size_t>(1 + std::max(0, random_starts)));
  summary.all_converged.reserve(static_cast<std::size_t>(1 + std::max(0, random_starts)));
  summary.all_iterations.reserve(static_cast<std::size_t>(1 + std::max(0, random_starts)));
  summary.all_start_indices.reserve(static_cast<std::size_t>(1 + std::max(0, random_starts)));
  summary.screen_values.reserve(static_cast<std::size_t>(std::max(0, random_starts)));
  summary.screen_start_indices.reserve(static_cast<std::size_t>(std::max(0, random_starts)));

  summary.all_values.push_back(best_fit.value);
  summary.all_converged.push_back(best_fit.convergence ? 1 : 0);
  summary.all_iterations.push_back(best_fit.iterations);
  summary.all_start_indices.push_back(1);

  int best_start_index = 1;
  int n_triaged = 0;
  int n_fully_optimized = 1;

  if (random_starts > 0) {
    std::vector<Candidate> candidates;
    candidates.reserve(static_cast<std::size_t>(random_starts));

    for (int s = 0; s < random_starts; ++s) {
      arma::mat Tstart = random_orthogonal_start_cpp(k);
      State screen_state = compute_state_kxk(Tstart, S, C, BtB);
      Candidate cand;
      cand.start_index = s + 2;
      cand.Tstart = Tstart;
      cand.screen_value = screen_state.f;
      candidates.push_back(cand);
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& a, const Candidate& b) {
                return a.screen_value < b.screen_value;
              });

    for (const Candidate& cand : candidates) {
      summary.screen_values.push_back(cand.screen_value);
      summary.screen_start_indices.push_back(cand.start_index);
    }

    int keep = std::min(screen_keep, random_starts);
    keep = std::max(keep, 0);

    for (int i = 0; i < keep; ++i) {
      const Candidate& cand = candidates[static_cast<std::size_t>(i)];

      FitResult triage_fit = run_single_oblique_fit(
        S, C, BtB, cand.Tstart, eps, triage_maxit, max_line_search, step0
      );
      ++n_triaged;

      summary.all_values.push_back(triage_fit.value);
      summary.all_converged.push_back(triage_fit.convergence ? 1 : 0);
      summary.all_iterations.push_back(triage_fit.iterations);
      summary.all_start_indices.push_back(cand.start_index);

      double rel_improve;
      if (std::isfinite(incumbent_value)) {
        rel_improve = (incumbent_value - triage_fit.value) /
          (std::abs(incumbent_value) + 1e-12);
      } else {
        rel_improve = triage_fit.valid ? std::numeric_limits<double>::infinity() :
          -std::numeric_limits<double>::infinity();
      }

      if (triage_fit.valid && rel_improve >= triage_improve_tol) {
        FitResult full_fit = run_single_oblique_fit(
          S, C, BtB, triage_fit.T, eps, maxit, max_line_search, step0
        );
        ++n_fully_optimized;

        summary.all_values.back() = full_fit.value;
        summary.all_converged.back() = full_fit.convergence ? 1 : 0;
        summary.all_iterations.back() = full_fit.iterations;

        if (fit_is_better_cpp(full_fit, best_fit)) {
          best_fit = full_fit;
          incumbent_value = full_fit.value;
          best_start_index = cand.start_index;
        }
      }
    }
  }

  summary.best_fit = best_fit;
  summary.best_start_index = best_start_index;
  summary.n_triaged = n_triaged;
  summary.n_fully_optimized = n_fully_optimized;
  return summary;
}

// Validate the shared oblique-solver control scalars. Used by both the single and
// batched entry points so the contract and its messages live in one place.
static void validate_oblique_scalars(double eps, int maxit, int max_line_search,
                                     double step0, int random_starts, int screen_keep,
                                     int triage_maxit, double triage_improve_tol) {
  if (!is_valid_scalar_cpp(eps) || eps <= 0.0) {
    Rcpp::stop("eps must be a positive finite scalar.");
  }
  if (maxit < 0) {
    Rcpp::stop("maxit must be non-negative.");
  }
  if (max_line_search < 0) {
    Rcpp::stop("max_line_search must be non-negative.");
  }
  if (!is_valid_scalar_cpp(step0) || step0 <= 0.0) {
    Rcpp::stop("step0 must be a positive finite scalar.");
  }
  if (random_starts < 0) {
    Rcpp::stop("random_starts must be non-negative.");
  }
  if (screen_keep < 0) {
    Rcpp::stop("screen_keep must be non-negative.");
  }
  if (triage_maxit < 0) {
    Rcpp::stop("triage_maxit must be non-negative.");
  }
  if (!is_valid_scalar_cpp(triage_improve_tol) || triage_improve_tol < 0.0) {
    Rcpp::stop("triage_improve_tol must be a non-negative finite scalar.");
  }
}

// Kaiser row normalization shared by both entry points: scale each row of A_work
// and B_work by the inverse row norm of A_work, returning the weights so the
// rotated loadings can be back-transformed. A zero/non-finite row norm is floored
// to 1 to leave that row unchanged.
static arma::vec kaiser_normalize_rows(arma::mat& A_work, arma::mat& B_work) {
  arma::vec W = sqrt(sum(square(A_work), 1));
  for (arma::uword i = 0; i < W.n_elem; ++i) {
    if (!std::isfinite(W(i)) || W(i) < 1e-15) {
      W(i) = 1.0;
    }
  }
  A_work.each_col() /= W;
  B_work.each_col() /= W;
  return W;
}

// Assemble the rotated loadings and factor correlations from a fitted
// transformation: L = A_work U (un-normalized when Kaiser weights were applied),
// Phi = symmetrized T'T with unit diagonal. Shared by both entry points.
static void finalize_oblique(const FitResult& best_fit, const arma::mat& A_work,
                             const arma::vec& W, bool normalize,
                             arma::mat& Lrot, arma::mat& Phi) {
  Lrot = A_work * best_fit.U;
  if (normalize) {
    Lrot.each_col() %= W;
  }
  Phi = best_fit.T.t() * best_fit.T;
  Phi = (Phi + Phi.t()) / 2.0;
  Phi.diag().fill(1.0);
}

//' Oblique Procrustes target rotation using a k x k inner objective
//'
//' Compute an oblique target rotation for a loading matrix using a
//' `targetQ`-compatible parameterization and a `k x k` objective.
//'
//' The rotated loading matrix is defined as
//' `L = A %*% solve(t(T))`, and the corresponding factor correlation matrix is
//' `Phi = t(T) %*% T`. The optimization is carried out over the transformation
//' matrix `T` under the oblique normalization constraint `diag(t(T) %*% T) = 1`.
//'
//' The line search is monotone: a candidate is accepted only if it satisfies the
//' sufficient-decrease condition, or, as a numerical fallback, if it at least
//' decreases the objective after all step halvings are exhausted. Non-invertible
//' candidate transformations are rejected rather than evaluated through a
//' pseudo-inverse.
//'
//' Additional random starts may be requested. To reduce runtime, the solver uses
//' a two-stage strategy for extra starts: cheap objective screening, followed by
//' short triage optimization, followed by full optimization only for starts that
//' improve on the current incumbent by at least `triage_improve_tol`.
//'
//' @param A Numeric matrix. Loading matrix to be rotated.
//' @param B Numeric matrix. Target loading matrix with the same dimensions as
//'   `A`.
//' @param S_r Optional numeric `k x k` matrix containing `crossprod(A)`.
//'   Supplying this is useful when the same `A` is rotated repeatedly. Ignored
//'   when `normalize = TRUE` because the normalized cross-product is different.
//' @param T_init_r Optional numeric `k x k` starting transformation matrix.
//'   If `NULL`, the identity matrix is used for the primary start.
//' @param eps Numeric scalar. Convergence tolerance for the projected-gradient
//'   norm.
//' @param maxit Integer scalar. Maximum number of full projected-gradient
//'   updates.
//' @param max_line_search Integer scalar. Maximum number of step-halving
//'   attempts after the initial trial step in each line-search phase.
//' @param step0 Numeric scalar. Initial step size used in the projected-gradient
//'   update.
//' @param normalize Logical scalar. If `TRUE`, apply Kaiser normalization before
//'   rotation and reverse it after rotation.
//' @param random_starts Integer scalar. Number of additional random starts.
//' @param screen_keep Integer scalar. Number of screened random starts retained
//'   for triage optimization.
//' @param triage_maxit Integer scalar. Number of short optimization iterations
//'   used in the triage stage.
//' @param triage_improve_tol Numeric scalar. Relative improvement required for a
//'   triaged start to be promoted to full optimization.
//'
//' @returns A named list containing the rotated loadings, transformation matrix,
//'   factor correlation matrix, target criterion value, convergence diagnostics,
//'   line-search diagnostics, and multi-start summaries.
//'
//' @details
//' The routine is intended for repeated oblique target rotations in workflows
//' such as bootstrap alignment or consensus alignment of exploratory factor
//' solutions across multiply imputed datasets. It follows the same oblique
//' transformation convention as `GPArotation::targetQ()`.
//'
//' @references
//' Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms
//' and software for arbitrary rotation criteria in factor analysis.
//' *Educational and Psychological Measurement*, 65, 676-696.
//'
//' Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika*, 40,
//' 33-51.
//'
// [[Rcpp::export(.oblique_procrustes)]]
Rcpp::List oblique_procrustes(const arma::mat& A,
                              const arma::mat& B,
                              Rcpp::Nullable<Rcpp::NumericMatrix> S_r = R_NilValue,
                              Rcpp::Nullable<Rcpp::NumericMatrix> T_init_r = R_NilValue,
                              double eps = 1e-5,
                              int maxit = 1000,
                              int max_line_search = 10,
                              double step0 = 1.0,
                              bool normalize = false,
                              int random_starts = 0,
                              int screen_keep = 2,
                              int triage_maxit = 25,
                              double triage_improve_tol = 0.0) {
  if (A.n_rows != B.n_rows || A.n_cols != B.n_cols) {
    Rcpp::stop("A and B must have identical dimensions.");
  }
  if (A.n_rows == 0 || A.n_cols == 0) {
    Rcpp::stop("A and B must be non-empty matrices.");
  }
  if (A.n_cols <= 1) {
    Rcpp::stop("Use the R wrapper for one-factor Procrustes alignment; oblique optimization is not needed.");
  }
  if (!all_finite_cpp(A) || !all_finite_cpp(B)) {
    Rcpp::stop("A and B must contain only finite values.");
  }
  validate_oblique_scalars(eps, maxit, max_line_search, step0, random_starts,
                           screen_keep, triage_maxit, triage_improve_tol);

  const unsigned int k = A.n_cols;

  arma::mat A_work = A;
  arma::mat B_work = B;
  arma::vec W;

  if (normalize) {
    W = kaiser_normalize_rows(A_work, B_work);
  }

  arma::mat S;
  if (normalize || S_r.isNull()) {
    S = A_work.t() * A_work;
  } else {
    S = Rcpp::as<arma::mat>(S_r);
    if (S.n_rows != k || S.n_cols != k) {
      Rcpp::stop("S must be a k x k matrix.");
    }
    if (!all_finite_cpp(S)) {
      Rcpp::stop("S must contain only finite values.");
    }
  }

  arma::mat C = A_work.t() * B_work;
  double BtB = accu(B_work % B_work);
  if (!all_finite_cpp(S) || !all_finite_cpp(C) || !std::isfinite(BtB)) {
    Rcpp::stop("The Procrustes objective could not be formed from finite values.");
  }

  arma::mat T_primary(k, k, fill::eye);
  if (!T_init_r.isNull()) {
    T_primary = Rcpp::as<arma::mat>(T_init_r);
    if (T_primary.n_rows != k || T_primary.n_cols != k) {
      Rcpp::stop("T_init must be a k x k matrix.");
    }
    if (!all_finite_cpp(T_primary)) {
      Rcpp::stop("T_init must contain only finite values.");
    }
    T_primary = normalize_cols_cpp(T_primary);

    arma::mat tmp_inv;
    if (!inverse_checked_cpp(T_primary, tmp_inv)) {
      Rcpp::stop("T_init must be nonsingular after column normalization.");
    }
  }

  ObliqueFitSummary summary = run_oblique_multistart(
    S, C, BtB, T_primary, eps, maxit, max_line_search, step0,
    random_starts, screen_keep, triage_maxit, triage_improve_tol
  );
  const FitResult& best_fit = summary.best_fit;

  arma::mat Lrot, Phi;
  finalize_oblique(best_fit, A_work, W, normalize, Lrot, Phi);

  Rcpp::NumericMatrix table_r = Rcpp::wrap(best_fit.Table);
  colnames(table_r) = Rcpp::CharacterVector::create("iter", "f", "log10_s", "step");

  return Rcpp::List::create(
    Rcpp::Named("loadings") = Lrot,
    Rcpp::Named("T") = best_fit.T,
    Rcpp::Named("Phi") = Phi,
    Rcpp::Named("value") = best_fit.value,
    Rcpp::Named("convergence") = best_fit.convergence,
    Rcpp::Named("valid") = best_fit.valid,
    Rcpp::Named("iterations") = best_fit.iterations,
    Rcpp::Named("kappa_T") = best_fit.kappa_T,
    Rcpp::Named("Table") = table_r,
    Rcpp::Named("line_search_failed") = best_fit.line_search_failed,
    Rcpp::Named("method") = "oblique_procrustes",
    Rcpp::Named("best_start_index") = summary.best_start_index,
    Rcpp::Named("all_start_indices") = Rcpp::wrap(summary.all_start_indices),
    Rcpp::Named("all_values") = Rcpp::wrap(summary.all_values),
    Rcpp::Named("all_converged") = Rcpp::wrap(summary.all_converged),
    Rcpp::Named("all_iterations") = Rcpp::wrap(summary.all_iterations),
    Rcpp::Named("screen_start_indices") = Rcpp::wrap(summary.screen_start_indices),
    Rcpp::Named("screen_values") = Rcpp::wrap(summary.screen_values),
    Rcpp::Named("n_random_starts") = random_starts,
    Rcpp::Named("n_screened") = static_cast<int>(summary.screen_values.size()),
    Rcpp::Named("n_triaged") = summary.n_triaged,
    Rcpp::Named("n_fully_optimized") = summary.n_fully_optimized
  );
}

//' Batched oblique Procrustes target rotation over a cube of loading matrices
//'
//' Align each slice of a loading-matrix cube to a single shared target using the
//' same oblique target rotation as `.oblique_procrustes()`, in one call. This
//' removes the per-replicate marshalling overhead of looping `PROCRUSTES()` in R
//' over bootstrap or multiple-imputation arrays.
//'
//' Each slice `A[, , i]` is aligned to `B`. For a single-factor cube the alignment
//' reduces to the closed-form sign match `T = sign(crossprod(A_i, B))` with factor
//' correlation `1`, matching the one-factor short-circuit in `PROCRUSTES()`. For
//' two or more factors the slice is warm-started from the closed-form orthogonal
//' Procrustes solution (mirroring `PROCRUSTES()`) and optimized with the same
//' multi-start oblique solver as `.oblique_procrustes()`. Random starts are drawn
//' serially with `R::rnorm` in the calling process.
//'
//' Slices are aligned independently. A slice that cannot be aligned (a non-finite
//' loading matrix, a failed warm-start decomposition, an invalid fit, or any
//' linear-algebra exception) is reported with `valid = FALSE` and `NA` for the
//' loadings, factor correlations, and all other per-slice diagnostics, rather than
//' aborting the whole call, so one degenerate replicate does not discard the rest.
//'
//' @param A Numeric array of dimension `n x m x b`: the `b` loading matrices to
//'   align.
//' @param B Numeric `n x m` target loading matrix shared across all slices.
//' @param eps Numeric scalar. Convergence tolerance for the projected-gradient
//'   norm.
//' @param maxit Integer scalar. Maximum number of full projected-gradient updates.
//' @param max_line_search Integer scalar. Maximum number of step-halving attempts
//'   after the initial trial step in each line-search phase.
//' @param step0 Numeric scalar. Initial step size used in the projected-gradient
//'   update.
//' @param normalize Logical scalar. If `TRUE`, apply Kaiser normalization before
//'   rotation and reverse it afterwards (ignored for single-factor slices).
//' @param random_starts Integer scalar. Number of additional random starts per
//'   slice.
//' @param screen_keep Integer scalar. Number of screened random starts retained
//'   for triage optimization.
//' @param triage_maxit Integer scalar. Number of short optimization iterations
//'   used in the triage stage.
//' @param triage_improve_tol Numeric scalar. Relative improvement required for a
//'   triaged start to be promoted to full optimization.
//'
//' @returns A named list with the aligned-loadings array `loadings` (`n x m x b`),
//'   the factor-correlation array `Phi` (`m x m x b`), and the per-slice
//'   diagnostics `valid`, `convergence`, `value`, `iterations`, and
//'   `line_search_failed`.
//'
// [[Rcpp::export(.oblique_procrustes_batch)]]
Rcpp::List oblique_procrustes_batch(const arma::cube& A,
                                    const arma::mat& B,
                                    double eps = 1e-5,
                                    int maxit = 1000,
                                    int max_line_search = 10,
                                    double step0 = 1.0,
                                    bool normalize = false,
                                    int random_starts = 0,
                                    int screen_keep = 2,
                                    int triage_maxit = 25,
                                    double triage_improve_tol = 0.0) {
  const arma::uword n = A.n_rows;
  const arma::uword m = A.n_cols;
  const arma::uword b = A.n_slices;

  if (n == 0 || m == 0 || b == 0) {
    Rcpp::stop("A must be a non-empty n x m x b array.");
  }
  if (B.n_rows != n || B.n_cols != m) {
    Rcpp::stop("B must have the same dimensions as each slice of A.");
  }
  if (!all_finite_cpp(B)) {
    Rcpp::stop("B must contain only finite values.");
  }
  validate_oblique_scalars(eps, maxit, max_line_search, step0, random_starts,
                           screen_keep, triage_maxit, triage_improve_tol);

  // Each slice is aligned independently. A slice that cannot be aligned -- a
  // non-finite loading matrix, a failed warm-start decomposition, an invalid fit,
  // or any unexpected linear-algebra exception -- is marked invalid and left as
  // NA, rather than aborting the whole batch. This preserves the per-replicate
  // failure isolation of the per-replicate PROCRUSTES() loop it replaces: one
  // degenerate bootstrap replicate drops out instead of losing the entire run.
  // An invalid slice is reported as valid = FALSE with everything else NA, so a
  // skipped replicate is never mistaken for a fit that ran but failed to converge.
  arma::cube loadings(n, m, b);
  arma::cube Phi(m, m, b);
  loadings.fill(NA_REAL);
  Phi.fill(NA_REAL);
  Rcpp::LogicalVector valid(b);
  Rcpp::LogicalVector convergence(b, NA_LOGICAL);
  Rcpp::NumericVector value(b, NA_REAL);
  Rcpp::IntegerVector iterations(b, NA_INTEGER);
  Rcpp::LogicalVector line_search_failed(b, NA_LOGICAL);

  if (m == 1) {
    // Single factor: oblique alignment reduces to a sign match (Phi = 1), with no
    // optimization or random starts, mirroring the one-factor short-circuit in
    // PROCRUSTES(). Kaiser normalization does not affect a sign and is skipped.
    // The factor count is a whole-cube constant, so this branch is taken once.
    for (arma::uword i = 0; i < b; ++i) {
      Rcpp::checkUserInterrupt();
      const R_xlen_t out_i = static_cast<R_xlen_t>(i);
      const arma::mat A_i = A.slice(i);
      if (!A_i.is_finite()) {
        continue;
      }
      double t = (arma::accu(A_i % B) >= 0.0) ? 1.0 : -1.0;
      arma::mat L = A_i * t;
      loadings.slice(i) = L;
      Phi.slice(i).fill(1.0);
      valid[out_i] = true;
      convergence[out_i] = true;
      value[out_i] = 0.5 * arma::accu(arma::square(L - B));
      iterations[out_i] = 0;
      line_search_failed[out_i] = false;
    }
  } else {
    // The target B (and, when not normalizing, crossprod(B)) is shared across all
    // slices, so hoist the invariant rather than rebuilding it per replicate.
    const double BtB_shared = normalize ? 0.0 : accu(B % B);

    for (arma::uword i = 0; i < b; ++i) {
      Rcpp::checkUserInterrupt();
      const R_xlen_t out_i = static_cast<R_xlen_t>(i);
      const arma::mat A_i = A.slice(i);
      if (!A_i.is_finite()) {
        continue;
      }

      try {
        arma::mat A_work;  // only materialized (and modified) under normalize
        arma::vec W;
        arma::mat S, C, T_primary;
        double BtB;
        bool warm_ok;

        if (normalize) {
          A_work = A_i;
          arma::mat B_work = B;
          W = kaiser_normalize_rows(A_work, B_work);
          S = A_work.t() * A_work;
          C = A_work.t() * B_work;
          BtB = accu(B_work % B_work);
          warm_ok = orthogonal_procrustes_T_cpp(A_work, B_work, T_primary);
        } else {
          S = A_i.t() * A_i;
          C = A_i.t() * B;
          BtB = BtB_shared;
          warm_ok = orthogonal_procrustes_T_cpp(A_i, B, T_primary);
        }
        if (!warm_ok) {
          continue;
        }

        ObliqueFitSummary summary = run_oblique_multistart(
          S, C, BtB, T_primary, eps, maxit, max_line_search, step0,
          random_starts, screen_keep, triage_maxit, triage_improve_tol
        );
        const FitResult& best_fit = summary.best_fit;
        if (!best_fit.valid) {
          continue;
        }

        // Form the loadings from the (normalized) matrix actually fitted, without
        // copying A_i again when no normalization is applied.
        const arma::mat& A_fit = normalize ? A_work : A_i;
        arma::mat Lrot, Phi_i;
        finalize_oblique(best_fit, A_fit, W, normalize, Lrot, Phi_i);

        loadings.slice(i) = Lrot;
        Phi.slice(i) = Phi_i;
        valid[out_i] = true;
        convergence[out_i] = best_fit.convergence;
        value[out_i] = best_fit.value;
        iterations[out_i] = best_fit.iterations;
        line_search_failed[out_i] = best_fit.line_search_failed;
      } catch (const std::exception&) {
        // Leave this slice as NA / invalid and continue with the remaining ones.
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("Phi") = Phi,
    Rcpp::Named("valid") = valid,
    Rcpp::Named("convergence") = convergence,
    Rcpp::Named("value") = value,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("line_search_failed") = line_search_failed
  );
}
