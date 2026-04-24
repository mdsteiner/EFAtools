#include <RcppArmadillo.h>
#include <algorithm>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

using namespace Rcpp;
using namespace arma;

static arma::mat normalize_cols_cpp(const arma::mat& X, double eps = 1e-15) {
  arma::rowvec norms = sqrt(sum(square(X), 0));
  norms.transform([&](double val) { return (val < eps) ? eps : val; });
  return X * diagmat(1.0 / norms.t());
}

static inline arma::mat safe_inv(const arma::mat& X) {
  arma::mat Y;
  bool ok = arma::solve(Y, X, arma::eye<arma::mat>(X.n_rows, X.n_cols), arma::solve_opts::fast);
  if (!ok) {
    Y = arma::pinv(X);
  }
  return Y;
}

struct State {
  double f;
  double s;
  arma::mat U;
  arma::mat G;
  arma::mat Gp;
};

struct FitResult {
  arma::mat T;
  arma::mat U;
  arma::mat Table;
  double value;
  bool convergence;
  int iterations;
  double kappa_T;
};

struct Candidate {
  int start_index;
  arma::mat Tstart;
  double screen_value;
};

static State compute_state_kxk(const arma::mat& Tmat,
                               const arma::mat& S,
                               const arma::mat& C,
                               double BtB) {
  arma::mat invT = safe_inv(Tmat);
  arma::mat U = invT.t();
  arma::mat SU = S * U;
  arma::mat GU = SU - C;

  double f = 0.5 * (accu(U % SU) - 2.0 * accu(U % C) + BtB);
  arma::mat G = - U * GU.t() * U;

  arma::rowvec proj = sum(Tmat % G, 0);
  arma::mat Gp = G - Tmat * diagmat(proj.t());
  double s = norm(vectorise(Gp), 2);

  State out;
  out.f = f;
  out.s = s;
  out.U = U;
  out.G = G;
  out.Gp = Gp;
  return out;
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

static FitResult run_single_oblique_fit(const arma::mat& S,
                                        const arma::mat& C,
                                        double BtB,
                                        arma::mat Tmat,
                                        double eps,
                                        int maxit,
                                        int max_line_search,
                                        double step0) {
  Tmat = normalize_cols_cpp(Tmat);

  arma::mat table(maxit + 1, 4, fill::none);
  int filled = 0;
  bool convergence = false;
  double al = step0;

  State state = compute_state_kxk(Tmat, S, C, BtB);

  for (int iter = 0; iter <= maxit; ++iter) {
    table(filled, 0) = static_cast<double>(iter);
    table(filled, 1) = state.f;
    table(filled, 2) = std::log10(std::max(state.s, std::numeric_limits<double>::epsilon()));
    table(filled, 3) = al;
    ++filled;

    if (state.s < eps) {
      convergence = true;
      break;
    }

    al = 2.0 * al;
    bool improved = false;
    arma::mat T_candidate;
    State cand;

    for (int i = 0; i <= max_line_search; ++i) {
      arma::mat X = Tmat - al * state.Gp;
      T_candidate = normalize_cols_cpp(X);
      cand = compute_state_kxk(T_candidate, S, C, BtB);

      if ((state.f - cand.f) > 0.5 * state.s * state.s * al) {
        Tmat = T_candidate;
        state = cand;
        improved = true;
        break;
      }
      al = al / 2.0;
    }

    if (!improved) {
      Tmat = T_candidate;
      state = cand;
    }
  }

  FitResult out;
  out.T = Tmat;
  out.U = state.U;
  out.Table = table.rows(0, filled - 1);
  out.value = state.f;
  out.convergence = convergence;
  out.iterations = filled;
  out.kappa_T = arma::cond(Tmat);
  return out;
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
//' Additional random starts may be requested. To reduce runtime, the solver uses
//' a two-stage strategy for extra starts: cheap objective screening, followed by
//' short triage optimization, followed by full optimization only for starts that
//' improve on the current incumbent by at least `triage_improve_tol`.
//'
//' @param A Numeric matrix. Loading matrix to be rotated.
//' @param B Numeric matrix. Target loading matrix with the same dimensions as
//'   `A`.
//' @param S_r Optional numeric `k x k` matrix containing `crossprod(A)`.
//'   Supplying this is useful when the same `A` is rotated repeatedly.
//' @param T_init_r Optional numeric `k x k` starting transformation matrix.
//'   If `NULL`, the identity matrix is used for the primary start.
//' @param eps Numeric scalar. Convergence tolerance for the projected gradient
//'   norm.
//' @param maxit Integer scalar. Maximum number of full projected-gradient
//'   iterations.
//' @param max_line_search Integer scalar. Maximum number of step-halving
//'   attempts in each line-search phase.
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
//'   and multi-start summaries.
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
  if (A.n_cols <= 1) {
    Rcpp::stop("Oblique rotation does not make sense for single-factor models.");
  }

  const unsigned int k = A.n_cols;

  arma::mat A_work = A;
  arma::mat B_work = B;
  arma::vec W;

  if (normalize) {
    W = sqrt(sum(square(A_work), 1));
    for (unsigned int i = 0; i < W.n_elem; ++i) {
      if (W(i) < 1e-15) {
        W(i) = 1.0;
      }
    }
    A_work.each_col() /= W;
    B_work.each_col() /= W;
  }

  arma::mat S;
  if (normalize || S_r.isNull()) {
    S = A_work.t() * A_work;
  } else {
    S = Rcpp::as<arma::mat>(S_r);
    if (S.n_rows != k || S.n_cols != k) {
      Rcpp::stop("S must be a k x k matrix.");
    }
  }

  arma::mat C = A_work.t() * B_work;
  double BtB = accu(B_work % B_work);

  arma::mat T_primary;
  if (T_init_r.isNull()) {
    T_primary.eye(k, k);
  } else {
    T_primary = Rcpp::as<arma::mat>(T_init_r);
    if (T_primary.n_rows != k || T_primary.n_cols != k) {
      Rcpp::stop("T_init must be a k x k matrix.");
    }
    T_primary = normalize_cols_cpp(T_primary);
  }

  FitResult best_fit = run_single_oblique_fit(S, C, BtB, T_primary, eps, maxit, max_line_search, step0);
  double incumbent_value = best_fit.value;

  std::vector<double> all_values;
  std::vector<int> all_converged;
  std::vector<int> all_iterations;
  all_values.reserve(static_cast<std::size_t>(1 + std::max(0, random_starts)));
  all_converged.reserve(static_cast<std::size_t>(1 + std::max(0, random_starts)));
  all_iterations.reserve(static_cast<std::size_t>(1 + std::max(0, random_starts)));

  all_values.push_back(best_fit.value);
  all_converged.push_back(best_fit.convergence ? 1 : 0);
  all_iterations.push_back(best_fit.iterations);

  int best_start_index = 1;

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

    int keep = std::min(screen_keep, random_starts);
    keep = std::max(keep, 0);

    for (int i = 0; i < keep; ++i) {
      const Candidate& cand = candidates[static_cast<std::size_t>(i)];

      FitResult triage_fit = run_single_oblique_fit(
        S, C, BtB, cand.Tstart, eps, triage_maxit, max_line_search, step0
      );

      all_values.push_back(triage_fit.value);
      all_converged.push_back(triage_fit.convergence ? 1 : 0);
      all_iterations.push_back(triage_fit.iterations);

      double rel_improve = (incumbent_value - triage_fit.value) /
        (std::abs(incumbent_value) + 1e-12);

      if (rel_improve >= triage_improve_tol) {
        FitResult full_fit = run_single_oblique_fit(
          S, C, BtB, triage_fit.T, eps, maxit, max_line_search, step0
        );

        all_values.back() = full_fit.value;
        all_converged.back() = full_fit.convergence ? 1 : 0;
        all_iterations.back() = full_fit.iterations;

        bool better = false;
        if (full_fit.convergence && !best_fit.convergence) {
          better = true;
        } else if ((full_fit.convergence == best_fit.convergence) &&
                   (full_fit.value < best_fit.value)) {
          better = true;
        }

        if (better) {
          best_fit = full_fit;
          incumbent_value = full_fit.value;
          best_start_index = cand.start_index;
        }
      }
    }
  }

  arma::mat Lrot = A_work * best_fit.U;
  if (normalize) {
    Lrot.each_col() %= W;
  }

  arma::mat Phi = best_fit.T.t() * best_fit.T;

  Rcpp::NumericMatrix table_r = Rcpp::wrap(best_fit.Table);
  colnames(table_r) = Rcpp::CharacterVector::create("iter", "f", "log10_s", "step");

  return Rcpp::List::create(
    Rcpp::Named("loadings") = Lrot,
    Rcpp::Named("T") = best_fit.T,
    Rcpp::Named("Phi") = Phi,
    Rcpp::Named("value") = best_fit.value,
    Rcpp::Named("convergence") = best_fit.convergence,
    Rcpp::Named("iterations") = best_fit.iterations,
    Rcpp::Named("kappa_T") = best_fit.kappa_T,
    Rcpp::Named("Table") = table_r,
    Rcpp::Named("best_start_index") = best_start_index,
    Rcpp::Named("all_values") = Rcpp::wrap(all_values),
    Rcpp::Named("all_converged") = Rcpp::wrap(all_converged),
    Rcpp::Named("all_iterations") = Rcpp::wrap(all_iterations)
  );
}
