#ifndef EFATOOLS_GPF_ENGINE_H
#define EFATOOLS_GPF_ENGINE_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "gpf_common.h"

// Gradient-projection (GP) manifold optimizer shared by the rotation solvers. The
// solver structure -- the projected-gradient iteration with a sufficient-decrease
// line search, and the screen -> triage -> optimize multi-start orchestration -- is
// identical across solvers; only three operations differ between them and are
// supplied by a `Manifold` strategy:
//
//   GpfState compute(const arma::mat& T) const;            // objective + projected gradient at T
//   bool     retract(const arma::mat& X, arma::mat& Tout); // map a raw step back onto the manifold
//   arma::mat normalize_start(const arma::mat& T) const;   // condition a starting transformation
//
// `compute` returns the objective `f`, the projected-gradient norm `s`, the projected
// gradient `Gp`, and a validity flag. `retract` returns false to reject a candidate
// (e.g. a failed decomposition), which the line search skips. This lets the oblique
// Procrustes solver (k x k inverse objective, column-normalization manifold) and the
// orthogonal rotation solver (criterion objective, Stiefel manifold) share one engine.
//
// Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms and
// software for arbitrary rotation criteria in factor analysis. Educational and
// Psychological Measurement, 65, 676-696.

// State of the projected-gradient iteration at a transformation T.
struct GpfState {
  double f;
  double s;
  arma::mat Gp;
  bool valid;
};

// A converged (or terminated) single-start fit.
struct GpfFit {
  arma::mat T;
  arma::mat Table;
  double value;
  bool convergence;
  bool valid;
  bool line_search_failed;
  int iterations;
};

struct GpfCandidate {
  int start_index;
  arma::mat Tstart;
  double screen_value;
};

// The winning fit plus the per-start multi-start diagnostics.
struct GpfSummary {
  GpfFit best_fit;
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

inline GpfState gpf_invalid_state(const arma::uword k) {
  GpfState out;
  out.f = std::numeric_limits<double>::infinity();
  out.s = std::numeric_limits<double>::infinity();
  out.Gp.zeros(k, k);
  out.valid = false;
  return out;
}

// Selection across multi-start fits is driven by the objective value: the attained
// objective measures solution quality, whereas convergence is only a termination flag
// (the projected-gradient norm fell below eps within maxit). A lower value always wins;
// an exact tie breaks toward the converged fit; an invalid fit (infinite objective)
// can never win.
inline bool gpf_fit_is_better(const GpfFit& candidate, const GpfFit& incumbent) {
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

// Validate the shared solver control scalars so the contract and its messages live in
// one place.
inline void validate_gpf_scalars(double eps, int maxit, int max_line_search,
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

// One projected-gradient optimization from a starting transformation. The line search
// is monotone: a candidate is accepted only if it satisfies the sufficient-decrease
// (Armijo) condition, or, as a numerical fallback, if it at least decreases the
// objective after all step halvings are exhausted. The step is doubled between
// iterations and halved within the line search.
template <typename Manifold>
GpfFit run_single_gpf_fit(const Manifold& manifold,
                          arma::mat Tmat,
                          double eps,
                          int maxit,
                          int max_line_search,
                          double step0) {
  Tmat = manifold.normalize_start(Tmat);

  arma::mat table(maxit + 1, 4);
  table.fill(NA_REAL);
  int filled = 0;
  bool convergence = false;
  bool line_search_failed = false;
  double al = step0;

  GpfState state = manifold.compute(Tmat);

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
    GpfState best_decrease_state = state;

    double trial_step = al;
    for (int i = 0; i <= max_line_search; ++i) {
      arma::mat X = Tmat - trial_step * state.Gp;
      arma::mat T_candidate;

      if (manifold.retract(X, T_candidate)) {
        GpfState cand = manifold.compute(T_candidate);

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
      }

      trial_step = trial_step / 2.0;
    }

    if (!armijo_improved) {
      if (any_decrease) {
        // Preserve monotone descent even when the sufficient-decrease test is too
        // strict because of numerical noise near convergence.
        Tmat = best_decrease_T;
        state = best_decrease_state;
        al = best_decrease_step;
      } else {
        line_search_failed = true;
        break;
      }
    }
  }

  GpfFit out;
  out.T = Tmat;
  out.Table = table.rows(0, filled - 1);
  out.value = state.f;
  out.convergence = convergence;
  out.valid = state.valid;
  out.line_search_failed = line_search_failed;
  out.iterations = filled - 1;
  return out;
}

// Run the primary fit from T_primary, then optionally screen, triage, and fully
// optimize the requested random starts, keeping the best fit by objective value. To
// bound runtime, extra starts use a two-stage strategy: cheap objective screening,
// short triage optimization, then full optimization only for starts that improve on
// the incumbent by at least triage_improve_tol. Random orthogonal starts are drawn
// serially with R::rnorm.
template <typename Manifold>
GpfSummary run_gpf_multistart(const Manifold& manifold,
                              const arma::mat& T_primary,
                              double eps,
                              int maxit,
                              int max_line_search,
                              double step0,
                              int random_starts,
                              int screen_keep,
                              int triage_maxit,
                              double triage_improve_tol) {
  const arma::uword k = T_primary.n_cols;

  GpfSummary summary;
  GpfFit best_fit = run_single_gpf_fit(manifold, T_primary, eps, maxit, max_line_search, step0);
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
    std::vector<GpfCandidate> candidates;
    candidates.reserve(static_cast<std::size_t>(random_starts));

    for (int s = 0; s < random_starts; ++s) {
      arma::mat Tstart = random_orthogonal_start_cpp(k);
      GpfState screen_state = manifold.compute(Tstart);
      GpfCandidate cand;
      cand.start_index = s + 2;
      cand.Tstart = Tstart;
      cand.screen_value = screen_state.f;
      candidates.push_back(cand);
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const GpfCandidate& a, const GpfCandidate& b) {
                return a.screen_value < b.screen_value;
              });

    for (const GpfCandidate& cand : candidates) {
      summary.screen_values.push_back(cand.screen_value);
      summary.screen_start_indices.push_back(cand.start_index);
    }

    int keep = std::min(screen_keep, random_starts);
    keep = std::max(keep, 0);

    for (int i = 0; i < keep; ++i) {
      const GpfCandidate& cand = candidates[static_cast<std::size_t>(i)];

      GpfFit triage_fit = run_single_gpf_fit(
        manifold, cand.Tstart, eps, triage_maxit, max_line_search, step0
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
        GpfFit full_fit = run_single_gpf_fit(
          manifold, triage_fit.T, eps, maxit, max_line_search, step0
        );
        ++n_fully_optimized;

        summary.all_values.back() = full_fit.value;
        summary.all_converged.back() = full_fit.convergence ? 1 : 0;
        summary.all_iterations.back() = full_fit.iterations;

        if (gpf_fit_is_better(full_fit, best_fit)) {
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

#endif  // EFATOOLS_GPF_ENGINE_H
