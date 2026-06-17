#ifndef EFATOOLS_GPF_ENGINE_H
#define EFATOOLS_GPF_ENGINE_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
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

// An invalid state has an infinite objective so it can never win selection or pass the
// convergence test. Its projected gradient is never read (the line search only steps
// from valid states), so Gp is left empty rather than allocated.
inline GpfState gpf_invalid_state() {
  GpfState out;
  out.f = std::numeric_limits<double>::infinity();
  out.s = std::numeric_limits<double>::infinity();
  out.valid = false;
  return out;
}

// Assemble the rotated loadings and factor correlations from a fitted oblique
// transformation, shared by the oblique entry points (the Procrustes target rotation and
// the criterion rotations). The loadings transform U = solve(t(T)) is reconstructed from
// the winning T here rather than carried through the iteration state: L = A_work U
// (un-normalized by the Kaiser weights W when they were applied), Phi = symmetrized t(T) T
// with unit diagonal. An invalid winner yields a zero U (a valid fit always has an
// invertible T, so this only zeros genuinely degenerate results).
inline void finalize_oblique(const GpfFit& best_fit, const arma::mat& A_work,
                             const arma::vec& W, bool normalize,
                             arma::mat& Lrot, arma::mat& Phi) {
  arma::mat invT;
  arma::mat U;
  if (best_fit.valid && inverse_checked_cpp(best_fit.T, invT)) {
    U = invT.t();
  } else {
    U.zeros(best_fit.T.n_rows, best_fit.T.n_cols);
  }
  Lrot = A_work * U;
  if (normalize) {
    Lrot.each_col() %= W;
  }
  Phi = best_fit.T.t() * best_fit.T;
  Phi = (Phi + Phi.t()) / 2.0;
  Phi.diag().fill(1.0);
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

// Validate the loading matrix shared by every criterion rotation entry point, keeping the
// input contract and its messages with the scalar-control contract above. `family` is the
// rotation family name ("Orthogonal"/"Oblique") used in the at-least-two-factors message; the
// single-factor case is handled by the R wrapper, so the compiled entries require k >= 2.
inline void validate_gpf_input(const arma::mat& L, const char* family) {
  if (L.n_rows == 0 || L.n_cols == 0) {
    Rcpp::stop("L must be a non-empty matrix.");
  }
  if (L.n_cols < 2) {
    Rcpp::stop(std::string(family) + " rotation requires at least two factors; use the "
               "R wrapper for single-factor solutions.");
  }
  if (!all_finite_cpp(L)) {
    Rcpp::stop("L must contain only finite values.");
  }
}

// One projected-gradient optimization from a starting transformation. The step is doubled
// between iterations and halved within the line search.
//
// `fwindow` selects the line-search acceptance rule. With the default `fwindow == 0` the line
// search is monotone: a candidate is accepted only if it satisfies the sufficient-decrease
// (Armijo) condition, or, as a numerical fallback, if it at least decreases the objective after
// all step halvings are exhausted; otherwise the optimization terminates. The lowest-objective
// iterate is necessarily the final one, which is reported, and convergence is the projected
// gradient falling below `eps`. This is the rule used by every smooth rotation criterion (and the
// Procrustes objective), whose gradient vanishes at the optimum.
//
// With `fwindow >= 1` the acceptance is non-monotone (Grippo, Lampariello & Lucidi, 1986): a
// candidate is accepted if it sufficiently decreases the largest objective over the last
// `fwindow` iterations (rather than the current one), and when no step is accepted the
// smallest-step candidate is taken unconditionally instead of terminating. This lets the optimizer
// cross the kinks of a criterion whose gradient does not vanish at the optimum (the simplimax
// criterion reselects the k smallest squared loadings at every evaluation, so it is only piecewise
// smooth and a strictly monotone search stalls short of the minimum). Because such a search may end
// on an up-step and its projected-gradient norm need not reach `eps`, the non-monotone path instead
// (a) reports the lowest-objective iterate it visited and (b) reports convergence when that best
// objective has stopped improving for `fwindow` consecutive iterations by the end of the run (it
// has stalled at the optimum) -- the projected-gradient tolerance is not the right convergence
// signal there. The search still runs to `maxit` (its non-monotone steps can cross a kink to a
// deeper basin). The `fwindow == 0` path is exactly the monotone behavior above; `fwindow` changes
// the acceptance test, the reported iterate, and the convergence flag, not the search trajectory.
template <typename Manifold>
GpfFit run_single_gpf_fit(const Manifold& manifold,
                          arma::mat Tmat,
                          double eps,
                          int maxit,
                          int max_line_search,
                          double step0,
                          int fwindow = 0) {
  Tmat = manifold.normalize_start(Tmat);

  // Record per-iteration diagnostics in a growing buffer rather than preallocating a
  // (maxit + 1) x 4 table. maxit is a documented, user-forwarded control that can be
  // large (e.g. EFA_AVERAGE passes maxit = 5e4) while the optimization typically
  // converges in far fewer iterations, so an eager allocation would reserve -- and
  // NA-fill -- a table that is then trimmed away, and `maxit + 1` could overflow.
  std::vector<double> diag_rows;  // row-major: iter, f, log10(s), step per recorded row
  diag_rows.reserve(4 * 64);
  int filled = 0;
  bool convergence = false;
  bool line_search_failed = false;
  double al = step0;

  GpfState state = manifold.compute(Tmat);

  // Non-monotone path (fwindow >= 1) bookkeeping: the lowest-objective iterate visited (reported in
  // place of the final, possibly-up-step iterate) and a stall counter for the convergence flag. The
  // counter resets whenever the best objective drops by more than a relative tolerance; once it
  // stays flat for `fwindow` iterations the search has stalled at the (piecewise-smooth) optimum.
  // Only the winning transformation needs a matrix copy (the reported value/validity are scalars),
  // and only when fwindow >= 1, so the monotone path takes no extra copy.
  arma::mat best_T;
  double best_f = state.f;
  bool best_valid = state.valid;
  double stall_ref_f = state.f;
  int iters_since_improve = 0;
  if (fwindow >= 1) {
    best_T = Tmat;
  }

  for (int iter = 0; iter <= maxit; ++iter) {
    diag_rows.push_back(static_cast<double>(iter));
    diag_rows.push_back(state.f);
    diag_rows.push_back(std::log10(std::max(state.s, std::numeric_limits<double>::epsilon())));
    diag_rows.push_back(al);
    ++filled;

    if (fwindow >= 1 && state.valid) {
      // Track the lowest-objective iterate (for reporting) and stalling (for the convergence flag).
      if (state.f < best_f) {
        best_f = state.f;
        best_valid = state.valid;
        best_T = Tmat;
      }
      // The improvement threshold scales with the current best level (unit floor), so the stall
      // test is appropriate to the criterion's magnitude rather than pinned to the starting value.
      const double ftol = 1e-8 * std::max(1.0, std::abs(stall_ref_f));
      if (state.f < stall_ref_f - ftol) {
        stall_ref_f = state.f;
        iters_since_improve = 0;
      } else {
        ++iters_since_improve;
      }
    }

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

    // Non-monotone reference: the largest objective over the last `fwindow` recorded
    // iterations (the f-values live in column 1 of the diagnostics buffer). With the default
    // fwindow == 0 this is just the current objective, so the acceptance test below reduces
    // exactly to the monotone sufficient-decrease test.
    double target_f = state.f;
    if (fwindow >= 1) {
      const int lo = std::max(0, filled - fwindow);
      for (int r = lo; r < filled; ++r) {
        target_f = std::max(target_f, diag_rows[4 * static_cast<std::size_t>(r) + 1]);
      }
    }

    bool armijo_improved = false;

    // Monotone fallback bookkeeping (fwindow == 0 only): the best strict decrease tried. The
    // matrices are assigned only when a candidate strictly decreases the objective (any_decrease),
    // so they need no per-iteration initial copy.
    bool any_decrease = false;
    double best_decrease_value = state.f;
    double best_decrease_step = al;
    arma::mat best_decrease_T;
    GpfState best_decrease_state;

    // Non-monotone fallback bookkeeping (fwindow >= 1 only): the smallest-step valid candidate (the
    // line search halves the step, so the last valid candidate is the smallest-step one), accepted
    // unconditionally when no step satisfies the sufficient-decrease test.
    bool have_last_valid = false;
    double last_valid_step = al;
    arma::mat last_valid_T;
    GpfState last_valid_state;

    double trial_step = al;
    for (int i = 0; i <= max_line_search; ++i) {
      arma::mat X = Tmat - trial_step * state.Gp;
      arma::mat T_candidate;

      if (manifold.retract(X, T_candidate)) {
        GpfState cand = manifold.compute(T_candidate);

        if (cand.valid) {
          if (fwindow >= 1) {
            have_last_valid = true;
            last_valid_step = trial_step;
            last_valid_T = T_candidate;
            last_valid_state = cand;
          } else if (cand.f < best_decrease_value) {
            any_decrease = true;
            best_decrease_value = cand.f;
            best_decrease_step = trial_step;
            best_decrease_T = T_candidate;
            best_decrease_state = cand;
          }

          if ((target_f - cand.f) > 0.5 * state.s * state.s * trial_step) {
            Tmat = T_candidate;
            state = cand;
            al = trial_step;
            armijo_improved = true;
            break;
          }
        }
      }

      trial_step = trial_step / 2.0;
    }

    if (!armijo_improved) {
      if (fwindow >= 1) {
        // GPArotation's GPFoblq/GPForth take the smallest-step candidate even when it does not
        // decrease the (non-monotone) objective, so the optimizer can step across a kink rather
        // than stop. Terminate only if the line search produced no valid candidate at all.
        if (have_last_valid) {
          Tmat = last_valid_T;
          state = last_valid_state;
          al = last_valid_step;
        } else {
          line_search_failed = true;
          break;
        }
      } else if (any_decrease) {
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

  // For the non-monotone path, a search that exhausted maxit while its best objective stopped
  // improving has stalled at the (piecewise-smooth) optimum -- report it as converged, since the
  // projected-gradient tolerance is not the right convergence signal there. The search itself still
  // runs to maxit (the non-monotone steps can cross a kink to a deeper basin, so terminating early
  // would forfeit the lower minima those crossings find). A genuine line-search failure is never
  // reclassified as converged.
  if (fwindow >= 1 && !line_search_failed && iters_since_improve >= fwindow) {
    convergence = true;
  }

  // filled >= 1 (iteration 0 is always recorded before any break), so the table has at
  // least one row. The buffer holds exactly `filled` rows of the four diagnostics.
  arma::mat table(filled, 4);
  for (int r = 0; r < filled; ++r) {
    const std::size_t base = 4 * static_cast<std::size_t>(r);
    table(r, 0) = diag_rows[base + 0];
    table(r, 1) = diag_rows[base + 1];
    table(r, 2) = diag_rows[base + 2];
    table(r, 3) = diag_rows[base + 3];
  }

  // The non-monotone path reports its lowest-objective iterate (best_*); the monotone path reports
  // its final iterate (which is also its lowest, by construction). For the non-monotone path the
  // reported iterate may therefore be earlier than the last row of the diagnostics `table` (and
  // than `iterations`), which always describe the full trajectory. The table is diagnostic only and
  // is not surfaced by the rotation entry points -- those read T/value/convergence/valid.
  const bool nonmono = (fwindow >= 1);

  GpfFit out;
  out.T = nonmono ? best_T : Tmat;
  out.Table = table;
  out.value = nonmono ? best_f : state.f;
  out.convergence = convergence;
  out.valid = nonmono ? best_valid : state.valid;
  out.line_search_failed = line_search_failed;
  out.iterations = filled - 1;
  return out;
}

// Run the primary fit from T_primary, then optimize the requested random starts and keep the best
// fit by objective value. With the default `full_multistart == false`, extra starts use a
// two-stage strategy to bound runtime: cheap objective screening, short triage optimization, then
// full optimization only for starts that improve on the incumbent by at least triage_improve_tol
// (appropriate when the rational start usually lies in the global basin, as for the smooth
// criteria). With `full_multistart == true`, every random start is fully optimized once and the
// screen/triage controls are ignored -- the right strategy for a strongly multimodal criterion
// whose early-iteration screen is uninformative (the simplimax criterion). Random orthogonal
// starts are drawn serially with R::rnorm.
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
                              double triage_improve_tol,
                              int fwindow = 0,
                              bool full_multistart = false) {
  const arma::uword k = T_primary.n_cols;

  GpfSummary summary;
  GpfFit best_fit = run_single_gpf_fit(manifold, T_primary, eps, maxit, max_line_search,
                                       step0, fwindow);
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

  if (random_starts > 0 && full_multistart) {
    // Full multistart: fully optimize every random start once and keep the best. This is the
    // standard remedy for criteria whose rational (identity) start is not in the global basin
    // (Kiers, 1994; Browne, 2001); the screen-and-triage heuristic below assumes it is, so its
    // early-iteration objective screen is uninformative for them. Random orthogonal starts are
    // drawn serially with R::rnorm in the same order as the screen loop, so the RNG stream is
    // consumed identically. Only the per-start all_* diagnostics are recorded here; the screen
    // diagnostics (screen_values, screen_start_indices, n_triaged) describe the screen-and-triage
    // path and stay empty/zero in this mode (the rotation entry points read only best_fit).
    for (int s = 0; s < random_starts; ++s) {
      arma::mat Tstart = random_orthogonal_start_cpp(k);
      GpfFit fit = run_single_gpf_fit(manifold, Tstart, eps, maxit, max_line_search, step0,
                                      fwindow);
      ++n_fully_optimized;
      summary.all_values.push_back(fit.value);
      summary.all_converged.push_back(fit.convergence ? 1 : 0);
      summary.all_iterations.push_back(fit.iterations);
      summary.all_start_indices.push_back(s + 2);
      if (gpf_fit_is_better(fit, best_fit)) {
        best_fit = fit;
        best_start_index = s + 2;
      }
    }
  } else if (random_starts > 0) {
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
        manifold, cand.Tstart, eps, triage_maxit, max_line_search, step0, fwindow
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
          manifold, triage_fit.T, eps, maxit, max_line_search, step0, fwindow
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
