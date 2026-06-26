#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include "gpf_engine.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

using namespace Rcpp;
using namespace arma;

// `normalize_cols_cpp` (column normalization onto diag(t(T) %*% T) = 1) and
// `inverse_checked_cpp` (rejecting checked inverse) are shared leaf helpers defined in
// gpf_common.h, included via gpf_engine.h above.

// Oblique Procrustes manifold for the shared gradient-projection engine. The
// transformation T lives on the column-normalized manifold diag(t(T) %*% T) = 1; the
// rotated loadings are L = A %*% solve(t(T)). The objective is the k x k form
// 0.5 * (tr(U' S U) - 2 tr(U' C) + tr(B'B)) of 0.5 * ||A solve(t(T)) - B||_F^2, with
// U = solve(t(T)), S = A'A, C = A'B. A non-invertible T yields an infinite (invalid)
// objective and is rejected rather than evaluated through a pseudo-inverse. U is not
// carried in the iteration state; it is reconstructed from the winning T in
// finalize_oblique().
struct ProcrustesManifold {
  const arma::mat& S;
  const arma::mat& C;
  double BtB;

  GpfState compute(const arma::mat& Tmat) const {
    arma::mat invT;
    if (!inverse_checked_cpp(Tmat, invT)) {
      return gpf_invalid_state();
    }

    arma::mat U = invT.t();
    arma::mat SU = S * U;
    arma::mat GU = SU - C;

    double f = 0.5 * (accu(U % SU) - 2.0 * accu(U % C) + BtB);
    arma::mat G = -U * GU.t() * U;

    arma::mat Gp;
    double s;
    oblique_tangent_project(Tmat, G, Gp, s);

    GpfState out;
    out.f = f;
    out.s = s;
    out.Gp = Gp;
    out.valid = std::isfinite(f) && std::isfinite(s) && all_finite_cpp(U) &&
      all_finite_cpp(G) && all_finite_cpp(Gp);

    if (!out.valid) {
      return gpf_invalid_state();
    }
    return out;
  }

  // The column-normalization projection always returns a candidate on the manifold;
  // a singular candidate is rejected later by the inverse check in compute().
  bool retract(const arma::mat& X, arma::mat& T_out) const {
    T_out = normalize_cols_cpp(X);
    return true;
  }

  arma::mat normalize_start(const arma::mat& Tmat) const {
    return normalize_cols_cpp(Tmat);
  }
};

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

// Kaiser row normalization shared by both entry points: scale each row of the
// loadings A_work by its inverse row norm, returning the weights so the rotated
// loadings can be back-transformed. The target is left unnormalized, matching
// GPArotation::targetQ(normalize = TRUE), which normalizes only the loadings. A
// zero/non-finite row norm is floored to 1 to leave that row unchanged.
static arma::vec kaiser_normalize_rows(arma::mat& A_work) {
  arma::vec W = sqrt(sum(square(A_work), 1));
  for (arma::uword i = 0; i < W.n_elem; ++i) {
    if (!std::isfinite(W(i)) || W(i) < 1e-15) {
      W(i) = 1.0;
    }
  }
  A_work.each_col() /= W;
  return W;
}

// `finalize_oblique()` (reconstruct the rotated loadings and factor correlations from a
// fitted transformation) is a shared inline in gpf_engine.h, used by this Procrustes path
// and the criterion rotations in rotate.cpp.

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
//' Non-invertible candidate transformations are rejected rather than evaluated
//' through a pseudo-inverse.
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
//' @param normalize Logical scalar. If `TRUE`, apply Kaiser normalization to the
//'   loadings (only) before rotation and reverse it afterwards; the target is left
//'   unnormalized, matching `GPArotation::targetQ(normalize = TRUE)`.
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
//' Browne, M. W. (2001). An overview of analytic rotation in exploratory factor
//' analysis. *Multivariate Behavioral Research*, 36, 111-150.
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
  validate_gpf_scalars(eps, maxit, max_line_search, step0, random_starts,
                       screen_keep, triage_maxit, triage_improve_tol);

  const unsigned int k = A.n_cols;

  arma::mat A_work = A;
  arma::vec W;

  if (normalize) {
    // normalize only the loadings; the target B stays raw (targetQ parity)
    W = kaiser_normalize_rows(A_work);
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

  arma::mat C = A_work.t() * B;
  double BtB = accu(B % B);
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

  ProcrustesManifold manifold{S, C, BtB};
  GpfSummary summary = run_gpf_multistart(
    manifold, T_primary, eps, maxit, max_line_search, step0,
    random_starts, screen_keep, triage_maxit, triage_improve_tol
  );
  const GpfFit& best_fit = summary.best_fit;

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
    Rcpp::Named("kappa_T") = cond_checked_cpp(best_fit.T),
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
//' @param normalize Logical scalar. If `TRUE`, apply Kaiser normalization to the
//'   loadings (only) before rotation and reverse it afterwards, leaving the target
//'   unnormalized (ignored for single-factor slices).
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
  validate_gpf_scalars(eps, maxit, max_line_search, step0, random_starts,
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
    // The target B is raw in both the normalized and unnormalized objective (only
    // the loadings are Kaiser-normalized), so B'B is shared across all slices;
    // hoist it rather than rebuilding it per replicate.
    const double BtB_shared = accu(B % B);

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
        bool warm_ok;

        // The target B is raw in both branches, so BtB = BtB_shared throughout.
        if (normalize) {
          A_work = A_i;
          W = kaiser_normalize_rows(A_work);   // normalize only the loadings
          S = A_work.t() * A_work;
          C = A_work.t() * B;
          warm_ok = orthogonal_procrustes_T_cpp(A_work, B, T_primary);
        } else {
          S = A_i.t() * A_i;
          C = A_i.t() * B;
          warm_ok = orthogonal_procrustes_T_cpp(A_i, B, T_primary);
        }
        if (!warm_ok) {
          continue;
        }

        ProcrustesManifold manifold{S, C, BtB_shared};
        GpfSummary summary = run_gpf_multistart(
          manifold, T_primary, eps, maxit, max_line_search, step0,
          random_starts, screen_keep, triage_maxit, triage_improve_tol
        );
        const GpfFit& best_fit = summary.best_fit;
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
