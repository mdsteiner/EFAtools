// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <limits>

// Native two-step (conditional) polychoric correlation estimator. Thresholds are
// fixed once per variable from the marginal cumulative proportions; each pairwise
// correlation then maximises the bivariate-ordinal log-likelihood over rho with
// those thresholds held fixed (Olsson, 1979, two-step estimator; matches
// polycor::polychor(ML = FALSE) and the psych two-step).
//
// Each ordinal cell probability is the integral of the bivariate normal density over
// its (possibly half-infinite) threshold rectangle. Rather than forming it by
// differencing the bivariate CDF at the four corners - which suffers catastrophic
// cancellation for the near-impossible tail cells of highly correlated items and can
// even return negatives - the rectangle is computed by conditioning on X:
//   P(a0<X<=a1, b0<Y<=b1) = int_{a0}^{a1} phi(x) [Phi((b1-rho x)/s) - Phi((b0-rho x)/s)] dx,
// with s = sqrt(1-rho^2) (Plackett, 1954; Genz & Bretz, 2009). The 1-D integral is a
// sum of non-negative terms, so it keeps full relative accuracy for tiny cells with no
// cancellation. It is evaluated by Gauss-Legendre quadrature, infinite X-cuts clamped to
// a finite range where phi is negligible, and infinite Y-cuts handled directly by the
// normal CDF (Phi(+/-Inf) = 1/0). The per-pair correlation is found by Brent's (1973)
// method. Every per-pair quantity is pure C++/Armadillo (std::erfc for the normal CDF),
// so the pair loop is allocation-light and safe to run under OpenMP.
//
// References:
//   Olsson, U. (1979). Maximum likelihood estimation of the polychoric correlation
//     coefficient. Psychometrika, 44, 443-460.
//   Plackett, R. L. (1954). A reduction formula for normal multivariate integrals.
//     Biometrika, 41, 351-360.
//   Genz, A., & Bretz, F. (2009). Computation of Multivariate Normal and t Probabilities.
//     Springer.
//   Brent, R. P. (1973). Algorithms for Minimization without Derivatives. Prentice-Hall.
//   Higham, N. J. (2002). Computing the nearest correlation matrix - a problem from
//     finance. IMA Journal of Numerical Analysis, 22, 329-343.

static const double POLY_INF = std::numeric_limits<double>::infinity();
static const double POLY_INV_SQRT_2PI = 0.39894228040143267794;  // 1/sqrt(2*pi)
static const double POLY_INV_SQRT_2 = 0.70710678118654752440;    // 1/sqrt(2)
// Beyond this |x| the standard normal density is below ~1e-16; used to clamp the
// infinite outer (X) integration limits to a finite range without loss of accuracy.
static const double POLY_CLAMP = 8.5;
// Gauss-Legendre node count for the rectangle quadrature. The integrand is a smooth
// Gaussian over each threshold band, so this converges quickly; 16 nodes reproduce the
// reference estimators to within ~1e-5 and higher orders do not change the result.
static const int POLY_GL_N = 16;

// Standard normal CDF via the complementary error function (thread-safe, allocation-free,
// accurate to floating-point precision). Phi(+/-Inf) evaluates cleanly to 1 / 0.
static inline double std_norm_cdf(double z) {
  return 0.5 * std::erfc(-z * POLY_INV_SQRT_2);
}
static inline double std_norm_pdf(double z) {
  return POLY_INV_SQRT_2PI * std::exp(-0.5 * z * z);
}

// Gauss-Legendre nodes/weights on [-1, 1] (Golub-Welsch via Newton on the Legendre roots).
// Computed once per call, used read-only inside the parallel pair loop.
static void gauss_legendre(int n, std::vector<double>& x, std::vector<double>& w) {
  const double pi = 3.14159265358979323846;
  x.assign(n, 0.0);
  w.assign(n, 0.0);
  int m = (n + 1) / 2;
  for (int i = 0; i < m; i++) {
    double z = std::cos(pi * (i + 0.75) / (n + 0.5));
    double pp = 0.0, z1;
    do {
      double p1 = 1.0, p2 = 0.0;
      for (int j = 1; j <= n; j++) {
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
      }
      pp = n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    } while (std::abs(z - z1) > 1e-15);
    double ww = 2.0 / ((1.0 - z * z) * pp * pp);
    x[i] = -z;       w[i] = ww;
    x[n - 1 - i] = z; w[n - 1 - i] = ww;
  }
}

// Rectangle probability P(a0 < X <= a1, b0 < Y <= b1; rho) via 1-D conditioning on X.
// Infinite Y-cuts are handled by std_norm_cdf directly; infinite X-cuts are clamped to
// +/-POLY_CLAMP. Cancellation-free and non-negative by construction.
static double bvn_rect(double a0, double a1, double b0, double b1, double rho,
                       const std::vector<double>& gx, const std::vector<double>& gw) {
  double lo = a0 < -POLY_CLAMP ? -POLY_CLAMP : a0;
  double hi = a1 >  POLY_CLAMP ?  POLY_CLAMP : a1;
  if (hi <= lo) return 0.0;
  double s = std::sqrt(1.0 - rho * rho);
  double mid = 0.5 * (lo + hi), half = 0.5 * (hi - lo);
  double acc = 0.0;
  for (std::size_t k = 0; k < gx.size(); k++) {
    double x = mid + half * gx[k];
    double inner = std_norm_cdf((b1 - rho * x) / s) - std_norm_cdf((b0 - rho * x) / s);
    acc += gw[k] * std_norm_pdf(x) * inner;
  }
  return half * acc;
}

// Brent's (1973) derivative-free minimiser on [lo, hi] (the parabolic-interpolation /
// golden-section method also used by R's optimise()). Pure C++ and deterministic.
template <typename F>
static double brent_min(F f, double lo, double hi, double tol, int max_iter) {
  const double gold = 0.3819660112501051;  // (3 - sqrt(5)) / 2
  const double zeps = 1e-12;
  double a = lo, b = hi;
  double x = a + gold * (b - a);
  double w = x, v = x;
  double fx = f(x), fw = fx, fv = fx;
  double d = 0.0, e = 0.0;

  for (int iter = 0; iter < max_iter; iter++) {
    double xm = 0.5 * (a + b);
    double tol1 = tol * std::abs(x) + zeps;
    double tol2 = 2.0 * tol1;
    if (std::abs(x - xm) <= (tol2 - 0.5 * (b - a))) break;

    bool use_golden = true;
    if (std::abs(e) > tol1) {
      double rr = (x - w) * (fx - fv);
      double q = (x - v) * (fx - fw);
      double pp = (x - v) * q - (x - w) * rr;
      q = 2.0 * (q - rr);
      if (q > 0.0) pp = -pp;
      q = std::abs(q);
      double etemp = e;
      e = d;
      if (std::abs(pp) < std::abs(0.5 * q * etemp) &&
          pp > q * (a - x) && pp < q * (b - x)) {
        d = pp / q;
        double u = x + d;
        if (u - a < tol2 || b - u < tol2) {
          d = (xm - x >= 0.0) ? std::abs(tol1) : -std::abs(tol1);
        }
        use_golden = false;
      }
    }
    if (use_golden) {
      e = (x >= xm) ? (a - x) : (b - x);
      d = gold * e;
    }

    double u = (std::abs(d) >= tol1) ? (x + d)
                                     : (x + ((d >= 0.0) ? std::abs(tol1) : -std::abs(tol1)));
    double fu = f(u);
    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }
  return x;
}

// Two-step polychoric correlation for one variable pair, thresholds fixed. Builds the
// contingency table from the pairwise-complete rows, optionally applies the empty-cell
// continuity correction, then minimises the negative log-likelihood over rho. Pure C++
// (OpenMP-safe); reuses caller-owned scratch buffers so the hot path does not allocate.
// Returns NA_REAL when the pair has no overlapping complete cases (the R wrapper turns a
// resulting NA into a classed condition).
static double polychoric_pair(const int* xp, int nrow, int ci, int cj,
                              const std::vector<double>& tau_i,
                              const std::vector<double>& tau_j,
                              int Ki, int Kj, double correct,
                              const std::vector<double>& gx, const std::vector<double>& gw,
                              std::vector<double>& tab, std::vector<double>& rowp,
                              std::vector<double>& colcdf, std::vector<double>& rcut,
                              std::vector<double>& ccut, std::vector<double>& xnode,
                              std::vector<double>& wpdf) {
  const int n = (int) gx.size();
  std::fill(tab.begin(), tab.begin() + (std::size_t) Ki * Kj, 0.0);
  double total = 0.0;
  const int* coli = xp + (std::size_t) ci * nrow;
  const int* colj = xp + (std::size_t) cj * nrow;
  for (int r = 0; r < nrow; r++) {
    int a = coli[r];
    int b = colj[r];
    if (a == NA_INTEGER || b == NA_INTEGER) continue;  // pairwise-complete
    tab[(std::size_t) a * Kj + b] += 1.0;
    total += 1.0;
  }
  if (total <= 0.0) return NA_REAL;  // no overlapping complete cases for this pair

  // Empty-cell continuity correction: add `correct` pseudo-counts to zero cells (mirrors
  // psych's correct= for sparse tables). correct = 0 leaves the table untouched.
  if (correct > 0.0) {
    std::size_t ncell = (std::size_t) Ki * Kj;
    for (std::size_t idx = 0; idx < ncell; idx++) {
      if (tab[idx] == 0.0) tab[idx] = correct;
    }
  }

  // Pad the fixed thresholds with +/-Inf so each ordinal cell is a (half-infinite) rectangle.
  rcut[0] = -POLY_INF; rcut[Ki] = POLY_INF;
  for (int k = 0; k < Ki - 1; k++) rcut[k + 1] = tau_i[k];
  ccut[0] = -POLY_INF; ccut[Kj] = POLY_INF;
  for (int k = 0; k < Kj - 1; k++) ccut[k + 1] = tau_j[k];

  // Precompute the rho-INDEPENDENT part of the quadrature once per pair: for each X (row)
  // band, the clamped limits, the abscissae, and the weight*phi(x) factors depend only on the
  // fixed thresholds, so they are reused across every Brent evaluation instead of being
  // rebuilt each time. A band lying entirely beyond +/-POLY_CLAMP is marked with zero weight.
  for (int a = 0; a < Ki; a++) {
    double lo = rcut[a]     < -POLY_CLAMP ? -POLY_CLAMP : rcut[a];
    double hi = rcut[a + 1] >  POLY_CLAMP ?  POLY_CLAMP : rcut[a + 1];
    std::size_t base = (std::size_t) a * n;
    if (hi > lo) {
      double mid = 0.5 * (lo + hi), half = 0.5 * (hi - lo);
      for (int k = 0; k < n; k++) {
        double x = mid + half * gx[k];
        xnode[base + k] = x;
        wpdf[base + k] = gw[k] * half * std_norm_pdf(x);
      }
    } else {
      for (int k = 0; k < n; k++) wpdf[base + k] = 0.0;
    }
  }

  // Smallest positive double, used only to keep log() finite for an utterly negligible
  // (underflowing) cell; the conditioning integral is otherwise non-negative and accurate.
  const double p_floor = std::numeric_limits<double>::min();

  // Negative log-likelihood at rho: only the inner (column) normal CDFs depend on rho. The
  // per-row inner CDFs are shared across the cells of that row.
  auto neg_loglik = [&](double rho) {
    double s = std::sqrt(1.0 - rho * rho);
    double val = 0.0;
    for (int a = 0; a < Ki; a++) {
      std::size_t base = (std::size_t) a * n;
      for (int b = 0; b < Kj; b++) rowp[b] = 0.0;
      for (int k = 0; k < n; k++) {
        double w = wpdf[base + k];
        if (w == 0.0) continue;  // band beyond the clamp contributes nothing
        double x = xnode[base + k];
        for (int b = 0; b <= Kj; b++) colcdf[b] = std_norm_cdf((ccut[b] - rho * x) / s);
        for (int b = 0; b < Kj; b++) rowp[b] += w * (colcdf[b + 1] - colcdf[b]);
      }
      for (int b = 0; b < Kj; b++) {
        double nab = tab[(std::size_t) a * Kj + b];
        if (nab == 0.0) continue;
        double p = rowp[b];
        if (p < p_floor) p = p_floor;
        val -= nab * std::log(p);
      }
    }
    return val;
  };

  // Optimise on (-1, 1) bounded away from the singular endpoints, as in polycor.
  return brent_min(neg_loglik, -0.9999, 0.9999, 1e-8, 200);
}

// Internal C++ backend for the two-step polychoric/tetrachoric matrix, called from
// .polychoric(), which owns the user-facing validation and classed conditions.
// [[Rcpp::export(.polychoric_cpp)]]
Rcpp::List polychoric_cpp(Rcpp::IntegerMatrix x, std::string acov,
                          double correct, bool nearest_pd, int n_threads) {
  // acov is accepted for forward compatibility; only the matrix is produced here.
  if (acov != "none") {
    Rcpp::stop("Only acov = \"none\" is currently supported.");
  }
  const int nrow = x.nrow();
  const int p = x.ncol();
  if (p < 2) {
    Rcpp::stop("At least two variables are required for a correlation matrix.");
  }
  if (n_threads < 1) n_threads = 1;

  const int* xp = x.begin();

  // Gauss-Legendre rule for the rectangle quadrature, computed once on the main thread.
  std::vector<double> gx, gw;
  gauss_legendre(POLY_GL_N, gx, gw);

  // Category counts and thresholds, computed once per variable on the main thread
  // (R::qnorm is not called inside the parallel region below). The codes are 0-based
  // and consecutive (the R wrapper recodes them), so the number of categories is
  // max(code) + 1 and no marginal proportion lands exactly on 0.
  std::vector<int> Kcat(p, 0);
  std::vector< std::vector<double> > tau(p);
  int Kmax = 2;
  for (int j = 0; j < p; j++) {
    const int* col = xp + (std::size_t) j * nrow;
    int maxc = -1;
    int n_obs = 0;  // bounded by nrow (an int)
    for (int r = 0; r < nrow; r++) {
      int v = col[r];
      if (v == NA_INTEGER) continue;
      if (v > maxc) maxc = v;
      n_obs++;
    }
    int K = maxc + 1;
    Kcat[j] = K;
    if (K > Kmax) Kmax = K;

    std::vector<double> freq(K > 0 ? K : 1, 0.0);
    for (int r = 0; r < nrow; r++) {
      int v = col[r];
      if (v == NA_INTEGER) continue;
      freq[v] += 1.0;
    }

    std::vector<double> th(K > 1 ? K - 1 : 0);
    double cum = 0.0;
    for (int k = 0; k < K - 1; k++) {
      cum += freq[k];
      th[k] = R::qnorm(cum / (double) n_obs, 0.0, 1.0, 1, 0);
    }
    tau[j] = th;
  }

  // Flat list of upper-triangle pairs, parallelised over OpenMP. Each pair writes a
  // distinct rho slot, so the result is independent of the thread count (no reduction
  // races) and therefore bit-for-bit deterministic.
  std::vector<int> pair_i, pair_j;
  pair_i.reserve((std::size_t) p * (p - 1) / 2);
  pair_j.reserve((std::size_t) p * (p - 1) / 2);
  for (int i = 0; i < p; i++) {
    for (int j = i + 1; j < p; j++) {
      pair_i.push_back(i);
      pair_j.push_back(j);
    }
  }
  const int npairs = (int) pair_i.size();
  std::vector<double> rho((std::size_t) npairs, 0.0);

  #pragma omp parallel num_threads(n_threads)
  {
    // Thread-local scratch reused across this thread's pairs (no allocation in the loop).
    std::vector<double> tab((std::size_t) Kmax * Kmax);
    std::vector<double> rowp(Kmax);
    std::vector<double> colcdf(Kmax + 1);
    std::vector<double> rcut(Kmax + 1), ccut(Kmax + 1);
    std::vector<double> xnode((std::size_t) Kmax * POLY_GL_N);
    std::vector<double> wpdf((std::size_t) Kmax * POLY_GL_N);

    #pragma omp for schedule(dynamic)
    for (int t = 0; t < npairs; t++) {
      int i = pair_i[t];
      int j = pair_j[t];
      rho[t] = polychoric_pair(xp, nrow, i, j, tau[i], tau[j], Kcat[i], Kcat[j], correct,
                               gx, gw, tab, rowp, colcdf, rcut, ccut, xnode, wpdf);
    }
  }

  arma::mat Rmat(p, p, arma::fill::eye);
  for (int t = 0; t < npairs; t++) {
    double v = rho[t];
    Rmat(pair_i[t], pair_j[t]) = v;
    Rmat(pair_j[t], pair_i[t]) = v;
  }

  // Optional nearest-PD projection, gated by `nearest_pd` so callers can choose how it
  // composes with downstream smoothing. The default path needs only the smallest eigenvalue,
  // so it uses the cheaper values-only eigendecomposition; the eigenvectors are computed only
  // on the rare projection branch. A non-finite matrix (a pair with no overlapping complete
  // cases) is left untouched here; the R wrapper raises a classed condition on the NA.
  bool pd_adjusted = false;
  if (Rmat.is_finite()) {
    arma::vec eigval;
    if (!arma::eig_sym(eigval, Rmat)) {
      Rcpp::stop("Eigendecomposition of the polychoric correlation matrix failed.");
    }
    if (nearest_pd && eigval.min() < std::numeric_limits<double>::epsilon()) {
      // Clip the eigenvalues to a small positive floor and rescale to a unit
      // diagonal (a simple nearest-correlation projection; Higham, 2002).
      arma::vec ev;
      arma::mat eigvec;
      if (!arma::eig_sym(ev, eigvec, Rmat)) {
        Rcpp::stop("Eigendecomposition of the polychoric correlation matrix failed.");
      }
      const double floor_eps = 1e-8;
      ev.transform([floor_eps](double e) { return e < floor_eps ? floor_eps : e; });
      arma::mat Rp = eigvec * arma::diagmat(ev) * eigvec.t();
      arma::vec d = arma::sqrt(Rp.diag());
      arma::mat Dinv = arma::diagmat(1.0 / d);
      Rmat = Dinv * Rp * Dinv;
      Rmat = 0.5 * (Rmat + Rmat.t());
      Rmat.diag().ones();
      pd_adjusted = true;
    }
  }

  Rcpp::List thr(p);
  for (int j = 0; j < p; j++) {
    thr[j] = Rcpp::NumericVector(tau[j].begin(), tau[j].end());
  }

  return Rcpp::List::create(
    Rcpp::Named("R") = Rmat,
    Rcpp::Named("thresholds") = thr,
    Rcpp::Named("pd_adjusted") = pd_adjusted,
    Rcpp::Named("acov") = R_NilValue);
}

// Bivariate normal rectangle probability P(a0 < X <= a1, b0 < Y <= b1; rho), exposed only
// so the test suite can cross-check the cell-probability integral against mnormt::sadmvn /
// mvtnorm::pmvnorm. Arguments are recycled to the longest length (standard R recycling).
// [[Rcpp::export(.bvn_rect_cpp)]]
Rcpp::NumericVector bvn_rect_cpp(Rcpp::NumericVector a0, Rcpp::NumericVector a1,
                                 Rcpp::NumericVector b0, Rcpp::NumericVector b1,
                                 Rcpp::NumericVector rho) {
  R_xlen_t na0 = a0.size(), na1 = a1.size(), nb0 = b0.size(),
           nb1 = b1.size(), nr = rho.size();
  if (na0 == 0 || na1 == 0 || nb0 == 0 || nb1 == 0 || nr == 0) {
    return Rcpp::NumericVector(0);
  }
  R_xlen_t n = na0;
  if (na1 > n) n = na1;
  if (nb0 > n) n = nb0;
  if (nb1 > n) n = nb1;
  if (nr  > n) n = nr;

  std::vector<double> gx, gw;
  gauss_legendre(POLY_GL_N, gx, gw);
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; i++) {
    out[i] = bvn_rect(a0[i % na0], a1[i % na1], b0[i % nb0], b1[i % nb1],
                      rho[i % nr], gx, gw);
  }
  return out;
}
