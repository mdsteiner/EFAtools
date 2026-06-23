// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <limits>

// Native two-step (conditional) polychoric correlation estimator. Thresholds are
// fixed from the marginal cumulative proportions; each pairwise correlation then
// maximises the bivariate-ordinal log-likelihood over rho with those thresholds held
// fixed (Olsson, 1979, two-step estimator; matches polycor::polychor(ML = FALSE) and
// the psych two-step). On complete data the thresholds are the per-variable (full-column)
// marginals; when the data have missing values each pairwise correlation instead uses the
// thresholds of its own pairwise-complete cases, so the thresholds and the contingency
// table always come from the same cases.
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
// normal CDF (Phi(+/-Inf) = 1/0). The per-pair correlation maximises the log-likelihood by
// Fisher scoring: the score is the closed-form derivative of each cell probability obtained by
// differentiating the same conditioning integral with respect to rho (so it is consistent with
// the quadrature used for the probabilities), and the step is normalised by the expected
// information of the cell-count multinomial (Olsson, 1979). The step is capped away from the
// singular endpoints and backtracked to keep the likelihood monotone, with Brent's (1973)
// method as a fallback. Every per-pair quantity is pure C++/Armadillo (std::erfc for the normal
// CDF), so the pair loop is allocation-light and safe to run under OpenMP.
//
// References:
//   Olsson, U. (1979). Maximum likelihood estimation of the polychoric correlation
//     coefficient. Psychometrika, 44, 443-460.
//   Plackett, R. L. (1954). A reduction formula for normal multivariate integrals.
//     Biometrika, 41, 351-360.
//   Genz, A., & Bretz, F. (2009). Computation of Multivariate Normal and t Probabilities.
//     Springer.
//   Brent, R. P. (1973). Algorithms for Minimization without Derivatives. Prentice-Hall.
//   Rebonato, R., & Jackel, P. (2000). The most general methodology to create a valid
//     correlation matrix for risk management and option pricing purposes. Journal of Risk,
//     2(2), 17-27.
//   Wichura, M. J. (1988). Algorithm AS 241: The percentage points of the normal
//     distribution. Applied Statistics, 37, 477-484.

static const double POLY_INF = std::numeric_limits<double>::infinity();
static const double POLY_INV_SQRT_2PI = 0.39894228040143267794;  // 1/sqrt(2*pi)
static const double POLY_INV_SQRT_2 = 0.70710678118654752440;    // 1/sqrt(2)
// Beyond this |x| the standard normal density is below ~1e-16; used to clamp the
// infinite outer (X) integration limits to a finite range without loss of accuracy.
static const double POLY_CLAMP = 8.5;
// Gauss-Legendre node count for the rectangle quadrature. The integrand is a smooth Gaussian
// over each threshold band whose conditional transition has width sqrt(1 - rho^2), so 12 nodes
// reproduce the reference estimators to within ~1e-5 for |rho| up to about 0.95. Near |rho| = 1
// the transition narrows and a fixed 12-node rule under-resolves it, biasing the estimate, so a
// pair whose estimate exceeds POLY_REFINE_RHO is re-estimated with the finer POLY_GL_N_HI rule.
// The common, moderate-correlation case keeps the cheap 12-node rule and is unchanged.
static const int POLY_GL_N = 12;
static const int POLY_GL_N_HI = 96;
static const double POLY_REFINE_RHO = 0.95;

// Standard normal CDF via the complementary error function (thread-safe, allocation-free,
// accurate to floating-point precision). Phi(+/-Inf) evaluates cleanly to 1 / 0.
static inline double std_norm_cdf(double z) {
  return 0.5 * std::erfc(-z * POLY_INV_SQRT_2);
}
static inline double std_norm_pdf(double z) {
  return POLY_INV_SQRT_2PI * std::exp(-0.5 * z * z);
}

// Standard normal quantile via Wichura's (1988) AS 241 (accurate to ~1e-16), thread-safe and
// allocation-free so it can be called inside the OpenMP pair loop, where R::qnorm cannot. Used
// for the per-pair (pairwise-complete) thresholds under missing data; returns +/-Inf at p = 1/0
// (an empty marginal category, handled by the +/-POLY_CLAMP clamping downstream).
static inline double std_norm_quantile(double p) {
  // Endpoints first: the AS241 rational evaluates Inf/Inf -> NaN at p = 0/1, so an empty
  // pairwise marginal category (cumulative proportion exactly 0 or 1) must be sent to the
  // infinite cut the threshold padding and +/-POLY_CLAMP clamping expect (as R::qnorm does).
  if (p <= 0.0) return -POLY_INF;
  if (p >= 1.0) return POLY_INF;
  double q = p - 0.5, r, val;
  if (std::fabs(q) <= 0.425) {
    r = 0.180625 - q * q;
    val = q * (((((((2509.0809287301226727 * r + 33430.575583588128105) * r +
            67265.770927008700853) * r + 45921.953931549871457) * r +
            13731.693765509461125) * r + 1971.5909503065514427) * r +
            133.14166789178437745) * r + 3.387132872796366608) /
          (((((((5226.495278852854561 * r + 28729.085735721942674) * r +
            39307.89580009271061) * r + 21213.794301586595867) * r +
            5394.1960214247511077) * r + 687.1870074920579083) * r +
            42.313330701600911252) * r + 1.0);
  } else {
    r = (q < 0.0) ? p : 1.0 - p;
    r = std::sqrt(-std::log(r));
    if (r <= 5.0) {
      r -= 1.6;
      val = (((((((7.7454501427834140764e-4 * r + 0.0227238449892691845833) * r +
              0.24178072517745061177) * r + 1.27045825245236838258) * r +
              3.64784832476320460504) * r + 5.7694972214606914055) * r +
              4.6303378461565452959) * r + 1.42343711074968357734) /
            (((((((1.05075007164441684324e-9 * r + 5.475938084995344946e-4) * r +
              0.0151986665636164571966) * r + 0.14810397642748007459) * r +
              0.68976733498510000455) * r + 1.6763848301838038494) * r +
              2.05319162663775882187) * r + 1.0);
    } else {
      r -= 5.0;
      val = (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r +
              0.0012426609473880784386) * r + 0.026532189526576123093) * r +
              0.29656057182850489123) * r + 1.7848265399172913358) * r +
              5.4637849111641143699) * r + 6.6579046435011037772) /
            (((((((2.04426310338993978564e-15 * r + 1.4215117583164458887e-7) * r +
              1.8463183175100546818e-5) * r + 7.868691311456132591e-4) * r +
              0.0148753612908506148525) * r + 0.13692988092273580531) * r +
              0.59983220655588793769) * r + 1.0);
    }
    if (q < 0.0) val = -val;
  }
  return val;
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

// Accumulate one Gauss-Legendre node's contribution to a row of cell probabilities and
// their rho-derivatives: P_ab += w * P(Y in column band b | X = x) and dP_ab/drho += the
// rho-derivative of the same. The conditional band CDF and its derivative (Plackett, 1954:
// d/drho Phi((c - rho x)/s) = phi((c - rho x)/s) (rho c - x)/s^3; an infinite column cut
// contributes zero) are shared by the likelihood evaluation in polychoric_pair and the ACOV
// cell probabilities in poly_cell_probs, so the formula lives in exactly one place. cut_j is
// padded with +/-Inf; colcdf/coldcdf are scratch of length >= Kj+1.
static inline void poly_accum_node(double x, double w, double rho, double s, double s3,
                                   const double* cut_j, int Kj, double* prow, double* dprow,
                                   double* colcdf, double* coldcdf) {
  for (int b = 0; b <= Kj; b++) {
    double c = cut_j[b];
    double z = (c - rho * x) / s;
    colcdf[b] = std_norm_cdf(z);
    coldcdf[b] = std::isfinite(c) ? std_norm_pdf(z) * (rho * c - x) / s3 : 0.0;
  }
  for (int b = 0; b < Kj; b++) {
    prow[b]  += w * (colcdf[b + 1]  - colcdf[b]);
    dprow[b] += w * (coldcdf[b + 1] - coldcdf[b]);
  }
}

// Two-step polychoric correlation for one variable pair. Builds the contingency table from the
// pairwise-complete rows, fixes the thresholds (the shared per-variable tau_i/tau_j, or - when
// use_local is set, i.e. the data have missing values - this pair's own marginal thresholds),
// optionally applies the empty-cell continuity correction, then minimises the negative
// log-likelihood over rho. Pure C++ (OpenMP-safe); reuses caller-owned scratch buffers so the
// hot path does not allocate. Returns NA_REAL when the pair has no overlapping complete cases
// (the R wrapper turns a resulting NA into a classed condition).
static double polychoric_pair(const int* xp, int nrow, int ci, int cj,
                              const std::vector<double>& tau_i,
                              const std::vector<double>& tau_j,
                              int Ki, int Kj, double correct, bool use_local,
                              const std::vector<double>& gx, const std::vector<double>& gw,
                              std::vector<double>& tab, std::vector<double>& pmat,
                              std::vector<double>& dpmat, std::vector<double>& colcdf,
                              std::vector<double>& coldcdf, std::vector<double>& rcut,
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

  // Pad the thresholds with +/-Inf so each ordinal cell is a (half-infinite) rectangle. With
  // complete data the per-variable (full-column) thresholds tau_i/tau_j are shared by every
  // pair. With missing data each pair instead uses the thresholds of its own pairwise-complete
  // cases (the marginals of this contingency table, taken before any continuity correction):
  // the thresholds and the table then come from the same cases, as in polycor / psych. An empty
  // marginal category gives a +/-Inf cut, which the +/-POLY_CLAMP clamping turns into a
  // zero-width (zero-probability, zero-count) band.
  rcut[0] = -POLY_INF; rcut[Ki] = POLY_INF;
  ccut[0] = -POLY_INF; ccut[Kj] = POLY_INF;
  if (use_local) {
    double cum = 0.0;
    for (int a = 0; a < Ki - 1; a++) {
      double rs = 0.0;
      for (int b = 0; b < Kj; b++) rs += tab[(std::size_t) a * Kj + b];
      cum += rs;
      rcut[a + 1] = std_norm_quantile(cum / total);
    }
    cum = 0.0;
    for (int b = 0; b < Kj - 1; b++) {
      double cs = 0.0;
      for (int a = 0; a < Ki; a++) cs += tab[(std::size_t) a * Kj + b];
      cum += cs;
      ccut[b + 1] = std_norm_quantile(cum / total);
    }
  } else {
    for (int k = 0; k < Ki - 1; k++) rcut[k + 1] = tau_i[k];
    for (int k = 0; k < Kj - 1; k++) ccut[k + 1] = tau_j[k];
  }

  // Empty-cell continuity correction: add `correct` pseudo-counts to zero cells (mirrors
  // psych's correct= for sparse tables). correct = 0 leaves the table untouched. ntab is the
  // total table mass (including any pseudo-counts) - the multinomial size used in the
  // expected information.
  double ntab = total;
  if (correct > 0.0) {
    std::size_t ncell = (std::size_t) Ki * Kj;
    for (std::size_t idx = 0; idx < ncell; idx++) {
      if (tab[idx] == 0.0) { tab[idx] = correct; ntab += correct; }
    }
  }

  // Precompute the rho-INDEPENDENT part of the quadrature once per pair: for each X (row)
  // band, the clamped limits, the abscissae, and the weight*phi(x) factors depend only on the
  // fixed thresholds, so they are reused across every likelihood evaluation instead of being
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

  // Smallest positive double, used only to keep log()/division finite for an utterly
  // negligible (underflowing) cell; the conditioning integral is otherwise non-negative.
  const double p_floor = std::numeric_limits<double>::min();

  // Evaluate the negative log-likelihood, the score dl/drho, and the expected (Fisher)
  // information at rho. The cell probabilities P_ab and their derivatives dP_ab/drho are built
  // from the SAME conditioning integral, so the score is exactly the derivative of the
  // quadrature objective (not of the unreachable exact probability) and the iteration
  // converges cleanly. Differentiating P_ab = int phi(x)[Phi((c1-rho x)/s) - Phi((c0-rho x)/s)] dx
  // under the integral sign gives d/drho Phi((c-rho x)/s) = phi((c-rho x)/s) (rho c - x) / s^3
  // (Plackett, 1954); an infinite column cut contributes zero. Only the inner column terms
  // depend on rho and are shared across the cells of a row. The information is the multinomial
  // expectation ntab * sum_ab (dP_ab/drho)^2 / P_ab (Olsson, 1979): positive, and needing no
  // second derivative.
  auto eval = [&](double rho, double& nll, double& score, double& info) {
    double s = std::sqrt(1.0 - rho * rho);
    double s3 = s * s * s;
    for (int a = 0; a < Ki; a++) {
      std::size_t base = (std::size_t) a * n;
      std::size_t pbase = (std::size_t) a * Kj;
      for (int b = 0; b < Kj; b++) { pmat[pbase + b] = 0.0; dpmat[pbase + b] = 0.0; }
      for (int k = 0; k < n; k++) {
        double w = wpdf[base + k];
        if (w == 0.0) continue;  // band beyond the clamp contributes nothing
        poly_accum_node(xnode[base + k], w, rho, s, s3, ccut.data(), Kj,
                        &pmat[pbase], &dpmat[pbase], colcdf.data(), coldcdf.data());
      }
    }
    // Accumulate the negative log-likelihood, score, and information over the cells.
    nll = 0.0;
    score = 0.0;
    double iacc = 0.0;
    for (int a = 0; a < Ki; a++) {
      std::size_t pbase = (std::size_t) a * Kj;
      for (int b = 0; b < Kj; b++) {
        double praw = pmat[pbase + b];
        double p = praw < p_floor ? p_floor : praw;
        double dP = dpmat[pbase + b];
        // Expected (Fisher) information sum_ab (dP)^2/P. A near-impossible cell carries no
        // information (dP^2/P -> 0 as P -> 0), but at |rho| near 1 its probability can underflow
        // while its rho-derivative, divided by s^3 = (1 - rho^2)^{3/2}, stays O(1); the floored
        // ratio would then explode and, through the Newton step score/info, fake a vanishing step
        // that stops the iteration at the warm start far from the optimum. Skip the underflowed
        // cells - their true contribution is negligible.
        if (praw > p_floor) iacc += dP * dP / praw;
        double nab = tab[pbase + b];
        if (nab != 0.0) {
          nll -= nab * std::log(p);
          // The score must be the gradient of this nll. When the cell probability has
          // underflowed (p held at p_floor), -nab*log(p) is locally constant in rho, so its
          // gradient is zero; gating the score on the same condition as the information keeps
          // the Newton step score/info consistent. Adding the unfloored derivative here would
          // inject a spurious large term that pulls a strongly-correlated pair (one with an
          // observed near-impossible off-diagonal cell) away from its optimum.
          if (praw > p_floor) score += nab * dP / p;
        }
      }
    }
    info = ntab * iacc;
  };

  // Warm start at the Pearson correlation of the category codes (lavaan / Olsson, 1979): it
  // is already close to rho, so the scoring iteration starts inside the basin and converges in
  // a few steps instead of overshooting from rho = 0 as a cold start would. Computed from the
  // contingency table with the category indices as scores; falls back to 0 for a degenerate
  // (zero-variance) margin.
  double sa = 0.0, sb = 0.0, saa = 0.0, sbb = 0.0, sab = 0.0;
  for (int a = 0; a < Ki; a++) {
    std::size_t pbase = (std::size_t) a * Kj;
    for (int b = 0; b < Kj; b++) {
      double nij = tab[pbase + b];
      if (nij == 0.0) continue;
      sa += nij * a; sb += nij * b;
      saa += nij * (double) a * a; sbb += nij * (double) b * b;
      sab += nij * (double) a * b;
    }
  }
  double cov = sab - sa * sb / ntab;
  double va = saa - sa * sa / ntab, vb = sbb - sb * sb / ntab;
  double rho0 = (va > 0.0 && vb > 0.0) ? cov / std::sqrt(va * vb) : 0.0;

  // Maximise the log-likelihood over (-maxcor, maxcor), bounded away from the singular
  // endpoints as in polycor, by damped Fisher scoring (step = score / information).
  // Convergence is judged on the full scoring step (the Newton decrement) at an interior
  // point. A scoring step can still overshoot the optimum into the near-singular region close
  // to rho = +/-1, where the cell probabilities underflow and the information blows up; to
  // prevent that, each step is capped to at most half the distance to the boundary, which
  // shrinks to zero as rho -> +/-1 yet never binds near an interior optimum (the score, and so
  // the step, vanishes there on its own). The capped step is still backtracked by halving to
  // keep the negative log-likelihood monotone. Any pair that fails to converge or yields
  // unusable information is finished by Brent's method on the same objective.
  //
  // The step tolerance is ~1e-5 in rho: near the optimum the log-likelihood is flat, so a
  // tighter target would chase changes below floating-point noise (a tolerance also used by
  // comparable scoring implementations); 1e-5 is far inside the accuracy of the estimate. The
  // same constant is reused as Brent's relative-x tolerance in the fallback (where it bounds
  // |rho - rho*| to ~1e-5 * |rho|), which is likewise far inside the estimate's accuracy.
  const double maxcor = 0.9999;
  const double tol = 1e-5;
  double rho = rho0 > maxcor ? maxcor : (rho0 < -maxcor ? -maxcor : rho0);
  double nll, score, info;
  eval(rho, nll, score, info);
  bool converged = false;
  for (int iter = 0; iter < 50; iter++) {
    // The expected information excludes underflowed cells (see eval), so a high-|rho| empty
    // tail cell can no longer inflate it; this guard is defensive, handing any pair that still
    // yields non-positive or non-finite information/score to the robust Brent fallback.
    if (!(info > 0.0) || !std::isfinite(info) || !std::isfinite(score)) break;
    double step = score / info;
    if (std::abs(step) < tol && std::abs(rho) < maxcor) {  // interior Newton decrement -> optimum
      converged = true;
      break;
    }
    double cap = 0.5 * (1.0 - std::abs(rho));  // never jump more than halfway to +/-1
    if (step >  cap) step =  cap;
    else if (step < -cap) step = -cap;
    double t = 1.0;
    bool improved = false;
    double cand = rho, nn = nll, sc = score, inf = info;
    for (int h = 0; h < 40; h++) {
      double c = rho + t * step;
      if (c < -maxcor) c = -maxcor;
      else if (c > maxcor) c = maxcor;
      double f, g, fi;
      eval(c, f, g, fi);
      if (f <= nll) { cand = c; nn = f; sc = g; inf = fi; improved = true; break; }
      t *= 0.5;
    }
    if (!improved || cand == rho) break;  // line search stalled / pinned at boundary -> Brent
    rho = cand; nll = nn; score = sc; info = inf;
  }
  if (!converged) {
    rho = brent_min([&](double r) { double f, g, fi; eval(r, f, g, fi); return f; },
                    -maxcor, maxcor, tol, 200);
  }
  return rho;
}

// Cell probabilities P_ab and their rho-derivatives dP_ab/drho at a fixed rho, for the
// ACOV. Same conditioning integral and rho-derivative as the likelihood evaluation in
// polychoric_pair, but at a single rho with no optimisation. cut_i/cut_j are the threshold
// vectors already padded with +/-Inf; pmat/dpmat (Ki*Kj) and colcdf/coldcdf (>= Kj+1) are
// caller-owned scratch so the per-pair call does not allocate.
static void poly_cell_probs(double rho, const std::vector<double>& cut_i,
                            const std::vector<double>& cut_j, int Ki, int Kj,
                            const std::vector<double>& gx, const std::vector<double>& gw,
                            std::vector<double>& pmat, std::vector<double>& dpmat,
                            std::vector<double>& colcdf, std::vector<double>& coldcdf) {
  const int n = (int) gx.size();
  double s = std::sqrt(1.0 - rho * rho);
  double s3 = s * s * s;
  for (int a = 0; a < Ki; a++) {
    double lo = cut_i[a]     < -POLY_CLAMP ? -POLY_CLAMP : cut_i[a];
    double hi = cut_i[a + 1] >  POLY_CLAMP ?  POLY_CLAMP : cut_i[a + 1];
    std::size_t pbase = (std::size_t) a * Kj;
    for (int b = 0; b < Kj; b++) { pmat[pbase + b] = 0.0; dpmat[pbase + b] = 0.0; }
    if (hi <= lo) continue;
    double mid = 0.5 * (lo + hi), half = 0.5 * (hi - lo);
    for (int k = 0; k < n; k++) {
      double x = mid + half * gx[k];
      double w = gw[k] * half * std_norm_pdf(x);
      poly_accum_node(x, w, rho, s, s3, cut_j.data(), Kj,
                      &pmat[pbase], &dpmat[pbase], colcdf.data(), coldcdf.data());
    }
  }
}

// Cross-information A21[k] = sum_cases (score_rho * score_tau_k) between rho and the
// thresholds of one variable of a pair, for the ACOV's threshold correction. Only the two
// contingency rows/columns adjacent to threshold k contribute; dpi_ab/dtau_k = phi(tau_k) *
// (conditional band over the OTHER variable) (Joreskog, 1994). The two variables differ only
// in the cell stride: with the threshold variable on the rows, cell (k,o) sits at index
// k*Kother + o and its lower-neighbour at +Kother; on the columns, cell (o,k) sits at
// o*Kself + k and its neighbour at +1. Passing (stride_self, stride_other) handles both, so
// the formula is written once. cut_self/cut_other are padded with +/-Inf.
static void poly_cross_info(int n_self_thr, int n_other, int stride_self, int stride_other,
                            const std::vector<double>& cut_self,
                            const std::vector<double>& phi_self,
                            const std::vector<double>& cut_other, double rho, double s,
                            const double* nab, const double* dxr, const double* pmat,
                            double p_floor, std::vector<double>& colb, double* A21out) {
  for (int k = 0; k < n_self_thr; k++) {
    double x = cut_self[k + 1];                       // tau_self[k]
    for (int o = 0; o <= n_other; o++) colb[o] = std_norm_cdf((cut_other[o] - rho * x) / s);
    double acc = 0.0;
    for (int o = 0; o < n_other; o++) {
      double band = colb[o + 1] - colb[o];
      std::size_t up = (std::size_t) k * stride_self + (std::size_t) o * stride_other;
      double pu = pmat[up] < p_floor ? p_floor : pmat[up];
      double term = nab[up] * dxr[up] / pu;
      std::size_t lo = up + stride_self;               // the adjacent row/column of band k
      double pl = pmat[lo] < p_floor ? p_floor : pmat[lo];
      term -= nab[lo] * dxr[lo] / pl;
      acc += band * term;
    }
    A21out[k] = phi_self[k] * acc;
  }
}

// Per-variable threshold "bread" for the polychoric ACOV. From the listwise category
// counts m (length K) and thresholds tau (length K-1), returns the per-category threshold
// influence IFth (K x (K-1)), IFth(a,.) = A11^{-1} SC.TH(a), where SC.TH(a) is the marginal
// (univariate) ML threshold score for a case in category a and A11 = sum_a m_a SC.TH(a)
// SC.TH(a)' is its outer-product (Fisher) information. Mirrors lavaan's lav_uvord_scores /
// the A11 block of muthen1984. An empty category (m_a = 0) contributes a zero score row.
static arma::mat poly_threshold_bread(const std::vector<double>& m, int K,
                                      const std::vector<double>& tau, double Nc) {
  int nth = K - 1;
  arma::mat S(K, nth, arma::fill::zeros);          // SC.TH(category a, threshold k)
  arma::vec mv(K);
  for (int a = 0; a < K; a++) {
    mv[a] = m[a];
    double inv_pm = (m[a] > 0.0) ? Nc / m[a] : 0.0;  // 1 / P(category a); 0 if empty
    if (a <= nth - 1) S(a, a)     +=  std_norm_pdf(tau[a]) * inv_pm;       // upper boundary
    if (a >= 1)       S(a, a - 1) += -std_norm_pdf(tau[a - 1]) * inv_pm;   // lower boundary
  }
  arma::mat A11 = S.t() * arma::diagmat(mv) * S;
  arma::mat A11inv;
  if (!arma::inv_sympd(A11inv, A11) && !arma::inv(A11inv, A11) &&
      !arma::pinv(A11inv, A11)) {
    A11inv = arma::zeros<arma::mat>(nth, nth);
  }
  return S * A11inv;                                // IFth (K x nth)
}

// Internal C++ backend for the two-step polychoric/tetrachoric matrix, called from
// .polychoric(), which owns the user-facing validation and classed conditions.
// [[Rcpp::export(.polychoric_cpp)]]
Rcpp::List polychoric_cpp(Rcpp::IntegerMatrix x, std::string acov,
                          double correct, bool nearest_pd, int n_threads) {
  if (acov != "none" && acov != "diag" && acov != "full") {
    Rcpp::stop("`acov` must be one of \"none\", \"diag\", or \"full\".");
  }
  const int nrow = x.nrow();
  const int p = x.ncol();
  if (p < 2) {
    Rcpp::stop("At least two variables are required for a correlation matrix.");
  }
  if (n_threads < 1) n_threads = 1;

  const int* xp = x.begin();

  // Gauss-Legendre rules for the rectangle quadrature, computed once on the main thread: the
  // base 12-node rule and the finer rule used to re-estimate near-collinear pairs (|rho| close
  // to 1), where the base rule under-resolves the narrow conditional transition.
  std::vector<double> gx, gw, gx_hi, gw_hi;
  gauss_legendre(POLY_GL_N, gx, gw);
  gauss_legendre(POLY_GL_N_HI, gx_hi, gw_hi);

  // Whether any cell is missing. With complete data every pair shares the per-variable
  // thresholds (the fast, exact path); with missing data each pair instead uses its own
  // pairwise-complete thresholds (computed inside the parallel loop via a thread-safe quantile).
  bool any_missing = false;
  for (std::size_t k = 0; k < (std::size_t) nrow * p; k++) {
    if (xp[k] == NA_INTEGER) { any_missing = true; break; }
  }

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
  // Records, per pair, whether the finer (POLY_GL_N_HI) rule produced the estimate, so the
  // ACOV below can reuse the same rule rather than re-deriving it from the rounded estimate
  // (a refined value can land just the other side of POLY_REFINE_RHO).
  std::vector<char> used_hi((std::size_t) npairs, 0);

  #pragma omp parallel num_threads(n_threads)
  {
    // Thread-local scratch reused across this thread's pairs (no allocation in the loop).
    std::vector<double> tab((std::size_t) Kmax * Kmax);
    std::vector<double> pmat((std::size_t) Kmax * Kmax);
    std::vector<double> dpmat((std::size_t) Kmax * Kmax);
    std::vector<double> colcdf(Kmax + 1);
    std::vector<double> coldcdf(Kmax + 1);
    std::vector<double> rcut(Kmax + 1), ccut(Kmax + 1);
    std::vector<double> xnode((std::size_t) Kmax * POLY_GL_N_HI);
    std::vector<double> wpdf((std::size_t) Kmax * POLY_GL_N_HI);

    #pragma omp for schedule(dynamic)
    for (int t = 0; t < npairs; t++) {
      int i = pair_i[t];
      int j = pair_j[t];
      double r = polychoric_pair(xp, nrow, i, j, tau[i], tau[j], Kcat[i], Kcat[j], correct,
                                 any_missing, gx, gw, tab, pmat, dpmat, colcdf, coldcdf,
                                 rcut, ccut, xnode, wpdf);
      // Near |rho| = 1 the 12-node rule under-resolves the conditional transition and biases
      // the estimate, so re-estimate the rare near-collinear pair with the finer rule.
      if (std::isfinite(r) && std::fabs(r) > POLY_REFINE_RHO) {
        r = polychoric_pair(xp, nrow, i, j, tau[i], tau[j], Kcat[i], Kcat[j], correct,
                            any_missing, gx_hi, gw_hi, tab, pmat, dpmat, colcdf, coldcdf,
                            rcut, ccut, xnode, wpdf);
        used_hi[t] = 1;
      }
      rho[t] = r;
    }
  }

  arma::mat Rmat(p, p, arma::fill::eye);
  for (int t = 0; t < npairs; t++) {
    double v = rho[t];
    Rmat(pair_i[t], pair_j[t]) = v;
    Rmat(pair_j[t], pair_i[t]) = v;
  }

  // Optional nearest-PD projection, gated by `nearest_pd` so callers can choose how it
  // composes with downstream smoothing. When requested it first takes the cheaper values-only
  // eigendecomposition to test definiteness and computes the eigenvectors only on the rare
  // projection branch; callers that do not ask for it (the default, including every bootstrap
  // replicate) skip the eigendecomposition entirely. A non-finite matrix (a pair with no
  // overlapping complete cases) is left untouched here; the R wrapper raises a classed
  // condition on the NA.
  bool pd_adjusted = false;
  if (nearest_pd && Rmat.is_finite()) {
    arma::vec eigval;
    if (!arma::eig_sym(eigval, Rmat)) {
      Rcpp::stop("Eigendecomposition of the polychoric correlation matrix failed.");
    }
    if (eigval.min() < std::numeric_limits<double>::epsilon()) {
      // Clip the eigenvalues to a small positive floor and rescale to a unit diagonal (the
      // spectral nearest-correlation projection of Rebonato & Jackel, 2000; an eigenvalue
      // clip-and-rescale, not Higham's iterative Frobenius-nearest algorithm).
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

  // Asymptotic covariance of the polychoric correlations (gated by `acov`). The two-step
  // estimator's ACOV must account for the estimated thresholds (Muthen, 1984; Joreskog,
  // 1994): each rho's influence function is its conditional score minus the threshold
  // influence carried through the implicit derivative drho/dtau. We follow lavaan's
  // empirical (outer-product) two-step construction (muthen1984): per case build the
  // correlation influence IF = (s_rho - A21 . A11^{-1} . s_tau) / A22 from the per-case
  // scores, then Gamma = crossprod(IF) -- the variance-scale covariance Var(rho-hat)
  // (lavaan's WLS.W; diag(WLS.W)^{-1} are the DWLS weights). Because every per-case score is
  // constant within a contingency cell, the diagonal reduces to a cell sum (the cheap path)
  // and equals diag(Gamma) by construction. When an ACOV is requested the R wrapper restricts
  // the input to the listwise-complete rows first, so the point estimate, thresholds, and
  // ACOV share one case set (a sandwich covariance must be that of the estimator that
  // produced the estimates); here every row is complete and the point-estimate thresholds
  // `tau` are reused directly. References:
  //   Muthen, B. (1984). A general structural equation model with dichotomous, ordered
  //     categorical, and continuous latent variable indicators. Psychometrika, 49, 115-132.
  //   Joreskog, K. G. (1994). On the estimation of polychoric correlations and their
  //     asymptotic covariance matrix. Psychometrika, 59, 381-389.
  Rcpp::RObject acov_out = R_NilValue;
  if (acov == "diag" || acov == "full") {
    const double p_floor = std::numeric_limits<double>::min();
    const int Nc = nrow;                  // listwise-complete (enforced by the R wrapper)
    const double Ncd = (double) Nc;

    // Per variable: padded threshold cuts and densities (reusing the point-estimate
    // thresholds, finite and exact on complete data), the marginal category counts, and the
    // per-category threshold influence (the A11^{-1} bread).
    std::vector< std::vector<double> > cut_v(p), phi_v(p);
    std::vector<arma::mat> IFth(p);
    for (int j = 0; j < p; j++) {
      int K = Kcat[j];
      const int* col = xp + (std::size_t) j * nrow;
      std::vector<double> mcnt(K, 0.0);
      for (int r = 0; r < nrow; r++) {
        if (col[r] != NA_INTEGER) mcnt[col[r]] += 1.0;  // complete here; guard misuse
      }
      std::vector<double> ph(K - 1), cut(K + 1);
      for (int k = 0; k < K - 1; k++) ph[k] = std_norm_pdf(tau[j][k]);
      cut[0] = -POLY_INF; cut[K] = POLY_INF;
      for (int k = 0; k < K - 1; k++) cut[k + 1] = tau[j][k];
      phi_v[j] = ph;
      cut_v[j] = cut;
      IFth[j] = poly_threshold_bread(mcnt, K, tau[j], Ncd);
    }

    std::vector<double> acov_diag((std::size_t) npairs, 0.0);
    arma::mat IF_cor;
    if (acov == "full") IF_cor.set_size(Nc, npairs);

    #pragma omp parallel num_threads(n_threads)
    {
      // Thread-local scratch reused across this thread's pairs.
      std::vector<double> pmat((std::size_t) Kmax * Kmax);
      std::vector<double> dpmat((std::size_t) Kmax * Kmax);
      std::vector<double> nab((std::size_t) Kmax * Kmax);
      std::vector<double> dxr((std::size_t) Kmax * Kmax);
      std::vector<double> cif((std::size_t) Kmax * Kmax);
      std::vector<double> colcdf(Kmax + 1), coldcdf(Kmax + 1), colb(Kmax + 1);
      std::vector<double> Ti(Kmax), Tj(Kmax), A21i(Kmax), A21j(Kmax);

      #pragma omp for schedule(dynamic)
      for (int t = 0; t < npairs; t++) {
        int i = pair_i[t], j = pair_j[t];
        int Ki = Kcat[i], Kj = Kcat[j];
        double r = rho[t];
        double s = std::sqrt(1.0 - r * r);
        const std::vector<double>& ci = cut_v[i];
        const std::vector<double>& cj = cut_v[j];

        // Use the same quadrature rule that estimated this pair (the finer rule for the rare
        // near-collinear pairs), so the cell probabilities are consistent with rho.
        bool hi = used_hi[t] != 0;
        poly_cell_probs(r, ci, cj, Ki, Kj, hi ? gx_hi : gx, hi ? gw_hi : gw,
                        pmat, dpmat, colcdf, coldcdf);

        // Cell counts (every row is listwise-complete; the NA guard only matters under direct
        // misuse of the C++ entry point).
        std::fill(nab.begin(), nab.begin() + (std::size_t) Ki * Kj, 0.0);
        const int* coli = xp + (std::size_t) i * nrow;
        const int* colj = xp + (std::size_t) j * nrow;
        for (int r0 = 0; r0 < nrow; r0++) {
          if (coli[r0] == NA_INTEGER || colj[r0] == NA_INTEGER) continue;
          nab[(std::size_t) coli[r0] * Kj + colj[r0]] += 1.0;
        }

        // dx.rho per cell and the outer-product rho-information A22 = sum_cases s_rho^2.
        double A22 = 0.0;
        for (int a = 0; a < Ki; a++) {
          for (int b = 0; b < Kj; b++) {
            std::size_t idx = (std::size_t) a * Kj + b;
            double pf = pmat[idx] < p_floor ? p_floor : pmat[idx];
            double d = dpmat[idx] / pf;
            dxr[idx] = d;
            A22 += nab[idx] * d * d;
          }
        }

        // Cross-information of rho with each variable's thresholds: variable i on the rows
        // (cell stride Kj, 1), variable j on the columns (cell stride 1, Kj).
        poly_cross_info(Ki - 1, Kj, Kj, 1, ci, phi_v[i], cj, r, s,
                        nab.data(), dxr.data(), pmat.data(), p_floor, colb, A21i.data());
        poly_cross_info(Kj - 1, Ki, 1, Kj, cj, phi_v[j], ci, r, s,
                        nab.data(), dxr.data(), pmat.data(), p_floor, colb, A21j.data());

        // Threshold-influence aggregates T_i(a) = A21_i . IFth_i(a), T_j(b) likewise.
        const arma::mat& Fi = IFth[i];
        const arma::mat& Fj = IFth[j];
        for (int a = 0; a < Ki; a++) {
          double tt = 0.0;
          for (int k = 0; k < Ki - 1; k++) tt += A21i[k] * Fi(a, k);
          Ti[a] = tt;
        }
        for (int b = 0; b < Kj; b++) {
          double tt = 0.0;
          for (int k = 0; k < Kj - 1; k++) tt += A21j[k] * Fj(b, k);
          Tj[b] = tt;
        }

        // Per-cell correlation influence and the diagonal (= cell sum of n_ab * IF^2).
        double invA22 = (A22 > 0.0) ? 1.0 / A22 : 0.0;
        double vsum = 0.0;
        for (int a = 0; a < Ki; a++) {
          for (int b = 0; b < Kj; b++) {
            std::size_t idx = (std::size_t) a * Kj + b;
            double v = (dxr[idx] - Ti[a] - Tj[b]) * invA22;
            cif[idx] = v;
            vsum += nab[idx] * v * v;
          }
        }
        acov_diag[t] = vsum;

        // Full Gamma: scatter the per-case influence into this pair's column (disjoint).
        if (acov == "full") {
          double* outcol = IF_cor.colptr(t);
          for (int r0 = 0; r0 < nrow; r0++) {
            outcol[r0] = (coli[r0] == NA_INTEGER || colj[r0] == NA_INTEGER)
                           ? 0.0 : cif[(std::size_t) coli[r0] * Kj + colj[r0]];
          }
        }
      }
    }

    if (acov == "diag") {
      acov_out = Rcpp::NumericVector(acov_diag.begin(), acov_diag.end());
    } else {
      arma::mat Gamma = IF_cor.t() * IF_cor;           // pstar x pstar, variance scale
      Gamma = 0.5 * (Gamma + Gamma.t());               // symmetrise away round-off
      acov_out = Rcpp::wrap(Gamma);
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
    Rcpp::Named("acov") = acov_out);
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

  std::vector<double> gx, gw, gx_hi, gw_hi;
  gauss_legendre(POLY_GL_N, gx, gw);
  gauss_legendre(POLY_GL_N_HI, gx_hi, gw_hi);
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; i++) {
    double rr = rho[i % nr];
    bool hi = std::fabs(rr) > POLY_REFINE_RHO;
    out[i] = bvn_rect(a0[i % na0], a1[i % na1], b0[i % nb0], b1[i % nb1],
                      rr, hi ? gx_hi : gx, hi ? gw_hi : gw);
  }
  return out;
}
