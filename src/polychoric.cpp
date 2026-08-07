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
// sum of non-negative terms, so the OUTER integration cannot cancel. The inner conditional
// band is still a difference of two normal CDFs, and in the far tails both of them round to
// the same value, so it is taken from whichever tail is nearer to zero (see
// std_norm_cdf_signed): the band, and hence the cell probability, then keeps full relative
// accuracy however deep in the tail it lies. It is evaluated by Gauss-Legendre quadrature,
// infinite X-cuts clamped to a finite range where phi is negligible, and infinite Y-cuts
// handled directly by the normal CDF (Phi(+/-Inf) = 1/0).
//
// The per-pair correlation maximises the log-likelihood by Fisher scoring: the score is the
// closed-form derivative of each cell probability obtained by
// differentiating the same conditioning integral with respect to rho (so it is consistent with
// the quadrature used for the probabilities), and the step is normalised by the expected
// information of the cell-count multinomial (Olsson, 1979). The step is capped away from the
// singular endpoints and backtracked to keep the likelihood monotone, with Brent's (1973)
// method as a fallback. Every per-pair quantity is pure C++/Armadillo (std::erfc for the normal
// CDF), so the pair loop is allocation-light and reuses caller-owned scratch.
//
// References:
//   Olsson, U. (1979). Maximum likelihood estimation of the polychoric correlation
//     coefficient. Psychometrika, 44, 443-460.
//   Plackett, R. L. (1954). A reduction formula for normal multivariate integrals.
//     Biometrika, 41, 351-360.
//   Genz, A., & Bretz, F. (2009). Computation of Multivariate Normal and t Probabilities.
//     Springer.
//   Higham, N. J. (2002). Accuracy and Stability of Numerical Algorithms (2nd ed.). SIAM.
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
// A pair with an empty contingency cell is re-estimated with the finer rule only once its cheap
// (12-node) estimate reaches this gate. A table with an empty cell biases the 12-node estimate
// only when it is strongly correlated (a sharp conditional transition over a wide tail band); at
// a low or moderate |rho| the conditional is smooth and the 12-node estimate is already accurate,
// so refining would be wasted work. The gate sits well below the worst under-shoot of a biased
// pair (empirically about 0.75, since the 12-node estimate of such a pair lands low), leaving
// margin while keeping the finer rule off the many weakly-correlated sparse pairs of heterogeneous
// ordinal data.
static const double POLY_EMPTY_REFINE_RHO = 0.65;
// Largest |rho| the estimator admits. The log-likelihood is maximised over the closed interval
// [-POLY_MAXCOR, POLY_MAXCOR], bounded away from the singular endpoints as in polycor (whose
// `maxcor` default is the same 0.9999), so the model never degenerates to a rank-one normal.
// It is therefore also the estimate reported for a pair whose table is at a Frechet bound, where
// the unconstrained maximiser is exactly +/-1 (see poly_frechet_bound).
static const double POLY_MAXCOR = 0.9999;
// Continuity correction added to the empty cell of a binary (2x2) pair before estimation. Every
// 2x2 table with an empty cell sits on a Frechet bound of its own marginals, so its uncorrected
// two-step estimate is the boundary +/-1 whatever the underlying correlation; nudging the empty
// cell moves the estimate back into the interior. Savalei (2011) recommends adding 0.5 for two
// categories (and nothing beyond two), which is also lavaan's default (`zero.add = c(0.5, 0)`) and
// psych's (`correct = 0.5`). Choi, Kim, Chen, & Zhang (2025) compare 0.1 / 0.25 / 0.5 against four
// ways of adding them and find the added value dominates while the manner is second-order, so the
// manner below is chosen for consistency with the shared thresholds rather than for accuracy.
// Changing this value also means changing the two literals in the `efa_cor_zero_cell` warning in
// .polychoric(), which state it to the user, and the documented convention in efa_fit()'s
// *Correlation methods* section.
static const double POLY_ZERO_ADD = 0.5;

static inline double std_norm_pdf(double z) {
  return POLY_INV_SQRT_2PI * std::exp(-0.5 * z * z);
}

// Signed-tail normal CDF: Phi(z) for z < 0, and -(1 - Phi(z)) for z >= 0. Each branch is one
// std::erfc call on a non-negative argument, so it costs the same as a plain Phi while always
// returning the tail NEARER to zero, computed to full relative accuracy (erfc keeps it; the
// same quantity formed as 1 - Phi(z) by subtraction would not). It is the form in which a
// probability band Phi(z1) - Phi(z0) can be differenced without cancellation -- see
// poly_signed_band, which does the differencing.
//
// Why relative accuracy is the requirement here, and not merely absolute: a band deep in a tail
// is a cell probability of the ordinal model, and the log-likelihood, its score dP/P and the
// asymptotic covariance are all assembled from that probability, so an error relative to its own
// magnitude propagates undiminished however small the cell is (Genz & Bretz, 2009, on tail
// accuracy for these integrals). Forming a small difference from two nearly equal quantities is
// the textbook way to lose exactly that (Higham, 2002).
//
// Concretely: the upper-tail band of a strongly correlated pair is Phi(13.7) - Phi(10.4), and
// both of those round to exactly 1.0, so the band -- and hence the cell probability -- comes out
// 0 although its true value (~1e-25 here, ~1e-45 for the extreme corner cells) is hundreds of
// orders of magnitude above the underflow floor. The cell's rho-derivative is built from
// DENSITIES rather than a CDF difference and keeps its accuracy, so the score dP/P would be
// inflated by those same hundreds of orders of magnitude.
static inline double std_norm_cdf_signed(double z) {
  return z < 0.0 ? 0.5 * std::erfc(-z * POLY_INV_SQRT_2)
                 : -0.5 * std::erfc(z * POLY_INV_SQRT_2);
}

// One probability band from the signed-tail CDF of its two conditional cuts, q0 = q(z0) and
// q1 = q(z1) (std_norm_cdf_signed). When both cuts lie on the same side of zero the band is just
// q1 - q0 -- Phi(z1) - Phi(z0) in the lower tail, Phic(z0) - Phic(z1) in the upper, neither of
// which cancels. When they straddle zero it is that difference PLUS ONE, since
//   Phi(z1) - Phi(z0) = 1 - Phic(z1) - Phi(z0),
// which subtracts from 1 two quantities that are each at most 0.5 and so loses at most one bit.
// No straddling band lies in a tail, which is where the accuracy matters.
//
// The straddling band identifies itself from the SIGN BIT of q, with no cut index to track: q is
// positive (or +0.0, when the lower tail underflows) exactly when z < 0, and negative (or -0.0,
// since IEEE gives -0.5 * 0.0 = -0.0 when the upper tail underflows or z = +Inf) exactly when
// z >= 0. Callers pass consecutive cuts of a sorted, +/-Inf-padded vector, so z is monotone and
// at most one band per row can straddle; keeping the test local to the band rather than shared
// across the row is what lets all three band consumers use one definition.
static inline double poly_signed_band(double q0, double q1) {
  double band = q1 - q0;
  if (std::signbit(q1) && !std::signbit(q0)) band += 1.0;
  return band;
}

// Standard normal quantile via Wichura's (1988) AS 241 (accurate to ~1e-16), allocation-free.
// Used for the per-pair (pairwise-complete) thresholds under missing data; returns +/-Inf at
// p = 1/0 (an empty marginal category, handled by the +/-POLY_CLAMP clamping downstream).
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
// Computed once per call, used read-only inside the pair loop.
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
// The conditional Y band comes from poly_signed_band, so it is cancellation-free in either
// tail (infinite Y-cuts included, q(-Inf) = +0 and q(+Inf) = -0); infinite X-cuts are clamped
// to +/-POLY_CLAMP. Non-negative by construction.
// NOT the estimator's integrator: the likelihood accumulates cell probabilities through
// poly_accum_node() below. This is reached only from .bvn_rect_cpp(), the entry point the
// test suite uses to cross-check that integral against an external bivariate-normal rule.
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
    double rx = rho * x;
    double inner = poly_signed_band(std_norm_cdf_signed((b0 - rx) / s),
                                    std_norm_cdf_signed((b1 - rx) / s));
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
//
// colcdf holds the SIGNED-TAIL CDF (std_norm_cdf_signed), which poly_signed_band differences
// without cancellation in either tail. The derivative row needs no such care: coldcdf holds
// densities, which decay to zero rather than to a common limit, so their difference does not
// cancel.
static inline void poly_accum_node(double x, double w, double rho, double s, double s3,
                                   const double* cut_j, int Kj, double* prow, double* dprow,
                                   double* colcdf, double* coldcdf) {
  double rx = rho * x;
  for (int b = 0; b <= Kj; b++) {
    double c = cut_j[b];
    double z = (c - rx) / s;
    colcdf[b] = std_norm_cdf_signed(z);
    coldcdf[b] = std::isfinite(c) ? std_norm_pdf(z) * (rho * c - x) / s3 : 0.0;
  }
  for (int b = 0; b < Kj; b++) {
    prow[b]  += w * poly_signed_band(colcdf[b], colcdf[b + 1]);
    dprow[b] += w * (coldcdf[b + 1] - coldcdf[b]);
  }
}

// Fix the interior thresholds of a contingency table from its cumulative marginals (divisor
// `denom`, the table mass): rcut[a + 1] = Phi^{-1}(cumulative row mass / denom) and ccut[b + 1]
// likewise for columns. The outer cuts (rcut[0]/rcut[Ki], ccut[0]/ccut[Kj]) are +/-Inf and set by
// the caller. Used by the pairwise-marginal (missing-data) pass.
static inline void poly_fix_thresholds(const std::vector<double>& tab, int Ki, int Kj,
                                       double denom, std::vector<double>& rcut,
                                       std::vector<double>& ccut) {
  double cum = 0.0;
  for (int a = 0; a < Ki - 1; a++) {
    double rs = 0.0;
    for (int b = 0; b < Kj; b++) rs += tab[(std::size_t) a * Kj + b];
    cum += rs;
    rcut[a + 1] = std_norm_quantile(cum / denom);
  }
  cum = 0.0;
  for (int b = 0; b < Kj - 1; b++) {
    double cs = 0.0;
    for (int a = 0; a < Ki; a++) cs += tab[(std::size_t) a * Kj + b];
    cum += cs;
    ccut[b + 1] = std_norm_quantile(cum / denom);
  }
}

// Whether a contingency table is exactly one of the two Frechet bounds of its own marginals:
// returns +1 for the upper (comonotone) bound, -1 for the lower (antimonotone) bound, 0 otherwise.
//
// The bounds are the extreme couplings of fixed marginals (Frechet, 1951): on the cumulative
// scale the upper bound is M(u, v) = min(u, v) and the lower W(u, v) = max(u + v - 1, 0). These
// are exactly the bivariate normal copulas at rho = +1 and rho = -1, and the two-step thresholds
// fix the model's marginals to this table's marginals, so at rho = +/-1 the model cell
// probabilities ARE the corresponding bound coupling of the observed marginals. A table equal to
// that coupling is therefore reproduced exactly at the boundary: the log-likelihood attains its
// saturated value there, the maximiser is exactly rho = +/-1, and the rho-information vanishes
// (the score is driven by probability flux into the structurally empty cells, which decays like
// exp(-c / (1 - rho^2))). Both the estimate and its asymptotic variance are then boundary
// quantities, which is why the caller reports +/-POLY_MAXCOR and refuses to compute a variance.
//
// Every quantity here is a count: the cumulative marginals and the differenced bound couplings
// are exact integers held in doubles (the table mass is at most the row count, far below 2^53),
// so the comparisons are exact and the verdict is bit-identical on every platform - no floating
// point enters the decision. Cost is O(Ki * Kj) per pair.
//
// A table with fewer than two populated rows or columns is reported as 0: there the two bounds
// coincide (both reduce to the table itself) so the sign would be arbitrary, and such a pair
// carries no information about rho in any case.
//
// Rc/Cc are caller-owned scratch of at least Ki + 1 / Kj + 1 elements, so the pair loop does not
// allocate; they are overwritten with the cumulative marginals on every call.
static int poly_frechet_bound(const std::vector<double>& tab, int Ki, int Kj,
                              std::vector<double>& Rc, std::vector<double>& Cc) {
  // Cumulative marginals R_a = sum_{a' <= a} row_a' and C_b likewise, with R_{-1} = C_{-1} = 0.
  // A marginal with fewer than two populated categories leaves the sign of the bound arbitrary,
  // so it is counted while the sums are formed rather than in a second pass.
  int nz_row = 0, nz_col = 0;
  Rc[0] = 0.0;
  for (int a = 0; a < Ki; a++) {
    double rs = 0.0;
    for (int b = 0; b < Kj; b++) rs += tab[(std::size_t) a * Kj + b];
    if (rs > 0.0) nz_row++;
    Rc[(std::size_t) a + 1] = Rc[(std::size_t) a] + rs;
  }
  Cc[0] = 0.0;
  for (int b = 0; b < Kj; b++) {
    double cs = 0.0;
    for (int a = 0; a < Ki; a++) cs += tab[(std::size_t) a * Kj + b];
    if (cs > 0.0) nz_col++;
    Cc[(std::size_t) b + 1] = Cc[(std::size_t) b] + cs;
  }
  if (nz_row < 2 || nz_col < 2) return 0;
  const double N = Rc[(std::size_t) Ki];

  // Recover each bound's cell counts from its cumulative table by inclusion-exclusion and compare.
  bool upper = true, lower = true;
  for (int a = 0; a < Ki && (upper || lower); a++) {
    for (int b = 0; b < Kj && (upper || lower); b++) {
      double obs = tab[(std::size_t) a * Kj + b];
      if (upper) {
        double cell = std::min(Rc[(std::size_t) a + 1], Cc[(std::size_t) b + 1])
                    - std::min(Rc[(std::size_t) a + 1], Cc[(std::size_t) b])
                    - std::min(Rc[(std::size_t) a],     Cc[(std::size_t) b + 1])
                    + std::min(Rc[(std::size_t) a],     Cc[(std::size_t) b]);
        if (cell != obs) upper = false;
      }
      if (lower) {
        double cell = std::max(Rc[(std::size_t) a + 1] + Cc[(std::size_t) b + 1] - N, 0.0)
                    - std::max(Rc[(std::size_t) a + 1] + Cc[(std::size_t) b]     - N, 0.0)
                    - std::max(Rc[(std::size_t) a]     + Cc[(std::size_t) b + 1] - N, 0.0)
                    + std::max(Rc[(std::size_t) a]     + Cc[(std::size_t) b]     - N, 0.0);
        if (cell != obs) lower = false;
      }
    }
  }
  // A table cannot be both once each marginal has two populated categories (guarded above).
  if (upper) return 1;
  if (lower) return -1;
  return 0;
}

// Continuity-correct a binary pair whose table has a single empty cell, in place; returns whether
// it applied. Only 2x2 tables with exactly one empty cell qualify - Savalei (2011) recommends the
// correction for two categories and against it beyond two, and with two or more empty cells the
// table carries too little information for a nudge to mean anything.
//
// The correction is marginal-preserving: POLY_ZERO_ADD goes into the empty cell AND into the cell
// diagonally opposite it, and comes back out of the two cells that share neither its row nor its
// column. Every row and column sum - and hence the table total - is left exactly as it was. That
// matters structurally here rather than only statistically: the thresholds are computed once per
// VARIABLE from its full-column marginals and shared by every pair it appears in, so a correction
// that moved this pair's margins would put the pair's own table out of step with the thresholds it
// is scored against, and with every other pair built on the same variable. Preserving the margins
// keeps the shared thresholds exactly valid, so nothing downstream needs recomputing. It is also
// what lavaan does (`zero.keep.margins`), and Choi et al. (2025) find the manner of adding matters
// far less than the value, so consistency is the right thing to optimise for.
//
// No feasibility check is needed: the table holds integer counts, so a 2x2 with exactly one empty
// cell has all three other cells >= 1 > POLY_ZERO_ADD, and none can be driven negative.
static bool poly_zero_correct_2x2(std::vector<double>& tab, int Ki, int Kj) {
  if (Ki != 2 || Kj != 2) return false;
  int zeros = 0, z = -1;
  for (int k = 0; k < 4; k++) {
    if (tab[(std::size_t) k] == 0.0) { zeros++; z = k; }
  }
  if (zeros != 1) return false;
  // Row-major c(a, b, c, d): the diagonal opposite of index z is 3 - z (a<->d, b<->c), and the
  // remaining two indices are the ones that lose the same amount.
  const int opp = 3 - z;
  for (int k = 0; k < 4; k++) {
    if (k == z || k == opp) tab[(std::size_t) k] += POLY_ZERO_ADD;
    else                    tab[(std::size_t) k] -= POLY_ZERO_ADD;
  }
  return true;
}

// Two-step polychoric correlation for one variable pair. Builds the contingency table from the
// pairwise-complete rows, fixes the thresholds (the shared per-variable tau_i/tau_j, or - when
// use_local is set, i.e. the data have missing values - this pair's own marginal thresholds),
// then minimises the negative log-likelihood over rho. Pure C++; reuses caller-owned scratch
// buffers so the hot path does not allocate. Returns NA_REAL when the pair has no overlapping
// complete cases (the R wrapper turns a resulting NA into a classed condition).
static double polychoric_pair(const int* xp, int nrow, int ci, int cj,
                              const std::vector<double>& tau_i,
                              const std::vector<double>& tau_j,
                              int Ki, int Kj, bool use_local,
                              const std::vector<double>& gx, const std::vector<double>& gw,
                              std::vector<double>& tab, std::vector<double>& pmat,
                              std::vector<double>& dpmat, std::vector<double>& colcdf,
                              std::vector<double>& coldcdf, std::vector<double>& rcut,
                              std::vector<double>& ccut, std::vector<double>& xnode,
                              std::vector<double>& wpdf, bool& had_empty, int& at_bound,
                              bool& zero_corrected) {
  const int n = (int) gx.size();
  at_bound = 0;
  zero_corrected = false;
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

  // Whether any contingency cell is empty, reported back through `had_empty`. A handful of
  // extreme cells then dominate the likelihood, and the fixed 12-node rule under-resolves the
  // wide (clamped) tail bands of such a table - mislocating the maximum, sometimes badly, even
  // at moderate rho - so the caller re-estimates these pairs with the finer quadrature rule,
  // exactly as it does for the near-collinear pairs. Dense tables (no empty cell) keep the
  // cheap 12-node rule.
  had_empty = false;
  for (std::size_t idx = 0, ncell = (std::size_t) Ki * Kj; idx < ncell; idx++) {
    if (tab[idx] == 0.0) { had_empty = true; break; }
  }

  // A table at one of the Frechet bounds of its own marginals is reproduced exactly at
  // rho = +/-1, so the maximiser over the estimator's domain is that domain's endpoint (see
  // poly_frechet_bound). Return it directly: there is nothing for the iteration to find, the
  // likelihood is numerically flat over the whole approach to the boundary - so an optimiser
  // stops at an arbitrary point there, which is what made the estimate platform-dependent - and
  // the quadrature would in any case under-resolve a conditional transition of width
  // sqrt(1 - POLY_MAXCOR^2) = 0.014. The verdict comes from integer cell counts alone, so this
  // branch is taken identically on every platform.
  //
  // colcdf/coldcdf serve as the detector's cumulative-marginal scratch: they are sized Kmax + 1,
  // and nothing has written or read them at this point -- the quadrature below fills them from
  // scratch on every node -- so borrowing them here keeps the pair loop allocation-free.
  at_bound = poly_frechet_bound(tab, Ki, Kj, colcdf, coldcdf);
  if (at_bound != 0) {
    // A binary pair is the one case where the bound is repaired rather than reported: every 2x2
    // table with an empty cell is at a bound, so without a correction no binary pair with a rare
    // response combination could ever be estimated in the interior. Nudge the empty cell (see
    // poly_zero_correct_2x2) and carry on into the ordinary iteration below; the correction
    // preserves the marginals, so the thresholds padded below are still this table's own. Larger
    // tables, and 2x2 tables with more than one empty cell, keep the boundary report.
    if (poly_zero_correct_2x2(tab, Ki, Kj)) {
      zero_corrected = true;
      at_bound = 0;
      // `had_empty` deliberately stays true: the corrected cell holds POLY_ZERO_ADD rather than a
      // count, so the table is still effectively empty there and the finer quadrature rule is
      // still needed once the estimate is at least moderately correlated.
    } else {
      return at_bound * POLY_MAXCOR;
    }
  }

  // Pad the thresholds with +/-Inf so each ordinal cell is a (half-infinite) rectangle. With
  // complete data the per-variable (full-column) thresholds tau_i/tau_j are shared by every
  // pair. With missing data each pair instead uses the thresholds of its own pairwise-complete
  // cases (the marginals of this contingency table): the thresholds and the table then come from
  // the same cases, as in polycor / psych. An empty marginal category gives a +/-Inf cut, which
  // the +/-POLY_CLAMP clamping turns into a zero-width (zero-probability, zero-count) band.
  rcut[0] = -POLY_INF; rcut[Ki] = POLY_INF;
  ccut[0] = -POLY_INF; ccut[Kj] = POLY_INF;
  if (use_local) {
    poly_fix_thresholds(tab, Ki, Kj, total, rcut, ccut);
  } else {
    for (int k = 0; k < Ki - 1; k++) rcut[k + 1] = tau_i[k];
    for (int k = 0; k < Kj - 1; k++) ccut[k + 1] = tau_j[k];
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

  // Smallest positive double, used only to keep log()/division finite for a cell whose
  // probability is genuinely below it - which, with the band taken from the nearer tail, means
  // a truly unrepresentable probability rather than one lost to cancellation. The conditioning
  // integral is otherwise non-negative.
  const double p_floor = std::numeric_limits<double>::min();

  // Evaluate the negative log-likelihood, the score dl/drho, and the expected (Fisher)
  // information at rho. The cell probabilities P_ab and their derivatives dP_ab/drho are built
  // from the SAME conditioning integral, so the score is exactly the derivative of the
  // quadrature objective (not of the unreachable exact probability) and the iteration
  // converges cleanly. Differentiating P_ab = int phi(x)[Phi((c1-rho x)/s) - Phi((c0-rho x)/s)] dx
  // under the integral sign gives d/drho Phi((c-rho x)/s) = phi((c-rho x)/s) (rho c - x) / s^3
  // (Plackett, 1954); an infinite column cut contributes zero. Only the inner column terms
  // depend on rho and are shared across the cells of a row. The information is the multinomial
  // expectation total * sum_ab (dP_ab/drho)^2 / P_ab (Olsson, 1979): positive, and needing no
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
        // Expected (Fisher) information sum_ab (dP)^2/P. A cell whose probability has
        // underflowed carries no information: the floored ratio would be meaningless (the
        // numerator is not floored with it) and, through the Newton step score/info, could
        // fake a vanishing step that stops the iteration far from the optimum. Skip those
        // cells - their true contribution is negligible. This bites only at |rho| within a
        // rounding unit of 1, where the true probability really is below the smallest
        // representable double; a merely tiny cell keeps its full relative accuracy and its
        // (correspondingly tiny) information.
        if (praw > p_floor) iacc += dP * dP / praw;
        double nab = tab[pbase + b];
        if (nab != 0.0) {
          nll -= nab * std::log(p);
          // The score must be the gradient of this nll. Once the probability underflows,
          // -nab*log(p) is pinned at -nab*log(p_floor) and so locally constant in rho, giving
          // it zero gradient; gating the score on the same condition as the information keeps
          // the Newton step score/info consistent with the objective being minimised.
          if (praw > p_floor) score += nab * dP / p;
        }
      }
    }
    info = total * iacc;
  };

  // Warm start at the Pearson correlation of the category codes (an inexpensive starting value
  // also used by lavaan): it is usually close to rho, so the scoring iteration starts inside the
  // basin and converges in a few steps instead of overshooting from rho = 0 as a cold start
  // would. Computed from the
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
  double cov = sab - sa * sb / total;
  double va = saa - sa * sa / total, vb = sbb - sb * sb / total;
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
  const double maxcor = POLY_MAXCOR;
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
    // Signed-tail CDF, differenced by poly_signed_band exactly as in poly_accum_node: these far
    // tails are the same near-impossible cells, and a band that cancelled to zero here would
    // silently drop a threshold's share of the cross-information.
    double rx = rho * x;
    for (int o = 0; o <= n_other; o++) colb[o] = std_norm_cdf_signed((cut_other[o] - rx) / s);
    double acc = 0.0;
    for (int o = 0; o < n_other; o++) {
      double band = poly_signed_band(colb[o], colb[o + 1]);
      std::size_t up = (std::size_t) k * stride_self + (std::size_t) o * stride_other;
      std::size_t lo = up + stride_self;               // the adjacent row/column of band k
      // Empty cells carry no cases and so contribute nothing; skipping them rather than
      // adding a 0 * ... term keeps A21 finite should such a cell sit at the probability
      // floor, where dx.rho (already dP/P, divided by P once more here) could overflow to
      // Inf and 0 * Inf = NaN -- same rationale as the A22 and diagonal accumulations.
      // Unlike those, a NaN here would not stay local: it spreads through the threshold
      // aggregates T into the influence of every cell of the pair, empty or not. Defensive
      // rather than load-bearing: dP/P is bounded once the band keeps its relative accuracy
      // (see std_norm_cdf_signed), because P and dP then underflow together.
      double term = 0.0;
      if (nab[up] != 0.0) {
        double pu = pmat[up] < p_floor ? p_floor : pmat[up];
        term += nab[up] * dxr[up] / pu;
      }
      if (nab[lo] != 0.0) {
        double pl = pmat[lo] < p_floor ? p_floor : pmat[lo];
        term -= nab[lo] * dxr[lo] / pl;
      }
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
                          bool nearest_pd) {
  if (acov != "none" && acov != "diag" && acov != "full") {
    Rcpp::stop("`acov` must be one of \"none\", \"diag\", or \"full\".");
  }
  const int nrow = x.nrow();
  const int p = x.ncol();
  // Unreachable defence: .polychoric() rejects p < 2 at the R level with a classed condition.
  if (p < 2) {
    Rcpp::stop("At least two variables are required for a correlation matrix.");
  }
  if (nrow < 2) {
    Rcpp::stop("At least two observations are required for polychoric correlations.");
  }

  // The R wrapper supplies 0-based consecutive codes. Keep that contract at the
  // registered native boundary as well: malformed negative or gapped codes would
  // otherwise index the frequency and contingency buffers out of bounds. A valid
  // maximum is necessarily below nrow because every consecutive category occurs.
  for (int j = 0; j < p; ++j) {
    std::vector<char> seen(static_cast<std::size_t>(nrow), 0);
    int maxc = -1;
    int n_seen = 0;
    for (int r = 0; r < nrow; ++r) {
      const int v = x(r, j);
      if (v == NA_INTEGER) continue;
      if (v < 0 || v >= nrow) {
        Rcpp::stop("Category codes must be non-negative, consecutive, and smaller than the number of observations.");
      }
      if (!seen[static_cast<std::size_t>(v)]) {
        seen[static_cast<std::size_t>(v)] = 1;
        ++n_seen;
      }
      if (v > maxc) maxc = v;
    }
    if (n_seen < 2) {
      Rcpp::stop("Every variable must have at least two observed categories.");
    }
    if (maxc + 1 != n_seen) {
      Rcpp::stop("Category codes must be consecutive from 0 within each variable.");
    }
  }

  const int* xp = x.begin();

  // Gauss-Legendre rules for the rectangle quadrature, computed once: the base 12-node rule and
  // the finer rule used to re-estimate near-collinear pairs (|rho| close to 1), where the base
  // rule under-resolves the narrow conditional transition.
  std::vector<double> gx, gw, gx_hi, gw_hi;
  gauss_legendre(POLY_GL_N, gx, gw);
  gauss_legendre(POLY_GL_N_HI, gx_hi, gw_hi);

  // Whether any cell is missing. With complete data every pair shares the per-variable
  // thresholds (the fast, exact path); with missing data each pair instead uses its own
  // pairwise-complete thresholds (computed inside the pair loop via the AS 241 quantile).
  bool any_missing = false;
  for (std::size_t k = 0; k < (std::size_t) nrow * p; k++) {
    if (xp[k] == NA_INTEGER) { any_missing = true; break; }
  }

  // Category counts and thresholds, computed once per variable. The codes are 0-based
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

  // Flat list of upper-triangle pairs. Each pair writes a distinct rho slot, so the
  // result is bit-for-bit deterministic.
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
  // Records, per pair, whether its contingency table is at a Frechet bound of its own marginals
  // (+1 upper / -1 lower, 0 otherwise): a boundary solution whose estimate is the domain endpoint
  // and whose asymptotic variance does not exist. Reported to R so the pairs can be named, and
  // consumed by the ACOV block below.
  std::vector<int> at_bound((std::size_t) npairs, 0);
  // Records, per pair, whether a binary pair's empty cell was continuity-corrected before
  // estimation (see poly_zero_correct_2x2). Such a pair IS estimated in the interior, so it is not
  // at_bound, but the estimate comes from a nudged table rather than from the observed counts, so
  // it has no asymptotic variance either and the ACOV block skips it alongside the bound pairs.
  std::vector<int> zero_corrected((std::size_t) npairs, 0);

  {
    // Scratch reused across pairs (no allocation in the loop).
    std::vector<double> tab((std::size_t) Kmax * Kmax);
    std::vector<double> pmat((std::size_t) Kmax * Kmax);
    std::vector<double> dpmat((std::size_t) Kmax * Kmax);
    std::vector<double> colcdf(Kmax + 1);
    std::vector<double> coldcdf(Kmax + 1);
    std::vector<double> rcut(Kmax + 1), ccut(Kmax + 1);
    std::vector<double> xnode((std::size_t) Kmax * POLY_GL_N_HI);
    std::vector<double> wpdf((std::size_t) Kmax * POLY_GL_N_HI);

    for (int t = 0; t < npairs; t++) {
      int i = pair_i[t];
      int j = pair_j[t];
      bool had_empty = false;
      int bnd = 0;
      bool zc = false;
      double r = polychoric_pair(xp, nrow, i, j, tau[i], tau[j], Kcat[i], Kcat[j],
                                 any_missing, gx, gw, tab, pmat, dpmat, colcdf, coldcdf,
                                 rcut, ccut, xnode, wpdf, had_empty, bnd, zc);
      at_bound[t] = bnd;
      zero_corrected[t] = zc ? 1 : 0;
      // The 12-node rule under-resolves the conditional transition near |rho| = 1, and the same
      // sharp transition over a wide tail band biases the estimate of a strongly-correlated table
      // with an empty cell -- where the 12-node estimate itself under-shoots, dropping below
      // POLY_REFINE_RHO so the near-collinear test alone would miss it. Re-estimate such pairs with
      // the finer rule: when the cheap estimate is near-collinear, or has an empty cell and is
      // already at least moderately correlated (POLY_EMPTY_REFINE_RHO). An empty cell at a low |rho|
      // carries no quadrature bias, so it keeps the cheap rule. A pair at a Frechet bound is
      // exempt: its estimate is the domain endpoint by construction, not a quadrature result.
      if (bnd == 0 && std::isfinite(r) &&
          (std::fabs(r) > POLY_REFINE_RHO ||
           (had_empty && std::fabs(r) > POLY_EMPTY_REFINE_RHO))) {
        // The second pass re-derives `had_empty`, the bound verdict, and the zero-cell correction
        // from the same table, so all three necessarily repeat what the first pass found; they are
        // received into throwaways rather than allowed to overwrite the values the branch was
        // chosen on.
        bool he2 = false;
        int bnd2 = 0;
        bool zc2 = false;
        r = polychoric_pair(xp, nrow, i, j, tau[i], tau[j], Kcat[i], Kcat[j],
                            any_missing, gx_hi, gw_hi, tab, pmat, dpmat, colcdf, coldcdf,
                            rcut, ccut, xnode, wpdf, he2, bnd2, zc2);
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
  // and equals diag(Gamma) by construction. This is the large-sample (per-observation, 1/N)
  // covariance: N * Gamma matches lavaan's per-case correlation NACOV, and lavaan's reported
  // vcov() differs only by the small-sample 1/(N - 1) versus 1/N convention (a factor of
  // (N - 1)/N that vanishes as N grows). When an ACOV is requested the R wrapper restricts
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
    // Zero-initialised, not merely sized: a pair at a Frechet bound is skipped below and never
    // writes its column, and an uninitialised column would feed garbage into the crossprod and so
    // into every OTHER pair's cross-covariances. Every estimable pair overwrites its own column in
    // full, so the fill costs one pass and changes nothing for them.
    if (acov == "full") IF_cor.zeros(Nc, npairs);

    {
      // Scratch reused across pairs.
      std::vector<double> pmat((std::size_t) Kmax * Kmax);
      std::vector<double> dpmat((std::size_t) Kmax * Kmax);
      std::vector<double> nab((std::size_t) Kmax * Kmax);
      std::vector<double> dxr((std::size_t) Kmax * Kmax);
      std::vector<double> cif((std::size_t) Kmax * Kmax);
      std::vector<double> colcdf(Kmax + 1), coldcdf(Kmax + 1), colb(Kmax + 1);
      std::vector<double> Ti(Kmax), Tj(Kmax), A21i(Kmax), A21j(Kmax);

      for (int t = 0; t < npairs; t++) {
        // A pair at a Frechet bound has no asymptotic variance to report. Its estimate is the
        // boundary of the parameter space, where the rho-information vanishes identically, so the
        // influence function (score - threshold correction) / A22 is 0/0: both parts underflow to
        // zero well before the boundary is reached, and the limit of the ratio is infinite. Any
        // finite value here would be an artefact of the probability floor rather than a variance -
        // which is exactly what differed by 30-odd orders of magnitude between platforms - and the
        // standard sandwich asymptotics do not apply at a boundary with singular information in any
        // case (Self & Liang, 1987; Rotnitzky, Cox, Bottai, & Robins, 2000). Fail closed with NA so
        // no consumer can silently weight or report noise; the DWLS weights refuse it and the
        // sandwich withholds the affected standard errors. Its influence column stays at the zero it
        // was initialised to, so the crossprod below stays finite for every OTHER pair, and this
        // pair's own row and column of Gamma are set to NA afterwards.
        //
        // A continuity-corrected binary pair is withheld for a different reason: its estimate is
        // interior, but it was computed from a nudged table rather than from the observed counts,
        // and the sandwich below would treat those nudged counts as if they were the data. That
        // approximation fails exactly where the correction is invoked -- simulated coverage of a
        // nominal 95% interval built this way ranges from 0 to 1 depending on the marginals -- so
        // the correction is a point-estimate device only and its variance is withheld too.
        if (at_bound[t] != 0 || zero_corrected[t] != 0) {
          acov_diag[t] = NA_REAL;
          continue;
        }
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
            // Empty cells contribute no cases, so skip them: at the probability floor d
            // could overflow to Inf, and 0 * Inf = NaN would silently zero the variance via
            // the invA22 guard below (same rationale as the diagonal accumulation).
            if (nab[idx] != 0.0) A22 += nab[idx] * d * d;
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
            // An empty cell carries no case and so contributes exactly nothing to the
            // sum over cases of IF^2. Skipping it rather than adding 0 * v * v also keeps
            // the diagonal finite should a structurally empty cell of a (near-)degenerate
            // pair have a model probability at the underflow floor, where dP/P would be
            // astronomically large and v * v would overflow to Inf -- 0 * Inf being NaN in
            // IEEE arithmetic. This is what the `full` scatter below already does implicitly
            // by indexing observed cases only, and is what makes diag == diag(Gamma) hold
            // everywhere.
            //
            // Keep this skip (and its A22 and A21 counterparts) even though no table now
            // reaches the overflow. With the band taken from the nearer tail, P and dP
            // underflow together and dP/P stays O(1e5) even in the deepest reachable tail, so
            // the skip is defensive: it is what makes the invariant hold by construction
            // rather than by the accident of a particular quadrature rule.
            if (nab[idx] != 0.0) vsum += nab[idx] * v * v;
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
      // A withheld pair's covariance with every other pair is as unavailable as its own variance,
      // so NA its whole row and column rather than let the zeroed influence column pass them off
      // as exact zeros. Done after the crossprod so the NAs cannot propagate into the pairs that
      // ARE estimable. This keeps diag(Gamma) equal to the "diag" vector, NA entries included.
      for (int t = 0; t < npairs; t++) {
        if (at_bound[t] != 0 || zero_corrected[t] != 0) {
          Gamma.row(t).fill(NA_REAL);
          Gamma.col(t).fill(NA_REAL);
        }
      }
      acov_out = Rcpp::wrap(Gamma);
    }
  }

  // `at_bound` and `zero_corrected` are reported as logicals over the pairs (utils::combn order,
  // i.e. this function's own pair loop) so the R wrapper can name the affected variable pairs in
  // its two diagnostics. They are mutually exclusive: a corrected pair is estimated in the interior
  // and so is not at a bound.
  //
  // `bound_rho` carries the boundary estimate itself, signed, for the pairs that are at one (0
  // elsewhere). The R warning has to state the value it reported, and the sign cannot be recovered
  // from the returned matrix: with `nearest_pd` the matrix is projected afterwards, so the entry is
  // no longer +/-POLY_MAXCOR. Carrying it here also keeps the constant in one place instead of
  // repeating it as a literal in the message.
  Rcpp::LogicalVector bound_out(npairs);
  Rcpp::LogicalVector corrected_out(npairs);
  Rcpp::NumericVector bound_rho_out(npairs);
  for (int t = 0; t < npairs; t++) {
    bound_out[t] = (at_bound[t] != 0);
    corrected_out[t] = (zero_corrected[t] != 0);
    bound_rho_out[t] = at_bound[t] * POLY_MAXCOR;
  }

  return Rcpp::List::create(
    Rcpp::Named("R") = Rmat,
    Rcpp::Named("pd_adjusted") = pd_adjusted,
    Rcpp::Named("acov") = acov_out,
    Rcpp::Named("at_bound") = bound_out,
    Rcpp::Named("zero_corrected") = corrected_out,
    Rcpp::Named("bound_rho") = bound_rho_out);
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
