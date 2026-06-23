// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(roptim)]]
#include <RcppArmadillo.h>
#include <roptim.h>

using namespace Rcpp;
using namespace arma;

// Symmetric eigendecomposition that fails loudly instead of silently leaving the
// outputs empty (which the callers would then index out of bounds). During the
// optimization the decomposition is taken on a rescaled matrix that can go
// non-finite or non-symmetric for a pathological correlation matrix, so turn a
// failed arma::eig_sym into a catchable R error rather than undefined behaviour.
static void eig_sym_checked(arma::vec& eigval, arma::mat& eigvec,
                            const arma::mat& X) {
  if (!arma::eig_sym(eigval, eigvec, X)) {
    Rcpp::stop("Eigendecomposition failed during factor extraction; the "
               "correlation matrix is not finite or not symmetric.");
  }
}

// Shared constructor guard: the eigen-extraction reads the largest n_fac eigenpairs and
// would index past the available eigenvalues if n_fac is out of range, so reject that up
// front with a clean message rather than an opaque bounds error.
static void check_n_fac(int n_fac, arma::uword n_cols) {
  if (n_fac < 1 || static_cast<arma::uword>(n_fac) >= n_cols) {
    Rcpp::stop("n_fac must be at least 1 and smaller than the number of "
               "variables (got n_fac = %d for %d variables).",
               n_fac, static_cast<int>(n_cols));
  }
}

// Shared constructor guard for the weighted estimators: the per-element weight matrix must
// match the correlation matrix.
static void check_weights(const arma::mat& W, const arma::mat& R) {
  if (W.n_rows != R.n_rows || W.n_cols != R.n_cols) {
    Rcpp::stop("The weight matrix must be square and match the correlation "
               "matrix dimensions.");
  }
}

// Shared scaffolding for the bounded L-BFGS-B factor estimators. Each estimator's
// objective and gradient are fused so that a single eigendecomposition of the
// rescaled correlation matrix per psi serves both, and a one-slot cache keyed on
// psi holds the last decomposition's objective, gradient, and factor loadings, so
// the optimiser's paired value/gradient evaluations at the same psi decompose only
// once. This base owns the cache, the value/gradient/loadings accessors, and the
// n_fac bound check; each estimator overrides compute() with its criterion math.
class FactorFunctor : public roptim::Functor {
public:
  FactorFunctor(const arma::mat& R, int n_fac) : R_(R), n_fac_(n_fac) {
    check_n_fac(n_fac, R.n_cols);
  }

  double operator()(const arma::vec& psi) override {
    refresh(psi);
    return cached_f_;
  }

  void Gradient(const arma::vec& psi, arma::vec& grad) override {
    refresh(psi);
    grad = cached_g_;
  }

  // Factor loadings of the most recently evaluated psi; call refresh() with the
  // optimum first so they correspond to the solution.
  const arma::mat& loadings() const { return cached_load_; }

  void refresh(const arma::vec& psi) {
    if (valid_ && psi.n_elem == cached_psi_.n_elem &&
        arma::all(psi == cached_psi_)) {
      return;
    }
    compute(psi, cached_f_, cached_g_, cached_load_);
    cached_psi_ = psi;
    valid_ = true;
  }

protected:
  // Criterion-specific objective f, gradient g, and factor loadings load from a
  // single eigendecomposition of the rescaled correlation matrix at psi.
  virtual void compute(const arma::vec& psi, double& f, arma::vec& g,
                       arma::mat& load) = 0;

  const arma::mat& R_;
  const int n_fac_;

private:
  arma::vec cached_psi_;
  double cached_f_ = 0.0;
  arma::vec cached_g_;
  arma::mat cached_load_;
  bool valid_ = false;
};

// Maximum-likelihood factor model. The objective is the discrepancy minimised by
// stats::factanal()/psych; the gradient is its derivative with respect to the
// uniquenesses psi (adapted from stats::factanal()).
class FaFunctor : public FactorFunctor {
public:
  using FactorFunctor::FactorFunctor;

protected:
  void compute(const arma::vec& psi, double& f, arma::vec& g,
               arma::mat& load) override {
    // One eigendecomposition of Rs = D^-1/2 R D^-1/2 (D = diag(psi)) feeds the
    // objective, gradient, and loadings.
    arma::mat sc = arma::diagmat(1.0 / sqrt(psi));
    arma::mat Rs = sc * R_ * sc;
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym_checked(eigval, eigvec, Rs);

    // descending order (largest eigenvalue first), matching R's eigen()
    arma::vec lambda = arma::flipud(eigval);
    arma::mat V = arma::fliplr(eigvec);

    // objective: -sum(log(e) - e) over the smallest p - n_fac eigenvalues
    arma::vec e = lambda.tail(lambda.n_elem - n_fac_);
    f = -(arma::accu(arma::log(e) - e) - n_fac_ + (int)R_.n_rows);

    // factor loadings from the leading n_fac eigenpairs (also used by the gradient)
    arma::vec top = lambda.head(n_fac_) - 1.0;
    top.elem(arma::find(top < 0)).zeros();
    load = arma::diagmat(arma::sqrt(psi)) *
      (V.head_cols(n_fac_) * arma::diagmat(arma::sqrt(top)));

    // gradient of the objective with respect to psi
    arma::mat gmat = load * load.t() + arma::diagmat(psi) - R_;
    g = gmat.diag() / arma::pow(psi, 2);
  }
};

// Unweighted-least-squares (minres) factor model. The objective is the sum of
// squared strictly-lower off-diagonal residuals of the reproduced correlation
// matrix (the diagonal is absorbed by the free uniquenesses); the gradient is the
// diagonal of the full residual (adapted from the psych package).
//
// The objective and gradient eigen-floor the shared decomposition differently (the
// gradient clips negative eigenvalues to 0; the objective lifts eigenvalues below
// machine epsilon to eps*100), so compute() forms a loading matrix for each. The
// gradient's floor-at-0 loadings are also the reported solution, the extraction
// stats::eigen() would give at the optimum.
class UlsFunctor : public FactorFunctor {
public:
  using FactorFunctor::FactorFunctor;

protected:
  void compute(const arma::vec& psi, double& f, arma::vec& g,
               arma::mat& load) override {
    // One eigendecomposition of Rs = R with diagonal replaced by 1 - psi
    // (= R - diag(psi) since R is a correlation matrix) feeds both criteria.
    arma::mat Rs = R_;
    Rs.diag() -= psi;
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym_checked(eigval, eigvec, Rs);

    // descending order (largest eigenvalue first), matching R's eigen(); keep the
    // leading n_fac eigenpairs
    arma::vec lambda = arma::vec(arma::flipud(eigval)).head(n_fac_);
    arma::mat V = arma::mat(arma::fliplr(eigvec)).head_cols(n_fac_);

    // gradient: clip negative eigenvalues to 0; gradient is diag(LL' + diag(psi) - R)
    arma::vec lam_g = lambda;
    lam_g.elem(arma::find(lam_g < 0)).zeros();
    arma::mat load_g = V * arma::diagmat(arma::sqrt(lam_g));
    arma::mat gmat = load_g * load_g.t() + arma::diagmat(psi) - R_;
    g = gmat.diag();

    // objective: lift eigenvalues below machine epsilon to eps*100, then sum the
    // squared strictly-lower off-diagonal residuals (the diagonal, absorbed by the
    // uniquenesses, is excluded so objective and gradient stay consistent)
    arma::vec lam_f = lambda;
    lam_f.elem(arma::find(lam_f < arma::datum::eps)).fill(arma::datum::eps * 100);
    arma::mat load_f = V * arma::diagmat(arma::sqrt(lam_f));
    arma::mat resid = arma::trimatl(Rs - load_f * load_f.t(), -1);
    f = arma::accu(arma::square(resid));

    load = std::move(load_g);
  }
};

// Warm start for the DWLS optimiser. This is the "weighted-ULS" solution: it minimises
// the same weighted off-diagonal objective but over the uniquenesses psi, reading the
// loadings off the top-n_fac eigendecomposition of R - diag(psi). The eigen-extraction
// only approximates the weighted low-rank fit, so this is a few per cent from the
// free-loadings optimum -- but, unlike the unweighted ULS solution, it lies in the global
// basin of the weighted objective, so polishing the free loadings from here reaches the
// DWLS optimum (starting the free-loadings optimiser at the unweighted ULS solution can
// instead converge to a higher local minimum). A numerical gradient (roptim's default)
// is adequate for a warm start.
class EigenPsiDwlsFunctor : public roptim::Functor {
public:
  EigenPsiDwlsFunctor(const arma::mat& R, const arma::mat& W, int n_fac)
      : R_(R), W_(W), n_fac_(n_fac) {
    check_n_fac(n_fac, R.n_cols);
    check_weights(W, R);
  }

  double operator()(const arma::vec& psi) override {
    arma::mat L = loadings_at(psi);
    arma::mat E = R_ - L * L.t();
    arma::mat WE = W_ % E;
    WE.diag().zeros();  // off-diagonal weighting only, regardless of W's diagonal
    return 0.5 * arma::accu(WE % E);
  }

  // Loadings from the leading n_fac eigenpairs of R - diag(psi), negative eigenvalues
  // clipped to zero (the standard ULS/minres extraction).
  arma::mat loadings_at(const arma::vec& psi) const {
    arma::mat A = R_;
    A.diag() -= psi;
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym_checked(eigval, eigvec, A);
    arma::vec lam = arma::vec(arma::flipud(eigval)).head(n_fac_);
    arma::mat V = arma::mat(arma::fliplr(eigvec)).head_cols(n_fac_);
    lam.elem(arma::find(lam < 0)).zeros();
    return V * arma::diagmat(arma::sqrt(lam));
  }

private:
  const arma::mat& R_;
  const arma::mat& W_;
  const int n_fac_;
};

// Diagonally-weighted least squares (DWLS) factor model. Unlike the ML and ULS
// functors -- which optimise the uniquenesses and read the loadings off an
// eigendecomposition of the reduced correlation matrix -- DWLS optimises the loadings
// L directly: under per-element weighting the best weighted low-rank fit of the
// reduced matrix has no eigendecomposition form (Srebro & Jaakkola, 2003), so an
// eigen-extraction would minimise the unweighted criterion and only approximate DWLS.
// The objective is the weighted strictly-lower off-diagonal residual sum of squares
//   F(L) = sum_{i<j} W_ij (R_ij - (L L')_ij)^2,
// with W the symmetric per-element weight matrix (W_ij = 1 / Var(rho_hat_ij)) and zero
// diagonal, so F = 0.5 * accu(W % (R - L L')^2). Its gradient with respect to L is
//   dF/dL = -2 (W % (R - L L')) L
// (the strictly-lower sum's factor; -4 would double-count the symmetric off-diagonals).
// A one-slot cache keyed on the parameter vector fuses the value and gradient so the
// optimiser's paired evaluations form the residual once.
class DwlsFunctor : public roptim::Functor {
public:
  DwlsFunctor(const arma::mat& R, const arma::mat& W, int n_fac)
      : R_(R), W_(W), p_(R.n_rows), n_fac_(n_fac) {
    check_n_fac(n_fac, R.n_cols);
    check_weights(W, R);
  }

  double operator()(const arma::vec& par) override {
    refresh(par);
    return cached_f_;
  }

  void Gradient(const arma::vec& par, arma::vec& grad) override {
    refresh(par);
    grad = cached_g_;
  }

private:
  void refresh(const arma::vec& par) {
    if (valid_ && par.n_elem == cached_par_.n_elem &&
        arma::all(par == cached_par_)) {
      return;
    }
    arma::mat L(par.memptr(), p_, n_fac_);
    arma::mat E = R_ - L * L.t();
    arma::mat WE = W_ % E;
    WE.diag().zeros();  // off-diagonal weighting only, regardless of W's diagonal
    cached_f_ = 0.5 * arma::accu(WE % E);
    cached_g_ = arma::vectorise(-2.0 * (WE * L));
    cached_par_ = par;
    valid_ = true;
  }

  const arma::mat& R_;
  const arma::mat& W_;
  const arma::uword p_;
  const int n_fac_;
  arma::vec cached_par_;
  double cached_f_ = 0.0;
  arma::vec cached_g_;
  bool valid_ = false;
};

// Raw result of the bounded optimiser, before the method-specific post-processing.
struct OptResult {
  arma::vec par;
  double value;
  int iter;
  int convergence;
};

// Bounded L-BFGS-B over the uniquenesses psi, shared by the estimators. roptim
// drives R's own lbfgsb routine, so the control mirrors stats::optim() exactly:
// parscale/fnscale plus the L-BFGS-B defaults factr/pgtol/lmm/maxit. The bounds
// and parscale are sized to the parameter vector. The functor supplies the
// method-specific objective and gradient.
template <typename Functor>
static OptResult run_lbfgsb(Functor& f, const arma::vec& start, double lower) {
  const arma::uword p = start.n_elem;

  roptim::Roptim<Functor> opt("L-BFGS-B");
  opt.set_lower(arma::vec(p).fill(lower));
  opt.set_upper(arma::vec(p).fill(1.0));
  opt.control.fnscale  = 1.0;
  opt.control.parscale = arma::vec(p).fill(0.01);
  opt.control.factr    = 1e7;
  opt.control.pgtol    = 0.0;
  opt.control.lmm      = 5;
  opt.control.maxit    = 100;

  arma::vec par = start;
  opt.minimize(f, par);

  return OptResult{par, opt.value(), opt.fncount(), opt.convergence()};
}

// Fit the maximum-likelihood factor model entirely in C++. Returns the loadings at
// the optimum together with the uniquenesses, objective value, evaluation count,
// and convergence code.
// [[Rcpp::export(.fit_ml_cpp)]]
Rcpp::List fit_ml_cpp(const arma::mat& R, const int n_fac, arma::vec start,
                      const double lower) {
  FaFunctor f(R, n_fac);
  OptResult res = run_lbfgsb(f, start, lower);

  // Evaluate at the returned optimum so the cached loadings match it (the
  // optimiser's last evaluation may be a rejected trial point).
  f.refresh(res.par);

  return Rcpp::List::create(
    Rcpp::Named("loadings")    = f.loadings(),
    Rcpp::Named("psi")         = res.par,
    Rcpp::Named("Fm")          = res.value,
    Rcpp::Named("iter")        = res.iter,
    Rcpp::Named("convergence") = res.convergence);
}

// Maximum-likelihood objective and gradient exposed for the estimator guard tests;
// thin wrappers over FaFunctor so the discrepancy/gradient math has a single
// source. The constructor enforces the n_fac bound and refresh() the finite-matrix
// bound, which is what these tests assert.
// [[Rcpp::export(.error_ml)]]
double error_ml(arma::vec psi, const arma::mat& R, const int n_fac) {
  FaFunctor f(R, n_fac);
  return f(psi);
}

// [[Rcpp::export(.grad_ml)]]
arma::vec grad_ml(arma::vec psi, const arma::mat& R, const int n_fac) {
  FaFunctor f(R, n_fac);
  arma::vec grad;
  f.Gradient(psi, grad);
  return grad;
}

// Fit the unweighted-least-squares (minres) factor model entirely in C++. Returns
// the loadings at the optimum together with the uniquenesses, objective value,
// evaluation count, and convergence code.
// [[Rcpp::export(.fit_uls_cpp)]]
Rcpp::List fit_uls_cpp(const arma::mat& R, const int n_fac, arma::vec start,
                       const double lower) {
  UlsFunctor f(R, n_fac);
  OptResult res = run_lbfgsb(f, start, lower);

  // Evaluate at the returned optimum so the cached loadings match it (the
  // optimiser's last evaluation may be a rejected trial point).
  f.refresh(res.par);

  return Rcpp::List::create(
    Rcpp::Named("loadings")    = f.loadings(),
    Rcpp::Named("psi")         = res.par,
    Rcpp::Named("Fm")          = res.value,
    Rcpp::Named("iter")        = res.iter,
    Rcpp::Named("convergence") = res.convergence);
}

// Unweighted-least-squares objective and gradient exposed for the estimator guard
// tests; thin wrappers over UlsFunctor so the criterion math has a single source.
// The constructor enforces the n_fac bound and refresh() the finite-matrix bound,
// which is what these tests assert.
// [[Rcpp::export(.uls_residuals)]]
double uls_residuals(arma::vec psi, const arma::mat& R, const int n_fac) {
  UlsFunctor f(R, n_fac);
  return f(psi);
}

// [[Rcpp::export(.grad_uls)]]
arma::vec grad_uls(arma::vec psi, const arma::mat& R, const int n_fac) {
  UlsFunctor f(R, n_fac);
  arma::vec grad;
  f.Gradient(psi, grad);
  return grad;
}

// Fit the diagonally-weighted least squares factor model entirely in C++. A weighted-ULS
// warm start minimises the weighted off-diagonal objective over the uniquenesses psi
// (loadings from the eigendecomposition), which lands in the global basin; the free
// loadings are then polished with an unconstrained quasi-Newton (BFGS) run to the true
// DWLS optimum. The objective is coercive (large loadings inflate the residual), so the
// unconstrained polish stays bounded; improper (Heywood) solutions are left for the caller
// to detect, as in lavaan. The optimum is reoriented to the principal-axis frame (L'L
// diagonal, descending) via a thin SVD so the unrotated solution matches the ML/ULS
// extraction convention; this is a pure reorientation (L L' is unchanged) and so leaves the
// objective untouched. Returns the loadings at the optimum together with the weighted
// objective value, polish evaluation count, and convergence code.
// [[Rcpp::export(.fit_dwls_cpp)]]
Rcpp::List fit_dwls_cpp(const arma::mat& R, const int n_fac, const arma::mat& W) {
  // weighted-ULS warm start over psi (numerical gradient; box [floor, 1])
  EigenPsiDwlsFunctor wf(R, W, n_fac);
  roptim::Roptim<EigenPsiDwlsFunctor> wopt("L-BFGS-B");
  wopt.set_lower(arma::vec(R.n_rows).fill(0.005));
  wopt.set_upper(arma::vec(R.n_rows).fill(1.0));
  wopt.control.maxit = 200;
  // squared-multiple-correlation start (uniqueness = 1 / diag(R^-1)), clamped off the box.
  // Use a non-throwing inverse so an indefinite (e.g. an unsmoothed bootstrap resample) but
  // invertible R still yields a warm start, matching how ULS/ML tolerate such matrices
  // (their eigen extraction clips negative eigenvalues); the start only needs to be roughly
  // right. Fall back to a flat 0.5 start if R is singular.
  arma::mat Rinv;
  bool inv_ok = arma::inv_sympd(Rinv, R) || arma::inv(Rinv, R);
  arma::vec psi = inv_ok
    ? arma::clamp(1.0 / arma::diagvec(Rinv), 0.05, 0.95)
    : arma::vec(R.n_rows).fill(0.5);
  wopt.minimize(wf, psi);
  arma::mat start_L = wf.loadings_at(psi);

  // polish the free loadings to the DWLS optimum
  DwlsFunctor f(R, W, n_fac);
  roptim::Roptim<DwlsFunctor> opt("BFGS");
  opt.control.maxit  = 500;
  opt.control.reltol = 1e-12;
  arma::vec par = arma::vectorise(start_L);
  opt.minimize(f, par);

  arma::mat L(par.memptr(), R.n_rows, n_fac);
  arma::mat U;
  arma::vec s;
  arma::mat V;
  arma::svd_econ(U, s, V, L);
  arma::mat Lc = U * arma::diagmat(s);

  arma::mat E  = R - Lc * Lc.t();
  arma::mat WE = W % E;
  WE.diag().zeros();  // off-diagonal weighting only, matching DwlsFunctor (the minimised objective)
  double Fm = 0.5 * arma::accu(WE % E);

  return Rcpp::List::create(
    Rcpp::Named("loadings")    = Lc,
    Rcpp::Named("Fm")          = Fm,
    Rcpp::Named("iter")        = opt.fncount(),
    Rcpp::Named("convergence") = opt.convergence());
}

// DWLS objective and gradient exposed for the estimator guard tests; thin wrappers over
// DwlsFunctor so the criterion math has a single source. `par` is the vectorised loading
// matrix. The constructor enforces the n_fac and weight-dimension bounds.
// [[Rcpp::export(.dwls_residuals)]]
double dwls_residuals(arma::vec par, const arma::mat& R, const int n_fac,
                      const arma::mat& W) {
  DwlsFunctor f(R, W, n_fac);
  return f(par);
}

// [[Rcpp::export(.grad_dwls)]]
arma::vec grad_dwls(arma::vec par, const arma::mat& R, const int n_fac,
                    const arma::mat& W) {
  DwlsFunctor f(R, W, n_fac);
  arma::vec grad;
  f.Gradient(par, grad);
  return grad;
}
