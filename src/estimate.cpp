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

// Objective and gradient of the maximum-likelihood factor model, fused so that a
// single eigendecomposition of the rescaled correlation matrix serves both. The
// objective is the discrepancy minimised by stats::factanal()/psych; the gradient
// is its derivative with respect to the uniquenesses psi (adapted from
// stats::factanal()). A one-slot cache keyed on psi holds the last decomposition's
// objective, gradient, and factor loadings, so the optimiser's paired
// value/gradient evaluations at the same psi decompose only once.
class FaFunctor : public roptim::Functor {
public:
  FaFunctor(const arma::mat& R, int n_fac) : R_(R), n_fac_(n_fac) {
    // The eigen-extraction reads the largest n_fac eigenpairs and would index past
    // the available eigenvalues if n_fac >= ncol(R); reject that up front.
    if (n_fac < 1 || static_cast<arma::uword>(n_fac) >= R.n_cols) {
      Rcpp::stop("n_fac must be at least 1 and smaller than the number of "
                 "variables (got n_fac = %d for %d variables).",
                 n_fac, static_cast<int>(R.n_cols));
    }
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
    cached_f_ = -(arma::accu(arma::log(e) - e) - n_fac_ + (int)R_.n_rows);

    // factor loadings from the leading n_fac eigenpairs (also used by the gradient)
    arma::vec top = lambda.head(n_fac_) - 1.0;
    top.elem(arma::find(top < 0)).zeros();
    arma::mat load = arma::diagmat(arma::sqrt(psi)) *
      (V.head_cols(n_fac_) * arma::diagmat(arma::sqrt(top)));

    // gradient of the objective with respect to psi
    arma::mat g = load * load.t() + arma::diagmat(psi) - R_;
    cached_g_ = g.diag() / arma::pow(psi, 2);

    cached_load_ = std::move(load);
    cached_psi_ = psi;
    valid_ = true;
  }

private:
  const arma::mat& R_;
  const int n_fac_;

  arma::vec cached_psi_;
  double cached_f_ = 0.0;
  arma::vec cached_g_;
  arma::mat cached_load_;
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
