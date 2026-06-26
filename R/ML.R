## Maximum Likelihood Estimation of Factor Loadings
## Thin fitter: returns the raw ML results; shared post-processing happens in
## .finalize_fit() / .estimate_model().
.ML <- function(x, n_factors, start_method = c("psych", "factanal")) {

  # Get correlation matrix entered or created in EFA
  R <- x

  ml <- .fit_ml(R, n_factors, start_method)

  L <- ml$loadings
  orig_R <- R
  h2 <- rowSums(L^2)
  diag(R) <- h2

  # raw fit, finalized by .estimate_model()
  list(
    L = L,
    h2 = h2,
    psi = ml$res$par,
    Fm = ml$res$value,
    iter = ml$res$counts[1],
    convergence = ml$res$convergence,
    orig_R = orig_R,
    R_final = R,
    settings = list(start_method = start_method)
  )
}

# function to obtain the ML fit; adapted from the psych package
.fit_ml <- function(R, n_fac, start_method) {

  if (start_method == "psych") {
    R.smc <- .smc_start(R)
    # psych::fa()-style start: the uniqueness start diag(R) - smc, or a flat 0.5
    # fallback sized to the variables (the optimiser parameter is the length-p
    # uniqueness vector).
    if((sum(R.smc) == n_fac) && (n_fac > 1)) {
      start <- rep(.5, ncol(R))
    }  else {
      start <- diag(R)- R.smc
    }
  } else if (start_method == "factanal") {
    start <- (1 - 0.5 * n_fac / ncol(R)) / diag(solve(R))
  }

  # Bounded L-BFGS-B over the uniquenesses, run entirely in C++ (roptim drives
  # R's lbfgsb routine, so the box constraints and control match stats::optim).
  ml <- .fit_ml_cpp(R, n_fac, start, .uniqueness_floor)

  res <- list(par = ml$psi, value = ml$Fm, counts = c(ml$iter, NA_integer_),
              convergence = ml$convergence)

  list(loadings = ml$loadings, res = res, R = R)
}
