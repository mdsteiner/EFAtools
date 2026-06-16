## Unweighted Least Squares Estimation of Factor Loadings
## Thin fitter: returns the raw ULS results; shared post-processing happens in
## .finalize_fit() / .estimate_model().
.ULS <- function(x, n_factors) {

  # Get correlation matrix entered or created in EFA
  R <- x

  uls <- .fit_uls(R, n_factors)

  L <- uls$loadings
  orig_R <- R
  h2 <- diag(L %*% t(L))
  diag(R) <- h2

  # ULS objective: sum of squared off-diagonal residuals
  Fm <- orig_R - (L %*% t(L))
  Fm <- sum(Fm[upper.tri(Fm)] ^ 2)

  # raw fit, finalized by .estimate_model()
  list(
    L = L,
    h2 = h2,
    psi = uls$res$par,
    Fm = Fm,
    iter = uls$res$counts[1],
    convergence = uls$res$convergence,
    orig_R = orig_R,
    R_final = R
  )
}

# function to obtain the uls fit; adapted from the psych package
.fit_uls <- function(R, n_fac) {

  start <- diag(R) - .smc_start(R)

  # Bounded L-BFGS-B over the uniquenesses, run entirely in C++ (roptim drives
  # R's lbfgsb routine, so the box constraints and control match stats::optim).
  uls <- .fit_uls_cpp(R, n_fac, start, .uniqueness_floor)

  res <- list(par = uls$psi, value = uls$Fm, counts = c(uls$iter, NA_integer_),
              convergence = uls$convergence)

  list(loadings = uls$loadings, res = res, R = R)
}
