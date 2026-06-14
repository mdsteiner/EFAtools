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

  res <- stats::optim(start, .uls_residuals, gr = .grad_uls, method = "L-BFGS-B",
               lower = .uniqueness_floor, upper = 1,
               control = c(list(fnscale = 1, parscale = rep(0.01, length(start)))),
               R = R, n_fac = n_fac)

  Lambda <- .FAout_wls(res$par, R, n_fac)

  result <- list(loadings = Lambda, res = res, R = R)
}

.FAout_wls <-  function(psi, R, n_fac) {
  diag(R) <- 1 - psi
  E <- eigen(R, symmetric = TRUE)

  L <- E$vectors[,1L:n_fac,drop=FALSE] %*%
    diag(sqrt(pmax(E$values[1L:n_fac,drop=FALSE], 0)), n_fac)
  return(L)
}
