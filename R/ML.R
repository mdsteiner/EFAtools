## Maximum Likelihood Estimation of Factor Loadings
## Thin fitter: returns the raw ML results; shared post-processing happens in
## .finalize_fit() / .estimate_model().
.ML <- function(x, n_factors, start_method = c("psych", "factanal")) {

  # Get correlation matrix entered or created in EFA
  R <- x

  ml <- .fit_ml(R, n_factors, start_method)

  L <- ml$loadings
  orig_R <- R
  h2 <- diag(L %*% t(L))
  diag(R) <- h2

  # raw fit, finalized by .estimate_model()
  list(
    L = L,
    h2 = h2,
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
    R.smc <- (1 - 1 / diag(solve(R)))
    if((sum(R.smc) == n_fac) && (n_fac > 1)) {
      start <- rep(.5, n_fac)
    }  else {
      start <- diag(R)- R.smc
    }
  } else if (start_method == "factanal") {
    start <- (1 - 0.5 * n_fac / ncol(R)) / diag(solve(R))
  }

  res <- stats::optim(start, .error_ml, gr = .grad_ml, method = "L-BFGS-B",
                      lower = .005, upper = 1,
                      control = c(list(fnscale = 1,
                                       parscale = rep(0.01, length(start)))),
                      R = R, n_fac = n_fac)

  Lambda <- .FAout(res$par, R, n_fac)

  result <- list(loadings = Lambda, res = res, R = R)

  result
}

# taken from factanal
.FAout <- function(psi, R, n_fac) {
  sc <- diag(1 / sqrt(psi))
  Rs <- sc %*% R %*% sc
  E <- eigen(Rs, symmetric = TRUE)
  L <- E$vectors[, seq_len(n_fac), drop = FALSE]
  load <- L %*% diag(sqrt(pmax(E$values[seq_len(n_fac)] - 1, 0)),
                     n_fac)
  diag(sqrt(psi)) %*% load
}
