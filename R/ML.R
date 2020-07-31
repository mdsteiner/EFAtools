## Maximum Likelihood Estimation of Factor Loadings
.ML <- function(x, n_factors, N = NA, start_method = c("psych", "factanal")) {

  # Get correlation matrix entered or created in EFA
  R <- x

  ml <- .fit_ml(R, n_factors, start_method)

  L <- ml$loadings
  orig_R <- R
  h2 <- diag(L %*% t(L))
  diag(R) <- h2

  # reverse the sign of loadings
  if (n_factors > 1) {
    signs <- sign(colSums(L))
    signs[signs == 0] <- 1
    L <- L %*% diag(signs)
  } else {
    if (sum(L) < 0) {
      L <- -as.matrix(L)
    } else {
      L <- as.matrix(L)
    }

  }

  if (!is.null(colnames(orig_R))) {
    # name the loading matrix so the variables can be identified
    rownames(L) <- colnames(orig_R)
  } else {
    varnames <- paste0("V", seq_len(ncol(orig_R)))
    colnames(orig_R) <- varnames
    rownames(orig_R) <- varnames
    rownames(L) <- varnames
  }

  colnames(L) <- paste0("F", seq_len(n_factors))

  vars_accounted <- .compute_vars(L_unrot = L, L_rot = L)

  colnames(vars_accounted) <- colnames(L)

  # compute fit indices
  fit_ind <- .gof(L, orig_R, N, "ML", ml$res$value)

  # create the output object
  class(L) <- "LOADINGS"

  # store the settings used:
  settings <- list(
    start_method = start_method
  )

  # Create and name communalities
  h2 <-  diag(L %*% t(L))
  names(h2) <- colnames(orig_R)

  # Create output
  output <- list(
    orig_R = orig_R,
    h2 = h2,
    orig_eigen = eigen(orig_R, symmetric = TRUE)$values,
    final_eigen = eigen(R, symmetric = TRUE)$values,
    iter = ml$res$counts[1],
    convergence = ml$res$convergence,
    unrot_loadings = L,
    vars_accounted = vars_accounted,
    fit_indices = fit_ind,
    settings = settings
  )

  output
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

# .error_ml2 <- function(psi, R, n_fac)
# {
#   sc <- diag(1/sqrt(psi))
#   Rs <- sc %*% R %*% sc
#   E <- eigen(Rs, symmetric = TRUE, only.values = TRUE)
#   e <- E$values[-(1:n_fac)]
#   e <- sum(log(e) - e) - n_fac + nrow(R)
#   -e
# }
# .grad_ml2 <- function(psi, R, n_fac) {
#   sc <- diag(1 / sqrt(psi))
#   Rs <- sc %*% R %*% sc
#   E <- eigen(Rs, symmetric = TRUE)
#   L <- E$vectors[, 1:n_fac, drop = FALSE]
#   load <- L %*% diag(sqrt(pmax(E$values[1:n_fac] - 1, 0)), n_fac)
#   load <- diag(sqrt(psi)) %*% load
#   g <- load %*% t(load) + diag(psi) - R     # g <- model - data
#   diag(g) / psi^2                             #normalized
# }

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
