## Unweighted Least Squares Estimation of Factor Loadings
.ULS <- function(x, n_factors, N = NA) {

  # Get correlation matrix entered or created in EFA
  R <- x

  uls <- .fit_uls(R, n_factors)

  L <- uls$loadings
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

  Fm <- orig_R - (L %*% t(L))
  Fm <- sum(Fm[upper.tri(Fm)] ^ 2)

  # compute fit indices
  fit_ind <- .gof(L, orig_R, N, "ULS", Fm)


  # create the output object
  class(L) <- "LOADINGS"

  # Create and name communalities
  h2 <-  diag(L %*% t(L))
  names(h2) <- colnames(orig_R)

  # Create output
  output <- list(
    orig_R = orig_R,
    h2 = h2,
    orig_eigen = eigen(orig_R, symmetric = TRUE)$values,
    final_eigen = eigen(R, symmetric = TRUE)$values,
    iter = uls$res$counts[1],
    convergence = uls$res$convergence,
    unrot_loadings = L,
    vars_accounted = vars_accounted,
    fit_indices = fit_ind
  )

  output
}

# function to obtain the uls fit; adapted from the psych package
.fit_uls <- function(R, n_fac) {

  start <- diag(R) - (1 - 1 / diag(solve(R)))

  res <- stats::optim(start, .uls_residuals, gr = .grad_uls, method = "L-BFGS-B",
               lower = .005, upper = 1,
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
