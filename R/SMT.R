# Sequential Model Test
# Chi sq & RMSEA(see grieder&grob)
SMT <- function(x, N){

  max_fac <- min(ncol(x)-2, ncol(x)/2 + 2)
  ps <- rep(0, max_fac)

  # see how this can be done
  zeromod <- EFA(x, n_factors = 1, method = "ML", rotation = "none", N = N)

  if (stats::pchisq(zeromod$null.chisq, zeromod$null.dof, lower.tail = F) < 0.05){

  # end done

    n_fac_chi <- 0

  } else {

    for (i in 1:max_fac) {

      temp <- EFA(x, n_factors = i, method = "ML", rotate ="none", N = N)
      ps[i] <- stats::pchisq(temp$fit_indices$chi, temp$fit_indices$df,
                             lower.tail = F)

    }

    if(any(ps > 0.05, na.rm = T)) {

      n_fac_chi <- which(ps > 0.05)[1]

    } else {

      n_fac_chi <- 0
    }

  }

  list(n_fac_chi = n_fac_chi)

}
