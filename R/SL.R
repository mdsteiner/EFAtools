SL <- function(L1, Phi1, cors = TRUE, criterion = .001,
                         criterion_type = "max_individual",
                         init_comm = "smc", paf_type = "norm") {

  # determine number of factors
  n_factors <- ncol(Phi1)

  # perform a factor analysis on the intercorrelation matrix of the first order
  # factors
  paf_p1 <- paf(Phi1, n_factors = 1, cors = cors, criterion = criterion,
                criterion_type = criterion_type, init_comm = init_comm,
                paf_type = paf_type)

  # extract second order loadings
  L2 <- paf_p1$loadings

  # Schmid-Leiman solution, direct loadings of second order factor
  L_sls_2 <- L1 %*% L2

  # compute uniqueness of higher order factor
  U2 <- sqrt(1 - diag(L2 %*% t(L2)))

  # Schmid-Leiman solution, residualized first order factor loadings
  L_sls_1 <- L1 %*% diag(U2)

  output <- list(
    L_sls_1 = L_sls_1,
    L_sls_2 = L_sls_2,
    L_2 = L2)

}
