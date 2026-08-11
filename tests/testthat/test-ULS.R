ULS_test <- .estimate_model(method = "ULS",test_models$baseline$cormat, n_factors = 3, N = 500)
ULS_test_1 <- .estimate_model(method = "ULS",test_models$baseline$cormat, n_factors = 1, N = 500)

test_that("output class and dimensions are correct", {
  expect_s3_class(ULS_test$unrot_loadings, "LOADINGS")
  expect_output(str(ULS_test), "List of 12")
  expect_s3_class(ULS_test_1$unrot_loadings, "LOADINGS")
  expect_output(str(ULS_test_1), "List of 12")

})

test_that("outputs are correct", {
  expect_equal(ULS_test$orig_R, test_models$baseline$cormat)
  expect_equal(sum(ULS_test$orig_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(ULS_test$final_eigen), ncol(test_models$baseline$cormat))
  expect_equal(ULS_test$convergence, 0)

  expect_equal(ULS_test_1$orig_R, test_models$baseline$cormat)
  expect_equal(sum(ULS_test_1$orig_eigen), ncol(test_models$baseline$cormat))
  expect_lt(sum(ULS_test_1$final_eigen), ncol(test_models$baseline$cormat))
  expect_equal(ULS_test_1$convergence, 0)

  # The reported eigenvalues are those of the original and of the reduced correlation
  # matrix, in decreasing order. Checked against the eigenvector-returning LAPACK
  # driver, which is a different routine from the values-only one used to compute them,
  # so the tolerance is set well above their disagreement at the noise floor (~6e-15)
  # and well below anything a real regression could hide in. Both factor counts are
  # checked: the single-factor path is where the extraction degenerates.
  for (fit in list(ULS_test, ULS_test_1)) {
    R_final <- test_models$baseline$cormat
    diag(R_final) <- fit$h2
    expect_equal(fit$orig_eigen,
                 eigen(test_models$baseline$cormat, symmetric = TRUE)$values,
                 tolerance = 1e-10)
    expect_equal(fit$final_eigen, eigen(R_final, symmetric = TRUE)$values,
                 tolerance = 1e-10)
    expect_equal(fit$orig_eigen, sort(fit$orig_eigen, decreasing = TRUE))
    expect_equal(fit$final_eigen, sort(fit$final_eigen, decreasing = TRUE))
  }
})

test_that("fit indices are returned correctly", {
  expect_output(str(ULS_test$fit_indices), "List of 18")
  expect_output(str(ULS_test_1$fit_indices), "List of 18")

  expect_type(ULS_test$fit_indices$chi, "double")
  expect_type(ULS_test$fit_indices$df, "double")
  expect_type(ULS_test$fit_indices$p_chi, "double")
  expect_type(ULS_test$fit_indices$CAF, "double")
  expect_type(ULS_test$fit_indices$CFI, "double")
  expect_type(ULS_test$fit_indices$RMSEA, "double")
  expect_type(ULS_test$fit_indices$RMSEA_LB, "double")
  expect_type(ULS_test$fit_indices$RMSEA_UB, "double")
  expect_type(ULS_test$fit_indices$AIC, "double")
  expect_type(ULS_test$fit_indices$BIC, "double")
  expect_type(ULS_test$fit_indices$Fm, "double")
  expect_type(ULS_test$fit_indices$chi_null, "double")
  expect_type(ULS_test$fit_indices$df_null, "double")
  expect_type(ULS_test$fit_indices$p_null, "double")

  expect_type(ULS_test_1$fit_indices$chi, "double")
  expect_type(ULS_test_1$fit_indices$df, "double")
  expect_type(ULS_test_1$fit_indices$p_chi, "double")
  expect_type(ULS_test_1$fit_indices$CAF, "double")
  expect_type(ULS_test_1$fit_indices$CFI, "double")
  expect_type(ULS_test_1$fit_indices$RMSEA, "double")
  expect_type(ULS_test_1$fit_indices$RMSEA_LB, "double")
  expect_type(ULS_test_1$fit_indices$RMSEA_UB, "double")
  expect_type(ULS_test_1$fit_indices$AIC, "double")
  expect_type(ULS_test_1$fit_indices$BIC, "double")
  expect_type(ULS_test_1$fit_indices$Fm, "double")
  expect_type(ULS_test_1$fit_indices$chi_null, "double")
  expect_type(ULS_test_1$fit_indices$df_null, "double")
  expect_type(ULS_test_1$fit_indices$p_null, "double")
})

test_that("the ULS objective is half the squared Frobenius norm of the reduced residual", {
  # The criterion is 0.5 * ||R - diag(psi) - LL'||_F^2 with L the top-n_fac eigen
  # extraction of the reduced matrix (negative eigenvalues clipped to zero), and its
  # gradient is the negated diagonal of that residual.
  set.seed(4)
  p <- 7L; m <- 2L
  A <- matrix(stats::rnorm(p * p), p, p)
  Rg <- stats::cov2cor(A %*% t(A))
  psi <- stats::runif(p, 0.25, 0.75)

  Rs <- Rg
  diag(Rs) <- diag(Rs) - psi
  ev <- eigen(Rs, symmetric = TRUE)
  lam <- pmax(ev$values[seq_len(m)], 0)
  L <- ev$vectors[, seq_len(m), drop = FALSE] %*% diag(sqrt(lam), nrow = m)
  E <- Rs - tcrossprod(L)

  expect_equal(.uls_residuals(psi, Rg, m), 0.5 * sum(E^2))
  expect_equal(as.vector(.grad_uls(psi, Rg, m)), -diag(E))
})

test_that("the ULS gradient matches finite differences", {
  # value and gradient must differentiate the same function: checked against central
  # differences at a non-trivial point away from the optimum and off the box boundary
  set.seed(3)
  p <- 9L; m <- 3L
  A <- matrix(stats::rnorm(p * p), p, p)
  Rg <- stats::cov2cor(A %*% t(A))
  psi <- stats::runif(p, 0.25, 0.75)

  ga <- as.vector(.grad_uls(psi, Rg, m))
  h <- 1e-6
  gn <- vapply(seq_along(psi), function(k) {
    e <- numeric(length(psi)); e[k] <- h
    (.uls_residuals(psi + e, Rg, m) - .uls_residuals(psi - e, Rg, m)) / (2 * h)
  }, numeric(1))

  expect_equal(ga, gn, tolerance = 1e-5)
})

test_that("MINRES is accepted as a synonym for ULS", {
  uls <- suppressWarnings(EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                              method = "ULS"))
  minres <- suppressWarnings(EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                                 method = "MINRES"))

  # minimum residual and unweighted least squares are the same estimator
  expect_equal(minres$unrot_loadings, uls$unrot_loadings)
  expect_equal(minres$h2, uls$h2)
  expect_equal(minres$fit_indices, uls$fit_indices)
  # the alias resolves to the canonical method name
  expect_identical(minres$settings$estimator, "ULS")
})

test_that("SL() and EFA_AVERAGE() accept method = 'MINRES'", {
  EFA_mod <- suppressWarnings(EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
                                  method = "PAF", rotation = "promax"))
  expect_no_error(suppressWarnings(SL(EFA_mod, method = "MINRES")))
  expect_no_error(suppressWarnings(
    EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                method = "MINRES", show_progress = FALSE)))
})

rm(ULS_test, ULS_test_1)
