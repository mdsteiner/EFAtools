# Stage-1 FIML engine: the EM estimate of the saturated MVN mean/covariance from raw data
# with missing values, and the saturated observed-data log-likelihood. Numerics are
# cross-checked against lavaan's two-stage / FIML moments (missing = "ml").

# A missing-at-random fixture: column 1 is fully observed and drives the missingness in the
# others, so the mechanism depends only on observed data.
make_mar_data <- function(n = 600, seed = 456) {
  set.seed(seed)
  p <- 5
  Sig <- 0.5 ^ abs(outer(seq_len(p), seq_len(p), "-"))     # AR(1)-type positive-definite cov
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sig)
  colnames(X) <- paste0("V", seq_len(p))
  X[X[, 1] >  0.7, 2] <- NA
  X[X[, 1] < -0.7, 3] <- NA
  X[X[, 1] >  1.0, 4] <- NA
  X[X[, 1] < -1.0, 5] <- NA
  X
}

test_that("complete data reproduces the MLE covariance and the Pearson correlation", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  set.seed(123)
  p <- 5
  Sig <- 0.5 ^ abs(outer(seq_len(p), seq_len(p), "-"))
  X <- MASS::mvrnorm(400, mu = rep(0, p), Sigma = Sig)
  colnames(X) <- paste0("V", seq_len(p))

  em <- .fiml_em_moments(X)
  n <- nrow(X)

  expect_true(em$converged)
  expect_equal(em$n, n)
  expect_equal(em$n_patterns, 1L)
  # With no missing data the EM fixed point is the ML covariance (1/n scaling) and its
  # standardisation is exactly the Pearson correlation.
  expect_equal(em$sigma, (n - 1) / n * stats::cov(X), tolerance = 1e-8)
  expect_equal(stats::cov2cor(em$sigma), stats::cor(X), tolerance = 1e-8)
})

test_that("FIML moments match lavaan two-stage (missing = 'ml') under MAR", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")

  X <- make_mar_data()
  df <- as.data.frame(X)
  em <- .fiml_em_moments(X)
  expect_true(em$converged)
  expect_gt(em$n_patterns, 1L)

  fit <- lavaan::lavCor(df, missing = "ml", output = "fit")
  ss <- lavaan::lavInspect(fit, "sampstat")
  expect_equal(em$sigma, ss$cov, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(em$mu, ss$mean, tolerance = 1e-4, ignore_attr = TRUE)

  Rml <- lavaan::lavCor(df, missing = "ml", output = "cor")
  expect_equal(stats::cov2cor(em$sigma), as.matrix(Rml),
               tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("the saturated log-likelihood matches a direct per-case computation", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- make_mar_data()
  em <- .fiml_em_moments(X)

  # Transparent reference: sum the observed-subvector MVN log-density case by case (the
  # pattern-grouped helper must reproduce it).
  ref <- 0
  for (i in seq_len(nrow(X))) {
    o <- which(!is.na(X[i, ]))
    So <- em$sigma[o, o, drop = FALSE]
    d <- X[i, o] - em$mu[o]
    ref <- ref - 0.5 * (length(o) * log(2 * pi) +
                          as.numeric(determinant(So, logarithm = TRUE)$modulus) +
                          drop(crossprod(d, solve(So, d))))
  }
  expect_equal(.fiml_loglik(X, em$mu, em$sigma), ref, tolerance = 1e-8)
})

test_that("structural guards raise classed errors", {
  # A variable pair that is never jointly observed (zero coverage).
  zero_cov <- matrix(c(1, 2, NA, NA,
                       NA, NA, 3, 4), ncol = 2)
  expect_error(.fiml_em_moments(zero_cov), class = "efa_fiml_zero_coverage")

  # A non-numeric column.
  chr <- data.frame(a = c("x", "y", "z"), b = c(1, 2, 3))
  expect_error(.fiml_em_moments(chr), class = "efa_fiml_not_numeric")

  # A variable with no observed values at all.
  no_obs <- cbind(c(1, 2, 3, 4), c(NA, NA, NA, NA))
  expect_error(.fiml_em_moments(no_obs), class = "efa_fiml_no_observed")
})

test_that("hitting max_iter warns and reports non-convergence", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- make_mar_data()
  expect_warning(em <- .fiml_em_moments(X, max_iter = 1L),
                 class = "efa_fiml_em_nonconvergence")
  expect_false(em$converged)
  expect_equal(em$iter, 1L)
})

test_that("the missing-data EM fixed point is stationary under an independent E-M step", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  # The complete-data and log-likelihood tests do not pin the missing-data E-step (the
  # conditional-mean completion and the missing-by-missing covariance correction) to an
  # independent reference without lavaan. Here an independently coded per-case E-M step
  # (using solve(), not the engine's pattern-grouped chol2inv) must reproduce the converged
  # moments: a sign error or omission in the correction would move the engine's fixed point
  # away from this correct operator's, failing the check.
  X <- make_mar_data()
  em <- .fiml_em_moments(X, tol = 1e-10)
  mu <- unname(em$mu)
  sigma <- unname(em$sigma)
  n <- em$n
  p <- ncol(X)

  T1 <- numeric(p)
  T2 <- matrix(0, p, p)
  for (i in seq_len(nrow(X))) {
    o <- which(!is.na(X[i, ]))
    m <- which(is.na(X[i, ]))
    yc <- numeric(p)
    yc[o] <- X[i, o]
    extra <- matrix(0, p, p)
    if (length(m)) {
      B <- sigma[m, o, drop = FALSE] %*% solve(sigma[o, o, drop = FALSE])
      yc[m] <- mu[m] + B %*% (X[i, o] - mu[o])
      extra[m, m] <- sigma[m, m, drop = FALSE] - B %*% sigma[o, m, drop = FALSE]
    }
    T1 <- T1 + yc
    T2 <- T2 + tcrossprod(yc) + extra
  }
  mu_ref <- T1 / n
  sigma_ref <- T2 / n - tcrossprod(mu_ref)

  expect_equal(mu, mu_ref, tolerance = 1e-7)
  expect_equal(sigma, sigma_ref, tolerance = 1e-7)
})

test_that("the FIML correlation is invariant to per-column rescaling", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  # FIML correlations must not depend on the variables' measurement scale; an absolute
  # convergence tolerance would stop early on small-variance columns and return a wrong
  # correlation. Rescaling each column must leave cov2cor(sigma) unchanged.
  X <- make_mar_data()
  em1 <- .fiml_em_moments(X)
  Xs <- sweep(X, 2L, c(1e-3, 10, 1, 1e2, 0.5), "*")
  em2 <- .fiml_em_moments(Xs)

  expect_true(em1$converged && em2$converged)
  expect_equal(stats::cov2cor(em1$sigma), stats::cov2cor(em2$sigma),
               tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("a non-positive-definite covariance aborts with a classed error", {
  set.seed(11)
  degenerate <- cbind(stats::rnorm(60), rep(3, 60), stats::rnorm(60))   # column 2 is constant
  expect_error(.fiml_em_moments(degenerate), class = "efa_fiml_not_posdef")
})

test_that("the diagonal-init fallback (too few complete cases) still converges", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  set.seed(202)
  p <- 4
  n <- 300
  Z <- MASS::mvrnorm(n, rep(0, p), 0.4 ^ abs(outer(seq_len(p), seq_len(p), "-")))
  colnames(Z) <- paste0("V", seq_len(p))
  # One NA per row (rotating column) leaves only two complete rows (< p + 1), forcing the
  # diagonal-init fallback, while every pair keeps ample joint coverage.
  for (i in seq_len(n - 2L)) Z[i, (i %% p) + 1L] <- NA

  expect_lt(sum(stats::complete.cases(Z)), p + 1L)       # the fallback path is exercised
  em <- .fiml_em_moments(Z)
  expect_true(em$converged)
  expect_true(all(eigen(em$sigma, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("fully-missing rows are dropped from n and the log-likelihood", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- make_mar_data()
  Xd <- rbind(X, NA)                                      # append an all-NA row
  em <- .fiml_em_moments(Xd)

  expect_equal(em$n, nrow(X))
  expect_lt(em$n, nrow(Xd))
  # .fiml_loglik ignores the all-NA row too, so it returns the same value with or without it.
  expect_equal(.fiml_loglik(Xd, em$mu, em$sigma),
               .fiml_loglik(X, em$mu, em$sigma), tolerance = 1e-10)
})

test_that("a non-positive or non-integer max_iter is rejected", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
  expect_error(.fiml_em_moments(X, max_iter = 0L), class = "efa_fiml_bad_max_iter")
  expect_error(.fiml_em_moments(X, max_iter = 2.5), class = "efa_fiml_bad_max_iter")
})

test_that("a non-positive tol is rejected", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
  expect_error(.fiml_em_moments(X, tol = 0), class = "efa_fiml_bad_tol")
  expect_error(.fiml_em_moments(X, tol = -1), class = "efa_fiml_bad_tol")
})
