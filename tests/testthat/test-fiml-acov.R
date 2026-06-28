# Asymptotic covariance of the Stage-1 saturated FIML estimates: the inverse observed
# information (negative Hessian) of theta = (mu, vech(sigma)) and its correlation-scale variant.
# The analytic score is cross-checked against a numerical gradient of the saturated
# log-likelihood, the theta-scale covariance against lavaan's saturated-model vcov, and the
# correlation-scale variant against a finite-difference standardisation and a Monte-Carlo
# sampling covariance.

# MAR fixture: column 1 is fully observed and drives the missingness in the others, so the
# mechanism depends only on observed data (mirrors the Stage-1 engine tests).
.acov_mar_data <- function(n = 500, seed = 321) {
  set.seed(seed)
  p <- 4
  Sig <- 0.5 ^ abs(outer(seq_len(p), seq_len(p), "-"))     # AR(1)-type positive-definite cov
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sig)
  colnames(X) <- paste0("V", seq_len(p))
  X[X[, 1] >  0.8, 2] <- NA
  X[X[, 1] < -0.8, 3] <- NA
  X[X[, 1] >  1.1, 4] <- NA
  X
}

test_that("the analytic score matches a numerical gradient of the log-likelihood", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- .acov_mar_data()
  em <- .fiml_em_moments(X)
  p <- ncol(X)
  pstar <- p * (p + 1L) / 2L

  loglik <- function(theta) {
    sig <- .efa_pooled_unvech(theta[p + seq_len(pstar)], p)
    .fiml_loglik(X, theta[seq_len(p)], sig)
  }

  # Evaluate away from the EM fixed point (where the score is ~0), keeping sigma posdef.
  set.seed(9)
  theta <- c(em$mu, .efa_pooled_vech(em$sigma)) + 0.02 * stats::rnorm(p + pstar)
  sig <- .efa_pooled_unvech(theta[p + seq_len(pstar)], p)
  skip_if(any(eigen(sig, symmetric = TRUE, only.values = TRUE)$values < 1e-6))

  ana <- .fiml_saturated_score(X, theta[seq_len(p)], sig)

  h <- 1e-5
  num <- vapply(seq_along(theta), function(k) {
    tp <- theta; tp[k] <- tp[k] + h
    tm <- theta; tm[k] <- tm[k] - h
    (loglik(tp) - loglik(tm)) / (2 * h)
  }, numeric(1))

  expect_equal(ana, num, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("the score vanishes at the EM fixed point", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- .acov_mar_data()
  em <- .fiml_em_moments(X, tol = 1e-10, max_iter = 1000L)
  sc <- .fiml_saturated_score(X, em$mu, em$sigma)
  # The score is O(n); a near-zero value relative to n confirms the fixed point is the MLE.
  expect_lt(max(abs(sc)) / em$n, 1e-5)
})

test_that("the information and asymptotic covariances are well-formed", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- .acov_mar_data()
  em <- .fiml_em_moments(X)
  p <- ncol(X)
  pstar <- p * (p + 1L) / 2L
  npair <- p * (p - 1L) / 2L

  ac <- .fiml_saturated_acov(X, em$mu, em$sigma)

  # theta-scale: square (p + p*), symmetric, positive definite, labelled.
  expect_equal(dim(ac$theta), c(p + pstar, p + pstar))
  expect_equal(rownames(ac$theta)[seq_len(p)], colnames(X))
  expect_true(isSymmetric(unname(ac$theta), tol = 1e-8))
  expect_gt(min(eigen(ac$theta, symmetric = TRUE, only.values = TRUE)$values), 0)

  # cor-scale: off-diagonal pairs, combn order/labels (same layout as .adf_gamma / .pair_labels).
  expect_equal(dim(ac$cor), c(npair, npair))
  expect_equal(rownames(ac$cor), .pair_labels(colnames(X)))
  expect_true(isSymmetric(unname(ac$cor), tol = 1e-8))
  expect_gt(min(eigen(ac$cor, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
  # each correlation carries strictly positive sampling variance (guards against a vacuously
  # PSD all-zero block from a mis-indexed Jacobian or sigma-block extraction).
  expect_true(all(diag(ac$cor) > 0))
})

test_that("a non-positive-definite covariance aborts with a classed error", {
  # A deterministic indefinite covariance (one negative eigenvalue); both the score and the ACOV
  # must refuse with the shared classed condition rather than leaking a bare chol() error.
  X <- matrix(c(1, 2, 3, 2, 1, 4, 3, 5, 2, 4, 1, 6), ncol = 3)
  colnames(X) <- paste0("V", seq_len(3))
  mu <- colMeans(X)
  bad <- matrix(c(1, 0.99, 0.99, 0.99, 1, -0.99, 0.99, -0.99, 1), 3)
  expect_lt(min(eigen(bad, symmetric = TRUE, only.values = TRUE)$values), 0)

  expect_error(.fiml_saturated_score(X, mu, bad), class = "efa_fiml_not_posdef")
  expect_error(.fiml_saturated_acov(X, mu, bad), class = "efa_fiml_not_posdef")
})

test_that("the theta-scale asymptotic covariance matches lavaan's saturated vcov", {
  skip_on_cran()
  skip_if_not_installed("lavaan")
  skip_if_not_installed("MASS")

  X <- .acov_mar_data()
  df <- as.data.frame(X)
  p <- ncol(X)
  vn <- colnames(X)
  em <- .fiml_em_moments(X)
  ours <- .fiml_saturated_acov(X, em$mu, em$sigma)$theta

  # Explicit saturated model: free every variance/covariance and mean, fit by FIML with the
  # observed information so the parameter vcov is the inverse observed information we compute.
  ij <- which(lower.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  cov_lines <- sprintf("%s ~~ %s", vn[ij[, 2L]], vn[ij[, 1L]])
  mean_lines <- sprintf("%s ~ 1", vn)
  mod <- paste(c(cov_lines, mean_lines), collapse = "\n")
  fit <- suppressWarnings(
    lavaan::lavaan(mod, data = df, missing = "ml",
                   information = "observed", se = "standard"))
  V <- lavaan::vcov(fit)

  # Map lavaan's parameters to our theta order by canonicalised labels (sorted endpoints for a
  # covariance; the orientation lavaan stores a "~~" pair in is not guaranteed).
  canon <- function(lab) {
    lab <- gsub("\\s", "", lab)
    if (grepl("~~", lab, fixed = TRUE)) {
      vv <- sort(strsplit(lab, "~~", fixed = TRUE)[[1]])
      paste0("cov:", vv[1L], ":", vv[2L])
    } else {
      paste0("mean:", sub("~1$", "", lab))
    }
  }
  lav_keys <- vapply(rownames(V), canon, character(1))
  cov_keys <- vapply(seq_len(nrow(ij)), function(t) {
    vv <- sort(c(vn[ij[t, 1L]], vn[ij[t, 2L]]))
    paste0("cov:", vv[1L], ":", vv[2L])
  }, character(1))
  our_keys <- c(paste0("mean:", vn), cov_keys)

  ord <- match(our_keys, lav_keys)
  expect_false(anyNA(ord))
  lav <- V[ord, ord, drop = FALSE]

  # Split tolerance as the polychoric ACOV tests do: variances relatively, the tiny off-diagonal
  # covariances absolutely (scaled to the variance magnitude). The residual is our EM vs
  # lavaan's optimiser plus the finite-difference Jacobian; measured agreement is far tighter.
  expect_lt(max(abs(diag(ours) - diag(lav)) / diag(lav)), 0.05)
  off <- upper.tri(ours)
  expect_lt(max(abs(ours[off] - lav[off])), 0.05 * max(diag(lav)))
})

test_that("the correlation-scale ACOV equals a finite-difference standardisation", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  X <- .acov_mar_data()
  em <- .fiml_em_moments(X)
  p <- ncol(X)
  pstar <- p * (p + 1L) / 2L
  ac <- .fiml_saturated_acov(X, em$mu, em$sigma)

  # Off-diagonal correlations of cov2cor(sigma), in combn order, as a function of vech(sigma).
  pairs <- utils::combn(p, 2L)
  cor_off <- function(vs) {
    R <- stats::cov2cor(.efa_pooled_unvech(vs, p))
    R[cbind(pairs[1L, ], pairs[2L, ])]
  }
  vs0 <- .efa_pooled_vech(em$sigma)
  h <- 1e-6
  Jn <- vapply(seq_len(pstar), function(k) {
    vp <- vs0; vp[k] <- vp[k] + h
    vm <- vs0; vm[k] <- vm[k] - h
    (cor_off(vp) - cor_off(vm)) / (2 * h)
  }, numeric(ncol(pairs)))

  sig_idx <- p + seq_len(pstar)
  Omega_sig <- ac$theta[sig_idx, sig_idx]
  Omega_R_fd <- Jn %*% Omega_sig %*% t(Jn)

  expect_equal(unname(ac$cor), unname(Omega_R_fd), tolerance = 1e-5)
})

test_that("the correlation-scale ACOV tracks the Monte-Carlo sampling covariance", {
  skip_on_cran()
  skip_if_not_slow()
  skip_if_not_installed("MASS")

  # The analytic asymptotic variance of each off-diagonal FIML correlation must track the
  # empirical sampling variance over many fresh MAR samples from the same population. The
  # analytic variance is averaged over the samples to remove single-sample noise (it is itself
  # estimated from data), so the check isolates the method's calibration; the tolerance then
  # covers the Monte-Carlo error in the empirical variances.
  p <- 4
  pairs <- utils::combn(p, 2L)
  R <- 1000L
  rmat <- matrix(NA_real_, R, ncol(pairs))
  amat <- matrix(NA_real_, R, ncol(pairs))
  for (r in seq_len(R)) {
    Xr <- .acov_mar_data(seed = 1000L + r)
    emr <- .fiml_em_moments(Xr)
    rmat[r, ] <- stats::cov2cor(emr$sigma)[cbind(pairs[1L, ], pairs[2L, ])]
    amat[r, ] <- diag(.fiml_saturated_acov(Xr, emr$mu, emr$sigma)$cor)
  }
  empirical <- apply(rmat, 2L, stats::var)
  analytic <- colMeans(amat)

  expect_lt(max(abs(analytic / empirical - 1)), 0.15)
})
