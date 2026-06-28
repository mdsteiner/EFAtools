# Stage-1 of two-stage / full-information maximum likelihood (FIML) correlation estimation:
# from raw data with missing values, EM-estimate the saturated multivariate-normal mean and
# covariance assuming multivariate normality and ignorable (MAR) missingness, plus the
# saturated observed-data log-likelihood that the fit-index and standard-error paths reuse.
# Two-stage EFA with missing data: Yuan, Marshall & Bentler (2002, Psychometrika 67:95-121).
# MAR/ignorability: Little & Rubin (2002). EM reference implementation: lavaan's
# lav_mvn_mi_h1_est_moments() / lav_mvnorm_missing_*().
#
# The asymptotic covariance of the saturated estimates (the meat for corrected two-stage
# standard errors) is the inverse observed information of the saturated log-likelihood. Under
# FIML the observed -- not the expected -- information is used: the expected information would
# have to be averaged over the random missingness distribution and is inconsistent under MAR
# (Yuan & Bentler, 2000; Savalei & Bentler, 2009, SEM). lavaan/Mplus make the same choice.

# Group row indices by missingness pattern (unique sets of observed columns); returns, per
# pattern, the observed/missing column indices and the rows that share it, so the EM and the
# log-likelihood can accumulate per pattern rather than per case.
.fiml_patterns <- function(obs) {
  # One key per row via a vectorised column-wise paste of the 0/1 mask (avoids a per-row
  # apply() and stays exact for any number of variables).
  keys <- do.call(paste0, asplit(matrix(as.integer(obs), nrow(obs)), 2L))
  lapply(split(seq_len(nrow(obs)), keys), function(rows) {
    pat <- obs[rows[1L], ]
    list(rows = rows, o = which(pat), m = which(!pat), freq = length(rows))
  })
}

# EM estimate of the saturated MVN mean and covariance from raw data with missing values
# (the FIML / two-stage Stage-1 engine; Yuan, Marshall & Bentler, 2002).
.fiml_em_moments <- function(data, max_iter = 500L, tol = 1e-5) {

  if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1 ||
      max_iter != as.integer(max_iter)) {
    cli::cli_abort("{.arg max_iter} must be a single positive integer.",
                   class = "efa_fiml_bad_max_iter")
  }

  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    cli::cli_abort("{.arg tol} must be a single positive number.",
                   class = "efa_fiml_bad_tol")
  }

  Y <- as.matrix(data)
  if (!is.numeric(Y)) {
    cli::cli_abort(
      c("{.arg data} must be numeric to estimate FIML moments.",
        "x" = "You supplied {.obj_type_friendly {data}}."),
      class = "efa_fiml_not_numeric"
    )
  }

  p <- ncol(Y)
  nms <- colnames(Y)
  if (is.null(nms)) nms <- paste0("V", seq_len(p))

  obs <- !is.na(Y)

  # Fully-missing rows carry no information for the saturated moments; drop them so n counts
  # the cases actually used.
  keep <- rowSums(obs) > 0L
  Y <- Y[keep, , drop = FALSE]
  obs <- obs[keep, , drop = FALSE]
  n <- nrow(Y)

  # Coverage: diagonal = observed cases per variable, off-diagonal = jointly-observed cases
  # per pair. A variable or a pair with zero coverage cannot be estimated.
  coverage <- crossprod(obs)
  if (any(diag(coverage) == 0)) {
    empty <- nms[diag(coverage) == 0]
    n_empty <- length(empty)
    cli::cli_abort(
      c("Every variable must have at least one observed value.",
        "x" = "{n_empty} variable{?s} ({.val {empty}}) {?has/have} no observed values."),
      class = "efa_fiml_no_observed"
    )
  }
  zero_pairs <- which(coverage == 0 & upper.tri(coverage), arr.ind = TRUE)
  if (nrow(zero_pairs) > 0L) {
    pair <- paste(nms[zero_pairs[1L, 1L]], nms[zero_pairs[1L, 2L]], sep = " & ")
    cli::cli_abort(
      c("Every variable pair must have at least one jointly observed case.",
        "x" = "Variables {.val {pair}} are never observed together (zero coverage)."),
      class = "efa_fiml_zero_coverage"
    )
  }

  patterns <- .fiml_patterns(obs)
  n_patterns <- length(patterns)

  # Init: mu from the available data; sigma from the complete cases when there are enough of
  # them and their covariance is well-conditioned, else a diagonal of column variances
  # (floored to stay positive definite). A near-singular complete-case covariance would seed
  # the EM with an ill-conditioned start, so it is rejected by an rcond threshold rather than
  # by chol() alone (which accepts a tiny positive pivot).
  mu <- colMeans(Y, na.rm = TRUE)
  cc <- stats::complete.cases(Y)
  sigma <- NULL
  if (sum(cc) > p) {
    s_cc <- stats::cov(Y[cc, , drop = FALSE])
    if (!inherits(try(chol(s_cc), silent = TRUE), "try-error") &&
        rcond(s_cc) > .Machine$double.eps) {
      sigma <- s_cc
    }
  }
  if (is.null(sigma)) {
    v <- apply(Y, 2L, stats::var, na.rm = TRUE)
    v[!is.finite(v) | v <= 0] <- 1
    sigma <- diag(v, nrow = p)
  }

  # Convergence is judged on the scale-invariant standardised moments (the standardised means
  # mu/sd, the log-variances, and the off-diagonal correlations), so the criterion does not
  # depend on the variables' measurement scale; an absolute tolerance on the raw moments would
  # stop early on small-variance data and over-iterate on large-variance data. lavaan's
  # lav_mvnorm_missing_* uses a likelihood-based criterion for the same scale-invariance.
  lt_off <- lower.tri(sigma)
  std_par <- function(mu, sigma) {
    v <- diag(sigma)
    s <- sqrt(v)
    c(mu / s, log(v), (sigma / tcrossprod(s))[lt_off])
  }
  par_old <- std_par(mu, sigma)

  converged <- FALSE
  for (iter in seq_len(max_iter)) {

    # E-step accumulators: T1 = sum of completed rows, T2 = sum of their outer products.
    T1 <- numeric(p)
    T2 <- matrix(0, p, p)

    for (pat in patterns) {
      o <- pat$o
      m <- pat$m
      Yo <- Y[pat$rows, o, drop = FALSE]

      if (length(m) == 0L) {
        # Complete pattern: the observed block contributes directly.
        T1[o] <- T1[o] + colSums(Yo)
        T2[o, o] <- T2[o, o] + crossprod(Yo)
        next
      }

      # Solve each pattern's Sigma_oo directly (lavaan instead inverts the full Sigma once and
      # downdates each pattern from it).
      Soo_inv <- chol2inv(chol(sigma[o, o, drop = FALSE]))
      B <- sigma[m, o, drop = FALSE] %*% Soo_inv          # conditional-regression block

      # Complete each row: observed entries as-is, missing entries -> conditional means
      # yhat_m = mu_m + B (y_o - mu_o).
      Yc <- matrix(0, length(pat$rows), p)
      Yc[, o] <- Yo
      Yc[, m] <- sweep(sweep(Yo, 2L, mu[o], "-") %*% t(B), 2L, mu[m], "+")

      T1 <- T1 + colSums(Yc)
      T2 <- T2 + crossprod(Yc)

      # Conditional-covariance correction for the missing x missing block, added once per case
      # (constant within the pattern). Omitting it biases sigma downward.
      corr <- sigma[m, m, drop = FALSE] - B %*% sigma[o, m, drop = FALSE]
      T2[m, m] <- T2[m, m] + pat$freq * corr
    }

    # M-step.
    mu <- T1 / n
    sigma <- T2 / n - tcrossprod(mu)
    sigma <- (sigma + t(sigma)) / 2                       # symmetrise away round-off

    # A loss of positive-definiteness signals a (near-)degenerate / collinear problem.
    if (inherits(try(chol(sigma), silent = TRUE), "try-error")) {
      cli::cli_abort(
        c("The FIML covariance estimate stopped being positive definite during the EM iteration.",
          "i" = "A variable may be (near-)constant or collinear; check the data."),
        class = "efa_fiml_not_posdef"
      )
    }

    par_new <- std_par(mu, sigma)
    if (max(abs(par_new - par_old)) < tol) {
      converged <- TRUE
      break
    }
    par_old <- par_new
  }

  if (!converged) {
    cli::cli_warn(
      c("FIML EM did not converge in {max_iter} iteration{?s}.",
        "i" = "Increase {.arg max_iter} or relax {.arg tol}; the last iterate is returned."),
      class = "efa_fiml_em_nonconvergence"
    )
  }

  names(mu) <- nms
  dimnames(sigma) <- list(nms, nms)

  # The saturated log-likelihood; this call also re-validates positive-definiteness, so a
  # near-singular fixed point that slipped past the per-iteration chol() guard (chol accepts a
  # sub-eps eigenvalue) aborts here with efa_fiml_not_posdef rather than being returned.
  logl <- .fiml_loglik(Y, mu, sigma, patterns = patterns)

  list(mu = mu, sigma = sigma, logl = logl, n = n,
       n_patterns = n_patterns, iter = iter, converged = converged)
}

# Saturated (h1) observed-data multivariate-normal log-likelihood for data with missing
# values: the sum over cases of the observed-variable MVN log-density, grouped by missingness
# pattern. Reused by the FIML fit-index and standard-error paths.
.fiml_loglik <- function(data, mu, sigma, patterns = NULL) {

  # A non-positive-definite sigma has no valid MVN density; refuse with a classed error rather
  # than letting chol() raise a bare base error (this helper is reused with externally supplied
  # moments). A principal submatrix of a positive-definite matrix is itself positive definite,
  # so this single check covers Sigma_oo for every pattern. The threshold matches the
  # non-positive-definite guard in .prepare_cor_input().
  if (any(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values < .Machine$double.eps)) {
    cli::cli_abort(
      c("The covariance matrix is not positive definite, so the FIML log-likelihood is undefined.",
        "i" = "A variable may be (near-)constant or collinear."),
      class = "efa_fiml_not_posdef"
    )
  }

  Y <- as.matrix(data)

  # The EM engine already grouped the rows it passes (over its fully-missing-row-filtered data),
  # so accept those patterns directly; otherwise group here so the helper stays self-contained
  # for external callers that pass raw, unfiltered data (fully-missing rows carry no density).
  if (is.null(patterns)) {
    obs <- !is.na(Y)
    keep <- rowSums(obs) > 0L
    Y <- Y[keep, , drop = FALSE]
    patterns <- .fiml_patterns(obs[keep, , drop = FALSE])
  }

  log2pi <- log(2 * pi)
  ll <- 0
  for (pat in patterns) {
    o <- pat$o
    Yo <- Y[pat$rows, o, drop = FALSE]
    ch <- chol(sigma[o, o, drop = FALSE])
    logdet <- 2 * sum(log(diag(ch)))                     # log|Sigma_oo|
    C <- sweep(Yo, 2L, mu[o], "-")
    quad <- rowSums((C %*% chol2inv(ch)) * C)            # (y_o - mu_o)' Sigma_oo^{-1} (.)
    ll <- ll + sum(-0.5 * (length(o) * log2pi + logdet + quad))
  }
  ll
}

# Analytic gradient (score) of the saturated observed-data log-likelihood with respect to
# theta = (mu, vech(sigma)), accumulated per missingness pattern; vech is the column-major
# lower triangle including the diagonal (the .efa_pooled_vech() convention). The covariance
# block is returned with respect to the distinct free parameters: an off-diagonal sigma_ij
# appears in both Sigma[i, j] and Sigma[j, i], so its derivative carries the factor 2 from the
# duplication-matrix chain rule. Differentiating this score gives the observed information.
.fiml_saturated_score <- function(data, mu, sigma, patterns = NULL) {

  # No valid score without a positive-definite covariance (mirrors .fiml_loglik() /
  # .fiml_saturated_acov(); a principal submatrix of a positive-definite matrix is positive
  # definite, so every pattern's Sigma_oo is then safe to factor below). This also classifies the
  # rare case where a finite-difference perturbation in .fiml_saturated_acov() tips a near-
  # singular Sigma indefinite, so that path raises the classed condition rather than a bare chol().
  if (any(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values < .Machine$double.eps)) {
    cli::cli_abort(
      c("The covariance matrix is not positive definite, so the FIML score is undefined.",
        "i" = "A variable may be (near-)constant or collinear."),
      class = "efa_fiml_not_posdef"
    )
  }

  Y <- as.matrix(data)
  p <- ncol(Y)

  # The EM passes its (fully-missing-row-filtered) patterns; otherwise group here, matching
  # .fiml_loglik() so an external caller can pass raw, unfiltered data.
  if (is.null(patterns)) {
    obs <- !is.na(Y)
    keep <- rowSums(obs) > 0L
    Y <- Y[keep, , drop = FALSE]
    patterns <- .fiml_patterns(obs[keep, , drop = FALSE])
  }

  mu_grad <- numeric(p)
  Gfull <- matrix(0, p, p)                                # symmetric-matrix gradient w.r.t. Sigma

  for (pat in patterns) {
    o <- pat$o
    C <- sweep(Y[pat$rows, o, drop = FALSE], 2L, mu[o], "-")
    Soo_inv <- chol2inv(chol(sigma[o, o, drop = FALSE]))
    S1 <- colSums(C)                                      # sum of (y_o - mu_o)
    S2 <- crossprod(C)                                    # sum of (y_o - mu_o)(y_o - mu_o)'

    mu_grad[o] <- mu_grad[o] + as.numeric(Soo_inv %*% S1)
    Gfull[o, o] <- Gfull[o, o] -
      0.5 * (pat$freq * Soo_inv - Soo_inv %*% S2 %*% Soo_inv)
  }

  # Map the symmetric matrix gradient to the distinct vech parameters: double the off-diagonals
  # (chain rule), leave the diagonal as is.
  Gpar <- 2 * Gfull
  diag(Gpar) <- diag(Gfull)
  c(mu_grad, .efa_pooled_vech(Gpar))
}

# Asymptotic covariance of the Stage-1 saturated FIML estimates: the inverse observed
# information (negative Hessian of the saturated log-likelihood) of theta = (mu, vech(sigma)),
# and its correlation-scale variant for the off-diagonal correlations of cov2cor(sigma). The
# information is the negative numerical Jacobian of the analytic score .fiml_saturated_score();
# the correlation variant standardises the covariance block by the delta method. Returns a list
# with `theta` (the full ACOV on the (mu, vech(sigma)) scale) and `cor` (the off-diagonal
# correlation ACOV in utils::combn(p, 2) order, the meat the corrected two-stage standard errors
# reuse, on the same scale and layout as the continuous ADF covariance .adf_gamma()).
.fiml_saturated_acov <- function(data, mu, sigma, patterns = NULL) {

  # No valid information without a positive-definite covariance (mirrors .fiml_loglik(); a
  # principal submatrix of a positive-definite matrix is positive definite, so every pattern's
  # Sigma_oo is then safe).
  if (any(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values < .Machine$double.eps)) {
    cli::cli_abort(
      c("The covariance matrix is not positive definite, so the FIML information is undefined.",
        "i" = "A variable may be (near-)constant or collinear."),
      class = "efa_fiml_not_posdef"
    )
  }

  Y <- as.matrix(data)
  p <- ncol(Y)
  nms <- colnames(Y)
  if (is.null(nms)) nms <- names(mu)
  if (is.null(nms)) nms <- paste0("V", seq_len(p))

  # Group the rows once so every score evaluation in the finite-difference loop reuses the same
  # missingness patterns (and the same fully-missing-row filtering).
  if (is.null(patterns)) {
    obs <- !is.na(Y)
    keep <- rowSums(obs) > 0L
    Y <- Y[keep, , drop = FALSE]
    patterns <- .fiml_patterns(obs[keep, , drop = FALSE])
  }

  pstar <- p * (p + 1L) / 2L
  theta0 <- c(mu, .efa_pooled_vech(sigma))
  score_at <- function(theta) {
    sig <- .efa_pooled_unvech(theta[p + seq_len(pstar)], p)
    .fiml_saturated_score(Y, theta[seq_len(p)], sig, patterns = patterns)
  }

  # Index pairs (row >= col) of vech(sigma) in column-major order, shared below by the finite-
  # difference step scale, the theta labels, and the standardisation Jacobian.
  sds <- sqrt(diag(sigma))
  ij <- which(lower.tri(diag(p), diag = TRUE), arr.ind = TRUE)

  # Observed information = -d(score)/d(theta) by central differences. A per-parameter relative
  # step (eps^(1/3) is the central-difference optimum) keeps the Jacobian accurate across the
  # differing scales of the means, variances, and covariances; each parameter's own scale (the
  # variable SD for a mean, sd_i sd_j for sigma_ij) floors a near-zero parameter so the step never
  # overshoots its variable's spread (which could otherwise drive a small variance negative).
  d <- p + pstar
  scale_k <- c(sds, sds[ij[, 1L]] * sds[ij[, 2L]])
  h <- .Machine$double.eps^(1 / 3) * pmax(abs(theta0), scale_k)

  info <- matrix(0, d, d)
  for (k in seq_len(d)) {
    tp <- theta0; tp[k] <- tp[k] + h[k]
    tm <- theta0; tm[k] <- tm[k] - h[k]
    info[, k] <- -(score_at(tp) - score_at(tm)) / (2 * h[k])
  }
  info <- (info + t(info)) / 2

  # A singular information signals genuine near-degeneracy (coverage is validated upstream in
  # .fiml_em_moments()); refuse with its own classed condition -- distinct from the non-positive-
  # definite-covariance abort above, which is a different failure -- rather than a bare solve()
  # error.
  ch <- try(chol(info), silent = TRUE)
  if (inherits(ch, "try-error") || rcond(info) < .Machine$double.eps) {
    cli::cli_abort(
      c("The FIML information matrix is singular, so the saturated asymptotic covariance is undefined.",
        "i" = "A variable or variable pair may be (near-)degenerate or have very low coverage."),
      class = "efa_fiml_singular_information"
    )
  }
  Omega <- chol2inv(ch)
  Omega <- (Omega + t(Omega)) / 2

  # Label theta: means by variable name, vech(sigma) by "name_i-name_j" in the column-major
  # lower-triangle order (ij computed above).
  theta_lab <- c(nms, paste(nms[ij[, 1L]], nms[ij[, 2L]], sep = "-"))
  dimnames(Omega) <- list(theta_lab, theta_lab)

  # Correlation-scale variant. Take the vech(sigma) block of Omega (the MARGINAL covariance of
  # the covariance estimates -- not the inverse of the information's sigma-block -- so the
  # missingness-induced mean/covariance dependence is carried through), then standardise it by
  # the delta method to the off-diagonal correlations r_ij = sigma_ij / sqrt(sigma_ii sigma_jj).
  sig_idx <- p + seq_len(pstar)
  Omega_sig <- Omega[sig_idx, sig_idx, drop = FALSE]

  # Standardisation Jacobian J = d r_ij / d vech(sigma): rows are the off-diagonal pairs in
  # utils::combn() order, columns the vech(sigma) positions. vech_pos maps an element (row, col),
  # row >= col, to its column-major vech index (ij computed above). The two variance terms are the
  # delta-method correction (the continuous analogue of .adf_gamma()'s -1/2 r (w_i^2 + w_j^2)).
  vech_pos <- matrix(0L, p, p)
  vech_pos[ij] <- seq_len(pstar)
  pairs <- utils::combn(p, 2L)
  J <- matrix(0, ncol(pairs), pstar)
  for (t in seq_len(ncol(pairs))) {
    i <- pairs[1L, t]; j <- pairs[2L, t]                 # i < j
    rij <- sigma[i, j] / (sds[i] * sds[j])
    J[t, vech_pos[j, i]] <- 1 / (sds[i] * sds[j])        # d r / d sigma_ij  (lower-tri (j, i))
    J[t, vech_pos[i, i]] <- -0.5 * rij / sigma[i, i]
    J[t, vech_pos[j, j]] <- -0.5 * rij / sigma[j, j]
  }
  Omega_R <- J %*% Omega_sig %*% t(J)
  Omega_R <- (Omega_R + t(Omega_R)) / 2
  pair_lab <- .pair_labels(nms)
  dimnames(Omega_R) <- list(pair_lab, pair_lab)

  list(theta = Omega, cor = Omega_R)
}
