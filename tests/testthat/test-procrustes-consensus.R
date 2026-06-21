# Property tests for the internal GPA-consensus engine that backs
# EFA_POOLED(target_method = "consensus"): symmetry across the m inputs,
# monotone descent under the standard Gower 1975 GPA update, agreement with
# a hand-coded GPA reference loop, stationary-point convergence, and the
# oblique-rotation guard.

gpa_consensus_target <- EFAtools:::.gpa_consensus_target

# Synthetic m-list around a simple-structure L_true with random per-imputation
# orthogonal unrotations + small noise. Inputs are kept small to minimise
# test runtime.
.make_imp_list <- function(m, p, k, sigma = 0, seed = 1L) {
  set.seed(seed)
  L_true <- matrix(0, p, k)
  per <- p %/% k
  for (j in seq_len(k)) {
    L_true[((j - 1L) * per + 1L):(j * per), j] <- runif(per, 0.6, 0.8)
  }
  unrot <- lapply(seq_len(m), function(d) {
    L_true %*% qr.Q(qr(matrix(rnorm(k * k), k, k))) +
      matrix(rnorm(p * k, sd = sigma), p, k)
  })
  list(L_true = L_true, unrot = unrot)
}


test_that(".gpa_consensus_target is symmetric in its m inputs (orthogonal)", {
  d <- .make_imp_list(m = 3L, p = 8L, k = 2L, sigma = 0.02, seed = 11)
  fit1 <- gpa_consensus_target(
    unrotated_list = d$unrot,
    init_targets = list(d$L_true), start = 1L,
    rotation = "orthogonal",
    tol = 1e-8, loss_tol = 1e-10, max_iter = 100, min_iter = 2L
  )
  perm <- c(3L, 1L, 2L)
  fit2 <- gpa_consensus_target(
    unrotated_list = d$unrot[perm],
    init_targets = list(d$L_true), start = 1L,
    rotation = "orthogonal",
    tol = 1e-8, loss_tol = 1e-10, max_iter = 100, min_iter = 2L
  )
  expect_equal(unname(fit1$target), unname(fit2$target), tolerance = 1e-7)
})


test_that("orthogonal consensus loss is monotone non-increasing (Gower 1975 GPA)", {
  # match_target = FALSE, alpha = 1, rotation = orthogonal: block-coordinate
  # GPA, guaranteed monotone non-increase of the mean squared discrepancy.
  d <- .make_imp_list(m = 4L, p = 9L, k = 3L, sigma = 0.05, seed = 23)
  fit <- gpa_consensus_target(
    unrotated_list = d$unrot,
    rotation = "orthogonal",
    match_target = FALSE, alpha = 1,
    tol = 1e-10, loss_tol = 1e-12, max_iter = 100, min_iter = 2L
  )
  expect_true(all(diff(fit$history$loss) <= 1e-12))
})


test_that("consensus target is a fixed point of one more centroid update", {
  d <- .make_imp_list(m = 3L, p = 8L, k = 2L, sigma = 0.03, seed = 31)
  tol <- 1e-8
  fit <- gpa_consensus_target(
    unrotated_list = d$unrot,
    rotation = "orthogonal",
    tol = tol, loss_tol = 1e-12, max_iter = 200, min_iter = 2L
  )
  next_aligned <- lapply(d$unrot, function(L) {
    PROCRUSTES(L, fit$target, rotation = "orthogonal")$loadings
  })
  centroid_next <- Reduce(`+`, next_aligned) / length(next_aligned)
  rel <- norm(unname(centroid_next) - unname(fit$target), type = "F") /
    (norm(unname(fit$target), type = "F") + 1e-12)
  expect_lt(rel, 100 * tol)
})


test_that("orthogonal consensus matches a hand-coded Gower 1975 GPA loop (m = 2)", {
  # Independent reference: alternate orthogonal Procrustes per d, then
  # arithmetic-mean centroid. Same limit at the convergence tolerance.
  d <- .make_imp_list(m = 2L, p = 8L, k = 2L, sigma = 0.05, seed = 37)
  fit <- gpa_consensus_target(
    unrotated_list = d$unrot,
    init_targets = list(d$unrot[[1L]]), start = 1L,
    rotation = "orthogonal",
    match_target = FALSE, alpha = 1,
    tol = 1e-12, loss_tol = 1e-14, max_iter = 200, min_iter = 2L
  )
  target <- d$unrot[[1L]]
  for (it in seq_len(200L)) {
    aligned <- lapply(d$unrot, function(L) {
      s <- svd(crossprod(L, target))
      L %*% (s$u %*% t(s$v))
    })
    new_target <- Reduce(`+`, aligned) / length(aligned)
    rel <- norm(new_target - target, type = "F") /
      (norm(target, type = "F") + 1e-12)
    target <- new_target
    if (rel < 1e-12) break
  }
  expect_equal(unname(fit$target), unname(target), tolerance = 1e-8)
})


test_that("oblique GPA-consensus aborts with a classed condition for k > 1", {
  # The naive iterated-oblique-Procrustes design has degenerate fixed points
  # for k > 1; Lorenzo-Seva & Van Ginkel (2016) instead apply a Promin step
  # on top of the centroid, which this engine does not implement.
  d <- .make_imp_list(m = 2L, p = 6L, k = 2L, sigma = 0.02, seed = 43)
  expect_error(
    gpa_consensus_target(
      unrotated_list = d$unrot,
      init_targets = list(d$L_true), start = 1L,
      rotation = "oblique"
    ),
    class = "efa_consensus_oblique_unsupported"
  )
})
