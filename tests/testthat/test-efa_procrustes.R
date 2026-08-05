
# EFAs of ten similar datasets (bootstrap samples of DOSPERT_raw). Only the first two are
# read by the blocks that run everywhere; samples 3..10 feed the warm-start comparison
# alone, which is skip_on_cran(). Each is a six-factor ML fit with promax on 3,123 x 30
# rows, so they are fitted on first request rather than up front.
#
# The resample row indices are drawn eagerly, in the original order and under the original
# seed, and only the fits are deferred. That keeps every fixture byte-identical to the
# eager version: an ML extraction with promax draws no random numbers, so interleaving the
# draws with the fits never mattered (verified against the eager construction).
set.seed(42)
boot_rows <- lapply(1:10, function(x) sample(1:nrow(DOSPERT_raw), replace = TRUE))

efa_at <- local({
  cache <- vector("list", 10)
  function(d) {
    if (is.null(cache[[d]])) {
      cache[[d]] <<- EFA(DOSPERT_raw[boot_rows[[d]], ], n_factors = 6, method = "ML",
                         rotation = "promax")
    }
    cache[[d]]
  }
})

unrot_loadings <- function(d) efa_at(d)$unrot_loadings
rot_loadings <- function(d) efa_at(d)$rot_loadings


test_that("efa_procrustes matches psych::Procrustes and GPArotation::targetQ outputs", {
  expect_equal(
    psych::Procrustes(unrot_loadings(2), rot_loadings(1))$loadings,
    efa_procrustes(unrot_loadings(2), rot_loadings(1))$loadings,
    ignore_attr = TRUE
    )
  # The default single start of GPArotation::targetQ() is itself trapped in a
  # local minimum on this matrix; random starts are needed for it to reach the
  # global oblique target-rotation optimum that efa_procrustes() now finds.
  skip_if_not_installed("GPArotation")
  set.seed(42)
  oracle <- suppressWarnings(
    GPArotation::targetQ(unrot_loadings(2), Target = rot_loadings(1),
                         randomStarts = 20)
  )
  expect_equal(
    oracle$loadings,
    efa_procrustes(unrot_loadings(2), rot_loadings(1), rotation = "oblique")$loadings,
    tolerance = 1e-3,
    ignore_attr = TRUE
  )
})


test_that("oblique efa_procrustes recovers a known oblique transform from its default start", {
  skip_on_cran()

  # Self-recovery: with A = L0 %*% t(T0) and target L0, the oblique optimum is T0
  # with objective 0. The default warm orthogonal-Procrustes start lands in the
  # global basin and recovers the generating loadings and factor correlations.
  L0 <- matrix(0, 9, 3)
  L0[1:3, 1] <- c(.75, .70, .65)
  L0[4:6, 2] <- c(.72, .68, .62)
  L0[7:9, 3] <- c(.80, .74, .66)
  L0[1, 2] <- .15; L0[5, 3] <- .12; L0[8, 1] <- .10
  Phi0 <- matrix(c(1, .3, .2, .3, 1, .25, .2, .25, 1), 3)
  T0 <- chol(Phi0)                       # unit-norm columns; t(T0) %*% T0 == Phi0
  A <- L0 %*% t(T0)

  fit <- efa_procrustes(A, Target = L0, rotation = "oblique")

  expect_lt(fit$value, 1e-6)
  expect_true(fit$convergence)
  expect_equal(unclass(fit$loadings), L0, tolerance = 1e-4, ignore_attr = TRUE)
  expect_equal(fit$Phi, Phi0, tolerance = 1e-4, ignore_attr = TRUE)

  # matches the GPArotation target-rotation oracle on the same clean problem
  skip_if_not_installed("GPArotation")
  oracle <- GPArotation::targetQ(A, Target = L0)
  expect_equal(unclass(fit$loadings), unclass(oracle$loadings),
               tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("oblique efa_procrustes warm start escapes identity-start local minima", {
  skip_on_cran()

  # On the 6-factor bootstrap loadings the identity start is trapped in poor local
  # minima (still reporting convergence) for some matrices; the default warm
  # orthogonal-Procrustes start reaches a far lower target-criterion value.
  target <- rot_loadings(1)
  warm <- vapply(2:10, function(d)
    efa_procrustes(unrot_loadings(d), Target = target, rotation = "oblique")$value,
    numeric(1))
  identity_start <- vapply(2:10, function(d)
    efa_procrustes(unrot_loadings(d), Target = target, rotation = "oblique",
               T_init = diag(6))$value,
    numeric(1))

  # The default start is never worse than the identity start. Where the two starts land in
  # the same basin the objectives differ only by the solver's own stopping slack: the
  # projected-gradient norm is driven below oblique_eps = 1e-5, and at a flat oblique optimum
  # the corresponding objective gap is ||g||^2 / (2 * lambda), which for a small curvature
  # lambda can exceed 1e-16-scale rounding by several orders. The slack is therefore stated
  # at the solver's own resolution rather than at machine precision -- this is an all() over
  # nine bootstrap matrices, so a single such tie would otherwise fail it. The contract is
  # "never worse", which the >1-unit rescue asserted next shows is not vacuous.
  expect_true(all(warm <= identity_start + 1e-4))
  # ... and rescues at least one matrix the identity start traps badly
  expect_gt(max(identity_start - warm), 1)
})

test_that("oblique multi-start selection keeps the lowest objective, not the converged one", {
  # A converged suboptimal local minimum must not be preferred over a strictly
  # better but not-yet-converged fit. Use a reflected-factor stationary point
  # (a genuine converged local minimum) as the primary start, then let the random
  # starts reach the global optimum under a budget too small to converge there.
  set.seed(6)
  k <- 2
  p <- 8
  L0 <- matrix(runif(p * k, -.1, .1), p, k)
  for (j in seq_len(k)) L0[((j - 1) * 4 + 1):(j * 4), j] <- runif(4, .5, .85)
  T0 <- chol(matrix(c(1, .7, .7, 1), 2))
  A <- L0 %*% t(T0)

  # a converged suboptimal stationary point (objective ~2.13 vs 0 at the global)
  trapped <- efa_procrustes(A, Target = L0, rotation = "oblique",
                        T_init = T0 %*% diag(c(-1, 1)), oblique_maxit = 1000)$T

  set.seed(2)
  fit <- efa_procrustes(A, Target = L0, rotation = "oblique", T_init = trapped,
                    oblique_random_starts = 12, oblique_screen_keep = 12,
                    oblique_maxit = 3, oblique_triage_maxit = 2)

  # the returned fit is the minimum-objective start ...
  expect_equal(fit$value, min(fit$all_values))
  expect_lt(fit$value, 0.5)
  # ... even though the primary start (index 1) converged at a much higher
  # objective: a convergence-first rule would have returned that worse solution.
  expect_true(as.logical(fit$all_converged[1]))
  expect_gt(fit$all_values[1], fit$value + 0.5)
})


test_that("efa_procrustes rejects invalid matrices and oblique controls", {
  A <- rot_loadings(1)

  expect_error(efa_procrustes("a", A), class = "efa_not_matrix")
  expect_error(efa_procrustes(A, replace(A, 1, NA)), class = "efa_nonfinite")
  expect_error(efa_procrustes(matrix(numeric(0), 0, 0), A), class = "efa_empty_matrix")

  # the oblique control checks run only on the oblique branch
  ob <- function(...) efa_procrustes(A, A, rotation = "oblique", ...)
  expect_error(ob(oblique_maxit = 2.5), class = "efa_not_integer")
  expect_error(ob(oblique_eps = -1), class = "efa_out_of_range")
  expect_error(ob(oblique_eps = c(1, 2)), class = "efa_not_scalar")
  expect_error(ob(oblique_normalize = "y"), class = "efa_not_flag")

  k <- ncol(A)
  T_zero <- diag(k)
  T_zero[, 1] <- 0
  expect_error(ob(T_init = T_zero), class = "efa_zero_column")
  # nonzero columns throughout, but rank k - 1
  T_sing <- diag(k)
  T_sing[, k] <- T_sing[, 1]
  expect_error(ob(T_init = T_sing), class = "efa_singular")
})


test_that("an asymmetric S warns but still rotates", {
  A <- unrot_loadings(2)
  S <- crossprod(A)
  S[1, 2] <- S[1, 2] + 1

  # S should normally be crossprod(A); an asymmetric one is almost always a mistake,
  # but it is a legal weight matrix, so the rotation proceeds after the warning.
  expect_warning(fit <- efa_procrustes(A, rot_loadings(1), rotation = "oblique", S = S),
                 class = "efa_not_symmetric")
  expect_identical(dim(unclass(fit$loadings)), dim(unclass(A)))
})


rm(boot_rows, efa_at, unrot_loadings, rot_loadings)
