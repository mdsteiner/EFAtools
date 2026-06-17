
# list of EFAs from similar datasets (10 bootstrap samples)
set.seed(42)
efa_list <- lapply(1:10,
                   function (x) {
                     EFA(DOSPERT_raw[sample(1:nrow(DOSPERT_raw), replace = TRUE),],
                         n_factors = 6, method = "ML", rotation = "promax")
                   })

# loadings lists
unrot_loadings <- lapply(efa_list, `[[`, "unrot_loadings")
rot_loadings <- lapply(efa_list, `[[`, "rot_loadings")


test_that("PROCRUSTES matches psych::Procrustes and GPArotation::targetQ outputs", {
  expect_equal(
    psych::Procrustes(unrot_loadings[[2]], rot_loadings[[1]])$loadings,
    PROCRUSTES(unrot_loadings[[2]], rot_loadings[[1]])$loadings,
    ignore_attr = TRUE
    )
  # The default single start of GPArotation::targetQ() is itself trapped in a
  # local minimum on this matrix; random starts are needed for it to reach the
  # global oblique target-rotation optimum that PROCRUSTES() now finds.
  skip_if_not_installed("GPArotation")
  set.seed(42)
  oracle <- suppressWarnings(
    GPArotation::targetQ(unrot_loadings[[2]], Target = rot_loadings[[1]],
                         randomStarts = 20)
  )
  expect_equal(
    oracle$loadings,
    PROCRUSTES(unrot_loadings[[2]], rot_loadings[[1]], rotation = "oblique")$loadings,
    tolerance = 1e-3,
    ignore_attr = TRUE
  )
})


test_that("oblique PROCRUSTES recovers a known oblique transform from its default start", {
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

  fit <- PROCRUSTES(A, Target = L0, rotation = "oblique")

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

test_that("oblique PROCRUSTES warm start escapes identity-start local minima", {
  skip_on_cran()

  # On the 6-factor bootstrap loadings the identity start is trapped in poor local
  # minima (still reporting convergence) for some matrices; the default warm
  # orthogonal-Procrustes start reaches a far lower target-criterion value.
  target <- rot_loadings[[1]]
  warm <- vapply(2:10, function(d)
    PROCRUSTES(unrot_loadings[[d]], Target = target, rotation = "oblique")$value,
    numeric(1))
  identity_start <- vapply(2:10, function(d)
    PROCRUSTES(unrot_loadings[[d]], Target = target, rotation = "oblique",
               T_init = diag(6))$value,
    numeric(1))

  # the default start is never worse than the identity start ...
  expect_true(all(warm <= identity_start + 1e-6))
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
  trapped <- PROCRUSTES(A, Target = L0, rotation = "oblique",
                        T_init = T0 %*% diag(c(-1, 1)), oblique_maxit = 1000)$T

  set.seed(2)
  fit <- PROCRUSTES(A, Target = L0, rotation = "oblique", T_init = trapped,
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


rm(efa_list, unrot_loadings, rot_loadings)
