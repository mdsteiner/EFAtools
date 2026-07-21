# A three-factor population built from a shipped loading pattern and moderate
# factor intercorrelations; reused across the recovery and reproducibility tests.
Lambda_pop <- population_models$loadings$baseline
Phi_pop <- population_models$phis_3$moderate
R_pop <- efa_simulate(Lambda = Lambda_pop, Phi = Phi_pop, return_pop = TRUE)$population

test_that("return_pop yields the population correlation matrix", {
  expect_true(is.matrix(R_pop))
  expect_equal(dim(R_pop), c(nrow(Lambda_pop), nrow(Lambda_pop)))
  expect_equal(diag(R_pop), rep(1, nrow(Lambda_pop)), ignore_attr = TRUE)
  expect_true(isSymmetric(R_pop))
  expect_equal(rownames(R_pop), rownames(Lambda_pop))
  # Assembled from Lambda * Phi * Lambda' + Psi and standardized: the within-factor
  # entries are .6 * 1 * .6 = .36, the between-factor entries .6 * .3 * .6 = .108.
  expect_equal(R_pop[1, 2], 0.36)
  expect_equal(R_pop[1, 7], 0.108)
  # Building from Lambda/Phi and passing the resulting matrix as R agree.
  expect_equal(efa_simulate(R = R_pop, return_pop = TRUE)$population, R_pop)
})

test_that("simulated data recover the population correlation matrix", {
  # At a large N the sample correlations of normal data converge to the population
  # correlation; check both input paths (model components and a supplied R).
  dat_model <- efa_simulate(N = 1e5, Lambda = Lambda_pop, Phi = Phi_pop, seed = 42)$data
  dat_R <- efa_simulate(N = 1e5, R = R_pop, seed = 42)$data

  expect_true(is.matrix(dat_model))
  expect_false(is.data.frame(dat_model))
  expect_null(attr(dat_model, "class"))
  expect_equal(dim(dat_model), c(1e5, nrow(Lambda_pop)))
  expect_equal(colnames(dat_model), rownames(Lambda_pop))

  expect_lt(max(abs(stats::cor(dat_model) - R_pop)), 0.03)
  expect_lt(max(abs(stats::cor(dat_R) - R_pop)), 0.03)
})

test_that("n_datasets > 1 returns a list of datasets", {
  sims <- efa_simulate(N = 200, R = R_pop, n_datasets = 3, seed = 1)$data
  expect_type(sims, "list")
  expect_length(sims, 3)
  expect_true(all(vapply(sims, is.matrix, logical(1))))
  expect_equal(dim(sims[[1]]), c(200, nrow(R_pop)))
  # future.seed = TRUE binds each replicate to its own stream, so the datasets differ.
  expect_false(isTRUE(all.equal(sims[[1]], sims[[2]])))
})

test_that("a positive-semidefinite but singular population uses the eigen fallback", {
  # Two perfectly correlated variables: the correlation matrix is singular (a zero
  # eigenvalue), so the Cholesky fails and the symmetric eigen square root is used.
  R_sing <- matrix(1, 2, 2)
  dat <- efa_simulate(N = 500, R = R_sing, seed = 7)$data
  expect_equal(dim(dat), c(500, 2))
  expect_equal(stats::cor(dat)[1, 2], 1, tolerance = 1e-8)
})

test_that("a fixed seed is reproducible and leaves the RNG stream unchanged", {
  a <- efa_simulate(N = 100, R = R_pop, seed = 123)
  b <- efa_simulate(N = 100, R = R_pop, seed = 123)
  expect_identical(a, b)

  # A supplied seed must not have a lasting side effect on the caller's stream.
  set.seed(1)
  state_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  invisible(efa_simulate(N = 100, R = R_pop, seed = 99))
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   state_before)

  # When no stream existed yet, the one set.seed() creates is removed again.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  invisible(efa_simulate(N = 100, R = R_pop, seed = 5))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("a fixed seed is reproducible at 1 vs 2 workers", {
  skip_on_cran()
  # Each dataset is drawn in its own future, and future.seed = TRUE binds every
  # replicate's RNG stream to its index. With a fixed `seed` the draws must
  # therefore be identical regardless of the number of workers. The multisession
  # workers are fresh R processes that load the installed package, so run this
  # under devtools::check() / after devtools::install() for the worker code to
  # match the main process. (multicore is unavailable on Windows.)
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  run <- function() efa_simulate(N = 200, R = R_pop, n_datasets = 4, seed = 2024)

  future::plan(future::sequential)
  one <- run()

  future::plan(future::multisession, workers = 2)
  two <- run()

  expect_equal(one, two, tolerance = 1e-10)
})

test_that("invalid population specifications raise a classed error", {
  # Neither Lambda nor R.
  expect_error(efa_simulate(N = 10), class = "efa_simulate_input")
  # Both Lambda and R.
  expect_error(efa_simulate(N = 10, Lambda = Lambda_pop, R = R_pop),
               class = "efa_simulate_input")
  # Phi/Psi supplied alongside a ready R.
  expect_error(efa_simulate(N = 10, R = R_pop, Phi = Phi_pop),
               class = "efa_simulate_input")
  # A non-symmetric R.
  R_asym <- R_pop
  R_asym[1, 2] <- R_asym[1, 2] + 0.1
  expect_error(efa_simulate(N = 10, R = R_asym), class = "efa_simulate_input")
})

test_that("an indefinite population is rejected", {
  # A symmetric matrix with a negative eigenvalue is not a valid covariance.
  R_indef <- matrix(c(1, 1, 0,
                      1, 1, 1,
                      0, 1, 1), 3, 3)
  expect_lt(min(eigen(R_indef, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_error(efa_simulate(N = 10, R = R_indef), class = "efa_simulate_not_pd")
})

test_that("a Heywood case in the factor model is rejected", {
  # A single-factor loading above 1 implies a communality above 1, leaving no
  # unique variance to simulate from.
  Lambda_hey <- matrix(1.2, 3, 1)
  expect_error(efa_simulate(N = 10, Lambda = Lambda_hey),
               class = "efa_simulate_heywood")
})

test_that("non-finite or non-symmetric inputs raise a classed error", {
  # any.missing = FALSE rejects NA/NaN but not Inf, so finiteness is checked
  # explicitly on every population input.
  L_inf <- Lambda_pop
  L_inf[1, 1] <- Inf
  expect_error(efa_simulate(N = 10, Lambda = L_inf), class = "efa_simulate_input")

  R_inf <- R_pop
  R_inf[1, 1] <- Inf
  expect_error(efa_simulate(N = 10, R = R_inf), class = "efa_simulate_input")

  # A matrix-valued Psi must be finite and symmetric, like Phi and R.
  p <- nrow(Lambda_pop)
  Psi_asym <- diag(0.5, p)
  Psi_asym[1, 2] <- 0.3
  expect_error(efa_simulate(N = 10, Lambda = Lambda_pop, Psi = Psi_asym),
               class = "efa_simulate_input")

  Psi_inf <- diag(0.5, p)
  Psi_inf[1, 1] <- Inf
  expect_error(efa_simulate(N = 10, Lambda = Lambda_pop, Psi = Psi_inf),
               class = "efa_simulate_input")
})

test_that("a valid matrix-valued Psi matches the default unique variances", {
  # A diagonal Psi equal to the standardizing unique variances reproduces the
  # default population exactly.
  Psi <- diag(1 - diag(Lambda_pop %*% Phi_pop %*% t(Lambda_pop)))
  R_from_mat <- efa_simulate(Lambda = Lambda_pop, Phi = Phi_pop, Psi = Psi,
                             return_pop = TRUE)$population
  expect_equal(R_from_mat, R_pop)
})

test_that("return_pop rejects an indefinite population", {
  # The positive-semidefinite screen runs before the return_pop exit, so an
  # invalid population is never returned as 'the population correlation matrix'.
  R_indef <- matrix(c(1, 1, 0,
                      1, 1, 1,
                      0, 1, 1), 3, 3)
  expect_error(efa_simulate(R = R_indef, return_pop = TRUE),
               class = "efa_simulate_not_pd")
})

test_that("a communality of exactly 1 gives a drawable singular population", {
  # psi = 0 (communality exactly 1) is allowed: the population is positive
  # semi-definite but singular and drawn via the eigen fallback, not rejected.
  L1 <- matrix(1, 3, 1)
  R1 <- efa_simulate(Lambda = L1, return_pop = TRUE)$population
  expect_equal(R1, matrix(1, 3, 3), ignore_attr = TRUE)
  dat <- efa_simulate(N = 200, Lambda = L1, seed = 3)$data
  expect_equal(dim(dat), c(200, 3))
})

test_that("variable names fall back to the columns of R when rows are unnamed", {
  R2 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  colnames(R2) <- c("item1", "item2")
  Rp <- efa_simulate(R = R2, return_pop = TRUE)$population
  expect_equal(colnames(Rp), c("item1", "item2"))
  dat <- efa_simulate(N = 50, R = R2, seed = 1)$data
  expect_equal(colnames(dat), c("item1", "item2"))
})

test_that("empirical marginals reproduce the target R and are preserved", {
  # Ruscio-Kaczetow reproduces the population correlation whatever the supplied
  # marginals are; check recovery at a large N and that each drawn column's
  # values come only from the matching (skewed) empirical source column.
  p <- nrow(Lambda_pop)
  set.seed(2024)
  x_skew <- matrix(stats::rexp(400 * p), nrow = 400, ncol = p)

  dat <- efa_simulate(N = 1e4, Lambda = Lambda_pop, Phi = Phi_pop,
                      marginals = "empirical", marginal_data = x_skew,
                      n_factors = 3, seed = 42)$data

  expect_true(is.matrix(dat))
  expect_equal(dim(dat), c(1e4, p))
  expect_equal(colnames(dat), rownames(Lambda_pop))
  expect_lt(max(abs(stats::cor(dat) - R_pop)), 0.05)
  # Each drawn column carries only values from its empirical marginal.
  expect_true(all(vapply(seq_len(p),
                         function(i) all(dat[, i] %in% x_skew[, i]),
                         logical(1))))
})

test_that("empirical n_factors defaults to the model's factor count", {
  p <- nrow(Lambda_pop)
  set.seed(3)
  x_emp <- matrix(stats::rnorm(200 * p), nrow = 200, ncol = p)

  # Built from Lambda: n_factors defaults to ncol(Lambda) and the call runs.
  dat <- efa_simulate(N = 500, Lambda = Lambda_pop, Phi = Phi_pop,
                      marginals = "empirical", marginal_data = x_emp, seed = 1)$data
  expect_equal(dim(dat), c(500, p))

  # Supplied via R: no factor count is implied, so n_factors is required.
  expect_error(
    efa_simulate(N = 500, R = R_pop, marginals = "empirical",
                 marginal_data = x_emp, seed = 1),
    class = "efa_simulate_input")
})

test_that("empirical-path inputs are validated", {
  p <- nrow(Lambda_pop)
  set.seed(5)
  x_emp <- matrix(stats::rnorm(100 * p), nrow = 100, ncol = p)

  # marginal_data is required for the empirical path.
  expect_error(
    efa_simulate(N = 100, Lambda = Lambda_pop, Phi = Phi_pop, marginals = "empirical"),
    class = "efa_simulate_input")
  # Wrong number of columns.
  expect_error(
    efa_simulate(N = 100, Lambda = Lambda_pop, Phi = Phi_pop, marginals = "empirical",
                 marginal_data = x_emp[, 1:5]),
    class = "efa_simulate_input")
  # n_factors must be smaller than the number of variables.
  expect_error(
    efa_simulate(N = 100, Lambda = Lambda_pop, Phi = Phi_pop, marginals = "empirical",
                 marginal_data = x_emp, n_factors = p),
    class = "efa_simulate_input")
  # Non-finite marginal_data.
  x_bad <- x_emp
  x_bad[1, 1] <- Inf
  expect_error(
    efa_simulate(N = 100, Lambda = Lambda_pop, Phi = Phi_pop, marginals = "empirical",
                 marginal_data = x_bad),
    class = "efa_simulate_input")
  # A single-row marginal source would resample a length-one vector (tripping
  # sample()'s length-1 behaviour), so it is rejected.
  expect_error(
    efa_simulate(N = 100, Lambda = Lambda_pop, Phi = Phi_pop, marginals = "empirical",
                 marginal_data = x_emp[1, , drop = FALSE]),
    class = "efa_simulate_input")
  # A constant (zero-variance) column cannot reproduce a correlation.
  x_const <- x_emp
  x_const[, 2] <- 4
  expect_error(
    efa_simulate(N = 100, Lambda = Lambda_pop, Phi = Phi_pop, marginals = "empirical",
                 marginal_data = x_const),
    class = "efa_simulate_input")
  # The reproduction needs at least two cases to draw.
  expect_error(
    efa_simulate(N = 1, Lambda = Lambda_pop, Phi = Phi_pop, marginals = "empirical",
                 marginal_data = x_emp),
    class = "efa_simulate_input")
  # marginal_data / n_factors only apply to the empirical path.
  expect_error(
    efa_simulate(N = 100, Lambda = Lambda_pop, Phi = Phi_pop, marginal_data = x_emp),
    class = "efa_simulate_input")
  expect_error(
    efa_simulate(N = 100, Lambda = Lambda_pop, Phi = Phi_pop, n_factors = 2),
    class = "efa_simulate_input")
})

test_that("empirical n_datasets > 1 returns a reproducible list", {
  p <- nrow(Lambda_pop)
  set.seed(7)
  x_emp <- matrix(stats::rexp(150 * p), nrow = 150, ncol = p)

  run <- function() {
    efa_simulate(N = 300, Lambda = Lambda_pop, Phi = Phi_pop, marginals = "empirical",
                 marginal_data = x_emp, n_factors = 3, n_datasets = 2, seed = 11)
  }
  sims1 <- run()
  sims2 <- run()

  expect_type(sims1$data, "list")
  expect_length(sims1$data, 2)
  expect_true(all(vapply(sims1$data, is.matrix, logical(1))))
  expect_equal(dim(sims1$data[[1]]), c(300, p))
  # A fixed seed is reproducible; each replicate has its own stream, so the two
  # datasets differ.
  expect_identical(sims1, sims2)
  expect_false(isTRUE(all.equal(sims1$data[[1]], sims1$data[[2]])))
})

# Population and empirical-moment estimators (divisor n, matching the theoretical
# targets) shared by the Vale-Maurelli / independent-generator moment tests.
.pop_skew <- function(x) {
  x <- x - mean(x)
  mean(x^3) / mean(x^2)^1.5
}
.pop_exkurt <- function(x) {
  x <- x - mean(x)
  mean(x^4) / mean(x^2)^2 - 3
}
# A compact two-factor population used for moment recovery (six variables keeps
# the large-N draws affordable while the fourth-moment estimator settles).
Lambda_vm <- matrix(0, 6, 2)
Lambda_vm[1:3, 1] <- 0.7
Lambda_vm[4:6, 2] <- 0.7
Phi_vm <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
R_vm <- efa_simulate(Lambda = Lambda_vm, Phi = Phi_vm, return_pop = TRUE)$population

test_that("Vale-Maurelli data carry the target moments and correlation", {
  skip_on_cran()
  sk <- 0.8
  ku <- 1.5
  dat <- efa_simulate(N = 2e5, Lambda = Lambda_vm, Phi = Phi_vm, marginals = "VM",
                      skewness = sk, kurtosis = ku, seed = 42)$data

  expect_true(is.matrix(dat))
  expect_equal(dim(dat), c(2e5, 6))
  # Fleishman marginals are centred with unit variance by construction.
  expect_lt(max(abs(colMeans(dat))), 0.02)
  expect_lt(max(abs(apply(dat, 2, stats::var) - 1)), 0.02)
  # Achieved skewness/(excess) kurtosis match the targets; the fourth moment is a
  # high-variance estimator, so it carries a looser tolerance than skewness.
  expect_lt(max(abs(apply(dat, 2, .pop_skew) - sk)), 0.05)
  expect_lt(max(abs(apply(dat, 2, .pop_exkurt) - ku)), 0.3)
  expect_lt(max(abs(stats::cor(dat) - R_vm)), 0.02)
})

test_that("independent-generator data carry the target moments and correlation", {
  skip_on_cran()
  sk <- 0.8
  ku <- 1.5
  dat <- efa_simulate(N = 2e5, R = R_vm, marginals = "IG",
                      skewness = sk, kurtosis = ku, seed = 42)$data

  expect_equal(dim(dat), c(2e5, 6))
  expect_lt(max(abs(colMeans(dat))), 0.02)
  expect_lt(max(abs(apply(dat, 2, stats::var) - 1)), 0.02)
  expect_lt(max(abs(apply(dat, 2, .pop_skew) - sk)), 0.05)
  expect_lt(max(abs(apply(dat, 2, .pop_exkurt) - ku)), 0.3)
  expect_lt(max(abs(stats::cor(dat) - R_vm)), 0.02)
})

test_that("per-variable skewness and kurtosis are honoured", {
  skip_on_cran()
  # The first three variables are near-normal, the last three markedly skewed.
  sk <- c(rep(0, 3), rep(1, 3))
  ku <- c(rep(0, 3), rep(2, 3))
  dat <- efa_simulate(N = 2e5, R = R_vm, marginals = "VM",
                      skewness = sk, kurtosis = ku, seed = 7)$data
  expect_lt(max(abs(apply(dat, 2, .pop_skew) - sk)), 0.05)
  expect_lt(max(abs(apply(dat, 2, .pop_exkurt) - ku)), 0.3)
})

test_that("unattainable target moments raise a classed error", {
  # Below the universal bound (excess kurtosis >= skewness^2 - 2): skew 3 needs
  # excess kurtosis of at least 7.
  expect_error(
    efa_simulate(N = 10, R = R_vm, marginals = "VM", skewness = 3, kurtosis = 0),
    class = "efa_simulate_infeasible_moments")
  # Inside the universal bound but outside Fleishman's tighter region (its lower
  # kurtosis boundary at skewness 0 is about -1.15): the coefficient solve fails.
  expect_error(
    efa_simulate(N = 10, R = R_vm, marginals = "VM", skewness = 0, kurtosis = -1.4),
    class = "efa_simulate_infeasible_moments")
  # The independent-generator moments can be infeasible even when the requested
  # marginal moments are attainable, because they depend on the population.
  expect_error(
    efa_simulate(N = 10, R = R_vm, marginals = "IG", skewness = 3, kurtosis = 0),
    class = "efa_simulate_infeasible_moments")
  # A modest, attainable marginal request can still be infeasible via the
  # population: strong correlations blow up the implied generator moments past
  # Fleishman's region. The abort is reported for the request, not the internals.
  R_hi <- matrix(0.98, 4, 4)
  diag(R_hi) <- 1
  expect_error(
    efa_simulate(N = 10, R = R_hi, marginals = "IG", skewness = 1, kurtosis = 2),
    class = "efa_simulate_infeasible_moments")
})

test_that("a non-positive-definite Vale-Maurelli intermediate matrix is guarded", {
  # This population is (barely) positive definite, but the strong non-normality
  # inflates the required intermediate correlations enough to make the pairwise
  # intermediate matrix indefinite.
  R_bad <- matrix(c(1, 0.857, 0.918,
                    0.857, 1, 0.585,
                    0.918, 0.585, 1), 3, 3)
  expect_gt(min(eigen(R_bad, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_error(
    efa_simulate(N = 10, R = R_bad, marginals = "VM", skewness = 2.5, kurtosis = 10),
    class = "efa_simulate_intermediate_not_pd")
  # force_pd projects the intermediate matrix to the nearest correlation matrix
  # and warns, rather than aborting.
  expect_warning(
    dat <- efa_simulate(N = 500, R = R_bad, marginals = "VM", skewness = 2.5,
                        kurtosis = 10, force_pd = TRUE, seed = 1),
    class = "efa_simulate_pd_forced")
  expect_equal(dim(dat$data), c(500, 3))
})

test_that("the independent generator requires a positive-definite population", {
  # Two perfectly correlated variables give a singular population, which has no
  # Cholesky factor for the independent-generator mixture.
  R_sing <- matrix(1, 3, 3)
  expect_error(
    efa_simulate(N = 10, R = R_sing, marginals = "IG", skewness = 1, kurtosis = 2),
    class = "efa_simulate_input")
})

test_that("VM/IG moment arguments are validated and gated", {
  # skewness/kurtosis only apply to the VM/IG paths.
  expect_error(efa_simulate(N = 10, R = R_vm, skewness = 1), class = "efa_simulate_input")
  expect_error(efa_simulate(N = 10, R = R_vm, kurtosis = 1), class = "efa_simulate_input")
  # force_pd only applies to the Vale-Maurelli path, not normal or IG.
  expect_error(efa_simulate(N = 10, R = R_vm, force_pd = TRUE), class = "efa_simulate_input")
  expect_error(
    efa_simulate(N = 10, R = R_vm, marginals = "IG", skewness = 1, kurtosis = 2,
                 force_pd = TRUE),
    class = "efa_simulate_input")
  expect_error(
    efa_simulate(N = 10, R = R_vm, marginals = "empirical", marginal_data = R_vm,
                 skewness = 1),
    class = "efa_simulate_input")
  # VM/IG need at least one target moment.
  expect_error(efa_simulate(N = 10, R = R_vm, marginals = "VM"),
               class = "efa_simulate_input")
  # A moment vector must have length 1 or p, and must be a vector (not a matrix).
  expect_error(
    efa_simulate(N = 10, R = R_vm, marginals = "IG", skewness = c(1, 2)),
    class = "efa_simulate_input")
  expect_error(
    efa_simulate(N = 10, R = R_vm, marginals = "VM", skewness = matrix(0.5, 6, 1),
                 kurtosis = 1),
    class = "efa_simulate_input")
})

test_that("VM/IG draws are reproducible and support multiple datasets", {
  a <- efa_simulate(N = 200, R = R_vm, marginals = "VM", skewness = 1, kurtosis = 2,
                    seed = 11)
  b <- efa_simulate(N = 200, R = R_vm, marginals = "VM", skewness = 1, kurtosis = 2,
                    seed = 11)
  expect_identical(a, b)
  expect_equal(colnames(a$data), colnames(R_vm))

  ig_a <- efa_simulate(N = 200, R = R_vm, marginals = "IG", skewness = 1, kurtosis = 2,
                       seed = 11)
  ig_b <- efa_simulate(N = 200, R = R_vm, marginals = "IG", skewness = 1, kurtosis = 2,
                       seed = 11)
  expect_identical(ig_a, ig_b)

  sims <- efa_simulate(N = 200, R = R_vm, marginals = "VM", skewness = 1,
                       kurtosis = 2, n_datasets = 3, seed = 11)$data
  expect_type(sims, "list")
  expect_length(sims, 3)
  # Each replicate has its own stream, so the datasets differ.
  expect_false(isTRUE(all.equal(sims[[1]], sims[[2]])))
})

test_that("categories discretize the draw into ordered category codes", {
  # A five-category cut of normal marginals: the output is an integer matrix of
  # codes 1..5, every category is filled at a moderate N, and categorization
  # attenuates the product-moment correlation below the latent value.
  dat <- efa_simulate(N = 2000, R = R_vm, categories = 5, seed = 1)$data
  expect_true(is.matrix(dat))
  expect_type(dat, "integer")
  expect_equal(dim(dat), c(2000, ncol(R_vm)))
  expect_equal(colnames(dat), colnames(R_vm))
  expect_true(all(dat %in% 1:5))
  expect_true(all(vapply(seq_len(ncol(dat)),
                         function(j) length(unique(dat[, j])) == 5L, logical(1))))
  # Within-factor pair (both load on the first factor): the ordinal Pearson
  # correlation is smaller in magnitude than the latent .49.
  expect_lt(stats::cor(dat)[1, 2], R_vm[1, 2])
})

test_that("match = 'polychoric' round-trips through .polychoric() to the target R", {
  skip_on_cran()
  # Normal latents thresholded into ordered categories: the two-step polychoric
  # estimator is consistent for the latent correlation (Olsson, 1979), so at a
  # large N its estimate returns the target correlation. Both equally probable
  # categories and asymmetric category proportions must round-trip.
  dat_eq <- efa_simulate(N = 1e5, R = R_vm, categories = 5, match = "polychoric",
                         seed = 42)$data
  R_eq <- .polychoric(dat_eq)$R
  expect_lt(max(abs(R_eq - R_vm)), 0.02)

  props <- rep(list(c(0.1, 0.2, 0.4, 0.3)), ncol(R_vm))
  dat_sk <- efa_simulate(N = 1e5, R = R_vm, categories = props, match = "polychoric",
                         seed = 42)$data
  R_sk <- .polychoric(dat_sk)$R
  expect_lt(max(abs(R_sk - R_vm)), 0.02)
})

test_that("polychoric matching recovers the requested category proportions", {
  skip_on_cran()
  # With normal latents, cutting at qnorm(cumsum(p)) gives population category
  # proportions equal to `p`; at a large N the realized proportions match.
  props <- c(0.1, 0.2, 0.4, 0.3)
  dat <- efa_simulate(N = 1e5, R = R_vm, categories = rep(list(props), ncol(R_vm)),
                      match = "polychoric", seed = 7)$data
  realized <- vapply(seq_len(ncol(dat)),
                     function(j) tabulate(dat[, j], nbins = 4L) / nrow(dat),
                     numeric(4L))
  expect_lt(max(abs(realized - props)), 0.01)
})

test_that("Vale-Maurelli thresholds recover the requested category proportions", {
  skip_on_cran()
  # A VM draw is X = f(Z) with f the monotone Fleishman cubic, so P(X <= f(tau)) =
  # P(Z <= tau): cutting at the mapped thresholds f(tau) reproduces the proportions
  # tau was built for, while the raw normal-scale tau would cut the wrong quantiles.
  sk <- 1.5
  ku <- 4
  dat <- efa_simulate(N = 1e5, R = R_vm, marginals = "VM", skewness = sk, kurtosis = ku,
                      categories = 5, seed = 7)$data
  realized <- vapply(seq_len(ncol(dat)),
                     function(j) tabulate(dat[, j], nbins = 5L) / nrow(dat),
                     numeric(5L))
  expect_lt(max(abs(realized - 0.2)), 0.01)

  # Asymmetric requested proportions likewise.
  props <- c(0.1, 0.2, 0.4, 0.3)
  dat_p <- efa_simulate(N = 1e5, R = R_vm, marginals = "VM", skewness = sk, kurtosis = ku,
                        categories = rep(list(props), ncol(R_vm)), seed = 7)$data
  realized_p <- vapply(seq_len(ncol(dat_p)),
                       function(j) tabulate(dat_p[, j], nbins = 4L) / nrow(dat_p),
                       numeric(4L))
  expect_lt(max(abs(realized_p - props)), 0.01)

  # Platykurtic marginals (excess kurtosis below 0) are the case a global-monotonicity
  # test would wrongly reject: d < 0 there, yet the cubic increases across the
  # thresholds, so the proportions must come back just as exactly.
  dat_f <- efa_simulate(N = 1e5, R = R_vm, marginals = "VM", skewness = 0, kurtosis = -1,
                        categories = 5, seed = 7)$data
  realized_f <- vapply(seq_len(ncol(dat_f)),
                       function(j) tabulate(dat_f[, j], nbins = 5L) / nrow(dat_f),
                       numeric(5L))
  expect_lt(max(abs(realized_f - 0.2)), 0.01)

  # The mapping is what buys this: cutting the same draw at the unmapped normal
  # thresholds misses the request by far more than the Monte Carlo error above.
  cont <- efa_simulate(N = 1e5, R = R_vm, marginals = "VM", skewness = sk,
                       kurtosis = ku, seed = 7)$data
  naive <- vapply(seq_len(ncol(cont)),
                  function(j) tabulate(findInterval(cont[, j], stats::qnorm((1:4) / 5)) + 1L,
                                       nbins = 5L) / nrow(cont),
                  numeric(5L))
  expect_gt(max(abs(naive - 0.2)), 0.05)
})

test_that(".vm_thresholds maps a cubic that increases across the thresholds", {
  tau <- stats::qnorm((1:4) / 5)
  # The normal marginal has coefficients (0, 1, 0, 0): an identity map that leaves the
  # thresholds where they are, with no fallback reported.
  norm_map <- .vm_thresholds(list(tau), matrix(c(0, 1, 0, 0), 1, 4))
  expect_equal(norm_map$thresholds[[1]], tau)
  expect_length(norm_map$non_monotone, 0L)

  ft <- .fleishman_table(1.5, 4)
  expect_equal(.vm_thresholds(list(tau), ft)$thresholds[[1]],
               ft[1, 1L] + ft[1, 2L] * tau + ft[1, 3L] * tau^2 + ft[1, 4L] * tau^3)

  # A platykurtic marginal has d < 0, so f' has real roots and the cubic is not monotone
  # on the whole line -- but its turning points are far out in the tails, so it does
  # increase across ordinary thresholds and must still be mapped.
  ft_flat <- .fleishman_table(0, -1)
  expect_lt(ft_flat[1, 4L], 0)
  flat_map <- .vm_thresholds(list(tau), ft_flat)
  expect_length(flat_map$non_monotone, 0L)
  expect_false(is.unsorted(flat_map$thresholds[[1]]))

  # Pushing a threshold past the turning point (|z| = 2.25 here) is what forces the
  # fallback: the mapping would no longer preserve the category order.
  far <- c(-2.6, tau)
  fallback <- .vm_thresholds(list(far), ft_flat)
  expect_equal(fallback$thresholds[[1]], far)
  expect_identical(fallback$non_monotone, 1L)
})

test_that(".vm_thresholds falls back when the tail folds mass back below an outer cut", {
  # A threshold just INSIDE the turning point passes the derivative screen, but the cubic
  # re-crosses the mapped cut a little further out and the tail mass beyond that point
  # lands in the wrong category: at excess kurtosis -1 the turning point is |z| = 2.25 and
  # a top category of .014 (threshold 2.197) would come out near .0035 -- a factor of four
  # off, silently, under the derivative screen alone. The mass screen catches it.
  ft_flat <- .fleishman_table(0, -1)
  tau_tail <- stats::qnorm(0.986)
  expect_gt(ft_flat[1, 2L] + 2 * ft_flat[1, 3L] * tau_tail +
              3 * ft_flat[1, 4L] * tau_tail^2, 0)   # derivative screen alone passes
  fb <- .vm_thresholds(list(tau_tail), ft_flat)
  expect_identical(fb$non_monotone, 1L)
  expect_equal(fb$thresholds[[1]], tau_tail)

  # Not a platykurtosis-only phenomenon: a skewed leptokurtic marginal dips on one side,
  # and a 5% bottom category (threshold -1.645, turning point -1.69) folds back 84% of
  # the requested mass. It must fall back too.
  ft_skew <- .fleishman_table(1.5, 3)
  fb_skew <- .vm_thresholds(list(stats::qnorm(0.05)), ft_skew)
  expect_identical(fb_skew$non_monotone, 1L)

  # The screen is a tolerance, not a monotonicity dogma: ordinary interior thresholds on
  # the same platykurtic marginal still map (their re-crossings carry negligible mass).
  ok <- .vm_thresholds(list(stats::qnorm((1:4) / 5)), ft_flat)
  expect_length(ok$non_monotone, 0L)
})

test_that("thresholds that outrun the Fleishman cubic warn and keep the normal scale", {
  # Same construction end to end: a strongly platykurtic marginal with a category
  # boundary beyond the cubic's turning point falls back with a classed warning.
  props <- rep(list(c(0.005, 0.7, 0.28, 0.015)), ncol(R_vm))
  # The fallback empties that category too, and necessarily so: the untransformed
  # threshold (-2.58) lies outside the range of the platykurtic draw, which is precisely
  # why the normal-scale cut points do not reproduce the requested proportions.
  expect_warning(
    expect_warning(
      efa_simulate(N = 5000, R = R_vm, marginals = "VM", skewness = 0, kurtosis = -1,
                   categories = props, seed = 11),
      class = "efa_simulate_threshold_fallback"),
    class = "efa_simulate_empty_category")
})

test_that("match = 'thresholds' cuts non-normal marginals without error", {
  # The naive mode works with the standardized non-normal marginals: the
  # Vale-Maurelli draw is cut as is, yielding a valid integer ordinal matrix (its
  # polychoric correlation departs from the target, the documented cost of non-normal
  # marginals).
  dat <- efa_simulate(N = 2000, R = R_vm, marginals = "VM", skewness = 1,
                      kurtosis = 2, categories = 4, match = "thresholds", seed = 3)$data
  expect_type(dat, "integer")
  expect_equal(dim(dat), c(2000, ncol(R_vm)))
  expect_true(all(dat %in% 1:4))
  dat2 <- efa_simulate(N = 2000, R = R_vm, marginals = "VM", skewness = 1,
                       kurtosis = 2, categories = 4, match = "thresholds", seed = 3)$data
  expect_identical(dat, dat2)
})

test_that("match = 'polychoric' rejects non-normal marginals", {
  # The polychoric model assumes a normal latent, so pairing it with non-normal
  # marginals is a contradiction rather than an approximation.
  expect_error(
    efa_simulate(N = 10, R = R_vm, marginals = "VM", skewness = 1, kurtosis = 2,
                 categories = 5, match = "polychoric"),
    class = "efa_simulate_match_conflict")
  expect_error(
    efa_simulate(N = 10, R = R_vm, marginals = "IG", skewness = 1, kurtosis = 2,
                 categories = 5, match = "polychoric"),
    class = "efa_simulate_match_conflict")
  expect_error(
    efa_simulate(N = 10, R = R_vm, marginals = "empirical", categories = 5,
                 match = "polychoric"),
    class = "efa_simulate_match_conflict")
})

test_that("categories are not supported with empirical marginals", {
  # Ordinal thresholds are on the standard-normal scale, but the empirical draw
  # carries the arbitrary scale of marginal_data, so the combination is rejected
  # (rather than silently returning mis-proportioned or degenerate categories).
  set.seed(1)
  md <- matrix(stats::rexp(200 * ncol(R_vm)), ncol = ncol(R_vm))
  expect_error(
    efa_simulate(N = 100, R = R_vm, marginals = "empirical", marginal_data = md,
                 n_factors = 2, categories = 5),
    class = "efa_simulate_input")
})

test_that("match only applies when categories requests ordinal output", {
  expect_error(efa_simulate(N = 10, R = R_vm, match = "polychoric"),
               class = "efa_simulate_input")
  expect_error(efa_simulate(N = 10, R = R_vm, match = "thresholds"),
               class = "efa_simulate_input")
})

test_that("invalid category specifications are rejected", {
  p <- ncol(R_vm)
  # Counts must be whole numbers of at least two categories, of length 1 or p.
  expect_error(efa_simulate(N = 10, R = R_vm, categories = 1),
               class = "efa_simulate_input")
  expect_error(efa_simulate(N = 10, R = R_vm, categories = 2.5),
               class = "efa_simulate_input")
  expect_error(efa_simulate(N = 10, R = R_vm, categories = rep(3L, p + 1L)),
               class = "efa_simulate_input")
  # A list of proportions must have one element per variable, each a valid
  # probability vector of at least two strictly-positive proportions summing to 1.
  expect_error(efa_simulate(N = 10, R = R_vm, categories = list(c(0.5, 0.5))),
               class = "efa_simulate_input")
  expect_error(
    efa_simulate(N = 10, R = R_vm, categories = rep(list(c(0.2, 0.2, 0.2)), p)),
    class = "efa_simulate_input")
  expect_error(
    efa_simulate(N = 10, R = R_vm, categories = rep(list(c(0, 0.5, 0.5)), p)),
    class = "efa_simulate_input")
  expect_error(efa_simulate(N = 10, R = R_vm, categories = rep(list(1), p)),
               class = "efa_simulate_input")
  # A dimensioned object (matrix/array) is rejected rather than flattened
  # column-major into per-variable counts (a p-entry matrix passes the length check).
  expect_error(efa_simulate(N = 10, R = R_vm, categories = matrix(2:7, 2, 3)),
               class = "efa_simulate_input")
  # ... and likewise for a matrix supplied as a proportion list element.
  expect_error(
    efa_simulate(N = 10, R = R_vm,
                 categories = rep(list(matrix(c(0.25, 0.25, 0.25, 0.25), 2, 2)), p)),
    class = "efa_simulate_input")
  # Proportions that pass the sum tolerance but whose interior cumulative sum would
  # cross 1 are renormalized rather than driving qnorm to NaN (which would crash
  # findInterval with an unclassed error); the tiny last category is left empty. The
  # total sits just inside the sum tolerance (1 + 9e-7 < 1 + 1e-6) rather than exactly
  # on it, so the check is not decided by the last-bit rounding of the sum (which
  # differs across BLAS/platforms).
  props_edge <- rep(list(c(0.5, 0.5000005, 4e-7)), p)
  expect_warning(
    dat_edge <- efa_simulate(N = 500, R = R_vm, categories = props_edge, seed = 1),
    class = "efa_simulate_empty_category")
  expect_type(dat_edge$data, "integer")
})

test_that("empty categories are reported with a warning", {
  # Three cases cannot fill five categories, so some are left empty.
  expect_warning(
    efa_simulate(N = 3, R = R_vm, categories = 5, seed = 1),
    class = "efa_simulate_empty_category")
})

test_that("ordinal draws are reproducible and support multiple datasets", {
  a <- efa_simulate(N = 500, R = R_vm, categories = 4, seed = 11)
  b <- efa_simulate(N = 500, R = R_vm, categories = 4, seed = 11)
  expect_identical(a, b)

  sims <- efa_simulate(N = 500, R = R_vm, categories = 4, n_datasets = 3, seed = 11)$data
  expect_type(sims, "list")
  expect_length(sims, 3)
  expect_true(all(vapply(sims, function(d) is.matrix(d) && is.integer(d), logical(1))))
  expect_equal(dim(sims[[1]]), c(500, ncol(R_vm)))
  expect_false(isTRUE(all.equal(sims[[1]], sims[[2]])))
})

test_that("MCAR holes each variable at the target expected rate", {
  # MCAR removes values independently of the data, so the achieved marginal rate per
  # variable is the target up to sampling error (expected rate, not an exact count).
  dat <- efa_simulate(N = 1e4, R = R_vm, missing = "MCAR", missing_prop = 0.2, seed = 1)$data
  expect_true(anyNA(dat))
  expect_lt(max(abs(colMeans(is.na(dat)) - 0.2)), 0.02)
})

test_that("MAR and MNAR achieve the target expected missing rate", {
  d_mar <- efa_simulate(N = 1e4, R = R_vm, missing = "MAR", missing_prop = 0.25,
                        seed = 2)$data
  d_mnar <- efa_simulate(N = 1e4, R = R_vm, missing = "MNAR", missing_prop = 0.25,
                         seed = 2)$data
  # Calibrating the logistic intercept per variable targets the mean missing
  # probability, so the achieved marginal rate matches the target.
  expect_lt(max(abs(colMeans(is.na(d_mar)) - 0.25)), 0.02)
  expect_lt(max(abs(colMeans(is.na(d_mnar)) - 0.25)), 0.02)
})

test_that("missingness is orthogonal to the marginal distribution", {
  # Missingness applies to any marginal type (only ordinal cutting is gated to the
  # standard-normal-scale marginals); empirical and VM draws accept it unchanged.
  p <- nrow(Lambda_pop)
  set.seed(30)
  x_emp <- matrix(stats::rexp(300 * p), nrow = 300, ncol = p)
  d_emp <- efa_simulate(N = 5000, Lambda = Lambda_pop, Phi = Phi_pop,
                        marginals = "empirical", marginal_data = x_emp, n_factors = 3,
                        missing = "MAR", missing_prop = 0.2, seed = 31)$data
  expect_lt(max(abs(colMeans(is.na(d_emp)) - 0.2)), 0.03)

  d_vm <- efa_simulate(N = 5000, R = R_vm, marginals = "VM", skewness = 1, kurtosis = 2,
                       missing = "MNAR", missing_prop = 0.2, seed = 32)$data
  expect_lt(max(abs(colMeans(is.na(d_vm)) - 0.2)), 0.03)
})

test_that(".calibrate_missing_intercept hits the target mean probability", {
  set.seed(14)
  z <- as.vector(scale(stats::rnorm(1000)))
  for (b in c(0, 0.5, 1.5, 3)) {
    for (target in c(0.05, 0.2, 0.5, 0.8)) {
      a <- .calibrate_missing_intercept(z, b, target)
      expect_equal(mean(stats::plogis(a + b * z)), target, tolerance = 1e-6)
    }
  }
  # A constant predictor (z = 0) reduces to plogis(a) = target.
  a0 <- .calibrate_missing_intercept(rep(0, 100), 2, 0.3)
  expect_equal(stats::plogis(a0), 0.3, tolerance = 1e-6)
})

test_that("the MAR mask depends on the predictor, not the variable itself", {
  # Independent columns, so a dependence on the predictor cannot leak through
  # correlation with the variable being holed.
  set.seed(11)
  n <- 2e4
  x <- matrix(stats::rnorm(n * 3), n, 3)
  # pred_idx = c(2, 3, 1): variable 1's missingness is driven by variable 2.
  out <- .apply_missing(x, "MAR", prop = 0.3, strength = 3, pred_idx = c(2L, 3L, 1L))
  miss1 <- as.numeric(is.na(out[, 1]))
  expect_gt(stats::cor(miss1, x[, 2]), 0.15)      # depends on the predictor (var 2)
  expect_lt(abs(stats::cor(miss1, x[, 1])), 0.05) # not on variable 1 itself
  expect_equal(mean(miss1), 0.3, tolerance = 0.02)
  # The non-holed values are left byte-identical to the input.
  keep <- !is.na(out[, 1])
  expect_identical(out[keep, 1], x[keep, 1])
})

test_that("the MNAR mask depends on the variable's own value", {
  set.seed(12)
  n <- 2e4
  x <- matrix(stats::rnorm(n * 3), n, 3)
  out <- .apply_missing(x, "MNAR", prop = 0.3, strength = 3, pred_idx = NULL)
  miss1 <- as.numeric(is.na(out[, 1]))
  expect_gt(stats::cor(miss1, x[, 1]), 0.15)
  # A positive slope preferentially removes high values, so the observed mean shifts low.
  expect_lt(mean(out[!is.na(out[, 1]), 1]), -0.1)
})

test_that(".apply_missing reads predictors from the complete draw", {
  # Column 2's predictor is column 1, which is holed first. The mask must use the
  # complete values of column 1, not the NA-holed ones (else the probabilities would be
  # NA and the intended dependence would break).
  set.seed(13)
  n <- 8000
  x <- matrix(stats::rnorm(n * 2), n, 2)
  out <- .apply_missing(x, "MAR", prop = 0.3, strength = 3, pred_idx = c(2L, 1L))
  miss2 <- as.numeric(is.na(out[, 2]))
  expect_false(anyNA(miss2))
  expect_gt(stats::cor(miss2, x[, 1]), 0.15)
  expect_equal(mean(miss2), 0.3, tolerance = 0.03)
})

test_that(".resolve_missing_predictor resolves and validates specs", {
  p <- 5L
  vn <- paste0("v", seq_len(p))
  # NULL gives the cyclic next neighbour.
  expect_identical(.resolve_missing_predictor(NULL, p, vn), c(2L, 3L, 4L, 5L, 1L))
  # Names resolve to indices.
  expect_identical(.resolve_missing_predictor(vn[c(2:5, 1)], p, vn),
                   c(2L, 3L, 4L, 5L, 1L))
  # A self-predicting index errors (that would be MNAR).
  expect_error(.resolve_missing_predictor(seq_len(p), p, vn),
               class = "efa_simulate_input")
  # A single shared predictor is not allowed (wrong length).
  expect_error(.resolve_missing_predictor(2L, p, vn), class = "efa_simulate_input")
  # Out-of-range index, unknown name, and p < 2 error.
  expect_error(.resolve_missing_predictor(c(2L, 3L, 4L, 5L, 99L), p, vn),
               class = "efa_simulate_input")
  expect_error(.resolve_missing_predictor("nope", p, vn), class = "efa_simulate_input")
  expect_error(.resolve_missing_predictor(NULL, 1L, "v1"), class = "efa_simulate_input")
})

test_that("missingness is reproducible and supports multiple datasets", {
  a <- efa_simulate(N = 300, R = R_vm, missing = "MAR", missing_prop = 0.2, seed = 5)
  b <- efa_simulate(N = 300, R = R_vm, missing = "MAR", missing_prop = 0.2, seed = 5)
  expect_identical(a, b)

  sims <- efa_simulate(N = 300, R = R_vm, missing = "MCAR", missing_prop = 0.2,
                       n_datasets = 3, seed = 5)$data
  expect_type(sims, "list")
  expect_length(sims, 3)
  expect_true(all(vapply(sims, anyNA, logical(1))))
  # Each replicate has its own stream, so the NA patterns differ.
  expect_false(identical(is.na(sims[[1]]), is.na(sims[[2]])))
})

test_that("missingness combines with ordinal categories", {
  # The mask acts on the latent draw and the NAs propagate into the codes, so the NA
  # positions match a continuous draw at the same seed and the codes stay valid.
  dat <- efa_simulate(N = 3000, R = R_vm, categories = 4, missing = "MCAR",
                      missing_prop = 0.15, seed = 6)$data
  cont <- efa_simulate(N = 3000, R = R_vm, missing = "MCAR", missing_prop = 0.15,
                       seed = 6)$data
  expect_type(dat, "integer")
  expect_true(anyNA(dat))
  expect_true(all(dat[!is.na(dat)] %in% 1:4))
  expect_lt(max(abs(colMeans(is.na(dat)) - 0.15)), 0.03)
  # Discretizing does not move the holes.
  expect_identical(is.na(dat), is.na(cont))
})

test_that("MCAR and MNAR work with a single variable, MAR does not", {
  R1 <- matrix(1, 1, 1)
  expect_true(anyNA(efa_simulate(N = 500, R = R1, missing = "MCAR",
                                 missing_prop = 0.2, seed = 1)$data))
  expect_true(anyNA(efa_simulate(N = 500, R = R1, missing = "MNAR",
                                 missing_prop = 0.2, seed = 1)$data))
  # MAR needs another variable to predict from.
  expect_error(
    efa_simulate(N = 500, R = R1, missing = "MAR", missing_prop = 0.2),
    class = "efa_simulate_input")
})

test_that("MCAR draws its runif per column, leaving the pre-mask draw intact", {
  # MCAR consumes one runif(n) per column just like MAR/MNAR, so the non-missing values
  # equal the corresponding complete draw (the mask only overwrites, never reorders).
  full <- efa_simulate(N = 1000, R = R_vm, seed = 21)$data
  holed <- efa_simulate(N = 1000, R = R_vm, missing = "MCAR", missing_prop = 0.3,
                        seed = 21)$data
  keep <- !is.na(holed)
  expect_equal(holed[keep], full[keep])
})

test_that("MAR predictor can be given by name or index", {
  R_named <- R_vm
  nm <- paste0("v", seq_len(ncol(R_vm)))
  dimnames(R_named) <- list(nm, nm)
  a <- efa_simulate(N = 400, R = R_named, missing = "MAR", missing_prop = 0.2,
                    missing_predictor = nm[c(2:6, 1)], seed = 8)$data
  b <- efa_simulate(N = 400, R = R_named, missing = "MAR", missing_prop = 0.2,
                    missing_predictor = c(2:6, 1), seed = 8)$data
  expect_identical(is.na(a), is.na(b))
})

test_that("missing-data arguments are validated and gated", {
  p <- ncol(R_vm)
  # Companion settings require a mechanism.
  expect_error(efa_simulate(N = 100, R = R_vm, missing_prop = 0.1),
               class = "efa_simulate_input")
  expect_error(efa_simulate(N = 100, R = R_vm, missing_strength = 1),
               class = "efa_simulate_input")
  expect_error(efa_simulate(N = 100, R = R_vm, missing_predictor = 2),
               class = "efa_simulate_input")
  # missing_prop is required and must be a single proportion in (0, 1).
  expect_error(efa_simulate(N = 100, R = R_vm, missing = "MCAR"),
               class = "efa_simulate_input")
  expect_error(efa_simulate(N = 100, R = R_vm, missing = "MCAR", missing_prop = 0),
               class = "efa_simulate_input")
  expect_error(efa_simulate(N = 100, R = R_vm, missing = "MCAR", missing_prop = 1),
               class = "efa_simulate_input")
  expect_error(efa_simulate(N = 100, R = R_vm, missing = "MCAR", missing_prop = 1.2),
               class = "efa_simulate_input")
  expect_error(efa_simulate(N = 100, R = R_vm, missing = "MCAR",
                            missing_prop = c(0.1, 0.2)),
               class = "efa_simulate_input")
  # strength/predictor do not apply to MCAR.
  expect_error(
    efa_simulate(N = 100, R = R_vm, missing = "MCAR", missing_prop = 0.1,
                 missing_strength = 1),
    class = "efa_simulate_input")
  expect_error(
    efa_simulate(N = 100, R = R_vm, missing = "MCAR", missing_prop = 0.1,
                 missing_predictor = 2),
    class = "efa_simulate_input")
  # A non-finite strength is rejected on the value-dependent paths.
  expect_error(
    efa_simulate(N = 100, R = R_vm, missing = "MNAR", missing_prop = 0.1,
                 missing_strength = Inf),
    class = "efa_simulate_input")
  # MNAR has no separate predictor.
  expect_error(
    efa_simulate(N = 100, R = R_vm, missing = "MNAR", missing_prop = 0.1,
                 missing_predictor = 2),
    class = "efa_simulate_input")
  # A MAR predictor cannot be the variable itself, out of range, or unknown.
  expect_error(
    efa_simulate(N = 100, R = R_vm, missing = "MAR", missing_prop = 0.1,
                 missing_predictor = seq_len(p)),
    class = "efa_simulate_input")
  expect_error(
    efa_simulate(N = 100, R = R_vm, missing = "MAR", missing_prop = 0.1,
                 missing_predictor = c(2:6, 99)),
    class = "efa_simulate_input")
})

# ---- model error (none / CB / TKL / WB) --------------------------------------

# The 18-variable, three-factor baseline population gives ample residual degrees of
# freedom (df = 102) for injecting model error.
Lambda_me <- population_models$loadings$baseline
Phi_me <- population_models$phis_3$moderate
p_me <- nrow(Lambda_me)
q_me <- ncol(Lambda_me)
df_me <- ((p_me - q_me)^2 - (p_me + q_me)) / 2

test_that("model error is off by default and when no target is given", {
  # A bare Lambda call returns the exact population (backward-compatible: the default
  # method is CB, but nothing is perturbed without a target).
  sim <- efa_simulate(N = 200, Lambda = Lambda_me, Phi = Phi_me, seed = 1)
  expect_s3_class(sim, "efa_simulated")
  expect_null(sim$model_error)
  rp <- efa_simulate(Lambda = Lambda_me, Phi = Phi_me, return_pop = TRUE)
  expect_equal(rp$population[1, 2], 0.36)
  expect_null(rp$model_error)
})

test_that("CB model error achieves the target RMSEA to numerical precision", {
  skip_on_cran()
  for (tr in c(0.03, 0.05, 0.08)) {
    sim <- efa_simulate(N = 100, Lambda = Lambda_me, Phi = Phi_me,
                        target_rmsea = tr, seed = 1)
    expect_equal(sim$model_error$method, "CB")
    # Self-consistent: the achieved RMSEA equals the target (the model stays the minimizer).
    expect_equal(sim$model_error$rmsea, tr, tolerance = 1e-6)
    # The perturbed population is a proper, positive-definite correlation matrix.
    expect_equal(diag(sim$population), rep(1, p_me), ignore_attr = TRUE)
    expect_gt(min(eigen(sim$population, symmetric = TRUE, only.values = TRUE)$values), 0)
    expect_true(is.finite(sim$model_error$kappa))
  }
})

test_that("CB achieved RMSEA matches an independent factanal fit", {
  skip_on_cran()
  # factanal fits the best-fitting q-factor model; CB makes the specified model that
  # minimizer, so the minimized discrepancy reproduces the target RMSEA.
  rp <- efa_simulate(Lambda = Lambda_me, Phi = Phi_me, target_rmsea = 0.05,
                     return_pop = TRUE, seed = 2)
  fit <- stats::factanal(covmat = rp$population, factors = q_me, n.obs = 1e6)
  rmsea_min <- sqrt(fit$criteria[["objective"]] / df_me)
  expect_equal(rmsea_min, 0.05, tolerance = 1e-4)
  expect_equal(rp$model_error$rmsea, 0.05, tolerance = 1e-6)
})

test_that("TKL matches a single RMSEA target closely and both targets approximately", {
  skip_on_cran()
  s1 <- efa_simulate(Lambda = Lambda_me, Phi = Phi_me, model_error = "TKL",
                     target_rmsea = 0.05, return_pop = TRUE, seed = 3)
  expect_equal(s1$model_error$method, "TKL")
  expect_equal(s1$model_error$rmsea, 0.05, tolerance = 1e-3)
  # Two targets, two knobs: a close compromise rather than an exact match (absolute band).
  s2 <- efa_simulate(Lambda = Lambda_me, Phi = Phi_me, model_error = "TKL",
                     target_rmsea = 0.05, target_cfi = 0.95, return_pop = TRUE, seed = 4)
  expect_lt(abs(s2$model_error$rmsea - 0.05), 0.02)
  expect_lt(abs(s2$model_error$cfi - 0.95), 0.02)
})

test_that("TKL can target the CFI alone", {
  skip_on_cran()
  s <- efa_simulate(Lambda = Lambda_me, Phi = Phi_me, model_error = "TKL",
                    target_cfi = 0.95, return_pop = TRUE, seed = 5)
  expect_lt(abs(s$model_error$cfi - 0.95), 0.02)
})

test_that("WB produces a realized RMSEA in a relaxed band around the target", {
  skip_on_cran()
  # A single inverse-Wishart draw: the realized RMSEA varies around (and tends to exceed)
  # the target, so check a broad band across a few seeds rather than a tight tolerance.
  ach <- vapply(1:8, function(s) {
    efa_simulate(Lambda = Lambda_me, Phi = Phi_me, model_error = "WB",
                 target_rmsea = 0.05, return_pop = TRUE, seed = s)$model_error$rmsea
  }, numeric(1))
  expect_true(all(ach > 0.02 & ach < 0.12))
  # The precision parameter is m = 1 / target_rmsea^2.
  m <- efa_simulate(Lambda = Lambda_me, Phi = Phi_me, model_error = "WB",
                    target_rmsea = 0.05, return_pop = TRUE, seed = 1)$model_error$m
  expect_equal(m, 1 / 0.05^2)
})

test_that("WB rejects a target RMSEA too large for the number of variables", {
  # m = 1 / rmsea^2 must be at least p; rmsea = 0.3 -> m ~ 11 < 18.
  expect_error(
    efa_simulate(Lambda = Lambda_me, Phi = Phi_me, model_error = "WB",
                 target_rmsea = 0.3, return_pop = TRUE),
    class = "efa_simulate_model_error")
})

test_that("model error is reproducible and shared across replicates", {
  a <- efa_simulate(N = 150, Lambda = Lambda_me, Phi = Phi_me,
                    target_rmsea = 0.05, n_datasets = 2, seed = 7)
  b <- efa_simulate(N = 150, Lambda = Lambda_me, Phi = Phi_me,
                    target_rmsea = 0.05, n_datasets = 2, seed = 7)
  # Compared numerically rather than bit for bit: every step from the population
  # perturbation to the draw runs through the BLAS, and a threaded BLAS need not return
  # identical bits for two identical calls within one session.
  expect_equal(a, b)
  # One population perturbation is shared by all replicates, which still differ.
  expect_equal(a$model_error$rmsea, 0.05, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(a$data[[1]], a$data[[2]])))
})

test_that("model error requires a factor model and residual degrees of freedom", {
  R <- efa_simulate(Lambda = Lambda_me, Phi = Phi_me, return_pop = TRUE)$population
  # A bare R has no specified factor model to misfit against.
  expect_error(efa_simulate(N = 10, R = R, target_rmsea = 0.05),
               class = "efa_simulate_input")
  # A model with no residual df (q too large) cannot carry model error: p = 3, q = 1 -> df = 0.
  L_sat <- matrix(0.5, 3, 1)
  expect_error(efa_simulate(N = 10, Lambda = L_sat, target_rmsea = 0.05),
               class = "efa_simulate_input")
})

test_that("model-error arguments are validated and gated", {
  # A target with model_error = "none".
  expect_error(efa_simulate(N = 10, Lambda = Lambda_me, Phi = Phi_me,
                            model_error = "none", target_rmsea = 0.05),
               class = "efa_simulate_input")
  # CB / WB target the RMSEA only, not the CFI.
  expect_error(efa_simulate(N = 10, Lambda = Lambda_me, Phi = Phi_me, target_cfi = 0.95),
               class = "efa_simulate_input")
  expect_error(efa_simulate(N = 10, Lambda = Lambda_me, Phi = Phi_me,
                            model_error = "WB", target_rmsea = 0.05, target_cfi = 0.95),
               class = "efa_simulate_input")
  # Out-of-range targets (the range is the open interval (0, 1); a target of 0 misfit or a CFI
  # of exactly 1 means an exact population and is rejected -- omit the target instead).
  expect_error(efa_simulate(N = 10, Lambda = Lambda_me, Phi = Phi_me, target_rmsea = 1.5),
               class = "efa_simulate_input")
  expect_error(efa_simulate(N = 10, Lambda = Lambda_me, Phi = Phi_me, target_rmsea = 0),
               class = "efa_simulate_input")
  expect_error(efa_simulate(N = 10, Lambda = Lambda_me, Phi = Phi_me,
                            model_error = "TKL", target_cfi = 1.5),
               class = "efa_simulate_input")
  expect_error(efa_simulate(N = 10, Lambda = Lambda_me, Phi = Phi_me,
                            model_error = "TKL", target_cfi = 1),
               class = "efa_simulate_input")
  # force_pd needs VM or a model-error target.
  expect_error(efa_simulate(N = 10, Lambda = Lambda_me, Phi = Phi_me, force_pd = TRUE),
               class = "efa_simulate_input")
})

test_that("model error rejects singular and non-factor-structure populations", {
  # A positive-semidefinite but SINGULAR population (two items with a communality of 1 are
  # perfectly correlated) is drawable but has no inverse, so the model-error solvers reject it
  # with a classed error rather than an opaque base-R linear-algebra crash.
  L_sing <- matrix(0, 6, 2)
  L_sing[1:3, 1] <- c(1, 1, 0.5)      # items 1 and 2 both fully explained -> R singular
  L_sing[4:6, 2] <- c(0.7, 0.6, 0.5)
  expect_error(efa_simulate(N = 10, Lambda = L_sing, target_rmsea = 0.05),
               class = "efa_simulate_model_error")

  # A correlated-residual (off-diagonal) Psi makes the population not an exact q-factor
  # structure, which the eigen-recovered loadings cannot represent, so it is rejected.
  Psi_corr <- diag(1 - diag(Lambda_me %*% Phi_me %*% t(Lambda_me)))
  Psi_corr[1, 2] <- Psi_corr[2, 1] <- 0.1
  expect_error(efa_simulate(N = 10, Lambda = Lambda_me, Phi = Phi_me, Psi = Psi_corr,
                            target_rmsea = 0.05),
               class = "efa_simulate_model_error")
  # ... and likewise for TKL (which would otherwise take the square root of a negative
  # uniqueness when a negative Psi diagonal pushes a communality above 1).
  Psi_neg <- diag(c(-0.02, rep(0.75, nrow(Lambda_me) - 1L)))
  expect_error(efa_simulate(N = 10, Lambda = matrix(0.5, nrow(Lambda_me), 1),
                            Psi = Psi_neg, model_error = "TKL", target_rmsea = 0.05),
               class = "efa_simulate_model_error")
})

test_that("model error works with an equivalent diagonal-matrix Psi", {
  skip_on_cran()
  # A diagonal matrix Psi is a proper factor structure, so model error still hits the target.
  Psi_diag <- diag(1 - diag(Lambda_me %*% Phi_me %*% t(Lambda_me)))
  sim <- efa_simulate(Lambda = Lambda_me, Phi = Phi_me, Psi = Psi_diag,
                      target_rmsea = 0.05, return_pop = TRUE, seed = 1)
  expect_equal(sim$model_error$rmsea, 0.05, tolerance = 1e-6)
})

test_that("model error combines with non-normal marginals and missingness", {
  skip_on_cran()
  # Model error is orthogonal to the marginal, ordinal, and missing-data options.
  sim <- efa_simulate(N = 2000, Lambda = Lambda_me, Phi = Phi_me, target_rmsea = 0.05,
                      marginals = "VM", skewness = 1, kurtosis = 2,
                      missing = "MCAR", missing_prop = 0.1, seed = 8)
  expect_equal(sim$model_error$rmsea, 0.05, tolerance = 1e-6)
  expect_true(anyNA(sim$data))
  expect_lt(max(abs(colMeans(is.na(sim$data)) - 0.1)), 0.03)
})

test_that(".efa_population_fit reports the population-limit fit reused from EFA()", {
  # An exact population fits its own model perfectly (RMSEA 0, CFI 1).
  R <- efa_simulate(Lambda = Lambda_me, Phi = Phi_me, return_pop = TRUE)$population
  h2 <- diag(Lambda_me %*% Phi_me %*% t(Lambda_me))
  Cstd <- R
  diag(Cstd) <- h2
  ev <- eigen(Cstd, symmetric = TRUE)
  L <- ev$vectors[, seq_len(q_me)] %*% diag(sqrt(ev$values[seq_len(q_me)]))
  pf <- .efa_population_fit(L, R)
  expect_equal(pf$rmsea, 0, tolerance = 1e-8)
  expect_equal(pf$cfi, 1, tolerance = 1e-8)
  expect_equal(pf$df, df_me)
})

test_that("the efa_simulated object has the documented shape", {
  sim <- efa_simulate(N = 50, Lambda = Lambda_me, Phi = Phi_me, target_rmsea = 0.05, seed = 1)
  expect_s3_class(sim, "efa_simulated")
  expect_named(sim, c("data", "population", "model_error", "settings"))
  expect_true(is.matrix(sim$data))
  expect_equal(dim(sim$population), c(p_me, p_me))
  expect_named(sim$model_error,
               c("method", "target_rmsea", "target_cfi", "rmsea", "cfi", "df", "kappa"))
  # return_pop drops the data but keeps the population and model-error record.
  rp <- efa_simulate(Lambda = Lambda_me, Phi = Phi_me, target_rmsea = 0.05, return_pop = TRUE)
  expect_null(rp$data)
  expect_equal(dim(rp$population), c(p_me, p_me))
  expect_equal(rp$model_error$rmsea, 0.05, tolerance = 1e-6)
})

test_that("print.efa_simulated summarises the object", {
  skip_on_cran()
  testthat::local_reproducible_output()
  # Scrub decimals: the achieved CFI is a BLAS-sensitive derived value, so pin the
  # wording and layout (and the integer dimensions) but not the computed decimals.
  sim <- efa_simulate(N = 100, Lambda = Lambda_me, Phi = Phi_me, target_rmsea = 0.05, seed = 1)
  expect_snapshot(print(sim), transform = scrub_num)
  rp <- efa_simulate(Lambda = Lambda_me, Phi = Phi_me, return_pop = TRUE)
  expect_snapshot(print(rp), transform = scrub_num)
})

rm(Lambda_pop, Phi_pop, R_pop, Lambda_vm, Phi_vm, R_vm,
   Lambda_me, Phi_me, p_me, q_me, df_me)
