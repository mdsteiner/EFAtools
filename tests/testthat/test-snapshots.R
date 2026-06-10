# Captured-output baselines for the print methods, recorded with
# local_reproducible_output() so the text is plain (no ANSI colour) and fixed
# width. Only methods whose printed output is integer factor counts or fixed
# text are snapshotted; methods that print computed decimals (e.g. KMO,
# Bartlett's chi-square) are omitted, since those digits can drift across BLAS
# and platforms and would make the baselines flaky.

test_that("print.efa_retention output is stable for SMT", {
  local_reproducible_output()

  smt <- SMT(test_models$baseline$cormat, N = 500)
  expect_snapshot(print(smt))

  # Identity matrix: zero-factor solution.
  smt_id <- SMT(diag(5), N = 500)
  expect_snapshot(print(smt_id))
})

test_that("print.efa_retention output is stable for SCREE", {
  local_reproducible_output()

  scree <- SCREE(test_models$baseline$cormat)
  expect_snapshot(print(scree))

  scree_smc <- SCREE(test_models$baseline$cormat, eigen_type = "SMC")
  expect_snapshot(print(scree_smc))
})

test_that("print.efa_retention output is stable for CD", {
  local_reproducible_output()

  # CD is stochastic; seed and use a small simulation for a stable, fast snapshot
  set.seed(123)
  cd <- CD(GRiPS_raw, N_pop = 1000, N_samples = 100)
  expect_snapshot(print(cd))
})

test_that("print.efa_retention output is stable for PARALLEL", {
  local_reproducible_output()

  # seed the simulated eigenvalues so the printed suggestion cannot drift
  set.seed(123)
  pa <- PARALLEL(test_models$baseline$cormat, N = 500, eigen_type = "PCA")
  expect_snapshot(print(pa))

  # without real data there is no numeric suggestion (visual output only)
  pa_nodat <- PARALLEL(N = 20, n_vars = 5)
  expect_snapshot(print(pa_nodat))
})

test_that("print.efa_retention output is stable for KGC", {
  local_reproducible_output()

  kgc <- KGC(test_models$baseline$cormat)
  expect_snapshot(print(kgc))

  kgc_smc <- KGC(test_models$baseline$cormat, eigen_type = "SMC")
  expect_snapshot(print(kgc_smc))
})

test_that("print.efa_retention output is stable for NEST", {
  local_reproducible_output()

  nest <- NEST(test_models$baseline$cormat, N = 500)
  expect_snapshot(print(nest))
})

test_that("print.efa_retention output is stable for EKC", {
  local_reproducible_output()

  ekc <- EKC(test_models$baseline$cormat, N = 500)
  expect_snapshot(print(ekc))

  ekc_both <- EKC(test_models$baseline$cormat, N = 500,
                  type = c("BvA2017", "AM2019"))
  expect_snapshot(print(ekc_both))
})

test_that("print.efa_retention output is stable for HULL", {
  local_reproducible_output()

  hull <- HULL(test_models$baseline$cormat, N = 500, method = "ML")
  expect_snapshot(print(hull))

  hull_paf <- suppressMessages(HULL(test_models$baseline$cormat, N = 500))
  expect_snapshot(print(hull_paf))
})

test_that("print.efa_retention output is stable for MAP", {
  local_reproducible_output()

  map <- MAP(test_models$baseline$cormat)
  expect_snapshot(print(map))
})
