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

test_that("print.SCREE output is stable", {
  local_reproducible_output()

  scree <- SCREE(test_models$baseline$cormat)
  expect_snapshot(print(scree, plot = FALSE))

  scree_smc <- SCREE(test_models$baseline$cormat, eigen_type = "SMC")
  expect_snapshot(print(scree_smc, plot = FALSE))
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
