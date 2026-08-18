# Captured-output baselines for the print methods, recorded with
# local_reproducible_output() so the text is plain (no ANSI colour) and fixed
# width. Only methods whose printed output is integer factor counts or fixed
# text are snapshotted; methods that print computed decimals (e.g. KMO,
# Bartlett's chi-square) are omitted, since those digits can drift across BLAS
# and platforms and would make the baselines flaky.

test_that("print.efa_retention output is stable for SMT", {
  skip_if_not_slow()
  local_reproducible_output()

  smt <- SMT(test_models$baseline$cormat, N = 500)
  expect_snapshot(print(smt))

  # Identity matrix: zero-factor solution.
  smt_id <- SMT(diag(5), N = 500)
  expect_snapshot(print(smt_id))
})

test_that("print.efa_retention output is stable for SCREE", {
  skip_if_not_slow()
  local_reproducible_output()

  scree <- SCREE(test_models$baseline$cormat)
  expect_snapshot(print(scree))

  scree_smc <- SCREE(test_models$baseline$cormat, eigen_type = "SMC")
  expect_snapshot(print(scree_smc))
})

test_that("print.efa_retention output is stable for CD", {
  skip_if_not_slow()
  local_reproducible_output()

  # CD is stochastic; seed and use a small simulation for a stable, fast snapshot
  set.seed(123)
  cd <- CD(GRiPS_raw, N_pop = 1000, N_samples = 100)
  expect_snapshot(print(cd))
})

test_that("print.efa_retention output is stable for PARALLEL", {
  skip_if_not_slow()
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
  skip_if_not_slow()
  local_reproducible_output()

  kgc <- KGC(test_models$baseline$cormat)
  expect_snapshot(print(kgc))

  kgc_smc <- KGC(test_models$baseline$cormat, eigen_type = "SMC")
  expect_snapshot(print(kgc_smc))
})

test_that("print.efa_retention output is stable for NEST", {
  skip_if_not_slow()
  local_reproducible_output()

  nest <- NEST(test_models$baseline$cormat, N = 500)
  expect_snapshot(print(nest))
})

test_that("print.efa_retention output is stable for EKC", {
  skip_if_not_slow()
  local_reproducible_output()

  ekc <- EKC(test_models$baseline$cormat, N = 500)
  expect_snapshot(print(ekc))
})

test_that("print.efa_retention output is stable for HULL", {
  skip_if_not_slow()
  local_reproducible_output()

  hull <- HULL(test_models$baseline$cormat, N = 500, method = "ML")
  expect_snapshot(print(hull))

  hull_paf <- suppressMessages(HULL(test_models$baseline$cormat, N = 500))
  expect_snapshot(print(hull_paf))
})

test_that("print.efa_retention output is stable for MAP", {
  skip_if_not_slow()
  local_reproducible_output()

  map <- MAP(test_models$baseline$cormat)
  expect_snapshot(print(map))
})

test_that("print.efa_retain output is stable", {
  skip_if_not_slow()
  local_reproducible_output()

  # deterministic criteria only (no simulations); suitability = FALSE keeps
  # BLAS-sensitive decimals (chi-square, KMO) out of the snapshot
  nf <- efa_retain(test_models$baseline$cormat, N = 500, suitability = FALSE,
                   criteria = c("EKC", "KGC", "MAP", "SMT"),
                   eigen_type_other = c("PCA", "SMC"))
  expect_snapshot(print(nf))
})

test_that("every print method is exactly cat(format(x), sep = \"\\n\")", {
  skip_on_cran()
  local_reproducible_output()

  cm <- test_models$baseline$cormat
  fit <- efa_fit(cm, n_factors = 3, N = 500, estimator = "PAF", rotation = "promax")
  fit_ml <- efa_fit(cm, n_factors = 3, N = 500, estimator = "ML", rotation = "promax")
  sl <- efa_schmid_leiman(fit, estimator = "PAF")
  # correlation-matrix input, so efa_scores() notes that it returns weights only
  scores <- suppressMessages(efa_scores(cm, f = fit))

  # One object per printable class whose fixture is cheap to build. The three that are
  # not - efa_average, efa_group, and efa_mi - assert the same identity alongside their
  # own fixtures. The two loading tables print a trailing blank line for console spacing,
  # which is the only documented departure from the rule.
  fixtures <- list(
    efa                  = fit,
    summary.efa          = summary(fit),
    efa_loadings         = fit$rot_loadings,
    efa_sl_loadings      = sl$sl,
    efa_schmid_leiman    = sl,
    OMEGA                = OMEGA(sl, type = "EFAtools",
                                 factor_corres = sl$sl[, c("F1", "F2", "F3")] >= .2),
    efa_reliability      = efa_reliability(fit),
    efa_scores           = scores,
    summary.efa_scores   = summary(scores),
    efa_compare          = efa_compare(fit$unrot_loadings, fit_ml$unrot_loadings),
    efa_retention        = efa_ekc(cm, N = 500),
    efa_retain           = efa_retain(cm, N = 500, suitability = FALSE,
                                      criteria = c("EKC", "KGC")),
    efa_screen           = efa_screen(cm, N = 500),
    efa_kmo              = efa_kmo(cm),
    efa_bartlett         = efa_bartlett(cm, N = 500),
    efa_power            = efa_power(df = 100, N = 200),
    efa_simulated        = efa_simulate(N = 100, Lambda = population_models$loadings$baseline,
                                        seed = 42),
    efa_estimate_control = estimate_control(type = "SPSS"),
    efa_rotate_control   = rotate_control(type = "SPSS")
  )
  # trailing blank line, see above
  spacer <- c("efa_loadings", "efa_sl_loadings")

  for (nm in names(fixtures)) {
    obj <- fixtures[[nm]]
    expected <- format(obj)
    if (nm %in% spacer) expected <- c(expected, "")
    expect_identical(utils::capture.output(print(obj)), expected, info = nm)
  }
})

test_that("print.efa_retain summarises the criteria's suggestions", {
  local_reproducible_output()

  # A cheap deterministic set whose criteria disagree, pinning the summary line under
  # the section rule. The modal value is counted one vote per criterion: the empirical
  # Kaiser criterion and the Kaiser-Guttman criterion both suggest 3, while the
  # internally tied minimum average partial (1 and 3) is undecided and casts no vote
  # rather than casting two.
  nf <- efa_retain(test_models$baseline$cormat, N = 500, suitability = FALSE,
                   criteria = c("EKC", "KGC", "MAP"), eigen_type_other = "PCA")
  expect_snapshot(print(nf))
})
