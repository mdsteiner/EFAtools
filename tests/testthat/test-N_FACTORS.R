# N_FACTORS aggregates the heavy retention criteria (CD's bootstrap simulation,
# PARALLEL's eigen Monte Carlo, HULL's per-factor refits), so the file-top
# `nf_grips` fixture dominates this file (~28s). Skipped by default; opt in with
# `Sys.setenv(EFATOOLS_TEST_SLOW = "true")`. See helper-slow.R.
if (is_slow_test()) {
set.seed(42)
nf_grips <- suppressMessages(suppressWarnings(N_FACTORS(GRiPS_raw)))
}  # is_slow_test()

test_that("output class and dimensions are correct", {
  skip_if_not_slow()
  expect_s3_class(nf_grips, "N_FACTORS")
  expect_named(nf_grips, c("suitability", "outputs", "n_factors", "not_run",
                           "settings"))
  expect_type(nf_grips$outputs, "list")
  expect_type(nf_grips$settings, "list")
  expect_type(nf_grips$n_factors, "double")

  # default criteria, in registry order, each a unified retention object
  expect_named(nf_grips$outputs, c("CD", "EKC", "HULL", "MAP", "NEST",
                                   "PARALLEL"))
  for (out in nf_grips$outputs) {
    expect_s3_class(out, "efa_retention")
  }

  # all requested criteria ran (GRiPS is raw data, so CD is not skipped)
  expect_null(nf_grips$not_run)

  expect_s3_class(nf_grips$suitability$bartlett, "BARTLETT")
  expect_s3_class(nf_grips$suitability$kmo, "KMO")

  expect_named(nf_grips$n_factors,
               c("CD", "EKC_BvA2017", "HULL_CAF", "HULL_CFI", "HULL_RMSEA",
                 "MAP_TR2", "MAP_TR4", "NEST", "PARALLEL_SMC"))

  expect_named(nf_grips$settings, c("criteria", "suitability", "N", "use",
                                    "n_factors_max", "N_pop", "N_samples", "alpha",
                                    "cor_method", "max_iter_CD", "n_fac_theor",
                                    "method", "gof", "eigen_type_HULL",
                                    "eigen_type_other", "n_factors", "n_datasets",
                                    "percent", "decision_rule",
                                    "ekc_type", "n_datasets_nest",
                                    "alpha_nest"))
})

test_that("print and plot work on the aggregate object", {
  skip_if_not_slow()
  # the suitability block is rendered, including the Bartlett verdict and the
  # KMO label ladder (GRiPS has a marvellous KMO and a significant Bartlett test)
  txt <- paste(format(nf_grips), collapse = "\n")
  expect_match(txt, "Bartlett", fixed = TRUE)
  expect_match(txt, "significant", fixed = TRUE)
  expect_match(txt, "Kaiser-Meyer-Olkin", fixed = TRUE)
  expect_match(txt, "marvellous", fixed = TRUE)

  # one ggplot per plottable criterion that was run (MAP/NEST have no plot)
  p <- plot(nf_grips)
  expect_named(p, c("CD", "EKC", "HULL", "PARALLEL"))
  for (p_i in p) {
    expect_s3_class(p_i, "ggplot")
  }
})

test_that("format.N_FACTORS is the source of truth and honours the colour state", {
  skip_if_not_slow()
  # print() is exactly cat(format(x), sep = "\n"), so the two agree line for line.
  expect_identical(utils::capture.output(print(nf_grips)), format(nf_grips))

  old <- options(cli.num_colors = 256)
  on.exit(options(old), add = TRUE)

  # With colours on the report embeds ANSI ...
  expect_true(any(grepl("\033", format(nf_grips), fixed = TRUE)))

  # ... and with colours off it is plain.
  options(cli.num_colors = 1)
  expect_false(any(grepl("\033", format(nf_grips), fixed = TRUE)))
})

test_that("visual criteria and suitability = FALSE are handled", {
  skip_if_not_slow()
  nf_scree <- N_FACTORS(test_models$baseline$cormat, suitability = FALSE,
                        criteria = c("MAP", "SCREE"))

  expect_null(nf_scree$suitability)
  expect_named(nf_scree$outputs, c("MAP", "SCREE"))
  # the visual scree plot contributes no numeric suggestion
  expect_named(nf_scree$n_factors, c("MAP_TR2", "MAP_TR4"))

  txt <- paste(format(nf_scree), collapse = "\n")
  expect_match(txt, "no numeric suggestion", fixed = TRUE)
  expect_false(grepl("suitability of the data", txt, fixed = TRUE))
})

test_that("missing N skips Bartlett's test but still runs N-free criteria", {
  skip_if_not_slow()
  # A correlation matrix without N: Bartlett's test needs N and is skipped with a
  # classed warning, while KMO and the N-free retention criteria still run -
  # instead of the suitability check aborting the whole call.
  expect_warning(
    nf_noN <- N_FACTORS(test_models$baseline$cormat,
                        criteria = c("KGC", "MAP", "SCREE")),
    class = "efa_suitability_no_n"
  )
  expect_null(nf_noN$suitability$bartlett)
  expect_s3_class(nf_noN$suitability$kmo, "KMO")
  expect_named(nf_noN$outputs, c("KGC", "MAP", "SCREE"), ignore.order = TRUE)
})

test_that("a failing criterion warns, is reported, and the others still run", {
  skip_if_not_slow()
  # EKC needs N, which is missing for a correlation matrix; MAP still runs
  expect_warning(
    nf_fail <- N_FACTORS(test_models$baseline$cormat, suitability = FALSE,
                         criteria = c("EKC", "MAP")),
    class = "efa_criterion_failed"
  )
  expect_named(nf_fail$outputs, "MAP")
  expect_named(nf_fail$n_factors, c("MAP_TR2", "MAP_TR4"))
  expect_named(nf_fail$not_run, "EKC")
  expect_match(paste(format(nf_fail), collapse = "\n"), "could not be run",
               fixed = TRUE)
})

test_that("a skipped criterion is reported in not_run", {
  skip_if_not_slow()
  expect_warning(
    nf_skip <- N_FACTORS(test_models$baseline$cormat, N = 500),
    class = "efa_criterion_skipped"
  )
  # CD needs raw data and is skipped on a correlation matrix
  expect_named(nf_skip$not_run, "CD")
  expect_false("CD" %in% names(nf_skip$outputs))
})

test_that("an all-failed run is a hard error", {
  skip_if_not_slow()
  # CD is the only requested criterion and is skipped on a correlation matrix
  expect_error(
    suppressWarnings(N_FACTORS(test_models$baseline$cormat, suitability = FALSE,
                               criteria = "CD")),
    class = "efa_no_criteria"
  )
})

x <- rnorm(100)
y <- rnorm(100)
z <- x + y



burt <- matrix(c(1.00,  0.83,  0.81,  0.80,   0.71, 0.70, 0.54, 0.53,  0.59,  0.24, 0.13,
                 0.83,  1.00,  0.87,  0.62,   0.59, 0.44, 0.58, 0.44,  0.23,  0.45,  0.21,
                 0.81,  0.87,  1.00,  0.63,   0.37, 0.31, 0.30, 0.12,  0.33,  0.33,  0.36,
                 0.80,  0.62,  0.63,  1.00,   0.49, 0.54, 0.30, 0.28,  0.42,  0.29, -0.06,
                 0.71,  0.59,  0.37,  0.49,   1.00, 0.54, 0.34, 0.55,  0.40,  0.19, -0.10,
                 0.70,  0.44,  0.31,  0.54,   0.54, 1.00, 0.50, 0.51,  0.31,  0.11,  0.10,
                 0.54,  0.58,  0.30,  0.30,   0.34, 0.50, 1.00, 0.38,  0.29,  0.21,  0.08,
                 0.53,  0.44,  0.12,  0.28,   0.55, 0.51, 0.38, 1.00,  0.53,  0.10, -0.16,
                 0.59,  0.23,  0.33,  0.42,   0.40, 0.31, 0.29, 0.53,  1.00, -0.09, -0.10,
                 0.24,  0.45,  0.33,  0.29,   0.19, 0.11, 0.21, 0.10, -0.09,  1.00,  0.41,
                 0.13,  0.21,  0.36, -0.06,  -0.10, 0.10, 0.08, -0.16, -0.10, 0.41,  1.00),
               nrow = 11, ncol = 11)

test_that("errors etc. are thrown correctly", {
  skip_if_not_slow()
  expect_error(N_FACTORS(1:10), class = "efa_input_not_matrix")
  expect_warning(N_FACTORS(GRiPS_raw, N = 10, criteria = "MAP",
                           suitability = FALSE), class = "efa_n_from_data")
  expect_error(N_FACTORS(cbind(x, y, z, z + 1, y + 1, x + 1)), class = "efa_cor_singular")
  expect_warning(N_FACTORS(test_models$baseline$cormat, N = 500), class = "efa_criterion_skipped")
  # burt is near-singular: it is smoothed, and PARALLEL finds every real SMC
  # eigenvalue above its reference (no crossing), so both warnings are expected
  expect_warning(
    expect_warning(N_FACTORS(burt, N = 170, criteria = c("PARALLEL", "EKC")),
                   class = "efa_parallel_no_crossing"),
    class = "efa_cor_smoothed")
})

if (is_slow_test()) rm(nf_grips)
rm(x, y, z, burt)
