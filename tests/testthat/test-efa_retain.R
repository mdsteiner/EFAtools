# efa_retain aggregates the heavy retention criteria (CD's bootstrap simulation,
# PARALLEL's eigen Monte Carlo, HULL's per-factor refits), so the file-top
# `nf_grips` fixture dominates this file (~28s). Skipped by default; opt in with
# `Sys.setenv(EFATOOLS_TEST_SLOW = "true")`. See helper-slow.R.
if (is_slow_test()) {
set.seed(42)
nf_grips <- suppressMessages(suppressWarnings(efa_retain(GRiPS_raw)))
}  # is_slow_test()

test_that("output class and dimensions are correct", {
  skip_if_not_slow()
  expect_s3_class(nf_grips, c("efa_retain", "N_FACTORS"), exact = TRUE)
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

  expect_s3_class(nf_grips$suitability$bartlett, "efa_bartlett")
  expect_s3_class(nf_grips$suitability$kmo, "efa_kmo")

  expect_named(nf_grips$n_factors,
               c("CD", "EKC_BvA2017", "HULL_CAF", "HULL_CFI", "HULL_RMSEA",
                 "MAP_TR2", "MAP_TR4", "NEST", "PARALLEL_SMC"))

  expect_named(nf_grips$settings, c("criteria", "suitability", "N", "use",
                                    "n_factors_max", "N_pop", "N_samples", "alpha",
                                    "cor_method", "max_iter_CD", "n_fac_theor",
                                    "estimator", "method", "gof", "gof_used",
                                    "eigen_type_HULL",
                                    "eigen_type_other", "n_factors", "n_datasets",
                                    "percent", "decision_rule",
                                    "ekc_type", "n_datasets_nest",
                                    "alpha_nest", "estimate_control"))
})

test_that("the summary line counts one vote per criterion", {
  # `outputs` elements are efa_retention objects; only their $n_factors is read
  crit <- function(...) list(n_factors = c(...))

  # nothing to summarise: no criterion determined a number, or a lone suggestion is
  # already carried by the bullet printed below the line
  expect_null(.retention_summary(list(A = crit(A = NA_real_))))
  expect_null(.retention_summary(list(A = crit(A = 3))))

  expect_equal(.retention_summary(list(A = crit(A = 3), B = crit(B = 3))),
               "2 suggestions from 2 criteria, all suggesting 3 factors.")

  # a three-variant criterion casts one vote, so it cannot outvote two single-variant
  # criteria that agree
  expect_equal(
    .retention_summary(list(HULL = crit(HULL_CAF = 4, HULL_CFI = 4, HULL_RMSEA = 4),
                            EKC = crit(EKC = 2), MAP = crit(MAP = 2))),
    "5 suggestions from 3 criteria, ranging from 2 to 4 factors (most common: 2).")

  # that one vote is the criterion's own modal variant
  expect_equal(
    .retention_summary(list(HULL = crit(HULL_CAF = 2, HULL_CFI = 2, HULL_RMSEA = 5),
                            EKC = crit(EKC = 2))),
    "4 suggestions from 2 criteria, ranging from 2 to 5 factors (most common: 2).")

  # a criterion whose variants tie has no modal value and abstains: here 3 is the mode
  # of the per-variant vector, but only SMT decides, so no value wins more than one vote
  expect_equal(
    .retention_summary(list(EKC = crit(EKC_BvA2017 = 3, EKC_AM2019 = 2),
                            KGC = crit(KGC_PCA = 3, KGC_SMC = 1),
                            SMT = crit(SMT_chi = 3, SMT_RMSEA = 3))),
    "6 suggestions from 3 criteria, ranging from 1 to 3 factors.")

  # every deciding criterion picked a different number: the range speaks for itself
  expect_equal(
    .retention_summary(list(A = crit(A = 1), B = crit(B = 2), C = crit(C = 3))),
    "3 suggestions from 3 criteria, ranging from 1 to 3 factors.")

  # tied winners among the voters are both reported
  expect_equal(
    .retention_summary(list(A = crit(A = 2), B = crit(B = 2),
                            C = crit(C = 4), D = crit(D = 4))),
    "4 suggestions from 4 criteria, ranging from 2 to 4 factors (most common: 2 and 4).")

  # a variant that could not determine a number contributes nothing to either count
  expect_equal(
    .retention_summary(list(A = crit(A_x = 2, A_y = NA_real_), B = crit(B = 3))),
    "2 suggestions from 2 criteria, ranging from 2 to 3 factors.")
})

test_that("settings record the Hull goodness-of-fit indices that actually ran", {
  skip_if_not_slow()
  cormat <- test_models$baseline$cormat

  # PAF supports only the CAF, so the Hull method reduces the requested set to it
  r_paf <- suppressMessages(
    efa_retain(cormat, N = 500, suitability = FALSE, criteria = "HULL", estimator = "PAF"))
  expect_equal(r_paf$settings$gof, c("CAF", "CFI", "RMSEA"))
  expect_equal(r_paf$settings$gof_used, "CAF")
  expect_named(r_paf$n_factors, "HULL_CAF")

  # with an estimator that supports all three, requested and used coincide
  r_ml <- efa_retain(cormat, N = 500, suitability = FALSE, criteria = "HULL",
                     estimator = "ML")
  expect_equal(r_ml$settings$gof_used, r_ml$settings$gof)

  # no HULL result to read it off
  r_no_hull <- efa_retain(cormat, N = 500, suitability = FALSE, criteria = "EKC")
  expect_identical(r_no_hull$settings$gof_used, NA_character_)
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

  # one ggplot per plottable criterion that was run (MAP has no plot)
  p <- plot(nf_grips)
  expect_named(p, c("CD", "EKC", "HULL", "NEST", "PARALLEL"))
  for (p_i in p) {
    expect_s3_class(p_i, "ggplot")
  }
})

test_that("the aggregate report carries each criterion's own context line", {
  # The lines are compared verbatim, so pin the console width: the subtitle is longer
  # than a narrow console and would otherwise be wrapped across two of them.
  local_reproducible_output()

  # The bullets are keyed by eigenvalue type / fit index, so the group has to repeat
  # the criterion's subtitle -- otherwise the aggregate report shows a bare "PCA: 3"
  # with nothing saying what that key is. Criteria without a subtitle add no line.
  nf <- efa_retain(test_models$baseline$cormat, N = 500, suitability = FALSE,
                   criteria = c("EKC", "KGC"), eigen_type_other = c("PCA", "SMC"))
  txt <- format(nf)

  expect_true(any(txt == nf$outputs$KGC$subtitle))
  # and it sits between the criterion heading and that criterion's bullets
  expect_equal(which(txt == nf$outputs$KGC$subtitle),
               which(txt == "Kaiser-Guttman criterion") + 1)
  # EKC has no subtitle, so its heading is followed straight by its bullet
  expect_null(nf$outputs$EKC$subtitle)
})

test_that("format.efa_retain is the source of truth and honours the colour state", {
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
  nf_scree <- efa_retain(test_models$baseline$cormat, suitability = FALSE,
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
  # A correlation matrix without N: Bartlett's test needs N and is skipped with a
  # classed warning, while KMO and the N-free retention criteria still run -
  # instead of the suitability check aborting the whole call.
  expect_warning(
    nf_noN <- efa_retain(test_models$baseline$cormat,
                        criteria = c("KGC", "MAP", "SCREE")),
    class = "efa_suitability_no_n"
  )
  expect_null(nf_noN$suitability$bartlett)
  expect_s3_class(nf_noN$suitability$kmo, "efa_kmo")
  expect_named(nf_noN$outputs, c("KGC", "MAP", "SCREE"), ignore.order = TRUE)
})

test_that("a failing criterion warns, is reported, and the others still run", {
  # EKC needs N, which is missing for a correlation matrix; MAP still runs
  expect_warning(
    nf_fail <- efa_retain(test_models$baseline$cormat, suitability = FALSE,
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
  # CD needs raw data and is skipped on a correlation matrix; pairing it with the cheap MAP
  # exercises the skip path without the simulating default criteria, so this runs on CRAN.
  expect_warning(
    nf_skip <- efa_retain(test_models$baseline$cormat, N = 500,
                          criteria = c("CD", "MAP")),
    class = "efa_criterion_skipped"
  )
  # CD needs raw data and is skipped on a correlation matrix
  expect_named(nf_skip$not_run, "CD")
  expect_false("CD" %in% names(nf_skip$outputs))
})

test_that("an all-failed run is a hard error", {
  # CD is the only requested criterion and is skipped on a correlation matrix
  expect_error(
    suppressWarnings(efa_retain(test_models$baseline$cormat, suitability = FALSE,
                               criteria = "CD")),
    class = "efa_no_criteria"
  )

  # when a criterion failed with an error rather than being skipped, its condition is
  # chained onto the abort, so the cause is reported with it
  err <- tryCatch(
    suppressWarnings(efa_retain(test_models$baseline$cormat, suitability = FALSE,
                                criteria = "EKC")),
    error = function(e) e
  )
  expect_s3_class(err, "efa_no_criteria")
  expect_s3_class(err$parent, "efa_n_required")
})

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
  expect_error(efa_retain(1:10), class = "efa_input_not_matrix")
  expect_warning(efa_retain(GRiPS_raw, N = 10, criteria = "MAP",
                           suitability = FALSE), class = "efa_n_from_data")
  expect_error(efa_retain(sing_raw), class = "efa_cor_singular")
  # CD is skipped on a correlation matrix; a small criteria set keeps the skip-warning check
  # off the simulating default criteria so it runs on CRAN.
  expect_warning(efa_retain(test_models$baseline$cormat, N = 500,
                            criteria = c("CD", "MAP")), class = "efa_criterion_skipped")
  # burt is near-singular: it is smoothed, and PARALLEL finds every real SMC
  # eigenvalue above its reference (no crossing), so both warnings are expected.
  # PARALLEL's verdict comes out of a simulation, so seed it rather than inheriting
  # whatever state the preceding tests left (which differs between the default and the
  # slow gate, because the gated bodies return before drawing).
  set.seed(42)
  expect_warning(
    expect_warning(efa_retain(burt, N = 170, criteria = c("PARALLEL", "EKC")),
                   class = "efa_parallel_no_crossing"),
    class = "efa_cor_smoothed")
})

if (is_slow_test()) rm(nf_grips)
rm(burt)
