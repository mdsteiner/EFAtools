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

  # the defaults: PCA eigenvalues for parallel analysis, the CAF for the Hull method
  expect_named(nf_grips$n_factors,
               c("CD", "EKC_BvA2017", "HULL_CAF", "MAP_TR2", "MAP_TR4", "NEST",
                 "PARALLEL_PCA"))

  expect_named(nf_grips$settings, c("criteria", "suitability", "N", "use",
                                    "n_factors_max", "N_pop", "N_samples", "alpha_CD",
                                    "alpha", "cor_method", "max_iter_CD", "n_fac_theor",
                                    "estimator", "method", "gof", "gof_used",
                                    "eigen_type_HULL", "decision_rule_HULL",
                                    "eigen_type_other", "n_factors", "n_datasets",
                                    "percent", "decision_rule", "n_datasets_nest",
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
    .retention_summary(list(MAP = crit(MAP_TR2 = 3, MAP_TR4 = 2),
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

test_that("the summary line renders its counts whatever options(scipen) is", {
  # An efa_retention object holds its counts as doubles, and a double is rendered in
  # scientific notation under a negative options(scipen) ("3e+00" instead of "3"). Both
  # branches of the line carry counts, so both are pinned.
  crit <- function(...) list(n_factors = c(...))
  withr::local_options(scipen = -5)

  expect_equal(.retention_summary(list(A = crit(A = 3), B = crit(B = 3))),
               "2 suggestions from 2 criteria, all suggesting 3 factors.")
  expect_equal(
    .retention_summary(list(A = crit(A = 2), B = crit(B = 2), C = crit(C = 4))),
    "3 suggestions from 3 criteria, ranging from 2 to 4 factors (most common: 2).")
})

test_that("settings record the Hull goodness-of-fit indices that actually ran", {
  skip_if_not_slow()
  cormat <- test_models$baseline$cormat

  # PAF supports only the CAF, so the Hull method reduces the requested set to it
  expect_message(
    r_paf <- efa_retain(cormat, N = 500, suitability = FALSE, criteria = "HULL",
                        estimator = "PAF", gof = c("CAF", "CFI", "RMSEA")),
    class = "efa_hull_gof_caf")
  expect_equal(r_paf$settings$gof, c("CAF", "CFI", "RMSEA"))
  expect_equal(r_paf$settings$gof_used, "CAF")
  expect_named(r_paf$n_factors, "HULL_CAF")

  # with an estimator that supports all three, requested and used coincide
  r_ml <- efa_retain(cormat, N = 500, suitability = FALSE, criteria = "HULL",
                     estimator = "ML", gof = c("CAF", "CFI", "RMSEA"))
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
  expect_match(txt, "KMO", fixed = TRUE)
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
  expect_false(grepl("Sampling adequacy", txt, fixed = TRUE))
})

# The classes of the package warnings that `expr` gives, one element per warning. Warnings of
# other packages (e.g. on loading one) are left out, so that they cannot fail a test.
warning_classes <- function(expr) {
  cls <- character(0)
  withCallingHandlers(expr, warning = function(w) {
    cls <<- c(cls, class(w)[1])
    invokeRestart("muffleWarning")
  })
  cls[startsWith(cls, "efa_")]
}

test_that("missing N skips Bartlett's test but still runs N-free criteria", {
  # A correlation matrix without N: Bartlett's test needs N and is not computed (the report
  # says so), while KMO and the N-free retention criteria still run. KGC with EFA
  # eigenvalues needs no N either.
  expect_length(warning_classes(
    nf_noN <- efa_retain(test_models$baseline$cormat,
                         criteria = c("KGC", "MAP", "SCREE"),
                         eigen_type_other = c("PCA", "SMC", "EFA"))
  ), 0L)
  expect_null(nf_noN$suitability$bartlett)
  expect_s3_class(nf_noN$suitability$kmo, "efa_kmo")
  expect_named(nf_noN$outputs, c("KGC", "MAP", "SCREE"))
  expect_null(nf_noN$not_run)
})

test_that("criteria that need N are skipped, with a warning also by default", {
  cormat <- test_models$baseline$cormat
  needs_n <- c("EKC", "HULL", "NEST", "PARALLEL", "SMT")
  expect_setequal(
    names(Filter(function(e) isTRUE(e$needs_N), .retention_registry)), needs_n)

  # each is skipped before it runs, not reported as a failure
  expect_identical(
    warning_classes(nf <- efa_retain(cormat, suitability = FALSE,
                                     criteria = c(needs_n, "MAP"))),
    "efa_criterion_skipped")
  expect_named(nf$not_run, needs_n)
  expect_named(nf$outputs, "MAP")

  # a missing N is most likely a mistake, so the default criteria report it too; CD, which
  # cannot use a correlation matrix at all, is skipped without a warning
  expect_identical(warning_classes(nf_def <- efa_retain(cormat, suitability = FALSE)),
                   "efa_criterion_skipped")
  expect_named(nf_def$not_run, c("CD", "EKC", "HULL", "NEST", "PARALLEL"))
})

test_that("default criteria that cannot use the input are skipped without a warning", {
  cormat <- test_models$baseline$cormat
  set.seed(1)
  expect_length(warning_classes(
    nf <- efa_retain(cormat, N = 500, n_datasets = 20, n_datasets_nest = 20)), 0L)
  expect_named(nf$not_run, "CD")
  expect_named(nf$outputs, c("EKC", "HULL", "MAP", "NEST", "PARALLEL"))

  # NULL selects the default criteria, so it skips silently too
  expect_length(warning_classes(
    nf_null <- efa_retain(cormat, N = 500, criteria = NULL, suitability = FALSE,
                          n_datasets = 20, n_datasets_nest = 20)), 0L)
  expect_named(nf_null$not_run, "CD")
})

test_that("with polychoric correlations, CD and Bartlett's test use Pearson correlations", {
  # HULL, NEST, and PARALLEL compare against continuous reference data and are skipped
  # without a warning; CD resamples the observed values, so it runs on the Pearson
  # correlations of the same data, as in the original method (Ruscio & Roche, 2012)
  set.seed(1)
  expect_length(warning_classes(
    nf <- efa_retain(GRiPS_raw, cor_method = "poly", N_pop = 500, N_samples = 20)), 0L)
  expect_named(nf$outputs, c("CD", "EKC", "MAP"))
  expect_named(nf$not_run, c("HULL", "NEST", "PARALLEL"))
  expect_identical(nf$outputs$CD$settings$cor_method, "pearson")
  expect_false(is.null(nf$outputs$CD$subtitle))
  set.seed(1)
  expect_equal(nf$outputs$CD$results,
               efa_cd(GRiPS_raw, N_pop = 500, N_samples = 20)$results)

  # Bartlett's test uses the Pearson correlations of the raw data, the KMO the polychoric ones
  expect_identical(nf$suitability$bartlett$settings$cor_method, "pearson")
  expect_equal(nf$suitability$bartlett$chisq,
               suppressMessages(efa_bartlett(GRiPS_raw))$chisq)
  expect_identical(nf$suitability$kmo$settings$cor_method, "poly")

  # ordered factors are coded by their levels, as for the polychoric correlations
  ord <- as.data.frame(lapply(as.data.frame(GRiPS_raw), ordered))
  set.seed(1)
  expect_length(warning_classes(
    nf_ord <- efa_retain(ord, cor_method = "poly", criteria = c("CD", "MAP"),
                         N_pop = 500, N_samples = 20)), 0L)
  expect_equal(nf_ord$suitability$bartlett$chisq, nf$suitability$bartlett$chisq)
  expect_equal(nf_ord$outputs$CD$results, nf$outputs$CD$results)

  # without raw data there is no Pearson matrix: the matrix supplied is tested
  expect_identical(
    efa_retain(stats::cor(GRiPS_raw), N = 810, cor_method = "poly",
               criteria = "MAP")$suitability$bartlett$settings$cor_method,
    "poly")

  # a criterion that you request is still reported
  expect_warning(efa_retain(GRiPS_raw, cor_method = "poly", suitability = FALSE,
                            criteria = c("MAP", "SMT")),
                 class = "efa_criterion_skipped")

  # pairwise-complete ordinal data with inconsistent pairs: both the polychoric and the
  # Pearson matrix need smoothing, but only the analysed (polychoric) one is reported
  set.seed(3)
  n <- 150
  cut5 <- function(v) as.integer(cut(v, c(-Inf, -1, -.3, .3, 1, Inf)))
  a <- rnorm(n); b <- rnorm(n); s <- rnorm(n)
  inc <- data.frame(x1 = c(cut5(a), cut5(b), rep(NA, n)),
                    x2 = c(cut5(a + rnorm(n, sd = .2)), rep(NA, n), cut5(s)),
                    x3 = c(rep(NA, n), cut5(b + rnorm(n, sd = .2)),
                           cut5(-s + rnorm(n, sd = .2))))
  expect_identical(
    warning_classes(efa_retain(inc, cor_method = "poly", criteria = "MAP")),
    "efa_cor_smoothed")
})

test_that("skipped and failed criteria give at most one warning of each kind", {
  # Without N, CD (needs raw data), EKC, and SMT (need N) are skipped; KGC fails, as its ML
  # fit has no starting values; MAP runs
  cls <- warning_classes(
    nf <- efa_retain(test_models$baseline$cormat, suitability = FALSE,
                     criteria = c("CD", "EKC", "KGC", "MAP", "SMT"),
                     eigen_type_other = "EFA",
                     estimate_control = estimate_control(start_method = NA)))
  expect_setequal(cls, c("efa_criterion_skipped", "efa_criterion_failed"))
  expect_length(cls, 2L)
  expect_named(nf$not_run, c("CD", "EKC", "KGC", "SMT"))
  expect_named(nf$outputs, "MAP")
  expect_named(nf$n_factors, c("MAP_TR2", "MAP_TR4"))
  expect_match(paste(format(nf), collapse = "\n"), "Criteria not run", fixed = TRUE)
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

test_that("an all-failed run is a hard error with no warning before it", {
  # CD is the only requested criterion and is skipped on a correlation matrix
  expect_length(warning_classes(expect_error(
    efa_retain(test_models$baseline$cormat, suitability = FALSE, criteria = "CD"),
    class = "efa_no_criteria"
  )), 0L)

  # when a criterion failed with an error rather than being skipped, its condition is
  # chained onto the abort, so the cause is reported with it
  expect_length(warning_classes(err <- tryCatch(
    efa_retain(test_models$baseline$cormat, N = 500, suitability = FALSE,
               criteria = "SMT", estimate_control = estimate_control(start_method = NA)),
    error = function(e) e
  )), 0L)
  expect_s3_class(err, "efa_no_criteria")
  expect_s3_class(err$parent, "efa_ml_start_missing")
})

test_that("a failed Hull method reports its own error", {
  # too few indicators: the abort chains the error of the Hull method itself
  cnd <- expect_error(
    efa_retain(test_models$baseline$cormat[1:5, 1:5], N = 500, criteria = "HULL",
               suitability = FALSE),
    class = "efa_no_criteria")
  expect_s3_class(cnd$parent, "efa_hull_min_indicators")

  # with the default criteria, the reason is the headline of the error that efa_hull()
  # gives for the same data
  set.seed(1)
  expect_identical(
    warning_classes(nf <- efa_retain(GRiPS_raw[, 1:5], N_pop = 500, N_samples = 20,
                                     n_datasets = 20, n_datasets_nest = 20)),
    "efa_criterion_failed")
  hull_cnd <- expect_error(efa_hull(GRiPS_raw[, 1:5]), class = "efa_hull_min_indicators")
  expect_identical(nf$not_run[["HULL"]], unname(hull_cnd$message))
})

test_that("invalid arguments abort before any criterion runs", {
  cormat <- test_models$baseline$cormat
  bad <- list(list(n_datasets = 0), list(N_pop = -1), list(N_pop = 1),
              list(N_samples = 0), list(alpha_CD = 2), list(max_iter_CD = 0),
              list(n_factors_max = 0), list(n_fac_theor = -1), list(n_factors = 0),
              list(percent = 101), list(n_datasets_nest = 0), list(alpha_nest = -1),
              list(show_progress = NA), list(suitability = "yes"), list(N = 0))
  for (b in bad) {
    expect_error(do.call(efa_retain, c(list(cormat, criteria = "MAP"), b)),
                 class = "efa_invalid_argument", info = names(b))
  }
  # also when the criterion that uses the argument is requested: an error, not a warning
  expect_length(warning_classes(expect_error(
    efa_retain(cormat, N = 500, criteria = c("MAP", "PARALLEL"), n_datasets = 0),
    class = "efa_invalid_argument")), 0L)
})

test_that("parallel analysis and the Hull method use PCA eigenvalues and the CAF", {
  cormat <- test_models$baseline$cormat
  # CAF works with every estimator, so PAF gives no message about the indices
  set.seed(1)
  expect_no_message(
    nf <- efa_retain(cormat, N = 500, suitability = FALSE, criteria = c("HULL", "PARALLEL"),
                     estimator = "PAF", n_datasets = 20),
    class = "efa_hull_gof_caf")
  expect_named(nf$n_factors, c("HULL_CAF", "PARALLEL_PCA"))
  expect_identical(nf$settings$gof_used, "CAF")
  expect_identical(nf$outputs$HULL$settings$eigen_type, "PCA")
})

test_that("the Hull method reuses the simulated eigenvalues of the parallel analysis", {
  # With this seed, the mean rule suggests 4 factors and the 95th percentile 3, so the
  # two rules give different upper bounds (5 and 4). Each eigenvalue is at least 0.5%
  # away from the reference value that decides it.
  cormat <- WJIV_ages_3_5$cormat
  calls <- 0L
  pa <- efa_parallel
  local_mocked_bindings(efa_parallel = function(...) {
    calls <<- calls + 1L
    pa(...)
  })

  set.seed(1)
  nf <- efa_retain(cormat, N = 300, suitability = FALSE, criteria = c("HULL", "PARALLEL"),
                   n_datasets = 20)
  expect_identical(calls, 1L)
  # the upper bound is the suggestion of the percentile rule plus one, and at least 3,
  # read off the same reference eigenvalues that the PARALLEL criterion compares with
  rec <- nf$outputs$PARALLEL$results[[1]]
  n_pct <- which(rec$y <= rec$references[["95 Percentile"]])[1] - 1
  expect_lt(n_pct, nf$n_factors[["PARALLEL_PCA"]])
  expect_equal(nf$outputs$HULL$settings$n_fac_max, max(n_pct + 1, 3))
  # the parallel analysis ran first, but the outputs keep the registry order
  expect_named(nf$outputs, c("HULL", "PARALLEL"))
  # efa_hull() gives the same result from the same reference draws
  set.seed(1)
  expect_equal(nf$outputs$HULL, efa_hull(cormat, N = 300, n_datasets = 20))

  # with the same rule for both, the Hull method takes the suggestion of PARALLEL
  calls <- 0L
  set.seed(1)
  nf_means <- efa_retain(cormat, N = 300, suitability = FALSE,
                         criteria = c("HULL", "PARALLEL"), n_datasets = 20,
                         decision_rule_HULL = "means")
  expect_identical(calls, 1L)
  expect_equal(nf_means$outputs$HULL$settings$n_fac_max,
               max(nf_means$n_factors[["PARALLEL_PCA"]] + 1, 3))

  # without a shared eigenvalue type, the Hull method runs its own parallel analysis
  calls <- 0L
  efa_retain(test_models$baseline$cormat, N = 500, suitability = FALSE,
             criteria = c("HULL", "PARALLEL"), n_datasets = 20, eigen_type_HULL = "SMC")
  expect_identical(calls, 2L)
})

test_that("alpha is deprecated in favour of alpha_CD", {
  cormat <- test_models$baseline$cormat
  expect_warning(
    nf <- efa_retain(cormat, criteria = "MAP", suitability = FALSE, alpha = .2),
    class = "lifecycle_warning_deprecated")
  expect_identical(nf$settings$alpha_CD, .2)
  # the former settings key stays as an alias
  expect_identical(nf$settings$alpha, .2)
  expect_error(efa_retain(cormat, criteria = "MAP", alpha = .2, alpha_CD = .1),
               class = "efa_renamed_arg")

  # alpha_CD reaches the comparison data
  set.seed(1)
  nf_cd <- efa_retain(GRiPS_raw, criteria = "CD", suitability = FALSE, alpha_CD = .1,
                      N_pop = 500, N_samples = 20)
  expect_identical(nf_cd$outputs$CD$settings$alpha, .1)
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
    expect_warning(efa_retain(burt, N = 170, criteria = c("PARALLEL", "EKC"),
                              eigen_type_other = "SMC"),
                   class = "efa_parallel_no_crossing"),
    class = "efa_cor_smoothed")
})

if (is_slow_test()) rm(nf_grips)
rm(burt)
