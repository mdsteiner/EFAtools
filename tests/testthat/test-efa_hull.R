# HULL refits the model at every candidate factor count for every (method, gof)
# combination, which makes the file-top fixture block the single heaviest in the
# suite (~60s). Skip the fixtures + the tests below by default; opt in with
# `Sys.setenv(EFATOOLS_TEST_SLOW = "true")` to run the full HULL coverage. See
# helper-slow.R for the convention.
if (is_slow_test()) {
hull_cor_paf <- suppressMessages(efa_hull(test_models$baseline$cormat, N = 500))
hull_cor_ml <- efa_hull(test_models$baseline$cormat, N = 500, estimator = "ML")
hull_cor_uls <- efa_hull(test_models$baseline$cormat, N = 500, estimator = "ULS")
hull_cor_uls_CFI <- efa_hull(test_models$baseline$cormat, N = 500, estimator = "ULS",
                         gof = "CFI")
hull_cor_ml_nf <- suppressWarnings(efa_hull(test_models$baseline$cormat, N = 500,
                                        estimator = "ML", n_fac_theor = 12))

hull_PCA <- efa_hull(test_models$baseline$cormat, N = 500, estimator = "ML",
                 eigen_type = "PCA")
hull_EFA <- efa_hull(test_models$baseline$cormat, N = 500, estimator = "ML",
                 eigen_type = "EFA")

hull_raw_paf <- suppressMessages(suppressWarnings(efa_hull(GRiPS_raw)))
hull_raw_ml <- suppressMessages(suppressWarnings(efa_hull(GRiPS_raw, estimator = "ML")))
hull_raw_uls <- suppressMessages(suppressWarnings(efa_hull(GRiPS_raw, estimator = "ULS")))
hull_raw_uls_CFI <- suppressMessages(suppressWarnings(efa_hull(GRiPS_raw,
                                                           estimator = "ULS",
                                                           gof = "CFI")))
hull_raw_ml_nf <- suppressMessages(suppressWarnings(efa_hull(GRiPS_raw, N = 500,
                                                         estimator = "ML",
                                                         n_fac_theor = 7)))
}  # is_slow_test()


test_that("output class and dimensions are correct", {
  skip_if_not_slow()
  objs <- list(hull_cor_paf, hull_cor_ml, hull_cor_uls, hull_cor_uls_CFI,
               hull_PCA, hull_EFA, hull_raw_paf, hull_raw_ml, hull_raw_uls,
               hull_raw_uls_CFI)
  for (obj in objs) {
    expect_s3_class(obj, "efa_retention")
    expect_length(obj, 9)
    expect_length(obj$settings, 9)
  }
})

test_that("n_fac_max is correctly specified", {
  skip_if_not_slow()
  expect_lte(hull_cor_paf$settings$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))
  expect_lte(hull_cor_ml$settings$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))
  expect_lte(hull_cor_uls$settings$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))
  expect_lte(hull_cor_uls_CFI$settings$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))
  expect_lte(hull_PCA$settings$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))
  expect_lte(hull_EFA$settings$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))

  expect_lte(hull_raw_paf$settings$n_fac_max,
             .det_max_factors(ncol(GRiPS_raw)))
  expect_lte(hull_raw_ml$settings$n_fac_max,
             .det_max_factors(ncol(GRiPS_raw)))
  expect_lte(hull_raw_uls$settings$n_fac_max,
             .det_max_factors(ncol(GRiPS_raw)))
  expect_lte(hull_raw_uls_CFI$settings$n_fac_max,
             .det_max_factors(ncol(GRiPS_raw)))

  expect_equal(hull_cor_ml_nf$settings$n_fac_max,
             .det_max_factors(ncol(test_models$baseline$cormat)))
  expect_equal(hull_raw_ml_nf$settings$n_fac_max,
               .det_max_factors(ncol(GRiPS_raw)))
})

test_that("the factor-search bound stays over-identified with the minimum number of indicators", {
  skip_if_not_slow()
  # With 6 indicators the largest over-identified model has 2 factors (a 3-factor
  # solution is just-identified, df = 0). The search bound J must be capped at that
  # maximum and never forced up to 3, so the Hull method cannot select a df = 0 model.
  R6 <- test_models$baseline$cormat[1:6, 1:6]
  out <- suppressWarnings(suppressMessages(efa_hull(R6, N = 500, estimator = "ULS")))
  expect_equal(out$settings$n_fac_max, .det_max_factors(6))
  expect_lte(max(out$n_factors, na.rm = TRUE), .det_max_factors(6))
})

test_that("records are correctly returned", {
  skip_if_not_slow()
  expect_named(hull_cor_paf$n_factors, "CAF")
  expect_equal(hull_cor_paf$results[[1]]$plot_type, "hull")
  checkmate::expect_numeric(hull_cor_paf$results[[1]]$y)
  checkmate::expect_numeric(hull_cor_paf$results[[1]]$x)
  checkmate::expect_logical(hull_cor_paf$results[[1]]$on_hull)

  expect_named(hull_cor_uls$n_factors, c("CAF", "CFI", "RMSEA"))
  expect_named(hull_cor_uls_CFI$n_factors, "CFI")

  expect_named(hull_raw_paf$n_factors, "CAF")
  expect_named(hull_raw_uls$n_factors, c("CAF", "CFI", "RMSEA"))
  expect_named(hull_raw_uls_CFI$n_factors, "CFI")
})

test_that("n_factors are correctly returned", {
  skip_if_not_slow()
  expect_equal(hull_cor_paf$n_factors[["CAF"]], 3)
  expect_false("CFI" %in% names(hull_cor_paf$n_factors))
  expect_false("RMSEA" %in% names(hull_cor_paf$n_factors))

  expect_equal(hull_cor_uls$n_factors[["CAF"]], 3)
  expect_equal(hull_cor_uls$n_factors[["CFI"]], 1)
  expect_equal(hull_cor_uls$n_factors[["RMSEA"]], 1)
  expect_equal(hull_cor_uls$n_fac_CAF, hull_cor_uls$n_factors[["CAF"]])
  expect_equal(hull_cor_uls$n_fac_CFI, hull_cor_uls$n_factors[["CFI"]])
  expect_equal(hull_cor_uls$n_fac_RMSEA, hull_cor_uls$n_factors[["RMSEA"]])

  expect_equal(hull_cor_uls_CFI$n_factors[["CFI"]], 1)
  expect_false("CAF" %in% names(hull_cor_uls_CFI$n_factors))
  expect_false("RMSEA" %in% names(hull_cor_uls_CFI$n_factors))
  expect_equal(hull_cor_uls_CFI$n_fac_CAF, NA_real_)
  expect_equal(hull_cor_uls_CFI$n_fac_CFI, hull_cor_uls_CFI$n_factors[["CFI"]])
  expect_equal(hull_cor_uls_CFI$n_fac_RMSEA, NA_real_)

  expect_equal(hull_raw_paf$n_factors[["CAF"]], 1)
  expect_false("CFI" %in% names(hull_raw_paf$n_factors))
  expect_false("RMSEA" %in% names(hull_raw_paf$n_factors))

  expect_equal(hull_raw_uls$n_factors[["CAF"]], 1)
  expect_equal(hull_raw_uls$n_factors[["CFI"]], 1)
  expect_equal(hull_raw_uls$n_factors[["RMSEA"]], 1)

  expect_equal(hull_raw_uls_CFI$n_factors[["CFI"]], 1)
  expect_false("CAF" %in% names(hull_raw_uls_CFI$n_factors))
  expect_false("RMSEA" %in% names(hull_raw_uls_CFI$n_factors))
})

test_that("the convex-hull elimination tests every interior triplet", {
  # Hand-built fit table (columns: nfactors, goodness-of-fit, df, st). The fit
  # values increase monotonically, so the boundary step removes nothing; the
  # geometry is set so the last interior solution (two factors) lies below the
  # line connecting its neighbours and must be dropped from the hull.
  s <- cbind(0:3, c(0.00, 0.80, 0.85, 1.00), c(10, 6, 3, 1), 0)
  out <- .hull_calc(s, J = 3, gof_t = "CAF")

  # the below-chord last interior solution is off the hull (NA st) ...
  expect_true(is.na(out$s_complete[3, 4]))
  # ... while the supporting interior solution remains on it
  expect_false(is.na(out$s_complete[2, 4]))
  expect_equal(out$retain, 1)

  # a last interior solution lying above its neighbours' chord is kept
  s2 <- cbind(0:3, c(0.00, 0.85, 0.95, 1.00), c(10, 6, 3, 1), 0)
  out2 <- .hull_calc(s2, J = 3, gof_t = "CAF")
  expect_false(is.na(out2$s_complete[3, 4]))

  # a fully convex fit curve collapses the hull below three solutions; the
  # max-fit fallback is used instead of indexing past the remaining solutions
  s3 <- cbind(0:4, c(0.00, 0.10, 0.25, 0.55, 1.00), c(10, 8, 6, 3, 1), 0)
  expect_warning(out3 <- .hull_calc(s3, J = 4, gof_t = "CAF"),
                 class = "efa_hull_few_solutions")
  expect_equal(out3$retain, 4)
})

test_that("a non-finite goodness-of-fit value is dropped with a classed warning", {
  # A model with an undefined fit value (e.g. an NA CFI for a Heywood or
  # near-singular solution) cannot lie on the hull; it must be excluded with a
  # classed warning instead of crashing the elimination comparisons.
  s <- cbind(0:4, c(0.00, 0.80, NA, 0.90, 1.00), c(10, 8, 6, 3, 1), 0)
  expect_warning(out <- .hull_calc(s, J = 4, gof_t = "CFI"),
                 class = "efa_hull_na_fit")
  expect_false(is.na(out$retain))
  # the dropped (NA-fit) solution carries no hull membership
  expect_true(is.na(out$s_complete[3, 4]))
})


burt <- .burt_cormat()

test_that("errors etc are thrown correctly", {
  skip_if_not_slow()
  expect_error(efa_hull(1:5), class = "efa_input_not_matrix")
  expect_message(suppressWarnings(efa_hull(GRiPS_raw)), class = "efa_cor_from_data")
  expect_warning(suppressMessages(efa_hull(GRiPS_raw)), class = "efa_hull_min_factors")
  expect_warning(
    suppressWarnings(
      suppressMessages(efa_hull(IDS2_R, N = 20)),
      classes = "efa_hull_min_factors"
    ),
    class = "efa_hull_few_solutions"
  )
  expect_warning(
    suppressWarnings(
      suppressMessages(efa_hull(GRiPS_raw, N = 20)),
      classes = "efa_hull_min_factors"
    ),
    class = "efa_n_from_data"
  )
  expect_error(efa_hull(test_models$baseline$cormat), class = "efa_n_required")
  expect_error(efa_hull(test_models$baseline$cormat, estimator = "ML"), class = "efa_n_required")
  expect_error(efa_hull(test_models$baseline$cormat, estimator = "ULS"), class = "efa_n_required")

  expect_error(efa_hull(sing_raw, estimator = "ML"), class = "efa_cor_singular")
  expect_error(efa_hull(sing_cor, N = sing_N, estimator = "ML"),
               class = "efa_cor_singular")

  expect_error(efa_hull(matrix(rnorm(50), ncol = 5)), class = "efa_hull_min_indicators")

  expect_message(suppressWarnings(efa_hull(GRiPS_raw)), class = "efa_hull_gof_caf")

  expect_warning(efa_hull(test_models$baseline$cormat, n_fac_theor = 13, N = 500), class = "efa_hull_max_factors")
})

test_that("a non-positive-definite input is smoothed and its hull flagged", {
  # burt has to be smoothed, its hull collapses below three solutions, and the
  # selected solution is inadmissible. Several warnings fire, so collect their
  # classes rather than stopping at the first one (or muffling them wholesale).
  # An explicit n_fac_theor is not varied here: the internal parallel analysis
  # already suggests more factors than any small theoretical value, so the search
  # bound J -- and with it the whole result -- is the same either way.
  #
  # efa_hull() derives J from a simulated parallel analysis, so seed it rather than
  # inheriting whatever state the preceding tests left (which differs between the
  # default and the slow gate).
  set.seed(42)
  classes <- .warn_classes(efa_hull(burt, N = 20, estimator = "ML"))

  expect_true("efa_cor_smoothed" %in% classes)
  expect_true("efa_hull_few_solutions" %in% classes)
  expect_true("efa_hull_inadmissible" %in% classes)
})

test_that("a hull with no finite goodness-of-fit value aborts", {
  # Every solution has an undefined fit value, so nothing can lie on the hull; the
  # elimination must abort with a classed condition instead of indexing an empty
  # matrix. The excluded-solutions warning fires first.
  s <- cbind(0:3, rep(NA_real_, 4), c(10, 6, 3, 1), 0)
  expect_error(
    suppressWarnings(.hull_calc(s, J = 3, gof_t = "CFI")),
    class = "efa_hull_no_fit"
  )
})

if (is_slow_test()) {
  rm(hull_cor_paf, hull_cor_ml, hull_cor_uls, hull_cor_uls_CFI, hull_raw_paf,
     hull_raw_ml, hull_raw_uls, hull_raw_uls_CFI, hull_raw_ml_nf, hull_cor_ml_nf,
     hull_PCA, hull_EFA)
}
rm(burt)
