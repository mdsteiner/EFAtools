# HULL refits the model at every candidate factor count for every (method, gof)
# combination, which makes the file-top fixture block the single heaviest in the
# suite (~60s). Skip the fixtures + the tests below by default; opt in with
# `Sys.setenv(EFATOOLS_TEST_SLOW = "true")` to run the full HULL coverage. See
# helper-slow.R for the convention.
if (is_slow_test()) {
hull_cor_paf <- suppressMessages(HULL(test_models$baseline$cormat, N = 500))
hull_cor_ml <- HULL(test_models$baseline$cormat, N = 500, method = "ML")
hull_cor_uls <- HULL(test_models$baseline$cormat, N = 500, method = "ULS")
hull_cor_uls_CFI <- HULL(test_models$baseline$cormat, N = 500, method = "ULS",
                         gof = "CFI")
hull_cor_ml_nf <- suppressWarnings(HULL(test_models$baseline$cormat, N = 500,
                                        method = "ML", n_fac_theor = 12))

hull_PCA <- HULL(test_models$baseline$cormat, N = 500, method = "ML",
                 eigen_type = "PCA")
hull_EFA <- HULL(test_models$baseline$cormat, N = 500, method = "ML",
                 eigen_type = "EFA")

hull_raw_paf <- suppressMessages(suppressWarnings(HULL(GRiPS_raw)))
hull_raw_ml <- suppressMessages(suppressWarnings(HULL(GRiPS_raw, method = "ML")))
hull_raw_uls <- suppressMessages(suppressWarnings(HULL(GRiPS_raw, method = "ULS")))
hull_raw_uls_CFI <- suppressMessages(suppressWarnings(HULL(GRiPS_raw,
                                                           method = "ULS",
                                                           gof = "CFI")))
hull_raw_ml_nf <- suppressMessages(suppressWarnings(HULL(GRiPS_raw, N = 500,
                                                         method = "ML",
                                                         n_fac_theor = 7)))
}  # is_slow_test()


test_that("output class and dimensions are correct", {
  skip_if_not_slow()
  objs <- list(hull_cor_paf, hull_cor_ml, hull_cor_uls, hull_cor_uls_CFI,
               hull_PCA, hull_EFA, hull_raw_paf, hull_raw_ml, hull_raw_uls,
               hull_raw_uls_CFI)
  for (obj in objs) {
    expect_s3_class(obj, "efa_retention")
    expect_length(obj, 7)
    expect_length(obj$settings, 8)
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

  expect_equal(hull_cor_uls_CFI$n_factors[["CFI"]], 1)
  expect_false("CAF" %in% names(hull_cor_uls_CFI$n_factors))
  expect_false("RMSEA" %in% names(hull_cor_uls_CFI$n_factors))

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
  skip_if_not_slow()
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
  skip_if_not_slow()
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


# Create singular correlation matrix for tests
set.seed(7)
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z, rnorm(10), rnorm(10), rnorm(10)), ncol = 6)
cor_sing <- stats::cor(dat_sing)

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

test_that("errors etc are thrown correctly", {
  skip_if_not_slow()
  expect_error(HULL(1:5), class = "efa_input_not_matrix")
  expect_message(suppressWarnings(HULL(GRiPS_raw)), class = "efa_cor_from_data")
  expect_warning(suppressMessages(HULL(GRiPS_raw)), class = "efa_hull_min_factors")
  expect_warning(suppressMessages(HULL(IDS2_R, N = 20)), class = "efa_hull_few_solutions")
  expect_warning(suppressMessages(HULL(GRiPS_raw, N = 20)), class = "efa_n_from_data")
  expect_error(HULL(test_models$baseline$cormat), class = "efa_n_required")
  expect_error(HULL(test_models$baseline$cormat, method = "ML"), class = "efa_n_required")
  expect_error(HULL(test_models$baseline$cormat, method = "ULS"), class = "efa_n_required")

  expect_error(HULL(dat_sing, method = "ML"), class = "efa_cor_singular")
  expect_error(HULL(cor_sing, N = 20, method = "ML"), class = "efa_cor_singular")

  expect_error(HULL(matrix(rnorm(50), ncol = 5)), class = "efa_hull_min_indicators")

  expect_message(suppressWarnings(HULL(GRiPS_raw)), 'Only CAF can be used as gof if method "PAF" is used. Setting gof to "CAF"\n')

  expect_warning(HULL(test_models$baseline$cormat, n_fac_theor = 13, N = 500), class = "efa_hull_max_factors")
  # the burt matrix is smoothed and selects an inadmissible/few-solution hull;
  # assert the smoothing warning and muffle the other (incidental) warnings
  suppressWarnings(
    expect_warning(HULL(burt, N = 20, method = "ML"), class = "efa_cor_smoothed"))
  suppressWarnings(
    expect_warning(HULL(burt, N = 20, method = "ML", n_fac_theor = 1), class = "efa_cor_smoothed"))
})

if (is_slow_test()) {
  rm(hull_cor_paf, hull_cor_ml, hull_cor_uls, hull_cor_uls_CFI, hull_raw_paf,
     hull_raw_ml, hull_raw_uls, hull_raw_uls_CFI, hull_raw_ml_nf, hull_cor_ml_nf,
     hull_PCA, hull_EFA)
}
rm(x, y, z, dat_sing, cor_sing, burt)
