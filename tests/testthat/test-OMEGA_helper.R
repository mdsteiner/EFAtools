## Tests for .OMEGA_LAVAAN --------

if (requireNamespace("lavaan", quietly = TRUE)) {
lav_mod_1 <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
              F2 =~ V7 + V8 + V9 + V10 + V11 + V12
              F3 =~ V13 + V14 + V15 + V16 + V17 + V18
              g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
                   V13 + V14 + V15 + V16 + V17 + V18'
lav_fit_1 <- suppressWarnings(lavaan::cfa(lav_mod_1,
                                          sample.cov = test_models$baseline$cormat,
                                          sample.nobs = 500, estimator = "ml",
                                          orthogonal = TRUE))
om_lav_bi_add <- .OMEGA_LAVAAN(lav_fit_1, g_name = "g")
om_lav_bi_noadd <- .OMEGA_LAVAAN(lav_fit_1, g_name = "g", add_ind = FALSE)

lav_mod_2 <- 'g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
                    V13 + V14 + V15 + V16 + V17 + V18'
lav_fit_2 <- suppressWarnings(lavaan::cfa(lav_mod_2,
                                         sample.cov = test_models$baseline$cormat,
                                         sample.nobs = 500, estimator = "ml",
                                         orthogonal = TRUE))
om_lav_1_add <- suppressMessages(.OMEGA_LAVAAN(lav_fit_2))
om_lav_1_noadd <- suppressMessages(.OMEGA_LAVAAN(lav_fit_2, add_ind = FALSE))

lav_fit_3 <- suppressWarnings(lavaan::cfa(lav_mod_1, sample.cov =
                                            list(test_models$baseline$cormat,
                                                 test_models$baseline$cormat),
                                          sample.nobs = c(500, 500),
                                          estimator = "ml", orthogonal = TRUE))
om_lav_gr_add <- .OMEGA_LAVAAN(lav_fit_3, g_name = "g", group_names = c("Some",
                                                                    "Others"))
om_lav_gr_noadd <- .OMEGA_LAVAAN(lav_fit_3, g_name = "g", add_ind = FALSE,
                                 group_names = c("Some", "Others"))

lav_fit_4 <- suppressWarnings(lavaan::cfa(lav_mod_2, sample.cov =
                                            list(test_models$baseline$cormat,
                                                 test_models$baseline$cormat),
                                          sample.nobs = c(500, 500),
                                          estimator = "ml", orthogonal = TRUE))
om_lav_1_gr_add <- suppressMessages(.OMEGA_LAVAAN(lav_fit_4,
                                                  group_names = c("Some",
                                                                  "Others")))
om_lav_1_gr_noadd <- suppressMessages(.OMEGA_LAVAAN(lav_fit_4, add_ind = FALSE,
                                                    group_names = c("Some",
                                                                    "Others")))

lav_mod_ho_1 <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
               F2 =~ V7 + V8 + V9 + V10 + V11 + V12
               F3 =~ V13 + V14 + V15 + V16 + V17 + V18
               g =~ F1 + F2 + F3'
lav_fit_ho_1 <- suppressWarnings(lavaan::cfa(lav_mod_ho_1,
                                           sample.cov = test_models$baseline$cormat,
                                           sample.nobs = 500, estimator = "ml"))
om_lav_ho_1 <- suppressMessages(.OMEGA_LAVAAN(lav_fit_ho_1, g_name = "g"))
}

test_that("output class and dimensions are correct", {
  skip_if_not_installed("lavaan")
  expect_s3_class(om_lav_bi_add, "OMEGA")
  expect_s3_class(om_lav_1_add, "OMEGA")
  expect_s3_class(om_lav_ho_1, "OMEGA")
  expect_s3_class(om_lav_gr_add, "OMEGA")
  expect_s3_class(om_lav_1_gr_add, "OMEGA")
  expect_s3_class(om_lav_bi_noadd, "OMEGA")
  expect_s3_class(om_lav_1_noadd, "OMEGA")
  expect_s3_class(om_lav_gr_noadd, "OMEGA")
  expect_s3_class(om_lav_1_gr_noadd, "OMEGA")

  expect_output(str(om_lav_bi_add), "List of 2")
  expect_output(str(om_lav_1_add), "OMEGA")
  expect_output(str(om_lav_ho_1), "List of 2")
  expect_output(str(om_lav_gr_add), "List of 2")
  expect_output(str(om_lav_1_gr_add), "List of 2")
  expect_output(str(om_lav_bi_add), "List of 2")
  expect_output(str(om_lav_1_add), "OMEGA")
  expect_output(str(om_lav_ho_1), "List of 2")
  expect_output(str(om_lav_gr_add), "List of 2")
  expect_output(str(om_lav_1_gr_add), "List of 2")
})

test_that("output is correct (including group names for multiple groups)", {
  skip_if_not_installed("lavaan")
  expect_equal(rowSums(om_lav_bi_add[, 2:3]), om_lav_bi_add[, 1], tolerance = 1e-3)
  expect_equal(om_lav_bi_add[1, 4], 0.849, tolerance = 1e-3)
  expect_equal(om_lav_bi_add[1, 5], 0.672, tolerance = 1e-3)
  expect_equal(om_lav_bi_add[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(unname(om_lav_1_add[1]), 0.868, tolerance = 1e-3)
  expect_equal(unname(om_lav_1_add[2]), 0.869, tolerance = 1e-3)
  expect_equal(rowSums(om_lav_ho_1[, 2:3]), om_lav_ho_1[, 1], tolerance = 1e-3)
  expect_equal(om_lav_ho_1[1, 4], 0.848, tolerance = 1e-3)
  expect_equal(om_lav_ho_1[1, 5], 0.684, tolerance = 1e-3)
  expect_equal(om_lav_ho_1[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(rowSums(om_lav_gr_add$Some[, 2:3]), om_lav_gr_add$Some[, 1],
               tolerance = 1e-3)
  expect_equal(rowSums(om_lav_gr_add$Others[, 2:3]), om_lav_gr_add$Others[, 1],
               tolerance = 1e-3)
})

# Preparations for error tests
if (requireNamespace("lavaan", quietly = TRUE)) {
lav_mod_NA <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6 + V17
               F2 =~ V7 + V8 + V9 + V10 + V11 + V12 + V2
               F3 =~ V13 + V14 + V15 + V16 + V17 + V18 + V10
               g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
                    V13 + V14 + V15 + V16 + V17 + V18'
lav_fit_NA <- suppressWarnings(lavaan::cfa(lav_mod_NA,
                                           sample.cov = test_models$baseline$cormat,
                                           sample.nobs = 500, estimator = "ml"))

lav_mod_inv <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
               F2 =~ V7 + V8 + V9 + V10 + V11 + V12
               F3 =~ V13 + V14 + V15 + V16 + V17 + V18'
lav_fit_inv <- suppressWarnings(lavaan::cfa(lav_mod_inv,
                                           sample.cov = test_models$baseline$cormat,
                                           sample.nobs = 500, estimator = "ml"))

lav_mod_bi_red <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
               F2 =~ V7 + V8 + V9 + V10 + V11 + V12
               F3 =~ V13 + V14 + V15 + V16 + V17 + V18
               g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 +
                    V13 + V14 + V15 + V16 + V17 + V18'
lav_fit_bi_red <- suppressWarnings(lavaan::cfa(lav_mod_bi_red,
                                              sample.cov = test_models$baseline$cormat,
                                              sample.nobs = 500, estimator = "ml",
                                              orthogonal = TRUE))

lav_mod_ho_2 <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
               F2 =~ V7 + V8 + V9
               F3 =~ V10 + V11 + V12
               F4 =~ V13 + V14 + V15 + V16 + V17 + V18
               g =~ F1 + F2
               h =~ F3 + F4'
lav_fit_ho_2 <- suppressWarnings(lavaan::cfa(lav_mod_ho_2,
                                           sample.cov = test_models$baseline$cormat,
                                           sample.nobs = 500, estimator = "ml"))
}

test_that("errors are thrown correctly", {
  skip_if_not_installed("lavaan")
  # OMEGA's lavaan path now delegates spec-building (and its validation) to the
  # shared reliability adapter, so the aborts carry the adapter's
  # efa_reliability_* classes; the informational notices remain OMEGA's own
  # efa_omega_* conditions.
  expect_error(.OMEGA_LAVAAN(lav_fit_NA, g_name = "g"), class = "efa_reliability_na_loadings")
  expect_error(.OMEGA_LAVAAN(lav_fit_1, g_name = "fu"), class = "efa_reliability_g_name")
  expect_message(.OMEGA_LAVAAN(lav_fit_2, add_ind = FALSE), class = "efa_omega_single_factor")
  expect_message(.OMEGA_LAVAAN(lav_fit_2), class = "efa_omega_single_factor")
  expect_message(.OMEGA_LAVAAN(lav_fit_ho_1, g_name = "g"), class = "efa_omega_g_second_order")
  expect_error(.OMEGA_LAVAAN(lav_fit_inv, g_name = "F3"), class = "efa_reliability_invalid_lavaan")
  expect_message(.OMEGA_LAVAAN(lav_fit_bi_red, g_name = "g"), class = "efa_omega_few_loadings")
  expect_error(.OMEGA_LAVAAN(lav_fit_ho_2, g_name = "g"), class = "efa_reliability_higher_order")
})

test_that("a bifactor model fitted with correlated factors is rejected", {
  skip_if_not_installed("lavaan")
  # The coefficients split a composite's variance into a general and a per-group-factor
  # part, which needs uncorrelated factors. lavaan::cfa() does not impose that, and the
  # same syntax fitted without orthogonal = TRUE returns factor correlations the
  # coefficients would silently drop -- here up to .24 on the covariance scale, worth
  # .016 on the whole-scale omega total.
  lav_fit_obl <- suppressWarnings(lavaan::cfa(lav_mod_1,
                                              sample.cov = test_models$baseline$cormat,
                                              sample.nobs = 500, estimator = "ml",
                                              orthogonal = FALSE))
  expect_error(.rel_adapt_lavaan(lav_fit_obl, g_name = "g"),
               class = "efa_reliability_correlated_factors")
  expect_error(.OMEGA_LAVAAN(lav_fit_obl, g_name = "g"),
               class = "efa_reliability_correlated_factors")
  expect_error(efa_reliability(lav_fit_obl, g_name = "g"),
               class = "efa_reliability_correlated_factors")
  # The supported shapes are unaffected: an orthogonal bifactor, a second-order model
  # and a single factor all leave psi diagonal.
  expect_no_error(.rel_adapt_lavaan(lav_fit_1, g_name = "g"))
  expect_no_error(suppressMessages(.rel_adapt_lavaan(lav_fit_ho_1, g_name = "g")))
  expect_no_error(.rel_adapt_lavaan(lav_fit_2))
})


## Tests for .OMEGA_FLEX -------

## Use with an output from the SL function, with type EFAtools
efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
sl_mod <- SL(efa_mod, type = "EFAtools", method = "PAF")
om_sl_add <- .OMEGA_FLEX(sl_mod, type = "EFAtools",
                     factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                     variance = "correlation")
om_sl_noadd <- .OMEGA_FLEX(sl_mod, type = "EFAtools",
                     factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                     variance = "correlation", add_ind = FALSE)

# Explicit factor names
om_sl_named_add <- .OMEGA_FLEX(sl_mod, type = "EFAtools",
                           fac_names = c("Fa1", "Fa2", "Fa3"),
                           factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                           variance = "correlation")
om_sl_named_noadd <- .OMEGA_FLEX(sl_mod, type = "EFAtools",
                           fac_names = c("Fa1", "Fa2", "Fa3"),
                           factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                           variance = "correlation", add_ind = FALSE)

## Use with an output from the psych::schmid function, with type psych
schmid_mod <- psych::schmid(test_models$baseline$cormat, nfactors = 3,
                            n.obs = 500, fm = "pa", rotate = "Promax")

# Omega also with type "psych"
# Find correlation matrix from phi and pattern matrix from psych::schmid output
om_schmid_1_add <- .OMEGA_FLEX(schmid_mod, type = "psych", variance = "correlation")
om_schmid_1_noadd <- .OMEGA_FLEX(schmid_mod, type = "psych", variance = "correlation",
                                add_ind = FALSE)
# Enter correlation matrix manually
om_schmid_2 <- .OMEGA_FLEX(schmid_mod, type = "psych",
                           cormat = test_models$baseline$cormat,
                           variance = "correlation")

## Manually specify components (here with type "EFAtools")
om_man_1_add <- .OMEGA_FLEX(model = NULL, type = "EFAtools",
                            var_names = rownames(sl_mod$sl),
                  g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                  u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
                  factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                  variance = "correlation")
om_man_1_noadd <- .OMEGA_FLEX(model = NULL, type = "EFAtools",
                              var_names = rownames(sl_mod$sl),
                              g_load = sl_mod$sl[, "g"],
                              s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                              u2 = sl_mod$sl[, "u2"],
                              cormat = test_models$baseline$cormat,
                              factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                              variance = "correlation", add_ind = FALSE)
# Now with other type of variance (model-based instead of based on corrmat)
om_man_2 <- .OMEGA_FLEX(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                  g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                  u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
                  factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                  variance = "sums_load")

test_that("output class and dimensions are correct", {
  expect_s3_class(om_sl_add, "OMEGA")
  expect_s3_class(om_sl_named_add, "OMEGA")
  expect_s3_class(om_schmid_1_add, "OMEGA")
  expect_s3_class(om_schmid_2, "OMEGA")
  expect_s3_class(om_man_1_add, "OMEGA")
  expect_s3_class(om_man_2, "OMEGA")
  expect_s3_class(om_sl_noadd, "OMEGA")
  expect_s3_class(om_sl_named_noadd, "OMEGA")
  expect_s3_class(om_schmid_1_noadd, "OMEGA")
  expect_s3_class(om_man_1_noadd, "OMEGA")

  expect_output(str(om_sl_add), "List of 2")
  expect_output(str(om_sl_named_add), "List of 2")
  expect_output(str(om_schmid_1_add), "List of 2")
  expect_output(str(om_schmid_2), "List of 2")
  expect_output(str(om_man_1_add), "List of 2")
  expect_output(str(om_man_2), "List of 2")
  expect_output(str(om_sl_noadd), "List of 2")
  expect_output(str(om_sl_named_noadd), "List of 2")
  expect_output(str(om_schmid_1_noadd), "List of 2")
  expect_output(str(om_man_1_noadd), "List of 2")
})

test_that("output is correct", {
  # On a subscale row, omega total counts every factor the composite's variables load
  # on while omega subscale counts only that factor's own contribution, so the three
  # coefficients no longer add up: the gap is what the composite receives from the
  # other group factors through its cross-loadings. Every fixture here has an estimated
  # (dense) group-loading matrix, so that gap is strictly positive -- asserting only
  # `tot >= hier + sub` would pass just as well if the cross terms were dropped again,
  # which makes the gap exactly zero.
  tol <- sqrt(.Machine$double.eps)
  for (om in list(om_sl_add, om_sl_named_add, om_schmid_1_add, om_schmid_2,
                  om_man_1_add, om_man_2)) {
    expect_true(all(om[2:4, 1] - rowSums(om[2:4, 2:3]) > tol))
  }
  # And the gap is exactly the mass the composite draws from the other group factors,
  # 1' L_other L_other' 1 over its observed variance.
  s_all <- unclass(sl_mod$sl[, c("F1", "F2", "F3")])
  map_sl <- s_all >= .2
  for (j in 1:3) {
    mem <- which(map_sl[, j])
    expect_equal(unname(om_sl_add[j + 1, 1] - sum(om_sl_add[j + 1, 2:3])),
                 sum(colSums(s_all[mem, -j, drop = FALSE])^2) /
                   sum(test_models$baseline$cormat[mem, mem]))
  }
  expect_equal(om_sl_add[1, 4], 0.842, tolerance = 1e-3)
  expect_equal(om_sl_add[1, 5], 0.652, tolerance = 1e-3)
  expect_equal(om_sl_add[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(om_sl_named_add[1, 4], 0.842, tolerance = 1e-3)
  expect_equal(om_sl_named_add[1, 5], 0.652, tolerance = 1e-3)
  expect_equal(om_sl_named_add[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(om_schmid_1_add[1, 4], 0.845, tolerance = 1e-3)
  expect_equal(om_schmid_1_add[1, 5], 0.661, tolerance = 1e-3)
  expect_equal(om_schmid_1_add[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(om_schmid_2[1, 4], 0.845, tolerance = 1e-3)
  expect_equal(om_schmid_2[1, 5], 0.661, tolerance = 1e-3)
  expect_equal(om_schmid_2[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(om_man_1_add[1, 4], 0.842, tolerance = 1e-3)
  expect_equal(om_man_1_add[1, 5], 0.652, tolerance = 1e-3)
  expect_equal(om_man_1_add[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(om_man_2[1, 4], 0.842, tolerance = 1e-3)
  expect_equal(om_man_2[1, 5], 0.652, tolerance = 1e-3)
  expect_equal(om_man_2[1, 6], 0.706, tolerance = 1e-3)

  expect_gte(om_sl_add[1, 1] + tol, sum(om_sl_add[1, 2:3]))
  expect_gte(om_sl_named_add[1, 1] + tol, sum(om_sl_named_add[1, 2:3]))
  expect_gte(om_schmid_1_add[1, 1] + tol, sum(om_schmid_1_add[1, 2:3]))
  expect_gte(om_schmid_2[1, 1] + tol, sum(om_schmid_2[1, 2:3]))
  expect_gte(om_man_1_add[1, 1] + tol, sum(om_man_1_add[1, 2:3]))
  expect_gte(om_man_2[1, 1] + tol, sum(om_man_2[1, 2:3]))

  # With variance = "sums_load", the general-factor omega total equals
  # hierarchical + subscale exactly (shared total-variance denominator).
  expect_equal(unname(om_man_2[1, 1]), unname(sum(om_man_2[1, 2:3])))
})

test_that("variance = 'sums_load' whole-scale omega total uses the full loading columns", {
  # A dense group-loading matrix: item 1 is assigned to F1 but also loads on F2.
  # The whole-scale composite contains all items, so its model-implied omega total
  # counts every group loading (full columns), giving 1 - sum(u2)/V, with the
  # general and group variances partitioning it exactly.
  g_load <- rep(0.5, 6)
  s_load <- matrix(0, 6, 2); s_load[1:3, 1] <- 0.5; s_load[4:6, 2] <- 0.5
  s_load[1, 2] <- 0.3
  u2 <- 1 - g_load^2 - rowSums(s_load^2)
  fc <- matrix(0, 6, 2); fc[1:3, 1] <- 1; fc[4:6, 2] <- 1   # item 1 assigned to F1 only

  om <- .OMEGA_FLEX(model = NULL, type = "EFAtools", var_names = paste0("V", 1:6),
                    g_load = g_load, s_load = s_load, u2 = u2,
                    factor_corres = fc, variance = "sums_load")

  full_s <- colSums(s_load)                       # full columns (all items)
  V <- sum(g_load)^2 + sum(full_s^2) + sum(u2)    # model-implied composite variance

  # omega total is the model-implied total, using the full loading columns, and the
  # whole-scale subscale omega is the full group-factor variance (not the
  # assigned-only mass), so tot = hier + sub exactly.
  expect_equal(unname(om[1, "tot"]), 1 - sum(u2) / V)
  expect_equal(unname(om[1, "sub"]), sum(full_s^2) / V)
  expect_equal(unname(om[1, "tot"]), unname(om[1, "hier"] + om[1, "sub"]))
  # the off-assignment loading makes the full mass strictly exceed the assigned one
  expect_gt(sum(full_s^2), sum(colSums(s_load * fc)^2))
})

test_that("type = 'psych' reproduces psych::omega (g omegas and ECV)", {
  skip_on_cran()
  skip_if_not_installed("psych")

  po <- psych::omega(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                     fm = "pa", rotate = "Promax", plot = FALSE)

  # The psych-type path derives the same variable-to-factor correspondences as
  # psych::omega and bases the ECV on the (unzeroed) group loadings, so the ECV and
  # omega hierarchical are the same quantity and match exactly. Omega total is not: it
  # is the model-implied common variance of all variables here and the observed total
  # variance less the unique variances in psych, which agree only to the extent that the
  # model reproduces the correlations. This fixture is a three-factor population matrix
  # fitted with three factors, so the residual mass is about 6e-5 and the tolerance below
  # holds with room to spare -- a worse-fitting fixture would separate the two on purpose.
  expect_equal(unname(om_schmid_1_add[1, "ECV"]), unname(po$ECV), tolerance = 1e-4)
  expect_equal(unname(om_schmid_1_add[1, "tot"]), po$omega.tot, tolerance = 1e-3)
  expect_equal(unname(om_schmid_1_add[1, "hier"]), po$omega_h, tolerance = 1e-3)
})

test_that("errors are thrown correctly", {
  expect_error(.OMEGA_FLEX(schmid_mod, type = "psych",
                           cormat = matrix(rnorm(50), ncol = 5),
                           variance = "correlation"), class = "efa_omega_not_cormat")
  expect_error(.OMEGA_FLEX(model = NULL, type = "EFAtools",
                           var_names = rownames(sl_mod$sl),
                           g_load = sl_mod$sl[, "g"],
                           s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                           u2 = sl_mod$sl[, "u2"],
                           factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                           variance = "correlation"), class = "efa_omega_need_cormat")
  expect_error(.OMEGA_FLEX(model = NULL, type = "EFAtools",
                           var_names = rownames(sl_mod$sl),
                           g_load = sl_mod$sl[, "g"],
                           s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                           u2 = sl_mod$sl[, "u2"], cormat = matrix(rnorm(50), ncol = 5),
                           factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                           variance = "correlation"), class = "efa_omega_not_cormat")
  expect_error(.OMEGA_FLEX(schmid_mod, type = "EFAtools",
                           variance = "correlation"), class = "efa_omega_need_corres")
  expect_warning(.OMEGA_FLEX(schmid_mod, type = "psych",
                             variance = "sums_load"), class = "efa_omega_variance_override")
  expect_warning(.OMEGA_FLEX(schmid_mod, type = "psych",
                             factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                             variance = "correlation"), class = "efa_omega_corres_override")
  # A factor_corres with the wrong number of rows (one fewer than the number of
  # items) must be rejected, not silently used to corrupt the omega outputs.
  expect_error(.OMEGA_FLEX(sl_mod, type = "EFAtools",
                           factor_corres = (sl_mod$sl[, c("F1", "F2", "F3")] >= .2)[-1, ],
                           variance = "correlation"))
})

test_that("group-factor H index uses only the factor's assigned items", {
  # Item 6 has a sizeable loading on F2 but is assigned to F1; that cross-loading
  # must not enter F2's H index (Hancock & Mueller, 2001).
  g_load <- rep(0.5, 6)
  s_load <- matrix(0, 6, 2)
  s_load[1:3, 1] <- 0.6
  s_load[4:5, 2] <- 0.6
  s_load[6, 1] <- 0.5
  s_load[6, 2] <- 0.5
  fc <- matrix(0, 6, 2)
  fc[1:3, 1] <- 1
  fc[6, 1] <- 1
  fc[4:5, 2] <- 1
  u2 <- 1 - g_load^2 - rowSums(s_load^2)
  om <- .OMEGA_FLEX(model = NULL, type = "EFAtools", var_names = paste0("V", 1:6),
                    g_load = g_load, s_load = s_load, u2 = u2,
                    factor_corres = fc, variance = "sums_load")
  s2 <- s_load[fc[, 2] == 1, 2]
  expect_equal(unname(om["2", "H"]), 1 / (1 + 1 / sum(s2^2 / (1 - s2^2))))
})

test_that("PUC counts each contaminated pair once (no double counting)", {
  # Items 3 and 4 belong to BOTH group factors, so the pair {3, 4} is a single
  # contaminated correlation even though it sits within two group factors.
  fc <- matrix(c(1, 1, 1, 1, 0, 0,
                 0, 0, 1, 1, 1, 1), ncol = 2)
  s_load <- 0.4 * fc
  g_load <- rep(0.5, 6)
  u2 <- 1 - g_load^2 - rowSums(s_load^2)
  om <- .OMEGA_FLEX(model = NULL, type = "EFAtools", var_names = paste0("V", 1:6),
                    g_load = g_load, s_load = s_load, u2 = u2,
                    factor_corres = fc, variance = "sums_load")
  # C(6, 2) = 15 pairs; 11 unique contaminated -> PUC = 1 - 11/15.
  expect_equal(unname(om[1, "PUC"]), 1 - 11 / 15)
})

test_that("PUC for an incomplete-g bifactor depends only on the group factors", {
  skip_if_not_installed("lavaan")
  # Bifactor where the general factor does not load on V18; PUC must still be the
  # group-factor-only quantity (Reise et al., 2013): 1 - 3 * C(6,2) / C(18,2).
  mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
          F2 =~ V7 + V8 + V9 + V10 + V11 + V12
          F3 =~ V13 + V14 + V15 + V16 + V17 + V18
          g  =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
                V13 + V14 + V15 + V16 + V17'
  fit <- suppressWarnings(lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
                                      sample.nobs = 500, estimator = "ml",
                                      orthogonal = TRUE))
  om <- suppressMessages(suppressWarnings(.OMEGA_LAVAAN(fit, g_name = "g")))
  expect_equal(unname(om[1, "PUC"]), 1 - 3 * choose(6, 2) / choose(18, 2))
})

test_that("an inconsistent cormat triggers an out-of-range warning", {
  # Strong general loadings combined with a weakly correlated matrix push
  # omega hierarchical above 1 (and above omega total); this must not pass silently.
  g_load <- rep(0.7, 6)
  s_load <- matrix(0, 6, 2)
  s_load[1:3, 1] <- 0.1
  s_load[4:6, 2] <- 0.1
  u2 <- 1 - g_load^2 - rowSums(s_load^2)
  Rm <- matrix(0.3, 6, 6)
  diag(Rm) <- 1
  rownames(Rm) <- colnames(Rm) <- paste0("V", 1:6)
  expect_warning(
    .OMEGA_FLEX(model = NULL, type = "EFAtools", var_names = paste0("V", 1:6),
                g_load = g_load, s_load = s_load, u2 = u2, cormat = Rm,
                factor_corres = s_load > 0, variance = "correlation"),
    class = "efa_omega_out_of_range")

  # And the coefficients stay in range on solutions whose loadings came from the
  # correlation matrix they are scored against. Every omega total is now the composite's
  # model-implied common variance, which is weakly larger than the general-plus-own-column
  # term it replaced, so this is the direction in which the guard could start firing on
  # ordinary input. Asserted on the values rather than on the absence of the warning,
  # because an over-extracted solution legitimately raises the empty-factor one.
  tol_range <- sqrt(.Machine$double.eps)
  for (k in 2:5) {
    mod_k <- suppressWarnings(EFA(test_models$baseline$cormat, N = 500, n_factors = k,
                                  type = "EFAtools", method = "PAF",
                                  rotation = "promax"))
    sl_k <- suppressWarnings(SL(mod_k, type = "psych", method = "PAF"))
    om_k <- suppressWarnings(
      .OMEGA_FLEX(sl_k, type = "psych", cormat = test_models$baseline$cormat,
                  variance = "correlation"))
    expect_true(all(om_k[, c("tot", "hier", "sub")] <= 1 + tol_range, na.rm = TRUE))
    expect_true(all(om_k[, "hier"] <= om_k[, "tot"] + tol_range, na.rm = TRUE))
  }
})

test_that("the out-of-range guard runs in both variance conventions", {
  # A share of variance above 1 is inadmissible whichever total variance it was divided by.
  # One way there in the model-implied convention is an improper solution, whose unique
  # variances sum to something negative for a composite; the Heywood check reports that from
  # the uniquenesses, but only when the indices are computed, and never says that the
  # returned coefficients are themselves out of range.
  g <- rep(0.6, 6)
  s <- matrix(0, 6, 2)
  s[1:3, 1] <- 0.5
  s[4:6, 2] <- 0.5
  spec <- list(g_load = g, s_load = s, u2 = rep(-0.2, 6),
               map = cbind(rep(1:0, each = 3), rep(0:1, each = 3)), cormat = NULL,
               var_names = paste0("V", 1:6), fac_names = c("F1", "F2"))
  # Each of these solutions raises both conditions, so the one an expectation is not about is
  # muffled to keep it from leaking out of the test.
  muffle <- function(expr, cls) {
    withCallingHandlers(expr,
      warning = function(w) if (inherits(w, cls)) invokeRestart("muffleWarning"))
  }
  expect_warning(om <- muffle(.reliability_core(spec, "sums_load"),
                              "efa_reliability_heywood"),
                 class = "efa_omega_out_of_range")
  expect_gt(unname(om["g", "tot"]), 1)

  # A proper solution with uncorrelated group factors stays silent, as it always has been.
  spec$u2 <- 1 - g^2 - rowSums(s^2)
  expect_no_warning(.reliability_core(spec, "sums_load"),
                    class = "efa_omega_out_of_range")

  # The second route needs no improper solution at all, which is what makes the guard worth
  # running in this convention: omega subscale divides a factor's own squared loading sum by
  # its composite's whole variance, and once the group factors are correlated that sum is not
  # a part of it -- the composite's cross-loadings on the other factor contribute negative
  # common variance and take its total below what its own factor puts in. Every uniqueness
  # here is positive and Phi is positive definite, so nothing else in the function signals it.
  S <- matrix(0, 8, 2)
  S[1:4, 1] <- .85; S[1:4, 2] <- -.30
  S[5:8, 2] <- .85; S[5:8, 1] <- -.30
  Phi <- matrix(c(1, .6, .6, 1), 2, 2)
  u2 <- 1 - diag(S %*% Phi %*% t(S))
  proper <- list(g_load = rep(0, 8), s_load = S, u2 = u2, Phi = Phi,
                 map = cbind(rep(1:0, each = 4), rep(0:1, each = 4)), cormat = NULL,
                 var_names = paste0("V", 1:8), fac_names = c("F1", "F2"))
  expect_true(all(u2 > 0))
  expect_gt(min(eigen(Phi, symmetric = TRUE, only.values = TRUE)$values), 0)

  expect_warning(om <- .reliability_core(proper, "sums_load"),
                 class = "efa_omega_out_of_range")
  expect_gt(unname(om["F1", "sub"]), 1)
  # and the improper-solution report stays silent, so this is the only signal there is
  expect_no_warning(muffle(.reliability_core(proper, "sums_load"),
                           "efa_omega_out_of_range"),
                    class = "efa_reliability_heywood")
})

test_that("a single group-factor column is handled without error", {
  g_load <- rep(0.7, 6)
  s_load <- matrix(0.2, 6, 1)
  u2 <- 1 - g_load^2 - rowSums(s_load^2)
  fc <- matrix(1, 6, 1)
  expect_s3_class(
    .OMEGA_FLEX(model = NULL, type = "EFAtools", var_names = paste0("V", 1:6),
                g_load = g_load, s_load = s_load, u2 = u2,
                factor_corres = fc, variance = "sums_load"), "OMEGA")
  Rm <- matrix(0.6, 6, 6)
  diag(Rm) <- 1
  rownames(Rm) <- colnames(Rm) <- paste0("V", 1:6)
  expect_s3_class(
    suppressWarnings(.OMEGA_FLEX(model = NULL, type = "EFAtools",
                var_names = paste0("V", 1:6), g_load = g_load, s_load = s_load,
                u2 = u2, cormat = Rm, factor_corres = fc,
                variance = "correlation")), "OMEGA")
})

test_that("a group factor with no assigned items reports NA coefficients and warns", {
  # type = 'psych' assigns each item to its highest-loading factor; F2 never wins,
  # so its correspondence column is all zero and its coefficients are undefined.
  g <- rep(0.5, 6)
  s <- matrix(c(rep(0.6, 6), rep(0.3, 6)), ncol = 2)
  u2 <- 1 - g^2 - rowSums(s^2)
  # Model-implied (consistent) correlation matrix so only the empty-factor warning
  # fires (no out-of-range, no variance override).
  L <- cbind(g, s)
  Rm <- L %*% t(L); diag(Rm) <- 1
  rownames(Rm) <- colnames(Rm) <- paste0("V", 1:6)
  expect_warning(
    om <- .OMEGA_FLEX(model = NULL, type = "psych", var_names = paste0("V", 1:6),
                      g_load = g, s_load = s, u2 = u2, cormat = Rm,
                      variance = "correlation"),
    class = "efa_omega_empty_factor")
  expect_true(all(is.na(om["2", c("tot", "hier", "sub", "H")])))
})

test_that("PUC membership follows factor_corres, not the loading magnitudes", {
  # Item 2 is assigned to F1 by factor_corres but has a zero loading on F1; it must
  # still count as an F1 member for PUC, consistent with the H index.
  g <- rep(0.5, 4)
  s <- matrix(0, 4, 2); s[1, 1] <- 0.4; s[3, 2] <- 0.4; s[4, 2] <- 0.4
  fc <- matrix(0, 4, 2); fc[1:2, 1] <- 1; fc[3:4, 2] <- 1
  u2 <- 1 - g^2 - rowSums(s^2)
  om <- .OMEGA_FLEX(model = NULL, type = "EFAtools", var_names = paste0("V", 1:4),
                    g_load = g, s_load = s, u2 = u2, factor_corres = fc,
                    variance = "sums_load")
  # F1 = {1, 2}, F2 = {3, 4}: 2 contaminated pairs of C(4, 2) = 6 -> PUC = 1 - 2/6.
  expect_equal(unname(om[1, "PUC"]), 1 - 2 / 6)
})

test_that("the sums_load whole-scale row honours the factor correlations", {
  # The whole-scale composite's model-implied common variance is 1' S Phi S' 1, the same
  # quantity the subscale rows of this branch read off `common`. Scoring the whole-scale
  # row from the squared loading-column sums alone drops the factor correlations, which
  # would put an orthogonal model and a correlated one in one table.
  S <- matrix(0, 6, 2)
  S[1:3, 1] <- c(.75, .70, .65)
  S[4:6, 2] <- c(.60, .65, .70)
  Phi <- matrix(c(1, .5, .5, 1), 2, 2)
  common <- S %*% Phi %*% t(S)
  u2 <- 1 - diag(common)
  spec <- list(g_load = rep(0, 6), s_load = S, u2 = u2, Phi = Phi,
               map = cbind(rep(1:0, each = 3), rep(0:1, each = 3)), cormat = NULL,
               var_names = paste0("V", 1:6), fac_names = c("F1", "F2"))

  om <- suppressWarnings(.reliability_core(spec, "sums_load"))
  expect_equal(unname(om["g", "tot"]), sum(common) / (sum(common) + sum(u2)))

  # Not a vacuous identity: the factor correlations move it well away from the value the
  # column sums alone give.
  s <- colSums(S)
  expect_gt(unname(om["g", "tot"]), sum(s^2) / (sum(s^2) + sum(u2)))

  # The same spec without factor correlations keeps the loading-column-sum arithmetic
  # exactly, which is every solution the shipped adapters can produce.
  spec$Phi <- NULL
  om_orth <- suppressWarnings(.reliability_core(spec, "sums_load"))
  expect_identical(unname(om_orth["g", "tot"]), sum(s^2) / (sum(s^2) + sum(u2)))
})

test_that("variables keyed against each other are flagged, in both variance modes", {
  cls <- "efa_reliability_mixed_keying"
  make <- function(g, s) {
    u2 <- 1 - g^2 - rowSums(s^2)
    L <- cbind(g, s)
    Rm <- L %*% t(L)
    diag(Rm) <- 1
    rownames(Rm) <- colnames(Rm) <- paste0("V", seq_along(g))
    list(g_load = g, s_load = s, u2 = u2, map = abs(s) > 0, cormat = Rm,
         var_names = paste0("V", seq_along(g)), fac_names = c("F1", "F2"))
  }

  # Half the variables keyed the other way, so the true score variance of their
  # unit-weighted sum cancels instead of accumulating.
  s <- matrix(0, 6, 2)
  s[1:3, 1] <- 0.2
  s[4:6, 2] <- -0.2
  mixed <- make(c(rep(0.7, 3), rep(-0.7, 3)), s)
  for (v in c("correlation", "sums_load")) {
    expect_warning(.reliability_core(mixed, v), class = cls)
  }

  # The same components all keyed alike are silent.
  aligned <- make(rep(0.7, 6), abs(s))
  for (v in c("correlation", "sums_load")) {
    expect_no_warning(.reliability_core(aligned, v), class = cls)
  }

  # Nothing to judge: no common variance at all, or a malformed component the callers
  # report separately.
  expect_no_warning(.rel_check_keying(matrix(0, 4, 4)), class = cls)
  expect_no_warning(.rel_check_keying(matrix(NA_real_, 4, 4)), class = cls)

  # The variables the message names are those subtracting from the composite, which is
  # guidance rather than the trigger: a solution can fall below the cutoff with no single
  # variable's own common variance turning negative, and it must still warn.
  opposed <- matrix(c(1, -0.4, -0.4, 1), 2, 2)
  expect_true(all(rowSums(opposed) > 0))
  expect_warning(.rel_check_keying(opposed), class = cls)

  # A single factor reaches the core on a spec with no group factors, whose common
  # variance comes from the general column alone; the check must fire there too.
  single <- function(g) .rel_single_factor_spec(
    list(g_load = rep(0, length(g)), s_load = matrix(g, ncol = 1), u2 = rep(0.51, 6),
         cormat = tcrossprod(g) + diag(0.51, 6),
         var_names = paste0("V", seq_along(g))))
  for (v in c("correlation", "sums_load")) {
    expect_warning(.reliability_core(single(c(rep(0.7, 3), rep(-0.7, 3))), v),
                   class = cls)
    expect_no_warning(.reliability_core(single(rep(0.7, 6)), v), class = cls)
  }
})

test_that(".rel_single_factor_label takes a name where there is one, F1 where there is not", {
  # A name is used as given ...
  expect_identical(.rel_single_factor_label("Ability"), "Ability")
  expect_identical(.rel_single_factor_label("g"), "g")
  # ... including one that happens to read as a digit.
  expect_identical(.rel_single_factor_label("1"), "1")

  # ... and everything that is not a name takes the default first-factor label: no name at
  # all, a bare column position (what .rel_adapt_manual() falls back to), a missing one, and
  # an empty one. Not "g", which would claim a general factor over group factors a
  # single-factor solution does not have.
  expect_identical(.rel_single_factor_label(NULL), "F1")
  expect_identical(.rel_single_factor_label(character(0)), "F1")
  expect_identical(.rel_single_factor_label(1L), "F1")
  expect_identical(.rel_single_factor_label(seq_len(1)), "F1")
  expect_identical(.rel_single_factor_label(NA_character_), "F1")
  expect_identical(.rel_single_factor_label(""), "F1")

  # More than one name cannot be the label of one factor.
  expect_identical(.rel_single_factor_label(c("F1", "F2")), "F1")
})

rm(efa_mod, sl_mod,
   om_sl_add, om_sl_noadd, om_sl_named_add, om_sl_named_noadd, schmid_mod,
   om_schmid_1_add, om_schmid_1_noadd, om_schmid_2, om_man_1_add, om_man_1_noadd,
   om_man_2)
if (requireNamespace("lavaan", quietly = TRUE)) {
  rm(lav_mod_1, lav_fit_1, om_lav_bi_add, om_lav_bi_noadd, lav_mod_2, lav_fit_2,
     om_lav_1_add, om_lav_1_noadd, lav_fit_3, lav_fit_4, om_lav_gr_add, om_lav_gr_noadd,
     om_lav_1_gr_add, om_lav_1_gr_noadd, lav_mod_ho_1, lav_mod_ho_2, lav_fit_ho_1,
     lav_fit_ho_2, om_lav_ho_1, lav_mod_NA, lav_fit_NA,
     lav_mod_inv, lav_fit_inv, lav_mod_bi_red, lav_fit_bi_red)
}
