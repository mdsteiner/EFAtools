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
  expect_error(.OMEGA_LAVAAN(lav_fit_NA, g_name = "g"), class = "efa_omega_na_loadings")
  expect_error(.OMEGA_LAVAAN(lav_fit_1, g_name = "fu"), class = "efa_omega_g_name")
  expect_message(.OMEGA_LAVAAN(lav_fit_2, add_ind = FALSE), class = "efa_omega_single_factor")
  expect_message(.OMEGA_LAVAAN(lav_fit_2), class = "efa_omega_single_factor")
  expect_message(.OMEGA_LAVAAN(lav_fit_ho_1, g_name = "g"), class = "efa_omega_g_second_order")
  expect_error(.OMEGA_LAVAAN(lav_fit_inv, g_name = "F3"), class = "efa_omega_invalid_lavaan")
  expect_message(.OMEGA_LAVAAN(lav_fit_bi_red, g_name = "g"), class = "efa_omega_few_loadings")
  expect_error(.OMEGA_LAVAAN(lav_fit_ho_2, g_name = "g"), class = "efa_omega_higher_order")
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
  expect_equal(rowSums(om_sl_add[2:4, 2:3]), om_sl_add[2:4, 1], tolerance = 1e-3)
  expect_equal(om_sl_add[1, 4], 0.842, tolerance = 1e-3)
  expect_equal(om_sl_add[1, 5], 0.652, tolerance = 1e-3)
  expect_equal(om_sl_add[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(rowSums(om_sl_named_add[2:4, 2:3]), om_sl_named_add[2:4, 1],
               tolerance = 1e-3)
  expect_equal(om_sl_named_add[1, 4], 0.842, tolerance = 1e-3)
  expect_equal(om_sl_named_add[1, 5], 0.652, tolerance = 1e-3)
  expect_equal(om_sl_named_add[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(rowSums(om_schmid_1_add[2:4, 2:3]), om_schmid_1_add[2:4, 1],
               tolerance = 1e-3)
  expect_equal(om_schmid_1_add[1, 4], 0.845, tolerance = 1e-3)
  expect_equal(om_schmid_1_add[1, 5], 0.661, tolerance = 1e-3)
  expect_equal(om_schmid_1_add[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(rowSums(om_schmid_2[2:4, 2:3]), om_schmid_2[2:4, 1],
               tolerance = 1e-3)
  expect_equal(om_schmid_2[1, 4], 0.845, tolerance = 1e-3)
  expect_equal(om_schmid_2[1, 5], 0.661, tolerance = 1e-3)
  expect_equal(om_schmid_2[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(rowSums(om_man_1_add[2:4, 2:3]), om_man_1_add[2:4, 1],
               tolerance = 1e-3)
  expect_equal(om_man_1_add[1, 4], 0.842, tolerance = 1e-3)
  expect_equal(om_man_1_add[1, 5], 0.652, tolerance = 1e-3)
  expect_equal(om_man_1_add[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(rowSums(om_man_2[2:4, 2:3]), om_man_2[2:4, 1], tolerance = 1e-3)
  expect_equal(om_man_2[1, 4], 0.842, tolerance = 1e-3)
  expect_equal(om_man_2[1, 5], 0.652, tolerance = 1e-3)
  expect_equal(om_man_2[1, 6], 0.706, tolerance = 1e-3)

  expect_gte(om_sl_add[1, 1], sum(om_sl_add[1, 2:3]))
  expect_gte(om_sl_named_add[1, 1], sum(om_sl_named_add[1, 2:3]))
  expect_gte(om_schmid_1_add[1, 1], sum(om_schmid_1_add[1, 2:3]))
  expect_gte(om_schmid_2[1, 1], sum(om_schmid_2[1, 2:3]))
  expect_gte(om_man_1_add[1, 1], sum(om_man_1_add[1, 2:3]))
  expect_gte(om_man_2[1, 1], sum(om_man_2[1, 2:3]))

  # With variance = "sums_load", the general-factor omega total equals
  # hierarchical + subscale exactly (shared total-variance denominator).
  expect_equal(unname(om_man_2[1, 1]), unname(sum(om_man_2[1, 2:3])))
})

test_that("type = 'psych' reproduces psych::omega (g omegas and ECV)", {
  skip_on_cran()
  skip_if_not_installed("psych")

  po <- psych::omega(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                     fm = "pa", rotate = "Promax", plot = FALSE)

  # The psych-type path derives the same variable-to-factor correspondences as
  # psych::omega and bases the ECV on the (unzeroed) group loadings, so the
  # general-factor omegas and the ECV reproduce psych::omega.
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
