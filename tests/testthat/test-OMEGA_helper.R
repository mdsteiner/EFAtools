## Tests for .OMEGA_LAVAAN --------

lav_mod_1 <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
              F2 =~ V7 + V8 + V9 + V10 + V11 + V12
              F3 =~ V13 + V14 + V15 + V16 + V17 + V18
              g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
                   V13 + V14 + V15 + V16 + V17 + V18'
lav_fit_1 <- suppressWarnings(lavaan::cfa(lav_mod_1,
                                          sample.cov = test_models$baseline$cormat,
                                          sample.nobs = 500, estimator = "ml",
                                          orthogonal = TRUE))
om_lav_bi <- .OMEGA_LAVAAN(lav_fit_1, g_name = "g")

lav_mod_2 <- 'g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
                    V13 + V14 + V15 + V16 + V17 + V18'
lav_fit_2 <- suppressWarnings(lavaan::cfa(lav_mod_2,
                                         sample.cov = test_models$baseline$cormat,
                                         sample.nobs = 500, estimator = "ml",
                                         orthogonal = TRUE))
om_lav_1 <- suppressMessages(.OMEGA_LAVAAN(lav_fit_2))

lav_fit_3 <- suppressWarnings(lavaan::cfa(lav_mod_1, sample.cov =
                                            list(test_models$baseline$cormat,
                                                 test_models$baseline$cormat),
                                          sample.nobs = c(500, 500),
                                          estimator = "ml", orthogonal = TRUE))
om_lav_gr <- .OMEGA_LAVAAN(lav_fit_3, g_name = "g", group_names = c("Some",
                                                                    "Others"))

lav_mod_ho <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
               F2 =~ V7 + V8 + V9 + V10 + V11 + V12
               F3 =~ V13 + V14 + V15 + V16 + V17 + V18
               g =~ F1 + F2 + F3'
lav_fit_ho <- suppressWarnings(lavaan::cfa(lav_mod_ho,
                                           sample.cov = test_models$baseline$cormat,
                                           sample.nobs = 500, estimator = "ml"))
om_lav_ho <- suppressMessages(.OMEGA_LAVAAN(lav_fit_ho, g_name = "g"))

test_that("output class and dimensions are correct", {
  expect_is(om_lav_bi, "OMEGA")
  expect_is(om_lav_1, "OMEGA")
  expect_is(om_lav_ho, "OMEGA")
  expect_is(om_lav_gr, "OMEGA")

  expect_output(str(om_lav_bi), "List of 2")
  expect_output(str(om_lav_1), "OMEGA")
  expect_output(str(om_lav_ho), "List of 2")
  expect_output(str(om_lav_gr), "List of 2")
})

test_that("output is correct (including group names for multiple groups)", {
  expect_equal(rowSums(om_lav_bi[, 2:3]), om_lav_bi[, 1], tolerance = 1e-3)
  expect_equal(om_lav_1[1], 0.868, tolerance = 1e-3)
  expect_equal(rowSums(om_lav_ho[, 2:3]), om_lav_ho[, 1], tolerance = 1e-3)
  expect_equal(rowSums(om_lav_gr$Some[, 2:3]), om_lav_gr$Some[, 1], tolerance = 1e-3)
  expect_equal(rowSums(om_lav_gr$Others[, 2:3]), om_lav_gr$Others[, 1], tolerance = 1e-3)
})

# Preparations for error tests
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

test_that("errors are thrown correctly", {
  expect_error(.OMEGA_LAVAAN(lav_fit_NA, g_name = "g"), " Some loadings are NA or NaN. No omegas are computed.\n")
  expect_error(.OMEGA_LAVAAN(lav_fit_1, g_name = "fu"), " Could not find the specified name of the general factor in the entered lavaan solution. Please check the spelling.\n")
  expect_message(.OMEGA_LAVAAN(lav_fit_2), " Model contained a single factor. Only omega total is returned.\n")
  expect_message(.OMEGA_LAVAAN(lav_fit_ho, g_name = "g"), " The general factor you specified is a second-order factor. Omegas are found on the Schmid-Leiman transformed second-order solution.\n")
  expect_error(.OMEGA_LAVAAN(lav_fit_inv, g_name = "F3"), " Your lavaan input is invalid, no omegas are computed. Either provide a bifactor model, a second-order model, or a model with a single factor.\n")
  expect_message(.OMEGA_LAVAAN(lav_fit_bi_red, g_name = "g"), " Some variables have less than two loadings. Did you really enter a bifactor model? Either provide a bifactor model, a second-order model, or a model with a single factor.\n", fixed = TRUE)
})

## Tests for .OMEGA_FLEX -------

## Use with an output from the SL function, with type EFAtools
efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
sl_mod <- SL(efa_mod, type = "EFAtools", method = "PAF")
om_sl <- .OMEGA_FLEX(sl_mod, type = "EFAtools",
                     factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                     variance = "correlation")

# Explicit factor names
om_sl_named <- .OMEGA_FLEX(sl_mod, type = "EFAtools",
                           fac_names = c("Fa1", "Fa2", "Fa3"),
                           factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                           variance = "correlation")

## Use with an output from the psych::schmid function, with type psych
schmid_mod <- psych::schmid(test_models$baseline$cormat, nfactors = 3,
                            n.obs = 500, fm = "pa", rotate = "Promax")
# Omega also with type "psych"
# Find correlation matrix from phi and pattern matrix from psych::schmid output
om_schmid_1 <- .OMEGA_FLEX(schmid_mod, type = "psych", variance = "correlation")
# Enter correlation matrix manually
om_schmid_2 <- .OMEGA_FLEX(schmid_mod, type = "psych",
                           cormat = test_models$baseline$cormat,
                           variance = "correlation")

## Manually specify components (here with type "EFAtools")
om_man_1 <- .OMEGA_FLEX(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                  g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                  u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
                  factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                  variance = "correlation")
# Now with other type of variance (model-based instead of based on corrmat)
om_man_2 <- .OMEGA_FLEX(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                  g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                  u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
                  factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                  variance = "sums_load")

test_that("output class and dimensions are correct", {
  expect_is(om_sl, "OMEGA")
  expect_is(om_sl_named, "OMEGA")
  expect_is(om_schmid_1, "OMEGA")
  expect_is(om_schmid_2, "OMEGA")
  expect_is(om_man_1, "OMEGA")
  expect_is(om_man_2, "OMEGA")

  expect_output(str(om_sl), "List of 2")
  expect_output(str(om_sl_named), "List of 2")
  expect_output(str(om_schmid_1), "List of 2")
  expect_output(str(om_schmid_2), "List of 2")
  expect_output(str(om_man_1), "List of 2")
  expect_output(str(om_man_2), "List of 2")
})

test_that("output is correct", {
  expect_equal(rowSums(om_sl[2:4, 2:3]), om_sl[2:4, 1], tolerance = 1e-3)
  expect_equal(rowSums(om_sl_named[2:4, 2:3]), om_sl_named[2:4, 1], tolerance = 1e-3)
  expect_equal(rowSums(om_schmid_1[2:4, 2:3]), om_schmid_1[2:4, 1], tolerance = 1e-3)
  expect_equal(rowSums(om_schmid_2[2:4, 2:3]), om_schmid_2[2:4, 1], tolerance = 1e-3)
  expect_equal(rowSums(om_man_1[2:4, 2:3]), om_man_1[2:4, 1], tolerance = 1e-3)
  expect_equal(rowSums(om_man_2[2:4, 2:3]), om_man_2[2:4, 1], tolerance = 1e-3)

  expect_gte(om_sl[1, 1], sum(om_sl[1, 2:3]))
  expect_gte(om_sl_named[1, 1], sum(om_sl_named[1, 2:3]))
  expect_gte(om_schmid_1[1, 1], sum(om_schmid_1[1, 2:3]))
  expect_gte(om_schmid_2[1, 1], sum(om_schmid_2[1, 2:3]))
  expect_gte(om_man_1[1, 1], sum(om_man_1[1, 2:3]))
  expect_gte(om_man_2[1, 1], sum(om_man_2[1, 2:3]))
})

test_that("errors are thrown correctly", {
  expect_error(.OMEGA_FLEX(schmid_mod, type = "psych",
                           cormat = matrix(rnorm(50), ncol = 5),
                           variance = "correlation"), " 'x' was not a correlation matrix. Check the cormat input, specify the Phi and pattern arguments instead, or set variance to 'sums_load'\n")
  expect_error(.OMEGA_FLEX(model = NULL, type = "EFAtools",
                           var_names = rownames(sl_mod$sl),
                           g_load = sl_mod$sl[, "g"],
                           s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                           u2 = sl_mod$sl[, "u2"],
                           factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                           variance = "correlation"), " Either specify the cormat argument or the Phi and pattern arguments, or set variance to 'sums_load'\n")
  expect_error(.OMEGA_FLEX(model = NULL, type = "EFAtools",
                           var_names = rownames(sl_mod$sl),
                           g_load = sl_mod$sl[, "g"],
                           s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                           u2 = sl_mod$sl[, "u2"], cormat = matrix(rnorm(50), ncol = 5),
                           factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                           variance = "correlation"), " 'x' was not a correlation matrix. Check the cormat input, specify the Phi and pattern arguments instead, or set variance to 'sums_load'\n")
  expect_error(.OMEGA_FLEX(schmid_mod, type = "EFAtools",
                           variance = "correlation"), " Either specify the factor_corres argument or set type = 'psych' to find variable-to-factor correspondences using the highest group factor loading per variable.\n")
  expect_warning(.OMEGA_FLEX(schmid_mod, type = "psych",
                             variance = "sums_load"), " Variance is specified. Variance is used with value ' sums_load '. Results may differ from the specified type\n")
  expect_warning(.OMEGA_FLEX(schmid_mod, type = "psych",
                             factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                             variance = "correlation"), " Argument factor_corres is specified. Specified variable-to-factor correspondences are taken. To compute factor correspondences as done in psych, leave factor_corres = NULL.\n")
})

rm(lav_mod_1, lav_fit_1, om_lav_bi, lav_mod_2, lav_fit_2, om_lav_1, lav_fit_3,
   om_lav_gr, lav_mod_ho, lav_fit_ho, om_lav_ho, lav_mod_NA, lav_fit_NA,
   lav_mod_inv, lav_fit_inv, lav_mod_bi_red, lav_fit_bi_red, efa_mod, sl_mod,
   om_sl, om_sl_named, schmid_mod, om_schmid_1, om_schmid_2, om_man_1, om_man_2)
