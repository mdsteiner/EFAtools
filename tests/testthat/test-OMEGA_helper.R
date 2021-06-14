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

test_that("output class and dimensions are correct", {
  expect_is(om_lav_bi_add, "OMEGA")
  expect_is(om_lav_1_add, "OMEGA")
  expect_is(om_lav_ho_1, "OMEGA")
  expect_is(om_lav_gr_add, "OMEGA")
  expect_is(om_lav_1_gr_add, "OMEGA")
  expect_is(om_lav_bi_noadd, "OMEGA")
  expect_is(om_lav_1_noadd, "OMEGA")
  expect_is(om_lav_gr_noadd, "OMEGA")
  expect_is(om_lav_1_gr_noadd, "OMEGA")

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

test_that("errors are thrown correctly", {
  expect_error(.OMEGA_LAVAAN(lav_fit_NA, g_name = "g"), " Some loadings are NA or NaN. No omegas are computed.\n")
  expect_error(.OMEGA_LAVAAN(lav_fit_1, g_name = "fu"), " Could not find the specified name of the general factor in the entered lavaan solution. Please check the spelling.\n")
  expect_message(.OMEGA_LAVAAN(lav_fit_2, add_ind = FALSE), " Model contained a single factor. Only omega total is returned.\n")
  expect_message(.OMEGA_LAVAAN(lav_fit_2), " Model contained a single factor. Only omega total and H index are returned.\n")
  expect_message(.OMEGA_LAVAAN(lav_fit_ho_1, g_name = "g"), " The general factor you specified is a second-order factor. Omegas are found on the Schmid-Leiman transformed second-order solution.\n")
  expect_error(.OMEGA_LAVAAN(lav_fit_inv, g_name = "F3"), " Your lavaan input is invalid, no omegas are computed. Either provide a bifactor model, a second-order model, or a model with a single factor.\n")
  expect_message(.OMEGA_LAVAAN(lav_fit_bi_red, g_name = "g"), " Some variables have less than two loadings. Did you really enter a bifactor model? Either provide a bifactor model, a second-order model, or a model with a single factor.\n", fixed = TRUE)
  expect_error(.OMEGA_LAVAAN(lav_fit_ho_2, g_name = "g"), " Your higher-order model had either more than two latent strata or more than 1 second-order factor. This function only works for second-order models with 1 second-order factor.\n")
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
  expect_is(om_sl_add, "OMEGA")
  expect_is(om_sl_named_add, "OMEGA")
  expect_is(om_schmid_1_add, "OMEGA")
  expect_is(om_schmid_2, "OMEGA")
  expect_is(om_man_1_add, "OMEGA")
  expect_is(om_man_2, "OMEGA")
  expect_is(om_sl_noadd, "OMEGA")
  expect_is(om_sl_named_noadd, "OMEGA")
  expect_is(om_schmid_1_noadd, "OMEGA")
  expect_is(om_man_1_noadd, "OMEGA")

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
  expect_equal(om_sl_add[1, 5], 0.659, tolerance = 1e-3)
  expect_equal(om_sl_add[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(rowSums(om_sl_named_add[2:4, 2:3]), om_sl_named_add[2:4, 1],
               tolerance = 1e-3)
  expect_equal(om_sl_named_add[1, 4], 0.842, tolerance = 1e-3)
  expect_equal(om_sl_named_add[1, 5], 0.659, tolerance = 1e-3)
  expect_equal(om_sl_named_add[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(rowSums(om_schmid_1_add[2:4, 2:3]), om_schmid_1_add[2:4, 1],
               tolerance = 1e-3)
  expect_equal(om_schmid_1_add[1, 4], 0.845, tolerance = 1e-3)
  expect_equal(om_schmid_1_add[1, 5], 0.668, tolerance = 1e-3)
  expect_equal(om_schmid_1_add[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(rowSums(om_schmid_2[2:4, 2:3]), om_schmid_2[2:4, 1],
               tolerance = 1e-3)
  expect_equal(om_schmid_2[1, 4], 0.845, tolerance = 1e-3)
  expect_equal(om_schmid_2[1, 5], 0.668, tolerance = 1e-3)
  expect_equal(om_schmid_2[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(rowSums(om_man_1_add[2:4, 2:3]), om_man_1_add[2:4, 1],
               tolerance = 1e-3)
  expect_equal(om_man_1_add[1, 4], 0.842, tolerance = 1e-3)
  expect_equal(om_man_1_add[1, 5], 0.659, tolerance = 1e-3)
  expect_equal(om_man_1_add[1, 6], 0.706, tolerance = 1e-3)
  expect_equal(rowSums(om_man_2[2:4, 2:3]), om_man_2[2:4, 1], tolerance = 1e-3)
  expect_equal(om_man_2[1, 4], 0.842, tolerance = 1e-3)
  expect_equal(om_man_2[1, 5], 0.659, tolerance = 1e-3)
  expect_equal(om_man_2[1, 6], 0.706, tolerance = 1e-3)

  expect_gte(om_sl_add[1, 1], sum(om_sl_add[1, 2:3]))
  expect_gte(om_sl_named_add[1, 1], sum(om_sl_named_add[1, 2:3]))
  expect_gte(om_schmid_1_add[1, 1], sum(om_schmid_1_add[1, 2:3]))
  expect_gte(om_schmid_2[1, 1], sum(om_schmid_2[1, 2:3]))
  expect_gte(om_man_1_add[1, 1], sum(om_man_1_add[1, 2:3]))
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

rm(lav_mod_1, lav_fit_1, om_lav_bi_add, om_lav_bi_noadd, lav_mod_2, lav_fit_2,
   om_lav_1_add, om_lav_1_noadd, lav_fit_3, lav_fit_4, om_lav_gr_add, om_lav_gr_noadd,
   om_lav_1_gr_add, om_lav_1_gr_noadd, lav_mod_ho_1, lav_mod_ho_2, lav_fit_ho_1,
   lav_fit_ho_2, om_lav_ho_1, lav_mod_NA, lav_fit_NA,
   lav_mod_inv, lav_fit_inv, lav_mod_bi_red, lav_fit_bi_red, efa_mod, sl_mod,
   om_sl_add, om_sl_noadd, om_sl_named_add, om_sl_named_noadd, schmid_mod,
   om_schmid_1_add, om_schmid_1_noadd, om_schmid_2, om_man_1_add, om_man_1_noadd,
   om_man_2)
