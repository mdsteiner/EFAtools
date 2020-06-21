#TEST OUTPUT CONTENT AND
# SPECIFIC DETAILS ABOUT INPUTS (also: specific variants for model) IN OMEGA_HELPER

## Tests for .OMEGA_LAVAAN --------

lav_mod_1 <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
        F2 =~ V7 + V8 + V9 + V10 + V11 + V12
        F3 =~ V13 + V14 + V15 + V16 + V17 + V18
        g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
             V13 + V14 + V15 + V16 + V17 + V18'
lav_fit_1 <- lavaan::cfa(lav_mod_1, sample.cov = test_models$baseline$cormat,
                       sample.nobs = 500, estimator = "ml", orthogonal = TRUE)
om_lav_bi <- .OMEGA_LAVAAN(lav_fit_1, g_name = "g")

lav_mod_2 <- 'g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
                    V13 + V14 + V15 + V16 + V17 + V18'
lav_fit_2 <- lavaan::cfa(lav_mod_2, sample.cov = test_models$baseline$cormat,
                       sample.nobs = 500, estimator = "ml", orthogonal = TRUE)
om_lav_1 <- .OMEGA_LAVAAN(lav_fit_2)

lav_fit_3 <- lavaan::cfa(lav_mod_1, sample.cov =
                           list(test_models$baseline$cormat,
                                test_models$baseline$cormat),
                         sample.nobs = c(500, 500), estimator = "ml",
                         orthogonal = TRUE)
om_lav_gr <- .OMEGA_LAVAAN(lav_fit_3, g_name = "g", group_names = c("Some",
                                                                    "Others"))

test_that("output class and dimensions are correct", {
  expect_is(om_lav_bi, "OMEGA")
  expect_is(om_lav_1, "OMEGA")
  expect_is(om_lav_gr, "OMEGA")

  expect_output(str(om_lav_bi), "List of 2")
  expect_output(str(om_lav_1), "OMEGA")
  expect_output(str(om_lav_gr), "List of 2")
})

test_that("output is correct (including group names for multiple groups)", {
  expect_equal(rowSums(om_lav_bi[, 2:3]), om_lav_bi[, 1], tolerance = 1e-3)
  expect_equal(om_lav_1[1], 0.868, tolerance = 1e-3)
  expect_equal(rowSums(om_lav_gr$Some[, 2:3]), om_lav_gr$Some[, 1], tolerance = 1e-3)
  expect_equal(rowSums(om_lav_gr$Others[, 2:3]), om_lav_gr$Others[, 1], tolerance = 1e-3)
})

# Preparations for error tests
lav_mod_NA <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6 + V17
        F2 =~ V7 + V8 + V9 + V10 + V11 + V12 + V2
        F3 =~ V13 + V14 + V15 + V16 + V17 + V18 + V10
        g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
             V13 + V14 + V15 + V16 + V17 + V18'
lav_fit_NA <- lavaan::cfa(lav_mod_NA, sample.cov = test_models$baseline$cormat,
                         sample.nobs = 500, estimator = "ml")

x <- rnorm(10)
y <- rnorm(10)
z <- x + y
a <- rnorm(10)
dat_sing <- matrix(c(a, x, y, z), ncol = 4)
cor_sing <- stats::cor(dat_sing)
colnames(cor_sing) <- c("V1", "V2", "V3", "V4")

lav_mod_conv <- 'F1 =~ V1 + V2
                 F2 =~ V3 + V4
                 g =~ V1 + V2 + V3 + V4'
lav_fit_conv <- lavaan::cfa(lav_mod_conv, sample.cov = cor_sing,
                            sample.nobs = 10, estimator = "ml")

lav_mod_hier <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
               F2 =~ V7 + V8 + V9 + V10 + V11 + V12
               F3 =~ V13 + V14 + V15 + V16 + V17 + V18
               g =~ F1 + F2 + F3'
lav_fit_hier <- lavaan::cfa(lav_mod_hier, sample.cov = test_models$baseline$cormat,
                            sample.nobs = 500, estimator = "ml")

lav_mod_hier2 <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
               F2 =~ V7 + V8 + V9 + V10 + V11 + V12 + V13
               F3 =~ V13 + V14 + V15 + V16 + V17 + V18
               g =~ F1 + F2 + F3'
lav_fit_hier2 <- lavaan::cfa(lav_mod_hier2, sample.cov = test_models$baseline$cormat,
                            sample.nobs = 500, estimator = "ml")

test_that("errors are thrown correctly", {
  expect_error(.OMEGA_LAVAAN(lav_fit_conv, g_name = "g"), " Model did not converge. No omegas are computed.")
  expect_error(.OMEGA_LAVAAN(lav_fit_NA, g_name = "g"), " Some loadings are NA or NaN. No omegas are computed.")
  expect_message(.OMEGA_LAVAAN(lav_fit_2), " Model contained a single factor. Only omega total is returned.")
  expect_error(.OMEGA_LAVAAN(lav_fit_hier, g_name = "g"), " You did not fit a bifactor model. Omegas cannot be computed. Either provide a bifactor model or a model with a single factor.")
  expect_message(.OMEGA_LAVAAN(lav_fit_hier2, g_name = "g"), " Some variables have less than two loadings. Did you really enter a bifactor model? Either provide a bifactor model or a model with a single factor.")
})

## Tests for .OMEGA_FLEX -------

## Use with an output from the SL function, with type EFAtools
efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
sl_mod <- SL(efa_mod, type = "EFAtools", method = "PAF")
om_sl <- OMEGA(sl_mod, type = "EFAtools", factor_corres = rep(c(3, 2, 1),
                                                              each = 6))

## Use with an output from the psych::schmid function, with type psych
schmid_mod <- psych::schmid(test_models$baseline$cormat, nfactors = 3,
                            n.obs = 500, fm = "pa", rotate = "Promax")
# Find correlation matrix from phi and pattern matrix from psych::schmid outpu
om_schmid <- OMEGA(schmid_mod, type = "psych")

## Manually specify components
om_man <- OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
                factor_corres = rep(c(3, 2, 1), each = 6))


rm(lav_mod, lav_fit, om_lav, efa_mod, sl_mod, om_sl, schmid_mod, om_schmid, om_man)
