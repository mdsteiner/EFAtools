## Use with a lavaan output
lav_mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
        F2 =~ V7 + V8 + V9 + V10 + V11 + V12
        F3 =~ V13 + V14 + V15 + V16 + V17 + V18
        g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
             V13 + V14 + V15 + V16 + V17 + V18'
lav_fit <- lavaan::cfa(lav_mod, sample.cov = test_models$baseline$cormat,
                   sample.nobs = 500, estimator = "ml", orthogonal = TRUE)
om_lav <- OMEGA(lav_fit, g_name = "g")

## Use with an output from the SL function, with type EFAtools
efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
sl_mod <- SL(efa_mod, type = "EFAtools", method = "PAF")
om_sl <- OMEGA(sl_mod, type = "EFAtools",
               factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2)

## Use with an output from the psych::schmid function, with type psych
schmid_mod <- psych::schmid(test_models$baseline$cormat, nfactors = 3,
                            n.obs = 500, fm = "pa", rotate = "Promax")
# Find correlation matrix from phi and pattern matrix from psych::schmid outpu
om_schmid <- OMEGA(schmid_mod, type = "psych")

## Manually specify components
om_man <- OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
                factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2)

test_that("output class and dimensions are correct", {
  expect_is(om_lav, "OMEGA")
  expect_is(om_sl, "OMEGA")
  expect_is(om_schmid, "OMEGA")
  expect_is(om_man, "OMEGA")

  expect_output(str(om_lav), "List of 2")
  expect_output(str(om_sl), "List of 2")
  expect_output(str(om_schmid), "List of 2")
  expect_output(str(om_man), "List of 2")
})

test_that("errors are thrown correctly", {
  expect_error(OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                     g_load = sl_mod$sl[, "g"], s_load = 1:7,
                     u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
                     factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2), " Specification of 's_load' was invalid. Please either leave this 'NULL' if you enter a model input or specify a matrix of loadings from a Schmid-Leiman solution of class matrix or SLLOADINGS.\n")
  expect_error(OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                     g_load = sl_mod$sl[, "g"],
                     s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                     u2 = sl_mod$sl[, "u2"],
                     factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                     pattern = 1:5), " Specification of 'pattern' was invalid. Please either leave this NULL or specify a matrix of pattern coefficients form an oblique factor solution of class matrix, loadings, or LOADINGS.\n")
  expect_warning(OMEGA(sl_mod, type = "EFAtools", g_load = sl_mod$sl[, "g"],
                       s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                       u2 = sl_mod$sl[, "u2"],
                       factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2), " You entered a model and specified at least one of the arguments 'var_names', 'g_load', 's_load', or 'u2'. These arguments are ignored. To use specific values for these, leave model = NULL and specify all arguments separately.\n")
  expect_error(OMEGA(model = 1:4), " Invalid input for model. Either enter a lavaan, psych::schmid or SL object or specify the arguments 'var_names', 'g_load', and 's_load'.\n")
  expect_error(OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                     s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                     u2 = sl_mod$sl[, "u2"],
                     factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2), " Please specify all of the following arguments: 'var_names', 'g_load', 's_load', 'u2'\n")
})

rm(lav_mod, lav_fit, om_lav, efa_mod, sl_mod, om_sl, schmid_mod, om_schmid, om_man)
