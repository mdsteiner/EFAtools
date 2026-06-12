## Use with a lavaan output
if (requireNamespace("lavaan", quietly = TRUE)) {
lav_mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
        F2 =~ V7 + V8 + V9 + V10 + V11 + V12
        F3 =~ V13 + V14 + V15 + V16 + V17 + V18
        g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
             V13 + V14 + V15 + V16 + V17 + V18'
lav_fit <- lavaan::cfa(lav_mod, sample.cov = test_models$baseline$cormat,
                   sample.nobs = 500, estimator = "ml", orthogonal = TRUE)
om_lav <- OMEGA(lav_fit, g_name = "g")
}

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
  skip_if_not_installed("lavaan")
  expect_s3_class(om_lav, "OMEGA")
  expect_s3_class(om_sl, "OMEGA")
  expect_s3_class(om_schmid, "OMEGA")
  expect_s3_class(om_man, "OMEGA")

  expect_output(str(om_lav), "List of 2")
  expect_output(str(om_sl), "List of 2")
  expect_output(str(om_schmid), "List of 2")
  expect_output(str(om_man), "List of 2")
})

test_that("errors are thrown correctly", {
  expect_error(OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                     g_load = sl_mod$sl[, "g"], s_load = 1:7,
                     u2 = sl_mod$sl[, "u2"], cormat = test_models$baseline$cormat,
                     factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2), class = "efa_omega_bad_s_load")
  expect_error(OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                     g_load = sl_mod$sl[, "g"],
                     s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                     u2 = sl_mod$sl[, "u2"],
                     factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2,
                     pattern = 1:5), class = "efa_omega_bad_pattern")
  expect_warning(OMEGA(sl_mod, type = "EFAtools", g_load = sl_mod$sl[, "g"],
                       s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                       u2 = sl_mod$sl[, "u2"],
                       factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2), class = "efa_omega_model_args_ignored")
  expect_error(OMEGA(model = 1:4), class = "efa_omega_invalid_model")
  expect_error(OMEGA(model = NULL, type = "EFAtools", var_names = rownames(sl_mod$sl),
                     s_load = sl_mod$sl[, c("F1", "F2", "F3")],
                     u2 = sl_mod$sl[, "u2"],
                     factor_corres = sl_mod$sl[, c("F1", "F2", "F3")] >= .2), class = "efa_omega_missing_args")
})

test_that("print output is stable", {
  local_reproducible_output()

  # single group, full coefficient matrix (general factor plus group factors)
  expect_snapshot(print(om_sl), transform = scrub_num)

  # single group, single factor
  expect_snapshot(print(structure(0.85, class = "OMEGA")), transform = scrub_num)
  expect_snapshot(print(structure(c(0.85, 2.10), class = "OMEGA")), transform = scrub_num)

  # multiple groups, full coefficient matrix (long header with ECV and PUC)
  om_mg6 <- structure(list(GroupA = unclass(om_sl), GroupB = unclass(om_sl)),
                      class = "OMEGA")
  expect_snapshot(print(om_mg6), transform = scrub_num)

  # multiple groups, three-column matrix (omega total, hierarchical, subscale)
  om_mg3 <- structure(list(GroupA = unclass(om_sl)[, 1:3],
                           GroupB = unclass(om_sl)[, 1:3]), class = "OMEGA")
  expect_snapshot(print(om_mg3), transform = scrub_num)

  # multiple groups, single factor
  expect_snapshot(print(structure(list(GroupA = 0.85, GroupB = 0.80), class = "OMEGA")),
                  transform = scrub_num)
  expect_snapshot(print(structure(list(GroupA = c(0.85, 2.10), GroupB = c(0.80, 1.90)),
                                  class = "OMEGA")), transform = scrub_num)
})

rm(efa_mod, sl_mod, om_sl, schmid_mod, om_schmid, om_man)
if (requireNamespace("lavaan", quietly = TRUE)) {
  rm(lav_mod, lav_fit, om_lav)
}
