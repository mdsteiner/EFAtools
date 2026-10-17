## Use with an output from the EFAtools::EFA function, both with type EFAtools
EFA_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
SL_EFAtools <- efa_schmid_leiman(EFA_mod, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF")

# with type SPSS and method ULS. The second-order fit notes that only PAF is validated
# against SPSS; that note is asserted where it belongs (test-superseded.R) and is muffled
# here so the fixture does not raise a warning outside any test.
SL_SPSS <- suppressWarnings(
  efa_schmid_leiman(EFA_mod, estimate_control = estimate_control(type = "SPSS"),
                    estimator = "ULS"))

## Use with an output from the psych::fa function with type psych
fa_mod <- psych::fa(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                    fm = "pa", rotate = "Promax", n.rotations = 1)
SL_psych <- efa_schmid_leiman(fa_mod, estimate_control = estimate_control(type = "psych"), estimator = "PAF")

## Use more flexibly by entering a pattern matrix and phi directly, with method
## ML
SL_flex <- efa_schmid_leiman(EFA_mod$rot_loadings, Phi = EFA_mod$Phi, estimator = "ML",
                             estimate_control = estimate_control(type = "EFAtools"))

## Use with a second-order lavaan solution
if (requireNamespace("lavaan", quietly = TRUE)) {
lav_mod_ho <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
               F2 =~ V7 + V8 + V9 + V10 + V11 + V12
               F3 =~ V13 + V14 + V15 + V16 + V17 + V18
               g =~ F1 + F2 + F3'
lav_fit_ho <- suppressWarnings(lavaan::cfa(lav_mod_ho,
                                           sample.cov = test_models$baseline$cormat,
                                           sample.nobs = 500, estimator = "ml"))
SL_lav <- efa_schmid_leiman(lav_fit_ho, g_name = "g")
}

test_that("output class and dimensions are correct", {
  skip_if_not_installed("lavaan")
  expect_identical(class(SL_EFAtools), c("efa_schmid_leiman", "SL"))
  expect_identical(class(SL_EFAtools$sl), c("efa_sl_loadings", "SLLOADINGS"))
  expect_s3_class(SL_SPSS, "efa_schmid_leiman")
  expect_s3_class(SL_psych, "efa_schmid_leiman")
  expect_s3_class(SL_flex, "efa_schmid_leiman")
  expect_s3_class(SL_lav, "efa_schmid_leiman")

  expect_output(str(SL_EFAtools), "List of 7")
  expect_output(str(SL_SPSS), "List of 7")
  expect_output(str(SL_psych), "List of 7")
  expect_output(str(SL_flex), "List of 7")
  expect_output(str(SL_lav), "List of 7")
})

test_that("non-convergence of the second-order EFA is recorded and surfaced", {
  testthat::local_reproducible_output()

  # The convergence code of the second-order EFA is carried on the SL object.
  expect_equal(SL_SPSS$convergence, 0)

  # A non-zero code shows the banner: the optimiser estimators use the generic
  # message, PAF the iteration-specific one. Set the code explicitly so the test
  # does not depend on the second-order fit's actual convergence.
  sl_uls <- SL_SPSS
  sl_uls$convergence <- 1L
  expect_true(any(grepl("The optimiser did not converge", format(sl_uls), fixed = TRUE)))

  sl_paf <- SL_EFAtools
  sl_paf$convergence <- 1L
  expect_true(any(grepl("Maximum number of iterations reached without convergence",
                        format(sl_paf), fixed = TRUE)))

  # A converged fit shows no banner.
  expect_false(any(grepl("did not converge", format(SL_SPSS), fixed = TRUE)))
})

test_that("original correlation is correct", {
  skip_if_not_installed("lavaan")
  expect_equal(SL_EFAtools$orig_R, test_models$baseline$cormat)
  expect_equal(SL_SPSS$orig_R, test_models$baseline$cormat)
  expect_equal(SL_psych$orig_R, test_models$baseline$cormat)
  expect_equal(SL_flex$orig_R, NA)
  expect_equal(SL_lav$orig_R, test_models$baseline$cormat)
})

test_that("sl solution is correct", {
  skip_if_not_installed("lavaan")
  expect_equal(unname(SL_EFAtools$sl[, "h2"]) + unname(SL_EFAtools$sl[, "u2"]),
               rep(1, 18))
  expect_equal(unname(SL_SPSS$sl[, "h2"]) + unname(SL_SPSS$sl[, "u2"]),
               rep(1, 18))
  expect_equal(unname(SL_psych$sl[, "h2"]) + unname(SL_psych$sl[, "u2"]),
               rep(1, 18))
  expect_equal(unname(SL_flex$sl[, "h2"]) + unname(SL_flex$sl[, "u2"]),
               rep(1, 18))
  expect_equal(unname(SL_lav$sl[, "h2"]) + unname(SL_lav$sl[, "u2"]),
               rep(1, 18))

  expect_equal(unname(SL_EFAtools$sl[, "g"]) >= .20, rep(TRUE, 18))
  expect_equal(unname(SL_SPSS$sl[, "g"]) >= .20, rep(TRUE, 18))
  expect_equal(unname(SL_psych$sl[, "g"]) >= .20, rep(TRUE, 18))
  expect_equal(unname(SL_flex$sl[, "g"]) >= .20, rep(TRUE, 18))
  expect_equal(unname(SL_lav$sl[, "g"]) >= .20, rep(TRUE, 18))

  expect_equal(unname(SL_EFAtools$sl[13:18, "F1"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_SPSS$sl[13:18, "F1"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_psych$sl[13:18, "PA1"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_flex$sl[13:18, "F1"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_lav$sl[1:6, "F1"]) >= .20, rep(TRUE, 6))

  expect_equal(unname(SL_EFAtools$sl[1:12, "F1"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_SPSS$sl[1:12, "F1"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_psych$sl[1:12, "PA1"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_flex$sl[1:12, "F1"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_lav$sl[7:18, "F1"]) < .20, rep(TRUE, 12))

  expect_equal(unname(SL_EFAtools$sl[7:12, "F2"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_SPSS$sl[7:12, "F2"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_psych$sl[7:12, "PA2"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_flex$sl[7:12, "F2"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_lav$sl[7:12, "F2"]) >= .20, rep(TRUE, 6))

  expect_equal(unname(SL_EFAtools$sl[c(1:6, 13:18), "F2"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_SPSS$sl[c(1:6, 13:18), "F2"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_psych$sl[c(1:6, 13:18), "PA2"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_flex$sl[c(1:6, 13:18), "F2"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_lav$sl[c(1:6, 13:18), "F2"]) < .20, rep(TRUE, 12))

  expect_equal(unname(SL_EFAtools$sl[1:6, "F3"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_SPSS$sl[1:6, "F3"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_psych$sl[1:6, "PA3"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_flex$sl[1:6, "F3"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_lav$sl[13:18, "F3"]) >= .20, rep(TRUE, 6))

  expect_equal(unname(SL_EFAtools$sl[7:18, "F3"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_SPSS$sl[7:18, "F3"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_psych$sl[7:18, "PA3"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_flex$sl[7:18, "F3"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_lav$sl[1:12, "F3"]) < .20, rep(TRUE, 12))
})

test_that("sl loadings reproduce psych::schmid", {
  skip_on_cran()
  skip_if_not_installed("psych")

  # The Schmid-Leiman transformation of the same psych Promax solution must match
  # psych::schmid()'s g and group loadings (its columns are g, F1*, F2*, F3*).
  sch <- suppressMessages(suppressWarnings(
    psych::schmid(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                  fm = "pa", rotate = "Promax")))

  expect_equal(unname(SL_psych$sl[, c("g", "PA1", "PA2", "PA3")]),
               unname(sch$sl[, c("g", "F1*", "F2*", "F3*")]),
               tolerance = 1e-4)
})

test_that("the group factors keep the order and the names of a psych::fa() solution", {
  skip_on_cran()
  skip_if_not_installed("psych")

  # psych::fa() orders its columns by explained variance but keeps each factor's number in
  # its label, so the labels arrive out of order. The columns must stay where they are, so
  # that they match the fa() solution and psych::schmid() column for column.
  fa_def <- psych::fa(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                      fm = "pa", rotate = "Promax")
  sch <- suppressMessages(suppressWarnings(
    psych::schmid(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                  fm = "pa", rotate = "Promax")))
  sl <- efa_schmid_leiman(fa_def, estimate_control = estimate_control(type = "psych"))

  expect_identical(colnames(sl$sl), c("g", colnames(fa_def$loadings), "h2", "u2"))
  expect_identical(rownames(sl$L2), colnames(fa_def$loadings))
  expect_equal(unname(unclass(sl$sl)[, 1:4]), unname(sch$sl[, 1:4]), tolerance = 1e-4)
})

test_that("settings are returned correctly", {
  skip_if_not_installed("lavaan")
  # only the estimator and the estimation settings: the sample size and the inference
  # settings of the second-order EFA are placeholders
  paf_knobs <- c("max_iter", "init_comm", "criterion", "criterion_type", "abs_eigen")
  expect_named(SL_EFAtools$settings, c("estimator", "type", paf_knobs))
  expect_named(SL_SPSS$settings, c("estimator", "type"))
  expect_named(SL_psych$settings, c("estimator", "type", paf_knobs))
  expect_named(SL_flex$settings, c("estimator", "type", "start_method"))
  expect_equal(SL_lav$settings, NA)

  expect_equal(SL_EFAtools$settings$estimator, "PAF")
  expect_equal(SL_SPSS$settings$estimator, "ULS")
  expect_equal(SL_psych$settings$estimator, "PAF")
  expect_equal(SL_flex$settings$estimator, "ML")

  expect_equal(SL_EFAtools$settings$type, "EFAtools")
  expect_equal(SL_SPSS$settings$type, "SPSS")
  expect_equal(SL_psych$settings$type, "psych")
  expect_equal(SL_flex$settings$type, "EFAtools")

  expect_equal(SL_EFAtools$settings$max_iter, 300)
  expect_equal(SL_psych$settings$max_iter, 50)

  expect_equal(SL_EFAtools$settings$init_comm, "smc")
  expect_equal(SL_psych$settings$init_comm, "smc")

  expect_equal(SL_EFAtools$settings$criterion, 0.001)
  expect_equal(SL_psych$settings$criterion,  0.001)

  expect_equal(SL_EFAtools$settings$criterion_type, "sum")
  expect_equal(SL_psych$settings$criterion_type, "sum")

  expect_equal(SL_EFAtools$settings$abs_eigen, TRUE)
  expect_equal(SL_psych$settings$abs_eigen, FALSE)

  expect_equal(SL_flex$settings$start_method, "psych")
})


EFA_mod_unrot <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
                     type = "EFAtools", method = "PAF", rotation = "none")
EFA_mod_orth <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
                     type = "EFAtools", method = "PAF", rotation = "varimax")
fa_mod_unrot <- psych::fa(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                          fm = "pa", rotate = "none")
fa_mod_orth <- psych::fa(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                         fm = "pa", rotate = "Varimax", n.rotations = 1)

if (requireNamespace("lavaan", quietly = TRUE)) {
lav_mod_NA <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6 + V17
               F2 =~ V7 + V8 + V9 + V10 + V11 + V12 + V2
               F3 =~ V13 + V14 + V15 + V16 + V17 + V18 + V10
               g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 +
                    V13 + V14 + V15 + V16 + V17 + V18'
lav_fit_NA <- suppressWarnings(lavaan::cfa(lav_mod_NA,
                                           sample.cov = test_models$baseline$cormat,
                                           sample.nobs = 500, estimator = "ml"))

lav_mod_ho_inv <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
                   F2 =~ V7 + V8 + V9 + V10 + V11 + V12
                   g =~ F1 + F2 + V13 + V14 + V15 + V16 + V17 + V18'
lav_fit_ho_inv <- suppressWarnings(lavaan::cfa(lav_mod_ho_inv,
                                           sample.cov = test_models$baseline$cormat,
                                           sample.nobs = 500, estimator = "ml"))

# Improper second-order solution: a first-order factor residual variance turns
# negative while all loadings stay below 1 (a variance Heywood case).
lav_mod_var_hey <- 'F1 =~ V1 + V2
                    F2 =~ V3 + V4
                    F3 =~ V5 + V6
                    g =~ F1 + F2 + F3'
lav_fit_var_hey <- suppressWarnings(lavaan::cfa(lav_mod_var_hey,
                                           sample.cov = test_models$baseline$cormat,
                                           sample.nobs = 500, estimator = "ml"))

# Bifactor solution (general factor loads on the items, not the first-order
# factors): SL needs a second-order model, so this must abort cleanly.
lav_mod_bifactor <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
                     F2 =~ V7 + V8 + V9 + V10 + V11 + V12
                     F3 =~ V13 + V14 + V15 + V16 + V17 + V18
                     g =~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 +
                          V12 + V13 + V14 + V15 + V16 + V17 + V18'
lav_fit_bifactor <- suppressWarnings(lavaan::cfa(lav_mod_bifactor,
                                           sample.cov = test_models$baseline$cormat,
                                           sample.nobs = 500, estimator = "ml",
                                           orthogonal = TRUE))
}

test_that("errors are thrown correctly", {
  skip_if_not_installed("lavaan")
  expect_error(efa_schmid_leiman(1:5), class = "efa_sl_bad_input")
  expect_error(efa_schmid_leiman(EFA_mod, Phi = as.data.frame(EFA_mod$Phi)),
               class = "efa_invalid_argument")
  expect_error(efa_schmid_leiman(EFA_mod_unrot, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF"), class = "efa_sl_not_oblique")
  expect_error(efa_schmid_leiman(EFA_mod_orth, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF"), class = "efa_sl_not_oblique")
  expect_error(efa_schmid_leiman(fa_mod_unrot, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF"), class = "efa_sl_not_oblique")
  expect_error(efa_schmid_leiman(fa_mod_orth, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF"), class = "efa_sl_not_oblique")
  expect_error(efa_schmid_leiman(lav_fit_NA, g_name = "g"), class = "efa_sl_na_loadings")
  expect_error(efa_schmid_leiman(lav_fit_var_hey, g_name = "g"), class = "efa_sl_heywood")
  expect_error(efa_schmid_leiman(lav_fit_ho, g_name = "fu"), class = "efa_sl_g_name")
  expect_error(efa_schmid_leiman(EFA_mod$rot_loadings, estimate_control = estimate_control(type = "EFAtools"), estimator = "ML"), class = "efa_sl_phi_missing")
  expect_error(efa_schmid_leiman(lav_fit_bifactor, g_name = "g"), class = "efa_sl_not_second_order")
  expect_error(efa_schmid_leiman(lav_fit_bifactor), class = "efa_sl_not_second_order")

  # an oblique bifactor solution carries Phi, but it already has a general factor
  bq <- suppressWarnings(suppressMessages(
    efa_fit(test_models$baseline$cormat, N = 500, n_factors = 4, rotation = "bifactorQ",
            seed = 1)))
  expect_error(efa_schmid_leiman(bq), class = "efa_sl_bifactor")
})

test_that("a Phi given with a fitted solution replaces the fitted one, without a warning", {
  Phi_alt <- EFA_mod$Phi * .8
  diag(Phi_alt) <- 1

  expect_no_warning(sl_alt <- efa_schmid_leiman(EFA_mod, Phi = Phi_alt))
  sl_flex_alt <- efa_schmid_leiman(EFA_mod$rot_loadings, Phi = Phi_alt)
  expect_equal(sl_alt$sl, sl_flex_alt$sl)
  expect_equal(sl_alt$L2, sl_flex_alt$L2)
  expect_false(isTRUE(all.equal(sl_alt$sl, efa_schmid_leiman(EFA_mod)$sl)))
})

test_that("a lavaan general factor with direct item loadings keeps them", {
  skip_if_not_installed("lavaan")

  # g loads on F1, F2, and directly on V13 to V18. Its loadings are the direct path plus the
  # paths through the first-order factors, so every communality is the lavaan R-squared.
  expect_no_warning(sl <- efa_schmid_leiman(lav_fit_ho_inv))
  r2 <- lavaan::lavInspect(lav_fit_ho_inv, "r2")
  expect_equal(unname(sl$sl[, "h2"]), unname(r2[rownames(sl$sl)]), tolerance = 1e-8)
  std <- lavaan::lavInspect(lav_fit_ho_inv, "std")
  expect_equal(unname(sl$sl[13:18, "g"]), unname(std$lambda[13:18, "g"]), tolerance = 1e-8)

  # the reliability coefficients of the transformed solution then agree with those of the
  # lavaan fit, which efa_reliability() transforms itself
  res_sl <- efa_reliability(sl, variance = "sums_load")
  res_lav <- efa_reliability(lav_fit_ho_inv)
  omegas <- res_sl$coefficient %in% c("omega_total", "omega_hierarchical", "omega_subscale")
  expect_equal(res_sl$value[omegas], res_lav$value[omegas], tolerance = 1e-8)

  # and so do those of OMEGA(): the items that load only on g belong to no group factor
  expect_warning(om_sl <- OMEGA(sl, type = "psych", variance = "sums_load"),
                 class = "efa_omega_variance_override")
  om_lav <- suppressMessages(OMEGA(lav_fit_ho_inv, g_name = "g"))
  expect_equal(unclass(om_sl)[, c("tot", "hier", "sub")],
               unclass(om_lav)[, c("tot", "hier", "sub")], tolerance = 1e-8)
})

test_that("a lavaan solution gives its correlation matrix and finds its general factor", {
  skip_if_not_installed("lavaan")

  # the correlation matrix the model was fitted to is kept, so efa_reliability() works with
  # its default variance = "correlation"
  expect_no_error(efa_reliability(SL_lav))

  # the general factor is found without g_name, whatever its name
  mod_G <- sub("g =~", "G =~", lav_mod_ho)
  fit_G <- lavaan::cfa(mod_G, sample.cov = test_models$baseline$cormat, sample.nobs = 500)
  sl_G <- efa_schmid_leiman(fit_G)
  expect_equal(unclass(sl_G$sl), unclass(SL_lav$sl), tolerance = 1e-6)
  expect_identical(rownames(sl_G$L2), c("F1", "F2", "F3"))

  # the first-order factors keep their lavaan names
  mod_named <- gsub("F1", "Verbal", lav_mod_ho)
  fit_named <- lavaan::cfa(mod_named, sample.cov = test_models$baseline$cormat,
                           sample.nobs = 500)
  expect_identical(colnames(efa_schmid_leiman(fit_named)$sl),
                   c("g", "Verbal", "F2", "F3", "h2", "u2"))

  # The transformation ignores any other structure, so a model with other latent paths or
  # correlated disturbances is refused: here two second-order factors, a path between two
  # first-order factors, and a disturbance covariance.
  mod_two <- 'F1 =~ V1 + V2 + V3
              F2 =~ V4 + V5 + V6
              F3 =~ V7 + V8 + V9
              F4 =~ V10 + V11 + V12
              F5 =~ V13 + V14 + V15
              F6 =~ V16 + V17 + V18
              H1 =~ F1 + F2 + F3
              H2 =~ F4 + F5 + F6'
  fit_two <- suppressWarnings(
    lavaan::cfa(mod_two, sample.cov = test_models$baseline$cormat, sample.nobs = 500))
  expect_error(efa_schmid_leiman(fit_two), class = "efa_sl_higher_order")
  expect_error(efa_schmid_leiman(fit_two, g_name = "H1"), class = "efa_sl_higher_order")
  mod_six <- sub("H1 =~ F1 + F2 + F3\n              H2 =~ F4 + F5 + F6",
                 "g =~ F1 + F2 + F3 + F4 + F5 + F6", mod_two, fixed = TRUE)
  fit_with <- function(extra) {
    suppressWarnings(lavaan::cfa(paste(mod_six, extra, sep = "\n"),
                                 sample.cov = test_models$baseline$cormat,
                                 sample.nobs = 500))
  }
  expect_s3_class(efa_schmid_leiman(fit_with("")), "efa_schmid_leiman")
  expect_error(efa_schmid_leiman(fit_with("F1 ~ F3")), class = "efa_sl_higher_order")
  expect_error(efa_schmid_leiman(fit_with("F1 ~~ F3")), class = "efa_sl_correlated_factors")
  # a negative path counts as well, also for efa_reliability()
  fit_neg <- fit_with("F1 ~ -0.2*F3")
  expect_error(efa_schmid_leiman(fit_neg), class = "efa_sl_higher_order")
  expect_error(efa_reliability(fit_neg), class = "efa_reliability_higher_order")

  # a cross-loading item is assigned to every factor it loads on, as on the lavaan route
  fit_cross <- suppressWarnings(lavaan::cfa(sub("V6\n", "V6 + V7\n", lav_mod_ho, fixed = TRUE),
                                            sample.cov = test_models$baseline$cormat,
                                            sample.nobs = 500))
  map_sl <- attr(efa_reliability(efa_schmid_leiman(fit_cross)), "factor_map")
  map_lav <- attr(efa_reliability(fit_cross), "factor_map")
  expect_true(all(map_sl["V7", c("F1", "F2")]))
  expect_identical(unname(map_sl), unname(map_lav))

  # a multiple-group fit has one solution per group
  cm <- test_models$baseline$cormat
  fit_mg <- suppressWarnings(lavaan::cfa(lav_mod_ho, sample.cov = list(cm, cm),
                                         sample.nobs = c(250, 250)))
  expect_error(efa_schmid_leiman(fit_mg), class = "efa_sl_multigroup")
})

test_that("a standardized lavaan loading above 1 is not taken as a Heywood case", {
  skip_if_not_installed("lavaan")

  # With correlated factors, a standardized pattern coefficient is a regression weight and
  # can exceed 1 in a proper solution (Joreskog, 1999). Here x11 loads 1.2 on F1 and -0.5
  # on F2, and its residual variance is positive.
  Lam <- matrix(0, 11, 3)
  Lam[1:4, 1] <- Lam[5:8, 2] <- Lam[9:10, 3] <- .7
  Lam[11, 1:2] <- c(1.2, -.5)
  Phi_p <- tcrossprod(c(.9, .9, .6))
  diag(Phi_p) <- 1
  C <- Lam %*% Phi_p %*% t(Lam)
  diag(C) <- 1
  dimnames(C) <- list(paste0("x", 1:11), paste0("x", 1:11))
  mod <- 'F1 =~ x1 + x2 + x3 + x4 + x11
          F2 =~ x5 + x6 + x7 + x8 + x11
          F3 =~ x9 + x10
          g =~ F1 + F2 + F3'
  fit <- suppressWarnings(lavaan::cfa(mod, sample.cov = C, sample.nobs = 1000))
  expect_gt(lavaan::lavInspect(fit, "std")$lambda["x11", "F1"], 1)

  sl <- efa_schmid_leiman(fit)
  expect_equal(unname(sl$sl["x11", "h2"]),
               unname(lavaan::lavInspect(fit, "r2")["x11"]), tolerance = 1e-6)
})

test_that("the second-order fit's substantive warnings reach the user", {
  # A factor space with an analytic Heywood case: for one factor on three variables
  # lambda1^2 = r12 * r13 / r23, so these intercorrelations imply a communality of 1.44.
  # PAF and ULS abort on it; ML clamps the communality at the boundary and warns, and
  # that warning must not be swallowed -- the residualized loadings it leads to are
  # essentially zero and everything computed from them rests on it.
  Phi_hey <- matrix(c(1, .9, .8, .9, 1, .5, .8, .5, 1), 3, 3)
  L1 <- matrix(0, 9, 3)
  L1[1:3, 1] <- L1[4:6, 2] <- L1[7:9, 3] <- .7
  colnames(L1) <- c("F1", "F2", "F3")

  expect_warning(efa_schmid_leiman(L1, Phi = Phi_hey, estimator = "ML"),
                 class = "efa_heywood")
  # PAF aborts on the same input, after the fit has named the improper first-order
  # factor -- the abort says that there is a Heywood case, the warning says where.
  expect_warning(
    expect_error(efa_schmid_leiman(L1, Phi = Phi_hey, estimator = "PAF"),
                 class = "efa_sl_heywood"),
    class = "efa_heywood")

  # The identification warnings of the internal fit are structural -- one factor on a
  # k by k matrix of factor intercorrelations is just identified at k = 3 for every
  # input -- so a sound three-factor solution stays silent.
  expect_no_warning(efa_schmid_leiman(EFA_mod, estimator = "PAF"))
  expect_no_warning(efa_schmid_leiman(EFA_mod, estimator = "ML"))
})

test_that("too few first-order factors are reported on the user's solution", {
  # One first-order factor leaves nothing to orthogonalize. The diagnostic must name the
  # supplied solution rather than the internal one-factor fit on the 1 x 1 matrix of
  # factor intercorrelations, which would otherwise abort about a single "variable".
  expect_error(
    efa_schmid_leiman(matrix(c(.7, .6, .5, .4), 4, 1), Phi = matrix(1, 1, 1),
                      estimator = "PAF"),
    class = "efa_sl_too_few_factors")

  # Two first-order factors do not identify the second-order loadings. As in psych::schmid(),
  # they get the same absolute value, with the factor intercorrelation as their product, and
  # no second-order EFA is run.
  L2 <- matrix(0, 4, 2)
  L2[1:2, 1] <- 0.7
  L2[3:4, 2] <- 0.7
  colnames(L2) <- c("F1", "F2")
  expect_warning(
    sl2 <- efa_schmid_leiman(L2, Phi = matrix(c(1, .5, .5, 1), 2), estimator = "ML"),
    class = "efa_sl_underidentified")
  expect_equal(unname(sl2$L2[, 1]), rep(sqrt(.5), 2))
  expect_true(is.na(sl2$settings) && is.na(sl2$iter) && is.na(sl2$convergence))

  expect_warning(
    sl2_neg <- efa_schmid_leiman(L2, Phi = matrix(c(1, -.5, -.5, 1), 2)),
    class = "efa_sl_underidentified")
  expect_equal(unname(sl2_neg$L2[, 1]), c(sqrt(.5), -sqrt(.5)))
})

test_that("efa_fit()'s inference arguments are refused by SL()", {
  # The second-order EFA is internal: it runs against a placeholder N and keeps only its
  # loadings, so efa_schmid_leiman() takes no efa_fit() argument other than the estimator
  # and the estimation control.
  expect_false("..." %in% names(formals(efa_schmid_leiman)))

  # The frozen wrapper keeps its dots, and refuses the inference arguments there, also as an
  # abbreviation that R would partial-match to them.
  expect_error(SL(EFA_mod, type = "EFAtools", method = "ML", se = "information"),
               class = "efa_unused_dots")
  expect_error(SL(EFA_mod, type = "EFAtools", method = "ML", seed = 42),
               class = "efa_unused_dots")
  expect_error(SL(EFA_mod, type = "EFAtools", method = "ML", b_b = 10),
               class = "efa_unused_dots")

  # any other efa_fit() argument had no effect on the internal EFA and is not used
  expect_equal(SL(EFA_mod, type = "EFAtools", method = "ML", use = "complete.obs"),
               SL(EFA_mod, type = "EFAtools", method = "ML"))
})

test_that("a second-order Heywood case raises a classed error", {
  # A pattern matrix with a highly (and unevenly) intercorrelated factor space
  # makes the single second-order factor improper: a second-order communality
  # exceeds 1, so the residualized first-order loadings would be undefined.
  L1 <- matrix(0, 6, 3)
  L1[1:2, 1] <- 0.7
  L1[3:4, 2] <- 0.7
  L1[5:6, 3] <- 0.7
  colnames(L1) <- c("F1", "F2", "F3")
  Phi_heywood <- matrix(c(1, 0.8, 0.95,
                          0.8, 1, 0.7,
                          0.95, 0.7, 1), nrow = 3)

  # The second-order fit names the improper first-order factor before the abort.
  expect_warning(
    expect_error(efa_schmid_leiman(L1, Phi = Phi_heywood, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF"),
                 class = "efa_sl_heywood"),
    class = "efa_heywood")
})

test_that("oblique solutions with unlabelled factor columns keep all factors", {

  # A loading matrix whose columns carry no labels leaves the first-order factors to be
  # taken in the order they come in. Reordering them by a factor number parsed from
  # absent labels used to subset the loadings and Phi down to zero columns, which made
  # the second-order EFA see no variables. efa_fit() labels every rotation's columns, so
  # strip the labels to reach the path an externally built or bare-matrix input takes.
  cormat <- test_models$baseline$cormat
  cormat_list <- list(cormat, cormat, cormat)

  for (k in c(3L, 4L)) {

    fitted <- suppressMessages(efa_fit(cormat, N = 500, n_factors = k,
                                       estimator = "PAF", rotation = "oblimin"))
    pooled <- suppressMessages(efa_mi(cormat_list, N = 500, n_factors = k,
                                      estimator = "PAF", rotation = "oblimin"))

    # the fixtures only test what they are meant to if the columns are unlabelled
    colnames(fitted$rot_loadings) <- NULL
    colnames(pooled$rot_loadings) <- NULL
    dimnames(fitted$Phi) <- NULL
    dimnames(pooled$Phi) <- NULL

    for (obj in list(fitted, pooled)) {

      sl <- efa_schmid_leiman(obj, estimator = "PAF")

      expect_s3_class(sl, "efa_schmid_leiman")
      expect_identical(colnames(sl$sl),
                       c("g", paste0("F", seq_len(k)), "h2", "u2"))
      expect_identical(nrow(sl$sl), ncol(cormat))
      expect_false(anyNA(sl$sl))

      # the second-order EFA saw all k first-order factors: this is what the
      # zero-column subset destroyed, so it is the assertion that pins the bug
      expect_identical(dim(sl$L2), c(k, 1L))

      # loadings are real rather than the all-zero block a collapsed L1 produces
      expect_true(all(abs(sl$sl[, "g"]) > 0))
      expect_true(all(sl$sl[, "h2"] > 0 & sl$sl[, "h2"] < 1))

    }

  }

})

test_that("labelled factor columns keep their order and names", {

  # Columns whose labels are out of order are not sorted, so an indicator-to-factor map
  # built in the column order of the input still lines up with the result.
  perm <- c(2, 3, 1)
  L1 <- EFA_mod$rot_loadings[, perm]
  Phi_1 <- EFA_mod$Phi[perm, perm]

  sl_perm <- efa_schmid_leiman(L1, Phi = Phi_1, estimator = "PAF",
                               estimate_control = estimate_control(type = "EFAtools"))

  expect_identical(colnames(sl_perm$sl), c("g", "F2", "F3", "F1", "h2", "u2"))
  expect_equal(unclass(sl_perm$sl), unclass(SL_EFAtools$sl)[, c(1, perm + 1, 5, 6)])

  # the messages of the second-order EFA name the first-order factors, also when the
  # pattern matrix and Phi have no names
  sl_bare <- efa_schmid_leiman(unname(unclass(EFA_mod$rot_loadings)),
                               Phi = unname(EFA_mod$Phi))
  expect_identical(rownames(sl_bare$L2), c("F1", "F2", "F3"))
})

test_that("a named Phi is aligned to the loading factors", {
  L1 <- unclass(EFA_mod$rot_loadings)
  Phi <- EFA_mod$Phi
  base <- efa_schmid_leiman(
    L1, Phi = Phi, estimator = "PAF",
    estimate_control = estimate_control(type = "EFAtools")
  )

  perm <- c(3L, 1L, 2L)
  aligned <- efa_schmid_leiman(
    L1, Phi = Phi[perm, perm], estimator = "PAF",
    estimate_control = estimate_control(type = "EFAtools")
  )

  expect_equal(aligned$sl, base$sl)
  expect_equal(aligned$L2, base$L2)
})

test_that("print output is stable", {
  local_reproducible_output()

  expect_snapshot(print(SL_EFAtools), transform = scrub_num)
})

test_that("the second-order header follows the console width", {
  # Emitted verbatim so the "estimator = 'value'" token is never split, and packed to the
  # console around it -- the same treatment the efa_fit() header gets.
  header <- function(w) {
    out <- withr::with_options(list(cli.width = w, cli.num_colors = 1), format(SL_EFAtools))
    wrapped_item(out, "^EFA for second-order")
  }

  for (w in c(120L, 80L, 60L)) {
    lines <- header(w)
    expect_true(all(cli::ansi_nchar(lines, type = "width") <= w))
    expect_identical(paste(trimws(lines), collapse = " "),
                     "EFA for second-order loadings performed with estimator = 'PAF'")
  }
  expect_length(header(60L), 2L)
})

rm(EFA_mod, SL_EFAtools, SL_SPSS, fa_mod, SL_psych, SL_flex, EFA_mod_unrot,
   EFA_mod_orth, fa_mod_unrot, fa_mod_orth)
if (requireNamespace("lavaan", quietly = TRUE)) {
  rm(lav_mod_ho, lav_fit_ho, SL_lav, lav_mod_NA, lav_fit_NA, lav_mod_ho_inv,
     lav_fit_ho_inv, lav_mod_var_hey, lav_fit_var_hey, lav_mod_bifactor,
     lav_fit_bifactor)
}


test_that("the input error names the averaged matrix an efa_average object carries", {
  # An efa_average object is rejected while "a matrix" is listed as acceptable, and its
  # averaged loadings are exactly that. The message now says where they are, together with
  # the averaged factor correlations, which this transformation always needs.
  avg <- suppressWarnings(suppressMessages(
    efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                estimator = c("PAF", "ULS"), rotation = "promax", type = "EFAtools",
                start_method = "psych", show_progress = FALSE)))

  e <- tryCatch(efa_schmid_leiman(avg), efa_sl_bad_input = function(e) e)
  expect_s3_class(e, "efa_sl_bad_input")
  expect_snapshot(cat(conditionMessage(e)))

  # The components the message names are accepted.
  expect_s3_class(
    suppressWarnings(suppressMessages(
      efa_schmid_leiman(avg$loadings$average, Phi = avg$Phi$average))),
    "efa_schmid_leiman")
})
