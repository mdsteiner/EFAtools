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
  expect_equal(SL_lav$orig_R, NA)
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
  expect_equal(unname(SL_psych$sl[13:18, "F1"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_flex$sl[13:18, "F1"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_lav$sl[1:6, "F1"]) >= .20, rep(TRUE, 6))

  expect_equal(unname(SL_EFAtools$sl[1:12, "F1"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_SPSS$sl[1:12, "F1"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_psych$sl[1:12, "F1"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_flex$sl[1:12, "F1"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_lav$sl[7:18, "F1"]) < .20, rep(TRUE, 12))

  expect_equal(unname(SL_EFAtools$sl[7:12, "F2"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_SPSS$sl[7:12, "F2"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_psych$sl[7:12, "F2"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_flex$sl[7:12, "F2"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_lav$sl[7:12, "F2"]) >= .20, rep(TRUE, 6))

  expect_equal(unname(SL_EFAtools$sl[c(1:6, 13:18), "F2"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_SPSS$sl[c(1:6, 13:18), "F2"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_psych$sl[c(1:6, 13:18), "F2"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_flex$sl[c(1:6, 13:18), "F2"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_lav$sl[c(1:6, 13:18), "F2"]) < .20, rep(TRUE, 12))

  expect_equal(unname(SL_EFAtools$sl[1:6, "F3"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_SPSS$sl[1:6, "F3"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_psych$sl[1:6, "F3"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_flex$sl[1:6, "F3"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_lav$sl[13:18, "F3"]) >= .20, rep(TRUE, 6))

  expect_equal(unname(SL_EFAtools$sl[7:18, "F3"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_SPSS$sl[7:18, "F3"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_psych$sl[7:18, "F3"]) < .20, rep(TRUE, 12))
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

  expect_equal(unname(SL_psych$sl[, c("g", "F1", "F2", "F3")]),
               unname(sch$sl[, c("g", "F1*", "F2*", "F3*")]),
               tolerance = 1e-4)
})

test_that("settings are returned correctly", {
  skip_if_not_installed("lavaan")
  expect_named(SL_EFAtools$settings, c("estimator", "method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "input_type",
                                   "cor_method_used", "se", "b_boot", "ci", "seed", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen"))
  expect_named(SL_SPSS$settings, c("estimator", "method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "input_type",
                                   "cor_method_used", "se", "b_boot", "ci", "seed"))
  expect_named(SL_psych$settings, c("estimator", "method", "rotation", "type", "n_factors",
                                       "N", "use", "cor_method", "input_type",
                                       "cor_method_used", "se", "b_boot", "ci", "seed", "max_iter",
                                       "init_comm", "criterion", "criterion_type",
                                       "abs_eigen"))
  expect_named(SL_flex$settings, c("estimator", "method", "rotation", "type", "n_factors",
                                  "N", "use", "cor_method", "input_type",
                                  "cor_method_used", "se", "b_boot", "ci", "seed", "start_method"))
  expect_equal(SL_lav$settings, NA)

  expect_equal(SL_EFAtools$settings$estimator, "PAF")
  expect_equal(SL_SPSS$settings$estimator, "ULS")
  expect_equal(SL_psych$settings$estimator, "PAF")
  expect_equal(SL_flex$settings$estimator, "ML")

  expect_equal(SL_EFAtools$settings$rotation, "none")
  expect_equal(SL_SPSS$settings$rotation, "none")
  expect_equal(SL_psych$settings$rotation, "none")
  expect_equal(SL_flex$settings$rotation, "none")

  expect_equal(SL_EFAtools$settings$type, "EFAtools")
  expect_equal(SL_SPSS$settings$type, "SPSS")
  expect_equal(SL_psych$settings$type, "psych")
  expect_equal(SL_flex$settings$type, "EFAtools")

  expect_equal(SL_EFAtools$settings$n_factors, 1)
  expect_equal(SL_SPSS$settings$n_factors, 1)
  expect_equal(SL_psych$settings$n_factors, 1)
  expect_equal(SL_flex$settings$n_factors, 1)

  expect_equal(SL_EFAtools$settings$N, 100)
  expect_equal(SL_SPSS$settings$N, 100)
  expect_equal(SL_psych$settings$N, 100)
  expect_equal(SL_flex$settings$N, 100)

  expect_equal(SL_EFAtools$settings$use, "pairwise.complete.obs")
  expect_equal(SL_SPSS$settings$use, "pairwise.complete.obs")
  expect_equal(SL_psych$settings$use, "pairwise.complete.obs")
  expect_equal(SL_flex$settings$use, "pairwise.complete.obs")

  expect_equal(SL_EFAtools$settings$cor_method, "pearson")
  expect_equal(SL_SPSS$settings$cor_method, "pearson")
  expect_equal(SL_psych$settings$cor_method, "pearson")
  expect_equal(SL_flex$settings$cor_method, "pearson")

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
  expect_warning(efa_schmid_leiman(EFA_mod, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF", Phi = EFA_mod$Phi),
                 class = "efa_sl_phi_specified")
  expect_error(efa_schmid_leiman(EFA_mod_unrot, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF"), class = "efa_sl_not_oblique")
  expect_error(efa_schmid_leiman(EFA_mod_orth, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF"), class = "efa_sl_not_oblique")
  expect_warning(efa_schmid_leiman(fa_mod, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF", Phi = fa_mod$Phi), class = "efa_sl_phi_specified")
  expect_error(efa_schmid_leiman(fa_mod_unrot, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF"), class = "efa_sl_not_oblique")
  expect_error(efa_schmid_leiman(fa_mod_orth, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF"), class = "efa_sl_not_oblique")
  expect_error(efa_schmid_leiman(lav_fit_NA, g_name = "g"), class = "efa_sl_na_loadings")
  expect_error(efa_schmid_leiman(lav_fit_var_hey, g_name = "g"), class = "efa_sl_heywood")
  expect_error(efa_schmid_leiman(lav_fit_ho, g_name = "fu"), class = "efa_sl_g_name")
  expect_warning(efa_schmid_leiman(lav_fit_ho_inv, g_name = "g"), class = "efa_sl_second_order_loadings")
  expect_error(efa_schmid_leiman(EFA_mod$rot_loadings, estimate_control = estimate_control(type = "EFAtools"), estimator = "ML"), class = "efa_sl_phi_missing")
  expect_error(efa_schmid_leiman(lav_fit_bifactor, g_name = "g"), class = "efa_sl_not_second_order")
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

  # Two first-order factors do not identify the second-order fit either, which is a
  # warning and is raised before that fit runs.
  L2 <- matrix(0, 4, 2)
  L2[1:2, 1] <- 0.7
  L2[3:4, 2] <- 0.7
  colnames(L2) <- c("F1", "F2")
  expect_warning(
    efa_schmid_leiman(L2, Phi = matrix(c(1, .5, .5, 1), 2), estimator = "PAF"),
    class = "efa_sl_underidentified")
})

test_that("efa_fit()'s inference arguments are refused", {
  # They are efa_fit() formals, so they used to ride through the dots into the second-order
  # fit -- which runs against a placeholder N and reports none of them, so a standard error
  # was computed from a sample size that is not the data's and then thrown away, and `seed`
  # governed a fit that draws no random numbers. Only `settings` ever recorded them.
  args <- list(EFA_mod, estimate_control = estimate_control(type = "EFAtools"),
               estimator = "ML")
  expect_error(do.call(efa_schmid_leiman, c(args, list(se = "information"))),
               class = "efa_unused_dots")
  expect_error(do.call(efa_schmid_leiman, c(args, list(seed = 42))),
               class = "efa_unused_dots")
  expect_error(do.call(efa_schmid_leiman, c(args, list(b_boot = 10, ci = .90))),
               class = "efa_unused_dots")

  # the frozen wrapper forwards its dots into the same fit, so the refusal reaches it too
  expect_error(SL(EFA_mod, type = "EFAtools", method = "ML", se = "information"),
               class = "efa_unused_dots")

  # An abbreviation must be refused as well: the dots are spliced into the second-order fit
  # with do.call(), where R partial-matches them against efa_fit()'s formals, so `b_b` would
  # arrive there as `b_boot` without ever matching the refused name exactly. Only the
  # whitelist, which spells every accepted name in full, closes that route.
  expect_error(do.call(efa_schmid_leiman, c(args, list(b_b = 10))),
               class = "efa_unused_dots")
  expect_error(do.call(efa_schmid_leiman, c(args, list(see = 42))),
               class = "efa_unused_dots")
  expect_error(SL(EFA_mod, type = "EFAtools", method = "ML", b_b = 10),
               class = "efa_unused_dots")

  # a genuine efa_fit() argument still travels through the dots untouched -- `use` is an
  # efa_fit() formal and not one of efa_schmid_leiman()'s, so it really does go through `...`
  expect_no_error(suppressWarnings(suppressMessages(
    do.call(efa_schmid_leiman, c(args, list(use = "complete.obs"))))))
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

test_that("labelled factor columns are still ordered by factor number", {

  # The fallback for unlabelled columns must not disable the reordering when the
  # labels are present but out of order.
  L1 <- EFA_mod$rot_loadings[, c(2, 3, 1)]
  Phi_1 <- EFA_mod$Phi[c(2, 3, 1), c(2, 3, 1)]
  colnames(L1) <- c("F2", "F3", "F1")

  sl_scrambled <- efa_schmid_leiman(L1, Phi = Phi_1, estimator = "PAF",
                                    estimate_control = estimate_control(type = "EFAtools"))

  expect_equal(sl_scrambled$sl, SL_EFAtools$sl)

})

test_that("a named Phi is aligned to the loading factors before reordering", {
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

test_that("efa_schmid_leiman runs on an oblique fit whose column labels are permuted", {
  skip_if_not_installed("GPArotation")

  # The first-order columns are reordered by the number in their labels, which the
  # oblimin rotation returns permuted (F2, F1, F3). An unguarded reorder empties Phi,
  # and the second-order fit then aborts on a model with no variables.
  EFA_obl <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
                     estimator = "PAF", rotation = "oblimin")
  SL_obl <- efa_schmid_leiman(EFA_obl, estimator = "PAF",
                              estimate_control = estimate_control(type = "EFAtools"))

  expect_s3_class(SL_obl, "efa_schmid_leiman")
  expect_identical(class(SL_obl$sl), c("efa_sl_loadings", "SLLOADINGS"))
  expect_identical(dim(unclass(SL_obl$sl)), dim(unclass(SL_EFAtools$sl)))
  expect_true(all(is.finite(SL_obl$L2)))
  expect_length(SL_obl$L2, 3)
})

test_that("print output is stable", {
  local_reproducible_output()

  expect_snapshot(print(SL_EFAtools), transform = scrub_num)
})

rm(EFA_mod, SL_EFAtools, SL_SPSS, fa_mod, SL_psych, SL_flex, EFA_mod_unrot,
   EFA_mod_orth, fa_mod_unrot, fa_mod_orth)
if (requireNamespace("lavaan", quietly = TRUE)) {
  rm(lav_mod_ho, lav_fit_ho, SL_lav, lav_mod_NA, lav_fit_NA, lav_mod_ho_inv,
     lav_fit_ho_inv, lav_mod_var_hey, lav_fit_var_hey, lav_mod_bifactor,
     lav_fit_bifactor)
}

