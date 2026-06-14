## Use with an output from the EFAtools::EFA function, both with type EFAtools
EFA_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
SL_EFAtools <- SL(EFA_mod, type = "EFAtools", method = "PAF")

# with type SPSS and method ULS
SL_SPSS <- SL(EFA_mod, type = "SPSS", method = "ULS")

## Use with an output from the psych::fa function with type psych in SL
fa_mod <- psych::fa(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                    fm = "pa", rotate = "Promax", n.rotations = 1)
SL_psych <- SL(fa_mod, type = "psych", method = "PAF")

## Use more flexibly by entering a pattern matrix and phi directly, with method
## ML
SL_flex <- SL(EFA_mod$rot_loadings, Phi = EFA_mod$Phi, type = "EFAtools",
              method = "ML")

## Use with a second-order lavaan solution
if (requireNamespace("lavaan", quietly = TRUE)) {
lav_mod_ho <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
               F2 =~ V7 + V8 + V9 + V10 + V11 + V12
               F3 =~ V13 + V14 + V15 + V16 + V17 + V18
               g =~ F1 + F2 + F3'
lav_fit_ho <- suppressWarnings(lavaan::cfa(lav_mod_ho,
                                           sample.cov = test_models$baseline$cormat,
                                           sample.nobs = 500, estimator = "ml"))
SL_lav <- SL(lav_fit_ho, g_name = "g")
}

test_that("output class and dimensions are correct", {
  skip_if_not_installed("lavaan")
  expect_s3_class(SL_EFAtools, "SL")
  expect_s3_class(SL_SPSS, "SL")
  expect_s3_class(SL_psych, "SL")
  expect_s3_class(SL_flex, "SL")
  expect_s3_class(SL_lav, "SL")

  expect_output(str(SL_EFAtools), "List of 6")
  expect_output(str(SL_SPSS), "List of 6")
  expect_output(str(SL_psych), "List of 6")
  expect_output(str(SL_flex), "List of 6")
  expect_output(str(SL_lav), "List of 6")
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

test_that("settings are returned correctly", {
  skip_if_not_installed("lavaan")
  expect_named(SL_EFAtools$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "se", "b_boot", "ci", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen"))
  expect_named(SL_SPSS$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "se", "b_boot", "ci"))
  expect_named(SL_psych$settings, c("method", "rotation", "type", "n_factors",
                                       "N", "use", "cor_method", "se", "b_boot", "ci", "max_iter",
                                       "init_comm", "criterion", "criterion_type",
                                       "abs_eigen"))
  expect_named(SL_flex$settings, c("method", "rotation", "type", "n_factors",
                                  "N", "use", "cor_method", "se", "b_boot", "ci", "start_method"))
  expect_equal(SL_lav$settings, NA)

  expect_equal(SL_EFAtools$settings$method, "PAF")
  expect_equal(SL_SPSS$settings$method, "ULS")
  expect_equal(SL_psych$settings$method, "PAF")
  expect_equal(SL_flex$settings$method, "ML")

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
}

test_that("errors are thrown correctly", {
  skip_if_not_installed("lavaan")
  expect_error(SL(1:5), class = "efa_sl_bad_input")
  expect_warning(SL(EFA_mod, type = "EFAtools", method = "PAF", Phi = EFA_mod$Phi),
                 class = "efa_sl_phi_specified")
  expect_error(SL(EFA_mod_unrot, type = "EFAtools", method = "PAF"), class = "efa_sl_not_oblique")
  expect_error(SL(EFA_mod_orth, type = "EFAtools", method = "PAF"), class = "efa_sl_not_oblique")
  expect_warning(SL(fa_mod, type = "EFAtools", method = "PAF", Phi = fa_mod$Phi), class = "efa_sl_phi_specified")
  expect_error(SL(fa_mod_unrot, type = "EFAtools", method = "PAF"), class = "efa_sl_not_oblique")
  expect_error(SL(fa_mod_orth, type = "EFAtools", method = "PAF"), class = "efa_sl_not_oblique")
  expect_error(SL(lav_fit_NA, g_name = "g"), class = "efa_omega_na_loadings")
  expect_error(SL(lav_fit_var_hey, g_name = "g"), class = "efa_omega_heywood")
  expect_error(SL(lav_fit_ho, g_name = "fu"), class = "efa_omega_g_name")
  expect_warning(SL(lav_fit_ho_inv, g_name = "g"), class = "efa_sl_second_order_loadings")
  expect_error(SL(EFA_mod$rot_loadings, type = "EFAtools", method = "ML"), class = "efa_sl_phi_missing")
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

  expect_error(SL(L1, Phi = Phi_heywood, type = "EFAtools", method = "PAF"),
               class = "efa_omega_heywood")
})

test_that("print output is stable", {
  local_reproducible_output()

  expect_snapshot(print(SL_EFAtools), transform = scrub_num)
})

rm(EFA_mod, SL_EFAtools, SL_SPSS, fa_mod, SL_psych, SL_flex, EFA_mod_unrot,
   EFA_mod_orth, fa_mod_unrot, fa_mod_orth)
if (requireNamespace("lavaan", quietly = TRUE)) {
  rm(lav_mod_ho, lav_fit_ho, SL_lav, lav_mod_NA, lav_fit_NA, lav_mod_ho_inv,
     lav_fit_ho_inv, lav_mod_var_hey, lav_fit_var_hey)
}

