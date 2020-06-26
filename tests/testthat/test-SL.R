## Use with an output from the EFAtools::EFA function, both with type EFAtools
EFA_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
SL_EFAtools <- SL(EFA_mod, type = "EFAtools", method = "PAF")

# with type SPSS and method ULS
SL_SPSS <- SL(EFA_mod, type = "SPSS", method = "ULS")

## Use with an output from the psych::fa function with type psych in SL
fa_mod <- psych::fa(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                    fm = "pa", rotate = "Promax")
SL_psych <- SL(fa_mod, type = "psych", method = "PAF")

## Use more flexibly by entering a pattern matrix and phi directly, with method
## ML
SL_flex <- SL(EFA_mod$rot_loadings, Phi = EFA_mod$Phi, type = "EFAtools",
              method = "ML")

test_that("output class and dimensions are correct", {
  expect_is(SL_EFAtools, "SL")
  expect_is(SL_SPSS, "SL")
  expect_is(SL_psych, "SL")
  expect_is(SL_flex, "SL")

  expect_output(str(SL_EFAtools), "List of 6")
  expect_output(str(SL_SPSS), "List of 6")
  expect_output(str(SL_psych), "List of 6")
  expect_output(str(SL_flex), "List of 6")
})

test_that("original correlation is correct", {
  expect_equal(SL_EFAtools$orig_R, test_models$baseline$cormat)
  expect_equal(SL_SPSS$orig_R, test_models$baseline$cormat)
  expect_equal(SL_psych$orig_R, test_models$baseline$cormat)
  expect_null(SL_flex$orig_R)
})

test_that("sl solution is correct", {
  expect_equal(unname(SL_EFAtools$sl[, "h2"]) + unname(SL_EFAtools$sl[, "u2"]),
               rep(1, 18))
  expect_equal(unname(SL_SPSS$sl[, "h2"]) + unname(SL_SPSS$sl[, "u2"]),
               rep(1, 18))
  expect_equal(unname(SL_psych$sl[, "h2"]) + unname(SL_psych$sl[, "u2"]),
               rep(1, 18))
  expect_equal(unname(SL_flex$sl[, "h2"]) + unname(SL_flex$sl[, "u2"]),
               rep(1, 18))

  expect_equal(unname(SL_EFAtools$sl[, "g"]) >= .20, rep(TRUE, 18))
  expect_equal(unname(SL_SPSS$sl[, "g"]) >= .20, rep(TRUE, 18))
  expect_equal(unname(SL_psych$sl[, "g"]) >= .20, rep(TRUE, 18))
  expect_equal(unname(SL_flex$sl[, "g"]) >= .20, rep(TRUE, 18))

  expect_equal(unname(SL_EFAtools$sl[13:18, "F1"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_SPSS$sl[13:18, "F1"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_psych$sl[13:18, "F1"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_flex$sl[13:18, "F1"]) >= .20, rep(TRUE, 6))

  expect_equal(unname(SL_EFAtools$sl[1:12, "F1"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_SPSS$sl[1:12, "F1"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_psych$sl[1:12, "F1"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_flex$sl[1:12, "F1"]) < .20, rep(TRUE, 12))

  expect_equal(unname(SL_EFAtools$sl[7:12, "F2"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_SPSS$sl[7:12, "F2"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_psych$sl[7:12, "F2"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_flex$sl[7:12, "F2"]) >= .20, rep(TRUE, 6))

  expect_equal(unname(SL_EFAtools$sl[c(1:6, 13:18), "F2"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_SPSS$sl[c(1:6, 13:18), "F2"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_psych$sl[c(1:6, 13:18), "F2"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_flex$sl[c(1:6, 13:18), "F2"]) < .20, rep(TRUE, 12))

  expect_equal(unname(SL_EFAtools$sl[1:6, "F3"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_SPSS$sl[1:6, "F3"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_psych$sl[1:6, "F3"]) >= .20, rep(TRUE, 6))
  expect_equal(unname(SL_flex$sl[1:6, "F3"]) >= .20, rep(TRUE, 6))

  expect_equal(unname(SL_EFAtools$sl[7:18, "F3"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_SPSS$sl[7:18, "F3"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_psych$sl[7:18, "F3"]) < .20, rep(TRUE, 12))
  expect_equal(unname(SL_flex$sl[7:18, "F3"]) < .20, rep(TRUE, 12))
})

test_that("settings are returned correctly", {
  expect_named(SL_EFAtools$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method", "max_iter",
                                   "init_comm", "criterion", "criterion_type",
                                   "abs_eigen"))
  expect_named(SL_SPSS$settings, c("method", "rotation", "type", "n_factors",
                                   "N", "use", "cor_method"))
  expect_named(SL_psych$settings, c("method", "rotation", "type", "n_factors",
                                       "N", "use", "cor_method", "max_iter",
                                       "init_comm", "criterion", "criterion_type",
                                       "abs_eigen"))
  expect_named(SL_flex$settings, c("method", "rotation", "type", "n_factors",
                                  "N", "use", "cor_method", "start_method"))

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

  expect_equal(SL_EFAtools$settings$init_comm, "mac")
  expect_equal(SL_psych$settings$init_comm, "smc")

  expect_equal(SL_EFAtools$settings$criterion, 0.001)
  expect_equal(SL_psych$settings$criterion,  0.001)

  expect_equal(SL_EFAtools$settings$criterion_type, "sums")
  expect_equal(SL_psych$settings$criterion_type, "sums")

  expect_equal(SL_EFAtools$settings$abs_eigen, TRUE)
  expect_equal(SL_psych$settings$abs_eigen, FALSE)

  expect_equal(SL_flex$settings$start_method, "factanal")
})

# Prepare data for input
EFA_mod_unrot <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
                     type = "EFAtools", method = "PAF", rotation = "none")
EFA_mod_orth <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
                     type = "EFAtools", method = "PAF", rotation = "varimax")
fa_mod_unrot <- psych::fa(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                          fm = "pa", rotate = "none")
fa_mod_orth <- psych::fa(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
                         fm = "pa", rotate = "varimax")

test_that("errors are thrown correctly", {
  expect_error(SL(1:5), " 'x' is neither an object of class EFA or fa nor a matrix, nor of class LOADINGS or loadings.\n")
  expect_warning(SL(EFA_mod, type = "EFAtools", method = "PAF", Phi = EFA_mod$Phi),
                 " Phi argument is specified. Specified factor intercorrelations are taken. To take factor intercorrelations from the EFA output, leave Phi = NULL\n")
  expect_error(SL(EFA_mod_unrot, type = "EFAtools", method = "PAF"), " 'x' is either a non-rotated or orthogonal factor solution. SL needs an oblique factor solution\n")
  expect_error(SL(EFA_mod_orth, type = "EFAtools", method = "PAF"), " 'x' is either a non-rotated or orthogonal factor solution. SL needs an oblique factor solution\n")
  expect_warning(SL(fa_mod, type = "EFAtools", method = "PAF", Phi = fa_mod$Phi), " Phi argument is specified. Specified factor intercorrelations are taken. To take factor intercorrelations from the psych fa output, leave Phi = NULL\n")
  expect_error(SL(fa_mod_unrot, type = "EFAtools", method = "PAF"), " 'x' is either a non-rotated or orthogonal factor solution. SL needs an oblique factor solution\n")
  expect_error(SL(fa_mod_orth, type = "EFAtools", method = "PAF"), " 'x' is either a non-rotated or orthogonal factor solution. SL needs an oblique factor solution\n")
  expect_error(SL(EFA_mod$rot_loadings, type = "EFAtools", method = "ML"), " Phi not provided. Either enter an oblique factor solution from EFAtools::EFA or from psych::fa, or provide Phi\n")
})

rm(EFA_mod, SL_EFAtools, SL_SPSS, fa_mod, SL_psych, SL_flex, EFA_mod_unrot,
   EFA_mod_orth, fa_mod_unrot, fa_mod_orth)

