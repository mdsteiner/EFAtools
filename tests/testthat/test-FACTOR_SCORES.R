EFA_raw <- suppressMessages(EFA(DOSPERT_raw, n_factors = 10, type = "EFAtools",
                                method = "PAF", rotation = "oblimin"))
fac_scores_raw <- FACTOR_SCORES(DOSPERT_raw, f = EFA_raw, method = "Bartlett")

EFA_cor <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
               type = "EFAtools", method = "PAF", rotation = "oblimin")
fac_scores_cor <- suppressMessages(FACTOR_SCORES(test_models$baseline$cormat,
                                                 f = EFA_cor))

EFA_unrot <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                 type = "EFAtools", method = "PAF", rotation = "none")
fac_scores_unrot <- suppressMessages(FACTOR_SCORES(test_models$baseline$cormat,
                                                   f = EFA_unrot))

fac_scores_man <- suppressMessages(FACTOR_SCORES(test_models$baseline$cormat,
                                                 f = unclass(EFA_unrot$unrot_loadings)))


test_that("output is correct", {
  expect_s3_class(fac_scores_raw, "FACTOR_SCORES")
  expect_s3_class(fac_scores_cor, "FACTOR_SCORES")
  expect_s3_class(fac_scores_unrot, "FACTOR_SCORES")
  expect_s3_class(fac_scores_man, "FACTOR_SCORES")

  expect_named(fac_scores_raw, c("scores", "weights", "r.scores", "missing",
                                 "R2", "settings"))
  expect_named(fac_scores_cor, c("scores", "weights", "r.scores",
                                 "R2", "settings"))
  expect_named(fac_scores_unrot, c("scores", "weights", "r.scores", "R2",
                                   "settings"))
  expect_named(fac_scores_man, c("scores", "weights", "r.scores",
                                 "R2", "settings"))
})

test_that("factor scores have the expected shape and recover a known factor", {
  skip_on_cran()

  # raw-data scores: one row per observation, one column per factor (guards against
  # a transposed or weight-matrix mix-up in the returned object)
  expect_equal(dim(fac_scores_raw$scores), c(nrow(DOSPERT_raw), 10L))
  expect_false(anyNA(fac_scores_raw$scores))

  # On a strongly unidimensional scale the single estimated factor score must track
  # the factor it represents; the mean item response is a faithful proxy for that
  # latent factor, so the two correlate almost perfectly.
  efa_1f <- suppressMessages(EFA(GRiPS_raw, n_factors = 1, type = "EFAtools",
                                 method = "PAF"))
  fs_1f <- suppressMessages(FACTOR_SCORES(GRiPS_raw, f = efa_1f, method = "Bartlett"))
  expect_equal(dim(fs_1f$scores), c(nrow(GRiPS_raw), 1L))
  expect_gt(abs(stats::cor(fs_1f$scores[, 1], rowMeans(GRiPS_raw))), 0.95)
})

test_that("settings are returned correctly", {
  expect_named(fac_scores_raw$settings, "method")
  expect_named(fac_scores_cor$settings, "method")
  expect_named(fac_scores_unrot$settings, "method")
  expect_named(fac_scores_man$settings, "method")

  expect_equal(fac_scores_raw$settings$method, "Bartlett")
  expect_equal(fac_scores_cor$settings$method, "Thurstone")
  expect_equal(fac_scores_unrot$settings$method, "Thurstone")
  expect_equal(fac_scores_man$settings$method, "Thurstone")
})


test_that("warnings and errors are thrown correctly", {
  expect_error(FACTOR_SCORES(1:5), class = "efa_input_not_matrix")
  expect_error(FACTOR_SCORES(DOSPERT_raw, f = 1:5), class = "efa_scores_bad_f")
  expect_message(FACTOR_SCORES(test_models$baseline$cormat, f = EFA_cor), class = "efa_scores_needs_raw")
  expect_message(FACTOR_SCORES(test_models$baseline$cormat,
                               f = unclass(EFA_unrot$unrot_loadings)), class = "efa_scores_phi_null")

})


rm(EFA_raw, fac_scores_raw, EFA_cor, fac_scores_cor, EFA_unrot, fac_scores_unrot,
   fac_scores_man)
