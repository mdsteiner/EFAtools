

vec_s <- c("a" = 1, "b" = 2, "c" = 4)
vec_L <- c("A" = 1, "B" = 2, "C" = 4)

mat_s <- matrix(c(0, 0, 0, 1), ncol = 2)
colnames(mat_s) <- c("a", "b")
mat_L <- matrix(c(0, 0, 0, 1), ncol = 2)
colnames(mat_L) <- c("A", "B")

int <- COMPARE(1:10, 1:10)
dec <- COMPARE(c(.1, .2), c(.1, .1))
matr <- COMPARE(matrix(c(1,1,1,2), ncol = 2), matrix(c(1,1,1,1), ncol = 2))

SPSS_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                type = "SPSS", method = "PAF", rotation = "none")
psych_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                 type = "psych", method = "PAF", rotation = "none")
load <- COMPARE(SPSS_PAF$unrot_loadings, SPSS_PAF$unrot_loadings)
load_ro1 <- COMPARE(SPSS_PAF$unrot_loadings, SPSS_PAF$unrot_loadings,
                    reorder = "names")
load_ro2 <- COMPARE(SPSS_PAF$unrot_loadings, SPSS_PAF$unrot_loadings,
                    reorder = "none")

SPSS_PAF_1 <- EFA(test_models$baseline$cormat, n_factors = 1, N = 500,
                type = "SPSS", method = "PAF", rotation = "none")
psych_PAF_1 <- EFA(test_models$baseline$cormat, n_factors = 1, N = 500,
                 type = "psych", method = "PAF", rotation = "none")
load_F1 <- COMPARE(SPSS_PAF_1$unrot_loadings, psych_PAF_1$unrot_loadings)

test_that("output class and dimensions are correct", {
  expect_is(int, "COMPARE")
  expect_is(dec, "COMPARE")
  expect_is(matr, "COMPARE")
  expect_is(load, "COMPARE")
  expect_is(load_ro1, "COMPARE")
  expect_is(load_ro2, "COMPARE")
  expect_is(load_F1, "COMPARE")

  expect_named(int, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                      "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                      "diff_corres_cross", "g", "settings"))
  expect_named(dec, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                      "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                      "diff_corres_cross", "g", "settings"))
  expect_named(matr, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                      "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                      "diff_corres_cross", "g", "settings"))
  expect_named(load, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                       "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                       "diff_corres_cross", "g", "settings"))
  expect_named(load_ro1, c("diff", "mean_abs_diff", "median_abs_diff",
                           "min_abs_diff", "max_abs_diff", "max_dec", "are_equal",
                           "diff_corres", "diff_corres_cross", "g", "settings"))
  expect_named(load_ro2, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                       "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                       "diff_corres_cross", "g", "settings"))
  expect_named(load_F1, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                       "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                       "diff_corres_cross", "g", "settings"))

  expect_is(int$diff, "integer")
  expect_is(int$mean_abs_diff, "numeric")
  expect_is(int$median_abs_diff, "numeric")
  expect_is(int$min_abs_diff, "integer")
  expect_is(int$max_abs_diff, "integer")
  expect_is(int$max_dec, "numeric")
  expect_is(int$are_equal, "numeric")
  expect_equal(int$diff_corres, NA)
  expect_equal(int$diff_corres_cross, NA)
  expect_is(int$g, "numeric")
  expect_is(int$settings, "list")

  expect_is(dec$diff, "numeric")
  expect_is(dec$mean_abs_diff, "numeric")
  expect_is(dec$median_abs_diff, "numeric")
  expect_is(dec$min_abs_diff, "numeric")
  expect_is(dec$max_abs_diff, "numeric")
  expect_is(dec$max_dec, "integer")
  expect_is(dec$are_equal, "numeric")
  expect_equal(dec$diff_corres, NA)
  expect_equal(dec$diff_corres_cross, NA)
  expect_is(dec$g, "numeric")
  expect_is(dec$settings, "list")

  expect_is(matr$diff, "matrix")
  expect_is(matr$mean_abs_diff, "numeric")
  expect_is(matr$median_abs_diff, "numeric")
  expect_is(matr$min_abs_diff, "numeric")
  expect_is(matr$max_abs_diff, "numeric")
  expect_is(matr$max_dec, "numeric")
  expect_is(matr$are_equal, "numeric")
  expect_is(matr$diff_corres, "integer")
  expect_is(matr$diff_corres_cross, "integer")
  expect_is(matr$g, "numeric")
  expect_is(matr$settings, "list")

  expect_is(load$diff, "matrix")
  expect_is(load$mean_abs_diff, "numeric")
  expect_is(load$median_abs_diff, "numeric")
  expect_is(load$min_abs_diff, "numeric")
  expect_is(load$max_abs_diff, "numeric")
  expect_is(load$max_dec, "integer")
  expect_is(load$are_equal, "numeric")
  expect_is(load$diff_corres, "integer")
  expect_is(load$diff_corres_cross, "integer")
  expect_is(load$g, "numeric")
  expect_is(load$settings, "list")

  expect_is(load_ro1$diff, "matrix")
  expect_is(load_ro1$mean_abs_diff, "numeric")
  expect_is(load_ro1$median_abs_diff, "numeric")
  expect_is(load_ro1$min_abs_diff, "numeric")
  expect_is(load_ro1$max_abs_diff, "numeric")
  expect_is(load_ro1$max_dec, "integer")
  expect_is(load_ro1$are_equal, "numeric")
  expect_is(load_ro1$diff_corres, "integer")
  expect_is(load_ro1$diff_corres_cross, "integer")
  expect_is(load_ro1$g, "numeric")
  expect_is(load_ro1$settings, "list")

  expect_is(load_ro2$diff, "matrix")
  expect_is(load_ro2$mean_abs_diff, "numeric")
  expect_is(load_ro2$median_abs_diff, "numeric")
  expect_is(load_ro2$min_abs_diff, "numeric")
  expect_is(load_ro2$max_abs_diff, "numeric")
  expect_is(load_ro2$max_dec, "integer")
  expect_is(load_ro2$are_equal, "numeric")
  expect_is(load_ro2$diff_corres, "integer")
  expect_is(load_ro2$diff_corres_cross, "integer")
  expect_is(load_ro2$g, "numeric")
  expect_is(load_ro2$settings, "list")

  expect_is(load_F1$diff, "matrix")
  expect_is(load_F1$mean_abs_diff, "numeric")
  expect_is(load_F1$median_abs_diff, "numeric")
  expect_is(load_F1$min_abs_diff, "numeric")
  expect_is(load_F1$max_abs_diff, "numeric")
  expect_is(load_F1$max_dec, "integer")
  expect_is(load_F1$are_equal, "numeric")
  expect_is(load_F1$diff_corres, "numeric")
  expect_is(load_F1$diff_corres_cross, "numeric")
  expect_is(load_F1$g, "numeric")
  expect_is(load_F1$settings, "list")

})

test_that("COMPARE returns the correct values", {
  expect_equal(int$diff, rep(0, 10))
  expect_equal(int$mean_abs_diff, 0)
  expect_equal(int$median_abs_diff, 0)
  expect_equal(int$min_abs_diff, 0)
  expect_equal(int$max_abs_diff, 0)
  expect_equal(int$max_dec, 0)
  expect_equal(int$are_equal, 0)
  expect_equal(int$g, 0)

  expect_equal(dec$diff, c(0, 0.1))
  expect_equal(dec$mean_abs_diff, 0.05)
  expect_equal(dec$median_abs_diff, 0.05)
  expect_equal(dec$min_abs_diff, 0)
  expect_equal(dec$max_abs_diff, .1)
  expect_equal(dec$max_dec, 1)
  expect_equal(dec$are_equal, 0)
  expect_equal(dec$g, 0.0707, tolerance = .01)

  expect_equal(matr$diff, matrix(c(0, 0, 0, 1), ncol = 2))
  expect_equal(matr$mean_abs_diff, 0.25)
  expect_equal(matr$median_abs_diff, 0)
  expect_equal(matr$min_abs_diff, 0)
  expect_equal(matr$max_abs_diff, 1)
  expect_equal(matr$max_dec, 0)
  expect_equal(matr$are_equal, 0)
  expect_equal(matr$g, 0.5, tolerance = .01)
  expect_equal(matr$diff_corres, 1)
  expect_equal(matr$diff_corres_cross, 0)

  expect_equal(load$diff, matrix(rep(0, 54), ncol = 3,
                                 dimnames = list(c(paste0("V", seq_len(18))),
                                                 c(paste0("F", seq_len(3))))))
  expect_equal(load$mean_abs_diff, 0)
  expect_equal(load$median_abs_diff, 0)
  expect_equal(load$min_abs_diff, 0)
  expect_equal(load$max_abs_diff, 0)
  expect_equal(load$max_dec, 17)
  expect_equal(load$are_equal, 17)
  expect_equal(load$g, 0)
  expect_equal(load$diff_corres, 0)
  expect_equal(load$diff_corres_cross, 0)

  expect_equal(load_ro1$mean_abs_diff, 0)
  expect_equal(load_ro1$median_abs_diff, 0)
  expect_equal(load_ro1$min_abs_diff, 0)
  expect_equal(load_ro1$max_abs_diff, 0)
  expect_equal(load_ro1$max_dec, 17)
  expect_equal(load_ro1$are_equal, 17)
  expect_equal(load_ro1$g, 0)
  expect_equal(load_ro1$diff_corres, 0)
  expect_equal(load_ro1$diff_corres_cross, 0)

  expect_equal(load_ro2$mean_abs_diff, 0)
  expect_equal(load_ro2$median_abs_diff, 0)
  expect_equal(load_ro2$min_abs_diff, 0)
  expect_equal(load_ro2$max_abs_diff, 0)
  expect_equal(load_ro2$max_dec, 17)
  expect_equal(load_ro2$are_equal, 17)
  expect_equal(load_ro2$g, 0)
  expect_equal(load_ro2$diff_corres, 0)
  expect_equal(load_ro2$diff_corres_cross, 0)

  expect_equal(round(load_F1$diff, 4), matrix(rep(0, 18), ncol = 1,
                                              dimnames = list(c(paste0("V",
                                                                       seq_len(18))),
                                                              "F1")))
  expect_equal(round(load_F1$mean_abs_diff, 4), 0)
  expect_equal(round(load_F1$median_abs_diff, 4), 0)
  expect_equal(round(load_F1$min_abs_diff, 4), 0)
  expect_equal(round(load_F1$max_abs_diff, 4), 0)
  expect_equal(load_F1$max_dec, 15)
  expect_equal(load_F1$are_equal, 2)
  expect_equal(round(load_F1$g, 4), 0)
  expect_equal(load_F1$diff_corres, 0)
  expect_equal(load_F1$diff_corres_cross, 0)

  expect_equal(COMPARE(vec_s, vec_s[c(3, 1, 2)],
                       reorder = "names")$mean_abs_diff, 0)
  expect_equal(COMPARE(mat_s, mat_s[, c(2, 1)],
                       reorder = "names")$mean_abs_diff, 0)
  expect_equal(COMPARE(psych_PAF$unrot_loadings,
                       psych_PAF$unrot_loadings[, c(3, 1, 2)])$mean_abs_diff, 0)
})


test_that("errors etc. are thrown correctly", {
  expect_error(COMPARE(c(1, 2), 1), " 'x' and 'y' have different lengths Compare only works with identical dimensions.\n")
  expect_error(COMPARE(c(1, 2), c("1", "2")), " 'x' is of class numeric and 'y' is of class character but must be numeric vectors or matrices\n")
  expect_error(COMPARE(c(1, 2), data.frame(x = "1", y = "2")),
               " 'x' is of class numeric and 'y' is of class data.frame but must be numeric vectors or matrices\n")

  expect_error(COMPARE(matrix(c(0, 0, 0, 1), ncol = 2),
                       matrix(c(0, 0, 0, 1), ncol = 1)), " 'x' and 'y' have different dimensions. Can only compare matrices with identical dimensions.\n")

  expect_warning(COMPARE(vec_s, vec_s), " reorder was set to 'congruence', but this only works for matrices. To reorder vectors, set reorder = 'names'. Proceeding without reordering.\n")
  expect_warning(COMPARE(vec_s, vec_L, reorder = "names"),
                 " reorder = 'names' was used but names of x and y were not identical. Results might be inaccurate.\n")
  expect_warning(COMPARE(vec_s, 1:3, reorder = "names"), " reorder was set to 'names' but at least one of 'x' and 'y' was not named. Proceeding without reordering.\n")
  expect_warning(COMPARE(1:3, 1:3, reorder = "names"), " reorder was set to 'names' but at least one of 'x' and 'y' was not named. Proceeding without reordering.\n")

  expect_warning(COMPARE(matrix(c(0, 0, 0, 1), ncol = 2),
                         matrix(c(0, 0, 0, 1), ncol = 2),
                         reorder = "names"), " reorder was set to 'names' but at least one of 'x' and 'y' was not named. Proceeding without reordering.\n")
  expect_warning(COMPARE(mat_s, mat_L, reorder = "names"),
                 " reorder = 'names' was used but colnames of x and y were not identical. Results might be inaccurate.\n")

  expect_error(COMPARE(mat_s, mat_s),
               " Tucker's congruence coefficients contained NAs, cannot reorder columns based on congruence. Try another reordering method.\n")

})

rm(int, dec, matr, SPSS_PAF, psych_PAF, load, load_ro1, load_ro2, SPSS_PAF_1,
   psych_PAF_1, load_F1, vec_s, vec_L, mat_s, mat_L)
