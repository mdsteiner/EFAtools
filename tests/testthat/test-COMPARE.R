

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
  expect_s3_class(int, "COMPARE")
  expect_s3_class(dec, "COMPARE")
  expect_s3_class(matr, "COMPARE")
  expect_s3_class(load, "COMPARE")
  expect_s3_class(load_ro1, "COMPARE")
  expect_s3_class(load_ro2, "COMPARE")
  expect_s3_class(load_F1, "COMPARE")

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

  expect_type(int$diff, "integer")
  expect_type(int$mean_abs_diff, "double")
  expect_type(int$median_abs_diff, "double")
  expect_type(int$min_abs_diff, "integer")
  expect_type(int$max_abs_diff, "integer")
  expect_type(int$max_dec, "double")
  expect_type(int$are_equal, "double")
  expect_equal(int$diff_corres, NA)
  expect_equal(int$diff_corres_cross, NA)
  expect_type(int$g, "double")
  expect_type(int$settings, "list")

  expect_type(dec$diff, "double")
  expect_type(dec$mean_abs_diff, "double")
  expect_type(dec$median_abs_diff, "double")
  expect_type(dec$min_abs_diff, "double")
  expect_type(dec$max_abs_diff, "double")
  expect_type(dec$max_dec, "integer")
  expect_type(dec$are_equal, "double")
  expect_equal(dec$diff_corres, NA)
  expect_equal(dec$diff_corres_cross, NA)
  expect_type(dec$g, "double")
  expect_type(dec$settings, "list")

  checkmate::expect_matrix(matr$diff)
  expect_type(matr$mean_abs_diff, "double")
  expect_type(matr$median_abs_diff, "double")
  expect_type(matr$min_abs_diff, "double")
  expect_type(matr$max_abs_diff, "double")
  expect_type(matr$max_dec, "double")
  expect_type(matr$are_equal, "double")
  expect_type(matr$diff_corres, "integer")
  expect_type(matr$diff_corres_cross, "integer")
  expect_type(matr$g, "double")
  expect_type(matr$settings, "list")

  checkmate::expect_matrix(load$diff)
  expect_type(load$mean_abs_diff, "double")
  expect_type(load$median_abs_diff, "double")
  expect_type(load$min_abs_diff, "double")
  expect_type(load$max_abs_diff, "double")
  expect_type(load$max_dec, "integer")
  expect_type(load$are_equal, "double")
  expect_type(load$diff_corres, "integer")
  expect_type(load$diff_corres_cross, "integer")
  expect_type(load$g, "double")
  expect_type(load$settings, "list")

  checkmate::expect_matrix(load_ro1$diff)
  expect_type(load_ro1$mean_abs_diff, "double")
  expect_type(load_ro1$median_abs_diff, "double")
  expect_type(load_ro1$min_abs_diff, "double")
  expect_type(load_ro1$max_abs_diff, "double")
  expect_type(load_ro1$max_dec, "integer")
  expect_type(load_ro1$are_equal, "double")
  expect_type(load_ro1$diff_corres, "integer")
  expect_type(load_ro1$diff_corres_cross, "integer")
  expect_type(load_ro1$g, "double")
  expect_type(load_ro1$settings, "list")

  checkmate::expect_matrix(load_ro2$diff)
  expect_type(load_ro2$mean_abs_diff, "double")
  expect_type(load_ro2$median_abs_diff, "double")
  expect_type(load_ro2$min_abs_diff, "double")
  expect_type(load_ro2$max_abs_diff, "double")
  expect_type(load_ro2$max_dec, "integer")
  expect_type(load_ro2$are_equal, "double")
  expect_type(load_ro2$diff_corres, "integer")
  expect_type(load_ro2$diff_corres_cross, "integer")
  expect_type(load_ro2$g, "double")
  expect_type(load_ro2$settings, "list")

  checkmate::expect_matrix(load_F1$diff)
  expect_type(load_F1$mean_abs_diff, "double")
  expect_type(load_F1$median_abs_diff, "double")
  expect_type(load_F1$min_abs_diff, "double")
  expect_type(load_F1$max_abs_diff, "double")
  expect_type(load_F1$max_dec, "integer")
  expect_type(load_F1$are_equal, "double")
  expect_type(load_F1$diff_corres, "double")
  expect_type(load_F1$diff_corres_cross, "double")
  expect_type(load_F1$g, "double")
  expect_type(load_F1$settings, "list")

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
  expect_equal(load$max_dec, 17, tolerance = 1)
  expect_equal(load$are_equal, 17, tolerance = 1)
  expect_equal(load$g, 0)
  expect_equal(load$diff_corres, 0)
  expect_equal(load$diff_corres_cross, 0)

  expect_equal(load_ro1$mean_abs_diff, 0)
  expect_equal(load_ro1$median_abs_diff, 0)
  expect_equal(load_ro1$min_abs_diff, 0)
  expect_equal(load_ro1$max_abs_diff, 0)
  expect_equal(load_ro1$max_dec, 17, tolerance = 1)
  expect_equal(load_ro1$are_equal, 17, tolerance = 1)
  expect_equal(load_ro1$g, 0)
  expect_equal(load_ro1$diff_corres, 0)
  expect_equal(load_ro1$diff_corres_cross, 0)

  expect_equal(load_ro2$mean_abs_diff, 0)
  expect_equal(load_ro2$median_abs_diff, 0)
  expect_equal(load_ro2$min_abs_diff, 0)
  expect_equal(load_ro2$max_abs_diff, 0)
  expect_equal(load_ro2$max_dec, 17, tolerance = 1)
  expect_equal(load_ro2$are_equal, 17, tolerance = 1)
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


test_that("congruence reordering yields a true permutation when columns collide", {
  # x has orthonormal factors; y is built so that the greedy row-wise which.max
  # match assigns y-column 1 to BOTH x-factor 1 and x-factor 2 (a collision that
  # would duplicate one y column and drop another). The optimal one-to-one
  # assignment is the permutation c(1, 3, 2).
  x <- diag(3)
  y <- matrix(c(0.7, 0.7, 0.10,    # aligned with factors 1 and 2
                0.2, 0.1, 0.97,    # aligned with factor 3
                0.1, 0.2, 0.95),   # aligned with factor 3
              ncol = 3)

  cmp <- COMPARE(x, y, reorder = "congruence", corres = FALSE)

  # y as reordered (and sign-matched) inside COMPARE
  y_reordered <- x - cmp$diff

  # index in y that each reordered column came from
  recovered <- vapply(seq_len(ncol(y_reordered)), function(j) {
    which(vapply(seq_len(ncol(y)), function(k)
      isTRUE(all.equal(y_reordered[, j], y[, k])), logical(1)))
  }, integer(1))

  # a genuine permutation: every factor used exactly once, never duplicated or
  # dropped (the greedy which.max produced c(1, 1, 2))
  expect_setequal(recovered, seq_len(ncol(y)))
  expect_equal(recovered, c(1L, 3L, 2L))
})


test_that("errors etc. are thrown correctly", {
  expect_error(COMPARE(c(1, 2), 1), class = "efa_compare_dim_mismatch")
  expect_error(COMPARE(c(1, 2), c("1", "2")), class = "efa_compare_bad_input")
  expect_error(COMPARE(c(1, 2), data.frame(x = "1", y = "2")), class = "efa_compare_bad_input")

  expect_error(COMPARE(matrix(c(0, 0, 0, 1), ncol = 2),
                       matrix(c(0, 0, 0, 1), ncol = 1)), class = "efa_compare_dim_mismatch")

  expect_warning(COMPARE(vec_s, vec_s), class = "efa_compare_reorder_vectors")
  expect_warning(COMPARE(vec_s, vec_L, reorder = "names"), class = "efa_compare_reorder_mismatch")
  expect_warning(COMPARE(vec_s, 1:3, reorder = "names"), class = "efa_compare_reorder_unnamed")
  expect_warning(COMPARE(1:3, 1:3, reorder = "names"), class = "efa_compare_reorder_unnamed")

  expect_warning(COMPARE(matrix(c(0, 0, 0, 1), ncol = 2),
                         matrix(c(0, 0, 0, 1), ncol = 2),
                         reorder = "names"), class = "efa_compare_reorder_unnamed")
  expect_warning(COMPARE(mat_s, mat_L, reorder = "names"), class = "efa_compare_reorder_mismatch")

  expect_error(COMPARE(mat_s, mat_s), class = "efa_compare_congruence_na")

})

test_that("print output is stable", {
  local_reproducible_output()

  # matrix difference (variables x factors). matr/int hold exact 0/1 differences, so the
  # snapshots are recorded literally (no scrub) to pin the column alignment.
  expect_snapshot(print(matr))

  # vector difference (rendered as a single column, one value per row)
  expect_snapshot(print(int))
})

test_that("plot returns a ggplot and guards too-few differences", {
  # Smoke-test only (no vdiffr baseline) because geom_jitter draws random positions,
  # so the rendered plot is not reproducible across runs.
  p <- plot(matr)
  expect_s3_class(p, "ggplot")

  # too few differences to plot
  expect_error(plot(dec), class = "efa_compare_too_few_to_plot")
})

rm(int, dec, matr, SPSS_PAF, psych_PAF, load, load_ro1, load_ro2, SPSS_PAF_1,
   psych_PAF_1, load_F1, vec_s, vec_L, mat_s, mat_L)
