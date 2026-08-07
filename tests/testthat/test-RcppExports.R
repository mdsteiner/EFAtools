
x_mat <- matrix(c(.1, .1, .3,
                  .32, .4, .1,
                   0,  0,  0), ncol = 3, byrow = TRUE)
y_mat <- matrix(c(.1, .1, .3,
                  .4, .29, .1,
                   0,  0,  0), ncol = 3, byrow = TRUE)
test_that(".factor_corres works", {
  expect_equal(.factor_corres(x_mat, y_mat)$diff_corres, 1)
  expect_equal(.factor_corres(y_mat, x_mat)$diff_corres_cross, 1)

  expect_equal(.factor_corres(x_mat,
                              matrix(0, ncol = 3, nrow = 3))$diff_corres, 2)
  expect_equal(.factor_corres(x_mat,
                              matrix(0, ncol = 3, nrow = 3))$diff_corres_cross, 2)
  expect_equal(.factor_corres(matrix(0, ncol = 3, nrow = 3),
                              x_mat)$diff_corres, 2)
  expect_equal(.factor_corres(matrix(0, ncol = 3, nrow = 3),
                              x_mat)$diff_corres_cross, 2)

  expect_equal(.factor_corres(x_mat, x_mat)$diff_corres, 0)
  expect_equal(.factor_corres(x_mat, x_mat)$diff_corres_cross, 0)
})

test_that(".factor_corres cross-loading sets do not collide for >= 10 factors", {
  # One indicator loads on factors {1, 2}; the other on factor {12}. Without a
  # separator both collapse to "12" and the differing correspondence is missed.
  x12 <- matrix(0, nrow = 1, ncol = 12)
  x12[1, 1] <- x12[1, 2] <- 0.5
  y12 <- matrix(0, nrow = 1, ncol = 12)
  y12[1, 12] <- 0.5
  expect_equal(.factor_corres(x12, y12)$diff_corres_cross, 1)
})

test_that(".factor_corres ignores missing cells rather than the whole row", {
  x <- rbind(c(.8, NA), c(.1, .7))
  y <- rbind(c(.8, .1), c(.1, .7))
  out <- .factor_corres(x, y)
  expect_identical(out$x_corres, c(1L, 2L))
  expect_identical(out$diff_corres, 0L)
  expect_identical(out$diff_corres_cross, 0L)

  # A factor position missing on either side is excluded from both classifications.
  pairwise <- .factor_corres(matrix(c(.8, NA), nrow = 1),
                             matrix(c(.4, .9), nrow = 1))
  expect_identical(pairwise$x_corres, 1L)
  expect_identical(pairwise$y_corres, 1L)
  expect_identical(pairwise$diff_corres, 0L)

  # Rows without any jointly observed loading are unclassifiable and do not count
  # as a difference; if there are no comparable rows, the count is undefined.
  partly_empty <- .factor_corres(rbind(c(NA, NA), c(.7, .1)),
                                  rbind(c(.2, .8), c(.6, .2)))
  expect_true(is.na(partly_empty$x_corres[1]))
  expect_identical(partly_empty$diff_corres, 0L)
  fully_empty <- .factor_corres(matrix(c(NA, NA), nrow = 1),
                                 matrix(c(.2, .8), nrow = 1))
  expect_true(is.na(fully_empty$diff_corres))
})

test_that("native simulation and polychoric entry points reject malformed inputs", {
  expect_error(.parallel_sim(0, 3, 20, 1, 10), "n_datasets")
  expect_error(.parallel_sim(1, 3, 20, 9, 10), "eigen_type")
  expect_error(.parallel_sim(1, 3, 3, 2, 10), "larger than n_vars")

  bad_negative <- cbind(c(0L, 1L, -1L), c(0L, 1L, 0L))
  expect_error(.polychoric_cpp(bad_negative, "none", FALSE), "non-negative")
  bad_gap <- cbind(c(0L, 2L, 2L), c(0L, 1L, 0L))
  expect_error(.polychoric_cpp(bad_gap, "none", FALSE), "consecutive")
})

rm(x_mat, y_mat)
