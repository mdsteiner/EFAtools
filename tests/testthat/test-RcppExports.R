
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

rm(x_mat, y_mat)
