
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

rm(x_mat, y_mat)
