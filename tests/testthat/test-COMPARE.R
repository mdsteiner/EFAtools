
int <- COMPARE(1:10, 1:10)
dec <- COMPARE(c(.1, .2), c(.1, .1))
matr <- COMPARE(matrix(c(1,1,1,2), ncol = 2), matrix(c(1,1,1,1), ncol = 2))

test_that("output class and dimensions are correct", {
  expect_is(int, "COMPARE")
  expect_is(dec, "COMPARE")
  expect_is(matr, "COMPARE")

  expect_named(int, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                      "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                      "diff_corres_cross", "g", "settings"))
  expect_named(dec, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                      "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                      "diff_corres_cross", "g", "settings"))
  expect_named(matr, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
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
})


test_that("errors etc. are thrown correctly", {
  expect_error(COMPARE(c(1, 2), 1), " 'x' and 'y' have different lengths Compare only works with identical dimensions.\n")
  expect_error(COMPARE(c(1, 2), c("1", "2")), " 'x' is of class numeric and 'y' is of class character but must be numeric vectors or matrices\n")
  expect_error(COMPARE(c(1, 2), data.frame(x = "1", y = "2")),
               " 'x' is of class numeric and 'y' is of class data.frame but must be numeric vectors or matrices\n")

  expect_error(COMPARE(matrix(c(0, 0, 0, 1), ncol = 2),
                       matrix(c(0, 0, 0, 1), ncol = 1)), " 'x' and 'y' have different dimensions. Can only compare matrices with identical dimensions.\n")

})

rm(int, dec, matr)
