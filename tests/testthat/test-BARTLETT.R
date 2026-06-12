bart_cor <- BARTLETT(test_models$baseline$cormat, N = 500)
bart_raw <- BARTLETT(GRiPS_raw)

set.seed(500)
bart_rand <- BARTLETT(matrix(rnorm(100), ncol = 4))

test_that("output class and dimensions are correct", {
  expect_s3_class(bart_cor, "BARTLETT")
  expect_output(str(bart_cor), "List of 4")
  expect_s3_class(bart_raw, "BARTLETT")
  expect_output(str(bart_raw), "List of 4")
})


test_that("p-values and df are correct", {
  expect_lt(bart_cor$p_value, 0.0001)
  expect_lt(bart_raw$p_value, 0.0001)
  expect_gt(bart_rand$p_value, 0.05)

  expect_equal(bart_cor$df, 153)
  expect_equal(bart_raw$df, 28)
  expect_equal(bart_rand$df, 6)
})

test_that("settings are returned correctly", {
  expect_named(bart_cor$settings, c("N", "use", "cor_method"))
  expect_named(bart_raw$settings, c("N", "use", "cor_method"))
  expect_named(bart_rand$settings, c("N", "use", "cor_method"))

  expect_equal(bart_cor$settings$N, 500)
  expect_equal(bart_raw$settings$N, 810)
  expect_equal(bart_rand$settings$N, 25)

  expect_equal(bart_cor$settings$use, "pairwise.complete.obs")
  expect_equal(bart_raw$settings$use, "pairwise.complete.obs")
  expect_equal(bart_rand$settings$use, "pairwise.complete.obs")

  expect_equal(bart_cor$settings$cor_method, "pearson")
  expect_equal(bart_raw$settings$cor_method, "pearson")
  expect_equal(bart_rand$settings$cor_method, "pearson")
})

# Create singular correlation matrix for tests
x <- rnorm(10)
y <- rnorm(10)
z <- x + y

dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors are thrown correctly", {
  expect_error(BARTLETT(1:5), class = "efa_input_not_matrix")
  expect_error(BARTLETT(test_models$baseline$cormat), class = "efa_n_required")
  expect_message(BARTLETT(GRiPS_raw), class = "efa_cor_from_data")
  expect_warning(BARTLETT(GRiPS_raw, N = 20), class = "efa_n_from_data")
  expect_error(BARTLETT(dat_sing), class = "efa_cor_singular")
  expect_error(BARTLETT(cor_sing, N = 10), class = "efa_cor_singular")
  expect_warning(BARTLETT(cor_nposdef, N = 10), class = "efa_cor_smoothed")
})

test_that("print output is stable", {
  local_reproducible_output()

  # significant
  expect_snapshot(print(bart_cor), transform = scrub_num)

  # not significant
  expect_snapshot(print(bart_rand), transform = scrub_num)

  # test did not render a result
  bart_na <- structure(list(chisq = NA_real_, df = NA_integer_, p_value = NA_real_,
                            settings = list()), class = "BARTLETT")
  expect_snapshot(print(bart_na), transform = scrub_num)
})

rm(bart_cor, bart_raw, bart_rand, x, y, z, dat_sing, cor_sing, cor_nposdef)


