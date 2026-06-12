kmo_cor <- KMO(test_models$baseline$cormat)
kmo_raw <- KMO(GRiPS_raw)
dat_nonames <- test_models$baseline$cormat
colnames(dat_nonames) <- NULL
kmo_nona <- KMO(dat_nonames)

test_that("output class and dimensions are correct", {
  expect_s3_class(kmo_cor, "KMO")
  expect_output(str(kmo_cor), "List of 3")
  expect_s3_class(kmo_raw, "KMO")
  expect_output(str(kmo_raw), "List of 3")
  expect_s3_class(kmo_nona, "KMO")
  expect_output(str(kmo_nona), "List of 3")
})

test_that("KMO values are correct", {
  expect_equal(kmo_cor$KMO, 0.916, tolerance = 1e-3)
  expect_equal(kmo_raw$KMO, 0.955, tolerance = 1e-3)
  expect_equal(kmo_nona$KMO, 0.916, tolerance = 1e-3)

  expect_length(kmo_cor$KMO_i, ncol(test_models$baseline$cormat))
  expect_length(kmo_raw$KMO_i, ncol(GRiPS_raw))
  expect_length(kmo_nona$KMO_i, ncol(dat_nonames))
})

test_that("settings are returned correctly", {
  expect_named(kmo_cor$settings, c("use", "cor_method"))
  expect_named(kmo_raw$settings, c("use", "cor_method"))
  expect_named(kmo_nona$settings, c("use", "cor_method"))

  expect_equal(kmo_cor$settings$use, "pairwise.complete.obs")
  expect_equal(kmo_raw$settings$use, "pairwise.complete.obs")
  expect_equal(kmo_nona$settings$use, "pairwise.complete.obs")

  expect_equal(kmo_cor$settings$cor_method, "pearson")
  expect_equal(kmo_raw$settings$cor_method, "pearson")
  expect_equal(kmo_nona$settings$cor_method, "pearson")
})

# Create singular correlation matrix for tests
set.seed(42)
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)

cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors are thrown correctly", {
  expect_error(KMO(1:5), class = "efa_input_not_matrix")
  expect_message(KMO(GRiPS_raw), class = "efa_cor_from_data")
  expect_error(KMO(dat_sing), class = "efa_cor_singular")
  expect_error(KMO(cor_sing), class = "efa_cor_singular")
  expect_warning(KMO(cor_nposdef), class = "efa_cor_smoothed")
})

test_that("print output is stable", {
  local_reproducible_output()

  # high KMO: tick, verdict, and the per-variable values
  expect_snapshot(print(kmo_cor), transform = scrub_num)

  # low KMO: cross and the "not suitable" verdict
  kmo_low <- structure(list(KMO = 0.45, KMO_i = c(V1 = 0.40, V2 = 0.50, V3 = 0.45)),
                       class = "KMO")
  expect_snapshot(print(kmo_low), transform = scrub_num)

  # KMO value not available
  kmo_na <- structure(list(KMO = NA_real_, KMO_i = NA_real_), class = "KMO")
  expect_snapshot(print(kmo_na), transform = scrub_num)
})

rm(kmo_cor, kmo_raw, dat_nonames, kmo_nona, x, y, z, dat_sing, cor_sing, cor_nposdef)
