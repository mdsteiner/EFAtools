kmo_cor <- efa_kmo(test_models$baseline$cormat)
kmo_raw <- efa_kmo(GRiPS_raw)
dat_nonames <- test_models$baseline$cormat
colnames(dat_nonames) <- NULL
kmo_nona <- efa_kmo(dat_nonames)

test_that("output class and dimensions are correct", {
  expect_identical(class(kmo_cor), c("efa_kmo", "KMO"))
  expect_output(str(kmo_cor), "List of 3")
  expect_identical(class(kmo_raw), c("efa_kmo", "KMO"))
  expect_output(str(kmo_raw), "List of 3")
  expect_identical(class(kmo_nona), c("efa_kmo", "KMO"))
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

test_that("per-variable and overall KMO match psych", {
  skip_on_cran()
  skip_if_not_installed("psych")

  psych_cor <- psych::KMO(test_models$baseline$cormat)
  expect_equal(unname(kmo_cor$KMO_i), unname(psych_cor$MSAi), tolerance = 1e-4)
  expect_equal(kmo_cor$KMO, unname(psych_cor$MSA), tolerance = 1e-4)

  expect_equal(unname(kmo_raw$KMO_i),
               unname(psych::KMO(stats::cor(GRiPS_raw))$MSAi),
               tolerance = 1e-4)
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

test_that("errors are thrown correctly", {
  expect_error(efa_kmo(1:5), class = "efa_input_not_matrix")
  expect_message(efa_kmo(GRiPS_raw), class = "efa_cor_from_data")
  expect_error(efa_kmo(sing_raw), class = "efa_cor_singular")
  expect_error(efa_kmo(sing_cor), class = "efa_cor_singular")
  expect_warning(efa_kmo(cor_nposdef), class = "efa_cor_smoothed")
})

test_that("print output is stable", {
  local_reproducible_output()

  # high KMO: tick, verdict, and the per-variable values
  expect_snapshot(print(kmo_cor), transform = scrub_num)

  # low KMO: cross and the "not suitable" verdict
  kmo_low <- structure(list(KMO = 0.45, KMO_i = c(V1 = 0.40, V2 = 0.50, V3 = 0.45)),
                       class = c("efa_kmo", "KMO"))
  expect_snapshot(print(kmo_low), transform = scrub_num)

  # KMO value not available
  kmo_na <- structure(list(KMO = NA_real_, KMO_i = NA_real_),
                      class = c("efa_kmo", "KMO"))
  expect_snapshot(print(kmo_na), transform = scrub_num)
})

rm(kmo_cor, kmo_raw, dat_nonames, kmo_nona)
