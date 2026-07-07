map_cor <- MAP(test_models$baseline$cormat)
map_raw <- MAP(GRiPS_raw)

test_that("output class and dimensions are correct", {
  expect_s3_class(map_cor, "efa_retention")
  expect_length(map_cor, 7)
  expect_s3_class(map_raw, "efa_retention")
  expect_length(map_raw, 7)

  expect_named(map_cor$n_factors, c("TR2", "TR4"))
  expect_equal(.retention_record(map_cor, "TR2")$plot_type, "none")
})

test_that("identified number of factors is correct", {
  expect_equal(map_cor$n_factors[["TR2"]], 1)
  expect_equal(map_cor$n_factors[["TR4"]], 3)
})

test_that("criterion series are returned", {
  tr2 <- .retention_record(map_cor, "TR2")
  expect_length(tr2$y, ncol(test_models$baseline$cormat))
  expect_length(tr2$x, ncol(test_models$baseline$cormat))
  expect_equal(tr2$x[1], 0)
  # the suggested m minimizes the criterion
  expect_equal(tr2$x[which.min(tr2$y)], map_cor$n_factors[["TR2"]])
})

test_that("settings are returned correctly", {
  expect_named(map_cor$settings, c("use", "cor_method"))
  expect_equal(map_cor$settings$use, "pairwise.complete.obs")
  expect_equal(map_cor$settings$cor_method, "pearson")
})

test_that("errors are thrown correctly", {
  expect_error(MAP(1:5), class = "efa_input_not_matrix")
  expect_message(MAP(GRiPS_raw), class = "efa_cor_from_data")
})

rm(map_cor, map_raw)
