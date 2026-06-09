test_that("efa_retention plot methods return ggplot objects", {
  ekc <- EKC(test_models$baseline$cormat, N = 500)
  p_ekc <- plot(ekc)
  expect_s3_class(p_ekc, "ggplot")
  expect_equal(ekc$results[[1]]$plot_type, "eigen")

  hull <- HULL(test_models$baseline$cormat, N = 500, method = "ML")
  p_hull <- plot(hull)
  expect_s3_class(p_hull, "ggplot")
  expect_equal(hull$results[[1]]$plot_type, "hull")

  kgc <- KGC(test_models$baseline$cormat)
  p_kgc <- plot(kgc)
  expect_s3_class(p_kgc, "ggplot")
  expect_equal(.retention_record(kgc, "PCA")$plot_type, "eigen")
})

test_that("plot-less criteria return NULL with a message", {
  for (obj in list(NEST(test_models$baseline$cormat, N = 500),
                   MAP(test_models$baseline$cormat),
                   SMT(test_models$baseline$cormat, N = 500))) {
    expect_message(p <- plot(obj), "No plot is available")
    expect_null(p)
  }
})

test_that("EKC eigen plot is visually stable", {
  skip_if_not_installed("vdiffr")

  ekc <- EKC(test_models$baseline$cormat, N = 500)
  vdiffr::expect_doppelganger("EKC eigen plot", plot(ekc))

  ekc_both <- EKC(test_models$baseline$cormat, N = 500,
                  type = c("BvA2017", "AM2019"))
  vdiffr::expect_doppelganger("EKC eigen plot both types", plot(ekc_both))
})

test_that("HULL hull plot is visually stable", {
  skip_if_not_installed("vdiffr")

  hull <- HULL(test_models$baseline$cormat, N = 500, method = "ML")
  vdiffr::expect_doppelganger("HULL hull plot", plot(hull))
})

test_that("KGC eigen plot is visually stable", {
  skip_if_not_installed("vdiffr")

  kgc <- KGC(test_models$baseline$cormat)
  vdiffr::expect_doppelganger("KGC eigen plot", plot(kgc))
})
