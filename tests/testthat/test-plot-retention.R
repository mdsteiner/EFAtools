test_that("efa_retention plot methods return ggplot objects", {
  skip_if_not_slow()
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

  # CD is stochastic, so only smoke-test that the plot builds with its custom
  # y-axis label (no vdiffr baseline)
  set.seed(123)
  cd <- CD(GRiPS_raw, N_pop = 1000, N_samples = 100)
  p_cd <- plot(cd)
  expect_s3_class(p_cd, "ggplot")
  expect_equal(.retention_record(cd, "CD")$y_label, "RMSE eigenvalues")

  scree <- SCREE(test_models$baseline$cormat)
  p_scree <- plot(scree)
  expect_s3_class(p_scree, "ggplot")

  # PARALLEL: smoke-test only (no vdiffr baseline) because the simulated
  # reference lines vary with the RNG state and the future plan's chunking,
  # so an SVG baseline would not be portable
  pa <- PARALLEL(test_models$baseline$cormat, N = 500, eigen_type = c("PCA", "SMC"))
  p_pa <- plot(pa)
  expect_s3_class(p_pa, "ggplot")

  # no real data: the plot shows only the dashed reference series
  pa_nodat <- PARALLEL(N = 20, n_vars = 5, eigen_type = "PCA")
  expect_s3_class(plot(pa_nodat), "ggplot")
})

test_that("plot-less criteria return NULL with a message", {
  for (obj in list(NEST(test_models$baseline$cormat, N = 500),
                   MAP(test_models$baseline$cormat),
                   SMT(test_models$baseline$cormat, N = 500))) {
    expect_message(p <- plot(obj), "No plot is available")
    expect_null(p)
  }
})

test_that("eigen plot of an empty record returns NULL with a message", {
  # e.g. CD on a tiny dataset that suggests 0 factors -> empty x/y record
  obj <- .new_efa_retention(
    "CD",
    results = list(list(name = "CD", label = "Suggested number of factors",
                        n_factors = 0, plot_type = "eigen",
                        x = integer(0), y = numeric(0))),
    settings = list()
  )
  expect_message(p <- plot(obj), "No plot is available")
  expect_null(p)
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

test_that("SCREE eigen plot is visually stable", {
  skip_if_not_installed("vdiffr")

  scree <- SCREE(test_models$baseline$cormat)
  vdiffr::expect_doppelganger("SCREE eigen plot", plot(scree))
})
