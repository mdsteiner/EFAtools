test_that(".new_efa_retention rejects unknown criterion ids", {
  rec <- list(list(name = "X", label = "X", n_factors = 1, plot_type = "none"))
  expect_error(
    .new_efa_retention("NOPE", results = rec, settings = list()),
    class = "efa_unknown_criterion"
  )
})

test_that(".new_efa_retention builds the documented shape for a known id", {
  rec <- list(
    list(name = "BvA2017", label = "Original", n_factors = 3, plot_type = "eigen"),
    list(name = "AM2019", label = "Adapted", n_factors = 2, plot_type = "eigen")
  )
  out <- .new_efa_retention("EKC", results = rec, settings = list(a = 1))

  expect_s3_class(out, "efa_retention")
  expect_equal(unname(out$criterion["label"]), "Empirical Kaiser Criterion")
  expect_named(out$n_factors, c("BvA2017", "AM2019"))
  expect_equal(unname(out$n_factors), c(3, 2))
  expect_equal(out$status, "ok")
})

test_that("format.efa_retention returns plain text even when styling is enabled", {
  # Force colours on so a styled print() embeds ANSI; format() must not.
  old <- options(cli.num_colors = 256)
  on.exit(options(old), add = TRUE)

  ekc <- EKC(test_models$baseline$cormat, N = 500)
  styled <- utils::capture.output(print(ekc))
  plain <- format(ekc)

  expect_true(any(grepl("\033", styled, fixed = TRUE)))
  expect_false(any(grepl("\033", plain, fixed = TRUE)))
})
