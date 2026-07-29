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
})

test_that(".n_factors_ctl() defaults mirror efa_retain()'s", {
  # efa_power() calls .n_factors_ctl() relying on its defaults, so a default drifting from
  # efa_retain()'s would silently retain factors under undocumented settings.
  ctl <- .n_factors_ctl()
  shared <- intersect(names(ctl), names(formals(efa_retain)))
  ret <- lapply(formals(efa_retain)[shared], eval)
  # efa_retain() resolves these with a single-choice match.arg(), so its effective default
  # is the first element; every other default (including the several.ok ones, whose whole
  # vector is the default) has to match .n_factors_ctl()'s in full
  single <- c("use", "cor_method", "estimator", "eigen_type_HULL", "decision_rule")
  ret[single] <- lapply(ret[single], `[`, 1L)
  expect_equal(ctl[shared], ret[shared])
})

test_that("format.efa_retention is the source of truth and honours the colour state", {
  ekc <- EKC(test_models$baseline$cormat, N = 500)

  # print() is exactly cat(format(x), sep = "\n"), so the two agree line for line.
  expect_identical(utils::capture.output(print(ekc)), format(ekc))

  old <- options(cli.num_colors = 256)
  on.exit(options(old), add = TRUE)

  # With colours on the report embeds ANSI ...
  expect_true(any(grepl("\033", format(ekc), fixed = TRUE)))

  # ... and with colours off it is plain.
  options(cli.num_colors = 1)
  expect_false(any(grepl("\033", format(ekc), fixed = TRUE)))
})
