# The estimation settings must reach every factor-retention criterion that fits a model, and
# only through the control object: a criterion that accepts a control but never passes it on to
# its fit is indistinguishable, from the outside, from one that ignores it. Each criterion below
# is therefore shown to HONOUR the control, not merely to take it.
#
# Two ways to show it, one per estimator:
#   * PAF (what the criteria that fit use by default): `max_iter = 1` halts the iterations after
#     a single step, so the communalities -- and everything derived from them -- must change.
#   * ML: `start_method = NA` leaves the optimiser without starting values, which efa_fit()
#     rejects with a classed error. That the error surfaces proves the control reached the fit.
#     This is the only route for efa_smt(), whose models are maximum likelihood by definition,
#     so of the estimation knobs only `start_method` applies to it.
#
# CD, EKC and MAP fit no model (CD's `max_iter` caps its own comparison-data generation, not an
# EFA), so they take no control at all.

cormat <- test_models$baseline$cormat

pinned <- estimate_control(max_iter = 1)
no_start <- estimate_control(start_method = NA)

test_that("efa_hull() honours the estimation control", {
  args <- list(cormat, N = 500, method = "PAF", gof = "CAF", n_datasets = 10)

  set.seed(42L)
  default <- suppressWarnings(do.call(efa_hull, args))
  set.seed(42L)
  tuned <- suppressWarnings(do.call(efa_hull, c(args, list(estimate_control = pinned))))

  # the goodness-of-fit values the hull is built from come from the 0..J factor fits
  expect_false(isTRUE(all.equal(default$results, tuned$results)))

  # the control reaches those fits on the ML route too
  expect_error(
    efa_hull(cormat, N = 500, method = "ML", n_datasets = 10,
             estimate_control = no_start),
    class = "efa_ml_start_missing")
})

test_that("efa_kgc() honours the estimation control", {
  default <- efa_kgc(cormat, eigen_type = "EFA")
  tuned <- suppressWarnings(efa_kgc(cormat, eigen_type = "EFA",
                                    estimate_control = pinned))

  expect_false(isTRUE(all.equal(.retention_record(default, "EFA")$y,
                                .retention_record(tuned, "EFA")$y)))
})

test_that("efa_scree() honours the estimation control", {
  default <- efa_scree(cormat, eigen_type = "EFA")
  tuned <- suppressWarnings(efa_scree(cormat, eigen_type = "EFA",
                                      estimate_control = pinned))

  expect_false(isTRUE(all.equal(.retention_record(default, "EFA")$y,
                                .retention_record(tuned, "EFA")$y)))
})

test_that("efa_parallel() honours the estimation control", {
  args <- list(cormat, N = 500, eigen_type = "EFA", n_datasets = 2)

  set.seed(42L)
  default <- suppressWarnings(do.call(efa_parallel, args))
  set.seed(42L)
  tuned <- suppressWarnings(do.call(efa_parallel,
                                    c(args, list(estimate_control = pinned))))

  # the real-data EFA eigenvalues; the shared seed rules out a mere simulation difference
  expect_false(isTRUE(all.equal(.retention_record(default, "EFA")$y,
                                .retention_record(tuned, "EFA")$y)))
})

test_that("efa_nest() honours the estimation control", {
  args <- list(cormat, N = 500, n_datasets = 2)

  set.seed(42L)
  default <- suppressWarnings(do.call(efa_nest, args))
  set.seed(42L)
  tuned <- suppressWarnings(do.call(efa_nest, c(args, list(estimate_control = pinned))))

  # the reference eigenvalues are simulated from the (nf - 1)-factor PAF reference models
  expect_false(isTRUE(all.equal(.retention_record(default, "NEST")$reference,
                                .retention_record(tuned, "NEST")$reference)))

  expect_error(
    do.call(efa_nest, c(args, list(method = "ML", estimate_control = no_start))),
    class = "efa_ml_start_missing")
})

test_that("efa_smt() honours the estimation control", {
  # SMT could not be tuned at all before: it had no way to receive a control
  expect_error(efa_smt(cormat, N = 500, estimate_control = no_start),
               class = "efa_ml_start_missing")

  # the other ML starting values are accepted and still produce a valid result
  expect_s3_class(
    efa_smt(cormat, N = 500,
            estimate_control = estimate_control(start_method = "factanal")),
    "efa_retention")
})

test_that("efa_schmid_leiman() honours the estimation control", {
  efa_mod <- efa_fit(cormat, n_factors = 3, N = 500, method = "PAF",
                     rotation = "promax")

  # type = "none" resolves nothing from a preset, so the second-order fit only runs if the
  # knobs supplied in the control actually reach it
  expect_no_error(
    efa_schmid_leiman(efa_mod, method = "PAF",
                      estimate_control = estimate_control(
                        type = "none", init_comm = "smc", criterion = 1e-3,
                        criterion_type = "sum", max_iter = 300, abs_eigen = TRUE)))

  expect_error(
    efa_schmid_leiman(efa_mod, method = "ML", estimate_control = no_start),
    class = "efa_ml_start_missing")
})

test_that("efa_retain() threads the estimation control into every criterion that fits", {
  ids <- c("KGC", "NEST", "PARALLEL", "SCREE")
  args <- list(cormat, N = 500, suitability = FALSE, criteria = ids, method = "PAF",
               eigen_type_other = "EFA", n_datasets = 2, n_datasets_nest = 2)

  set.seed(42L)
  default <- suppressWarnings(suppressMessages(do.call(efa_retain, args)))
  set.seed(42L)
  tuned <- suppressWarnings(suppressMessages(
    do.call(efa_retain, c(args, list(estimate_control = pinned)))))

  for (id in ids) {
    expect_false(isTRUE(all.equal(default$outputs[[id]]$results,
                                  tuned$outputs[[id]]$results)),
                 info = id)
  }

  # the control reaches SMT as well: its models are fitted with ML, so the missing starting
  # values abort that criterion (efa_retain() excludes a failing criterion and carries on)
  smt <- suppressWarnings(suppressMessages(
    efa_retain(cormat, N = 500, suitability = FALSE, criteria = c("EKC", "SMT"),
               estimate_control = no_start)))
  expect_named(smt$not_run, "SMT")
  expect_named(smt$outputs, "EKC")

  expect_identical(tuned$settings$estimate_control, pinned)
})

test_that("efa_retain() rejects a flat knob that R would partial-match to a formal", {
  # `max_iter` is a unique prefix of the `max_iter_CD` formal, so R matched it there: the knob
  # silently capped the comparison-data iterations instead of tuning the fits, and never
  # reached the dots -- or the guard. `max_iter_CD` now sits behind `...`, where R matches it by
  # its full name only, so the knob lands in the dots and is named for what it is.
  expect_error(efa_retain(cormat, N = 500, criteria = "EKC", max_iter = 500),
               class = "efa_flat_knob_in_dots")

  # the same holds when the knob arrives through a caller's own `...`, where the partial match
  # would happen before efa_retain() ever ran and no inspection of the call could catch it
  forwarder <- function(...) efa_retain(cormat, N = 500, criteria = "EKC", ...)
  expect_error(forwarder(max_iter = 500), class = "efa_flat_knob_in_dots")

  # the comparison-data argument it was mistaken for still works, spelled out
  expect_no_error(suppressWarnings(suppressMessages(
    efa_retain(cormat, N = 500, criteria = "EKC", max_iter_CD = 5))))
  expect_identical(
    suppressWarnings(suppressMessages(
      efa_retain(cormat, N = 500, criteria = "EKC",
                 max_iter_CD = 5)))$settings$max_iter_CD,
    5)
})

test_that("a bare flat knob is rejected even when the call runs no fit", {
  # with eigen_type = "PCA" no model is fitted, so the knob never reached efa_fit()'s guard
  # and was dropped without a word
  expect_error(efa_kgc(cormat, eigen_type = "PCA", max_iter = 500),
               class = "efa_flat_knob_in_dots")
  expect_error(efa_scree(cormat, eigen_type = "PCA", init_comm = "mac"),
               class = "efa_flat_knob_in_dots")
  expect_error(efa_parallel(cormat, N = 500, eigen_type = "PCA", n_datasets = 2,
                            type = "SPSS"),
               class = "efa_flat_knob_in_dots")
  expect_error(efa_nest(cormat, N = 500, n_datasets = 2, criterion = 1e-4),
               class = "efa_flat_knob_in_dots")
  expect_error(efa_hull(cormat, N = 500, gof = "CAF", abs_eigen = TRUE),
               class = "efa_flat_knob_in_dots")
  expect_error(
    efa_schmid_leiman(efa_fit(cormat, n_factors = 3, N = 500, method = "PAF",
                              rotation = "promax"),
                      max_iter = 300),
    class = "efa_flat_knob_in_dots")
})

test_that("an estimation control that is not a control object is rejected", {
  expect_error(efa_smt(cormat, N = 500, estimate_control = list(max_iter = 1)),
               class = "efa_control_input")
  expect_error(efa_kgc(cormat, eigen_type = "PCA", estimate_control = "SPSS"),
               class = "efa_control_input")
  expect_error(efa_retain(cormat, N = 500, criteria = "EKC",
                          estimate_control = rotate_control()),
               class = "efa_control_input")
})

test_that("N_FACTORS() now tunes the criteria its dots always claimed to reach", {
  # NEST and SMT never received the repacked control (the registry did not pass it on), so a
  # legacy call silently ran the default preset there. A coarse `criterion` ends the PAF
  # iterations almost at once; unlike `max_iter`, which R partial-matches to the frozen
  # `max_iter_CD` formal, it really does travel through the wrapper's dots.
  args <- list(cormat, N = 500, suitability = FALSE, criteria = "NEST",
               method = "PAF", n_datasets_nest = 2)

  set.seed(42L)
  old <- suppressWarnings(suppressMessages(
    do.call(N_FACTORS, c(args, list(criterion = 0.5)))))
  set.seed(42L)
  new <- suppressWarnings(suppressMessages(
    do.call(efa_retain, c(args, list(estimate_control = estimate_control(criterion = 0.5))))))
  set.seed(42L)
  default <- suppressWarnings(suppressMessages(do.call(efa_retain, args)))

  expect_equal(old$outputs$NEST$results, new$outputs$NEST$results)
  expect_false(isTRUE(all.equal(old$outputs$NEST$results,
                                default$outputs$NEST$results)))
})
