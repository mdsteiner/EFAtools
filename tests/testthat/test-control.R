# Unit tests for the estimate_control() / rotate_control() control constructors.

test_that("constructors return classed lists with the documented fields", {
  ec <- estimate_control()
  expect_s3_class(ec, "efa_estimate_control")
  expect_equal(ec$type, "EFAtools")
  expect_true(is.na(ec$init_comm))
  expect_true(is.na(ec$criterion))
  expect_true(is.na(ec$abs_eigen))
  expect_equal(ec$start_method, "psych")
  # The two FIML knobs are not preset-driven, so they carry the EM's own defaults rather
  # than the NA "resolve from the type preset" sentinel.
  expect_equal(ec$fiml_max_iter, 500)
  expect_equal(ec$fiml_tol, 1e-5)
  expect_named(ec, c("type", "init_comm", "criterion", "criterion_type",
                     "max_iter", "abs_eigen", "start_method",
                     "fiml_max_iter", "fiml_tol"))

  rc <- rotate_control()
  expect_s3_class(rc, "efa_rotate_control")
  expect_equal(rc$type, "EFAtools")
  expect_true(rc$normalize)
  expect_equal(rc$precision, 1e-5)
  expect_true(is.na(rc$p_type))
  expect_equal(rc$random_starts, 100)
  expect_type(rc$extra_args, "list")
  expect_length(rc$extra_args, 0)
  expect_named(rc, c("type", "normalize", "precision", "order_type",
                     "varimax_type", "p_type", "k", "random_starts", "extra_args"))
})

test_that("rotate_control captures criterion-specific extras in extra_args", {
  rc <- rotate_control(type = "EFAtools", k = 3, gam = 0.5, delta = 0.02)
  expect_equal(rc$k, 3)
  expect_equal(rc$extra_args$gam, 0.5)
  expect_equal(rc$extra_args$delta, 0.02)
})

test_that("rotate_control rejects a misplaced tuning knob among its extras", {
  # the fit drops an extra whose name collides with one of its formals, so a knob
  # smuggled in here would be silently ignored -- refuse it at construction instead
  expect_error(rotate_control(max_iter = 2), class = "efa_control_input")
  expect_error(rotate_control(criterion = 1e-4), class = "efa_control_input")
  expect_error(rotate_control(P_type = "norm"), class = "efa_control_input")
  expect_error(rotate_control(randomStarts = 7), class = "efa_control_input")
})

test_that("bad estimate_control inputs abort with a classed condition", {
  expect_error(estimate_control(init_comm = "xyz"), class = "efa_control_input")
  expect_error(estimate_control(criterion = -1), class = "efa_control_input")
  expect_error(estimate_control(criterion = 0), class = "efa_control_input")
  expect_error(estimate_control(criterion_type = "foo"), class = "efa_control_input")
  expect_error(estimate_control(max_iter = 0), class = "efa_control_input")
  expect_error(estimate_control(max_iter = 2.5), class = "efa_control_input")
  expect_error(estimate_control(abs_eigen = "yes"), class = "efa_control_input")
  expect_error(estimate_control(start_method = "bad"), class = "efa_control_input")
  # NaN is not the unset sentinel; it must be rejected, not stored.
  expect_error(estimate_control(criterion = NaN), class = "efa_control_input")
  expect_error(estimate_control(init_comm = NaN), class = "efa_control_input")
  expect_error(estimate_control(abs_eigen = NaN), class = "efa_control_input")
  # a criterion at or above 1 is not a convergence tolerance; the bound matches the one the
  # fit enforces, so a control can never carry a value the fit would go on to reject.
  expect_error(estimate_control(criterion = 5), class = "efa_control_input")
  expect_error(estimate_control(criterion = 1), class = "efa_control_input")
  # the bound is exclusive at 1 only: anything below it is a legitimate tolerance
  expect_s3_class(estimate_control(criterion = 0.999), "efa_estimate_control")
  # an invalid `type` is rejected with the same class as the other knobs
  expect_error(estimate_control(type = "bogus"), class = "efa_control_input")
  expect_error(rotate_control(type = "bogus"), class = "efa_control_input")
  # the FIML EM knobs mirror the bounds the EM itself enforces, and -- not being
  # preset-driven -- do not accept the NA sentinel either
  expect_error(estimate_control(fiml_max_iter = 0), class = "efa_control_input")
  expect_error(estimate_control(fiml_max_iter = 2.5), class = "efa_control_input")
  expect_error(estimate_control(fiml_max_iter = NA), class = "efa_control_input")
  expect_error(estimate_control(fiml_tol = 0), class = "efa_control_input")
  expect_error(estimate_control(fiml_tol = -1e-4), class = "efa_control_input")
  expect_error(estimate_control(fiml_tol = NA), class = "efa_control_input")
  # A tolerance at or above 1 is met on the first EM iteration, which would return the
  # starting moments while reporting convergence; the EM itself only demands tol > 0, so the
  # control is the only place that can refuse it (as it does for `criterion`).
  expect_error(estimate_control(fiml_tol = 1), class = "efa_control_input")
  expect_error(estimate_control(fiml_tol = 5), class = "efa_control_input")
  expect_s3_class(estimate_control(fiml_max_iter = 2000, fiml_tol = 0.999),
                  "efa_estimate_control")
})

test_that("efa_fit rejects unused dots when no rotation runs", {
  # with rotation = "none" no rotation engine runs, so nothing can consume the dots;
  # this also catches a misspelled control argument instead of silently ignoring it
  expect_error(
    efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500, rotation = "none",
            estimate_controls = estimate_control(max_iter = 500)),
    class = "efa_unused_dots")
  expect_error(
    efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500, gam = 0.5),
    class = "efa_unused_dots")
})

test_that("efa_fit rejects a dot argument the selected rotation does not consume", {
  cm <- test_models$baseline$cormat

  # a misspelled criterion parameter would silently run the engine default (gam = 0)
  expect_error(
    efa_fit(cm, n_factors = 3, N = 500, rotation = "oblimin", gamma = 0.5),
    class = "efa_unused_dots")
  # another criterion's parameter is just as unconsumable as a misspelling
  expect_error(
    efa_fit(cm, n_factors = 3, N = 500, rotation = "geominT", gam = 0.5),
    class = "efa_unused_dots")
  # varimax and promax run outside the native engines and take no extras at all
  expect_error(
    efa_fit(cm, n_factors = 3, N = 500, rotation = "varimax", bogus_arg = 1),
    class = "efa_unused_dots")
  expect_error(
    efa_fit(cm, n_factors = 3, N = 500, rotation = "promax", maxit = 100),
    class = "efa_unused_dots")

  # the consumable extras still reach the engine: oblimin's gam changes the criterion.
  # gam = 0.5 rewards correlated factors strongly enough to trip the extreme-Phi check on
  # this fixture, which is beside the point here -- only that the extra reaches the engine.
  fit_default <- efa_fit(cm, n_factors = 3, N = 500, rotation = "oblimin")
  fit_gam <- suppressWarnings(
    efa_fit(cm, n_factors = 3, N = 500, rotation = "oblimin", gam = 0.5))
  expect_false(identical(fit_default$rot_loadings, fit_gam$rot_loadings))
  expect_no_error(efa_fit(cm, n_factors = 3, N = 500, rotation = "geominQ", delta = 0.05))

  # a control-carried extra stays lenient per rotation: one control may serve fits whose
  # rotation never consumes it (efa_average() carries maxit across its whole grid)
  expect_no_error(efa_fit(cm, n_factors = 3, N = 500, rotation = "varimax",
                          rotate_control = rotate_control(maxit = 5e4)))
  expect_no_error(efa_fit(cm, n_factors = 3, N = 500, rotation = "none",
                          rotate_control = rotate_control(maxit = 5e4)))
})

test_that("the fit stays the backstop for a control edited after construction", {
  # The constructor's bounds mirror the fit's, so a control it builds can never carry a value the
  # fit rejects. A control is a plain list, though, so it can be edited afterwards -- the fit
  # keeps its own guard for that, and it must not be mistaken for dead code.
  ec <- estimate_control(criterion = 1e-3)
  ec$criterion <- 1

  expect_error(
    efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500, estimator = "PAF",
            estimate_control = ec),
    class = "efa_criterion_too_large")
})

test_that("efa_fit's argument assertions abort with a classed condition", {
  cm <- test_models$baseline$cormat

  expect_error(efa_fit(cm, n_factors = 3, N = -5), class = "efa_invalid_argument")
  expect_error(efa_fit(cm, n_factors = -1, N = 500), class = "efa_invalid_argument")

  # a value routed through a control reaches the same guard
  ec <- estimate_control()
  ec$max_iter <- -1
  expect_error(
    efa_fit(cm, n_factors = 3, N = 500, estimator = "PAF", estimate_control = ec),
    class = "efa_invalid_argument")
})

test_that("estimate_control accepts the start_method spellings the flat interface took", {
  # start_method only governs the ML optimiser, so NA ("not needed here") is admissible; the
  # fit rejects it, and only when the estimator really is ML.
  expect_true(is.na(estimate_control(start_method = NA)$start_method))
  # abbreviations resolve, as checkmate::matchArg() did for the flat interface
  expect_identical(estimate_control(start_method = "fact")$start_method, "factanal")
  expect_identical(estimate_control(start_method = "psy")$start_method, "psych")
  expect_error(estimate_control(start_method = "bad"), class = "efa_control_input")
})

test_that("bad rotate_control inputs abort with a classed condition", {
  # an extra no rotation engine can consume is a misspelling; refuse to store it
  expect_error(rotate_control(gamma = 0.5), class = "efa_control_input")
  expect_error(rotate_control(eps = 1e-6), class = "efa_control_input")
  expect_error(rotate_control(order_type = "x"), class = "efa_control_input")
  expect_error(rotate_control(varimax_type = "x"), class = "efa_control_input")
  expect_error(rotate_control(p_type = "x"), class = "efa_control_input")
  expect_error(rotate_control(k = -1), class = "efa_control_input")
  expect_error(rotate_control(precision = 0), class = "efa_control_input")
  expect_error(rotate_control(random_starts = -1), class = "efa_control_input")
  expect_error(rotate_control(random_starts = 2.5), class = "efa_control_input")
  expect_error(rotate_control(k = NaN), class = "efa_control_input")
  expect_error(rotate_control(p_type = NaN), class = "efa_control_input")
  # a rotation tolerance above 1 is not a tolerance; the bound matches the fit's.
  expect_error(rotate_control(precision = 2), class = "efa_control_input")
  # normalize, precision and random_starts have no NA state.
  expect_error(rotate_control(normalize = NA), class = "efa_control_input")
  expect_error(rotate_control(precision = NA), class = "efa_control_input")
})

test_that("rotate_control accepts random_starts = 0 (warm start only)", {
  # 0 is a meaningful setting, not a missing one: it runs the rotation from its warm start
  # with no random restarts, which the flat interface always allowed.
  expect_equal(rotate_control(random_starts = 0)$random_starts, 0)
})

test_that("NA-defaulted knobs are accepted (they resolve from the preset later)", {
  expect_s3_class(estimate_control(init_comm = NA, criterion = NA,
                                   criterion_type = NA, max_iter = NA,
                                   abs_eigen = NA),
                  "efa_estimate_control")
  expect_s3_class(rotate_control(order_type = NA, varimax_type = NA,
                                 p_type = NA, k = NA),
                  "efa_rotate_control")
})

# The core requirement: a control object built from a given type resolves (through
# the same resolver the fitting cores use) to the same settings EFA() would use.

test_that("estimate_control(type) resolves to the settings EFA() uses for PAF", {
  fields <- c("init_comm", "criterion", "criterion_type", "max_iter", "abs_eigen")
  for (ty in c("EFAtools", "psych", "SPSS")) {
    ec <- estimate_control(type = ty)
    resolved <- .resolve_settings(
      type = ty,
      user = list(init_comm = ec$init_comm, criterion = ec$criterion,
                  criterion_type = ec$criterion_type, max_iter = ec$max_iter,
                  abs_eigen = ec$abs_eigen),
      preset = .efa_presets$PAF
    )
    efa <- suppressWarnings(
      EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
          method = "PAF", type = ty))
    for (nm in fields) expect_equal(resolved[[nm]], efa$settings[[nm]])
  }
})

test_that("rotate_control(type) resolves to the settings EFA() uses per rotation", {
  cormat <- test_models$baseline$cormat

  for (ty in c("EFAtools", "psych", "SPSS")) {
    rc <- rotate_control(type = ty)

    # promax: the richest preset (normalize, P_type, order_type, varimax_type, k).
    resolved_p <- .resolve_settings(
      type = ty,
      user = list(normalize = rc$normalize, P_type = rc$p_type,
                  order_type = rc$order_type, varimax_type = rc$varimax_type,
                  k = rc$k),
      preset = .efa_presets$PROMAX
    )
    efa_p <- suppressWarnings(
      EFA(cormat, n_factors = 3, N = 500, rotation = "promax", type = ty))
    for (nm in c("normalize", "P_type", "order_type", "varimax_type", "k")) {
      expect_equal(resolved_p[[nm]], efa_p$settings[[nm]])
    }

    # varimax: normalize, order_type, varimax_type.
    resolved_v <- .resolve_settings(
      type = ty,
      user = list(normalize = rc$normalize, order_type = rc$order_type,
                  varimax_type = rc$varimax_type),
      preset = .efa_presets$VARIMAX
    )
    efa_v <- suppressWarnings(
      EFA(cormat, n_factors = 3, N = 500, rotation = "varimax", type = ty))
    for (nm in c("normalize", "order_type", "varimax_type")) {
      expect_equal(resolved_v[[nm]], efa_v$settings[[nm]])
    }

    # oblimin: only normalize and order_type are preset-driven.
    resolved_o <- .resolve_settings(
      type = ty,
      user = list(normalize = rc$normalize, order_type = rc$order_type),
      preset = .efa_presets$ROTATE_OBLQ
    )
    efa_o <- suppressWarnings(
      EFA(cormat, n_factors = 3, N = 500, rotation = "oblimin", type = ty))
    for (nm in c("normalize", "order_type")) {
      expect_equal(resolved_o[[nm]], efa_o$settings[[nm]])
    }
  }
})

test_that("a knob pinned alongside a type is carried through and wins over the preset", {
  ec <- estimate_control(type = "SPSS", max_iter = 10)
  expect_equal(ec$max_iter, 10)

  resolved <- suppressWarnings(.resolve_settings(
    type = "SPSS",
    user = list(init_comm = ec$init_comm, criterion = ec$criterion,
                criterion_type = ec$criterion_type, max_iter = ec$max_iter,
                abs_eigen = ec$abs_eigen),
    preset = .efa_presets$PAF
  ))
  expect_equal(resolved$max_iter, 10)
  # the SPSS preset default (25) is overridden by the pin.
  expect_equal(resolved$criterion_type, "max_individual")
})

test_that("control objects print", {
  testthat::local_reproducible_output()

  expect_snapshot(print(estimate_control(type = "SPSS")))
  expect_snapshot(print(estimate_control(type = "none", init_comm = "smc",
                                         criterion = 1e-3, criterion_type = "sum",
                                         max_iter = 300, abs_eigen = TRUE)))
  expect_snapshot(print(rotate_control(type = "psych")))
  expect_snapshot(print(rotate_control(type = "none", normalize = FALSE,
                                       order_type = "eigen", varimax_type = "kaiser",
                                       p_type = "norm", k = 4, random_starts = 50,
                                       gam = 0.5)))
})
