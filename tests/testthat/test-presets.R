# Unit tests for the type-preset resolver shared by the estimation and rotation
# helpers.

test_that("preset defaults are filled when arguments are NA", {
  res <- .resolve_settings(
    type = "SPSS",
    user = list(init_comm = NA, criterion = NA, criterion_type = NA,
                max_iter = NA, abs_eigen = NA),
    preset = .efa_presets$PAF
  )
  expect_equal(res$init_comm, "smc")
  expect_equal(res$criterion_type, "max_individual")
  expect_equal(res$max_iter, 25)
  expect_equal(res$abs_eigen, TRUE)

  res_psych <- .resolve_settings(
    type = "psych",
    user = list(normalize = TRUE, order_type = NA, varimax_type = NA),
    preset = .efa_presets$VARIMAX
  )
  expect_equal(res_psych$order_type, "eigen")
  expect_equal(res_psych$varimax_type, "svd")
})

test_that("specified arguments are kept and untouched ones take the preset default", {
  res <- suppressWarnings(.resolve_settings(
    type = "EFAtools",
    user = list(init_comm = NA, criterion = NA, criterion_type = "max_individual",
                max_iter = NA, abs_eigen = NA),
    preset = .efa_presets$PAF
  ))
  expect_equal(res$criterion_type, "max_individual")
  expect_equal(res$max_iter, 300)
})

test_that("one consolidated warning is emitted regardless of how many arguments are pinned", {
  expect_warning(
    .resolve_settings(
      type = "EFAtools",
      user = list(init_comm = "smc", criterion = NA, criterion_type = NA,
                  max_iter = NA, abs_eigen = NA),
      preset = .efa_presets$PAF
    ),
    class = "efa_type_override"
  )

  warns <- capture_warnings(.resolve_settings(
    type = "EFAtools",
    user = list(init_comm = "smc", criterion = 0.01, criterion_type = "sum",
                max_iter = 100, abs_eigen = FALSE),
    preset = .efa_presets$PAF
  ))
  expect_length(warns, 1)
})

test_that("normalize is treated as pinned only when FALSE", {
  expect_warning(
    .resolve_settings(
      type = "EFAtools",
      user = list(normalize = FALSE, order_type = NA, varimax_type = NA),
      preset = .efa_presets$VARIMAX
    ),
    class = "efa_type_override"
  )

  expect_no_warning(.resolve_settings(
    type = "EFAtools",
    user = list(normalize = TRUE, order_type = NA, varimax_type = NA),
    preset = .efa_presets$VARIMAX
  ))
})

test_that("type = 'none' passes the user values through without a condition", {
  expect_no_warning(res <- .resolve_settings(
    type = "none",
    user = list(normalize = TRUE, order_type = "eigen", varimax_type = "svd"),
    preset = .efa_presets$VARIMAX
  ))
  expect_equal(res$order_type, "eigen")
  expect_equal(res$varimax_type, "svd")
})

test_that("type = 'none' errors when a required argument is missing", {
  expect_error(
    .resolve_settings(
      type = "none",
      user = list(normalize = TRUE, order_type = NA, varimax_type = "svd"),
      preset = .efa_presets$VARIMAX
    ),
    class = "efa_type_none"
  )

  # normalize has no NA state, so leaving it at its default must not error.
  expect_silent(.resolve_settings(
    type = "none",
    user = list(normalize = TRUE, order_type = "eigen"),
    preset = .efa_presets$ROTATE_ORTH
  ))
})

test_that("preset tables hold the documented defaults", {
  expect_equal(.efa_presets$PAF$SPSS$max_iter, 25)
  expect_equal(.efa_presets$PAF$psych$abs_eigen, FALSE)
  expect_equal(.efa_presets$PROMAX$psych$P_type, "unnorm")
  expect_equal(.efa_presets$ROTATE_ORTH$SPSS$order_type, "ss_factors")
})
