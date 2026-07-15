## Tests for .reliability_result() -- the tidy long output built from computed
## coefficient matrices -- and its format/print method.

efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
sl_mod <- efa_schmid_leiman(efa_mod, estimate_control = estimate_control(type = "EFAtools"), method = "PAF")
fc <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2
spec <- list(g_load = sl_mod$sl[, "g"], s_load = sl_mod$sl[, c("F1", "F2", "F3")],
             u2 = sl_mod$sl[, "u2"], map = fc, cormat = test_models$baseline$cormat,
             var_names = rownames(sl_mod$sl), fac_names = c("F1", "F2", "F3"))

# A full coefficient matrix (also carries the core's CR/AVE columns, which the result
# must drop).
mat <- .reliability_core(spec, "correlation", add_ind = TRUE, add_rel = TRUE)
res <- .reliability_result(mat, settings = list(variance = "correlation"))

# Multigroup: a named list of matrices, one entry per group.
res_mg <- .reliability_result(list(GroupA = mat, GroupB = mat),
                              settings = list(variance = "correlation"))

# Single factor: a one-row general-factor matrix carrying only omega total and H.
g_sf <- rep(0.7, 6)
u2_sf <- 1 - g_sf^2
mat_sf <- matrix(c(sum(g_sf)^2 / (sum(g_sf)^2 + sum(u2_sf)), .h_index(g_sf)),
                 nrow = 1, dimnames = list("g", c("tot", "H")))
res_sf <- .reliability_result(mat_sf, settings = list(variance = "sums_load"))

test_that("the result is a long data.frame with the documented schema", {
  expect_s3_class(res, "efa_reliability")
  expect_s3_class(res, "data.frame")
  expect_identical(names(res), c("coefficient", "level", "factor", "group", "value"))
  expect_type(res$coefficient, "character")
  expect_type(res$level, "character")
  expect_type(res$factor, "character")
  expect_type(res$group, "character")
  expect_type(res$value, "double")

  # level is general for the g row, group for the group factors, nothing else.
  expect_setequal(res$level, c("general", "group"))
  expect_true(all(res$level[res$factor == "g"] == "general"))
  expect_true(all(res$level[res$factor != "g"] == "group"))

  # A single unnamed group carries an NA group label.
  expect_true(all(is.na(res$group)))
})

test_that("only the surfaced menu is carried; CR/AVE/GLB/beta never appear", {
  # mat carries CR and AVE, but the result surfaces neither, nor any deferred surrogate.
  expect_true(all(c("CR", "AVE") %in% colnames(mat)))
  expect_false(any(c("CR", "AVE", "GLB", "beta") %in% res$coefficient))
  expect_setequal(res$coefficient,
                  c("omega_total", "omega_hierarchical", "omega_subscale",
                    "alpha", "H", "ECV", "PUC"))
})

test_that("the kind attribute tags each surfaced coefficient", {
  kind <- attr(res, "kind")
  expect_type(kind, "character")
  expect_identical(unname(kind[c("ECV", "PUC")]), c("common_variance", "common_variance"))
  expect_true(all(kind[c("omega_total", "omega_hierarchical", "omega_subscale",
                         "alpha", "H")] == "reliability"))
})

test_that("the settings attribute round-trips", {
  expect_identical(attr(res, "settings"), list(variance = "correlation"))
})

test_that("NA cells are dropped; ECV and PUC appear only for the general factor", {
  expect_false(anyNA(res$value))
  expect_true(all(res$level[res$coefficient %in% c("ECV", "PUC")] == "general"))
  # No group-factor row carries ECV or PUC.
  grp_rows <- res[res$level == "group", ]
  expect_false(any(grp_rows$coefficient %in% c("ECV", "PUC")))
})

test_that("a named list yields one labelled block per group", {
  expect_identical(unique(res_mg$group), c("GroupA", "GroupB"))
  # Each group holds the same coefficient set.
  expect_setequal(res_mg$coefficient[res_mg$group == "GroupA"],
                  res_mg$coefficient[res_mg$group == "GroupB"])
})

test_that("an unnamed list labels groups by position rather than collapsing them", {
  # A list is always multigroup; unnamed entries must get distinct positional labels so
  # each group renders as its own block (they must not all fall to a single NA label).
  res_un <- .reliability_result(list(mat, mat_sf))
  expect_identical(unique(res_un$group), c("1", "2"))
  expect_false(anyNA(res_un$group))
  # Both groups keep their own coefficients.
  expect_setequal(res_un$coefficient[res_un$group == "2"], c("omega_total", "H"))
})

test_that("a matrix with no surfaced columns yields an empty result, not an error", {
  # A CR/AVE-only matrix carries no registry column; the builder returns a 0-row
  # efa_reliability rather than aborting on a NULL assembly.
  cr_only <- matrix(c(0.9, 0.8), nrow = 1, dimnames = list("g", c("CR", "AVE")))
  expect_no_error(empty <- .reliability_result(cr_only))
  expect_s3_class(empty, "efa_reliability")
  expect_identical(nrow(empty), 0L)
})

test_that("a single-factor matrix carries only its available coefficients", {
  expect_setequal(res_sf$coefficient, c("omega_total", "H"))
  expect_true(all(res_sf$factor == "g"))
  expect_true(all(res_sf$level == "general"))
})

test_that("print output is stable", {
  local_reproducible_output()

  # single group, full menu (reliability and common-variance blocks)
  expect_snapshot(print(res), transform = scrub_num)

  # multiple groups (per-group headers)
  expect_snapshot(print(res_mg), transform = scrub_num)

  # single factor (reliability block only, omega total and H)
  expect_snapshot(print(res_sf), transform = scrub_num)
})

rm(efa_mod, sl_mod, fc, spec, mat, res, res_mg, g_sf, u2_sf, mat_sf, res_sf)
