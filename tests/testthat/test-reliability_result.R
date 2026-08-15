## Tests for .reliability_result() -- the tidy long output built from computed
## coefficient matrices -- and its format/print method.

efa_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
               type = "EFAtools", method = "PAF", rotation = "promax")
sl_mod <- efa_schmid_leiman(efa_mod, estimate_control = estimate_control(type = "EFAtools"), estimator = "PAF")
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

test_that("the general factor is the first row, whatever it is labelled", {
  # The general factor is row one of every matrix the core builds; the label on it varies
  # (a single-factor solution carries its own factor's name there), so the level cannot be
  # read off the label. Reading it off the label would also mislabel a group factor a user
  # names "g" as the whole-scale row.
  named_sf <- .reliability_result(
    matrix(c(0.86, 0.87), nrow = 1, dimnames = list("ability", c("tot", "H"))))
  expect_true(all(named_sf$factor == "ability"))
  expect_true(all(named_sf$level == "general"))

  collide <- .reliability_result(
    matrix(c(0.8, 0.7, 0.6, 0.5), nrow = 2,
           dimnames = list(c("F1", "g"), c("tot", "H"))))
  expect_identical(collide$level[collide$factor == "F1"], rep("general", 2L))
  expect_identical(collide$level[collide$factor == "g"], rep("group", 2L))
})

test_that("the correlated-factors note agrees with the table it introduces", {
  local_reproducible_output()
  # The note states whether the two omega columns coincide, so it has to be decided from
  # the values as printed: two coefficients less than one last decimal apart can still
  # straddle a rounding boundary and print as different numbers, which "equals" would
  # then contradict.
  build <- function(tot, sub) {
    m <- cbind(tot = c(0.9, tot), sub = c(NA_real_, sub))
    rownames(m) <- c("g", paste0("F", seq_along(tot)))
    .reliability_result(m, settings = list(variance = "correlation",
                                           no_general = TRUE))
  }
  txt <- function(x) gsub("\\s+", " ", paste(cli::ansi_strip(format(x)), collapse = " "))
  equal_note <- "subscale omega equals its total omega"

  expect_match(txt(build(c(0.769, 0.745), c(0.769, 0.745))), equal_note, fixed = TRUE)

  # Apart by 9e-4 -- inside the last printed decimal, but on either side of it.
  differ <- txt(build(c(0.7694, 0.745), c(0.7685, 0.745)))
  expect_false(grepl(equal_note, differ, fixed = TRUE))
  expect_match(differ, "F1 .769 .768", fixed = TRUE)
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
