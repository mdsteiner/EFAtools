

vec_s <- c("a" = 1, "b" = 2, "c" = 4)
vec_L <- c("A" = 1, "B" = 2, "C" = 4)

mat_s <- matrix(c(0, 0, 0, 1), ncol = 2)
colnames(mat_s) <- c("a", "b")
mat_L <- matrix(c(0, 0, 0, 1), ncol = 2)
colnames(mat_L) <- c("A", "B")

int <- efa_compare(1:10, 1:10)
dec <- efa_compare(c(.1, .2), c(.1, .1))
matr <- efa_compare(matrix(c(1,1,1,2), ncol = 2), matrix(c(1,1,1,1), ncol = 2))

SPSS_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                type = "SPSS", method = "PAF", rotation = "none")
psych_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                 type = "psych", method = "PAF", rotation = "none")
load <- efa_compare(SPSS_PAF$unrot_loadings, SPSS_PAF$unrot_loadings)
load_ro1 <- efa_compare(SPSS_PAF$unrot_loadings, SPSS_PAF$unrot_loadings,
                    reorder = "names")
load_ro2 <- efa_compare(SPSS_PAF$unrot_loadings, SPSS_PAF$unrot_loadings,
                    reorder = "none")

SPSS_PAF_1 <- EFA(test_models$baseline$cormat, n_factors = 1, N = 500,
                type = "SPSS", method = "PAF", rotation = "none")
psych_PAF_1 <- EFA(test_models$baseline$cormat, n_factors = 1, N = 500,
                 type = "psych", method = "PAF", rotation = "none")
load_F1 <- efa_compare(SPSS_PAF_1$unrot_loadings, psych_PAF_1$unrot_loadings)

test_that("output class and dimensions are correct", {
  # the old class string is kept alongside the new one, so `inherits(x, "COMPARE")`
  # in code written against the superseded name still resolves
  expect_identical(class(int), c("efa_compare", "COMPARE"))
  expect_identical(class(dec), c("efa_compare", "COMPARE"))
  expect_identical(class(matr), c("efa_compare", "COMPARE"))
  expect_s3_class(load, "efa_compare")
  expect_s3_class(load_ro1, "efa_compare")
  expect_s3_class(load_ro2, "efa_compare")
  expect_s3_class(load_F1, "efa_compare")

  expect_named(int, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                      "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                      "diff_corres_cross", "g", "settings"))
  expect_named(dec, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                      "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                      "diff_corres_cross", "g", "settings"))
  expect_named(matr, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                      "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                      "diff_corres_cross", "g", "settings"))
  expect_named(load, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                       "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                       "diff_corres_cross", "g", "settings"))
  expect_named(load_ro1, c("diff", "mean_abs_diff", "median_abs_diff",
                           "min_abs_diff", "max_abs_diff", "max_dec", "are_equal",
                           "diff_corres", "diff_corres_cross", "g", "settings"))
  expect_named(load_ro2, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                       "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                       "diff_corres_cross", "g", "settings"))
  expect_named(load_F1, c("diff", "mean_abs_diff", "median_abs_diff", "min_abs_diff",
                       "max_abs_diff", "max_dec", "are_equal", "diff_corres",
                       "diff_corres_cross", "g", "settings"))

  expect_type(int$diff, "integer")
  expect_type(int$mean_abs_diff, "double")
  expect_type(int$median_abs_diff, "double")
  expect_type(int$min_abs_diff, "integer")
  expect_type(int$max_abs_diff, "integer")
  expect_type(int$max_dec, "double")
  expect_type(int$are_equal, "double")
  expect_equal(int$diff_corres, NA)
  expect_equal(int$diff_corres_cross, NA)
  expect_type(int$g, "double")
  expect_type(int$settings, "list")

  expect_type(dec$diff, "double")
  expect_type(dec$mean_abs_diff, "double")
  expect_type(dec$median_abs_diff, "double")
  expect_type(dec$min_abs_diff, "double")
  expect_type(dec$max_abs_diff, "double")
  expect_type(dec$max_dec, "integer")
  expect_type(dec$are_equal, "double")
  expect_equal(dec$diff_corres, NA)
  expect_equal(dec$diff_corres_cross, NA)
  expect_type(dec$g, "double")
  expect_type(dec$settings, "list")

  checkmate::expect_matrix(matr$diff)
  expect_type(matr$mean_abs_diff, "double")
  expect_type(matr$median_abs_diff, "double")
  expect_type(matr$min_abs_diff, "double")
  expect_type(matr$max_abs_diff, "double")
  expect_type(matr$max_dec, "double")
  expect_type(matr$are_equal, "double")
  expect_type(matr$diff_corres, "integer")
  expect_type(matr$diff_corres_cross, "integer")
  expect_type(matr$g, "double")
  expect_type(matr$settings, "list")

  checkmate::expect_matrix(load$diff)
  expect_type(load$mean_abs_diff, "double")
  expect_type(load$median_abs_diff, "double")
  expect_type(load$min_abs_diff, "double")
  expect_type(load$max_abs_diff, "double")
  expect_type(load$max_dec, "integer")
  expect_type(load$are_equal, "double")
  expect_type(load$diff_corres, "integer")
  expect_type(load$diff_corres_cross, "integer")
  expect_type(load$g, "double")
  expect_type(load$settings, "list")

  checkmate::expect_matrix(load_ro1$diff)
  expect_type(load_ro1$mean_abs_diff, "double")
  expect_type(load_ro1$median_abs_diff, "double")
  expect_type(load_ro1$min_abs_diff, "double")
  expect_type(load_ro1$max_abs_diff, "double")
  expect_type(load_ro1$max_dec, "integer")
  expect_type(load_ro1$are_equal, "double")
  expect_type(load_ro1$diff_corres, "integer")
  expect_type(load_ro1$diff_corres_cross, "integer")
  expect_type(load_ro1$g, "double")
  expect_type(load_ro1$settings, "list")

  checkmate::expect_matrix(load_ro2$diff)
  expect_type(load_ro2$mean_abs_diff, "double")
  expect_type(load_ro2$median_abs_diff, "double")
  expect_type(load_ro2$min_abs_diff, "double")
  expect_type(load_ro2$max_abs_diff, "double")
  expect_type(load_ro2$max_dec, "integer")
  expect_type(load_ro2$are_equal, "double")
  expect_type(load_ro2$diff_corres, "integer")
  expect_type(load_ro2$diff_corres_cross, "integer")
  expect_type(load_ro2$g, "double")
  expect_type(load_ro2$settings, "list")

  checkmate::expect_matrix(load_F1$diff)
  expect_type(load_F1$mean_abs_diff, "double")
  expect_type(load_F1$median_abs_diff, "double")
  expect_type(load_F1$min_abs_diff, "double")
  expect_type(load_F1$max_abs_diff, "double")
  expect_type(load_F1$max_dec, "integer")
  expect_type(load_F1$are_equal, "double")
  expect_type(load_F1$diff_corres, "double")
  expect_type(load_F1$diff_corres_cross, "double")
  expect_type(load_F1$g, "double")
  expect_type(load_F1$settings, "list")

})

test_that("efa_compare returns the correct values", {
  expect_equal(int$diff, rep(0, 10))
  expect_equal(int$mean_abs_diff, 0)
  expect_equal(int$median_abs_diff, 0)
  expect_equal(int$min_abs_diff, 0)
  expect_equal(int$max_abs_diff, 0)
  expect_equal(int$max_dec, 0)
  expect_equal(int$are_equal, 0)
  expect_equal(int$g, 0)

  expect_equal(dec$diff, c(0, 0.1))
  expect_equal(dec$mean_abs_diff, 0.05)
  expect_equal(dec$median_abs_diff, 0.05)
  expect_equal(dec$min_abs_diff, 0)
  expect_equal(dec$max_abs_diff, .1)
  expect_equal(dec$max_dec, 1)
  expect_equal(dec$are_equal, 0)
  expect_equal(dec$g, 0.0707, tolerance = .01)

  expect_equal(matr$diff, matrix(c(0, 0, 0, 1), ncol = 2))
  expect_equal(matr$mean_abs_diff, 0.25)
  expect_equal(matr$median_abs_diff, 0)
  expect_equal(matr$min_abs_diff, 0)
  expect_equal(matr$max_abs_diff, 1)
  expect_equal(matr$max_dec, 0)
  # 1 vs 2 in one cell: the two differ already in the integer part, so there is no
  # decimal place they agree to -- NA, not the 0 that means "agree to zero decimals".
  expect_true(is.na(matr$are_equal))
  expect_equal(matr$g, 0.5, tolerance = .01)
  expect_equal(matr$diff_corres, 1)
  expect_equal(matr$diff_corres_cross, 0)

  expect_equal(load$diff, matrix(rep(0, 54), ncol = 3,
                                 dimnames = list(c(paste0("V", seq_len(18))),
                                                 c(paste0("F", seq_len(3))))))
  expect_equal(load$mean_abs_diff, 0)
  expect_equal(load$median_abs_diff, 0)
  expect_equal(load$min_abs_diff, 0)
  expect_equal(load$max_abs_diff, 0)
  expect_equal(load$max_dec, 17, tolerance = 1)
  expect_equal(load$are_equal, 17, tolerance = 1)
  expect_equal(load$g, 0)
  expect_equal(load$diff_corres, 0)
  expect_equal(load$diff_corres_cross, 0)

  expect_equal(load_ro1$mean_abs_diff, 0)
  expect_equal(load_ro1$median_abs_diff, 0)
  expect_equal(load_ro1$min_abs_diff, 0)
  expect_equal(load_ro1$max_abs_diff, 0)
  expect_equal(load_ro1$max_dec, 17, tolerance = 1)
  expect_equal(load_ro1$are_equal, 17, tolerance = 1)
  expect_equal(load_ro1$g, 0)
  expect_equal(load_ro1$diff_corres, 0)
  expect_equal(load_ro1$diff_corres_cross, 0)

  expect_equal(load_ro2$mean_abs_diff, 0)
  expect_equal(load_ro2$median_abs_diff, 0)
  expect_equal(load_ro2$min_abs_diff, 0)
  expect_equal(load_ro2$max_abs_diff, 0)
  expect_equal(load_ro2$max_dec, 17, tolerance = 1)
  expect_equal(load_ro2$are_equal, 17, tolerance = 1)
  expect_equal(load_ro2$g, 0)
  expect_equal(load_ro2$diff_corres, 0)
  expect_equal(load_ro2$diff_corres_cross, 0)

  expect_equal(round(load_F1$diff, 4), matrix(rep(0, 18), ncol = 1,
                                              dimnames = list(c(paste0("V",
                                                                       seq_len(18))),
                                                              "F1")))
  expect_equal(round(load_F1$mean_abs_diff, 4), 0)
  expect_equal(round(load_F1$median_abs_diff, 4), 0)
  expect_equal(round(load_F1$min_abs_diff, 4), 0)
  expect_equal(round(load_F1$max_abs_diff, 4), 0)
  expect_equal(load_F1$max_dec, 15)
  expect_equal(load_F1$are_equal, 2)
  expect_equal(round(load_F1$g, 4), 0)
  expect_equal(load_F1$diff_corres, 0)
  expect_equal(load_F1$diff_corres_cross, 0)

  expect_equal(efa_compare(vec_s, vec_s[c(3, 1, 2)],
                       reorder = "names")$mean_abs_diff, 0)
  expect_equal(efa_compare(mat_s, mat_s[, c(2, 1)],
                       reorder = "names")$mean_abs_diff, 0)
  expect_equal(efa_compare(psych_PAF$unrot_loadings,
                       psych_PAF$unrot_loadings[, c(3, 1, 2)])$mean_abs_diff, 0)
})


test_that("are_equal counts decimal agreement without floating-point miscounts", {
  # 0.57 is held as 0.5699999999999999; a plain trunc(x * 100) reads its second
  # decimal as 6 and would wrongly report agreement with 0.5699999 to two
  # decimals. They agree only to the first decimal.
  expect_equal(efa_compare(c(0.57, 0.1), c(0.5699999, 0.1),
                       reorder = "none")$are_equal, 1)
  # 0.6285 * 1e4 underflows to 6284.999...; a plain trunc drops a place and would
  # under-report. 0.6285 and 0.62851 agree to four decimals.
  expect_equal(efa_compare(c(0.6285, 0.1), c(0.62851, 0.1),
                       reorder = "none")$are_equal, 4)
})


test_that("are_equal separates disagreeing integer parts from zero-decimal agreement", {
  # Both pairs disagree from the first decimal on, but only the first pair already
  # disagrees in its integer part; are_equal must not report the same value for the two.
  expect_true(is.na(efa_compare(c(1.5, 2.5, 3.5), c(9.5, 8.5, 7.5),
                                reorder = "none")$are_equal))
  expect_equal(efa_compare(c(0.5, 0.5, 0.5), c(0.9, 0.9, 0.9),
                           reorder = "none")$are_equal, 0)
})


test_that("congruence reordering yields a true permutation when columns collide", {
  # x has orthonormal factors; y is built so that the greedy row-wise which.max
  # match assigns y-column 1 to BOTH x-factor 1 and x-factor 2 (a collision that
  # would duplicate one y column and drop another). The optimal one-to-one
  # assignment is the permutation c(1, 3, 2).
  x <- diag(3)
  y <- matrix(c(0.7, 0.7, 0.10,    # aligned with factors 1 and 2
                0.2, 0.1, 0.97,    # aligned with factor 3
                0.1, 0.2, 0.95),   # aligned with factor 3
              ncol = 3)

  cmp <- efa_compare(x, y, reorder = "congruence", corres = FALSE)

  # y as reordered (and sign-matched) inside efa_compare
  y_reordered <- x - cmp$diff

  # index in y that each reordered column came from
  recovered <- vapply(seq_len(ncol(y_reordered)), function(j) {
    which(vapply(seq_len(ncol(y)), function(k)
      isTRUE(all.equal(y_reordered[, j], y[, k])), logical(1)))
  }, integer(1))

  # a genuine permutation: every factor used exactly once, never duplicated or
  # dropped (the greedy which.max produced c(1, 1, 2))
  expect_setequal(recovered, seq_len(ncol(y)))
  expect_equal(recovered, c(1L, 3L, 2L))
})


test_that("errors etc. are thrown correctly", {
  expect_error(efa_compare(c(1, 2), 1), class = "efa_compare_dim_mismatch")
  expect_error(efa_compare(c(1, 2), c("1", "2")), class = "efa_compare_bad_input")
  expect_error(efa_compare(c(1, 2), data.frame(x = "1", y = "2")), class = "efa_compare_bad_input")

  expect_error(efa_compare(matrix(c(0, 0, 0, 1), ncol = 2),
                       matrix(c(0, 0, 0, 1), ncol = 1)), class = "efa_compare_dim_mismatch")

  expect_warning(efa_compare(vec_s, vec_s), class = "efa_compare_reorder_vectors")
  expect_warning(efa_compare(vec_s, vec_L, reorder = "names"), class = "efa_compare_reorder_mismatch")
  expect_warning(efa_compare(vec_s, 1:3, reorder = "names"), class = "efa_compare_reorder_unnamed")
  expect_warning(efa_compare(1:3, 1:3, reorder = "names"), class = "efa_compare_reorder_unnamed")

  expect_warning(efa_compare(matrix(c(0, 0, 0, 1), ncol = 2),
                         matrix(c(0, 0, 0, 1), ncol = 2),
                         reorder = "names"), class = "efa_compare_reorder_unnamed")
  expect_warning(efa_compare(mat_s, mat_L, reorder = "names"), class = "efa_compare_reorder_mismatch")

  expect_error(efa_compare(mat_s, mat_s), class = "efa_compare_congruence_na")

  # A column whose norm is below sqrt(.Machine$double.eps) but non-zero makes
  # Tucker's congruence undefined just as an exactly-zero column does; it must
  # surface the reorder error, not a low-level zero-column abort.
  nz <- matrix(c(1e-9, 1e-9, 0.8, 0.6), ncol = 2)
  expect_error(efa_compare(nz, nz, reorder = "congruence", corres = FALSE),
               class = "efa_compare_congruence_na")

  # Non-finite inputs, and finite values whose squared magnitudes overflow,
  # both leave Tucker's congruence undefined; each must surface the same reorder
  # error rather than a non-finite/solver abort from the alignment internals.
  m_inf <- matrix(c(Inf, 0.5, 0.5, 0.8, 0.6, 0.7), ncol = 2)
  expect_error(efa_compare(m_inf, m_inf, reorder = "congruence", corres = FALSE),
               class = "efa_compare_congruence_na")
  m_big <- matrix(c(1e200, 1e200, 1e200, 0.8, 0.6, 0.7), ncol = 2)
  expect_error(efa_compare(m_big, m_big, reorder = "congruence", corres = FALSE),
               class = "efa_compare_congruence_na")

})

test_that("efa_compare is NA-safe and honours na.rm (vector and matrix)", {
  # NA-containing input used to crash in .decimals ("missing value where
  # TRUE/FALSE needed") before efa_compare could return; the documented na.rm
  # argument must work for both vector and matrix input (reorder = "none"
  # avoids the congruence path, which aborts earlier on NAs).
  v_na <- c(1, 2, NA, 4)
  v_ok <- c(1, 2, 3, 4)
  expect_s3_class(efa_compare(v_na, v_ok, reorder = "none", na.rm = TRUE), "efa_compare")
  expect_s3_class(efa_compare(v_na, v_ok, reorder = "none", na.rm = FALSE), "efa_compare")

  m_na <- matrix(c(1, NA, 1, 2), ncol = 2)
  m_ok <- matrix(c(1, 1, 1, 1), ncol = 2)
  expect_s3_class(efa_compare(m_na, m_ok, reorder = "none", na.rm = TRUE), "efa_compare")
  expect_s3_class(efa_compare(m_na, m_ok, reorder = "none", na.rm = FALSE), "efa_compare")

  # na.rm = TRUE drops the missing entry from the RMSE statistic rather than
  # propagating NA into it.
  cmp_rm <- efa_compare(v_na, v_ok, reorder = "none", na.rm = TRUE)
  expect_false(is.na(cmp_rm$g))
  expect_equal(cmp_rm$g, 0)

  # are_equal mirrors the other statistics under missingness: NA (printed "none")
  # when na.rm = FALSE and an element is missing, and a number when na.rm = TRUE
  # (consistent with max_dec, which always drops missings).
  expect_true(is.na(efa_compare(c(0.5, NA), c(0.5, 0.5),
                            reorder = "none", na.rm = FALSE)$are_equal))
  expect_equal(efa_compare(c(0.5, NA), c(0.5, 0.5),
                       reorder = "none", na.rm = TRUE)$are_equal, 1)

  # Even with na.rm = TRUE, if every pair has a missing value there is nothing to
  # compare, so are_equal is NA rather than spuriously climbing to max_dec (the
  # empty all() would otherwise report agreement to every decimal place). The
  # all-missing overlap also makes min()/max() warn ("no non-missing arguments"),
  # which is unrelated to are_equal, so it is suppressed here.
  expect_true(is.na(suppressWarnings(
    efa_compare(c(1.234, NA), c(NA, 5.678), reorder = "none", na.rm = TRUE))$are_equal))

  # Near-zero values render in scientific notation, which used to crash the
  # decimal-count helper that drives max_dec; efa_compare must still return.
  expect_s3_class(efa_compare(c(0.5, 0.00003), c(0.5, 0.4), reorder = "none"),
                  "efa_compare")
})

test_that("print output is stable", {
  local_reproducible_output()

  # matrix difference (variables x factors). matr/int hold exact 0/1 differences, so the
  # snapshots are recorded literally (no scrub) to pin the column alignment.
  expect_snapshot(print(matr))

  # vector difference (rendered as a single column, one value per row)
  expect_snapshot(print(int))
})

test_that("display settings can be overridden at print time", {
  local_reproducible_output()

  # The recorded settings are unchanged; only this call's report differs.
  plain <- cli::ansi_strip(format(matr))
  two <- cli::ansi_strip(format(matr, digits = 2))
  expect_match(paste(plain, collapse = " "), "0.2500", fixed = TRUE)
  expect_match(paste(two, collapse = " "), "0.25", fixed = TRUE)
  expect_no_match(paste(two, collapse = " "), "0.2500", fixed = TRUE)
  expect_equal(matr$settings$digits, 4)

  # print() forwards its dots to format(), so the same override works there.
  expect_identical(cli::ansi_strip(utils::capture.output(print(matr, digits = 2))), two)

  # print_diff drops the difference table.
  expect_true(any(grepl("Elementwise differences", plain, fixed = TRUE)))
  expect_false(any(grepl("Elementwise differences",
                         cli::ansi_strip(format(matr, print_diff = FALSE)), fixed = TRUE)))

  # The three colouring thresholds change only the green/red styling, so they are
  # observable solely with colours on: each must flip the ANSI of its own line while
  # leaving the reported numbers untouched. matr's mean and max are .25 and 1, so a
  # threshold on either side of those is decisive; round_red is read against
  # are_equal, which is NA for matr (always red) and 0 for int.
  old <- options(cli.num_colors = 256)
  on.exit(options(old), add = TRUE)
  expect_false(identical(format(matr, m_red = 1), format(matr, m_red = 0)))
  expect_false(identical(format(matr, range_red = 2), format(matr, range_red = 0)))
  expect_false(identical(format(int, round_red = 0), format(int, round_red = 9)))
  # ... and none of them alters the numbers themselves.
  strip_matr <- function(...) cli::ansi_strip(format(matr, ...))
  expect_identical(strip_matr(m_red = 1), strip_matr(m_red = 0))
  expect_identical(strip_matr(range_red = 2), strip_matr(range_red = 0))
  options(cli.num_colors = 1)

  # An override is validated exactly as efa_compare() validates the argument.
  expect_error(format(matr, digits = -1))

  # omitting an argument leaves the recorded value in force
  expect_identical(format(matr), format(matr, digits = matr$settings$digits))
})

test_that("the minimum-decimals line is dropped when it carries no information", {
  local_reproducible_output()

  # Two ordinary doubles carry full double precision, so the count says nothing.
  full <- efa_compare(c(1/3, 2/7, 1/6), c(1/3, 2/7, 1/6 + 1e-4), reorder = "none")
  expect_gte(full$max_dec, 15)
  expect_no_match(paste(cli::ansi_strip(format(full)), collapse = " "),
                  "Minimum number of decimals provided", fixed = TRUE)

  # A rounded input bounds the comparison, so the line is informative and shown.
  expect_match(paste(cli::ansi_strip(format(dec)), collapse = " "),
               "Minimum number of decimals provided", fixed = TRUE)
})

test_that("plot returns a ggplot and guards too-few differences", {
  # Smoke-test only (no vdiffr baseline) because geom_jitter draws random positions,
  # so the rendered plot is not reproducible across runs.
  p <- plot(matr)
  expect_s3_class(p, "ggplot")

  # plot_red is a drawing setting, overridable per call without recomputing: it
  # decides which differences are flagged and is reported in the subtitle.
  expect_true("large difference" %in% p$data$color)
  p2 <- plot(matr, plot_red = 2)
  expect_setequal(p2$data$color, "acceptable difference")
  expect_match(p2$labels$subtitle, "difference coloring: 2", fixed = TRUE)
  expect_equal(matr$settings$plot_red, 0.01)

  # too few differences to plot
  expect_error(plot(dec), class = "efa_compare_too_few_to_plot")
})


# Reference implementation of the statistics and factor correspondences that
# .compare_loadings() returns. The test below pins the extracted core against
# accidental drift: it mirrors the helper's computation (a regression guard) rather
# than deriving it independently, so the statistic *values* are anchored by the
# value-based efa_compare tests above. Only the vector NA path is written differently
# here (an explicit numeric/integer test, not a catch-all else); the two agree for
# the matrix/numeric/integer inputs exercised.
ref_compare_loadings <- function(x, y, thresh = 0.3, na.rm = FALSE, corres = TRUE) {

  if (inherits(x, "matrix")) {
    if (ncol(x) > 1 && isTRUE(corres)) {
      corres_list <- .factor_corres(x, y, thresh = thresh)
      diff_corres <- corres_list$diff_corres
      diff_corres_cross <- corres_list$diff_corres_cross
    } else {
      diff_corres <- 0
      diff_corres_cross <- 0
    }
  } else if (inherits(x, c("numeric", "integer"))) {
    diff_corres <- NA
    diff_corres_cross <- NA
  }

  diff <- x - y

  sq <- diff ^ 2
  n_ok <- if (na.rm) sum(!is.na(diff)) else length(diff)
  g <- sqrt(sum(sq, na.rm = na.rm) / n_ok)

  mean_abs_diff <- mean(abs(diff), na.rm = na.rm)
  median_abs_diff <- stats::median(abs(diff), na.rm = na.rm)

  min_abs_diff <- min(abs(diff), na.rm = na.rm)
  max_abs_diff <- max(abs(diff), na.rm = na.rm)

  max_dec <- min(c(.decimals(x), .decimals(y)))

  ax <- abs(x)
  ay <- abs(y)
  if (anyNA(diff) && (!na.rm || all(is.na(diff)))) {
    are_equal <- NA_real_
  } else {
    are_equal <- NA_real_
    for (d in 0:max_dec) {
      if (!isTRUE(all(trunc(signif(ax * 10^d, 13)) == trunc(signif(ay * 10^d, 13)),
                      na.rm = na.rm))) {
        break
      }
      are_equal <- as.double(d)
    }
  }

  list(
    diff = diff,
    mean_abs_diff = mean_abs_diff,
    median_abs_diff = median_abs_diff,
    min_abs_diff = min_abs_diff,
    max_abs_diff = max_abs_diff,
    max_dec = max_dec,
    are_equal = are_equal,
    diff_corres = diff_corres,
    diff_corres_cross = diff_corres_cross,
    g = g
  )
}

test_that(".compare_loadings reproduces efa_compare's pre-extraction computation", {
  set_dn <- function(m) {
    colnames(m) <- paste0("F", seq_len(ncol(m)))
    rownames(m) <- paste0("V", seq_len(nrow(m)))
    m
  }

  m2x <- set_dn(matrix(c(0.8, 0.1, 0.2, 0.7, 0.15, 0.75, 0.1, 0.05), ncol = 2))
  m2y <- set_dn(matrix(c(0.75, 0.12, 0.25, 0.68, 0.2, 0.7, 0.08, 0.02), ncol = 2))
  m1x <- set_dn(matrix(c(0.5, 0.6, 0.7, 0.8), ncol = 1))
  m1y <- set_dn(matrix(c(0.52, 0.58, 0.72, 0.79), ncol = 1))
  # signed disagreement (0.57 vs -0.5699999) and float-noise magnitudes exercise
  # are_equal / max_dec.
  vx  <- c(a = 0.57, b = 0.6285, c = 0.1)
  vy  <- c(a = -0.5699999, b = 0.62851, c = 0.1)
  vix <- 1:6
  viy <- c(1L, 2L, 3L, 4L, 5L, 7L)
  # NA-containing matrix, corres kept on to mirror efa_compare's default path.
  mnax <- set_dn(matrix(c(0.8, NA, 0.2, 0.7, 0.15, 0.75, NA, 0.05), ncol = 2))
  mnay <- set_dn(matrix(c(0.75, 0.12, 0.25, 0.68, 0.2, 0.7, 0.08, 0.02), ncol = 2))
  vnax <- c(0.5, NA, 0.30003, 0.4)
  vnay <- c(0.5, 0.5, 0.3, 0.4)
  # near-zero value renders in scientific notation (drives the .decimals path).
  scix <- c(0.5, 0.00003, 0.4)
  sciy <- c(0.5, 0.4, 0.4)
  # every pair missing under na.rm = TRUE -> are_equal NA, min()/max() warn.
  allmiss_x <- c(1.234, NA)
  allmiss_y <- c(NA, 5.678)

  cases <- list(
    list(x = m2x, y = m2y, corres = TRUE),
    list(x = m2x, y = m2y, corres = FALSE),
    list(x = m2x, y = m2y, corres = TRUE, thresh = 0.5),
    list(x = m1x, y = m1y),
    list(x = vx, y = vy),
    list(x = vix, y = viy),
    list(x = mnax, y = mnay, na.rm = FALSE),
    list(x = mnax, y = mnay, na.rm = TRUE),
    list(x = vnax, y = vnay, na.rm = FALSE),
    list(x = vnax, y = vnay, na.rm = TRUE),
    list(x = scix, y = sciy),
    list(x = allmiss_x, y = allmiss_y, na.rm = TRUE)
  )

  for (cs in cases) {
    got <- suppressWarnings(do.call(.compare_loadings, cs))
    ref <- suppressWarnings(do.call(ref_compare_loadings, cs))
    expect_identical(got, ref)
  }
})

rm(int, dec, matr, SPSS_PAF, psych_PAF, load, load_ro1, load_ro2, SPSS_PAF_1,
   psych_PAF_1, load_F1, vec_s, vec_L, mat_s, mat_L)
