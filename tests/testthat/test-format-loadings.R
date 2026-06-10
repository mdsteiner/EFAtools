# Snapshot and structural tests for the shared loading-matrix renderer
# (.efa_format_matrix) behind format.LOADINGS / format.SLLOADINGS. The inputs are literal
# matrices with hand-chosen values (including negatives, a Heywood loading > 1, h2 > 1 and
# u2 < 0), so the printed decimals are platform-stable and the snapshots pin down vertical
# alignment, styling structure, the legend, and the block-wise wrapping directly. Snapshots
# are recorded under local_reproducible_output() (plain text, no ANSI).

make_loadings <- function() {
  L <- matrix(c( 0.82, -0.11,  0.05,
                 0.45,  0.60, -0.02,
                -0.30,  0.71,  0.123,
                 0.04, -0.05,  1.08),
              nrow = 4, byrow = TRUE,
              dimnames = list(c("fun", "friends_long_name", "enjoy", "hurt"),
                              c("F1", "F2", "F3")))
  class(L) <- "LOADINGS"
  L
}

make_wide_loadings <- function() {
  W <- matrix(round(seq(-0.9, 0.9, length.out = 4 * 8), 2), nrow = 4,
              dimnames = list(c("v1", "v2", "v3", "v4"), paste0("F", 1:8)))
  class(W) <- "LOADINGS"
  W
}

make_sl <- function(heywood = TRUE) {
  f1 <- if (heywood) 1.05 else 0.62
  h2_2 <- if (heywood) 1.10 else 0.70
  u2_2 <- if (heywood) -0.10 else 0.30
  sl <- matrix(c(0.70,  0.20,  0.10, 0.53, 0.47,
                 0.65,  f1,   -0.05, h2_2, u2_2,
                 0.40,  0.10,  0.55, 0.46, 0.54),
               nrow = 3, byrow = TRUE,
               dimnames = list(c("i1", "i2", "i3"),
                               c("g", "F1", "F2", "h2", "u2")))
  class(sl) <- "SLLOADINGS"
  sl
}

test_that("format.LOADINGS aligns decimals and renders a plain table", {
  local_reproducible_output()
  expect_snapshot(print(make_loadings()))
})

test_that("format.LOADINGS prints communalities and a legend", {
  local_reproducible_output()
  h2 <- c(0.70, 0.58, 0.63, 1.18)
  expect_snapshot(print(make_loadings(), h2 = h2, legend = TRUE))
})

test_that("format.LOADINGS sorts rows when requested", {
  local_reproducible_output()
  expect_snapshot(print(make_loadings(), sort_loadings = "clustered"))
})

test_that("color = FALSE never emits ANSI even when colours are on", {
  withr::local_options(cli.num_colors = 256)
  L <- make_loadings()
  expect_true(cli::ansi_has_any(paste(format(L, color = TRUE), collapse = "")))
  expect_false(cli::ansi_has_any(paste(format(L, color = FALSE), collapse = "")))
})

test_that("wide matrices wrap into stacked blocks via max_factors_per_block", {
  local_reproducible_output()
  expect_snapshot(print(make_wide_loadings(), max_factors_per_block = 3))
})

test_that("wide matrices wrap into stacked blocks on a narrow console", {
  local_reproducible_output()
  withr::local_options(cli.width = 40)
  expect_snapshot(print(make_wide_loadings()))
})

test_that("block wrapping stacks blocks and never folds a row", {
  local_reproducible_output()
  withr::local_options(cli.width = 40)
  out <- format(make_wide_loadings())
  # more than one block: each block carries its own "(block i/n)" label
  expect_gt(sum(grepl("(block ", out, fixed = TRUE)), 1L)
  # no emitted line exceeds the console width -> columns were split, lines not wrapped
  expect_true(all(cli::ansi_nchar(out) <= 40L))
})

test_that("format.SLLOADINGS flags a Heywood case", {
  local_reproducible_output()
  expect_snapshot(print(make_sl(heywood = TRUE)))
})

test_that("format.SLLOADINGS prints cleanly without Heywood cases", {
  local_reproducible_output()
  expect_snapshot(print(make_sl(heywood = FALSE)))
})
