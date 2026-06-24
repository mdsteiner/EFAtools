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

# A non-loading "corr"-role table (variances accounted for) with enough factor columns to
# overflow an 80-column console; values are deterministic so the printed decimals are stable.
make_wide_vars <- function(n = 12) {
  prop <- seq(0.30, 0.02, length.out = n)
  va <- rbind("SS loadings" = seq(3.5, 0.4, length.out = n),
              "Prop Var"    = prop,
              "Cum Var"     = cumsum(prop))
  colnames(va) <- paste0("F", seq_len(n))
  va
}

# A symmetric factor-intercorrelation matrix (Phi) wide enough to overflow an 80-column
# console; only the lower triangle (plus the unit diagonal) is ever shown.
make_wide_phi <- function(n = 12) {
  phi <- matrix(0, n, n)
  phi[lower.tri(phi)] <- seq(0.10, 0.65, length.out = n * (n - 1) / 2)
  phi <- phi + t(phi)
  diag(phi) <- 1
  dimnames(phi) <- list(paste0("F", seq_len(n)), paste0("F", seq_len(n)))
  phi
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

test_that("color controls print styling while format() stays plain", {
  withr::local_options(cli.num_colors = 256)
  L <- make_loadings()
  # styling is controlled by `color` in print(), which embeds ANSI when colours are on
  expect_true(cli::ansi_has_any(paste(
    utils::capture.output(print(L, color = TRUE)), collapse = "")))
  expect_false(cli::ansi_has_any(paste(
    utils::capture.output(print(L, color = FALSE)), collapse = "")))
  # format() is the plain-text representation regardless of `color`
  expect_false(cli::ansi_has_any(paste(format(L, color = TRUE), collapse = "")))
  expect_false(cli::ansi_has_any(paste(format(L, color = FALSE), collapse = "")))
  # format() returns only the rendered table lines, with no trailing blank element
  fmt <- format(L)
  expect_true(nzchar(fmt[length(fmt)]))
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

test_that("wide non-loading tables wrap into stacked blocks", {
  local_reproducible_output() # console width 80
  # a rectangular variances-accounted table and a symmetric Phi (lower triangle only)
  expect_snapshot(.print_efa_matrix(make_wide_vars(), role = "corr"))
  expect_snapshot(.print_efa_matrix(make_wide_phi(), role = "corr",
                                    lower_only = TRUE))
})

test_that("non-loading block wrapping splits columns without folding rows", {
  local_reproducible_output() # console width 80

  corr_lines <- function(values, lower_only = FALSE) {
    .efa_format_matrix(values,
                       row_labels = rownames(values),
                       col_labels = colnames(values),
                       col_roles = rep("corr", ncol(values)),
                       lower_only = lower_only)
  }

  vars <- corr_lines(make_wide_vars())
  # more than one block, and no emitted line exceeds the console width
  expect_gt(sum(grepl("(block ", vars, fixed = TRUE)), 1L)
  expect_true(all(cli::ansi_nchar(vars) <= 80L))

  phi <- corr_lines(make_wide_phi(), lower_only = TRUE)
  expect_gt(sum(grepl("(block ", phi, fixed = TRUE)), 1L)
  expect_true(all(cli::ansi_nchar(phi) <= 80L))
  # a later block's leading rows fall entirely in the blanked upper triangle and are dropped,
  # so no emitted line is a value-less row label (a phantom blank row would print as bare "Fk")
  expect_false(any(grepl("^F[0-9]+ *$", phi)))
})

test_that("format.SLLOADINGS flags a Heywood case", {
  local_reproducible_output()
  expect_snapshot(print(make_sl(heywood = TRUE)))
})

test_that("format.SLLOADINGS prints cleanly without Heywood cases", {
  local_reproducible_output()
  expect_snapshot(print(make_sl(heywood = FALSE)))
})

test_that("format.SLLOADINGS counts a communality-only Heywood item once", {
  # h2 > 1 (and u2 < 0) with no single loading > 1 must be flagged exactly once:
  # not missed (as a loading-only count would), nor multiplied across coupled cells.
  sl <- matrix(c(0.70, 0.20, 0.10, 0.53, 0.47,
                 0.80, 0.75, 0.10, 1.20, -0.20,
                 0.40, 0.10, 0.55, 0.46, 0.54),
               nrow = 3, byrow = TRUE,
               dimnames = list(c("i1", "i2", "i3"), c("g", "F1", "F2", "h2", "u2")))
  class(sl) <- "SLLOADINGS"
  out <- cli::ansi_strip(utils::capture.output(print(sl)))
  expect_true(any(grepl("Results contain a Heywood case!", out, fixed = TRUE)))
  expect_false(any(grepl("Heywood cases", out, fixed = TRUE)))
})

test_that("format.SLLOADINGS returns plain text even when styling is enabled", {
  old <- options(cli.num_colors = 256)
  on.exit(options(old), add = TRUE)

  styled <- utils::capture.output(print(make_sl(heywood = TRUE)))
  plain <- format(make_sl(heywood = TRUE))

  expect_true(any(grepl("\033", styled, fixed = TRUE)))
  expect_false(any(grepl("\033", plain, fixed = TRUE)))
  # no trailing blank element
  expect_true(nzchar(plain[length(plain)]))
})

test_that("print.LOADINGS argument validators raise classed conditions", {
  L <- make_loadings()

  expect_error(print(L, cutoff = -1), class = "efa_print_invalid_cutoff")
  expect_error(print(L, digits = 1.5), class = "efa_print_invalid_digits")
  expect_error(print(L, max_name_length = 0),
               class = "efa_print_invalid_max_name_length")
  expect_error(print(L, max_factor_name_length = 0),
               class = "efa_print_invalid_max_factor_name_length")
  expect_error(print(L, max_factors_per_block = 0),
               class = "efa_print_invalid_max_factors_per_block")
  expect_error(print(L, color = NA), class = "efa_print_invalid_color")
  expect_error(print(L, legend = NA), class = "efa_print_invalid_legend")
  # h2 with the wrong length
  expect_error(print(L, h2 = c(0.5, 0.5)), class = "efa_print_invalid_h2")

  # x that is not a valid numeric matrix is rejected by the validator
  valid_args <- list(cutoff = .3, digits = 3, max_name_length = 10, h2 = NULL,
                     color = TRUE, max_factor_name_length = NULL,
                     max_factors_per_block = NULL, legend = FALSE)
  expect_error(do.call(.validate_loadings_print_args,
                       c(list(x = "not a matrix"), valid_args)),
               class = "efa_print_invalid_x")
  expect_error(do.call(.validate_loadings_print_args,
                       c(list(x = matrix(numeric(0), 0, 0)), valid_args)),
               class = "efa_print_invalid_x")
})

test_that("a named h2 that omits a row name raises a classed condition", {
  # make_loadings() has row names; a named h2 missing one of them must abort in
  # .align_loadings_h2 rather than silently mismatching communalities to rows.
  L <- make_loadings()
  expect_error(
    print(L, h2 = c(fun = .7, friends_long_name = .6, enjoy = .6, wrong = .5)),
    class = "efa_print_invalid_h2"
  )
})
