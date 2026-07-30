# Snapshot and structural tests for the shared loading-matrix renderer
# (.efa_format_matrix) behind format.efa_loadings / format.efa_sl_loadings. The inputs are literal
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
  class(L) <- c("efa_loadings", "LOADINGS")
  L
}

make_wide_loadings <- function() {
  W <- matrix(round(seq(-0.9, 0.9, length.out = 4 * 8), 2), nrow = 4,
              dimnames = list(c("v1", "v2", "v3", "v4"), paste0("F", 1:8)))
  class(W) <- c("efa_loadings", "LOADINGS")
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
  class(sl) <- c("efa_sl_loadings", "SLLOADINGS")
  sl
}

test_that("format.efa_loadings aligns decimals and renders a plain table", {
  local_reproducible_output()
  expect_snapshot(print(make_loadings()))
})

test_that("format.efa_loadings prints communalities", {
  local_reproducible_output()
  # local_reproducible_output() turns colours off, so the requested legend is omitted:
  # it only describes styling, and none is rendered here (see the legend test below).
  h2 <- c(0.70, 0.58, 0.63, 1.18)
  expect_snapshot(print(make_loadings(), h2 = h2, legend = TRUE))
})

test_that("format.efa_loadings sorts rows when requested", {
  local_reproducible_output()
  expect_snapshot(print(make_loadings(), sort_loadings = "clustered"))
})

test_that("format.efa_loadings is the source of truth and honours the colour state", {
  L <- make_loadings()

  # print() is cat(format(x), sep = "\n") plus one blank line for console spacing
  expect_identical(utils::capture.output(print(L)), c(format(L), ""))

  withr::local_options(cli.num_colors = 256)
  # styling is controlled by `color`, in format() exactly as in print()
  expect_true(cli::ansi_has_any(paste(
    utils::capture.output(print(L, color = TRUE)), collapse = "")))
  expect_false(cli::ansi_has_any(paste(
    utils::capture.output(print(L, color = FALSE)), collapse = "")))
  expect_true(cli::ansi_has_any(paste(format(L, color = TRUE), collapse = "")))
  expect_false(cli::ansi_has_any(paste(format(L, color = FALSE), collapse = "")))

  # with colours off the table is plain whatever the caller asked for
  withr::local_options(cli.num_colors = 1)
  expect_false(cli::ansi_has_any(paste(format(L, color = TRUE), collapse = "")))

  # format() returns only the rendered table lines, with no trailing blank element
  fmt <- format(L)
  expect_true(nzchar(fmt[length(fmt)]))
})

test_that("the loading legend is printed only when its styling is rendered", {
  L <- make_loadings()
  legend_shown <- function(...) {
    any(grepl("Legend:", utils::capture.output(print(L, legend = TRUE, ...)),
              fixed = TRUE))
  }

  # The legend names bold/grey marks, so it is worth printing only when those marks are
  # actually rendered. That needs both the caller's `color` and a colour-capable target:
  # cli emits no escapes at one colour, where a legend would describe invisible styling.
  withr::with_options(list(cli.num_colors = 256), {
    expect_true(legend_shown(color = TRUE))
    expect_false(legend_shown(color = FALSE))
  })
  withr::with_options(list(cli.num_colors = 1), {
    expect_false(legend_shown(color = TRUE))
    expect_false(legend_shown(color = FALSE))
  })

  # Legend presence tracks ANSI presence exactly.
  withr::with_options(list(cli.num_colors = 256), {
    for (col in c(TRUE, FALSE)) {
      styled <- cli::ansi_has_any(paste(
        utils::capture.output(print(L, legend = TRUE, color = col)), collapse = ""))
      expect_equal(legend_shown(color = col), styled)
    }
  })

  # Suppressing the legend must not leave the spacer blank line behind: with no styling to
  # describe, `legend = TRUE` renders exactly as `legend = FALSE`.
  withr::with_options(list(cli.num_colors = 1), {
    expect_identical(utils::capture.output(print(L, legend = TRUE)),
                     utils::capture.output(print(L, legend = FALSE)))
  })

  # `format()` follows the same rule, since it is what `print()` renders: the legend
  # travels with the marks it names, and is dropped with them when colours are off.
  withr::with_options(list(cli.num_colors = 256), {
    out <- format(L, legend = TRUE)
    expect_true(any(grepl("Legend:", out, fixed = TRUE)))
    expect_true(cli::ansi_has_any(paste(out, collapse = "")))
  })
  withr::with_options(list(cli.num_colors = 1), {
    out <- format(L, legend = TRUE)
    expect_false(any(grepl("Legend:", out, fixed = TRUE)))
    expect_identical(out, format(L, legend = FALSE))
  })
})


test_that("the loading legend names the marks it describes", {
  local_reproducible_output()
  L <- make_loadings()
  withr::local_options(cli.num_colors = 8L)

  # The legend's wording is the contract; the exact escape sequences are cli's business,
  # so the assertion is on the stripped text. Nothing else pins these lines -- the printed
  # snapshot is taken with colours off, where the legend is correctly absent.
  out <- cli::ansi_strip(utils::capture.output(
    print(L, h2 = c(0.70, 0.58, 0.63, 1.18), legend = TRUE, cutoff = .3)))
  i <- which(out == "Legend:")
  expect_length(i, 1L)
  expect_false(nzchar(out[i - 1L]))                       # the spacer above it
  expect_identical(out[i + 1:3],
                   c("  bold = |loading| >= .300",
                     "  grey = below cutoff",
                     "  red h2/u2 = Heywood-relevant value"))
})


test_that("the efa report carries the loading legend through its nested render", {
  local_reproducible_output()

  # print.efa()/summary.efa() render the loading table through a format.efa_loadings()
  # call nested in their own cli::cli_format_method(), i.e. one cli_format_method inside
  # another. A sink makes cli report a single colour, so the legend would
  # be lost -- except that cli_format_method reads the colour count *before* installing its
  # sink and pins it for everything nested inside. Pinning `cli.num_colors` here would not
  # test that: it is the first thing cli consults, so it short-circuits the sink question
  # entirely and the assertion would hold with or without the outer cli_format_method.
  # `cli.default_num_colors` is read on the same side of the sink check as a real terminal's
  # isatty() detection, so it is the faithful stand-in. crayon.enabled and the environment
  # variables have to be cleared because each of them is consulted first.
  withr::local_envvar(c(NO_COLOR = NA, R_CLI_NUM_COLORS = NA))
  withr::local_options(cli.num_colors = NULL, crayon.enabled = NULL,
                       crayon.colors = NULL, knitr.in.progress = NULL,
                       cli.default_num_colors = 256L)
  skip_if_not(cli::num_ansi_colors() > 1L, "no colour-capable stand-in available")

  fit <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
                 estimator = "PAF", rotation = "promax")
  has_legend <- function(lines) any(grepl("Legend:", lines, fixed = TRUE))

  expect_true(has_legend(format(fit)))
  expect_true(has_legend(format(summary(fit))))
  expect_false(has_legend(format(fit, show_loading_legend = FALSE)))

  # The same table rendered on its own reaches the same verdict as it does nested in the
  # report: the legend is a property of the console theme, not of the call depth.
  expect_true(has_legend(
    format(fit$rot_loadings, h2 = fit$h2, legend = TRUE)))

  # The guard that makes the assertions above meaningful: nothing renders a legend once
  # the target cannot show the marks it names, whether the table is nested or standalone.
  withr::local_options(cli.num_colors = 1)
  expect_false(has_legend(format(fit)))
  expect_false(has_legend(
    format(fit$rot_loadings, h2 = fit$h2, legend = TRUE)))
})


test_that(".efa_styling_visible tracks both the caller and the target", {
  withr::with_options(list(cli.num_colors = 256), {
    expect_true(.efa_styling_visible(TRUE))
    expect_false(.efa_styling_visible(FALSE))
  })
  withr::with_options(list(cli.num_colors = 1), {
    expect_false(.efa_styling_visible(TRUE))
    expect_false(.efa_styling_visible(FALSE))
  })
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
  # a rectangular variances-accounted table and a symmetric Phi (lower triangle only),
  # through .efa_corr_lines(): the wrapper the reports actually use for these tables
  expect_snapshot(writeLines(.efa_corr_lines(make_wide_vars())))
  expect_snapshot(writeLines(.efa_corr_lines(make_wide_phi(), lower_only = TRUE)))
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

test_that("format.efa_sl_loadings flags a Heywood case", {
  local_reproducible_output()
  expect_snapshot(print(make_sl(heywood = TRUE)))
})

test_that("format.efa_sl_loadings prints cleanly without Heywood cases", {
  local_reproducible_output()
  expect_snapshot(print(make_sl(heywood = FALSE)))
})

test_that("format.efa_sl_loadings counts a communality-only Heywood item once", {
  # h2 > 1 (and u2 < 0) with no single loading > 1 must be flagged exactly once:
  # not missed (as a loading-only count would), nor multiplied across coupled cells.
  sl <- matrix(c(0.70, 0.20, 0.10, 0.53, 0.47,
                 0.80, 0.75, 0.10, 1.20, -0.20,
                 0.40, 0.10, 0.55, 0.46, 0.54),
               nrow = 3, byrow = TRUE,
               dimnames = list(c("i1", "i2", "i3"), c("g", "F1", "F2", "h2", "u2")))
  class(sl) <- c("efa_sl_loadings", "SLLOADINGS")
  out <- cli::ansi_strip(utils::capture.output(print(sl)))
  expect_true(any(grepl("Results contain a Heywood case!", out, fixed = TRUE)))
  expect_false(any(grepl("Heywood cases", out, fixed = TRUE)))
})

test_that("format.efa_sl_loadings is the source of truth and honours the colour state", {
  sl <- make_sl(heywood = TRUE)

  # print() is cat(format(x), sep = "\n") plus one blank line for console spacing
  expect_identical(utils::capture.output(print(sl)), c(format(sl), ""))

  old <- options(cli.num_colors = 256)
  on.exit(options(old), add = TRUE)
  expect_true(any(grepl("\033", format(sl), fixed = TRUE)))
  # `color = FALSE` governs the table cells, not the Heywood alert below them, which
  # cli styles from the active theme -- so it is not an ANSI switch for the whole block.
  # The table is the header row plus one row per item; the alert follows a blank line.
  table_lines <- seq_len(nrow(sl) + 1L)
  expect_false(any(grepl("\033", format(sl, color = FALSE)[table_lines], fixed = TRUE)))

  options(cli.num_colors = 1)
  plain <- format(sl)
  expect_false(any(grepl("\033", plain, fixed = TRUE)))
  # no trailing blank element
  expect_true(nzchar(plain[length(plain)]))
})

test_that("print.efa_loadings argument validators raise classed conditions", {
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

test_that("the shared argument checks keep their caller's frame and reject huge values", {
  L <- make_loadings()

  # The checks are shared helpers, so each abort has to report the validator that ran it
  # rather than the helper itself -- otherwise every message names the same internal
  # function whatever argument was wrong.
  cnd <- tryCatch(print(L, digits = -1), error = function(e) e)
  expect_s3_class(cnd, "efa_print_invalid_digits")
  expect_match(deparse(conditionCall(cnd))[1], ".validate_loadings_print_args",
               fixed = TRUE)

  # A whole number beyond the integer range must still reach the classed condition: an
  # as.integer() comparison would make it NA and abort with a bare "missing value where
  # TRUE/FALSE needed" instead.
  expect_error(print(L, digits = 1e10), class = "efa_print_invalid_digits")
  expect_error(print(L, max_name_length = 1e10),
               class = "efa_print_invalid_max_name_length")
  expect_error(print(L, max_factors_per_block = 3e9),
               class = "efa_print_invalid_max_factors_per_block")
  expect_error(format(efa_fit(test_models$baseline$cormat, n_factors = 2, N = 500,
                              estimator = "PAF"), digits = 1e10),
               class = "efa_print_invalid_digits")
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
