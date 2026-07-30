#' Print and format a multigroup factor analysis
#'
#' `print()` turns an [efa_group()] result into a sectioned report: a header recapping
#' the groups, the common number of factors, the estimator, the rotation, and the
#' alignment; a group-pair table of the matched Tucker congruences between the aligned
#' loadings; a per-pair summary of the cross-group loading differences (with the salient
#' and, when a bootstrap was run, the confidence-interval flags); and, when
#' `invariance = TRUE`, a group-pair by factor grid of the approximate-invariance verdicts.
#' `format()` assembles the same report and returns it as a character vector; `print()` is
#' `cat(format(x), sep = "\n")`. The lines follow the active console theme, so they are
#' plain when colours are disabled (for example when captured into a file or stripped with
#' [cli::ansi_strip()]). `print()` does not draw a plot; use [plot.efa_group()].
#'
#' @param x An object of class `efa_group` (output from [efa_group()]).
#' @param digits Integer. The number of decimal places the reported values are rounded to.
#'   Default is 3.
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a character
#'   vector with the report lines.
#'
#' @family factor analysis
#'
#' @export
#' @method print efa_group
#'
#' @examples
#' g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
#' mg <- efa_group(GRiPS_raw, groups = g, n_factors = 1)
#' mg
#'
#' # format() returns the same lines as a character vector:
#' writeLines(format(mg))
#'
print.efa_group <- function(x, digits = 3, ...) {
  cat(format(x, digits = digits, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.efa_group
#' @export
#' @method format efa_group
format.efa_group <- function(x, digits = 3, ...) {

  settings <- x$settings
  groups <- settings$groups
  m <- length(groups)
  k <- settings$n_factors
  # One row per group pair, in the same order as $diffs / $flags / $invariance.
  pairs <- paste(x$diffs$group_1, "vs", x$diffs$group_2)
  fac_names <- dimnames(x$congruence$matched)[[3L]]
  n_boot <- x$congruence$n_boot
  boot <- !is.null(n_boot)

  cli::cli_format_method({

    # -- Header ----------------------------------------------------------------
    cli::cli_rule(left = "{.strong Multigroup exploratory factor analysis}")
    cli::cli_text("")
    rot_txt <- if (is.null(settings$rotation) || identical(settings$rotation, "none")) {
      "unrotated"
    } else {
      paste0(settings$rotation, " rotation")
    }
    # Emitted verbatim so arbitrary group names are never interpreted as cli/glue markup.
    .efa_group_verbatim(paste0(
      m, " groups (", paste(groups, collapse = ", "), ") \u00b7 ",
      k, if (k == 1L) " factor" else " factors", " \u00b7 ",
      settings$estimator, " extraction \u00b7 ", rot_txt))
    .efa_group_verbatim(if (identical(settings$alignment, "consensus")) {
      "Aligned to a symmetric consensus target."
    } else {
      paste0("Aligned to reference group ", settings$reference_group, ".")
    })
    .efa_group_verbatim(paste0(
      "N: ", paste0(names(settings$N), " = ",
                    format(settings$N, trim = TRUE), collapse = ", ")))

    # -- Factor congruence -----------------------------------------------------
    .print_efa_rule("Factor congruence")
    cli::cli_text("Tucker's congruence between the aligned group loadings (matched factors):")
    cli::cli_text("")

    phi_mat <- matrix(NA_real_, nrow = length(pairs), ncol = k,
                      dimnames = list(pairs, fac_names))
    r <- 0L
    for (g in seq_len(m - 1L)) {
      for (h in (g + 1L):m) {
        r <- r + 1L
        phi_mat[r, ] <- x$congruence$matched[g, h, ]
      }
    }
    .efa_emit_lines(.efa_format_matrix(
      values = phi_mat, row_labels = pairs, col_labels = fac_names,
      col_roles = rep("corr", k), digits = digits))

    if (any(x$congruence$degenerate)) {
      cli::cli_text("")
      cli::cli_alert_info(
        "{.val NA}: a near-zero factor left the congruence undefined for that pair (see {.code $congruence$degenerate}).",
        wrap = TRUE)
    }
    if (boot) {
      cli::cli_text("")
      cli::cli_alert_info(
        "{round(settings$ci * 100)}% percentile bootstrap CIs ({n_boot} replicate{?s}) are in {.code $congruence$matched_ci}.",
        wrap = TRUE)
    }

    # -- Loading differences ---------------------------------------------------
    .print_efa_rule("Loading differences")
    dstr <- .efa_group_delta_str(settings$delta)
    cli::cli_text("Absolute loading differences by group pair (salience threshold |\u0394| \u2265 {dstr}):")
    cli::cli_text("")

    pk <- nrow(x$loadings[[1L]]) * k
    d <- x$diffs
    diff_cells <- cbind(
      .efa_num(d$mean_abs_diff, digits = digits, pad = FALSE),
      .efa_num(d$median_abs_diff, digits = digits, pad = FALSE),
      .efa_num(d$min_abs_diff, digits = digits, pad = FALSE),
      .efa_num(d$max_abs_diff, digits = digits, pad = FALSE),
      .efa_num(d$rmse, digits = digits, pad = FALSE),
      paste0(d$n_flagged, "/", pk))
    .efa_emit_lines(.efa_group_table(
      pairs, c("mean", "median", "min", "max", "rmse", "flagged"), diff_cells))

    if (boot) {
      n_exc <- sum(x$flags$ci_excludes_0, na.rm = TRUE)
      if (n_exc > 0L) {
        cli::cli_text("")
        cli::cli_alert_info(
          "{n_exc} loading{?s} {?has/have} a bootstrap CI excluding 0 (see {.code $flags}).",
          wrap = TRUE)
      }
    }

    # -- Approximate invariance (only when requested) --------------------------
    if (!is.null(x$invariance)) {
      .print_efa_rule("Approximate invariance")
      if (boot) {
        cli::cli_text("Congruence bands (Lorenzo-Seva & ten Berge, 2006), read off the CI lower bound:")
      } else {
        cli::cli_text("Congruence bands (Lorenzo-Seva & ten Berge, 2006: \u2265 .95 equal, \u2265 .85 fair):")
      }
      cli::cli_text("")

      # $invariance is enumerated in the same pair-major, factor-minor order as
      # `pairs` x `fac_names` (all follow the shared g < h enumeration over the groups),
      # so its verdicts fill the grid positionally -- as the congruence table above does --
      # rather than by re-joining and matching group-name strings.
      verdicts <- matrix(x$invariance$verdict, nrow = length(pairs), ncol = k,
                         byrow = TRUE, dimnames = list(pairs, fac_names))
      verdicts[is.na(verdicts)] <- "NA"
      .efa_emit_lines(.efa_group_table(
        pairs, fac_names, verdicts, style = .efa_group_verdict_style))

      # Tucker's congruence is invariant to a proportional rescaling of a factor's
      # loadings, so uniformly stronger loadings in one group are still graded "equal"
      # -- a verdict grid that then sits directly above a difference table contradicting
      # it. Point at the differences when that gap is actually present: every factor
      # graded equal, yet some pair's *average* absolute difference already reaches the
      # threshold set for a single salient cell. `delta = 0` flags every cell by
      # construction and so says nothing about the size of the gap.
      max_mad <- max(d$mean_abs_diff)
      if (settings$delta > 0 && !anyNA(x$invariance$verdict) &&
          all(x$invariance$verdict == "equal") && max_mad >= settings$delta) {
        # At least two decimals whatever `digits` is: a difference that rounds to "0"
        # would contradict the sentence it is quoted in.
        mad_str <- .efa_num(max_mad, digits = max(digits, 2L), pad = FALSE)
        cli::cli_text("")
        cli::cli_alert_info(
          "Every factor is graded {.val equal} while the aligned loadings differ by up to {mad_str} on average: Tucker's congruence is invariant to a proportional rescaling of a factor's loadings, so read the verdicts alongside {.code $diffs}.",
          wrap = TRUE)
      }
    }

  })
}

# Render a small labelled table (left-aligned row labels, centred bold headers,
# right-aligned cells) into console lines, following the column maths of
# .efa_format_matrix() without its role machinery. `cells` is an already-formatted
# character matrix; `style(padded, raw)` optionally restyles a cell (the raw, unpadded
# value drives the styling; the padded cell keeps the column aligned). A table wider
# than the console is split into vertically stacked column blocks by the same helper
# .efa_format_matrix() uses; there are no auxiliary columns here, so every block repeats
# only the row labels and its own header row.
.efa_group_table <- function(row_labels, headers, cells, style = NULL) {
  cells <- as.matrix(cells)
  n_col <- ncol(cells)

  row_width <- max(cli::ansi_nchar(row_labels, type = "width"), 1L)
  col_widths <- vapply(seq_len(n_col), function(j) {
    max(cli::ansi_nchar(headers[j], type = "width"),
        max(cli::ansi_nchar(cells[, j], type = "width")))
  }, integer(1L))

  blocks <- .efa_matrix_blocks(seq_len(n_col), integer(0), row_width, col_widths)

  out <- character(0)
  for (bb in seq_along(blocks)) {
    cols <- blocks[[bb]]

    header_cells <- vapply(cols, function(j) {
      as.character(cli::ansi_align(cli::style_bold(headers[j]),
                                   width = col_widths[j], align = "center"))
    }, character(1L))
    header <- paste0(strrep(" ", row_width), "  ",
                     paste(header_cells, collapse = "  "))

    body <- vapply(seq_len(nrow(cells)), function(i) {
      row_cells <- vapply(cols, function(j) {
        padded <- as.character(cli::ansi_align(cells[i, j], width = col_widths[j],
                                               align = "right"))
        if (is.null(style)) padded else style(padded, cells[i, j])
      }, character(1L))
      label <- as.character(cli::ansi_align(row_labels[i], width = row_width,
                                            align = "left"))
      paste0(label, "  ", paste(row_cells, collapse = "  "))
    }, character(1L))

    # A blank line separates stacked blocks; the repeated header row names their columns.
    out <- c(out, if (bb > 1L) "", header, body)
  }

  sub("[ ]+$", "", out)
}

# Emit a header line verbatim, wrapped to the console width with continuation lines
# indented: cli_verbatim() keeps arbitrary group names out of cli/glue markup, but it also
# does no wrapping of its own. strwrap() breaks *below* its `width`, so it is given one
# column more than the console has in order to fill it exactly. It wraps on display width
# (so wide characters are counted correctly) but it also normalises runs of whitespace and
# breaks only at spaces, so a group name carrying doubled spaces prints with single ones
# and a name longer than the console still overruns it.
.efa_group_verbatim <- function(txt) {
  cli::cli_verbatim(strwrap(txt, width = .efa_console_width() + 1L, exdent = 2L))
}

# Colour an approximate-invariance verdict cell by its band (a no-op when colour is off,
# since .efa_style's cli primitives are then no-ops).
.efa_group_verdict_style <- function(padded, raw) {
  colour <- switch(trimws(raw), equal = "green", fair = "yellow",
                   incongruent = "red", NULL)
  if (is.null(colour)) padded else .efa_style(padded, colour)
}

# Format the salience threshold for display: locale-independent (sprintf always uses ".",
# unlike format() under a comma OutDec) with the leading zero dropped, matching the
# correlation-scale numbers in the tables (which flow through .efa_num / sprintf). Shared by
# the print header and the differences-plot subtitle so the two show the threshold alike.
.efa_group_delta_str <- function(delta) sub("^(-?)0\\.", "\\1.", sprintf("%g", delta))
