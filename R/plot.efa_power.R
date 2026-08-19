#' Plot the RMSEA power curve
#'
#' @description
#' Draws the analytic RMSEA power (MacCallum, Browne, & Sugawara, 1996) of an
#' [efa_power()] result as a function of the total sample size, mirroring
#' `semTools::plotRMSEApower()` but returning a [ggplot2::ggplot] object rather than
#' drawing to the active device. The test, its null and alternative RMSEA, the
#' significance level, and the number of groups are taken from the object; only the
#' sample-size axis is swept, with an optional sweep of the degrees of freedom or the
#' alternative RMSEA to overlay several curves.
#'
#' @details
#' When the plotted curve is the object's own -- a single curve with neither `df` nor
#' `eps1` overridden -- it is annotated with the object's result: a dashed vertical line
#' at its sample size `x$N`, a dashed horizontal line at the reference power (the target
#' power when a sample size was solved for, otherwise the power achieved at `x$N`), and a
#' point at `x$N` and the achieved power. Overriding `df` or `eps1`, sweeping either as a
#' vector, or supplying an `n` that does not span `x$N` moves that point off the drawn
#' curve, so the marks are then omitted.
#'
#' @param x An object of class `efa_power` (output from [efa_power()]).
#' @param n numeric. The total sample sizes to evaluate. If `NULL` (the default) a
#'   sequence bracketing the object's sample size is chosen automatically.
#' @param df numeric. The model degrees of freedom (must be positive). Defaults to the
#'   object's `df`; a vector of length greater than one draws one curve per value.
#' @param eps1 numeric. The alternative-hypothesis RMSEA (must differ from the null
#'   `eps0`). Defaults to the object's `eps1`; a vector of length greater than one draws
#'   one curve per value. At most one of `df` and `eps1` may be a vector.
#' @param ... Not used; for consistency with the generic.
#'
#' @returns A [ggplot2::ggplot] object.
#'
#' @references
#' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis and
#' determination of sample size for covariance structure modeling. *Psychological
#' Methods, 1*(2), 130-149. \doi{10.1037/1082-989X.1.2.130}
#'
#' @family power analysis
#'
#' @importFrom rlang .data
#' @export
#' @method plot efa_power
#'
#' @examples
#' pw <- efa_power(df = 100, N = 200)
#'
#' # Power curve for the test of close fit, marking the object's own N
#' plot(pw)
#'
#' # Overlay several models by sweeping the degrees of freedom
#' plot(pw, df = c(50, 100, 200))
#'
#' # Sweep the alternative RMSEA instead
#' plot(pw, eps1 = c(0.06, 0.08, 0.10))
plot.efa_power <- function(x, n = NULL, df = NULL, eps1 = NULL, ...) {
  # Simulation-mode results are Monte-Carlo rates, not an analytic power curve, so
  # there is no sample-size axis to sweep.
  if (identical(x$settings$mode, "simulation")) {
    cli::cli_abort(
      c("There is no power curve to plot for a simulation-mode {.cls efa_power} object.",
        "i" = "Inspect the hit-rate, recovery, and convergence rates in the printed summary instead."),
      class = "efa_power_no_plot")
  }
  s <- x$settings
  if (is.null(df)) df <- s$df
  if (is.null(eps1)) eps1 <- s$eps1
  checkmate::assert_numeric(df, finite = TRUE, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(eps1, lower = 0, finite = TRUE, any.missing = FALSE,
                            min.len = 1)
  # An override must uphold the same invariants efa_power() enforces on the object: a
  # positive df, and an alternative RMSEA distinct from the null (else there is nothing
  # to detect). Reuse the constructor's condition classes so callers can catch either.
  if (any(df <= 0)) {
    cli::cli_abort("{.arg df} must be positive.", class = "efa_power_bad_df")
  }
  if (any(eps1 == s$eps0)) {
    cli::cli_abort(
      c("{.arg eps1} must differ from the null RMSEA {.arg eps0} ({s$eps0}).",
        "i" = "Equal null and alternative RMSEA leave no difference to detect."),
      class = "efa_power_unreachable")
  }

  # "Sweep df OR eps1": at most one dimension may vary, so a single colour aesthetic
  # carries the overlay. A 2-D grid of both would need faceting and is out of scope.
  df_swept <- length(df) > 1L
  eps1_swept <- length(eps1) > 1L
  if (df_swept && eps1_swept) {
    cli::cli_abort(
      c("Only one of {.arg df} and {.arg eps1} can be swept at a time.",
        "i" = "Give a vector to one of them and leave the other a single value."),
      class = "efa_power_plot_grid")
  }
  swept <- df_swept || eps1_swept

  # The swept argument supplies the colour levels; the other is held fixed. curve_df and
  # curve_eps1 pair element-wise with the levels for the per-curve power evaluation.
  if (eps1_swept) {
    levels_val <- eps1
    level_lab <- "RMSEA (H1)"
    curve_df <- rep(df[1L], length(levels_val))
    curve_eps1 <- eps1
  } else {
    levels_val <- df
    level_lab <- "df"
    curve_df <- df
    curve_eps1 <- rep(eps1[1L], length(levels_val))
  }

  # Sample-size axis: honour a supplied `n`, otherwise span from 1 up to where the
  # slowest curve saturates (its N for 0.99 power), so every curve reaches its plateau
  # and the object's N sits inside. `.efa_power_solve_N()` aborts (efa_power_unreached)
  # when a target is unreachable (e.g. wrong-side epsilons); fall back to twice the
  # object's N there, but let any other error surface.
  if (is.null(n)) {
    hi <- vapply(seq_along(levels_val), function(i) {
      tryCatch(
        .efa_power_solve_N(curve_df[i], s$eps0, curve_eps1[i], s$alpha, s$group,
                           s$type, target = 0.99),
        efa_power_unreached = function(e) 2 * x$N)
    }, numeric(1))
    n_hi <- max(hi, x$N)
    n <- unique(round(seq(1, n_hi, length.out = 100L)))
  } else {
    checkmate::assert_numeric(n, lower = 1, finite = TRUE, any.missing = FALSE,
                              min.len = 1)
  }

  # Power on the (N x curve) grid, one long-format row per point.
  rows <- lapply(seq_along(levels_val), function(i) {
    pw <- vapply(n, function(nn) {
      .efa_power_rmsea(nn, curve_df[i], s$eps0, curve_eps1[i], s$alpha, s$group,
                       s$type)$power
    }, numeric(1))
    data.frame(N = n, power = pw, level = levels_val[i])
  })
  dat <- do.call(rbind, rows)
  dat$level <- factor(dat$level, levels = unique(levels_val))

  fit_lbl <- if (s$type == "close") "close fit" else "not-close fit"
  n_lab <- if (s$group > 1) "Total sample size (N)" else "Sample size (N)"
  # The null bounds RMSEA from above for close fit and from below for not-close fit,
  # matching the comparator in format.efa_power().
  cmp <- if (s$type == "close") "\u2264" else "\u2265"
  sub <- paste0("alpha = ", .efa_num(s$alpha, pad = FALSE),
                " \u00b7 df = ", if (df_swept) "swept" else df[1L],
                " \u00b7 H0 RMSEA ", cmp, " ", .efa_num(s$eps0, pad = FALSE))
  if (!eps1_swept) {
    sub <- paste0(sub, " \u00b7 H1 RMSEA = ", .efa_num(eps1[1L], pad = FALSE))
  }

  # Annotate the object's own point only when the drawn curve is the object's -- a single
  # unmodified curve whose sample-size axis actually spans x$N. Any override or an `n`
  # that excludes x$N would place (x$N, x$power) off the curve.
  mark_object <- !swept && df[1L] == s$df && eps1[1L] == s$eps1 &&
    x$N >= min(n) && x$N <= max(n)

  if (swept) {
    p <- ggplot2::ggplot(dat, ggplot2::aes(.data$N, .data$power, colour = .data$level)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::scale_colour_viridis_d(end = 0.85) +
      ggplot2::labs(colour = level_lab)
  } else {
    p <- ggplot2::ggplot(dat, ggplot2::aes(.data$N, .data$power))
    if (mark_object) {
      # Reference power: the target when N was solved for, else the power at N.
      ref_power <- if (x$solve_for == "N") s$power else x$power
      p <- p +
        ggplot2::geom_hline(yintercept = ref_power, linetype = 2, colour = "grey50") +
        ggplot2::geom_vline(xintercept = x$N, linetype = 2, colour = "grey50")
    }
    p <- p + ggplot2::geom_line(linewidth = 0.8, colour = "#3B528B")
    if (mark_object) {
      p <- p + ggplot2::geom_point(data = data.frame(N = x$N, power = x$power),
                                   size = 2.5, colour = "#3B528B")
    }
  }

  p +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = n_lab, y = "Power",
                  title = paste0("Power for the RMSEA test of ", fit_lbl),
                  subtitle = sub) +
    .gg_theme()
}
