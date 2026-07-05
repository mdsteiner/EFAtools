#' Plot a multigroup factor analysis
#'
#' Two views of an [efa_group()] result, selected by `type`:
#'
#' - `"congruence"` (the default) plots the matched Tucker congruence of each factor
#'   between every group pair, with a percentile bootstrap confidence interval when one
#'   was computed (`b_boot > 0`). The Lorenzo-Seva and ten Berge (2006) reference bands
#'   (`.95` "equal", `.85` "fair") are drawn so a factor's cross-group similarity can be
#'   read against them at a glance.
#' - `"differences"` draws a heatmap of the signed cross-group loading differences (item
#'   by factor, one panel per group pair). Cells whose absolute difference reaches the
#'   salience threshold `delta` are outlined.
#'
#' @param x An object of class `efa_group` (output from [efa_group()]).
#' @param type character. Which plot to draw: `"congruence"` (per-factor congruence with
#'   confidence intervals) or `"differences"` (a per-item loading-difference heatmap).
#' @param ... Not used; for consistency with the generic.
#'
#' @returns A [ggplot2::ggplot] object.
#'
#' @references
#' Lorenzo-Seva, U., and ten Berge, J. M. F. (2006). Tucker's congruence coefficient as a
#' meaningful index of factor similarity. *Methodology*, 2, 57-64.
#' doi: 10.1027/1614-2241.2.2.57
#'
#' @family factor analysis
#'
#' @importFrom rlang .data
#' @export
#' @method plot efa_group
#'
#' @examples
#' g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
#' mg <- efa_group(GRiPS_raw, groups = g, n_factors = 1)
#'
#' # Per-factor congruence against the Lorenzo-Seva & ten Berge bands
#' plot(mg)
#'
#' # Per-item cross-group loading-difference heatmap
#' plot(mg, type = "differences")
plot.efa_group <- function(x, type = c("congruence", "differences"), ...) {
  type <- match.arg(type)
  if (type == "congruence") .gg_group_congruence(x) else .gg_group_differences(x)
}

# Per-factor matched Tucker congruence between every group pair, as points (or point
# ranges when a bootstrap supplied confidence intervals), against the Lorenzo-Seva & ten
# Berge (2006) reference bands. Returns a ggplot.
#' @importFrom rlang .data
.gg_group_congruence <- function(x) {
  M <- x$congruence$matched
  groups <- dimnames(M)[[1L]]
  fac_names <- dimnames(M)[[3L]]
  lower <- x$congruence$matched_ci$lower
  upper <- x$congruence$matched_ci$upper
  m <- length(groups)

  rows <- list()
  r <- 0L
  for (g in seq_len(m - 1L)) {
    for (h in (g + 1L):m) {
      for (f in seq_along(fac_names)) {
        r <- r + 1L
        rows[[r]] <- data.frame(
          pair = paste(groups[g], "vs", groups[h]),
          factor = fac_names[f],
          phi = M[g, h, f],
          lower = if (is.null(lower)) NA_real_ else lower[g, h, f],
          upper = if (is.null(upper)) NA_real_ else upper[g, h, f],
          stringsAsFactors = FALSE)
      }
    }
  }
  dat <- do.call(rbind, rows)
  dat$factor <- factor(dat$factor, levels = fac_names)
  dat$pair <- factor(dat$pair, levels = unique(dat$pair))
  has_ci <- any(!is.na(dat$lower))

  y_min <- suppressWarnings(min(c(dat$phi, dat$lower), na.rm = TRUE))
  if (!is.finite(y_min)) y_min <- 0.8

  # The .85 ("fair") and .95 ("equal") band edges mirror the cut-offs in `.invariance_band()`;
  # keep them in step if those are ever retuned.
  p <- ggplot2::ggplot(dat, ggplot2::aes(.data$factor, .data$phi, colour = .data$pair)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.85, ymax = 0.95,
                      fill = "grey80", alpha = 0.35) +
    ggplot2::geom_hline(yintercept = c(0.85, 0.95), linetype = 2, colour = "grey50")

  p <- if (has_ci) {
    p + ggplot2::geom_pointrange(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      position = ggplot2::position_dodge(width = 0.5), size = 0.5, na.rm = TRUE)
  } else {
    p + ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5),
                            size = 2.5, na.rm = TRUE)
  }

  p +
    ggplot2::coord_cartesian(ylim = c(min(0.8, y_min), 1)) +
    ggplot2::scale_colour_viridis_d(end = 0.8) +
    ggplot2::labs(x = "Factor", y = "Tucker's congruence",
                  title = "Factor congruence across groups", colour = "Group pair") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
}

# Heatmap of the signed cross-group loading differences (item by factor), one panel per
# group pair; cells reaching the salience threshold `delta` are outlined. Returns a ggplot.
#' @importFrom rlang .data
.gg_group_differences <- function(x) {
  d <- x$flags
  pair_labels <- paste(d$group_1, "vs", d$group_2)
  d$pair <- factor(pair_labels, levels = unique(pair_labels))
  d$indicator <- factor(d$indicator, levels = rev(unique(d$indicator)))
  d$factor <- factor(d$factor, levels = unique(d$factor))
  # which() (not the logical vector) so an NA `flagged`, should one ever arise, drops out
  # rather than injecting an all-NA row that geom_tile would try to outline.
  flagged <- d[which(d$flagged), , drop = FALSE]

  lim <- max(abs(d$diff), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1

  p <- ggplot2::ggplot(d, ggplot2::aes(.data$factor, .data$indicator)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$diff), colour = "grey92")
  if (nrow(flagged) > 0L) {
    p <- p + ggplot2::geom_tile(data = flagged, fill = NA, colour = "black",
                                linewidth = 0.7)
  }

  p +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-lim, lim)) +
    ggplot2::facet_wrap(ggplot2::vars(.data$pair)) +
    ggplot2::labs(x = "Factor", y = NULL, fill = "Loading\ndifference",
                  title = "Cross-group loading differences",
                  subtitle = paste0("Outlined cells: |difference| >= ",
                                    .efa_group_delta_str(x$settings$delta))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
}
