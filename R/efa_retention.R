# Unified result class for the factor-retention criteria. Every criterion
# returns an `efa_retention` object built by `.new_efa_retention()`, printed by
# one `print.efa_retention()`, and plotted by one `plot.efa_retention()` that
# dispatches to the two ggplot helpers below. The registry is the single source
# of truth mapping each criterion id to its display label.

# Registry of the factor-retention criteria (id -> metadata). Each `fun` runs
# its criterion from the N_FACTORS() control list `ctl`; `x` is the raw data
# when `needs_raw`, the prepared correlation matrix otherwise. Criteria that
# forward additional arguments to EFA() receive them via `ctl$dots`.
.retention_registry <- list(
  CD = list(
    label = "Comparison data", needs_raw = TRUE,
    fun = function(x, ctl) {
      CD(x, n_factors_max = ctl$n_factors_max, N_pop = ctl$N_pop,
         N_samples = ctl$N_samples, alpha = ctl$alpha, use = ctl$use,
         cor_method = ctl$cor_method, max_iter = ctl$max_iter_CD)
    }),
  EKC = list(
    label = "Empirical Kaiser Criterion", needs_raw = FALSE,
    fun = function(x, ctl) {
      EKC(x, N = ctl$N, use = ctl$use, cor_method = ctl$cor_method,
          type = ctl$ekc_type)
    }),
  HULL = list(
    label = "Hull method", needs_raw = FALSE,
    fun = function(x, ctl) {
      do.call(HULL, c(list(x, N = ctl$N, n_fac_theor = ctl$n_fac_theor,
                           method = ctl$method, gof = ctl$gof,
                           eigen_type = ctl$eigen_type_HULL, use = ctl$use,
                           cor_method = ctl$cor_method,
                           n_datasets = ctl$n_datasets, percent = ctl$percent,
                           decision_rule = ctl$decision_rule,
                           n_factors = ctl$n_factors),
                      ctl$dots))
    }),
  KGC = list(
    label = "Kaiser-Guttman criterion", needs_raw = FALSE,
    fun = function(x, ctl) {
      do.call(KGC, c(list(x, eigen_type = ctl$eigen_type_other, use = ctl$use,
                          cor_method = ctl$cor_method,
                          n_factors = ctl$n_factors, method = ctl$method),
                     ctl$dots))
    }),
  MAP = list(
    label = "Minimum average partial", needs_raw = FALSE,
    fun = function(x, ctl) {
      MAP(x, use = ctl$use, cor_method = ctl$cor_method)
    }),
  NEST = list(
    label = "Next Eigenvalue Sufficiency Test", needs_raw = FALSE,
    fun = function(x, ctl) {
      NEST(x, N = ctl$N, use = ctl$use, cor_method = ctl$cor_method,
           alpha = ctl$alpha_nest, n_datasets = ctl$n_datasets_nest,
           method = ctl$method)
    }),
  PARALLEL = list(
    label = "Parallel analysis", needs_raw = FALSE,
    fun = function(x, ctl) {
      do.call(PARALLEL, c(list(x, N = ctl$N, n_datasets = ctl$n_datasets,
                               percent = ctl$percent,
                               eigen_type = ctl$eigen_type_other, use = ctl$use,
                               cor_method = ctl$cor_method,
                               decision_rule = ctl$decision_rule,
                               n_factors = ctl$n_factors, method = ctl$method),
                          ctl$dots))
    }),
  SCREE = list(
    label = "Scree plot", needs_raw = FALSE, visual = TRUE,
    fun = function(x, ctl) {
      do.call(SCREE, c(list(x, eigen_type = ctl$eigen_type_other, use = ctl$use,
                            cor_method = ctl$cor_method,
                            n_factors = ctl$n_factors, method = ctl$method),
                       ctl$dots))
    }),
  SMT = list(
    label = "Sequential model tests", needs_raw = FALSE,
    fun = function(x, ctl) {
      SMT(x, N = ctl$N, use = ctl$use, cor_method = ctl$cor_method)
    })
)

# Construct an efa_retention object from a list of per-sub-variant records. The
# top-level `n_factors` named vector is derived from the records. `subtitle` is an
# optional one-line context string (e.g. the estimation method); `note` is an
# optional vector of cli info lines.
.new_efa_retention <- function(id, results, settings, status = "ok",
                               subtitle = NULL, note = NULL) {

  if (!id %in% names(.retention_registry)) {
    cli::cli_abort("Unknown factor-retention criterion id {.val {id}}.",
                   class = "efa_unknown_criterion")
  }

  label <- .retention_registry[[id]]$label

  n_factors <- vapply(results, function(r) as.numeric(r$n_factors), numeric(1))
  names(n_factors) <- vapply(results, function(r) r$name, character(1))

  structure(
    list(
      criterion = c(id = id, label = label),
      n_factors = n_factors,
      results = results,
      subtitle = subtitle,
      note = note,
      settings = settings,
      status = status
    ),
    class = "efa_retention"
  )
}

# Bullet lines (one per record) shared by format.efa_retention and
# format.N_FACTORS.
.retention_bullets <- function(results) {
  vapply(results, function(r) {
    value <- if (is.na(r$n_factors)) "not applicable" else as.character(r$n_factors)
    paste0(r$label, ": ", value)
  }, character(1))
}

#' Format method for efa_retention objects
#'
#' @param x an object of class efa_retention, returned by a factor-retention
#'   criterion (e.g. [EKC()] or [HULL()]).
#' @param ... not used.
#'
#' @returns A character vector with the formatted (plain, un-styled) output.
#'
#' @export
#' @method format efa_retention
#'
#' @examples
#' format(EKC(test_models$baseline$cormat, N = 500))
format.efa_retention <- function(x, ...) {
  cli::cli_format_method({
    cli::cli_rule(left = x$criterion[["label"]])
    if (!is.null(x$subtitle)) {
      cli::cli_text("{x$subtitle}")
    }
    # criteria with no numeric suggestion (e.g. the visual scree plot) skip the
    # bullets and rely on their subtitle/note
    if (any(!is.na(x$n_factors))) {
      cli::cli_text("")
      cli::cli_ul(.retention_bullets(x$results))
    }
    if (!is.null(x$note)) {
      cli::cli_text("")
      for (msg in x$note) {
        cli::cli_alert_info("{msg}")
      }
    }
  })
}

#' Print method for efa_retention objects
#'
#' @param x an object of class efa_retention, returned by a factor-retention
#'   criterion (e.g. [EKC()] or [HULL()]).
#' @param ... not used.
#'
#' @returns Invisibly returns `x`. Called for the side effect of printing the
#'   suggested number of factors.
#'
#' @export
#' @method print efa_retention
#'
#' @examples
#' EKC(test_models$baseline$cormat, N = 500)
print.efa_retention <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' Plot method for efa_retention objects
#'
#' Plots the result of a factor-retention criterion. Eigenvalue-based criteria
#' (e.g. [EKC()]) are shown as an eigenvalue plot, the Hull method ([HULL()]) as
#' a convex-hull plot. Criteria with more than one sub-variant are faceted.
#'
#' @param x an object of class efa_retention, returned by a factor-retention
#'   criterion (e.g. [EKC()] or [HULL()]).
#' @param ... not used.
#'
#' @returns A [ggplot2::ggplot] object, or invisibly `NULL` if the criterion has
#'   no plottable result.
#'
#' @export
#' @method plot efa_retention
#'
#' @examples
#' plot(EKC(test_models$baseline$cormat, N = 500))
plot.efa_retention <- function(x, ...) {

  plot_types <- unique(vapply(x$results, function(r) r$plot_type, character(1)))
  plot_types <- setdiff(plot_types, "none")

  if (length(plot_types) == 0) {
    cli::cli_inform("No plot is available for {x$criterion[['label']]}.")
    return(invisible(NULL))
  }

  if ("hull" %in% plot_types) {
    .gg_hull_plot(x)
  } else {
    .gg_eigen_plot(x)
  }

}

# Eigenvalue plot (factor index on x, eigenvalues on y) with an optional dashed
# reference series, an optional horizontal threshold, and an optional highlighted
# retained point. Covers the EKC/KGC/PARALLEL/SCREE/CD plots. Returns a ggplot.
#' @importFrom rlang .data
.gg_eigen_plot <- function(x) {

  # drop records with no plottable points (e.g. CD when it suggests 0 factors)
  records <- Filter(function(r) length(r$x) > 0, x$results)
  if (length(records) == 0) {
    cli::cli_inform("No plot is available for {x$criterion[['label']]}.")
    return(invisible(NULL))
  }

  variant_levels <- vapply(records, function(r) r$label, character(1))

  dat <- do.call(rbind, lapply(records, function(r) {
    data.frame(
      variant = r$label,
      factor = r$x,
      # a record may have no primary series (e.g. PARALLEL without real data
      # plots only its reference series)
      value = if (is.null(r$y)) NA_real_ else r$y,
      reference = if (is.null(r$reference)) NA_real_ else r$reference,
      stringsAsFactors = FALSE
    )
  }))
  dat$variant <- factor(dat$variant, levels = variant_levels)

  highlights <- do.call(rbind, lapply(records, function(r) {
    if (is.null(r$highlight) || is.na(r$highlight) || r$highlight < 1) return(NULL)
    data.frame(variant = r$label, factor = r$highlight,
               value = r$y[r$highlight], stringsAsFactors = FALSE)
  }))

  thresholds <- do.call(rbind, lapply(records, function(r) {
    if (is.null(r$threshold)) return(NULL)
    data.frame(variant = r$label, yintercept = r$threshold,
               stringsAsFactors = FALSE)
  }))

  p <- ggplot2::ggplot(dat, ggplot2::aes(.data$factor, .data$value)) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::geom_point(na.rm = TRUE)

  if (any(!is.na(dat$reference))) {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$reference),
                                linetype = 2, colour = "darkgray", na.rm = TRUE)
  }

  if (!is.null(thresholds)) {
    thresholds$variant <- factor(thresholds$variant, levels = variant_levels)
    p <- p + ggplot2::geom_hline(data = thresholds,
                                 ggplot2::aes(yintercept = .data$yintercept),
                                 linetype = 2, colour = "darkgray")
  }

  # records can carry several named reference series (e.g. PARALLEL's simulated
  # means and percentile eigenvalues), each drawn as a dashed coloured line with
  # a shared legend; criteria without `references` are unaffected
  ref_records <- Filter(function(r) !is.null(r$references), records)
  if (length(ref_records) > 0) {
    ref_dat <- do.call(rbind, lapply(ref_records, function(r) {
      data.frame(variant = r$label,
                 factor = rep(r$x, length(r$references)),
                 series = rep(names(r$references), each = length(r$x)),
                 value = unlist(r$references, use.names = FALSE),
                 stringsAsFactors = FALSE)
    }))
    ref_dat$variant <- factor(ref_dat$variant, levels = variant_levels)
    ref_dat$series <- factor(ref_dat$series, levels = unique(ref_dat$series))
    # label the primary series in the legend alongside the reference series
    # (drawn over the identical primary line, so only the legend key is added)
    ref_labels <- vapply(ref_records, function(r) r$label, character(1))
    real_dat <- dat[dat$variant %in% ref_labels & !is.na(dat$value), ,
                    drop = FALSE]
    if (nrow(real_dat) > 0) {
      p <- p +
        ggplot2::geom_line(data = real_dat,
                           ggplot2::aes(linetype = "Real Eigenvalues")) +
        ggplot2::scale_linetype_manual(values = c(`Real Eigenvalues` = 1))
    }
    p <- p +
      ggplot2::geom_line(data = ref_dat,
                         ggplot2::aes(.data$factor, .data$value,
                                      colour = .data$series),
                         linetype = 2, na.rm = TRUE) +
      ggplot2::scale_colour_viridis_d(end = 0.8) +
      ggplot2::labs(colour = NULL, linetype = NULL)
  }

  if (!is.null(highlights)) {
    highlights$variant <- factor(highlights$variant, levels = variant_levels)
    p <- p +
      ggplot2::geom_point(data = highlights, shape = 1, size = 4, colour = "red") +
      ggplot2::geom_text(data = highlights, ggplot2::aes(label = .data$factor),
                         colour = "red", vjust = -1)
  }

  if (length(variant_levels) > 1) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$variant), scales = "free_y")
  }

  # y-axis label defaults to "Eigenvalues" but a criterion can override it (e.g.
  # CD plots mean RMSE of the eigenvalues)
  y_label <- records[[1]]$y_label
  if (is.null(y_label)) y_label <- "Eigenvalues"

  p +
    ggplot2::scale_x_continuous(breaks = seq_len(max(dat$factor))) +
    ggplot2::labs(x = "Factor", y = y_label, title = x$criterion[["label"]]) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))

}

# Convex-hull plot (degrees of freedom on x, goodness-of-fit on y). Points on the
# hull are emphasised and connected; the retained solution is highlighted in red;
# each point is labelled with its number of factors. Records whose hull is empty
# (the degenerate fewer-than-three-solutions case) are dropped. Returns a ggplot.
#' @importFrom rlang .data
.gg_hull_plot <- function(x) {

  results <- Filter(function(r) any(r$on_hull), x$results)

  if (length(results) == 0) {
    cli::cli_inform("No hull plot is available for {x$criterion[['label']]}.")
    return(invisible(NULL))
  }

  variant_levels <- vapply(results, function(r) r$label, character(1))

  dat <- do.call(rbind, lapply(results, function(r) {
    data.frame(variant = r$label, df = r$x, fit = r$y,
               nfac = r$point_labels, on_hull = r$on_hull,
               stringsAsFactors = FALSE)
  }))
  dat$variant <- factor(dat$variant, levels = variant_levels)

  retained <- do.call(rbind, lapply(results, function(r) {
    idx <- which(r$point_labels == r$highlight)
    if (length(idx) == 0) return(NULL)
    data.frame(variant = r$label, df = r$x[idx], fit = r$y[idx],
               stringsAsFactors = FALSE)
  }))

  p <- ggplot2::ggplot(dat, ggplot2::aes(.data$df, .data$fit)) +
    ggplot2::geom_line(data = dat[dat$on_hull, , drop = FALSE]) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$on_hull)) +
    ggplot2::geom_text(ggplot2::aes(label = .data$nfac), vjust = -1, size = 3) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "darkgray"),
                                 guide = "none")

  if (!is.null(retained)) {
    retained$variant <- factor(retained$variant, levels = variant_levels)
    p <- p + ggplot2::geom_point(data = retained, shape = 1, size = 4,
                                 colour = "red")
  }

  if (length(variant_levels) > 1) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$variant), scales = "free_y")
  }

  p +
    ggplot2::labs(x = "Degrees of freedom", y = "Goodness of fit",
                  title = x$criterion[["label"]]) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))

}
