# Two-step polychoric / tetrachoric correlation matrix from raw ordinal data. Recodes
# each variable to consecutive 0-based categories, then defers the numerics to the C++
# backend (.polychoric_cpp); this wrapper owns the input validation and the classed
# conditions, which can only be raised at the R level.
.polychoric <- function(x, nearest_pd = FALSE,
                        binary_only = FALSE, acov = c("none", "diag", "full"),
                        label_acov = TRUE, error_call = rlang::caller_env()) {

  acov <- match.arg(acov)

  # Factor columns of a data frame become their integer codes; numeric data passes
  # through. The recoding below maps each column onto its own 0-based categories.
  x <- data.matrix(x)
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("V", seq_len(ncol(x)))
  }
  nms <- colnames(x)
  p <- ncol(x)

  # When an asymptotic covariance is requested, the point estimate, thresholds, and ACOV must
  # all come from the SAME cases: a sandwich covariance is only valid for the estimator that
  # produced the estimates, so a pairwise matrix paired with a listwise covariance is not the
  # covariance of any estimator (lavaan likewise computes categorical standard errors on a
  # single case set). Reduce to the listwise-complete rows up front, so the recode below also
  # collapses any category that survives only in incomplete rows. Without an ACOV the matrix
  # stays pairwise-complete for data efficiency (exploratory factor analysis needs no standard
  # errors from the matrix itself).
  if (acov != "none") {
    complete <- stats::complete.cases(x)
    if (sum(complete) < 2L) {
      cli::cli_abort(
        c("A polychoric asymptotic covariance needs at least two listwise-complete observations.",
          "x" = "{sum(complete)} row{?s} {?is/are} complete across all variables.",
          "i" = "Supply complete data, or impute the missing values, before requesting an {.arg acov}."),
        class = "efa_cor_no_complete_cases", call = error_call)
    }
    x <- x[complete, , drop = FALSE]
  }

  codes <- matrix(NA_integer_, nrow(x), p)
  n_cat <- integer(p)
  for (j in seq_len(p)) {
    col <- x[, j]
    lv <- sort(unique(col[!is.na(col)]))
    n_cat[j] <- length(lv)
    codes[, j] <- match(col, lv) - 1L
  }
  colnames(codes) <- nms

  # A constant (or empty) column has no thresholds and no information about any
  # correlation; abort rather than returning a degenerate matrix.
  if (any(n_cat < 2L)) {
    bad <- nms[n_cat < 2L]
    cli::cli_abort(
      c("Polychoric correlations need at least two response categories per variable.",
        "x" = "{cli::qty(bad)} Variable{?s} {.val {bad}} {?is/are} constant."),
      class = "efa_cor_constant_col", call = error_call)
  }

  # Tetrachoric correlations are the binary special case; reject any variable with
  # more than two categories so the request is not silently widened to polychoric.
  if (isTRUE(binary_only) && any(n_cat > 2L)) {
    bad <- nms[n_cat > 2L]
    cli::cli_abort(
      c("Tetrachoric correlations require binary variables (at most two categories).",
        "x" = "{cli::qty(bad)} Variable{?s} {.val {bad}} {?has/have} more than two categories."),
      class = "efa_cor_not_binary", call = error_call)
  }

  # Many categories usually means the variable is continuous, not ordinal; warn (rather than
  # abort, so a genuinely ordinal many-category variable still works) but proceed.
  if (any(n_cat >= 10L)) {
    many <- nms[n_cat >= 10L]
    cli::cli_warn(
      c("{cli::qty(many)} Variable{?s} {.val {many}} {?has/have} 10 or more response categories.",
        "i" = "Polychoric correlations assume a few ordered categories; check that {cli::qty(many)} {?this variable is/these variables are} not continuous."),
      class = "efa_cor_many_categories")
  }

  res <- .polychoric_cpp(codes, acov, isTRUE(nearest_pd))

  R <- res$R
  dimnames(R) <- list(nms, nms)
  thresholds <- res$thresholds
  names(thresholds) <- nms

  # A pair with no overlapping complete observations is uncomputable and comes back as
  # NA; abort with the same classed condition .prepare_cor_input() uses for an NA
  # correlation rather than returning a matrix EFA cannot use.
  if (anyNA(R)) {
    bad <- which(is.na(R) & upper.tri(R), arr.ind = TRUE)
    pairs <- paste0(nms[bad[, 1L]], "-", nms[bad[, 2L]])
    cli::cli_abort(
      c("The polychoric correlation could not be computed for {cli::qty(pairs)} variable pair{?s} {.val {pairs}}.",
        "i" = "{cli::qty(pairs)} {?This pair has/These pairs have} no overlapping complete observations."),
      class = "efa_cor_na", call = error_call)
  }

  # Surface the nearest-PD projection as the same classed warning .prepare_cor_input()
  # uses for smoothing, so downstream handling is consistent.
  if (isTRUE(res$pd_adjusted)) {
    cli::cli_warn(
      c("The polychoric correlation matrix was not positive definite; it has been projected to the nearest positive definite matrix.",
        "i" = "Inspect the results carefully."),
      class = "efa_cor_smoothed")
  }

  out <- list(R = R, thresholds = thresholds)

  # The asymptotic covariance of the off-diagonal correlations (Muthen, 1984; Joreskog,
  # 1994), at the requested level: a variance per element ("diag") or the full cross-pair
  # matrix ("full"), on the variance scale (Var(rho-hat)). It is the covariance of the
  # (un-projected) two-step estimates computed on the listwise-complete rows; `nearest_pd`
  # only smooths the returned R, it does not change the covariance. Pairs follow column order
  # of the upper triangle (i < j), so the labels match utils::combn(nms, 2).
  if (acov != "none") {
    av <- res$acov
    # The point estimate names the asymptotic covariance by variable pair and runs the
    # structural-sparsity diagnostic below; the per-replicate bootstrap recompute consumes the
    # covariance positionally and skips both (label_acov = FALSE), so neither is repeated per
    # replicate.
    if (label_acov) {
      # The two-step asymptotic covariance (and the DWLS weights / robust standard errors derived
      # from it) rely on large-sample theory that degrades for structurally empty contingency
      # cells, where the analytic variances can be far too small. Warn only when a cell is empty
      # despite a non-negligible expected count under independence (a structural rather than
      # incidental zero, by the usual >= 5 contingency-table rule); a rare corner combination in a
      # many-category table is not flagged. The covariance is still returned.
      n_obs <- nrow(codes)
      pairs_ij <- utils::combn(p, 2L)
      # Stop at the first structurally-sparse pair instead of building every pair's table.
      any_sparse <- FALSE
      for (k in seq_len(ncol(pairs_ij))) {
        i <- pairs_ij[1L, k]; j <- pairs_ij[2L, k]
        tab <- matrix(tabulate(codes[, i] * n_cat[j] + codes[, j] + 1L,
                               nbins = n_cat[i] * n_cat[j]),
                      n_cat[i], n_cat[j], byrow = TRUE)
        expected <- outer(rowSums(tab), colSums(tab)) / n_obs
        if (any(tab == 0L & expected >= 5)) { any_sparse <- TRUE; break }
      }
      if (any_sparse) {
        cli::cli_warn(
          c("Some response-category combinations are empty despite a non-negligible expected count.",
            "i" = "The polychoric asymptotic covariance (and any DWLS weights or robust standard errors derived from it) can be unreliable for such structurally sparse cells; interpret them with caution."),
          class = "efa_cor_sparse_cells")
      }

      labels <- .pair_labels(nms)
      if (acov == "diag") {
        names(av) <- labels
      } else {
        dimnames(av) <- list(labels, labels)
      }
    }
    out$acov <- av
  }

  out
}
