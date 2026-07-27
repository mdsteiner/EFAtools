# Two-step polychoric / tetrachoric correlation matrix from raw ordinal data. Recodes
# each variable to consecutive 0-based categories, then defers the numerics to the C++
# backend (.polychoric_cpp); this wrapper owns the input validation and the classed
# conditions, which can only be raised at the R level.
.polychoric <- function(x, nearest_pd = FALSE,
                        binary_only = FALSE, acov = c("none", "diag", "full"),
                        label_acov = TRUE, error_call = rlang::caller_env()) {

  acov <- match.arg(acov)

  # data.matrix() maps a factor to its LEVEL INDEX, so the level order silently becomes the
  # response order. That is correct for an ordered factor (and for numeric codes, which pass
  # through untouched), but an unordered factor built over character labels is ordered
  # alphabetically -- "always" < "never" < "often" < "rarely" -- which scrambles the responses
  # and attenuates the correlations badly (a true rho of 0.69 estimates as 0.33) without
  # producing anything a user could notice: the result is still a clean, unit-diagonal matrix.
  # Character columns take the same alphabetical route. Ordering is load-bearing here, so
  # refuse to guess it rather than accept the input on a silently arbitrary order.
  if (is.data.frame(x)) {
    bad <- vapply(x, function(col) is.character(col) ||
                    (is.factor(col) && !is.ordered(col)), logical(1))
    if (any(bad)) {
      nms_bad <- names(x)[bad]
      cli::cli_abort(
        c("Polychoric correlations need an explicit response order, but {cli::qty(nms_bad)} column{?s} {.val {nms_bad}} {?is/are} unordered.",
          "x" = "An unordered factor or character column would be ranked by its level order, which is alphabetical unless you set it.",
          "i" = "Supply numeric response codes, or convert with {.fun ordered} using the correct level order."),
        class = "efa_cor_unordered_factor", call = error_call)
    }
  } else if (is.character(x)) {
    cli::cli_abort(
      c("Polychoric correlations need numeric response codes, but {.arg x} is a character matrix.",
        "i" = "Supply numeric response codes, or convert the columns with {.fun ordered} using the correct level order."),
      class = "efa_cor_unordered_factor", call = error_call)
  }

  # Ordered-factor columns of a data frame become their level indices; numeric data passes
  # through. The recoding below maps each column onto its own 0-based categories.
  x <- data.matrix(x)
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("V", seq_len(ncol(x)))
  }
  nms <- colnames(x)
  p <- ncol(x)

  # A single variable has no correlation to estimate. Caught here so the request fails with a
  # classed condition; the C++ backend keeps the same check as unreachable defence.
  if (p < 2L) {
    cli::cli_abort(
      c("Polychoric correlations need at least two variables.",
        "x" = "{.arg x} has {p} column{?s}."),
      class = "efa_cor_too_few_vars", call = error_call)
  }

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

  out <- list(R = R)

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
      # Every pair is tabulated (rather than stopping at the first sparse one) so the warning
      # can name the variables the user has to act on. Each table is a single tabulate() and
      # this runs once per fit, not per bootstrap replicate (label_acov = FALSE skips it).
      # Only the pair indices are collected here; .pair_labels() below supplies the names, so
      # the warning and the covariance's dimnames cannot drift apart.
      is_sparse <- logical(ncol(pairs_ij))
      for (k in seq_len(ncol(pairs_ij))) {
        i <- pairs_ij[1L, k]; j <- pairs_ij[2L, k]
        tab <- matrix(tabulate(codes[, i] * n_cat[j] + codes[, j] + 1L,
                               nbins = n_cat[i] * n_cat[j]),
                      n_cat[i], n_cat[j], byrow = TRUE)
        expected <- outer(rowSums(tab), colSums(tab)) / n_obs
        is_sparse[k] <- any(tab == 0L & expected >= 5)
      }

      labels <- .pair_labels(nms)

      if (any(is_sparse)) {
        sparse <- labels[is_sparse]
        # Name the offending pairs, but cap the list: a heterogeneous ordinal set can have
        # dozens, and a warning that prints them all is unreadable.
        n_sparse <- length(sparse)
        shown <- sparse[seq_len(min(5L, n_sparse))]
        more <- n_sparse - length(shown)
        # When the list is truncated, drop cli's "and" before the last shown pair so the
        # trailing count closes the enumeration instead of adding a second conjunction.
        rest <- ""
        if (more > 0L) {
          shown <- cli::cli_vec(shown, style = list("vec-last" = ", "))
          rest <- paste0(", and ", more, " more")
        }
        cli::cli_warn(
          c("{n_sparse} variable pair{?s} {?has/have} an empty response-category combination despite a non-negligible expected count.",
            "x" = "Affected {cli::qty(n_sparse)}pair{?s}: {.val {shown}}{rest}.",
            "i" = "The polychoric asymptotic covariance (and any DWLS weights or robust standard errors derived from it) can be unreliable for such structurally sparse cells; interpret them with caution and consider collapsing rare response categories in these variables."),
          class = "efa_cor_sparse_cells")
      }

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
