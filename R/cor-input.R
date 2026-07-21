# Detect or build the analysis correlation matrix, then validate it once: assert the
# input type, smooth a non-positive-definite matrix, and resolve N. Shared by EFA and
# the suitability/retention functions so these checks live in one place.

# Unit diagonal on one shared tolerance: the correlation-vs-covariance classification
# hinges on .is_cormat() and .is_covmat() applying the SAME test, so any matrix falls on
# exactly one side of it.
.has_unit_diag <- function(x, na.rm = FALSE) {
  all(abs(diag(x) - 1) <= .Machine$double.eps * 100, na.rm = na.rm)
}

# Heuristically decide whether x is already a correlation matrix rather than raw data:
# square, symmetric, all entries in [-1, 1], and a unit diagonal. This cannot distinguish a
# correlation matrix from a (rare) raw-data matrix that happens to meet all four conditions,
# and it treats a covariance matrix as raw data (its diagonal is not all ones); callers that
# accept either input assume raw data unless these conditions hold.
.is_cormat <- function(x){

  if(nrow(x) == ncol(x) &&
     all(x >= (-1 + .Machine$double.eps * 100), na.rm = TRUE) &&
     all(x <= (1 + .Machine$double.eps * 100), na.rm = TRUE)){

    if (.has_unit_diag(x, na.rm = TRUE) &&
        isSymmetric(unclass(unname(x)))) {

      if (any(is.na(x))) {

        cli::cli_abort(
          c("{.arg x} looks like a correlation matrix but contains missing values.",
            "i" = "Please check the entered data."),
          class = "efa_cormat_has_na"
        )

      }

      TRUE

    } else {

      FALSE

    }


  } else {

    FALSE

  }

}

# Abort unless `x` is a matrix or data frame. `raw_only = TRUE` tailors the
# message for functions that accept only raw data (e.g. CD), which reject a
# correlation matrix downstream.
.assert_cor_input <- function(x, raw_only = FALSE, error_call = rlang::caller_env()) {
  if (!inherits(x, c("matrix", "data.frame"))) {
    lead <- if (raw_only) {
      "{.arg x} must be a data frame/matrix of raw data."
    } else {
      "{.arg x} must be a correlation matrix or a data frame/matrix of raw data."
    }
    cli::cli_abort(
      c(lead, "x" = "You supplied {.obj_type_friendly {x}}."),
      class = "efa_input_not_matrix",
      call = error_call
    )
  }
  if (.is_covmat(x)) {
    # cov2cor() is the fix only for a caller that accepts a correlation matrix; one that
    # needs raw data (raw_only) would just hit the correlation-matrix rejection next.
    cov_lead <- if (raw_only) {
      "{.arg x} looks like a covariance matrix, not a data frame/matrix of raw data."
    } else {
      "{.arg x} looks like a covariance matrix, not a correlation matrix or raw data."
    }
    cov_remedy <- if (raw_only) {
      "Supply the raw observations the covariance matrix was computed from."
    } else {
      "Standardise it first: {.code cov2cor(x)}."
    }
    cli::cli_abort(
      c(cov_lead,
        "x" = "It is square and symmetric, with variances rather than ones on the diagonal.",
        "i" = cov_remedy),
      class = "efa_input_is_covmat", call = error_call
    )
  }
  invisible(x)
}

# Recognise a covariance matrix, which is neither a correlation matrix (its diagonal is not
# all ones, so .is_cormat() returns FALSE) nor usable as raw data: passed on it would be fed
# to stats::cor(), i.e. the correlation of the covariance matrix's own p columns over p
# "cases", which is singular by construction and misdiagnoses the user's input. A square,
# symmetric, all-finite numeric matrix with positive non-unit diagonal entries is never raw
# data in practice -- symmetry alone is already decisive. Positive semi-definiteness is
# deliberately NOT required: the covariance matrices this guard exists for -- built with
# pairwise deletion, or transcribed from published tables at few decimals -- are routinely
# indefinite, and an indefinite correlation matrix is likewise accepted (and smoothed
# downstream), so demanding PSD of the covariance side only would send exactly the
# problematic inputs down the raw-data misdiagnosis. The positive-diagonal test keeps
# zero-diagonal distance/dissimilarity matrices out.
.is_covmat <- function(x) {

  # Cheapest test first, and before any coercion: raw data is almost never square, and
  # as.matrix() on a data frame copies the whole thing -- which would otherwise be paid on
  # every call for every raw-data input.
  if (ncol(x) < 2L || nrow(x) != ncol(x)) return(FALSE)

  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) return(FALSE)
    x <- as.matrix(x)
  }
  if (!is.numeric(x) || !all(is.finite(x))) return(FALSE)

  if (any(diag(x) <= 0)) return(FALSE)

  # A unit diagonal (on the shared tolerance, so the two classifiers are exactly
  # complementary) means this is a correlation matrix, handled by .is_cormat().
  if (.has_unit_diag(x)) return(FALSE)

  isSymmetric(unclass(unname(x)))
}

# Single source of truth for the correlation methods that are estimated from the
# raw data by .polychoric() rather than stats::cor(). Used both to route the
# computation and to reject criteria whose reference data are continuous.
.is_poly_cor <- function(cor_method) cor_method %in% c("poly", "tetra")

# Labels for the off-diagonal variable pairs in utils::combn(p, 2) column order
# ("name_i-name_j"), shared by every routine that names a pair-indexed vector or matrix
# (the polychoric and ADF asymptotic covariances and DWLS weights).
.pair_labels <- function(nms) {
  idx <- utils::combn(length(nms), 2L)
  paste(nms[idx[1L, ]], nms[idx[2L, ]], sep = "-")
}

# Build the symmetric per-element DWLS weight matrix W (W_ij = 1 / Var(rho_hat_ij), zero
# diagonal) from the diagonal polychoric asymptotic covariance returned by .polychoric():
# a length p(p - 1)/2 vector of off-diagonal variances ordered by the upper triangle
# (i < j) in utils::combn() column order. The combn order is row-major, which differs
# from R's column-major upper.tri(), so the pairs are placed via explicit (i, j) indexing
# rather than `W[upper.tri(W)] <- ...`. Shared by the DWLS point estimate and the per-
# replicate bootstrap recompute.
.poly_weight_matrix <- function(acov_diag, p) {
  # A non-positive or non-finite asymptotic variance means a (near-)degenerate variable pair
  # whose inverse-variance weight is undefined; refuse rather than emit Inf/NaN weights that
  # would silently corrupt the fit. (The bootstrap recompute catches this and drops the
  # replicate.)
  if (any(!is.finite(acov_diag)) || any(acov_diag <= 0)) {
    cli::cli_abort(
      c("A polychoric asymptotic variance was not positive, so the inverse-variance weights are undefined.",
        "i" = "A variable pair is (near-)degenerate (e.g. an empty or near-empty response category)."),
      class = "efa_dwls_degenerate_weight"
    )
  }
  W <- matrix(0, p, p)
  idx <- utils::combn(p, 2L)
  W[t(idx)] <- 1 / acov_diag
  W + t(W)
}

# One shared degeneracy check on the diagonal of the asymptotic covariance, applied where it
# is received so that every consumer -- the DWLS weights and the sandwich standard errors --
# takes the same branch on the same data. A (near-)degenerate variable pair, e.g. one whose
# contingency table is close to comonotone, is barely identified: its asymptotic variance
# explodes (values of 1e+26 are reachable), which leaves the DWLS weight ~0 so the pair drops
# out of the fit, but leaves the sandwich standard errors of that pair meaningless rather than
# merely large. Neither outcome should be silent, so warn once and name the pairs. Analysis
# proceeds: the variance is the honest answer for a pair the data cannot pin down. A variance
# above 1 already implies a +/- 1 SE interval as wide as the whole [-1, 1] range of a
# correlation, so it is the natural flag; non-finite and non-positive values are degenerate a
# fortiori.
.warn_acov_degenerate <- function(acov_diag, labels = names(acov_diag)) {
  bad <- !is.finite(acov_diag) | acov_diag <= 0 | acov_diag > 1
  if (!any(bad)) return(invisible(NULL))

  pairs <- if (is.null(labels)) paste0("pair ", which(bad)) else labels[bad]
  shown <- utils::head(pairs, 5L)
  more <- if (length(pairs) > length(shown)) {
    cli::format_inline(" and {length(pairs) - length(shown)} more")
  } else {
    ""
  }

  cli::cli_warn(
    c("The polychoric asymptotic covariance is degenerate for {cli::qty(pairs)} variable pair{?s} {.val {shown}}{more}.",
      "x" = "{cli::qty(pairs)} {?This pair is/These pairs are} barely identified (e.g. a near-empty or near-comonotone response table), so {?its/their} asymptotic variance is not usable.",
      "i" = "Any DWLS weight or robust/sandwich standard error involving {cli::qty(pairs)} {?it/them} is unreliable; a DWLS fit down-weights {cli::qty(pairs)} {?it/them} out of the solution.",
      "i" = "Consider collapsing rare response categories in the variables involved."),
    class = "efa_acov_degenerate"
  )
  invisible(NULL)
}

# Asymptotic-distribution-free (ADF; Browne, 1984) covariance of the off-diagonal sample
# correlations from raw continuous data, on the variance scale Var(rho-hat) and in
# utils::combn(p, 2) column order -- the continuous analogue of the polychoric asymptotic
# covariance, and the meat of the robust/sandwich standard errors. For case n and pair (i < j)
# the per-case influence is IF[n, (ij)] = (w_ni w_nj - 1/2 r_ij (w_ni^2 + w_nj^2)) / N, with w
# the columns standardised to unit variance; Gamma = crossprod(IF). The -1/2 r (w_ni^2 + w_nj^2)
# terms are the delta-method correction for estimating the marginal SDs (the standardisation
# nuisance, the continuous counterpart of the polychoric threshold influence). N * Gamma equals
# lavaan's correlation-structure NACOV (the basis of the MLM/MLR robust statistics). Kept in R:
# the cost is the BLAS crossprod (about 0.4 s at p = 40, N = 800); the influence build is cheap.
.adf_gamma <- function(x) {
  x <- as.matrix(x)
  N <- nrow(x)
  p <- ncol(x)

  # Centre each column and scale by its population SD, so colMeans(W^2) == 1 and crossprod(W)/N
  # is exactly the Pearson correlation matrix (the 1/N vs 1/(N-1) factors cancel).
  W <- sweep(x, 2L, colMeans(x), "-")
  W <- sweep(W, 2L, sqrt(colSums(W^2) / N), "/")
  Wsq <- W^2

  pairs <- utils::combn(p, 2L)
  pi <- pairs[1, ]
  pj <- pairs[2, ]
  r <- (crossprod(W) / N)[cbind(pi, pj)]

  # One influence column per off-diagonal pair (combn order), with the 1/N folded in so the
  # crossprod lands on the variance scale.
  IF <- (W[, pi, drop = FALSE] * W[, pj, drop = FALSE] -
           0.5 * (Wsq[, pi, drop = FALSE] + Wsq[, pj, drop = FALSE]) * rep(r, each = N)) / N
  Gamma <- crossprod(IF)
  Gamma <- (Gamma + t(Gamma)) / 2          # symmetrise away round-off

  if (!is.null(colnames(x))) {
    labels <- .pair_labels(colnames(x))
    dimnames(Gamma) <- list(labels, labels)
  }
  Gamma
}

# The `use` policies that listwise-delete incomplete rows (so N is the number of
# complete cases and resampling/correlation run on complete data).
.is_listwise_use <- function(use) use %in% c("complete.obs", "na.or.complete")

# Detect or compute the correlation matrix, check invertibility, and smooth a
# non-positive-definite matrix. Assumes `x` has already been validated as a
# matrix or data frame (see .assert_cor_input()). Shared by EFA and the
# suitability/retention functions so these checks live in one place. The flags
# reproduce each function's specifics (whether N is needed, the wording of the
# singular message, the SPSS positive-definite abort, etc.).
.prepare_cor_input <- function(x,
                               N = NA_real_,
                               use = "pairwise.complete.obs",
                               cor_method = "pearson",
                               N_policy = c("optional", "none", "required"),
                               acov = c("none", "diag", "full"),
                               dwls = FALSE,
                               inform_from_data = TRUE,
                               check_singular = TRUE,
                               posdef_abort = FALSE,
                               singular_tail = "no further analyses are performed",
                               N_required_msg = c(
                                 "{.arg N} is {.val NA} but a correlation matrix was entered.",
                                 "i" = "Provide {.arg N} or raw data."),
                               error_call = rlang::caller_env()) {

  N_policy <- match.arg(N_policy)
  acov <- match.arg(acov)

  # DWLS weights (1 / diag(Gamma)); populated only on the polychoric raw-data path when
  # an asymptotic covariance is requested, NULL otherwise.
  weights <- NULL
  # Full p(p-1)/2 x p(p-1)/2 asymptotic covariance of the off-diagonal correlations
  # (Var(rho-hat) scale); populated only when acov = "full" (the sandwich-SE path), NULL
  # otherwise. It is the meat of the robust/sandwich standard errors.
  Gamma <- NULL
  # The EM-estimated saturated mean/covariance and saturated log-likelihood; populated only
  # on the two-stage FIML raw-data path (cor_method = "fiml"), NULL otherwise.
  fiml <- NULL

  # TRUE exactly when an asymptotic covariance forces listwise deletion of incomplete rows:
  # the polychoric path deletes for any acov (inside .polychoric()), the Pearson/rank path only
  # for the full ADF covariance. Used both to report the override and to resolve N, so neither
  # claims a deletion that did not happen.
  acov_listwise <- acov == "full" || (.is_poly_cor(cor_method) && acov != "none")

  is_cormat <- .is_cormat(x)

  if (is_cormat) {

    R <- x

    if (N_policy == "required" && is.na(N)) {
      cli::cli_abort(N_required_msg, class = "efa_n_required", call = error_call)
    }

    # A correlation matrix carries no raw data, so an asymptotic covariance cannot be estimated
    # from it; report that the request is ignored rather than silently returning no weights.
    if (acov != "none") {
      cli::cli_warn(
        c("An asymptotic covariance was requested but {.arg x} is already a correlation matrix.",
          "i" = "{.arg acov} needs raw data and is ignored here."),
        class = "efa_acov_ignored"
      )
    }

  } else {

    if (inform_from_data) {
      cli::cli_inform(
        c("i" = "{.arg x} is not a correlation matrix; computing correlations from the raw data."),
        class = "efa_cor_from_data"
      )
    }

    if (N_policy != "none" && !is.na(N)) {
      cli::cli_warn(
        c("Both {.arg N} and raw data were supplied.",
          "i" = "Taking {.arg N} from the data."),
        class = "efa_n_from_data"
      )
    }

    # The asymptotic covariance is available from continuous (Pearson) data only as the full ADF
    # meat: a diagonal (DWLS-weight) covariance and any rank-correlation covariance are ordinal
    # constructs supplied on the polychoric path, not here. Reject the unsupported combinations
    # before reporting any listwise deletion below, so an override is not announced for a
    # computation that is then refused. EFA() gates these earlier; the guard keeps the helper's
    # contract local for any other caller.
    if (!.is_poly_cor(cor_method) && acov != "none" &&
        !(acov == "full" && cor_method == "pearson")) {
      cli::cli_abort(
        c("That asymptotic covariance is not available for {.code cor_method = {.val {cor_method}}}.",
          "x" = "Only {.code acov = \"full\"} with {.code cor_method = \"pearson\"} is supported for continuous data; {.code acov = \"diag\"} (DWLS weights) needs {.code cor_method = \"poly\"} or {.code \"tetra\"}."),
        class = "efa_acov_unsupported", call = error_call)
    }

    # An asymptotic covariance (the DWLS weights or the sandwich meat) must describe a single
    # set of complete cases, so incomplete rows are listwise-deleted below even though `use`
    # asks for pairwise-complete estimation. Surface that override rather than letting it change
    # the estimates silently.
    if (acov_listwise && use == "pairwise.complete.obs" && anyNA(x)) {
      cli::cli_inform(
        c("i" = "An asymptotic covariance requires complete cases; incomplete rows were dropped (listwise), overriding {.code use = \"pairwise.complete.obs\"}."),
        class = "efa_acov_listwise"
      )
    }

    if (cor_method == "fiml") {

      # Two-stage / full-information ML: EM-estimate the saturated multivariate-normal mean
      # and covariance from the raw data with missing values (assuming the data are missing
      # at random; Yuan, Marshall, & Bentler, 2002; Little & Rubin, 2002), then analyse the
      # standardised covariance. Every case contributes, so `use` does not apply and N is the
      # EM case count (resolved below). Routed apart from stats::cor(), which would instead
      # delete cases. The EM moments are carried forward in the returned list for downstream use.
      em <- .fiml_em_moments(x)
      R <- stats::cov2cor(em$sigma)
      # Mirror the other raw-data paths: label R by the input's column names (NULL for an
      # unnamed matrix) rather than the V1..Vp fallback the EM synthesises, so the analysis
      # matrix is named consistently across cor_method.
      dimnames(R) <- list(colnames(x), colnames(x))
      fiml <- list(mu = em$mu, sigma = em$sigma, logl = em$logl)
      # The saturated mean/covariance and saturated log-likelihood are carried forward so the
      # downstream fit indices can form the FIML likelihood-ratio chi-square (model and
      # independence-baseline log-likelihoods over the missingness patterns) rather than a naive
      # complete-case discrepancy on R; see .gof_fiml_chisq().

    } else if (.is_poly_cor(cor_method)) {

      # Polychoric/tetrachoric correlations come from the raw ordinal data.
      # .polychoric() handles missing data per pair, so reproduce the other `use`
      # policies on the data first to keep `use` meaning the same as for
      # stats::cor(): missing values abort under "all.obs"/"everything" (which do
      # not delete and so would yield an uncomputable result, matching how
      # stats::cor() errors or returns NAs there), listwise deletion applies under
      # "complete.obs"/"na.or.complete", and "pairwise.complete.obs" maps to the
      # per-pair handling. The wrapper owns the ordinal validation and classed
      # conditions (including the NA-pair abort); nearest_pd = FALSE so a
      # non-positive-definite result is smoothed once below by the shared
      # psych::cor.smooth() step rather than projected here as well. The returned
      # matrix already carries the variable names, so the colnames step that the
      # Pearson branch needs is not repeated.
      if (use %in% c("all.obs", "everything") && anyNA(x)) {
        cli::cli_abort(
          c("The correlation matrix could not be computed from the raw data because of missing values.",
            "i" = "Adjust {.arg use} (e.g. {.val pairwise.complete.obs}) or supply data with fewer missing values."),
          class = "efa_cor_na", call = error_call)
      }
      if (.is_listwise_use(use)) {
        x <- x[stats::complete.cases(x), , drop = FALSE]
        # Listwise deletion can remove every row; report that as missing data
        # rather than letting .polychoric() raise a misleading "constant variable".
        if (nrow(x) == 0L) {
          cli::cli_abort(
            c("The correlation matrix could not be computed from the raw data because of missing values.",
              "i" = "No complete cases remain under {.code use = {.val {use}}}; try {.val pairwise.complete.obs}."),
            class = "efa_cor_na", call = error_call)
        }
      }
      poly <- .polychoric(x, nearest_pd = FALSE,
                          binary_only = cor_method == "tetra",
                          acov = acov, error_call = error_call)
      R <- poly$R
      # When an asymptotic covariance is requested (the DWLS path), .polychoric()
      # estimates the matrix and the covariance on the listwise-complete rows, so the
      # per-element variances pair with the returned matrix; turn them into the symmetric
      # DWLS weight matrix. The weights are the asymptotic variances of the observed
      # (un-projected) polychoric correlations, as in lavaan, and are kept on that scale even
      # if R is subsequently projected to positive definiteness below.
      if (acov != "none") {
        # Screen the asymptotic variances once, before either consumer uses them, so a
        # degenerate pair is reported identically whether it reaches DWLS or the sandwich.
        # The diagonal is the same quantity at both levels (diag(Gamma) == the "diag" vector),
        # and .polychoric() has already labelled it by variable pair (diag() keeps those names).
        .warn_acov_degenerate(if (acov == "diag") poly$acov else diag(poly$acov))
      }
      if (acov == "diag") {
        weights <- .poly_weight_matrix(poly$acov, ncol(R))
      } else if (acov == "full") {
        # poly$acov is the full p(p-1)/2 x p(p-1)/2 asymptotic covariance: the sandwich meat.
        # Only DWLS needs the per-element inverse-variance weights, and its diagonal supplies
        # them (diag(Gamma) == the acov = "diag" vector); the ML / ULS sandwich path uses Gamma
        # alone, so the weights are built only when the estimator is DWLS.
        Gamma <- poly$acov
        if (dwls) weights <- .poly_weight_matrix(diag(Gamma), ncol(R))
      }

    } else {

      # A full ADF asymptotic covariance (the continuous sandwich meat) must describe the same
      # cases the correlation matrix was computed from -- a sandwich covariance is only valid for
      # the estimator that produced the estimates -- so reduce to the listwise-complete rows
      # first, mirroring the polychoric acov path. Without an ACOV the matrix stays pairwise for
      # data efficiency. (Only the full level is requested on the Pearson path; the diagonal DWLS
      # weights are an ordinal construct, so "diag" is left to the polychoric branch above.)
      if (acov == "full") {
        x <- x[stats::complete.cases(x), , drop = FALSE]
        if (nrow(x) < 2L) {
          cli::cli_abort(
            c("An asymptotic-distribution-free covariance needs at least two listwise-complete observations.",
              "x" = "{nrow(x)} row{?s} {?is/are} complete across all variables.",
              "i" = "Supply complete data, or impute the missing values, before requesting an {.arg acov}."),
            class = "efa_cor_no_complete_cases", call = error_call)
        }
      }

      # Missing values can make stats::cor() either throw a hard base error (e.g.
      # use = "all.obs", or "complete.obs" with no complete cases) or return NAs
      # (e.g. a column with no complete pairs under the chosen `use`). Catch both
      # instead of failing with an opaque base error here or later in
      # solve()/eigen(). A try-error without any NAs in the data has another cause
      # (e.g. a non-numeric or zero-variance column), so report that separately
      # rather than blaming missing values.
      R <- try(stats::cor(x, use = use, method = cor_method), silent = TRUE)
      if (inherits(R, "try-error") || anyNA(R)) {
        if (anyNA(x)) {
          cli::cli_abort(
            c("The correlation matrix could not be computed from the raw data because of missing values.",
              "i" = "Adjust {.arg use} (e.g. {.val pairwise.complete.obs}) or supply data with fewer missing values."),
            class = "efa_cor_na", call = error_call)
        }
        cli::cli_abort(
          c("The correlation matrix could not be computed from the raw data.",
            "i" = "Check that all columns are numeric and have non-zero variance."),
          class = "efa_cor_uncomputable", call = error_call)
      }
      colnames(R) <- colnames(x)

      # Full ADF covariance of the off-diagonal correlations (Browne, 1984): the meat of the
      # continuous robust/sandwich SEs, on the variance scale Var(rho-hat) and in utils::combn()
      # order, computed on the same listwise-complete rows as R.
      if (acov == "full") {
        Gamma <- .adf_gamma(x)
      }

    }

    if (N_policy != "none") {
      # Under listwise deletion stats::cor() drops incomplete rows, so N must be
      # the number of complete cases rather than the raw row count. Requesting a
      # polychoric ACOV (the DWLS path) likewise reduces to the listwise-complete rows
      # inside .polychoric(), regardless of `use`, so N follows the complete cases there
      # too. Under pairwise-complete estimation with missing data, N is the nominal row
      # count: each correlation is estimated from its own pairwise-complete subset, which
      # may be smaller, so fit statistics and analytic standard errors that scale with N
      # treat the data as if complete.
      N <- if (cor_method == "fiml") {
        # FIML analyses every case that carries at least one observed value (the EM case
        # count), independent of `use`.
        em$n
      } else if (.is_listwise_use(use) || acov_listwise) {
        sum(stats::complete.cases(x))
      } else {
        nrow(x)
      }
    }

  }

  # Check if the correlation matrix is invertible, if it is not, stop with message
  if (check_singular &&
      inherits(try(solve(R), silent = TRUE), "try-error")) {
    cli::cli_abort("The correlation matrix is singular; {singular_tail}.",
                   class = "efa_cor_singular", call = error_call)
  }

  # Check if correlation matrix is positive definite, if it is not, either stop
  # (SPSS type) or smooth the matrix and surface a single classed warning.
  # The threshold matches psych::cor.smooth()'s own trigger (smallest eigenvalue
  # below .Machine$double.eps), so a matrix that has already been smoothed - whose
  # eigenvalue floor sits well above this - is not re-flagged on each downstream
  # call (e.g. inside HULL -> PARALLEL -> EFA).
  if (any(eigen(R, symmetric = TRUE, only.values = TRUE)$values <
          .Machine$double.eps)) {

    if (posdef_abort) {
      cli::cli_abort(
        "The correlation matrix is not positive definite; no further analyses are performed.",
        class = "efa_cor_not_posdef", call = error_call)
    }

    # Muffle only cor.smooth's routine "smoothing was done" note (re-surfaced
    # below as a classed warning); let its serious "eigen values are NA" failure
    # warning propagate so an unrepairable matrix is not silently accepted.
    R <- withCallingHandlers(
      psych::cor.smooth(R),
      warning = function(w) {
        if (grepl("smoothing was done", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
    cli::cli_warn(
      c("The correlation matrix was not positive definite; it has been smoothed.",
        "i" = "Smoothing was applied via {.fun psych::cor.smooth}; inspect the results carefully."),
      class = "efa_cor_smoothed"
    )

  }

  list(R = R, N = N, is_cormat = is_cormat, weights = weights, Gamma = Gamma,
       fiml = fiml)
}

# Polychoric/tetrachoric correlations describe the observed data only. Criteria that cannot
# honour them reject the request with one shared classed condition rather than returning a
# silently invalid answer. `why` states which property fails: by default that the criterion
# compares the data against a separately generated continuous reference (efa_cd(),
# efa_parallel(), efa_nest(), efa_hull()), and for efa_smt() that its sequential tests rest on
# a normal-theory ML chi-square (and the RMSEA and AIC derived from it), which is not valid for
# such correlations. It is a cli template evaluated here, so it can interpolate `fn`.
#
# `fn` is the name the message shows, and every caller passes its own. Only a direct call
# reaches this message -- efa_retain() screens the same combination off the registry's
# `poly_ok` flag and reports it under the criterion id instead -- so naming the function the
# user actually called is what identifies the problem, including through a superseded wrapper,
# which forwards to the successor named here.
.reject_poly_reference <- function(cor_method, fn,
                                   why = "{.fn {fn}} compares the data against a reference computed with a continuous (Pearson) correlation.",
                                   error_call = rlang::caller_env()) {
  if (.is_poly_cor(cor_method)) {
    cli::cli_abort(
      c("{.val {cor_method}} correlations are not supported by {.fn {fn}}.",
        "x" = why,
        "i" = "Use {.val pearson}, {.val spearman}, or {.val kendall}."),
      class = "efa_cor_method_unsupported", call = error_call)
  }
  invisible(NULL)
}
