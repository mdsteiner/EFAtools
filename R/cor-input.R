# Detect or build the analysis correlation matrix, then validate it once: assert the
# input type, smooth a non-positive-definite matrix, and resolve N. Shared by EFA and
# the suitability/retention functions so these checks live in one place.

# Unit diagonal on one shared tolerance: the correlation-vs-covariance classification
# hinges on .is_cormat() and .is_covmat() applying the SAME test, so any matrix falls on
# exactly one side of it.
.has_unit_diag <- function(x, na.rm = FALSE) {
  all(abs(diag(x) - 1) <= .Machine$double.eps * 100, na.rm = na.rm)
}

# Everything that identifies a correlation matrix except symmetry: at least two columns,
# square, all-numeric, a unit diagonal, and every entry in [-1, 1]. Two callers need exactly
# this test -- .is_cormat(), which adds the symmetry requirement, and .assert_cor_input(),
# which rejects an input that meets it but is asymmetric -- so it lives in one place and the
# two cannot drift apart. Returns the matrix form of `x` when the shape matches and NULL
# otherwise, so the caller neither re-tests nor re-coerces.
.cormat_shape <- function(x) {

  # Cheapest test first, and before any coercion, for the reason .is_covmat() documents: raw
  # data is almost never square, and as.matrix() on a data frame copies the whole thing. The
  # two-column minimum matches .is_covmat()'s own guard and keeps a 1x1 input -- not an
  # analysable correlation matrix -- on the raw-data route, where it is rejected.
  if (ncol(x) < 2L || nrow(x) != ncol(x)) return(NULL)

  # A data frame is coerced here rather than tested in place: diag() on a data frame is a hard
  # base error, and isSymmetric() needs a real matrix too. The all-numeric test guards both the
  # coercion and the range comparisons that follow -- on a factor column the comparisons would
  # dispatch to Ops.factor and fail with an unclassed base error. Anything non-numeric is not a
  # correlation matrix, so it falls through to the raw-data path and is rejected there with a
  # classed condition.
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) return(NULL)
    x <- as.matrix(x)
  }
  if (!is.numeric(x)) return(NULL)

  if (!all(x >= (-1 - .Machine$double.eps * 100), na.rm = TRUE) ||
      !all(x <= (1 + .Machine$double.eps * 100), na.rm = TRUE)) return(NULL)

  if (!.has_unit_diag(x, na.rm = TRUE)) return(NULL)

  x
}

# Heuristically decide whether x is already a correlation matrix rather than raw data:
# .cormat_shape() plus symmetry. This cannot distinguish a correlation matrix from a (rare)
# raw-data matrix that happens to meet all of these conditions, and it treats a covariance
# matrix as raw data (its diagonal is not all ones); callers that accept either input assume
# raw data unless these conditions hold.
.is_cormat <- function(x){

  xm <- .cormat_shape(x)

  if (is.null(xm) || !isSymmetric(unclass(unname(xm)))) return(FALSE)

  if (any(is.na(xm))) {

    cli::cli_abort(
      c("{.arg x} looks like a correlation matrix but contains missing values.",
        "i" = "Please check the entered data."),
      class = "efa_cormat_has_na"
    )

  }

  TRUE

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
  # A square, unit-diagonal, in-range matrix that is not symmetric is a correlation matrix
  # transcribed one triangle at a time, or typed at too few decimals to mirror to machine
  # precision (isSymmetric()'s tolerance is 100 * .Machine$double.eps). Neither classifier
  # accepts it, so passed on it would be fed to stats::cor() as if its p columns were p cases
  # -- singular by construction, and reported as such, which sends the user looking for
  # collinearity that does not exist. It is deliberately not symmetrised here: with one
  # triangle left empty the other is authoritative, but an asymmetry from a transcription slip
  # leaves neither triangle authoritative, and repairing it either way would analyse a matrix
  # the user never entered.
  xm <- .cormat_shape(x)
  if (!is.null(xm) && !isSymmetric(unclass(unname(xm)))) {
    sym_lead <- if (raw_only) {
      "{.arg x} looks like a (non-symmetric) correlation matrix, not a data frame/matrix of raw data."
    } else {
      "{.arg x} looks like a correlation matrix but is not symmetric."
    }
    # The unit-diagonal and range tests run with na.rm = TRUE, so a matrix whose diagonal or
    # whose unused triangle was left empty reaches this branch; do not claim of missing cells
    # that they are in range.
    sym_detail <- if (anyNA(xm)) {
      "Every entry that is present lies in [-1, 1] and no diagonal entry departs from one, but part of the matrix is missing and {.code x[i, j]} and {.code x[j, i]} differ."
    } else {
      "Its diagonal is all ones and every entry lies in [-1, 1], but {.code x[i, j]} and {.code x[j, i]} differ."
    }
    # Which triangle to mirror has to be read off the data: the hint for a transcribed lower
    # triangle overwrites the entries of a transcribed upper one, which would leave a matrix
    # that passes every later check and analyses as the identity. An unused triangle is left
    # either empty (NA) or at zero. Mirroring is offered only when exactly one triangle was
    # filled in: when both carry entries neither is authoritative (the same reasoning that
    # keeps this branch from symmetrising on its own), and when neither does there is nothing
    # to mirror -- copying an empty triangle would also produce the identity.
    tri_l <- xm[lower.tri(xm)]
    tri_u <- xm[upper.tri(xm)]
    empty <- function(tri) all(is.na(tri) | tri == 0)
    sym_remedy <- if (raw_only) {
      "Supply the raw observations the correlation matrix was computed from."
    } else if (empty(tri_u) && !empty(tri_l)) {
      "Only the lower triangle carries entries; mirror it: {.code x[upper.tri(x)] <- t(x)[upper.tri(x)]}."
    } else if (empty(tri_l) && !empty(tri_u)) {
      "Only the upper triangle carries entries; mirror it: {.code x[lower.tri(x)] <- t(x)[lower.tri(x)]}."
    } else if (empty(tri_l) && empty(tri_u)) {
      "Neither triangle carries a correlation, so there is nothing to mirror; enter the off-diagonal correlations."
    } else {
      "Both triangles carry entries but they disagree, so neither can be mirrored onto the other; check the entries of the pairs that differ."
    }
    cli::cli_abort(
      c(sym_lead,
        "x" = sym_detail,
        "i" = sym_remedy),
      class = "efa_input_not_symmetric", call = error_call
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

# Cap a list of affected labels for a condition message at the first five plus a count of the
# rest: a heterogeneous ordinal set can have dozens of bad variable pairs, and a data set read
# in as character has one bad column per variable, so a message that prints them all is
# unreadable. Returns the pieces rather than a finished string so each caller keeps its own cli
# template. When the list is truncated, `shown` is a cli_vec whose "and" is dropped so the
# trailing count closes the enumeration instead of adding a second conjunction. Shared by every
# condition that enumerates affected variables -- the boundary and sparse-cell diagnostics, the
# asymptotic covariance screen, the DWLS weight abort, and the uncomputable-correlation column
# diagnosis -- so their wording cannot drift apart.
.cap_label_list <- function(labels) {
  shown <- utils::head(labels, 5L)
  more <- length(labels) - length(shown)
  if (more == 0L) return(list(shown = shown, rest = ""))
  list(shown = cli::cli_vec(shown, style = list("vec-last" = ", ")),
       rest = paste0(", and ", more, " more"))
}

# Build the symmetric per-element DWLS weight matrix W (W_ij = 1 / Var(rho_hat_ij), zero
# diagonal) from the diagonal polychoric asymptotic covariance returned by .polychoric():
# a length p(p - 1)/2 vector of off-diagonal variances ordered by the upper triangle
# (i < j) in utils::combn() column order. The combn order is row-major, which differs
# from R's column-major upper.tri(), so the pairs are placed via explicit (i, j) indexing
# rather than `W[upper.tri(W)] <- ...`. Shared by the DWLS point estimate and the per-
# replicate bootstrap recompute.
.poly_weight_matrix <- function(acov_diag, p) {
  # A non-positive or non-finite asymptotic variance means a variable pair whose inverse-variance
  # weight is undefined; refuse rather than emit Inf/NaN weights that would silently corrupt the
  # fit. NA arrives for a pair estimated at the boundary of the parameter space, which has no
  # asymptotic variance at all (see .polychoric()), so this refusal is deterministic rather than
  # dependent on which side of zero a floating-point artefact happened to land. (The bootstrap
  # recompute catches this and drops the replicate.)
  bad <- !is.finite(acov_diag) | acov_diag <= 0
  if (any(bad)) {
    n_bad <- sum(bad)
    pairs <- if (is.null(names(acov_diag))) paste0("pair ", which(bad)) else names(acov_diag)[bad]
    cap <- .cap_label_list(pairs)
    cli::cli_abort(
      c("DWLS needs an inverse-variance weight for every variable pair, but {n_bad} pair{?s} {?has/have} no usable asymptotic variance.",
        "x" = "Affected {cli::qty(n_bad)}pair{?s}: {.val {cap$shown}}{cap$rest}.",
        "i" = "The most common cause with binary variables is a response combination that never occurs: the correlation of such a pair is estimated with a continuity correction, which repairs the estimate but leaves it without a variance.",
        "i" = "Otherwise the pair carries too little information about its correlation to have one: its responses are perfectly ordered (never observed in reversed order), so the correlation sits on the boundary of the parameter space, or a response category is (near-)empty.",
        "i" = "Fit with {.code estimator = \"ULS\"} instead, which needs no such weights. Collapsing rare response categories in the variables involved, or dropping one item of a redundant pair, also resolves the cases that are genuinely under-identified."),
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
# takes the same branch on the same data. A variable pair whose asymptotic variance is unusable
# leaves the DWLS weight meaningless and the sandwich standard errors of that pair meaningless
# rather than merely large, and neither outcome should be silent, so warn once and name the pairs.
#
# The main case is a pair estimated at the boundary of the parameter space (a contingency table at
# a Frechet bound, i.e. perfectly ordered), which has no asymptotic variance at all and reaches
# this screen as NA by construction -- see .polychoric(). That is caught by the non-finite clause,
# so the branch taken here does not depend on the platform. The `> 1` clause is a residual net for
# any remaining pair whose variance is finite but useless: a variance above 1 already implies a
# +/- 1 SE interval as wide as the whole [-1, 1] range of a correlation. (Searches over
# strongly-correlated tables found no table off a Frechet bound that reaches it -- the largest was
# ~0.008 -- so it is expected to be inert, but it costs nothing to keep.) Non-positive values are
# degenerate a fortiori.
.warn_acov_degenerate <- function(acov_diag, labels = names(acov_diag)) {
  bad <- !is.finite(acov_diag) | acov_diag <= 0 | acov_diag > 1
  if (!any(bad)) return(invisible(NULL))

  pairs <- if (is.null(labels)) paste0("pair ", which(bad)) else labels[bad]
  cap <- .cap_label_list(pairs)

  cli::cli_warn(
    c("The polychoric asymptotic covariance is unavailable for {cli::qty(pairs)} variable pair{?s} {.val {cap$shown}}{cap$rest}.",
      "x" = "{cli::qty(pairs)} {?This pair was/These pairs were} either estimated with a continuity correction for a response combination that never occurs -- which repairs the correlation but leaves it without a variance -- or {?is/are} not identified well enough to have a usable one (e.g. a perfectly ordered or near-empty response table).",
      "i" = "Any robust/sandwich standard error involving {cli::qty(pairs)} {?it/them} is withheld. A DWLS fit refuses the data outright when the variance is missing or non-positive, and down-weights {cli::qty(pairs)} {?it/them} out of the solution when it is merely far too large.",
      "i" = "Fitting with {.code estimator = \"ULS\"} avoids the weights entirely; collapsing rare response categories in the variables involved resolves the under-identified cases."),
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

# Project a symmetric matrix onto the positive definite correlation matrices: floor the
# eigenvalues, rescale the spectrum to keep the trace at p, rebuild, and standardise back to a
# unit diagonal. Same construction and same floor as psych::cor.smooth(), so the two agree
# wherever both act; used where cor.smooth() declines to act on a matrix that is still not
# positive definite (see .prepare_cor_input()).
.project_cor_pd <- function(R) {
  eig_tol <- 1e-12                       # cor.smooth()'s own `eig.tol`, and its 100x floor
  e <- eigen(R, symmetric = TRUE)
  ev <- e$values
  ev[ev < eig_tol] <- 100 * eig_tol
  ev <- ev * ncol(R) / sum(ev)
  # ev * t(V) scales row i of t(V) by ev[i], i.e. diag(ev) %*% t(V) without forming diag(ev).
  out <- stats::cov2cor(e$vectors %*% (ev * t(e$vectors)))
  dimnames(out) <- dimnames(R)
  out
}

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
                               # EM budget for the two-stage FIML moments; ignored by every other
                               # cor_method. Taken FROM the control constructor rather than
                               # restated, so the entry points that thread an `estimate_control()`
                               # (efa_fit(), and what forwards to it) and those that cannot
                               # (efa_average()) can never estimate different correlations from the
                               # same data because one default was changed and the other was not.
                               # Both are evaluated lazily, inside the FIML branch only.
                               fiml_max_iter = estimate_control()$fiml_max_iter,
                               fiml_tol = estimate_control()$fiml_tol,
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
  # TRUE once `x` has been reduced to its listwise-complete rows below, so that a failure
  # diagnosed afterwards describes the rows the correlation was actually attempted on.
  listwise_reduced <- FALSE

  # TRUE exactly when an asymptotic covariance forces listwise deletion of incomplete rows:
  # the polychoric path deletes for any acov (inside .polychoric()), the Pearson/rank path only
  # for the full ADF covariance. Used both to report the override and to resolve N, so neither
  # claims a deletion that did not happen.
  acov_listwise <- acov == "full" || (.is_poly_cor(cor_method) && acov != "none")

  is_cormat <- .is_cormat(x)

  if (is_cormat) {

    # A data frame that passed .is_cormat() is a correlation matrix; carry it forward as a
    # matrix so every downstream consumer gets the type it expects (determinant(), the
    # eigendecompositions and the compiled estimators all reject a data frame).
    R <- as.matrix(x)

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
      em <- .fiml_em_moments(x, max_iter = fiml_max_iter, tol = fiml_tol)
      R <- stats::cov2cor(em$sigma)
      # Mirror the other raw-data paths: label R by the input's column names (NULL for an
      # unnamed matrix) rather than the V1..Vp fallback the EM synthesises, so the analysis
      # matrix is named consistently across cor_method.
      dimnames(R) <- list(colnames(x), colnames(x))
      fiml <- list(mu = em$mu, sigma = em$sigma, logl = em$logl,
                   converged = em$converged, iter = em$iter,
                   n_patterns = em$n_patterns, n = em$n)
      # The saturated mean/covariance and saturated log-likelihood are carried forward so the
      # downstream fit indices can form the FIML likelihood-ratio chi-square (model and
      # independence-baseline log-likelihoods over the missingness patterns) rather than a naive
      # complete-case discrepancy on R; see .gof_fiml_chisq(). The Stage-1 diagnostics travel
      # with them so the caller can report that the analysed matrix is the EM's last iterate
      # rather than the converged FIML estimate -- the warning alone leaves no record on the
      # fitted object, and callers that suppress per-fit warnings would surface nothing at all.

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
        # Recorded so the diagnosis below can say that it judged the rows actually used: a
        # column can be constant among the listwise-complete rows and vary in the full data.
        listwise_reduced <- TRUE
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
      #
      # stats::cor() raises exactly one warning of its own, "the standard deviation is zero"
      # for a constant column, and it always leaves NAs in the result -- so it is always
      # followed by the classed abort below, which names that column. Muffle it rather than
      # letting an untranslated base condition print ahead of, and in competition with, the
      # package's own diagnosis.
      R <- withCallingHandlers(
        try(stats::cor(x, use = use, method = cor_method), silent = TRUE),
        warning = function(w) invokeRestart("muffleWarning"))
      if (inherits(R, "try-error") || anyNA(R)) {
        if (anyNA(x)) {
          cli::cli_abort(
            c("The correlation matrix could not be computed from the raw data because of missing values.",
              "i" = "Adjust {.arg use} (e.g. {.val pairwise.complete.obs}) or supply data with fewer missing values."),
            class = "efa_cor_na", call = error_call)
        }
        # Name the columns rather than leaving the user to find them in a wide data set. The
        # three causes have different remedies, so they are reported separately, and any of
        # them can apply to several columns at once (a whole data set read in as character,
        # say). Usability is the test stats::cor() itself applies, is.numeric() OR
        # is.logical(), so dichotomous items stored as TRUE/FALSE -- which cor() correlates
        # without complaint -- are not misreported as non-numeric.
        #
        # Both remaining tests are written to yield FALSE rather than NA on a column they
        # cannot judge: sd() is NA for fewer than two values and NaN for a column carrying an
        # infinity, and a single NA in `const` would make `any(const)` NA and abort with an
        # unclassed base error in place of the classed condition below. Such a column falls
        # through to the generic bullet instead, which is what this branch reported before.
        nms <- colnames(x)
        if (is.null(nms)) nms <- paste0("V", seq_len(ncol(x)))   # unnamed input: by position
        usable <- if (is.data.frame(x)) {
          vapply(x, function(v) is.numeric(v) || is.logical(v), logical(1))
        } else {
          rep(is.numeric(x) || is.logical(x), ncol(x))
        }
        nonfin <- logical(ncol(x))
        const <- logical(ncol(x))
        for (j in which(usable)) {
          col_j <- if (is.data.frame(x)) x[[j]] else x[, j]
          if (!all(is.finite(col_j))) {
            nonfin[j] <- TRUE
          } else {
            # The exact test stats::cor() applies internally: sd == 0, not a tolerance.
            const[j] <- isTRUE(stats::sd(col_j) == 0)
          }
        }

        # The name lists are capped through the same helper as every other condition in this
        # file that enumerates affected variables, so the two cannot render differently. The
        # cli::qty() is repeated after the capped list because the empty `rest` string is
        # itself a substitution and would otherwise set the pluralisation quantity to one.
        bullets <- "The correlation matrix could not be computed from the raw data."
        if (any(!usable)) {
          n_type <- sum(!usable)
          cap_type <- .cap_label_list(nms[!usable])
          bullets <- c(
            bullets,
            "x" = "{cli::qty(n_type)}Column{?s} {.val {cap_type$shown}}{cap_type$rest}{cli::qty(n_type)} {?is/are} not numeric.",
            "i" = "Ordinal items stored as factors or character strings are correlated with {.code cor_method = \"poly\"} ({.val tetra} for binary items); anything else has to be converted or dropped."
          )
        }
        if (any(nonfin)) {
          n_inf <- sum(nonfin)
          cap_inf <- .cap_label_list(nms[nonfin])
          bullets <- c(
            bullets,
            "x" = "{cli::qty(n_inf)}Column{?s} {.val {cap_inf$shown}}{cap_inf$rest}{cli::qty(n_inf)} {?contains/contain} infinite values.",
            "i" = "An infinite value has no correlation with anything; check for a division by zero or an out-of-range missing-value code."
          )
        }
        if (any(const)) {
          n_const <- sum(const)
          cap_const <- .cap_label_list(nms[const])
          # Under a listwise-reducing `acov` the variance was judged on the complete cases
          # only, so a named column can still vary in the data the user supplied; say which
          # rows the verdict describes rather than advising a drop that may be wrong.
          const_remedy <- if (listwise_reduced) {
            "The variance was computed on the listwise-complete rows an {.arg acov} needs, which can be far fewer than the data supplied; either supply more complete cases or drop the constant {cli::qty(n_const)}column{?s}."
          } else {
            "A constant variable correlates with nothing; drop {cli::qty(n_const)}{?it/them} before the analysis."
          }
          bullets <- c(
            bullets,
            "x" = "{cli::qty(n_const)}Column{?s} {.val {cap_const$shown}}{cap_const$rest}{cli::qty(n_const)} {?has/have} zero variance.",
            "i" = const_remedy
          )
        }
        # Derived from what was actually reported, so adding a fourth cause above cannot leave
        # the fallback firing alongside a named one.
        if (length(bullets) == 1L) {
          too_few_rows <- if (nrow(x) < 2L) {
            " Correlations also need at least two observations."
          } else {
            ""
          }
          bullets <- c(
            bullets,
            "i" = paste0("Check that all columns are numeric and have non-zero variance.",
                         too_few_rows)
          )
        }
        cli::cli_abort(bullets, class = "efa_cor_uncomputable", call = error_call)
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

  # One symmetric eigendecomposition answers both checks below, so no separate p x p inverse
  # is built and discarded as a singularity test. For a symmetric matrix
  # min(|lambda|) / max(|lambda|) is exactly the reciprocal 2-norm condition number, and
  # min(lambda) < .Machine$double.eps is the positive-definiteness test. max(|lambda|) >= 1
  # here because R always carries a unit diagonal, so the ratio is well defined.
  #
  # The singularity threshold is p * .Machine$double.eps, the standard numerical-rank
  # tolerance (as in MASS::ginv and LAPACK's own rank estimates), NOT a bare epsilon. A
  # backward-stable symmetric eigensolver returns eigenvalues with absolute error of order
  # eps * ||R||_2 (Golub & Van Loan, 2013, Matrix Computations, 4th ed., sec. 8.1; Higham,
  # 2002, Accuracy and Stability of Numerical Algorithms, 2nd ed., ch. 5), so for an exactly
  # rank-deficient matrix the computed min(|lambda|) IS that rounding error: testing it
  # against eps alone compares a quantity against the noise floor that produced it, and the
  # verdict becomes a coin flip. Measured on 495 rank-deficient correlation matrices
  # (p = 4 to 40; a duplicated item, an item equal to the sum of two others, two redundant
  # items), a bare epsilon missed 124 of them while p * eps missed none.
  #
  # This threshold refuses everything solve() refused, which is what the check did before:
  # kappa_1 <= p * kappa_2 for a symmetric matrix (the norm-equivalence bounds in Golub & Van
  # Loan, sec. 2.3), so solve()'s gate (reciprocal 1-norm
  # condition number below eps) implies min(|lambda|) / max(|lambda|) < p * eps. It also
  # refuses a narrow band just above it -- 47 of the 600 matrices in that sweep -- which is
  # the intended direction: a correlation matrix with a condition number above 1 / (p * eps)
  # carries no usable inverse. Well-conditioned data is nowhere near the threshold; every
  # correlation matrix shipped with the package sits at 1e13 to 1e15 times eps.
  #
  # The two checks are NOT interchangeable and both are kept, in this order: an indefinite
  # matrix can be perfectly well conditioned (it is smoothed, not refused -- only the
  # eigenvalue test sees it, which is why the singularity test takes the absolute value and
  # the smoothing test does not), and a positive definite matrix can be numerically singular
  # (it is refused before any smoothing is attempted). The smoothing threshold stays at a bare
  # eps because it is pinned to psych::cor.smooth()'s own trigger, not to a rank tolerance.
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values

  # Check if the correlation matrix is invertible, if it is not, stop with message
  if (check_singular &&
      min(abs(ev)) / max(abs(ev)) < ncol(R) * .Machine$double.eps) {
    # The offending matrix travels with the condition, so a caller that wants to name
    # the culprits (efa_screen() reports the perfectly correlated pairs) does not have
    # to recompute a correlation matrix this function has already built.
    cli::cli_abort("The correlation matrix is singular; {singular_tail}.",
                   class = "efa_cor_singular", R = R, call = error_call)
  }

  # Check if correlation matrix is positive definite, if it is not, either stop
  # (SPSS type) or smooth the matrix and surface a single classed warning.
  # The threshold matches psych::cor.smooth()'s own trigger (smallest eigenvalue
  # below .Machine$double.eps), so a matrix that has already been smoothed - whose
  # eigenvalue floor sits well above this - is not re-flagged on each downstream
  # call (e.g. inside HULL -> PARALLEL -> EFA).
  if (min(ev) < .Machine$double.eps) {

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

    # cor.smooth() decides for itself whether to act, and it decides from an eigendecomposition
    # that asks for the eigenvectors, whereas the eigenvalues above are values-only. LAPACK runs
    # a different driver for each, so their smallest eigenvalues differ by the rounding of the
    # decomposition itself (of order eps * ||R||). For a matrix that is singular to working
    # precision that rounding IS the eigenvalue, so the two verdicts disagree on a substantial
    # share of such matrices -- measured at two in five over a sweep of exactly rank-deficient
    # correlation matrices of order 5 to 40 -- and cor.smooth() then returns the matrix
    # untouched, leaving a matrix reported as smoothed still indefinite. Which side of the
    # disagreement a build lands on is a property of its LAPACK, so verify the result here and
    # project it on the same eigenvalue floor when it was left alone.
    if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <
        .Machine$double.eps) {
      R <- .project_cor_pd(R)
    }

    # The remedy is named as the eigenvalue floor rather than as cor.smooth(), because on the
    # branch above it is .project_cor_pd() that applies it.
    cli::cli_warn(
      c("The correlation matrix was not positive definite; it has been smoothed.",
        "i" = "Its smallest eigenvalues were floored, as in {.fun psych::cor.smooth}; inspect the results carefully."),
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
        # Kendall's tau is deliberately absent: it is not on the Pearson scale the factor
        # model parameterises (for bivariate-normal data tau = (2/pi) asin(rho)), so
        # offering it as the remedy to a user who correctly reached for a polychoric
        # correlation would steer an ordinal analysis onto the most attenuated metric of
        # the three. Spearman is kept: it stays a correlation coefficient.
        "i" = "Use {.val pearson} or {.val spearman}."),
      class = "efa_cor_method_unsupported", call = error_call)
  }
  invisible(NULL)
}
