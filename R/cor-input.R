# Detect or build the analysis correlation matrix, then validate it once: assert the
# input type, smooth a non-positive-definite matrix, and resolve N. Shared by EFA and
# the suitability/retention functions so these checks live in one place.

# Checks if x is a correlation matrix
.is_cormat <- function(x){

  if(nrow(x) == ncol(x) &&
     all(x >= (-1 + .Machine$double.eps * 100), na.rm = TRUE) &&
     all(x <= (1 + .Machine$double.eps * 100), na.rm = TRUE)){

    if (round(sum(diag(x), na.rm = TRUE)) == nrow(x) && isSymmetric(unclass(unname(x)))) {

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
  invisible(x)
}

# Single source of truth for the correlation methods that are estimated from the
# raw data by .polychoric() rather than stats::cor(). Used both to route the
# computation and to reject criteria whose reference data are continuous.
.is_poly_cor <- function(cor_method) cor_method %in% c("poly", "tetra")

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
      c("A polychoric asymptotic variance was not positive, so the DWLS weights are undefined.",
        "i" = "A variable pair is (near-)degenerate."),
      class = "efa_dwls_degenerate_weight"
    )
  }
  W <- matrix(0, p, p)
  idx <- utils::combn(p, 2L)
  W[t(idx)] <- 1 / acov_diag
  W + t(W)
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
                               acov = c("none", "diag"),
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

  is_cormat <- .is_cormat(x)

  if (is_cormat) {

    R <- x

    if (N_policy == "required" && is.na(N)) {
      cli::cli_abort(N_required_msg, class = "efa_n_required", call = error_call)
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

    if (.is_poly_cor(cor_method)) {

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
      poly <- .polychoric(x, n_threads = 1L, nearest_pd = FALSE,
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
        weights <- .poly_weight_matrix(poly$acov, ncol(R))
      }

    } else {

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

    }

    if (N_policy != "none") {
      # Under listwise deletion stats::cor() drops incomplete rows, so N must be
      # the number of complete cases rather than the raw row count. Requesting a
      # polychoric ACOV (the DWLS path) likewise reduces to the listwise-complete rows
      # inside .polychoric(), regardless of `use`, so N follows the complete cases there
      # too.
      N <- if (.is_listwise_use(use) || acov != "none") {
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

  list(R = R, N = N, is_cormat = is_cormat, weights = weights)
}

# Polychoric/tetrachoric correlations describe the observed data only. Criteria
# that compare the data against a separately generated continuous reference (CD,
# PARALLEL, NEST) cannot honour them, so they reject the request with one shared
# classed condition rather than silently mixing an ordinal observed matrix with a
# Pearson reference distribution.
.reject_poly_reference <- function(cor_method, fn, error_call = rlang::caller_env()) {
  if (.is_poly_cor(cor_method)) {
    cli::cli_abort(
      c("{.val {cor_method}} correlations are not supported by {.fn {fn}}.",
        "x" = "{.fn {fn}} compares the data against a reference computed with a continuous (Pearson) correlation.",
        "i" = "Use {.val pearson}, {.val spearman}, or {.val kendall}."),
      class = "efa_cor_method_unsupported", call = error_call)
  }
  invisible(NULL)
}
