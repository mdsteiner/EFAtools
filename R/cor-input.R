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
                               inform_from_data = TRUE,
                               check_singular = TRUE,
                               posdef_abort = FALSE,
                               singular_tail = "no further analyses are performed",
                               N_required_msg = c(
                                 "{.arg N} is {.val NA} but a correlation matrix was entered.",
                                 "i" = "Provide {.arg N} or raw data."),
                               error_call = rlang::caller_env()) {

  N_policy <- match.arg(N_policy)

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

    if (N_policy != "none") {
      # Under listwise deletion stats::cor() drops incomplete rows, so N must be
      # the number of complete cases rather than the raw row count.
      N <- if (use %in% c("complete.obs", "na.or.complete")) {
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

  list(R = R, N = N, is_cormat = is_cormat)
}
