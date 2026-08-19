#' Estimate factor scores and score-quality diagnostics for an EFA model
#'
#' Computes factor-score weights (and, from raw data, the factor scores
#' themselves) natively for an [efa_fit()] solution or a directly supplied loading
#' matrix, together with the score-quality diagnostics that describe how well the
#' estimated scores represent the factors: the score intercorrelations, the
#' determinacy (validity) and univocality of each score, and Guttman's
#' indeterminacy index. Factor scores are returned only when raw data are
#' supplied; a correlation matrix yields the weights and diagnostics alone.
#'
#' @details
#' The `p` by `m` weight matrix `W` (standardized scores are `scale(X) %*% W`) is
#' computed from the structure matrix `S = Lambda %*% Phi`, the model
#' uniquenesses `Psi = diag(1 - h2)`, and the scoring correlation matrix `R`
#' according to `method`:
#' \describe{
#'   \item{`"regression"`}{Thurstone's (1935) regression scores, `W = R^-1 S`.}
#'   \item{`"Bartlett"`}{Bartlett's (1937) conditionally unbiased scores.}
#'   \item{`"Anderson"`}{Anderson & Rubin's (1956) uncorrelated, unit-variance
#'     scores; defined for orthogonal factors only.}
#'   \item{`"tenBerge"`}{ten Berge, Krijnen, Wansbeek & Shapiro's (1999) scores,
#'     which preserve the factor intercorrelations.}
#'   \item{`"Harman"`}{Harman's (1976) idealized-variable scores.}
#'   \item{`"components"`}{component scores, `W = Lambda`. These are formed from
#'     the raw, uncentered data (`X %*% W`) rather than the standardized data, so
#'     unlike the other methods they are on the scale of the input variables. The
#'     diagnostics below describe the standardized combination `scale(X) %*% W`,
#'     and therefore differ from the realized correlations of the returned scores
#'     whenever the variables have unequal variances.}
#' }
#' The determinacy (validity) of a score is its correlation with the factor it
#' estimates, computed from the returned weights; for regression scores it is the
#' multiple correlation between the factor and the observed variables (Guttman,
#' 1955; Grice, 2001). The off-diagonal score-factor correlations give the
#' univocality (the correlation of a score with the *other* factors), and
#' `2 rho^2 - 1` is Guttman's (1955) indeterminacy index, the minimum correlation
#' between two equally valid sets of scores. For a `method` other than
#' `"regression"` both quantities are specific to those scores: the determinacy is
#' the method's own score-factor correlation (never larger than the regression
#' value), and the reported `guttman` follows from it.
#'
#' Determinacies close to 1 mean the scores stand in for the factor with little
#' loss; Grice (2001) regards values of about .90 and above as the level required
#' before scores are interpreted for individual cases, and treats lower values as
#' usable only for group-level research. The Guttman index makes the same point
#' more sharply, because a factor score is never the factor: at `rho = .90` two
#' equally valid sets of scores can still correlate as low as .62, and at
#' `rho = .80` as low as .28, so the rank order of cases is not unique.
#'
#' Which `method` to prefer follows from what the scores are for. Regression
#' scores correlate most highly with the factor, but they are biased towards it
#' and correlate across factors even when the model is orthogonal. Bartlett
#' scores are conditionally unbiased, which makes them the choice when the scores
#' stand in for the factor in a later model. `"tenBerge"` reproduces the factor
#' intercorrelations `Phi`, so it is the choice when the scores will be correlated
#' with each other or with external variables. `"Anderson"` forces the scores to
#' be uncorrelated with unit variance and is appropriate only when the factors
#' themselves are orthogonal. `"components"` is a weighted sum of the observed
#' variables rather than an estimate of a common factor.
#'
#' @param x data.frame or matrix. Raw data (needed to obtain factor scores) or a
#'   correlation matrix (yields weights and diagnostics only). A correlation matrix
#'   gives the correlations the weights are derived from only where `f` is a loading
#'   matrix. An [efa_fit()] object carries the correlations it was fitted on, and
#'   those are used instead, so supply `rho` to derive the weights from another
#'   matrix. `x` describes the model variables either way. When raw data carry
#'   column names, they are matched to the model variables by name (any extra
#'   columns are ignored, and a model variable missing from `x` is an error).
#'   A named correlation matrix is likewise matched to the loading rows by name;
#'   its row and column names must use the same order, and it must carry one row and
#'   column for each model variable. Unnamed input is matched
#'   by position.
#'
#'   Raw data are scored as supplied: no imputation is performed, so a case with a
#'   missing value on any model variable receives `NA` scores (and is reported as
#'   not scored). A model variable that carries no usable spread in `x` -- constant,
#'   infinite, or observed fewer than twice -- is an error.
#' @param f object of class [efa_fit()], an `efa_loadings` object, or a matrix of factor
#'   loadings.
#' @param Phi matrix. Factor intercorrelations. Only used when a loading matrix is
#'   supplied directly in `f`; taken from the `efa` object otherwise, in which
#'   case a supplied `Phi` is ignored with a warning. Named rows and columns are
#'   matched to the loading columns and must use the same order.
#'   Default is `NULL`, in which case the factors are assumed uncorrelated.
#' @param rho matrix. Correlation matrix used to derive the scoring weights.
#'   Defaults to `NULL`, in which case `f$orig_R` is used for an `efa` object; for
#'   a directly supplied loading matrix, `x` itself when it is a correlation
#'   matrix, and `cor(x, use = "pairwise")` otherwise. Pass a matrix here to score
#'   against a correlation other than the one implied by `f`/`x`. Named rows and
#'   columns are matched to the loading rows; row and column names must use the
#'   same order.
#' @param method character. The factor-score method: one of `"regression"`
#'   (default), `"Bartlett"`, `"Anderson"`, `"tenBerge"`, `"Harman"`, or
#'   `"components"`.
#'
#' @return An object of class `efa_scores`, a list containing:
#'
#' \item{weights}{The `p` by `m` factor-score weight matrix.}
#' \item{scores}{The factor scores (`n` by `m`), or `NULL` when a correlation
#'   matrix was supplied. A case with a missing value on any model variable is not
#'   scored and keeps `NA` in every column.}
#' \item{r.scores}{The `m` by `m` correlations of the factor-score estimates. Like
#'   `score_cor` and `determinacy`, these describe the standardized combination
#'   `scale(x) %*% W`, which for `method = "components"` is not the scale the
#'   returned `scores` are on.}
#' \item{score_cor}{The `m` by `m` score-factor correlation matrix; its diagonal
#'   is the determinacy (validity) of each score and its off-diagonals the
#'   univocality.}
#' \item{determinacy}{A data frame with, per factor, the determinacy `rho`, the
#'   squared determinacy `rho2`, and Guttman's indeterminacy index `guttman`.}
#' \item{settings}{A list of the settings used, including the number of supplied
#'   observations `n_obs` and the number of them that could be scored `n_scored`.}
#'
#' @source Thurstone, L. L. (1935). *The vectors of mind*. University of Chicago
#'   Press.
#'
#'   Bartlett, M. S. (1937). The statistical conception of mental factors.
#'   *British Journal of Psychology, 28*, 97-104.
#'
#'   Anderson, T. W., & Rubin, H. (1956). Statistical inference in factor
#'   analysis. In *Proceedings of the Third Berkeley Symposium on Mathematical
#'   Statistics and Probability* (Vol. 5, pp. 111-150). University of California
#'   Press.
#'
#'   Guttman, L. (1955). The determinacy of factor score matrices with
#'   implications for five other basic problems of common-factor theory.
#'   *British Journal of Statistical Psychology, 8*, 65-81.
#'
#'   ten Berge, J. M. F., Krijnen, W. P., Wansbeek, T., & Shapiro, A. (1999).
#'   Some new results on correlation-preserving factor scores prediction methods.
#'   *Linear Algebra and its Applications, 289*, 311-318.
#'
#'   Grice, J. W. (2001). Computing and evaluating factor scores. *Psychological
#'   Methods, 6*, 430-450.
#'
#' @seealso [efa_fit()] for the solution these are computed from.
#'
#' @family factor scoring
#'
#' @export
#'
#' @examples
#' # Weights and score diagnostics from an EFA on a correlation matrix
#' efa <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                estimator = "PAF", rotation = "oblimin")
#' fs <- efa_scores(test_models$baseline$cormat, f = efa)
#' fs
#' summary(fs)
#'
#' # Factor scores from raw data (Bartlett method)
#' \donttest{
#' efa_raw <- efa_fit(GRiPS_raw, n_factors = 1, estimator = "PAF")
#' efa_scores(GRiPS_raw, f = efa_raw, method = "Bartlett")
#' }
#'
#' # Loadings supplied directly, with the factor intercorrelations
#' efa_scores(test_models$baseline$cormat, f = efa$rot_loadings, Phi = efa$Phi)
#'
efa_scores <- function(x, f, Phi = NULL, rho = NULL,
                       method = c("regression", "Bartlett", "Anderson",
                                  "tenBerge", "Harman", "components")) {

  .assert_cor_input(x)

  # A correlation matrix carries no cases, so scores cannot be formed; only the
  # weights and diagnostics are returned.
  is_cmat <- .is_cormat(x)
  if (is_cmat) {
    cli::cli_inform(
      c("i" = "{.arg x} is a correlation matrix; factor scores cannot be computed. Only factor weights and score diagnostics are returned. Enter raw data to get factor scores."),
      class = "efa_scores_needs_raw"
    )
  }

  method <- .match_arg_ci(method)
  checkmate::assert_matrix(Phi, null.ok = TRUE)

  if (!inherits(f, c("EFA", "matrix", "LOADINGS"))) {
    cli::cli_abort(
      "{.arg f} must be an {.cls efa} object (from {.fn efa_fit}), a matrix, or an {.cls efa_loadings} object.",
      class = "efa_scores_bad_f"
    )
  }

  cor_method <- NULL
  rotation <- NULL
  orig_R <- NULL

  if (inherits(f, "EFA")) {

    # The factor intercorrelations come from the fitted solution, so an explicitly
    # supplied `Phi` cannot be honoured here; say so instead of discarding it silently.
    if (!is.null(Phi)) {
      cli::cli_warn(
        c("{.arg Phi} is ignored when {.arg f} is an {.cls efa} object.",
          "i" = "The factor intercorrelations of the fitted solution are used instead.",
          "i" = "Pass the loading matrix (for example {.code f$rot_loadings}) in {.arg f} to score against a different {.arg Phi}."),
        class = "efa_scores_phi_ignored"
      )
    }

    Phi <- f$Phi
    cor_method <- f$settings$cor_method
    rotation <- f$settings$rotation
    # The loadings and communalities that the weights must stay consistent with:
    # the rotated pattern for a rotated solution, the unrotated loadings otherwise.
    # h2 is the model communality diag(Lambda Phi Lambda') (rotation-invariant), the
    # correct uniqueness source for the Bartlett/Anderson weighting.
    if (!identical(rotation, "none")) {
      Lambda <- unclass(f$rot_loadings)
    } else {
      Lambda <- unclass(f$unrot_loadings)
    }
    h2 <- f$h2
    orig_R <- unclass(f$orig_R)

  } else {

    Lambda <- unclass(f)

    if (is.null(Phi)) {
      cli::cli_inform(
        c("i" = "{.arg Phi} was left {.code NULL} and loadings were entered directly in {.arg f}; assuming uncorrelated factors."),
        class = "efa_scores_phi_null"
      )
    }

  }

  if (!is.matrix(Lambda) || !is.numeric(Lambda) ||
      nrow(Lambda) < 1L || ncol(Lambda) < 1L || any(!is.finite(Lambda))) {
    cli::cli_abort(
      "{.arg f} must provide a non-empty numeric loading matrix with only finite values.",
      class = "efa_scores_bad_loadings"
    )
  }

  # Label the factor dimension (F1..Fm) when the loadings carry no usable factor
  # names, so the weights, diagnostics, scores, and determinacy table are
  # self-describing (the same factor naming the EFA print methods use). Duplicate
  # names are treated as unusable too: they would otherwise abort the determinacy
  # table below with an opaque base-R "duplicate row.names" error.
  factor_names <- colnames(Lambda)
  if (is.null(factor_names) || anyNA(factor_names) ||
      !all(nzchar(factor_names)) || anyDuplicated(factor_names) > 0L) {
    colnames(Lambda) <- paste0("F", seq_len(ncol(Lambda)))
  }

  if (!is_cmat) {
    raw_columns_ok <- if (is.data.frame(x)) {
      all(vapply(x, function(z) is.numeric(z) || is.logical(z), logical(1)))
    } else {
      is.numeric(x) || is.logical(x)
    }
    if (!raw_columns_ok) {
      cli::cli_abort(
        "Raw {.arg x} must contain only numeric or logical variables.",
        class = "efa_scores_bad_x"
      )
    }
    # Logical variables represent binary 0/1 observations and are valid inputs,
    # but matrix algebra and stats::cor() should see numeric storage explicitly.
    if (is.logical(x) || (is.data.frame(x) && any(vapply(x, is.logical, logical(1))))) {
      x <- as.matrix(x)
      storage.mode(x) <- "double"
    }
  }

  # Uncorrelated factors when no Phi is available (EFA with an orthogonal/unrotated
  # solution, or a directly supplied loading matrix without Phi).
  if (is.null(Phi)) {
    Phi <- diag(ncol(Lambda))
  }

  Phi <- .align_correlation_axis(
    Phi, n = ncol(Lambda), target_names = colnames(Lambda), arg = "Phi"
  )

  # A directly supplied loading matrix carries no stored communalities, so derive
  # them from the (now resolved) model: h2 = diag(Lambda Phi Lambda').
  if (!inherits(f, "EFA")) {
    h2 <- diag(Lambda %*% Phi %*% t(Lambda))
  }

  # Anderson-Rubin orthogonalisation is undefined with a single factor (there is
  # nothing to make uncorrelated), so guard it before any computation.
  if (method == "Anderson" && ncol(Lambda) == 1L) {
    cli::cli_abort(
      c("Anderson-Rubin scores are not defined for a single factor.",
        "i" = "There is nothing to orthogonalise with one factor; use {.code method = \"regression\"} or {.code method = \"Bartlett\"}."),
      class = "efa_scores_anderson_single"
    )
  }

  # Anderson-Rubin scores are defined for orthogonal factors: the weights force the
  # score intercorrelations to the identity, so for a model whose factors are
  # correlated they report uncorrelated scores for factors that are not. Warn when Phi
  # is not diagonal up to rounding -- an orthogonal or unrotated solution carries no
  # Phi and resolves to an exact identity above, so this fires only on a genuinely
  # oblique solution (or a hand-supplied Phi).
  if (method == "Anderson") {
    off <- abs(Phi[upper.tri(Phi)])
    if (any(off > sqrt(.Machine$double.eps), na.rm = TRUE)) {
      cli::cli_warn(
        c("Anderson-Rubin scores are defined for orthogonal factors, but the factors of {.arg f} are correlated (largest {.code |r|} = {round(max(off, na.rm = TRUE), 3)}).",
          "i" = "The scores are orthogonalised regardless, so they are uncorrelated even though the factors of the model are not.",
          "i" = "Use {.code method = \"tenBerge\"} to preserve the factor intercorrelations."),
        class = "efa_scores_anderson_oblique"
      )
    }
  }

  # Align raw `x` to the model variables. The weight matrix `W` (and the structure
  # products it is built from) combines variables by column position, but `W`, like
  # `Lambda`, carries the fitted model's variable names in the model's order. A raw
  # `x` whose columns are in a different order would otherwise be scored against the
  # wrong variables -- and in the loadings-matrix path also corrupt the `cor(x)` the
  # weights are derived from. When both sides carry usable, unique names, reorder `x`
  # to the model order by name (dropping any extra, non-model columns); a model
  # variable absent from `x` cannot be scored, so abort. Fall back to positional
  # matching when either side is unnamed or carries duplicate names.
  if (!is_cmat) {
    model_vars <- rownames(Lambda)
    x_names <- colnames(x)
    usable_model_vars <- !is.null(model_vars) && !anyNA(model_vars) &&
      all(nzchar(model_vars)) && anyDuplicated(model_vars) == 0L
    if (usable_model_vars && !is.null(x_names)) {
      if (anyNA(x_names) || !all(nzchar(x_names)) || anyDuplicated(x_names) > 0L) {
        cli::cli_abort(
          "The column names of {.arg x} must be complete, non-empty, and unique for name-based alignment.",
          class = "efa_scores_x_names"
        )
      }
      missing_vars <- setdiff(model_vars, x_names)
      if (length(missing_vars) > 0L) {
        cli::cli_abort(
          c("{.arg x} is missing {cli::qty(missing_vars)} model variable{?s}: {.val {missing_vars}}.",
            "i" = "Factor scores need every variable in the loading matrix of {.arg f}."),
          class = "efa_scores_missing_vars"
        )
      }
      x <- x[, model_vars, drop = FALSE]
    }
  } else {
    # A correlation matrix `x` describes the same model axis, and takes the checks the
    # scoring matrix takes, whether or not it ends up being that matrix: a fitted `f`
    # carries the correlations the model was estimated from, so `x` is not read for the
    # weights there, and left unchecked a matrix of the wrong size or of entirely
    # different variables returned the fitted solution's weights as though it had produced
    # them. Aligned as well as checked, so the branch below that does read it gets the
    # model's variable order, as the raw branch above does. Coerced first because `x` also
    # takes a data frame, which a correlation matrix is as readily given as.
    x <- .align_correlation_axis(as.matrix(x), n = nrow(Lambda),
                                 target_names = rownames(Lambda), arg = "x")
  }

  # A model variable that carries no usable spread cannot be scored: a constant column
  # divides by a zero standard deviation and becomes `NaN`, and because a score sums
  # over all variables, one such column destroys the scores of every case. The rule is
  # the same for every method -- a model variable that does not vary means the fitted
  # correlation structure does not apply to these data. Checked after the alignment
  # above, so only the variables actually used for scoring are examined, and before the
  # scoring correlation is resolved, where the same column would otherwise surface as an
  # opaque non-finite correlation.
  #
  # sd() is exactly 0 for a constant column, NaN for one carrying an infinity, and NA
  # for one with fewer than two observed values, so the causes are separated before they
  # are reported -- an infinite value is not an absent one -- the same split the
  # uncomputable-correlation diagnosis makes. A small but positive standard deviation is
  # ill-conditioned rather than undefined and is left to standardise as usual. Each test
  # yields FALSE rather than NA on a column it cannot judge, so `any()` cannot abort with
  # an unclassed base error in place of the classed condition.
  if (!is_cmat) {
    sds <- if (is.data.frame(x)) {
      vapply(x, stats::sd, numeric(1), na.rm = TRUE)
    } else {
      apply(x, 2L, stats::sd, na.rm = TRUE)
    }
    nonfinite <- is.nan(sds)
    unobserved <- is.na(sds) & !nonfinite
    constant <- !is.na(sds) & sds == 0
    unusable <- nonfinite | unobserved | constant
    if (any(unusable)) {
      # Name the offenders, falling back to their column positions as character labels
      # so that cli::qty() reads the number of them rather than a position value. Long
      # lists go through the shared cap, as in every other condition that enumerates
      # affected variables.
      nms <- colnames(x)
      if (is.null(nms)) nms <- as.character(seq_along(sds))
      n_bad <- sum(unusable)
      cap <- .cap_label_list(nms[unusable])
      why <- if (all(constant[unusable])) {
        "{cli::qty(n_bad)}{?It has/They have} no variance, so {?it carries/they carry} no information about any factor."
      } else if (all(nonfinite[unusable])) {
        "{cli::qty(n_bad)}{?It contains/They contain} infinite values, which combine with no factor."
      } else if (all(unobserved[unusable])) {
        "{cli::qty(n_bad)}{?It has/They have} fewer than two observed values, so {?its/their} spread is undefined."
      } else {
        "{cli::qty(n_bad)}{?It has/They have} no variance, infinite values, or fewer than two observed values."
      }
      cli::cli_abort(
        c("{cli::qty(n_bad)}Model variable{?s} {.val {cap$shown}}{cap$rest}{cli::qty(n_bad)} cannot be scored from {.arg x}.",
          "x" = why,
          "i" = "Correct the affected values, or drop {cli::qty(n_bad)}{?it/them} from the factor model."),
        class = "efa_scores_constant_var"
      )
    }
  }

  # Scoring correlation matrix R. A user-supplied `rho` wins; otherwise use the
  # matrix the EFA was fit on (f$orig_R), so the weights are consistent with the
  # loadings even for a non-Pearson fit; otherwise fall back to the correlation
  # of `x` (the matrix itself when `x` is already a correlation matrix).
  if (!is.null(rho)) {
    R <- unclass(rho)
  } else if (!is.null(orig_R)) {
    R <- orig_R
  } else if (is_cmat) {
    R <- unclass(as.matrix(x))
  } else {
    R <- stats::cor(x, use = "pairwise")
  }

  R_source <- if (!is.null(rho)) "rho" else if (!is.null(orig_R)) "f$orig_R" else "x"
  R <- .align_correlation_axis(
    R, n = nrow(Lambda), target_names = rownames(Lambda), arg = R_source
  )

  # A non-Pearson EFA scored on raw data with rho = NULL: the weights and
  # diagnostics correctly use the fitted correlation (f$orig_R), but the scores are
  # formed from the Pearson-standardised raw data, so realised score properties can
  # differ from the reported diagnostics. Warn (and point to `rho`) unless the fit
  # is Pearson, a matrix was supplied via `rho`, or `x` is a correlation matrix.
  if (!is.null(cor_method) && !identical(cor_method, "pearson") &&
      is.null(rho) && !is_cmat) {
    cli::cli_warn(
      c("{.arg f} was fit with {.code cor_method = {.val {cor_method}}}; the weights and diagnostics use that fitted correlation ({.code f$orig_R}).",
        "i" = "The factor scores are formed from the Pearson-standardised raw data, so their realised properties may differ from the reported diagnostics.",
        "i" = "Pass the {.val {cor_method}} correlation matrix (or another) via {.arg rho} to control the scoring correlation."),
      class = "efa_scores_cor_method"
    )
  }

  # psych's components path multiplies `x %*% W` without coercing; coerce a raw
  # data.frame up front (a no-op for the scale()-based methods).
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  W <- .factor_score_weights(Lambda, Phi, R, h2, method = method)
  S <- Lambda %*% Phi
  diag_out <- .factor_score_diagnostics(W, R, S)

  # Factor scores only from raw data. Component scores are the raw (uncentered)
  # data times the loadings; the other methods use the standardized data. A case with
  # a missing value on any model variable cannot be scored (no imputation is
  # performed), so its score row is NA -- count those cases and say so, rather than
  # returning NA rows the caller only discovers downstream.
  scores <- NULL
  n_scored <- NA_integer_
  if (!is_cmat) {
    scores <- if (method == "components") x %*% W else scale(x) %*% W
    n_na <- sum(!stats::complete.cases(scores))
    n_scored <- nrow(scores) - n_na
    if (n_na > 0L) {
      cli::cli_warn(
        c("{n_na} case{?s} {?has/have} a missing value on at least one model variable, so {?its/their} factor score{?s} {?is/are} {.code NA}.",
          "i" = "No imputation is performed: handle the missing values before scoring, or drop the affected cases."),
        class = "efa_scores_incomplete_cases"
      )
    }
  }

  determinacy <- data.frame(
    rho = diag_out$determinacy,
    rho2 = diag_out$determinacy^2,
    guttman = diag_out$guttman,
    row.names = colnames(Lambda)
  )

  output <- list(
    weights = W,
    scores = scores,
    r.scores = diag_out$r.scores,
    score_cor = diag_out$score_cor,
    determinacy = determinacy,
    settings = list(
      method = method,
      cor_method = cor_method,
      n_factors = ncol(Lambda),
      n_obs = if (is.null(scores)) NA_integer_ else nrow(scores),
      n_scored = n_scored
    )
  )

  class(output) <- "efa_scores"

  output

}

# Validate a factor- or variable-correlation matrix and, when both dimensions and
# the corresponding model axis carry unique names, align it to the model by name.
# A wholly unnamed matrix keeps the documented positional semantics; partially
# named or duplicate-named input cannot be aligned safely and is rejected.
.align_correlation_axis <- function(x, n, target_names, arg) {
  if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(n, n))) {
    cli::cli_abort(
      "{.arg {arg}} must be a numeric {n} by {n} correlation matrix.",
      class = "efa_scores_matrix_dim"
    )
  }
  if (any(!is.finite(x))) {
    cli::cli_abort(
      "{.arg {arg}} must contain only finite values.",
      class = "efa_scores_matrix_nonfinite"
    )
  }

  usable_target <- !is.null(target_names) && !anyNA(target_names) &&
    all(nzchar(target_names)) && anyDuplicated(target_names) == 0L
  if (usable_target) {
    rn <- rownames(x)
    cn <- colnames(x)
    wholly_unnamed <- is.null(rn) && is.null(cn)

    if (!wholly_unnamed) {
      usable_names <- !is.null(rn) && !is.null(cn) && identical(rn, cn) &&
        !anyNA(rn) && !anyNA(cn) && all(nzchar(rn)) && all(nzchar(cn)) &&
        anyDuplicated(rn) == 0L && anyDuplicated(cn) == 0L
      same_names <- usable_names && setequal(rn, target_names) &&
        setequal(cn, target_names)
      if (!same_names) {
        cli::cli_abort(
          c("The row and column names of {.arg {arg}} do not identify the model axis unambiguously.",
            "i" = "Use the same unique names as the corresponding loading matrix dimension, or remove both sets of names to match by position."),
          class = "efa_scores_matrix_names"
        )
      }
      x <- x[target_names, target_names, drop = FALSE]
    }
  }

  # Validate after name alignment. Requiring the two labelled axes to use the same
  # order above also prevents a numerically symmetric matrix with one mislabelled
  # axis from becoming invalid only after subsetting.
  tol <- sqrt(.Machine$double.eps)
  if (max(abs(x - t(x))) > tol * max(1, max(abs(x)))) {
    cli::cli_abort(
      "{.arg {arg}} must be symmetric.",
      class = "efa_scores_matrix_symmetric"
    )
  }
  # Canonicalize accepted round-off asymmetry so symmetric eigensolvers and general
  # matrix operations consume the same matrix on every BLAS/LAPACK implementation.
  x <- x / 2 + t(x) / 2
  range_tol <- 100 * .Machine$double.eps
  if (max(abs(diag(x) - 1)) > range_tol || any(abs(x) > 1 + range_tol)) {
    cli::cli_abort(
      "{.arg {arg}} must have a unit diagonal and entries between -1 and 1.",
      class = "efa_scores_matrix_correlation"
    )
  }
  diag(x) <- 1
  x[x > 1] <- 1
  x[x < -1] <- -1

  x
}


#' Print and format an efa_scores object
#'
#' `print()` shows a concise overview of an [efa_scores()] result: a header naming
#' the method and whether factor scores were computed, and the per-factor
#' determinacy table (determinacy, squared determinacy, and Guttman index).
#' `summary()` returns a `summary.efa_scores` object whose print method adds the
#' full factor-weight matrix, the score validity/univocality matrix, and the score
#' intercorrelations. `format()` assembles the same report and returns it as a
#' character vector; `print()` is `cat(format(x), sep = "\n")`. The lines follow
#' the active console theme, so they are plain when colours are disabled (for
#' example when captured into a file or stripped with [cli::ansi_strip()]).
#'
#' @param x,object An object of class `efa_scores`; for the `summary.efa_scores`
#'   methods, the object returned by `summary()`.
#' @param digits numeric. Number of decimal places for the printed tables.
#'   Default is 3.
#' @param ... Not used; for consistency with the generics.
#'
#' @returns `print()` and the print method for `summary.efa_scores` objects return
#'   their argument invisibly. `format()` returns a character vector with the
#'   report lines. `summary()` returns an object of class `summary.efa_scores`.
#'
#' @family factor scoring
#'
#' @export
#'
#' @method print efa_scores
#'
#' @examples
#' efa <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                estimator = "PAF", rotation = "oblimin")
#' fs <- efa_scores(test_models$baseline$cormat, f = efa)
#' fs
#' summary(fs)
#'
#' # format() returns the same lines as a character vector:
#' writeLines(format(fs))
#'
print.efa_scores <- function(x, digits = 3, ...) {
  cat(format(x, digits = digits, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.efa_scores
#' @export
#' @method format efa_scores
format.efa_scores <- function(x, digits = 3, ...) {
  .render_efa_scores(x, view = "brief", digits = digits)
}

#' @rdname print.efa_scores
#' @export
#' @method summary efa_scores
summary.efa_scores <- function(object, digits = 3, ...) {
  structure(list(scores = object, opts = list(digits = digits)),
            class = "summary.efa_scores")
}

#' @rdname print.efa_scores
#' @export
#' @method print summary.efa_scores
print.summary.efa_scores <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.efa_scores
#' @export
#' @method format summary.efa_scores
format.summary.efa_scores <- function(x, digits = x$opts$digits, ...) {
  .render_efa_scores(x$scores, view = "full", digits = digits)
}

# Render the factor-score report at the requested depth: "brief" (print) shows the
# header and the per-factor determinacy table; "full" (summary) additionally shows
# the factor-weight matrix, the score validity/univocality matrix, and the score
# intercorrelations. Assembled inside a single cli container so prose wraps to the
# console width; the aligned tables are emitted verbatim (no reflow).
.render_efa_scores <- function(x, view = c("brief", "full"), digits = 3) {
  view <- .match_arg_ci(view)
  full <- identical(view, "full")

  method <- x$settings$method
  title <- paste0("Factor scores (", method, ")")

  # Every table here holds bounded coefficients or correlations, so all four go through
  # the shared renderer and print in the package's number convention (leading zero
  # dropped), matching the coefficient tables of efa_reliability() and the loading
  # tables. The renderer also splits a table wider than the console into stacked column
  # blocks, which the weight matrix reaches on a wide solution.
  cli::cli_format_method({
    cli::cli_text("")
    cli::cli_rule(left = "{.strong {title}}")
    cli::cli_text("")

    if (is.null(x$scores)) {
      cli::cli_text("Weights and diagnostics only (correlation-matrix input; no scores).")
    } else {
      n_scored <- x$settings$n_scored
      n_obs <- x$settings$n_obs
      n_fac <- x$settings$n_factors
      # The pointer names a field of the printed object, which is what `print()` was
      # handed; `summary()` wraps that object in a field of the same name, so the
      # pointer would resolve to the wrapper there and is left off the fuller report.
      # Written out twice rather than interpolated, because a `{.code ...}` carried in
      # through a substitution is emitted as literal text instead of being styled.
      if (full) {
        cli::cli_text("Scored {n_scored} of {n_obs} observation{?s} on {n_fac} factor{?s}.")
      } else {
        cli::cli_text("Scored {n_scored} of {n_obs} observation{?s} on {n_fac} factor{?s} (see {.code $scores}).")
      }
    }

    cli::cli_text("")
    cli::cli_rule(left = "{.strong Score determinacy}")
    cli::cli_text("")
    .efa_emit_lines(.efa_corr_lines(as.matrix(x$determinacy), digits = digits))

    if (full) {
      cli::cli_text("")
      cli::cli_rule(left = "{.strong Factor weights}")
      cli::cli_text("")
      .efa_emit_lines(.efa_corr_lines(x$weights, digits = digits))

      cli::cli_text("")
      cli::cli_rule(left = "{.strong Score validity and univocality}")
      cli::cli_text("")
      cli::cli_text("Diagonal: validity (score-factor correlation). Off-diagonal: univocality.")
      cli::cli_text("")
      .efa_emit_lines(.efa_corr_lines(x$score_cor, digits = digits))

      cli::cli_text("")
      cli::cli_rule(left = "{.strong Score intercorrelations}")
      cli::cli_text("")
      .efa_emit_lines(.efa_corr_lines(x$r.scores, digits = digits))
    }
  })
}


# Native factor-score weight engine ------------------------------------------
#
# Computes the p x m factor-score weight matrix W (standardized scores are
# scale(X) %*% W) natively, without psych::factor.scores(). With the structure
# matrix S = Lambda %*% Phi and the model uniquenesses Psi = diag(1 - h2), the
# supported methods are:
#   Thurstone / regression : W = R^-1 S                            (Thurstone, 1935)
#   Bartlett               : W = Psi^-1 Lambda (Lambda' Psi^-1 Lambda)^-1
#                                                                   (Bartlett, 1937)
#   Anderson (orthogonal)  : W = Psi^-1 Lambda
#                                (Lambda' Psi^-1 R Psi^-1 Lambda)^-1/2
#                                                          (Anderson & Rubin, 1956)
#   tenBerge               : L = Lambda Phi^1/2;
#                            W = R^-1/2 [R^-1/2 L (L' R^-1 L)^-1/2] Phi^1/2
#                     (ten Berge, Krijnen, Wansbeek & Shapiro, 1999; Grice, 2001)
#   Harman                 : M = Lambda (Lambda Phi)', diag(M) = 1; W = M^-1 Lambda
#                                                                    (Harman, 1976)
#   components             : W = Lambda
# For orthogonal regression the diagonal of Lambda' R^-1 Lambda gives the
# factor-score validity coefficients rho^2, i.e. each factor's squared multiple
# correlation with the observed variables (Grice, 2001).

# Symmetric matrix square root A^1/2 = V diag(sqrt(lambda)) V' via the eigen-
# decomposition; negative eigenvalues are clamped to .Machine$double.eps (as in
# psych's matSqrt). `nrow=` keeps diag() a matrix for a single eigenvalue.
.mat_sqrt <- function(x) {
  e <- eigen(x, symmetric = TRUE)
  ev <- e$values
  ev[ev < 0] <- .Machine$double.eps
  e$vectors %*% diag(sqrt(ev), nrow = length(ev)) %*% t(e$vectors)
}

# Symmetric inverse matrix square root A^-1/2; eigenvalues below
# .Machine$double.eps are clamped to 100 * .Machine$double.eps (as in psych's
# invMatSqrt). `nrow=` keeps diag() a matrix for a single eigenvalue.
.inv_mat_sqrt <- function(x) {
  e <- eigen(x, symmetric = TRUE)
  ev <- e$values
  ev[ev < .Machine$double.eps] <- 100 * .Machine$double.eps
  e$vectors %*% diag(1 / sqrt(ev), nrow = length(ev)) %*% t(e$vectors)
}

# Moore-Penrose pseudo-inverse via the SVD, dropping singular values at or below
# tol * (largest singular value); fallback for a singular correlation or Gram
# matrix (mirrors psych's Pinv).
.pinv <- function(X, tol = sqrt(.Machine$double.eps)) {
  s <- svd(X)
  pos <- s$d > max(tol * s$d[1], 0)
  if (all(pos)) {
    s$v %*% (t(s$u) / s$d)
  } else {
    s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
  }
}

# Factor-score weight matrix W (p x m) for the requested method; see the block
# comment above for the formulas and references. `h2` must be the model-
# reproduced communalities diag(Lambda Phi Lambda') (not rowSums(Lambda^2),
# which differ for oblique solutions) so the Bartlett/Anderson uniquenesses are
# correct. Anderson-Rubin scores are defined for orthogonal factors only.
.factor_score_weights <- function(Lambda, Phi, R, h2,
                                  method = c("Thurstone", "regression", "tenBerge",
                                             "Anderson", "Bartlett", "Harman",
                                             "components")) {

  method <- .match_arg_ci(method)
  if (method == "regression") method <- "Thurstone"

  S <- Lambda %*% Phi

  # solve(A, B) with a Moore-Penrose fallback when A is singular.
  solve_or_pinv <- function(A, B) {
    tryCatch(solve(A, B), error = function(e) .pinv(A) %*% B)
  }

  # The model uniquenesses u2 = 1 - h2 that the Bartlett and Anderson methods weight
  # by. Both need Psi^-1 Lambda = diag(1 / u2) %*% Lambda, which is the row scaling
  # `Lambda / u2`, so the dense p x p diagonal matrix is never formed. A non-positive
  # uniqueness makes the weighting undefined; abort instead of dividing by zero or
  # returning finite-looking output from an invalid negative weight. A very small
  # positive uniqueness is defined but ill-conditioned, so retain the existing warning
  # for that distinct case.
  uniquenesses <- function() {
    u2 <- 1 - h2
    nonpositive <- u2 <= 0
    if (any(nonpositive)) {
      bad <- rownames(Lambda)[nonpositive]
      if (is.null(bad)) bad <- which(nonpositive)
      cli::cli_abort(
        c("Non-positive uniqueness values make the {method} factor-score weights undefined (affected variables: {.val {bad}}).",
          "i" = "Use a scoring method that does not invert uniquenesses, or revise the factor model."),
        class = "efa_scores_nonpositive_uniqueness"
      )
    }
    near_zero <- u2 < .Machine$double.eps
    if (any(near_zero)) {
      bad <- rownames(Lambda)[near_zero]
      if (is.null(bad)) bad <- which(near_zero)
      cli::cli_warn(
        c("Near-zero uniqueness values were detected for {.val {bad}}.",
          "i" = "The affected {method} factor-score weights are unstable."),
        class = "efa_scores_heywood"
      )
    }
    u2
  }

  W <- switch(
    method,
    Thurstone = solve_or_pinv(R, S),
    Bartlett = {
      Lu <- Lambda / uniquenesses()                 # Psi^-1 Lambda
      Lu %*% solve_or_pinv(crossprod(Lu, Lambda), diag(ncol(Lambda)))
    },
    Anderson = {
      Lu <- Lambda / uniquenesses()                 # Psi^-1 Lambda
      Lu %*% .inv_mat_sqrt(crossprod(Lu, R %*% Lu))
    },
    tenBerge = {
      phi_sqrt <- .mat_sqrt(Phi)
      L <- Lambda %*% phi_sqrt
      r_inv_sqrt <- .inv_mat_sqrt(R)
      C <- r_inv_sqrt %*% L %*% .inv_mat_sqrt(crossprod(L, solve_or_pinv(R, L)))
      r_inv_sqrt %*% C %*% phi_sqrt
    },
    Harman = {
      M <- Lambda %*% t(S)
      diag(M) <- 1
      solve_or_pinv(M, Lambda)
    },
    components = Lambda
  )

  dimnames(W) <- dimnames(Lambda)
  W

}


# Score-quality diagnostics for a factor-score weight matrix W (from
# .factor_score_weights). With the standardized scores Fhat = Z %*% W (Z the standardized
# data, cov(Z) = R) and the structure matrix S = Lambda %*% Phi, C = W' R W is the
# covariance of the score estimates. Reports:
#   r.scores    : cov2cor(C), the score intercorrelations; for a well-behaved solution
#                 these reproduce the factor correlations Phi (exactly for ten Berge, and
#                 the identity for Anderson-Rubin).
#   score_cor   : the full score-factor correlation matrix diag(sd.in) W' S,
#                 sd.in = 1/sqrt(diag(C)); its diagonal is each score's determinacy /
#                 validity rho, its off-diagonals the univocality (the correlation of a
#                 score with the *other* factors).
#   determinacy : the diagonal of score_cor, the correlation between each factor and its
#                 estimated score for the *returned* (method-specific) weights, capped at
#                 1 with non-positive values set to NA (psych's convention).
#   guttman     : Guttman's (1955) indeterminacy index 2 rho^2 - 1, the minimum
#                 correlation between two equally valid sets of factor scores.
# The factors carry unit variance (diag(Phi) == 1), so Phi does not enter the score-
# factor correlations and is not needed here.
# References: Guttman (1955); ten Berge, Krijnen, Wansbeek & Shapiro (1999); Grice (2001).
.factor_score_diagnostics <- function(W, R, S) {

  C <- crossprod(W, R %*% W)                 # t(W) %*% R %*% W  (m x m)
  r.scores <- stats::cov2cor(C)

  sd.in <- 1 / sqrt(diag(C))
  # corr(Fhat_i, F_j): the factors are standardised (unit variance), so only the score
  # standard deviations enter the score-factor correlations.
  score_cor <- diag(sd.in, nrow = length(sd.in)) %*% crossprod(W, S)
  dimnames(score_cor) <- list(colnames(W), colnames(W))
  dimnames(r.scores)  <- dimnames(score_cor)

  # Determinacy / validity: the diagonal. Following psych's convention, cap at 1 and set
  # any non-positive value to NA (an ill-conditioned score whose "correlation" leaves
  # [0, 1]); a non-finite entry (a NaN reachable only from a non-PD R) is treated as NA
  # too, so the reported determinacy and Guttman index never leak a NaN.
  determinacy <- diag(score_cor)
  determinacy[determinacy > 1] <- 1
  determinacy[!is.finite(determinacy) | determinacy <= 0] <- NA
  names(determinacy) <- colnames(W)

  guttman <- 2 * determinacy^2 - 1

  list(r.scores = r.scores, score_cor = score_cor,
       determinacy = determinacy, guttman = guttman)

}
