#' Estimate factor scores and score-quality diagnostics for an EFA model
#'
#' Computes factor-score weights (and, from raw data, the factor scores
#' themselves) natively for an [EFA()] solution or a directly supplied loading
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
#'   \item{`"components"`}{component scores, `W = Lambda`.}
#' }
#' The determinacy of a score is its correlation with the factor it estimates,
#' `rho_i = cor(F_i, Fhat_i)`, computed from the returned weights. For regression
#' scores this equals the multiple correlation between the factor and the observed
#' variables, `sqrt(diag(S' R^-1 S))` (Guttman, 1955; Grice, 2001). The
#' off-diagonal score-factor correlations give the univocality (the correlation of
#' a score with the *other* factors), and `2 rho^2 - 1` is Guttman's (1955)
#' indeterminacy index, the minimum correlation between two equally valid sets of
#' scores. For a `method` other than `"regression"` both quantities are specific to
#' those scores: the determinacy is the chosen method's own score-factor correlation
#' (never larger than the regression value), and the reported `guttman` follows from
#' it rather than being the model's maximal indeterminacy index.
#'
#' When `f` is an [EFA()] object the scoring correlation matrix defaults to the
#' matrix the model was fit on (`f$orig_R`), so the weights stay consistent with
#' the loadings even for a non-Pearson correlation (e.g. polychoric). A matrix
#' supplied via `rho` overrides this.
#'
#' @param x data.frame or matrix. Raw data (needed to obtain factor scores) or a
#'   correlation matrix (yields weights and diagnostics only). When raw data carry
#'   column names, they are matched to the model variables by name (any extra
#'   columns are ignored, and a model variable missing from `x` is an error);
#'   unnamed data are matched by position.
#' @param f object of class [EFA()], a `LOADINGS` object, or a matrix of factor
#'   loadings.
#' @param Phi matrix. Factor intercorrelations. Only used when a loading matrix is
#'   supplied directly in `f`; taken from the `EFA` object otherwise. Default is
#'   `NULL`, in which case the factors are assumed uncorrelated.
#' @param rho matrix. Correlation matrix used to derive the scoring weights.
#'   Defaults to `NULL`, in which case `f$orig_R` is used for an `EFA` object and
#'   `cor(x, use = "pairwise")` otherwise. Pass a matrix here to score against a
#'   correlation other than the one implied by `f`/`x`.
#' @param method character. The factor-score method: one of `"regression"`
#'   (default), `"Bartlett"`, `"Anderson"`, `"tenBerge"`, `"Harman"`, or
#'   `"components"`.
#'
#' @return An object of class `efa_scores`, a list containing:
#'
#' \item{weights}{The `p` by `m` factor-score weight matrix.}
#' \item{scores}{The factor scores (`n` by `m`), or `NULL` when a correlation
#'   matrix was supplied.}
#' \item{r.scores}{The `m` by `m` correlations of the factor-score estimates.}
#' \item{score_cor}{The `m` by `m` score-factor correlation matrix; its diagonal
#'   is the determinacy (validity) of each score and its off-diagonals the
#'   univocality.}
#' \item{determinacy}{A data frame with, per factor, the determinacy `rho`, the
#'   squared determinacy `rho2`, and Guttman's indeterminacy index `guttman`.}
#' \item{settings}{A list of the settings used.}
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
#' @family factor scoring
#'
#' @export
#'
#' @examples
#' # Weights and score diagnostics from an EFA on a correlation matrix
#' efa <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'            method = "PAF", rotation = "oblimin")
#' fs <- efa_scores(test_models$baseline$cormat, f = efa)
#' fs
#' summary(fs)
#'
#' # Factor scores from raw data (Bartlett method)
#' \donttest{
#' efa_raw <- EFA(GRiPS_raw, n_factors = 1, method = "PAF")
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

  method <- match.arg(method)
  checkmate::assert_matrix(Phi, null.ok = TRUE)

  if (!inherits(f, c("EFA", "matrix", "LOADINGS"))) {
    cli::cli_abort(
      "{.arg f} must be an {.cls EFA} object, a matrix, or a {.cls LOADINGS} object.",
      class = "efa_scores_bad_f"
    )
  }

  cor_method <- NULL
  rotation <- NULL
  orig_R <- NULL

  if (inherits(f, "EFA")) {

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

  # Label the factor dimension (F1..Fm) when the loadings carry no usable factor
  # names, so the weights, diagnostics, scores, and determinacy table are
  # self-describing (the same factor naming the EFA print methods use). Duplicate
  # names are treated as unusable too: they would otherwise abort the determinacy
  # table below with an opaque base-R "duplicate row.names" error.
  if (is.null(colnames(Lambda)) || anyDuplicated(colnames(Lambda)) > 0L) {
    colnames(Lambda) <- paste0("F", seq_len(ncol(Lambda)))
  }

  # Uncorrelated factors when no Phi is available (EFA with an orthogonal/unrotated
  # solution, or a directly supplied loading matrix without Phi).
  if (is.null(Phi)) {
    Phi <- diag(ncol(Lambda))
  }

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
    if (!is.null(model_vars) && !is.null(x_names) &&
        anyDuplicated(model_vars) == 0L && anyDuplicated(x_names) == 0L) {
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
  # data times the loadings; the other methods use the standardized data.
  scores <- NULL
  if (!is_cmat) {
    scores <- if (method == "components") x %*% W else scale(x) %*% W
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
      n_obs = if (is.null(scores)) NA_integer_ else nrow(scores)
    )
  )

  class(output) <- "efa_scores"

  output

}


#' Print and format an efa_scores object
#'
#' `print()` shows a concise overview of an [efa_scores()] result: a header naming
#' the method and whether factor scores were computed, and the per-factor
#' determinacy table (determinacy, squared determinacy, and Guttman index).
#' [summary()] returns a `summary.efa_scores` object whose print method adds the
#' full factor-weight matrix, the score validity/univocality matrix, and the score
#' intercorrelations. `format()` assembles the same report and returns it as a
#' character vector; `print()` is `cat(format(x), sep = "\n")`. The lines follow
#' the active console theme, so they are plain when colours are disabled (for
#' example when captured into a file or stripped with [cli::ansi_strip()]).
#'
#' @param x,object An object of class `efa_scores`; for the `summary.efa_scores`
#'   methods, the object returned by [summary()].
#' @param digits numeric. Number of decimal places for the printed tables.
#'   Default is 3.
#' @param ... Not used; for consistency with the generics.
#'
#' @returns `print()` and the print method for `summary.efa_scores` objects return
#'   their argument invisibly. `format()` returns a character vector with the
#'   report lines (styled to the active console theme; plain when colours are
#'   disabled). `summary()` returns an object of class `summary.efa_scores`.
#'
#' @family factor scoring
#'
#' @export
#'
#' @method print efa_scores
#'
#' @examples
#' efa <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'            method = "PAF", rotation = "oblimin")
#' fs <- efa_scores(test_models$baseline$cormat, f = efa)
#' fs
#' summary(fs)
#'
#' # format() returns the same lines as plain text:
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
  view <- match.arg(view)
  full <- identical(view, "full")

  method <- x$settings$method
  title <- paste0("Factor scores (", method, ")")

  # Pre-format the tables outside the cli container (base print, not reflowed).
  det_lines <- utils::capture.output(print(round(x$determinacy, digits)))

  cli::cli_format_method({
    cli::cli_text("")
    cli::cli_rule(left = "{.strong {title}}")
    cli::cli_text("")

    if (is.null(x$scores)) {
      cli::cli_text("Weights and diagnostics only (correlation-matrix input; no scores).")
    } else {
      n_obs <- x$settings$n_obs
      n_fac <- x$settings$n_factors
      cli::cli_text("Scored {n_obs} observation{?s} on {n_fac} factor{?s}.")
    }

    cli::cli_text("")
    cli::cli_rule(left = "{.strong Score determinacy}")
    cli::cli_text("")
    cli::cli_verbatim(det_lines)

    if (full) {
      cli::cli_text("")
      cli::cli_rule(left = "{.strong Factor weights}")
      cli::cli_text("")
      cli::cli_verbatim(utils::capture.output(print(round(x$weights, digits))))

      cli::cli_text("")
      cli::cli_rule(left = "{.strong Score validity and univocality}")
      cli::cli_text("")
      cli::cli_text("Diagonal: validity (score-factor correlation). Off-diagonal: univocality.")
      cli::cli_text("")
      cli::cli_verbatim(utils::capture.output(print(round(x$score_cor, digits))))

      cli::cli_text("")
      cli::cli_rule(left = "{.strong Score intercorrelations}")
      cli::cli_text("")
      cli::cli_verbatim(utils::capture.output(print(round(x$r.scores, digits))))
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

  method <- match.arg(method)
  if (method == "regression") method <- "Thurstone"

  S <- Lambda %*% Phi

  # solve(A, B) with a Moore-Penrose fallback when A is singular.
  solve_or_pinv <- function(A, B) {
    tryCatch(solve(A, B), error = function(e) .pinv(A) %*% B)
  }

  # Psi^-1 = diag(1 / (1 - h2)), the model-uniqueness weighting shared by the
  # Bartlett and Anderson methods. A Heywood case (a communality at or above 1)
  # leaves a non-positive uniqueness and makes the weighting undefined/unstable;
  # warn (psych proceeds silently) and continue.
  psi_inverse <- function() {
    u2 <- 1 - h2
    heywood <- u2 < .Machine$double.eps
    if (any(heywood)) {
      bad <- rownames(Lambda)[heywood]
      if (is.null(bad)) bad <- which(heywood)
      cli::cli_warn(
        c(paste("{cli::qty(bad)}Heywood case{?s} detected for {.val {bad}}:",
                "a communality at or above 1 leaves a non-positive uniqueness."),
          "i" = "The {method} factor-score weights for the affected variable{?s} are unstable."),
        class = "efa_scores_heywood"
      )
    }
    diag(1 / u2, nrow = length(u2))
  }

  W <- switch(
    method,
    Thurstone = solve_or_pinv(R, S),
    Bartlett = {
      psi_inv <- psi_inverse()
      psi_inv %*% Lambda %*%
        solve_or_pinv(t(Lambda) %*% psi_inv %*% Lambda, diag(ncol(Lambda)))
    },
    Anderson = {
      psi_inv <- psi_inverse()
      psi_inv %*% Lambda %*%
        .inv_mat_sqrt(t(Lambda) %*% psi_inv %*% R %*% psi_inv %*% Lambda)
    },
    tenBerge = {
      L <- Lambda %*% .mat_sqrt(Phi)
      r_inv_sqrt <- .inv_mat_sqrt(R)
      C <- r_inv_sqrt %*% L %*% .inv_mat_sqrt(crossprod(L, solve_or_pinv(R, L)))
      r_inv_sqrt %*% C %*% .mat_sqrt(Phi)
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
