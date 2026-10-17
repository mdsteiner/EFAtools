#' Schmid-Leiman transformation
#'
#' Transforms an oblique factor solution into a hierarchical solution with one general
#' factor and group factors that are uncorrelated with it and with each other (Schmid &
#' Leiman, 1957). With the same settings, the results agree with [psych::schmid()] and with
#' the SPSS syntax of Wolff and Preising (2005).
#'
#' @param x An oblique [efa_fit()] or [psych::fa()] solution, a pattern matrix (give its
#'   factor intercorrelations in `Phi`), or a second-order CFA solution from `lavaan`.
#' @param Phi matrix. The factor intercorrelations of the pattern matrix in `x`. With an
#'   [efa_fit()] or [psych::fa()] solution in `x`, `Phi` replaces the fitted factor
#'   intercorrelations. Not used for a `lavaan` solution.
#' @param estimator character. The estimator of the second-order [efa_fit()]: `"PAF"`
#'   (default), `"ML"`, or `"ULS"`. `"MINRES"` is a synonym for `"ULS"`. Not used for a
#'   `lavaan` solution and for two first-order factors, for which no second-order EFA is done.
#' @param g_name character. The name of the general factor in a `lavaan` solution. If
#'   `NULL` (default), the latent variable that loads on the other latent variables is used.
#' @param estimate_control An [estimate_control()] object with the settings of the
#'   second-order [efa_fit()]. `NULL` (default) uses the [efa_fit()] defaults. Not used for
#'   a `lavaan` solution and for two first-order factors.
#'
#' @details
#' The factor intercorrelations are factor analyzed to get one second-order (general)
#' factor. The general-factor loadings of the items are their first-order loadings times
#' the second-order loadings. The group-factor loadings are the first-order loadings times
#' the residual standard deviations of the first-order factors. The general factor and the
#' group factors thus divide the common variance, for example for McDonald's omegas (see
#' [efa_reliability()]).
#'
#' The transformation has proportionality constraints: in each group factor, the ratio of
#' the group-factor loading to the general-factor loading is the same for all items. Thus,
#' it recovers a bifactor structure exactly only if each item loads on one group factor and
#' these ratios are equal. If items have cross-loadings, it overestimates their
#' general-factor loadings, and the result depends on the oblique rotation of the
#' first-order solution (Mansolf & Reise, 2016). The orthogonal bifactor rotation of
#' [efa_fit()] (`rotation = "bifactorT"`) does not have these constraints.
#'
#' With two first-order factors, the second-order loadings are not identified. The two
#' loadings then get the same absolute value, and their product is the factor
#' intercorrelation, as in [psych::schmid()].
#'
#' A `lavaan` solution must have one second-order factor and no other paths or
#' covariances between the latent variables. Its general-factor loadings also include the
#' direct loadings of the general factor on the items, if the model has them.
#'
#' The group factors keep the column order and the names of the input solution. If the
#' columns do not have unique names (or a name is `"g"`, `"h2"`, or `"u2"`), all group
#' factors are labelled `"F1"` to `"Fk"`.
#'
#' @returns A list of class `c("efa_schmid_leiman", "SL")`:
#' \item{orig_R}{The correlation matrix of the input solution. `NA` for a pattern matrix.}
#' \item{sl}{The Schmid-Leiman loadings: the general factor (`g`), the group factors, the
#'   communalities (`h2`), and the uniquenesses (`u2`).}
#' \item{L2}{The second-order loadings, with one row for each first-order factor.}
#' \item{vars_accounted}{The variances explained by the general factor and the group
#'   factors.}
#' \item{iter}{The number of iterations of the second-order EFA.}
#' \item{convergence}{The convergence code of the second-order EFA (0 = converged; see
#'   [efa_fit()]).}
#' \item{settings}{The estimator and the estimation settings of the second-order EFA.}
#'
#' `iter`, `convergence`, and `settings` are `NA` if no second-order EFA was done, that is,
#' for a `lavaan` solution and for two first-order factors.
#'
#' @source Schmid, J. & Leiman, J. M. (1957). The development of hierarchical
#' factor solutions. Psychometrika, 22(1), 53–61. doi:10.1007/BF02289209
#' @source Wolff, H.-G., & Preising, K. (2005). Exploring item and higher order
#' factor structure with the Schmid-Leiman solution: Syntax codes for SPSS and
#' SAS. Behavior Research Methods, 37 , 48–58. doi:10.3758/BF03206397
#' @source Mansolf, M., & Reise, S. P. (2016). Exploratory bifactor analysis: The
#' Schmid-Leiman orthogonalization and Jennrich-Bentler analytic rotations. Multivariate
#' Behavioral Research, 51(5), 698–717. doi:10.1080/00273171.2016.1215898
#'
#' @family factor rotation
#' @family reliability coefficients
#'
#' @export
#'
#' @examples
#' efa_mod <- efa_fit(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                    estimator = "PAF", rotation = "promax")
#' efa_schmid_leiman(efa_mod)
#'
#' # A pattern matrix with its factor intercorrelations, for example from another
#' # program, gives the same solution
#' efa_schmid_leiman(efa_mod$rot_loadings, Phi = efa_mod$Phi)
#'
#' \donttest{
#' # A psych::fa() solution, with the psych settings for the second-order EFA
#' fa_mod <- psych::fa(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
#'                     fm = "pa", rotate = "Promax")
#' efa_schmid_leiman(fa_mod, estimate_control = estimate_control(type = "psych"))
#'
#' # A lavaan second-order CFA solution
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#'   mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
#'           F2 =~ V7 + V8 + V9 + V10 + V11 + V12
#'           F3 =~ V13 + V14 + V15 + V16 + V17 + V18
#'           g =~ F1 + F2 + F3'
#'   fit <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
#'                      sample.nobs = 500)
#'   efa_schmid_leiman(fit)
#' }
#' }
efa_schmid_leiman <- function(x, Phi = NULL,
                              estimator = c("PAF", "ML", "ULS", "MINRES"),
                              g_name = NULL, estimate_control = NULL) {

  .assert_args({
    checkmate::assert_matrix(Phi, null.ok = TRUE)
    checkmate::assert_string(g_name, null.ok = TRUE)
  })
  .assert_estimate_control(estimate_control)
  estimator <- .match_arg_ci(estimator)
  # "MINRES" is a synonym for "ULS" (same estimator); resolve to the canonical name.
  if (estimator == "MINRES") estimator <- "ULS"

  if (!inherits(x, c("EFA", "fa", "lavaan", "matrix", "LOADINGS", "loadings"))) {
    cli::cli_abort(
      c("{.arg x} must be an oblique {.fn efa_fit} or {.fn psych::fa} solution, a second-order {.cls lavaan} solution, or a pattern matrix with {.arg Phi}.",
        "i" = "From an {.cls efa_average} object, use the averaged loading matrix in
               {.code $loadings$average}, with the averaged factor correlations in
               {.code $Phi$average}."),
      class = "efa_sl_bad_input")
  }

  orig_R <- NA
  iter <- NA
  convergence <- NA
  settings <- NA

  if (inherits(x, "lavaan")) {

    .require_lavaan()

    if (isFALSE(lavaan::lavInspect(x, what = "converged"))) {
      cli::cli_abort("The model did not converge; the Schmid-Leiman transformation is not performed.",
                     class = "efa_sl_no_converge")
    }

    std_sol <- suppressWarnings(lavaan::lavInspect(x, what = "std"))

    # A multiple-group or multilevel fit gives one set of matrices for each group or level.
    if (!is.matrix(std_sol$lambda)) {
      cli::cli_abort(
        c("{.arg x} has more than one group or level.",
          "i" = "Fit the model to each group separately."),
        class = "efa_sl_multigroup")
    }

    if (anyNA(std_sol$lambda)) {
      cli::cli_abort("Some loadings are {.val NA} or {.val NaN}; the Schmid-Leiman transformation is not performed.",
                     class = "efa_sl_na_loadings")
    }

    if (any(diag(std_sol$theta) <= 0) || any(diag(std_sol$psi) <= 0)) {
      cli::cli_abort("A Heywood case was detected (a variance of 0 or negative); the Schmid-Leiman transformation is not performed.",
                     class = "efa_sl_heywood")
    }

    lv_names <- colnames(std_sol$lambda)
    if (!is.null(g_name) && !g_name %in% lv_names) {
      cli::cli_abort(
        c("The {.cls lavaan} solution has no latent variable {.val {g_name}}.",
          "i" = "Its latent variables are {.val {lv_names}}."),
        class = "efa_sl_g_name")
    }

    # The second-order factor is the latent variable that loads on other latent variables.
    # The transformation needs exactly one, and first-order factors whose disturbances are
    # uncorrelated, so that the group factors are uncorrelated.
    tol <- .Machine$double.eps * 100
    beta <- std_sol$beta
    higher <- if (is.null(beta)) character() else colnames(beta)[colSums(abs(beta) > tol, na.rm = TRUE) > 0]
    if (length(higher) > 1L) {
      cli::cli_abort(
        c("More than one latent variable loads on other latent variables: {.val {higher}}.",
          "i" = "{.fn efa_schmid_leiman} needs a second-order model with one second-order factor."),
        class = "efa_sl_higher_order")
    }
    if (is.null(g_name)) g_name <- higher

    if (!isTRUE(g_name %in% higher)) {
      cli::cli_abort(
        c("{.arg x} does not appear to be a second-order CFA solution.",
          "i" = "{.fn efa_schmid_leiman} needs a second-order model; pass a bifactor solution directly to {.fn efa_reliability}."),
        class = "efa_sl_not_second_order"
      )
    }

    psi <- std_sol$psi
    if (any(abs(psi[upper.tri(psi)]) > tol, na.rm = TRUE)) {
      cli::cli_abort(
        c("The {.cls lavaan} solution has correlated latent variables or disturbances.",
          "i" = "{.fn efa_schmid_leiman} needs a second-order model whose first-order factors are related only through the second-order factor."),
        class = "efa_sl_correlated_factors")
    }

    fac_names <- setdiff(lv_names, g_name)
    labels <- .sl_group_names(fac_names, length(fac_names))
    L1 <- std_sol$lambda[, fac_names, drop = FALSE]
    L2 <- std_sol$beta[fac_names, g_name, drop = FALSE]

    # The general factor can also load on items directly, so its loadings are that direct
    # path plus the paths through the first-order factors.
    L_g <- std_sol$lambda[, g_name] + L1 %*% L2
    L_s <- .sl_group_loadings(L1, std_sol$psi, fac_names)

    S <- lavaan::lavInspect(x, what = "sampstat")$cov
    if (!is.null(S)) orig_R <- stats::cov2cor(unclass(S))[rownames(L1), rownames(L1)]

  } else {

    if (inherits(x, c("EFA", "fa"))) {

      if (inherits(x, "EFA") && identical(x$settings$rotation, "bifactorQ")) {
        cli::cli_abort(
          c("{.arg x} is a bifactor solution ({.code rotation = \"bifactorQ\"}), which already has a general factor.",
            "i" = "Transform an oblique solution without a general factor, or refit with {.code rotation = \"bifactorT\"} for an orthogonal bifactor solution."),
          class = "efa_sl_bifactor")
      }

      if (is.null(x[["Phi"]])) {
        cli::cli_abort("{.arg x} is a non-rotated or orthogonal factor solution, but SL needs an oblique solution.",
                       class = "efa_sl_not_oblique")
      }

      if (inherits(x, "EFA")) {
        L1 <- x$rot_loadings
        orig_R <- x$orig_R
      } else {
        L1 <- x$loadings
        orig_R <- unclass(x$r)
      }
      if (is.null(Phi)) Phi <- x[["Phi"]]

    } else {

      if (is.null(Phi)) {
        cli::cli_abort(
          c("{.arg Phi} was not provided.",
            "i" = "Enter an oblique solution from {.fn efa_fit} or {.fn psych::fa}, a second-order CFA from lavaan, or provide {.arg Phi}."),
          class = "efa_sl_phi_missing"
        )
      }
      L1 <- x

    }

    L1 <- unclass(L1)
    k <- ncol(L1)

    # The transformation regresses the first-order factors on one second-order factor, so
    # it needs at least two of them.
    if (k < 2L) {
      cli::cli_abort(
        c("A Schmid-Leiman transformation needs at least two first-order factors.",
          "x" = "{.arg x} has {k} first-order factor{?s}."),
        class = "efa_sl_too_few_factors"
      )
    }

    labels <- .sl_group_names(colnames(L1), k)
    Phi <- .align_correlation_axis(Phi, n = k, target_names = colnames(L1), arg = "Phi")
    # so that the messages of the second-order EFA name the first-order factors
    dimnames(Phi) <- list(labels, labels)

    if (k == 2L) {

      # One factor on two variables is not identified. As in psych::schmid(), the two
      # loadings get the same absolute value, and their product is the intercorrelation.
      cli::cli_warn(
        c("With two first-order factors, the second-order loadings are not identified.",
          "i" = "They get the same absolute value, and their product is the factor intercorrelation. Interpret the general factor with caution."),
        class = "efa_sl_underidentified")
      r12 <- Phi[1, 2]
      L2 <- matrix(sqrt(abs(r12)) * c(1, sign(r12)), 2L, 1L)

    } else {

      # One factor on three intercorrelations is just identified for every input, so that
      # warning says nothing about the data; the other warnings of this EFA (e.g., a Heywood
      # case) reach the user. N only avoids the warning about a missing sample size.
      fit <- withCallingHandlers(
        efa_fit(Phi, n_factors = 1, N = 100, estimator = estimator, rotation = "none",
                estimate_control = estimate_control),
        efa_just_identified = function(w) invokeRestart("muffleWarning"))

      iter <- fit$iter
      convergence <- fit$convergence
      # The sample size and the inference settings of this EFA are placeholders.
      settings <- fit$settings[setdiff(
        names(fit$settings),
        c("method", "rotation", "n_factors", "N", "use", "cor_method", "input_type",
          "cor_method_used", .fit_inference_args))]
      L2 <- unclass(fit$unrot_loadings)

    }

    # A second-order communality at or above 1 is a Heywood case: the residual variance of
    # that first-order factor, and thus its group-factor loadings, are undefined.
    comm_h <- L2[, 1]^2
    if (any(comm_h >= 1 + .Machine$double.eps)) {
      cli::cli_abort(
        c("A Heywood case was detected in the second-order factor analysis; no Schmid-Leiman solution is computed.",
          "i" = "A second-order communality is 1 or larger, so the residualized first-order loadings are undefined."),
        class = "efa_sl_heywood"
      )
    }

    L_g <- L1 %*% L2
    # pmax() clamps a communality within rounding of 1, so sqrt() never sees a tiny
    # negative value
    L_s <- L1 %*% diag(sqrt(pmax(0, 1 - comm_h)), nrow = k)

  }

  dimnames(L2) <- list(labels, "g")
  sl_load <- cbind(L_g, L_s)
  dimnames(sl_load) <- list(rownames(L1), c("g", labels))
  h2 <- rowSums(sl_load^2)

  sl <- cbind(sl_load, h2 = h2, u2 = 1 - h2)
  class(sl) <- c("efa_sl_loadings", "SLLOADINGS")

  output <- list(
    orig_R = orig_R,
    sl = sl,
    L2 = L2,
    vars_accounted = .compute_vars(L_unrot = sl_load, L_rot = sl_load),
    iter = iter,
    convergence = convergence,
    settings = settings
  )

  class(output) <- c("efa_schmid_leiman", "SL")

  output

}

# Labels of the group factors: the input's column names if each factor has a unique,
# non-empty one that does not clash with the "g", "h2", and "u2" columns of the result;
# otherwise "F1" to "Fk".
.sl_group_names <- function(nms, k) {
  usable <- length(nms) == k && !anyNA(nms) && all(nzchar(nms)) &&
    anyDuplicated(nms) == 0L && !any(nms %in% c("g", "h2", "u2"))
  if (usable) nms else paste0("F", seq_len(k))
}
