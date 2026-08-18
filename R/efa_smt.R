# Number of factors a sequential SMT rule suggests. `values` runs over the 0, 1, ...,
# max_fac factor models and `rule` maps that series to a logical "this model satisfies
# the rule" (taken as a function rather than a second vector so the series cannot be
# passed twice inconsistently). The search stops at the first model that satisfies the
# rule, or at the first undefined value: an NA breaks the strictly sequential test, so
# the search stops there rather than skipping ahead to a later model that satisfies it,
# and no number can be suggested. Shared by the chi-square and the RMSEA rule so the
# two cannot drift apart.
.smt_sequential_stop <- function(values, rule) {
  stop_at <- which(rule(values) | is.na(values))[1]
  if (is.na(stop_at) || is.na(values[stop_at])) NA_integer_ else stop_at - 1L
}

#' Sequential chi square model tests, RMSEA lower bound, and AIC
#'
#' Sequential chi square model tests (SMT) are a factor retention method where
#' multiple
#' EFAs with increasing numbers of factors are fitted and the number of factors
#' for which the Chi Square value first becomes non-significant is taken as the
#' suggested number of factors.
#' Preacher, Zhang, Kim, & Mels (2013) suggested a similar approach with the
#' lower bound of the 90% confidence interval of the Root Mean Square Error of
#' Approximation (RMSEA; Browne & Cudeck, 1992; Steiger & Lind, 1980), and with
#' the Akaike Information Criterion (AIC). For the RMSEA, the
#' number of factors for which this lower bound first falls below .05 is the
#' suggested number of factors to retain. For the AIC, it is the number of factors
#' where the AIC is lowest.
#'
#' @inheritParams efa_kgc
#' @param cor_method character. One of `"pearson"`, `"spearman"`, or `"kendall"`,
#'   passed to [stats::cor()]. `"poly"` and `"tetra"` are not supported because
#'   `SMT` rests on a normal-theory chi-square test that is not valid for
#'   polychoric / tetrachoric correlations.
#'  Default is `"pearson"`.
#' @param N numeric. The number of observations. Needs only be specified if a
#' correlation matrix is used. Must be larger than the number of variables.
#' @param estimate_control an [estimate_control()] object with the estimation settings for the
#'  sequential [efa_fit()] fits. `NULL` (default) uses the [efa_fit()] defaults. The sequential
#'  models are fitted with maximum likelihood (the chi-square, RMSEA, and AIC the SMT is built on
#'  are defined for it), so of the estimation knobs only `start_method` takes effect; the ones
#'  governing principal axis factoring do not apply. The models are unrotated, so no rotation
#'  settings apply either.
#'
#' @details
#' As a first step in the procedure, a maximum number of factors to extract is
#' determined for which the model is still over-identified (df > 0).
#'
#' Then, EFAs with increasing numbers of factors from 1 to the maximum number are
#' fitted with maximum likelihood estimation.
#'
#' For the SMT, first the significance of the chi
#' square value for a model with 0 factors is determined. If this value is
#' not significant, 0 factors are suggested to retain. If it is significant,
#' a model with 1 factor is estimated and the significance of its chi square value
#' is determined, and so on, until a non-significant result is obtained. The
#' suggested number of factors is the number of factors for the model where the
#' chi square value first becomes non-significant.
#'
#' Regarding the RMSEA, the suggested number of factors is the number of factors
#' for the model where the lower bound of the 90% confidence interval of the
#' RMSEA first falls below the .05 threshold.
#'
#' Regarding the AIC, the suggested number of factors is the number of factors
#' for the model with the lowest AIC.
#'
#' The sequential models are fitted without inequality constraints, so a solution
#' can be inadmissible (a Heywood case, or a fit that did not converge). Only the
#' models the three rules actually select are checked for this; if one of them is
#' inadmissible a warning is raised and the corresponding suggestion should be
#' interpreted with caution.
#'
#' In comparison with other prominent factor retention criteria, SMT performed
#' well at determining the number of factors to extract in EFA (Auerswald &
#' Moshagen, 2019). The RMSEA lower bound also performed well at determining the true
#' number of factors, while the AIC performed well at determining the
#' most generalizable model (Preacher, Zhang, Kim, & Mels, 2013).
#'
#' @returns An object of class `efa_retention` (see [print.efa_retention()] for
#'   the print method). SMT has no plot; [plot.efa_retention()] returns `NULL`
#'   with a message for it. Its main fields are:
#' \item{n_factors}{A named numeric vector (`"chi"`, `"RMSEA"`, `"AIC"`) with the
#'   suggested number of factors from the sequential chi-square model tests, the
#'   RMSEA lower bound, and the AIC.}
#' \item{results}{A list with one record per criterion, each holding the criterion
#'   values for the null model (zero factors) through the maximum number of
#'   factors.}
#' \item{settings}{A list of the settings used.}
#'
#' @source Auerswald, M., & Moshagen, M. (2019). How to determine the number of
#' factors to retain in exploratory factor analysis: A comparison of extraction
#' methods under realistic conditions. Psychological Methods, 24(4), 468–491.
#' https://doi.org/10.1037/met0000200
#' @source Browne, M.W., & Cudeck, R. (1992). Alternative ways of assessing model
#' fit. Sociological Methods and Research, 21, 230–258.
#' @source Preacher, K. J., Zhang G., Kim, C., & Mels, G. (2013). Choosing the
#' Optimal Number of Factors in Exploratory Factor Analysis: A Model Selection
#' Perspective, Multivariate Behavioral Research, 48(1), 28-56,
#' doi:10.1080/00273171.2012.710386
#' @source Steiger, J. H., & Lind, J. C. (1980, May). Statistically based tests
#' for the number of common factors. Paper presented at the annual meeting of
#' the Psychometric Society, Iowa City, IA.
#'
#' @family factor retention criteria
#'
#' @seealso [efa_retain()] as a wrapper function for this and the other factor
#'   retention criteria.
#'
#' @export
#'
#' @examples
#' SMT_base <- efa_smt(test_models$baseline$cormat, N = 500)
#' SMT_base
#'
efa_smt <- function(x, N = NA, use = c("pairwise.complete.obs", "all.obs",
                                     "complete.obs", "everything",
                                     "na.or.complete"),
                cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                estimate_control = NULL){

  # Perform argument checks
  .assert_cor_input(x)

  checkmate::assert_count(N, na.ok = TRUE)
  use <- .match_arg_ci(use)
  cor_method <- .match_arg_ci(cor_method)
  .reject_poly_reference(
    cor_method, "efa_smt",
    why = "{.fn {fn}} rests on a normal-theory chi-square test that is not valid for such correlations.")
  .assert_estimate_control(estimate_control)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(
    x, N = N, use = use, cor_method = cor_method, N_policy = "required",
    N_required_msg = c("{.arg N} is {.val NA}.",
                       "i" = "Provide {.arg N} or raw data."))
  R <- prep$R
  N <- prep$N
  .assert_n_gt_vars(N, ncol(R))

  # Prepare objects for sequential tests
  max_fac <- .det_max_factors(ncol(R))

  if(max_fac <= 0){
    cli::cli_abort(
      c("The model is underidentified, or just identified with a single factor; SMTs cannot be performed.",
        "i" = "Provide more indicators."),
      class = "efa_smt_underidentified"
    )
  }

  ps <- vector("double", max_fac)
  RMSEA_LB <- vector("double", max_fac)
  AIC <- vector("double", max_fac)
  # admissibility of each fitted model; only the selected ones are reported on below
  n_heywood <- vector("integer", max_fac)
  converged <- vector("logical", max_fac)

    # sequentially perform EFAs with 1 to the maximum number of factors
    for (i in seq_len(max_fac)) {

      temp <- suppressWarnings(suppressMessages(
        efa_fit(R, n_factors = i, estimator = "ML", rotation = "none", N = N,
                estimate_control = estimate_control)))
      ps[i] <- temp$fit_indices$p_chi
      RMSEA_LB[i] <- temp$fit_indices$RMSEA_LB
      AIC[i] <- temp$fit_indices$AIC
      n_heywood[i] <- length(temp$heywood)
      converged[i] <- !isTRUE(temp$convergence != 0)

    }

  # With which number of factors does the chi square first become
  # non-significant?

  # Null-model (zero-factor) statistics. These depend only on R and N (the
  # model-implied matrix is the identity), so compute them directly rather than
  # reading them off the last fitted EFA, whose fit_indices come back NA when
  # that model is degenerate (e.g. a Heywood / non-positive-definite case).
  # Mirrors the null-model block in .gof().
  m <- ncol(R)
  # Reuse one log-determinant of R for both null chi-squares (the Bartlett-corrected statistic
  # below and the uncorrected RMSEA baseline further down) instead of factorising R twice.
  ld_R <- determinant(R, logarithm = TRUE)
  chi_null <- .null_chisq(R, N, ld = ld_R)
  df_null <- (m^2 - m) / 2
  p_null <- stats::pchisq(chi_null, df_null, lower.tail = FALSE)

  # Retain the number of factors of the first model whose chi square is
  # non-significant. An NA p (e.g. a non-converged / Heywood model) stops the search
  # with no suggestion, as does a sequence in which every model is significant.
  p_seq <- c(p_null, ps)
  nfac_chi <- .smt_sequential_stop(p_seq, \(p) p > 0.05)

  # Calculate RMSEA (incl. lower bound of 90% CI) and AIC for the null model. The
  # RMSEA noncentrality is built on the uncorrected (N - 1) discrepancy scale (as in
  # .chi_fit_indices()), so the null model uses the uncorrected baseline chi-square;
  # p_null and AIC_null keep the Bartlett-corrected statistic computed above. When N is
  # too small for the Bartlett correction (chi_null is NA) the null-model asymptotics
  # break down, so leave the RMSEA bound NA too, matching how .gof() drops the
  # fitted-model fit indices at such N.
  chi_null_rmsea <- .null_chisq(R, N, ld = ld_R, corrected = FALSE)
  RMSEA_LB_null <- if (is.na(chi_null)) {
    NA_real_
  } else {
    sqrt(.rmsea_lambda(chi_null_rmsea, df_null, .95) / (df_null * (N - 1)))
  }

  AIC_null <- chi_null - 2 * df_null

  # Retain the number of factors of the first model whose RMSEA lower bound falls
  # below .05, under the same sequential convention as the chi-square rule above.
  RMSEA_seq <- c(RMSEA_LB_null, RMSEA_LB)
  nfac_RMSEA <- .smt_sequential_stop(RMSEA_seq, \(r) r < .05)

  # With which number of factors is the AIC lowest? (which.min returns the first
  # minimum, so ties yield a single, well-defined suggestion). which.min()
  # returns integer(0) when every AIC is NA (a degenerate tiny-N case where the
  # null and all factor models are undefined); fall back to NA so the record
  # keeps a length-1 suggestion.
  AIC_all <- c(AIC_null, AIC)
  nfac_AIC <- if (all(is.na(AIC_all))) NA_integer_ else which.min(AIC_all) - 1

  # The models at the upper end of the sequence over-extract and may be
  # inadmissible, so flag only the solutions the three rules actually select (as
  # efa_hull() does). The null model is not fitted and cannot be inadmissible.
  suggested <- c(chi = nfac_chi, RMSEA = nfac_RMSEA, AIC = nfac_AIC)
  inadmissible <- character(0)
  for (rule in names(suggested)) {
    k <- suggested[[rule]]
    if (is.na(k) || k < 1) next
    issues <- c(if (n_heywood[k] > 0) "Heywood case",
                if (!converged[k]) "non-convergence")
    if (length(issues) > 0) {
      inadmissible <- c(
        inadmissible,
        paste0(rule, ": ", .retention_count(k), " factor",
               if (k != 1) "s" else "",
               " (", paste(issues, collapse = ", "), ")")
      )
    }
  }

  if (length(inadmissible) > 0) {
    cli::cli_warn(
      c("The sequential model tests selected an inadmissible solution: {inadmissible}.",
        "i" = "The selected solution has a Heywood case or did not converge, so the suggested number of factors may be unreliable; interpret it with caution and cross-check with other criteria."),
      class = "efa_smt_inadmissible"
    )
  }

  # one record per criterion (values for the null model through max_fac factors)
  results <- list(
    list(name = "chi", label = "Sequential chi-square model tests",
         n_factors = nfac_chi, plot_type = "none",
         x = 0:max_fac, y = c(p_null, ps)),
    list(name = "RMSEA", label = "Lower bound of RMSEA 90% CI",
         n_factors = nfac_RMSEA, plot_type = "none",
         x = 0:max_fac, y = c(RMSEA_LB_null, RMSEA_LB)),
    list(name = "AIC", label = "Akaike Information Criterion",
         n_factors = nfac_AIC, plot_type = "none",
         x = 0:max_fac, y = AIC_all)
  )

  output <- .new_efa_retention(
    "SMT",
    results = results,
    settings = list(N = N, use = use, cor_method = cor_method)
  )

  return(output)

}
