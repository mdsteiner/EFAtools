#' Screen data for exploratory factor analysis
#'
#' @description
#' Checks whether your data are suitable for exploratory factor analysis (EFA). From a
#' correlation matrix or raw data, it reports the Kaiser-Meyer-Olkin (KMO) measure of
#' sampling adequacy, Bartlett's test of sphericity, the determinant and condition number
#' of the correlation matrix, and each variable's squared multiple correlation (SMC). When
#' you supply raw data, it also reports each variable's variance and percentage of missing
#' values, category counts for categorical variables, tests of multivariate normality, and
#' multivariate outliers.
#'
#' @param x data.frame or matrix. Raw data, or a correlation matrix. Needs at least three
#'   variables, none of which is a perfect linear combination of the others.
#' @param N numeric. The number of observations. Set this only when you supply a
#'   correlation matrix; it is needed for Bartlett's test of sphericity and is taken from
#'   the data automatically when you supply raw data. Default is `NA`.
#' @param use character. How to handle missing values in raw data. For
#'   `cor_method = "pearson"`, `"spearman"`, or `"kendall"` this is passed to
#'   [stats::cor()]. For `"poly"` or `"tetra"` the same rule is applied to the raw data
#'   before the correlations are estimated; `"all.obs"` and `"everything"` then stop with
#'   an error on any missing value, instead of returning `NA` correlations. Default is
#'   `"pairwise.complete.obs"`.
#' @param cor_method character. How to compute correlations from raw data: `"pearson"`,
#'   `"spearman"`, or `"kendall"` (via [stats::cor()]), or `"poly"` / `"tetra"` for
#'   polychoric / tetrachoric correlations of ordinal or binary data. A Spearman or Kendall
#'   correlation matrix is screened on its own scale, not converted to look like a Pearson
#'   correlation matrix; Kendall's tau in particular measures something different from a
#'   Pearson correlation, not just a rescaled version of it. Default is `"pearson"`.
#' @param mcd_alpha numeric. The proportion of cases used to build the robust outlier
#'   estimate, between 0.5 and 1. The default, `0.5`, is the most robust choice; a larger
#'   value uses more of the data but resists outliers less well. Used only with raw data.
#' @param outlier_cutoff numeric. The probability used to set the cutoff for flagging a
#'   multivariate outlier, between 0.5 and 0.9999. Default is `0.975`. Used only with raw
#'   data.
#' @param seed integer. A seed for the random subsets used by the outlier detection, so
#'   the result is reproducible. Does not affect your random-number generator elsewhere.
#'   Default is `NULL`. Used only with raw data.
#'
#' @details
#' The diagnostics are computed from the analysis correlation matrix \eqn{R}:
#' \describe{
#'   \item{KMO}{The Kaiser-Meyer-Olkin measure of sampling adequacy (Kaiser, 1970; Kaiser
#'     & Rice, 1974), overall and for each variable; see [efa_kmo()]. It shows how much
#'     common variance your variables share. Higher values are better; a common rule of
#'     thumb treats values below .50 as unacceptable.}
#'   \item{Bartlett}{Bartlett's (1951) test of sphericity: the likelihood-ratio test of
#'     whether the correlation matrix is an identity matrix, i.e., whether your variables
#'     correlate with each other at all; see [efa_bartlett()]. A significant result
#'     supports doing a factor analysis. The test needs the sample size `N`; without it,
#'     this diagnostic is skipped with a warning and `$bartlett` is `NULL`. If `N` is too
#'     small relative to the number of variables, the statistic is `NA`, also with a
#'     warning.}
#'   \item{Determinant}{The determinant of \eqn{R}, reported as a number only. It falls as
#'     you add variables even when the variables are not collinear, so a fixed cut-off on
#'     it (such as the 0.00001 often quoted from Field, 2018) says more about how many
#'     variables you have than about your data. Use the condition number instead.}
#'   \item{Condition number}{The ratio of the largest to the smallest eigenvalue of
#'     \eqn{R}. Its square root, the condition index, is the collinearity diagnostic of
#'     Belsley, Kuh & Welsch (1980); it drives the printed report and its recommendation.
#'     An index of 10 or less is rarely of interest. An index above 30 flags a near linear
#'     dependency: two or more variables that together carry much the same information.
#'     An index between the two is not negligible, but it stays below the value that
#'     flags a dependency. Belsley (1991) gives 30 as one example value and calls the
#'     choice of a cut-off "somewhat of an art form", so the report grades an index above
#'     30 by its position on the scale 1, 3, 10, 30, 100, 300, 1000: moderate (30 to 100),
#'     strong (100 to 300), or very strong (above 300). These values come from regression
#'     diagnostics on data that are not centred, but a correlation matrix is centred, so
#'     use them as a guide and not as a test.}
#'   \item{SMC}{The squared multiple correlation of each variable with all the others. A
#'     low value flags a variable that has little in common with the rest of your set.}
#'   \item{Variance and missing data}{For raw data: each variable's variance (over its
#'     available values) and percentage of missing values, computed from every row you
#'     supplied. These missing-value percentages explain why the correlation matrix's
#'     sample size (`N`) can be smaller than the number of rows in your data. Ordered-factor
#'     columns are recoded to integer levels first, so `variance` reflects those codes.}
#'   \item{Categories}{For raw data: for each variable with fewer than ten distinct values
#'     (treated as categorical), the count of responses in each category. A category with
#'     fewer than five responses is flagged as sparse, and an unused category between the
#'     smallest and largest response is flagged as empty. As a rule of thumb, items with
#'     fewer than five response categories are better analysed with `cor_method = "poly"`
#'     or `"tetra"` than with Pearson correlations (Rhemtulla et al., 2012).}
#'   \item{Multivariate normality}{For raw data, using only complete cases: two tests of
#'     multivariate normality, Mardia's (1970) test of skewness and kurtosis and the
#'     Henze-Zirkler (1990) test. A small p-value on either test suggests your data depart
#'     from a multivariate normal distribution, a reason to prefer a robust or ordinal
#'     method over normal-theory maximum likelihood. In a very small sample the kurtosis
#'     statistic is `NA`. The Henze-Zirkler p-value is not available with more than about
#'     50 to 60 variables; its test statistic is still reported.}
#'   \item{Outliers}{For raw data, using only complete cases: multivariate outliers, found
#'     from a robust estimate of each case's distance from the centre of your data (the
#'     minimum covariance determinant method; Rousseeuw & Van Driessen, 1999). A flagged
#'     case is unusually far from the rest of your sample. When there are too few complete
#'     cases, the variables are too collinear, or too many cases share identical answers,
#'     a plain (non-robust) distance is used instead, with a warning explaining why.}
#' }
#'
#' @returns An object of class `efa_screen`, a list containing:
#' \item{kmo}{A list with the overall KMO (`KMO`) and the per-variable KMO (`KMO_i`).}
#' \item{bartlett}{A list with Bartlett's chi-square statistic (`chisq`), its `p_value`,
#'   and its degrees of freedom (`df`); `chisq` and `p_value` are `NA` when `N` was too
#'   small for the correction. `NULL` when `N` is unavailable.}
#' \item{determinant}{The determinant of the correlation matrix.}
#' \item{condition}{The condition number of the correlation matrix (largest eigenvalue
#'   over smallest).}
#' \item{smc}{The per-variable squared multiple correlations.}
#' \item{per_item}{A data frame with one row per variable (row names are the variable
#'   names): `variance`, `missing` (percentage), `smc`, `kmo_i`, and `flags` (any
#'   sparse/empty-category issues). `NULL` when a correlation matrix is supplied instead
#'   of raw data.}
#' \item{normality}{A list with `mardia` (skewness `skewness`, `skewness_df`,
#'   `skewness_p`, kurtosis `kurtosis` and `kurtosis_p`, and the underlying `b1p`/`b2p`),
#'   `hz` (the Henze-Zirkler `statistic` and its `p_value`), and `n_complete` (the number
#'   of complete cases used). `NULL` without raw data, or a note explaining why when the
#'   complete-case data cannot support the tests.}
#' \item{outliers}{A list with `distances` (each complete case's robust distance, named by
#'   its row number), `cutoff` (the flagging threshold, on the same scale as `distances`),
#'   `flagged` (the row numbers exceeding `cutoff`), `center` and `cov` (the robust
#'   location and scatter), `method` (`"mcd"` or the `"classical"` fallback),
#'   `fallback_reason` (why the robust estimate was unavailable, when `method` is
#'   `"classical"`), and `n_complete`. `NULL` without raw data, or a note explaining why
#'   when no covariance can be formed.}
#' \item{categories}{A named list with the response-category counts for each categorical
#'   variable (in category order); `NA` for a variable treated as continuous. `NULL`
#'   without raw data.}
#' \item{note}{Explains why the raw-data diagnostics (`per_item`, `normality`,
#'   `outliers`, `categories`) are missing, when a correlation matrix is supplied instead
#'   of raw data. `NULL` when raw data are supplied.}
#' \item{settings}{The settings used: `N`, `n_obs` (rows in the raw data supplied, `NA`
#'   for a correlation-matrix input), `use`, `cor_method`, `mcd_alpha`, `outlier_cutoff`,
#'   and `seed`.}
#'
#' @source Bartlett, M. S. (1951). The effect of standardization on a Chi-square
#'   approximation in factor analysis. Biometrika, 38, 337-344.
#' @source Belsley, D. A. (1991). A guide to using the collinearity diagnostics. Computer
#'   Science in Economics and Management, 4, 33-50.
#' @source Belsley, D. A., Kuh, E. & Welsch, R. E. (1980). Regression diagnostics:
#'   Identifying influential data and sources of collinearity. Wiley.
#' @source Cochran, W. G. (1954). Some methods for strengthening the common
#'   \eqn{\chi^2} tests. Biometrics, 10, 417-451.
#' @source Croux, C. & Haesbroeck, G. (1999). Influence function and efficiency of the
#'   minimum covariance determinant scatter matrix estimator. Journal of Multivariate
#'   Analysis, 71, 161-190.
#' @source Field, A. (2018). Discovering statistics using IBM SPSS statistics (5th
#'   ed.). Sage.
#' @source Henze, N. & Zirkler, B. (1990). A class of invariant consistent tests for
#'   multivariate normality. Communications in Statistics - Theory and Methods, 19,
#'   3595-3617.
#' @source Kaiser, H. F. (1970). A second generation little jiffy. Psychometrika,
#'   35, 401-415.
#' @source Kaiser, H. F. & Rice, J. (1974). Little jiffy, mark IV. Educational and
#'   Psychological Measurement, 34, 111-117.
#' @source Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with
#'   applications. Biometrika, 57, 519-530.
#' @source Mardia, K. V. (1974). Applications of some measures of multivariate skewness
#'   and kurtosis in testing normality and robustness studies. Sankhya B, 36, 115-128.
#' @source Pison, G., Van Aelst, S. & Willems, G. (2002). Small sample corrections for LTS
#'   and MCD. Metrika, 55, 111-123.
#' @source Rhemtulla, M., Brosseau-Liard, P. E. & Savalei, V. (2012). When can
#'   categorical variables be treated as continuous? A comparison of robust continuous
#'   and categorical SEM estimation methods under suboptimal conditions. Psychological
#'   Methods, 17, 354-373.
#' @source Rousseeuw, P. J. & Van Driessen, K. (1999). A fast algorithm for the minimum
#'   covariance determinant estimator. Technometrics, 41, 212-223.
#'
#' @seealso [efa_kmo()] and [efa_bartlett()] for the individual suitability measures, and
#'   [efa_retain()] for factor retention criteria.
#'
#' @family factor analysis suitability
#'
#' @export
#'
#' @examples
#' # From a correlation matrix (supply N for Bartlett's test of sphericity)
#' efa_screen(test_models$baseline$cormat, N = 500)
#'
#' # From raw data (N is taken from the data; the seed makes the outlier
#' # diagnostics reproducible)
#' efa_screen(GRiPS_raw, seed = 1)
#'
efa_screen <- function(x, N = NA,
                       use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                               "everything", "na.or.complete"),
                       cor_method = c("pearson", "spearman", "kendall", "poly",
                                      "tetra"),
                       mcd_alpha = 0.5, outlier_cutoff = 0.975, seed = NULL) {

  # Perform argument checks
  .assert_cor_input(x)

  # Below three variables the report is degenerate rather than merely uninformative: the
  # KMO is 0/0 at p = 1 and identically .5 at p = 2 whatever the correlation, and
  # Bartlett's test has no or one degree of freedom, so the sections would contradict
  # each other. Refuse it instead of issuing a confident verdict on nothing.
  if (ncol(x) < 3L) {
    cli::cli_abort(
      c("Screening data for factor analysis needs at least three variables.",
        "x" = "{.arg x} has {ncol(x)} variable{?s}.",
        "i" = "The sampling adequacy of one or two variables is not defined: the KMO is
               undefined for a single variable and exactly {.val {0.5}} for any pair."),
      class = "efa_screen_too_few_vars")
  }

  use <- .match_arg_ci(use)
  cor_method <- .match_arg_ci(cor_method)
  checkmate::assert_count(N, na.ok = TRUE)
  checkmate::assert_number(mcd_alpha, lower = 0.5, upper = 1)
  # The endpoints of the natural [0, 1] range are degenerate rather than extreme: a
  # cutoff of 0 flags every observation and a cutoff of 1 puts the threshold at
  # infinity, so neither yields an outlier diagnostic. Bound the probability to the
  # usable range instead.
  checkmate::assert_number(outlier_cutoff, lower = 0.5, upper = 0.9999)
  checkmate::assert_int(seed, null.ok = TRUE)

  # Retain the raw data (when supplied) for the raw-data screening diagnostics;
  # NULL for a correlation-matrix input. data.matrix() codes factor / character
  # columns to their integer category codes, matching how .polychoric() reads
  # ordinal data (as.matrix() would coerce the whole frame to character as soon as
  # any one column is a factor, corrupting every column). data.matrix() only recodes
  # the columns of a data frame, so a non-numeric matrix (e.g. string-coded ordinal
  # responses) is passed through as.data.frame() first to receive the same coding.
  xr <- if (.is_cormat(x)) {
    NULL
  } else if (is.matrix(x) && !is.numeric(x)) {
    data.matrix(as.data.frame(x))
  } else {
    data.matrix(x)
  }

  # The level labels behind that coding, so the category tabulation can name its
  # counts by the responses the user supplied rather than by their integer codes.
  xl <- if (is.null(xr)) NULL else .screen_levels(x)

  # Detect or compute the correlation matrix, check it, and smooth it if needed.
  # N_policy = "optional" keeps N as NA for a correlation matrix without N (only
  # Bartlett's test needs it), and inform_from_data = FALSE suppresses the
  # "computing correlations from the raw data" note.
  # A singular matrix stays a refusal -- every suitability measure here is built from
  # R^-1, so there is nothing trustworthy to report -- but redundant items are the most
  # common reason to screen in the first place, so name them instead of leaving the user
  # to find them. The near-singular case is untouched and still yields a full report with
  # its multicollinearity flags.
  # The handler runs from tryCatch()'s own frame, so the re-raise is given this call
  # explicitly; otherwise the abort would name the internal helper rather than
  # efa_screen(), as every other condition here does.
  screen_call <- rlang::current_env()
  prep <- tryCatch(
    .prepare_cor_input(x, N = N, use = use, cor_method = cor_method,
                       N_policy = "optional", inform_from_data = FALSE),
    efa_cor_singular = function(cnd) .screen_singular_abort(cnd, call = screen_call))
  R <- prep$R
  N <- prep$N
  p <- ncol(R)

  # Variable names for the per-variable diagnostics, with a V-fallback when the
  # input is unnamed (matching KMO()).
  nms <- colnames(R)
  if (is.null(nms)) nms <- rownames(R)
  if (is.null(nms)) nms <- paste0("V", seq_len(p))

  # Two factorisations of R serve every correlation-side measure below, each formed once
  # rather than once per measure: they cost O(p^3), which is the bulk of the work when a
  # large correlation matrix is screened.
  #
  # The eigenvalues give the condition number and the log-determinant, both taken in
  # absolute value: the singular values of a symmetric matrix are its absolute eigenvalues,
  # so max|lambda| / min|lambda| is exactly the 2-norm condition number, and
  # prod(sign(lambda)) with sum(log|lambda|) are exactly the sign and log-modulus of the
  # determinant. The raw ratio and the raw logarithm would be defined only for a positive
  # definite R.
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values

  # R^-1, needed by the KMO's anti-image matrix and by the squared multiple correlations,
  # from the Cholesky factor, which uses the symmetry and positive definiteness that
  # solve()'s LU ignores. .prepare_cor_input() above has refused a numerically singular R
  # and smoothed an indefinite one, so the factorisation exists. Should a Cholesky pivot
  # still turn non-positive, the inverse is left unset and each measure falls back to the
  # route it takes when it is given no inverse -- which keeps .compute_kmo()'s classed
  # singularity abort rather than replacing it with a bare solve() error here.
  R_inv <- tryCatch(chol2inv(chol(R)), error = function(e) NULL)

  # KMO measure of sampling adequacy, overall and per variable
  kmo <- .compute_kmo(R, R_inv = R_inv)
  names(kmo$KMO_i) <- nms

  # Log-determinant of R: passed to the Bartlett statistic (the log form keeps that
  # statistic finite even when |R| underflows to 0) and used here to report the
  # determinant itself. Shaped like determinant()'s return value, which is what
  # .null_chisq() reads.
  ld <- list(sign = prod(sign(ev)), modulus = sum(log(abs(ev))))
  det_R <- ld$sign * exp(ld$modulus)

  # Bartlett's test of sphericity (Bartlett, 1951): the null-model chi-square with
  # df = p(p - 1)/2. It needs N; without it the test is skipped with a warning and the
  # remaining diagnostics are still returned.
  df <- p * (p - 1) / 2
  bartlett <- if (is.na(N)) {
    cli::cli_warn(
      c("{.arg N} is {.val NA}; Bartlett's test of sphericity was skipped.",
        "i" = "Provide {.arg N} or raw data to include it."),
      class = "efa_suitability_no_n"
    )
    NULL
  } else {
    chisq <- .null_chisq(R, N, ld = ld)
    # Same reason, same condition and same wording as efa_bartlett(), so a user reads the
    # same explanation from either entry point.
    if (is.na(chisq)) .warn_bartlett_n_too_small(N, p)
    list(chisq = chisq,
         p_value = stats::pchisq(chisq, df, lower.tail = FALSE),
         df = df)
  }

  # Condition number of R: the ratio of its largest to its smallest eigenvalue in
  # absolute value, i.e. the exact 2-norm condition number. Large values flag
  # near-collinearity.
  condition <- max(abs(ev)) / min(abs(ev))

  # Squared multiple correlation of each variable with the others (1 - 1/diag(R^-1))
  smc <- .smc_start(R, R_inv = R_inv)
  names(smc) <- nms

  # Raw-data descriptive diagnostics: per-item variance, percentage missing, and a
  # category tabulation with sparse/empty-category flags, gathered into a per-item
  # data frame together with the correlation-side smc and per-variable KMO. These
  # need the raw data, so they are NULL for a correlation-matrix input, where a
  # classed note records why.
  if (is.null(xr)) {
    per_item <- NULL
    categories <- NULL
    normality <- NULL
    outliers <- NULL
    note <- structure(
      list(message = paste(
        "Per-item variance, missing-data, category, normality, and outlier",
        "diagnostics require raw data; only a correlation matrix was supplied.")),
      class = "efa_screen_no_raw"
    )
  } else {
    raw <- .screen_raw_diagnostics(xr, nms, xl)
    # The helper's vectors, smc, and kmo$KMO_i are all in column (nms) order, so the
    # columns are assembled positionally. make.unique() guards against duplicate
    # variable names (which base data.frame() row names would reject) and keys the
    # per-item table and the category list identically, so the two stay aligned.
    item_names <- make.unique(nms)
    per_item <- data.frame(
      variance = raw$variance,
      missing  = raw$missing,
      smc      = smc,
      kmo_i    = kmo$KMO_i,
      flags    = raw$flags,
      row.names = item_names,
      stringsAsFactors = FALSE
    )
    categories <- raw$categories
    names(categories) <- item_names
    # Multivariate-normality tests on the complete cases. These need an invertible
    # complete-case covariance; if the data have too few complete cases (heavy
    # missingness) or collinear variables the tests are skipped so the remaining
    # diagnostics are still returned.
    normality <- tryCatch(
      .screen_normality(xr),
      efa_screen_singular = function(cnd) {
        cli::cli_warn(
          c("Multivariate normality tests were skipped.",
            "i" = paste("The complete-case covariance is singular (too few complete",
                        "cases or collinear variables).")),
          class = "efa_screen_no_mvn"
        )
        structure(
          list(message = paste(
            "Multivariate normality tests require an invertible complete-case",
            "covariance; it was singular (too few complete cases or collinear",
            "variables).")),
          class = "efa_screen_no_mvn"
        )
      }
    )
    # Robust-Mahalanobis outlier diagnostics from the complete cases: a FAST-MCD robust
    # location and scatter, the resulting robust distances, and the flagged rows. The
    # helper degrades gracefully (a classed warning plus classical distances, or a classed
    # note) so the remaining diagnostics are still returned.
    outliers <- .screen_outliers(xr, mcd_alpha = mcd_alpha,
                                 outlier_cutoff = outlier_cutoff, seed = seed)
    note <- NULL
  }

  settings <- list(N = N,
                   n_obs = if (is.null(xr)) NA_integer_ else nrow(xr),
                   use = use,
                   cor_method = cor_method,
                   mcd_alpha = mcd_alpha,
                   outlier_cutoff = outlier_cutoff,
                   seed = seed)

  output <- list(
    kmo = kmo,
    bartlett = bartlett,
    determinant = det_R,
    condition = condition,
    smc = smc,
    # Raw-data screening diagnostics; NULL when a correlation matrix is analysed.
    per_item = per_item,
    normality = normality,
    outliers = outliers,
    categories = categories,
    note = note,
    settings = settings
  )

  class(output) <- "efa_screen"

  output

}

# Re-raise a singular-correlation-matrix condition with the redundancy that caused it.
# `cnd` carries the offending matrix, so the perfectly correlated pairs are read off it
# rather than recomputed. A singular matrix without such a pair is singular through a
# linear combination of several variables (a total score, say), which no pair names.
# The pairs are attached to the re-raised condition (`$pairs`) so a caller can act on
# them without parsing the message. Pairs are named as elsewhere in the package
# (.polychoric()'s pair labels), "first-second".
.screen_singular_abort <- function(cnd, call = rlang::caller_env()) {
  R0 <- cnd$R
  # Only .prepare_cor_input() raises this class here, and it attaches the matrix; keep
  # the shared message rather than failing inside the handler if that ever changes.
  if (!is.matrix(R0)) {
    cli::cli_abort(conditionMessage(cnd), class = "efa_cor_singular", call = call)
  }
  nms <- colnames(R0)
  if (is.null(nms)) nms <- rownames(R0)
  if (is.null(nms)) nms <- paste0("V", seq_len(ncol(R0)))
  hi <- which(abs(R0) > 1 - 1e-8 & upper.tri(R0), arr.ind = TRUE)
  # paste0() recycles its separator against an empty index, so the no-pair case has to
  # be taken before it, not after.
  culprits <- if (nrow(hi) > 0L) {
    paste0(nms[hi[, "row"]], "-", nms[hi[, "col"]])
  } else {
    character()
  }
  cli::cli_abort(
    c("The correlation matrix is singular; no further analyses are performed.",
      if (length(culprits)) {
        c("x" = "{cli::qty(culprits)}Perfectly correlated variable pair{?s}: {.val {culprits}}.",
          "i" = "Drop one variable of each pair and screen again.")
      } else {
        c("x" = "No pair is perfectly correlated, so one or more variables are an exact
                 combination of the others.",
          "i" = "Look for total or difference scores among the variables and drop them.")
      }),
    class = "efa_cor_singular", pairs = culprits, call = call)
}

# Level labels of each column of a raw-data input, in column order: the character
# labels behind the integer codes data.matrix() assigns to a factor or character
# column, and NULL for a column it leaves numeric. A character column is coded via
# factor(), so its labels are its sorted distinct values.
.screen_levels <- function(x) {
  if (is.matrix(x) && is.numeric(x)) return(vector("list", ncol(x)))
  d <- if (is.matrix(x)) as.data.frame(x) else x
  lapply(d, function(col) {
    if (is.factor(col)) levels(col)
    else if (is.character(col)) levels(factor(col))
    else NULL
  })
}

# Per-item raw-data descriptive diagnostics from the coded raw matrix `xr`
# (data.matrix() output, so factor columns are their integer level codes). Every
# statistic is computed column by column on the values as supplied - using each
# column's non-missing values - and named by `nms`, which shares the raw matrix's
# column order. `levels_list` (from .screen_levels()) supplies the labels behind
# those codes so the category counts can be named by the original responses.
# Returns variance, percentage missing, category counts, and the sparse/empty-category
# flags.
.screen_raw_diagnostics <- function(xr, nms, levels_list = NULL) {

  p <- ncol(xr)
  variance <- numeric(p)
  missing <- numeric(p)
  flags <- character(p)
  categories <- vector("list", p)
  if (is.null(levels_list)) levels_list <- vector("list", p)

  for (j in seq_len(p)) {
    col <- xr[, j]

    # Variance over the available values (NA when fewer than two are present) and
    # the percentage missing over all supplied rows.
    variance[j] <- stats::var(col, na.rm = TRUE)
    missing[j] <- 100 * mean(is.na(col))

    lv <- sort(unique(col[!is.na(col)]))
    n_cat <- length(lv)

    # Variables with many distinct values are treated as continuous (the >= 10 rule
    # .polychoric() uses to flag a probably-continuous variable): no tabulation and
    # no sparse/empty screening applies to them.
    if (n_cat >= 10L) {
      categories[[j]] <- NA_integer_
      flags[j] <- NA_character_
      next
    }

    # Marginal category counts in level order (tabulate over the recoded levels, as
    # in .polychoric(); table() would key by value-as-string and drop NA). A column
    # that came in as a factor or character is labelled by its original levels, so the
    # tabulation is readable without re-mapping the codes; otherwise the names use a
    # non-scientific format so a large code labels its category as the plain value
    # rather than, e.g., "1e+08".
    counts <- tabulate(match(col, lv), nbins = n_cat)
    labs <- levels_list[[j]]
    names(counts) <- if (!is.null(labs) && all(lv == floor(lv)) &&
                         min(lv) >= 1L && max(lv) <= length(labs)) {
      labs[lv]
    } else {
      format(lv, scientific = FALSE, trim = TRUE)
    }
    categories[[j]] <- counts

    # Sparse: an observed category with fewer than five responses (the >= 5
    # contingency-table rule; a low-frequency category destabilises polychoric /
    # chi-square estimation). Empty: for integer-coded variables only, an unused
    # integer strictly between the smallest and largest observed value (an ordinal
    # item with a skipped category). With `lv` the sorted distinct integer codes, an
    # interior gap exists exactly when the integers spanned (max - min + 1) outnumber
    # those observed - a constant-time test that handles negative and non-contiguous
    # codes without enumerating the (possibly enormous) span. The span is formed in
    # double precision so an integer-typed code range wider than the 32-bit integer
    # limit does not overflow to NA.
    sparse <- any(counts < 5L)
    empty <- n_cat >= 2L && all(lv == floor(lv)) &&
      (as.double(max(lv)) - min(lv) + 1) > n_cat
    flags[j] <- paste(c("sparse", "empty")[c(sparse, empty)], collapse = ", ")
  }

  names(variance) <- names(missing) <- names(flags) <- names(categories) <- nms
  list(variance = variance, missing = missing, flags = flags,
       categories = categories)
}

# Multivariate-normality diagnostics from the coded raw matrix `xr`. Computes the shared
# linear algebra once - the maximum-likelihood (divisor-n) covariance of the complete
# cases and its inverse - then dispatches to Mardia's and the Henze-Zirkler tests, which
# both reuse the same Mahalanobis geometry. Aborts on a correlation-matrix input or a
# singular complete-case covariance (efa_screen() catches the latter and skips the tests).
.screen_normality <- function(xr) {

  # Guard against a correlation matrix reaching the raw-data tests.
  if (.is_cormat(xr)) {
    cli::cli_abort(
      "Multivariate normality tests need raw data, not a correlation matrix.",
      class = "efa_screen_normality_cormat"
    )
  }

  X <- xr[stats::complete.cases(xr), , drop = FALSE]
  n <- nrow(X)
  p <- ncol(X)

  # Centred data and the biased (divisor-n) covariance S; S^-1 via a Cholesky.
  # Mardia's null moments (p(p + 2), 8p(p + 2)/n) and the Henze-Zirkler statistic are
  # both defined with this maximum-likelihood covariance. S is singular whenever there
  # are p or fewer complete cases: centring costs one degree of freedom, so crossprod(Xc)
  # has rank at most n - 1, and n <= p is checked explicitly because at the n == p
  # boundary a rounding-positive pivot can let chol() succeed on an effectively singular
  # S and return a nonsensical inverse (negative squared Mahalanobis distances). The
  # chol() itself then catches the remaining rank-deficient cases (e.g. collinear
  # variables). efa_screen() catches this abort and skips the tests.
  Xc <- scale(X, scale = FALSE)
  S <- crossprod(Xc) / n
  ch <- if (n <= p) NULL else tryCatch(chol(S), error = function(e) NULL)
  if (is.null(ch)) {
    cli::cli_abort(
      "The complete-case covariance is singular; multivariate normality tests need an invertible covariance.",
      class = "efa_screen_singular"
    )
  }
  Sinv <- chol2inv(ch)

  # Shared Mahalanobis geometry: XSinv %*% t(Xc) gives the Gram matrix D = Xc S^-1 Xc'
  # (reused block by block below), and its diagonal d2 holds the squared Mahalanobis
  # distances. Both tests need the full n x n D, so the accumulation is chunked into
  # row-blocks to cap memory (~ block * n doubles instead of n^2).
  XSinv <- Xc %*% Sinv
  d2 <- rowSums(XSinv * Xc)
  block <- max(1L, min(n, 4194304L %/% n))
  # Both tests form the Gram blocks as XSinv %*% t(Xc), so t(Xc) is transposed once here -
  # alongside d2 - rather than re-transposed inside each helper.
  tXc <- t(Xc)

  list(mardia = .mardia(XSinv, tXc, d2, n, p, block),
       hz = .hz(XSinv, tXc, d2, n, p, block),
       n_complete = n)
}

# Mardia's (1970) multivariate skewness and kurtosis. b1p = (1/n^2) sum_ij d_ij^3 with
# d_ij = x_i' S^-1 x_j (the Gram matrix D), accumulated in row-blocks; b2p = (1/n) sum_i
# d_ii^2. The skewness statistic n b1p / 6 is chi-square with p(p+1)(p+2)/6 df and always
# carries Mardia's (1974) correction; b2p is standardised to N(0, 1) with its exact null
# moments. Ref: Mardia (1970, 1974).
.mardia <- function(XSinv, tXc, d2, n, p, block) {

  s3 <- 0
  for (start in seq.int(1L, n, by = block)) {
    idx <- start:min(start + block - 1L, n)
    D_block <- XSinv[idx, , drop = FALSE] %*% tXc   # |idx| x n block of D
    s3 <- s3 + sum(D_block^3)
  }
  b1p <- s3 / n^2
  b2p <- sum(d2^2) / n

  chi_df <- p * (p + 1) * (p + 2) / 6
  # Mardia's (1974) correction k, eq. (5.5), which he gives "so that E(A') = f for all n"
  # - not only for small samples. From his eq. (5.1) the uncorrected n b1p / 6 has
  # expectation chi_df / k exactly, so without k the statistic is biased low, by a factor
  # that grows with the ratio of variables to observations. The bias is a small fraction
  # of chi_df, but the chi-square scale is sqrt(2 chi_df), so it reaches nearly 3 standard
  # deviations at p = 60 and n = 150 and leaves the test with no power at all there.
  # The correction fixes the mean exactly; the chi-square variance stays wrong by about
  # -20% to +30%, with no published exact variance of b1p to correct it with, so the test
  # keeps a rejection rate of roughly .04 to .08 against a nominal .05.
  # Ref: Mardia (1974), eq. (5.1) and (5.5), and section 6, recommendation (1).
  k <- (p + 1) * (n + 1) * (n + 3) / (n * ((n + 1) * (p + 1) - 6))
  skew <- k * n * b1p / 6

  # b2p is standardised with its EXACT null moments under normality,
  #   E[b2p]   = p(p + 2)(n - 1) / (n + 1)                  Mardia (1970), eq. (3.16)
  #   Var[b2p] = 8p(p + 2)(n - 3)(n - p - 1)(n - p + 1) /
  #              ((n + 1)^2 (n + 3)(n + 5))                 Mardia (1974), eq. (5.3)
  # giving Mardia's (1974) statistic B', eq. (5.6). Both moments become the asymptotic
  # pair p(p + 2) and 8p(p + 2)/n as n grows, which most other implementations use. That
  # pair overstates both, by an amount that grows with the ratio of variables to
  # observations. The overstated mean pushes the statistic down, so the two-sided test
  # rejects data from an exact multivariate normal far too often: more than half the time
  # at p = 40 and n = 200. The overstated variance works the opposite way, so the exact
  # mean with the asymptotic variance - Mardia (1970), eq. (3.20) - is conservative
  # instead, and rejects almost nothing at p = 50 and n = 100. At p = 1 the variance
  # reduces to Fisher's (1930) univariate result, 24n(n - 2)(n - 3) / ((n + 1)^2 (n + 3)
  # (n + 5)). b2p stays mildly right-skewed after standardisation, so the normal
  # approximation splits the two tails unevenly (about 1.7% and 3.0% against a nominal
  # 2.5% at p = 40 and n = 200, and closer with more observations); the two-sided level
  # that the p-value reports holds at .05.
  # The exact variance is 0 at n = p + 1, the smallest sample the caller admits, where
  # the centred data span the whole space orthogonal to the mean, so b2p equals (n - 1)^2
  # with probability 1 and the statistic does not exist.
  kurt_var <- 8 * p * (p + 2) * (n - 3) * (n - p - 1) * (n - p + 1) /
    ((n + 1)^2 * (n + 3) * (n + 5))
  kurt <- if (kurt_var > 0) {
    (b2p - p * (p + 2) * (n - 1) / (n + 1)) / sqrt(kurt_var)
  } else {
    NA_real_
  }

  list(skewness = skew,
       skewness_df = chi_df,
       skewness_p = stats::pchisq(skew, chi_df, lower.tail = FALSE),
       kurtosis = kurt,
       kurtosis_p = 2 * stats::pnorm(abs(kurt), lower.tail = FALSE),
       b1p = b1p, b2p = b2p)
}

# Henze-Zirkler (1990) omnibus test of multivariate normality: a weighted integral of the
# distance between the empirical and normal characteristic functions, evaluated with the
# smoothing parameter b. The three-term statistic reuses the pairwise Mahalanobis
# distances D_jk = d_jj + d_kk - 2 G_jk (G the Gram matrix, accumulated in row-blocks);
# under normality it is approximately lognormal with the parameters (mu, si2) below. That
# approximation degenerates at many variables, where the p-value is withheld (see below).
# Ref: Henze & Zirkler (1990).
.hz <- function(XSinv, tXc, d2, n, p, block) {

  b <- (1 / sqrt(2)) * ((2 * p + 1) * n / 4)^(1 / (p + 4))
  b2 <- b^2

  s1 <- 0
  for (start in seq.int(1L, n, by = block)) {
    idx <- start:min(start + block - 1L, n)
    D_block <- XSinv[idx, , drop = FALSE] %*% tXc
    Djk <- outer(d2[idx], d2, "+") - 2 * D_block   # pairwise Mahalanobis^2
    s1 <- s1 + sum(exp(-b2 / 2 * Djk))
  }
  term1 <- s1 / n^2
  term2 <- 2 * (1 + b2)^(-p / 2) * mean(exp(-(b2 / (2 * (1 + b2))) * d2))
  term3 <- (1 + 2 * b2)^(-p / 2)
  HZ <- n * (term1 - term2 + term3)

  # Lognormal null distribution of HZ (Henze & Zirkler, 1990): mean mu and variance si2,
  # reparameterised to the log scale (pmu, psi); large HZ signals non-normality, so the
  # p-value is the upper tail.
  a <- 1 + 2 * b2
  wb <- (1 + b2) * (1 + 3 * b2)
  mu <- 1 - a^(-p / 2) * (1 + p * b2 / a + p * (p + 2) * b2^2 / (2 * a^2))
  si2 <- 2 * (1 + 4 * b2)^(-p / 2) +
    2 * a^(-p) * (1 + 2 * p * b2^2 / a^2 + 3 * p * (p + 2) * b2^4 / (4 * a^4)) -
    4 * wb^(-p / 2) * (1 + 3 * p * b2^2 / (2 * wb) + p * (p + 2) * b2^4 / (2 * wb^2))
  pmu <- log(sqrt(mu^4 / (si2 + mu^2)))
  psi <- sqrt(log((si2 + mu^2) / mu^2))

  # The lognormal null approximation breaks down at many variables. As p grows si2 falls
  # geometrically while mu goes to 1, so si2 + mu^2 rounds to mu^2, psi rounds to 0, and
  # the z-score (log(HZ) - pmu) / psi becomes -Inf or +Inf: the p-value is then exactly 0
  # or exactly 1, decided by a rounding-level residual, even on data from an exact
  # multivariate normal. The p-value is thus withheld when si2 is no longer resolvable
  # against mu^2, the point at which the reparameterisation loses its input. This is a
  # property of the null moments alone, so where it applies follows n and p and is not a
  # fixed ratio of them. Writing the ratio in cancellation-free form moves the point but
  # does not remove it, because HZ saturates at mu as well. The statistic stays correct
  # and is still reported; only its p-value is unavailable.
  if (!is.finite(psi) || si2 <= mu^2 * .Machine$double.eps) {
    reason <- paste("The lognormal null distribution of the Henze-Zirkler statistic has",
                    "no resolvable spread at", p, "variables, so the statistic cannot be",
                    "turned into a p-value.")
    cli::cli_warn(
      c("The Henze-Zirkler p-value is not available.", "i" = reason),
      class = "efa_screen_no_hz"
    )
    return(structure(
      list(statistic = HZ, p_value = NA_real_, message = reason),
      class = "efa_screen_no_hz"
    ))
  }

  list(statistic = HZ,
       p_value = stats::pnorm((log(HZ) - pmu) / psi, lower.tail = FALSE))
}

# Robust-Mahalanobis outlier diagnostics from the coded raw matrix `xr`. Estimates a
# high-breakdown location and scatter from the robustly rescaled complete cases with
# FAST-MCD, scales it to consistency at the normal model with a small-sample correction,
# reweights it, returns it to the supplied units, and flags observations whose squared
# robust distance exceeds the chi-square cutoff. Degrades gracefully: too few complete
# cases, near-collinear variables, or an exact fit fall back to the classical Mahalanobis
# distance with a classed warning naming which of the three applies, and a fully singular
# covariance yields a classed note. The random subsets are reproducible via `seed`, and the
# caller's RNG state is preserved. Aborts on a correlation-matrix input.
.screen_outliers <- function(xr, mcd_alpha, outlier_cutoff, nsamp = 500L,
                             seed = NULL) {

  # Guard against a correlation matrix reaching the raw-data diagnostics.
  if (.is_cormat(xr)) {
    cli::cli_abort(
      "Outlier diagnostics need raw data, not a correlation matrix.",
      class = "efa_screen_outliers_cormat"
    )
  }

  # A Mahalanobis distance needs a complete row, so the diagnostics use the complete
  # cases; `rows` maps them back to the supplied data for reporting.
  cc <- stats::complete.cases(xr)
  X <- xr[cc, , drop = FALSE]
  rows <- which(cc)
  n <- nrow(X)
  p <- ncol(X)

  # The MCD estimator is affine equivariant -- det(A S A') = det(A)^2 det(S), so ranking
  # covering subsets by the determinant of their covariance, and the Mahalanobis distance
  # built from the winning one, are unchanged by any nonsingular A (Rousseeuw & Van Driessen,
  # 1999). The rcond() gates that decide whether a subset is usable are not: on a covariance
  # they measure how far apart the variables' units are as much as how close the data lie to
  # a lower-dimensional hyperplane, and on four well-conditioned variables a scale ratio of
  # about twenty already tips them. Running the search on columns divided by a robust scale
  # turns them into a test of the subset's *correlation* matrix instead, which is the rank
  # test they were always meant to be. The median absolute deviation is itself
  # high-breakdown, so the scale cannot be set by the very outliers it exists to find; a
  # column with no spread about its median falls back to the standard deviation, and one
  # with no spread at all to 1.
  sc <- apply(X, 2L, stats::mad)
  bad <- !is.finite(sc) | sc <= 0
  if (any(bad)) sc[bad] <- apply(X[, bad, drop = FALSE], 2L, stats::sd)
  bad <- !is.finite(sc) | sc <= 0
  sc[bad] <- 1
  Z <- sweep(X, 2L, sc, "/", check.margin = FALSE)

  # Reproducible random subsets without side effects on the caller's RNG: save the current
  # .Random.seed (or arrange to remove a freshly created one) and restore it on exit, then
  # seed the stream when a seed is supplied.
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
  }
  if (!is.null(seed)) set.seed(seed)

  # Squared-distance critical value used for flagging; the reported cutoff (below) is on
  # the distance scale so that it is directly comparable to the reported `distances`.
  crit <- stats::qchisq(outlier_cutoff, p)

  # Robust location/scatter via FAST-MCD followed by a reweighting step. Too few complete
  # cases (n <= 2p) or a collinear MCD subset fall back to the classical mean/covariance
  # with a classed warning; if even the classical covariance is singular, no distances can
  # be formed and a classed note is returned.
  robust <- tryCatch(
    {
      if (n <= 2L * p) {
        cli::cli_abort("Too few complete cases for a robust covariance.",
                       class = "efa_screen_mcd_unusable")
      }
      h <- .mcd_hsize(n, p, mcd_alpha)
      fit <- .fast_mcd(Z, h, nsamp = nsamp)
      # Consistency (Croux & Haesbroeck, 1999) and small-sample (Pison et al., 2002)
      # scaling of the raw MCD scatter. The consistency factor uses the realised coverage
      # h/n; the small-sample correction uses the nominal `mcd_alpha` (as in its calibration).
      # Both are scalars in (p, n, mcd_alpha) alone, so they commute with the rescaling and
      # can be applied here, before the estimate is put back on the supplied scale.
      raw_center <- fit$center
      raw_cov <- .mcd_consistency(p, h / n) *
        .mcd_cnp2(p, n, mcd_alpha, reweighted = FALSE) * fit$cov
      # Reweighting: keep observations within the 0.975 cutoff of the raw fit, recompute
      # the mean/scatter, and rescale it (consistency at the 0.975 retention plus the
      # reweighted small-sample correction). Fall back to the raw fit if the reweighted
      # scatter would be singular.
      d2_raw <- stats::mahalanobis(Z, raw_center, raw_cov)
      w <- d2_raw <= stats::qchisq(0.975, p)
      rew_cov_sub <- if (sum(w) > p) stats::cov(Z[w, , drop = FALSE]) else NULL
      if (!is.null(rew_cov_sub) && rcond(rew_cov_sub) >= 1e-12) {
        list(center = colMeans(Z[w, , drop = FALSE]),
             cov = .mcd_consistency(p, 0.975) *
               .mcd_cnp2(p, n, mcd_alpha, reweighted = TRUE) * rew_cov_sub,
             method = "mcd")
      } else {
        list(center = raw_center, cov = raw_cov, method = "mcd")
      }
    },
    efa_screen_mcd_unusable = function(cnd) {
      # Fall back to the classical mean/covariance when it is invertible; otherwise signal
      # a total failure so a single, outcome-appropriate condition is raised below.
      # Gate on the reciprocal condition number, matching the reweighting gate above:
      # `mahalanobis()` inverts through `solve()`, which refuses anything whose rcond is
      # below the double-precision epsilon, so the gate has to test the same quantity.
      # Positive definiteness alone is not enough -- a collinear covariance can pass a
      # Cholesky on the sign of a rounding-level pivot (and does on some BLAS
      # implementations) yet still abort in `solve()`. A non-finite covariance is screened
      # off first: it has no condition number, and `rcond()` would abort on it.
      cov_cl <- stats::cov(Z)
      if (n > p && all(is.finite(cov_cl)) && rcond(cov_cl) >= 1e-12) {
        list(center = colMeans(Z), cov = cov_cl, method = "classical")
      } else {
        NULL
      }
    }
  )

  if (is.null(robust)) {
    # Name the condition that actually defeated the classical covariance in the handler
    # above, which gates on all three of n > p, a finite scatter, and its conditioning.
    reason <- if (n <= p) {
      "There are too few complete cases (n <= p) to form a covariance of full rank."
    } else if (!all(is.finite(stats::cov(Z)))) {
      # the same quantity the handler's finiteness gate rejected, recomputed here rather
      # than threaded out of it; this is a cold path, so the second cov() costs nothing
      paste("The complete-case covariance is not finite; the data contain non-finite",
            "or extreme values.")
    } else {
      paste("The complete-case covariance is singular: the variables are exactly or",
            "near-exactly linearly dependent.")
    }
    cli::cli_warn(
      c("Outlier diagnostics were skipped.", "i" = reason),
      class = "efa_screen_no_outliers"
    )
    return(structure(
      list(message = paste("Robust outlier diagnostics require an invertible",
                           "complete-case covariance.", reason),
           reason = reason),
      class = "efa_screen_no_outliers"
    ))
  }

  # A classed warning when the robust estimate was unavailable and classical distances were
  # used instead, naming the reason. Which of the three ways the search can fail applies is
  # read off the data rather than off whichever internal gate happened to fire first, since
  # the same gate is reachable from more than one cause: variables that are near-collinear
  # as a set starve the search of well-conditioned starting subsets, which is the same abort
  # a genuine exact fit produces. Too few complete cases is settled before the search runs;
  # a complete-case covariance that is itself ill-conditioned means the variables are
  # near-collinear; and a well-conditioned one that still admits no usable covering subset
  # means at least h of the cases lie on a lower-dimensional hyperplane. That distinction is
  # the whole of what a reader acts on -- a redundant item to drop, against tied answers
  # that are no fault of the variables -- so it is worth deciding on the evidence. The
  # threshold is the one the robust scatter itself has to clear in .fast_mcd(), on the same
  # rescaled columns.
  fallback_reason <- NULL
  if (identical(robust$method, "classical")) {
    fallback_reason <- if (n <= 2L * p) {
      "There are too few complete cases (n <= 2p) for a robust covariance."
    } else if (rcond(stats::cov(Z)) < 1e-3) {
      paste("The complete-case covariance is ill-conditioned: the variables are so nearly",
            "linearly dependent that no covering subset has a covariance stable enough to",
            "measure distances with. Look for redundant items rather than for outlying",
            "respondents.")
    } else {
      paste("At least half the complete cases lie exactly on a lower-dimensional",
            "hyperplane (an \"exact fit\"). This is common with coarse discrete items,",
            "where many respondents give identical answers on an item pair; it does not",
            "mean the data are collinear at the correlation level.")
    }
    cli::cli_warn(
      c("A robust (MCD) covariance could not be computed; classical Mahalanobis distances were used.",
        "i" = fallback_reason),
      class = "efa_screen_mcd_fallback"
    )
  }

  # Robust distances and outlier flags (squared distance beyond the chi-square cutoff).
  # `distances` and `flagged` share the supplied-data row numbers as their identifier: the
  # distances are named by row number (of the complete cases) and `flagged` holds the row
  # numbers whose distance exceeds the cutoff, so the two align even with incomplete rows.
  # The distances are taken on the rescaled columns, which is where the conditioning gates
  # above did their work: the Mahalanobis distance is the same either way, but putting the
  # scatter back on the supplied scale first can spoil its conditioning by as much as the
  # square of the ratio of the column scales, and `solve()` inside `mahalanobis()` would
  # then refuse a scatter the gates had already passed as usable.
  d2 <- stats::mahalanobis(Z, robust$center, robust$cov)
  distances <- sqrt(d2)
  names(distances) <- as.character(rows)
  flagged <- rows[d2 > crit]

  # The reported centre and scatter go back onto the supplied scale, so that they are in the
  # variables' own units; they describe the distances above, which the rescaling leaves
  # unchanged.
  robust$center <- robust$center * sc
  robust$cov <- outer(sc, sc) * robust$cov

  list(distances = distances,
       cutoff = sqrt(crit),
       flagged = flagged,
       center = robust$center,
       cov = robust$cov,
       method = robust$method,
       fallback_reason = fallback_reason,
       n_complete = n)
}

# Size of the MCD subset covering a proportion `alpha` of the observations (Rousseeuw &
# Van Driessen, 1999): at alpha = 0.5 the maximum-breakdown floor((n + p + 1)/2), rising
# to n at alpha = 1.
.mcd_hsize <- function(n, p, alpha) {
  n_half <- floor((n + p + 1) / 2)
  floor(2 * n_half - n + 2 * (n - n_half) * alpha)
}

# Consistency factor (Croux & Haesbroeck, 1999) that makes the raw MCD scatter consistent
# for the covariance at the multivariate normal model, for a subset covering a proportion
# `alpha`. Equivalent to alpha / F_{chi^2_{p+2}}(qchisq(alpha, p)).
.mcd_consistency <- function(p, alpha) {
  q <- stats::qchisq(alpha, p)
  alpha / stats::pgamma(q / 2, shape = p / 2 + 1)
}

# Small-sample correction factor (Pison, Van Aelst & Willems, 2002) that removes the finite-
# sample bias of the consistency-corrected MCD scatter. `reweighted` selects the factor for
# the raw (FALSE) or the reweighted (TRUE) estimator. For p > 2 the p = 2 and p = 3 fitted
# coefficients are interpolated in log p; for p <= 2 the published closed forms are used.
# The correction is evaluated at the alpha = 0.5 and alpha = 0.875 anchors and interpolated
# linearly in alpha. Ref: Pison et al. (2002), Metrika 55, 111-123.
.mcd_cnp2 <- function(p, n, alpha, reweighted = FALSE) {
  if (reweighted) {
    coef_875 <- matrix(c(-0.544482443573914, 1.25994483222292, 2,
                         -0.343791072183285, 1.25159004257133, 3), ncol = 2)
    coef_500 <- matrix(c(-1.02842572724793, 1.67659883081926, 2,
                         -0.26800273450853, 1.35968562893582, 3), ncol = 2)
    small_2 <- c(3.11101712909049, 1.91401056721863, 0.79473550581058, 1.10081930350091)
    small_1 <- c(1.11098143415027, 1.5182890270453, -0.66046776772861, 0.88939595831888)
  } else {
    coef_875 <- matrix(c(-0.455179464070565, 1.11192541278794, 2,
                         -0.294241208320834, 1.09649329149811, 3), ncol = 2)
    coef_500 <- matrix(c(-1.42764571687802, 1.26263336932151, 2,
                         -1.06141115981725, 1.28907991440387, 3), ncol = 2)
    small_2 <- c(0.673292623522027, 0.691365864961895, 0.446537815635445, 1.06690782995919)
    small_1 <- c(0.262024211897096, 0.604756680630497, -0.351584646688712, 1.01646567502486)
  }
  if (p > 2) {
    y500 <- log(-coef_500[1, ] / p^coef_500[2, ])
    y875 <- log(-coef_875[1, ] / p^coef_875[2, ])
    b500 <- solve(cbind(1, -log(coef_500[3, ] * p^2)), y500)
    b875 <- solve(cbind(1, -log(coef_875[3, ] * p^2)), y875)
    fp_500 <- 1 - exp(b500[1]) / n^b500[2]
    fp_875 <- 1 - exp(b875[1]) / n^b875[2]
  } else if (p == 2) {
    fp_500 <- 1 - exp(small_2[1]) / n^small_2[2]
    fp_875 <- 1 - exp(small_2[3]) / n^small_2[4]
  } else {
    fp_500 <- 1 - exp(small_1[1]) / n^small_1[2]
    fp_875 <- 1 - exp(small_1[3]) / n^small_1[4]
  }
  fp <- if (alpha <= 0.875) {
    fp_500 + (fp_875 - fp_500) / 0.375 * (alpha - 0.5)
  } else {
    fp_875 + (1 - fp_875) / 0.125 * (alpha - 0.875)
  }
  unname(1 / fp)
}

# One concentration step (Rousseeuw & Van Driessen, 1999): given a location/scatter, take
# the `h` observations with the smallest Mahalanobis distances and return their mean and
# covariance; each step cannot increase the covariance determinant. Returns NULL if the
# scatter is singular (so the calling candidate is abandoned).
.mcd_cstep <- function(X, mu, S, h) {
  ch <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(ch)) return(NULL)
  Xc <- sweep(X, 2L, mu, "-", check.margin = FALSE)
  d2 <- rowSums((Xc %*% chol2inv(ch)) * Xc)
  idx <- order(d2)[seq_len(h)]
  Xsub <- X[idx, , drop = FALSE]
  cov_new <- stats::cov(Xsub)
  # Reject a (near-)singular subset: its covariance cannot give stable distances, and the
  # minimum-determinant objective would otherwise collapse onto such a degenerate subset
  # (e.g. a tight cluster on a near-hyperplane) rather than the intended robust scatter.
  if (rcond(cov_new) < 1e-8) return(NULL)
  list(center = colMeans(Xsub), cov = cov_new)
}

# FAST-MCD (Rousseeuw & Van Driessen, 1999): find the h-subset of X whose covariance has the
# smallest determinant. From `nsamp` random (p + 1)-subsets (enlarged until non-singular),
# take two concentration steps each, keep the ten lowest-determinant candidates, and iterate
# those to convergence; return the best subset's mean and covariance. Aborts (classed,
# triggering the classical fallback) when no non-singular starting subset can be formed, when
# the minimum-determinant subset turns out to be singular (an exact fit: at least h
# observations lie on a lower-dimensional hyperplane, so a usable robust covariance does not
# exist), or when the winning subset's covariance is too ill-conditioned to give stable
# distances. Expects columns already on a common scale, so that its rcond() gates test rank
# rather than the variables' units.
.fast_mcd <- function(X, h, nsamp = 500L, maxit = 100L) {
  n <- nrow(X)
  p <- ncol(X)

  logdet <- function(S) as.numeric(determinant(S, logarithm = TRUE)$modulus)

  starts <- vector("list", nsamp)
  dets <- rep_len(NA_real_, nsamp)
  for (i in seq_len(nsamp)) {
    # Draw a (p + 1)-subset with a well-conditioned covariance, enlarging it if needed.
    sub <- sample.int(n, p + 1L)
    S <- stats::cov(X[sub, , drop = FALSE])
    tries <- 0L
    while (rcond(S) < 1e-8 && tries < n) {
      sub <- unique(c(sub, sample.int(n, 1L)))
      S <- stats::cov(X[sub, , drop = FALSE])
      tries <- tries + 1L
    }
    if (rcond(S) < 1e-8) next
    # Two concentration steps. A subset that reaches a degenerate scatter is dropped by
    # .mcd_cstep; only starts that completed at least one step (so `st` is a genuine
    # h-subset, not the raw (p + 1)-point start) are recorded, otherwise the tiny-sample
    # covariance of an unconcentrated start would win the minimum-determinant selection.
    st <- list(center = colMeans(X[sub, , drop = FALSE]), cov = S)
    concentrated <- FALSE
    for (k in 1:2) {
      nxt <- .mcd_cstep(X, st$center, st$cov, h)
      if (is.null(nxt)) break
      st <- nxt
      concentrated <- TRUE
    }
    if (concentrated) {
      starts[[i]] <- st
      dets[i] <- logdet(st$cov)
    }
  }

  ran <- which(!is.na(dets))
  if (!length(ran)) {
    cli::cli_abort("No well-conditioned MCD subset could be formed (exact fit).",
                   class = "efa_screen_mcd_unusable")
  }

  # Iterate the ten lowest-determinant candidates and keep the overall best that reaches a
  # stable fixed point (its determinant stops decreasing, i.e. the subset stabilises).
  keep <- ran[order(dets[ran])][seq_len(min(10L, length(ran)))]
  best <- NULL
  best_ld <- Inf
  for (i in keep) {
    st <- starts[[i]]
    prev <- logdet(st$cov)
    stabilised <- FALSE
    for (it in seq_len(maxit)) {
      nxt <- .mcd_cstep(X, st$center, st$cov, h)
      # A concentration step that lands on a singular h-subset has found a covering subset
      # of determinant zero, and no subset can have a smaller determinant than that: the
      # objective is minimised, and its minimiser is degenerate (an exact fit). Reaching
      # one therefore settles the question whichever candidate got there, and the search
      # stops. Dropping that candidate instead and returning the best non-degenerate fixed
      # point would label an estimate "MCD" that is demonstrably not the minimum-covariance-
      # determinant one, and would make the verdict depend on which random starts happened
      # to be drawn (Rousseeuw & Van Driessen, 1999, sec. 3).
      if (is.null(nxt)) {
        cli::cli_abort("The minimum-determinant subset is singular (exact fit).",
                       class = "efa_screen_mcd_unusable")
      }
      cur <- logdet(nxt$cov)
      st <- nxt
      if (cur >= prev - 1e-12) { prev <- cur; stabilised <- TRUE; break }
      prev <- cur
    }
    if (stabilised && prev < best_ld) {
      best_ld <- prev
      best <- st
    }
  }

  # No candidate settled within `maxit` concentration steps. Each step lowers the
  # determinant and there are finitely many h-subsets, so the iteration must reach a fixed
  # point eventually; exhausting the limit is a numerical stalemate rather than a statement
  # about the data. Signal the classical fallback.
  if (is.null(best)) {
    cli::cli_abort("No MCD subset settled within the concentration-step limit.",
                   class = "efa_screen_mcd_unusable")
  }
  # A near-singular best subset likewise marks a (near-)exact fit, so the robust distances
  # would be dominated by the near-null direction; fall back to the classical distances.
  # Note this is a property of the h-subset, not of the variables as a whole: coarse
  # discrete items routinely put more than h respondents on a hyperplane x_i = x_j while
  # the correlation matrix stays well conditioned, so the reported condition number is not
  # a proxy for it.
  if (rcond(best$cov) < 1e-3) {
    cli::cli_abort("The MCD scatter is near-singular (collinear variables / exact fit).",
                   class = "efa_screen_mcd_unusable")
  }
  best
}
