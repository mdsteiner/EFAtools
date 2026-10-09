#' Screen data for exploratory factor analysis
#'
#' @description
#' Checks whether your data are suitable for exploratory factor analysis (EFA). From a
#' correlation matrix or raw data, it reports the Kaiser-Meyer-Olkin (KMO) measure of
#' sampling adequacy, Bartlett's test of sphericity, the condition index of the correlation
#' matrix, and each variable's squared multiple correlation (SMC). When you supply raw data,
#' it also reports each variable's variance and percentage of missing values, category
#' counts for categorical variables, Mardia's tests of multivariate normality, and
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
#'   polychoric / tetrachoric correlations of ordinal or binary data. Spearman and Kendall
#'   matrices are screened as they are, without conversion to the Pearson scale. Default is
#'   `"pearson"`.
#' @param mcd_alpha numeric. The proportion of cases used to build the robust outlier
#'   estimate, between 0.5 and 1. `0.5` resists the most outliers (up to half the cases).
#'   The default, `0.75`, resists up to a quarter but is more efficient (Hubert et al.,
#'   2018), and Leys et al. (2019) recommend it. Used only with raw data in which all
#'   variables are continuous.
#' @param outlier_cutoff numeric. The probability used to set the cutoff for flagging a
#'   multivariate outlier, between 0.5 and 0.9999. Default is `0.975`. Used only with raw
#'   data in which all variables are continuous.
#' @param seed `r lifecycle::badge("deprecated")` Accepted and ignored. The outlier
#'   detection is deterministic, so it needs no seed.
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
#'     correlate with each other at all; see [efa_bartlett()]. At usual sample sizes the
#'     test is almost always significant, so a significant result is weak evidence of
#'     suitability, but a non-significant result is a strong warning. The test needs the
#'     sample size `N`; without it, `$bartlett` is `NULL`. If `N` is too small relative to
#'     the number of variables, the statistic is `NA`.}
#'   \item{Determinant}{The determinant of \eqn{R}, returned but not used for a verdict. It
#'     falls as you add variables, also when they are not collinear, so a fixed cut-off on
#'     it (such as the 0.00001 often quoted from Field, 2018) mostly reflects the number of
#'     variables.}
#'   \item{Condition index}{The square root of the condition number (the ratio of the
#'     largest to the smallest eigenvalue of \eqn{R}): the collinearity diagnostic of
#'     Belsley, Kuh & Welsch (1980). An index of 10 or less is rarely of interest; an index
#'     above 30 flags a near linear dependency (Belsley, 1991). The report grades such an
#'     index as moderate (30 to 100), strong (100 to 300), or very strong (above 300), and
#'     names the variables involved: those whose variance-decomposition proportions, summed
#'     over all indices above 30, are above .5 (Belsley, Kuh & Welsch, 1980). These values
#'     come from regression on data that are not centred, but a correlation matrix is
#'     centred, so use them as a guide and not as a test.}
#'   \item{SMC}{The squared multiple correlation of each variable with all the others. A
#'     low value flags a variable that has little in common with the rest of your set.}
#'   \item{Variance and missing data}{For raw data: each variable's variance (over its
#'     available values) and percentage of missing values over all rows you supplied.
#'     Factor columns are recoded to integer levels first, so `variance` reflects those
#'     codes. With `use = "complete.obs"` or `"na.or.complete"`, `N` is the number of
#'     complete rows; otherwise it is the number of rows.}
#'   \item{Categories}{For raw data: for each variable with fewer than ten distinct values
#'     (treated as categorical), the count of responses in each category. A category with
#'     fewer than five responses is flagged as sparse, and an unused category between the
#'     smallest and largest response is flagged as empty. As a rule of thumb, items with
#'     fewer than five response categories are better analysed with `cor_method = "poly"`
#'     or `"tetra"` than with Pearson correlations (Rhemtulla et al., 2012).}
#'   \item{Multivariate normality}{For raw data, using only complete cases: Mardia's (1970)
#'     tests of multivariate skewness and kurtosis, with the small-sample corrections of
#'     Mardia (1974). A small p-value suggests your data depart from a multivariate normal
#'     distribution, a reason to prefer robust standard errors over normal-theory ones. In
#'     a very small sample the kurtosis statistic is `NA`. With `cor_method = "poly"` or
#'     `"tetra"` the tests are skipped: they test the observed codes, but a polychoric
#'     analysis assumes normality only of the underlying continuous variables.}
#'   \item{Outliers}{For raw data, using only complete cases: each case's distance from
#'     the centre of your data.
#'     \itemize{
#'       \item If all variables are continuous, the distance is a robust one: it uses the
#'         reweighted minimum covariance determinant (MCD; Rousseeuw & Van Driessen, 1999),
#'         found with the deterministic DetMCD algorithm (Hubert, Rousseeuw & Verdonck,
#'         2012) and corrected for consistency (Croux & Haesbroeck, 1999) and small
#'         samples (Pison et al., 2002). A case is flagged when its squared distance is
#'         above the `outlier_cutoff` quantile of the chi-square distribution. When there
#'         are too few complete cases or the variables are nearly collinear, the classical
#'         (non-robust) distance is used instead, and `$outliers$fallback_reason` tells
#'         why.
#'       \item If any variable is categorical (fewer than ten distinct values, e.g. a
#'         rating-scale item), no case is flagged, and the report names the cases with the
#'         largest classical distances. The MCD theory needs continuous data (Hubert,
#'         Debruyne & Rousseeuw, 2018): with tied answers, many cases can share a
#'         hyperplane, which makes the MCD degenerate. Also, the distances of rating-scale
#'         data do not
#'         follow the chi-square distribution, so a chi-square cutoff flags many clean
#'         cases when items are skewed (Zijlstra, van der Ark & Sijtsma, 2007). To screen
#'         for careless responding, use indices made for that purpose.
#'     }
#'     Distance-based screening does not detect straight-lining (the same answer to all
#'     items).}
#' }
#'
#' @returns An object of class `efa_screen`, a list containing:
#' \item{kmo}{A list with the overall KMO (`KMO`) and the per-variable KMO (`KMO_i`).}
#' \item{bartlett}{A list with Bartlett's chi-square statistic (`chisq`), its `p_value`,
#'   and its degrees of freedom (`df`); `chisq` and `p_value` are `NA` when `N` was too
#'   small for the correction. `NULL` when `N` is unavailable.}
#' \item{determinant}{The determinant of the correlation matrix.}
#' \item{condition}{The condition number of the correlation matrix (largest eigenvalue
#'   over smallest). Its square root is the condition index.}
#' \item{collinear}{The names of the variables involved in a near linear dependency (see
#'   Details). `NULL` when no condition index is above 30.}
#' \item{smc}{The per-variable squared multiple correlations.}
#' \item{per_item}{A data frame with one row per variable (row names are the variable
#'   names): `variance`, `missing` (percentage), `smc`, `kmo_i`, and `flags` (any
#'   sparse/empty-category issues). `NULL` when a correlation matrix is supplied instead
#'   of raw data.}
#' \item{normality}{A list with `mardia` (skewness `skewness`, `skewness_df`,
#'   `skewness_p`, kurtosis `kurtosis` and `kurtosis_p`, and the underlying `b1p`/`b2p`)
#'   and `n_complete` (the number of complete cases used). `NULL` without raw data, or a
#'   note (`$message`) that tells why the tests were skipped.}
#' \item{outliers}{A list with `distances` (each complete case's distance, named by its
#'   row number), `cutoff` (the flagging threshold, on the same scale as `distances`),
#'   `flagged` (the row numbers above `cutoff`), `center` and `cov` (the location and
#'   scatter used), `method`, `fallback_reason`, and `n_complete`. `method` is `"mcd"`,
#'   `"classical"` (the fallback for continuous data; `fallback_reason` tells why), or
#'   `"ranking"` (categorical data: classical distances, `cutoff` is `NA`, and `flagged`
#'   is empty). `NULL` without raw data, or a note (`$message`) that tells why no
#'   distances could be computed.}
#' \item{categories}{A named list with the response-category counts for each categorical
#'   variable (in category order); `NA` for a variable treated as continuous. `NULL`
#'   without raw data.}
#' \item{note}{A note (`$message`) that tells why the raw-data diagnostics (`per_item`,
#'   `normality`, `outliers`, `categories`) are missing, when a correlation matrix is
#'   supplied instead of raw data. `NULL` when raw data are supplied.}
#' \item{settings}{The settings used: `N`, `n_obs` (rows in the raw data supplied, `NA`
#'   for a correlation-matrix input), `use`, `cor_method`, `mcd_alpha`, and
#'   `outlier_cutoff`.}
#'
#' @source Bartlett, M. S. (1951). The effect of standardization on a Chi-square
#'   approximation in factor analysis. Biometrika, 38, 337-344.
#' @source Belsley, D. A. (1991). A guide to using the collinearity diagnostics. Computer
#'   Science in Economics and Management, 4, 33-50.
#' @source Belsley, D. A., Kuh, E. & Welsch, R. E. (1980). Regression diagnostics:
#'   Identifying influential data and sources of collinearity. Wiley.
#' @source Croux, C. & Haesbroeck, G. (1999). Influence function and efficiency of the
#'   minimum covariance determinant scatter matrix estimator. Journal of Multivariate
#'   Analysis, 71, 161-190.
#' @source Field, A. (2018). Discovering statistics using IBM SPSS statistics (5th
#'   ed.). Sage.
#' @source Hubert, M., Debruyne, M. & Rousseeuw, P. J. (2018). Minimum covariance
#'   determinant and extensions. WIREs Computational Statistics, 10, e1421.
#' @source Hubert, M., Rousseeuw, P. J. & Verdonck, T. (2012). A deterministic algorithm
#'   for robust location and scatter. Journal of Computational and Graphical Statistics,
#'   21, 618-637.
#' @source Kaiser, H. F. (1970). A second generation little jiffy. Psychometrika,
#'   35, 401-415.
#' @source Kaiser, H. F. & Rice, J. (1974). Little jiffy, mark IV. Educational and
#'   Psychological Measurement, 34, 111-117.
#' @source Leys, C., Delacre, M., Mora, Y. L., Lakens, D. & Ley, C. (2019). How to
#'   classify, detect, and manage univariate and multivariate outliers, with emphasis on
#'   pre-registration. International Review of Social Psychology, 32, 5.
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
#' @source Zijlstra, W. P., van der Ark, L. A. & Sijtsma, K. (2007). Outlier detection in
#'   test and questionnaire data. Multivariate Behavioral Research, 42, 531-555.
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
#' # From raw data (N is taken from the data)
#' efa_screen(GRiPS_raw)
#'
efa_screen <- function(x, N = NA,
                       use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                               "everything", "na.or.complete"),
                       cor_method = c("pearson", "spearman", "kendall", "poly",
                                      "tetra"),
                       mcd_alpha = 0.75, outlier_cutoff = 0.975,
                       seed = lifecycle::deprecated()) {

  # Perform argument checks
  .assert_cor_input(x)

  # Below three variables the KMO is undefined (p = 1) or exactly .5 (p = 2), and
  # Bartlett's test has at most one degree of freedom.
  if (ncol(x) < 3L) {
    cli::cli_abort(
      c("Screening data for factor analysis needs at least three variables.",
        "x" = "{.arg x} has {ncol(x)} variable{?s}.",
        "i" = "The sampling adequacy of one or two variables is not defined: the KMO is
               undefined for a single variable and exactly {.val {0.5}} for any pair."),
      class = "efa_screen_too_few_vars")
  }

  if (lifecycle::is_present(seed)) {
    lifecycle::deprecate_warn(
      when = "1.2.0",
      what = "efa_screen(seed)",
      details = c("i" = "The outlier detection is deterministic, so it needs no seed.")
    )
  }

  use <- .match_arg_ci(use)
  cor_method <- .match_arg_ci(cor_method)
  checkmate::assert_count(N, na.ok = TRUE)
  checkmate::assert_number(mcd_alpha, lower = 0.5, upper = 1)
  # A cutoff probability of 0 flags every case and 1 puts the threshold at infinity.
  checkmate::assert_number(outlier_cutoff, lower = 0.5, upper = 0.9999)

  # The raw data coded as .polychoric() reads them: data.matrix() turns factor and character
  # columns into integer level codes (as.matrix() would turn every column into character).
  # A non-numeric matrix goes through as.data.frame() to get the same coding. NULL for a
  # correlation matrix.
  xr <- if (.is_cormat(x)) {
    NULL
  } else if (is.matrix(x) && !is.numeric(x)) {
    data.matrix(as.data.frame(x))
  } else {
    data.matrix(x)
  }
  # The level labels behind that coding, to name the category counts.
  xl <- if (is.null(xr)) NULL else .screen_levels(x)

  # N stays NA for a correlation matrix without N; only Bartlett's test needs it. A singular
  # matrix is refused, because every measure here needs R^-1, and the abort names the
  # redundant variables. The handler runs in tryCatch()'s frame, so it gets this call.
  screen_call <- rlang::current_env()
  prep <- tryCatch(
    .prepare_cor_input(x, N = N, use = use, cor_method = cor_method,
                       N_policy = "optional", inform_from_data = FALSE),
    efa_cor_singular = function(cnd) .screen_singular_abort(cnd, call = screen_call))
  R <- prep$R
  N <- prep$N
  p <- ncol(R)

  nms <- colnames(R)
  if (is.null(nms)) nms <- rownames(R)
  if (is.null(nms)) nms <- paste0("V", seq_len(p))

  # One eigendecomposition and one inverse serve all correlation-side measures. The absolute
  # eigenvalues of a symmetric matrix are its singular values, so their ratio is the 2-norm
  # condition number, and prod(sign) with sum(log|lambda|) is the log-determinant.
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  # If the Cholesky factor fails, each measure inverts R itself, which keeps
  # .compute_kmo()'s classed singularity abort.
  R_inv <- tryCatch(chol2inv(chol(R)), error = function(e) NULL)

  kmo <- .compute_kmo(R, R_inv = R_inv)
  names(kmo$KMO_i) <- nms

  # In determinant()'s shape, which .null_chisq() reads; the log form keeps Bartlett's
  # statistic finite when |R| underflows to 0.
  ld <- list(sign = prod(sign(ev)), modulus = sum(log(abs(ev))))
  det_R <- ld$sign * exp(ld$modulus)

  # Bartlett's test of sphericity (Bartlett, 1951) with df = p(p - 1)/2; NULL without N.
  df <- p * (p - 1) / 2
  bartlett <- if (is.na(N)) {
    NULL
  } else {
    chisq <- .null_chisq(R, N, ld = ld)
    list(chisq = chisq,
         p_value = stats::pchisq(chisq, df, lower.tail = FALSE),
         df = df)
  }

  condition <- max(abs(ev)) / min(abs(ev))
  collinear <- .screen_collinear(R, condition, make.unique(nms))

  # Squared multiple correlations, 1 - 1/diag(R^-1)
  smc <- .smc_start(R, R_inv = R_inv)
  names(smc) <- nms

  if (is.null(xr)) {
    per_item <- NULL
    categories <- NULL
    normality <- NULL
    outliers <- NULL
    note <- .screen_note("efa_screen_no_raw", paste(
      "Per-item variance, missing-data, category, normality, and outlier diagnostics",
      "require raw data; only a correlation matrix was supplied."))
  } else {
    raw <- .screen_raw_diagnostics(xr, nms, xl)
    # All columns are in the order of `nms`. make.unique() keeps duplicate variable names
    # valid as row names, and keys the per-item table and the category list identically.
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
    # Mardia's tests test the observed codes, but a polychoric analysis assumes normality
    # only of the underlying continuous variables, so the tests do not apply there.
    normality <- if (.is_poly_cor(cor_method)) {
      .screen_note("efa_screen_no_mvn", paste(
        "Multivariate normality tests were not computed: they do not apply to a",
        "polychoric or tetrachoric analysis."))
    } else {
      .screen_normality(xr) %||% .screen_note("efa_screen_no_mvn", paste(
        "Multivariate normality tests were not computed: the complete-case covariance is",
        "singular (too few complete cases or collinear variables)."))
    }
    # Categorical: fewer than ten distinct values, the rule of the category screening.
    categorical <- any(!is.na(raw$flags))
    outliers <- .screen_outliers(xr, categorical = categorical, mcd_alpha = mcd_alpha,
                                 outlier_cutoff = outlier_cutoff)
    note <- NULL
  }

  settings <- list(N = N,
                   n_obs = if (is.null(xr)) NA_integer_ else nrow(xr),
                   use = use,
                   cor_method = cor_method,
                   mcd_alpha = mcd_alpha,
                   outlier_cutoff = outlier_cutoff)

  output <- list(
    kmo = kmo,
    bartlett = bartlett,
    determinant = det_R,
    condition = condition,
    collinear = collinear,
    smc = smc,
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

# A classed note that tells why a diagnostic is missing.
.screen_note <- function(class, message) structure(list(message = message), class = class)

# Variables involved in the near linear dependencies that a condition index above 30 flags:
# those whose variance-decomposition proportions, summed over all such indices, are above
# .5. The sum handles competing dependencies (Belsley, Kuh & Welsch, 1980, ch. 3). NULL when
# no index is above 30.
.screen_collinear <- function(R, condition, nms) {
  if (!isTRUE(sqrt(condition) > 30)) return(NULL)
  e <- eigen(R, symmetric = TRUE)
  lam <- abs(e$values)
  phi <- sweep(e$vectors^2, 2L, lam, "/")        # phi[j, k] = v_jk^2 / lambda_k
  vdp <- phi / rowSums(phi)
  high <- sqrt(max(lam) / lam) > 30
  nms[rowSums(vdp[, high, drop = FALSE]) > .5]
}

# Re-raise a singular-correlation-matrix condition with the perfectly correlated pairs
# ("first-second", as .polychoric() labels pairs), read off the matrix the condition
# carries and attached as `$pairs`. Without such a pair, one or more variables are an
# exact linear combination of the others.
.screen_singular_abort <- function(cnd, call = rlang::caller_env()) {
  R0 <- cnd$R
  if (!is.matrix(R0)) {
    cli::cli_abort(conditionMessage(cnd), class = "efa_cor_singular", call = call)
  }
  nms <- colnames(R0)
  if (is.null(nms)) nms <- rownames(R0)
  if (is.null(nms)) nms <- paste0("V", seq_len(ncol(R0)))
  hi <- which(abs(R0) > 1 - 1e-8 & upper.tri(R0), arr.ind = TRUE)
  # paste0() recycles its separator against an empty index, so test for no pair first.
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

# Level labels of each raw-data column, in column order: the labels behind the integer codes
# that data.matrix() gives a factor or character column (for a character column, its sorted
# distinct values), and NULL for a numeric column.
.screen_levels <- function(x) {
  if (is.matrix(x) && is.numeric(x)) return(vector("list", ncol(x)))
  d <- if (is.matrix(x)) as.data.frame(x) else x
  lapply(d, function(col) {
    if (is.factor(col)) levels(col)
    else if (is.character(col)) levels(factor(col))
    else NULL
  })
}

# Per-item descriptive diagnostics from the coded raw matrix `xr`, column by column over each
# column's non-missing values: variance, percentage missing, category counts (named by the
# labels in `levels_list`, from .screen_levels()), and the sparse/empty-category flags.
.screen_raw_diagnostics <- function(xr, nms, levels_list = NULL) {

  p <- ncol(xr)
  variance <- numeric(p)
  missing <- numeric(p)
  flags <- character(p)
  categories <- vector("list", p)
  if (is.null(levels_list)) levels_list <- vector("list", p)

  for (j in seq_len(p)) {
    col <- xr[, j]

    variance[j] <- stats::var(col, na.rm = TRUE)
    missing[j] <- 100 * mean(is.na(col))

    lv <- sort(unique(col[!is.na(col)]))
    n_cat <- length(lv)

    # Ten or more distinct values: treated as continuous, as .polychoric() does.
    if (n_cat >= 10L) {
      categories[[j]] <- NA_integer_
      flags[j] <- NA_character_
      next
    }

    # Counts in level order. A factor or character column is labelled by its own levels;
    # other codes are labelled as plain values (not, e.g., "1e+08").
    counts <- tabulate(match(col, lv), nbins = n_cat)
    labs <- levels_list[[j]]
    names(counts) <- if (!is.null(labs) && all(lv == floor(lv)) &&
                         min(lv) >= 1L && max(lv) <= length(labs)) {
      labs[lv]
    } else {
      format(lv, scientific = FALSE, trim = TRUE)
    }
    categories[[j]] <- counts

    # Sparse: an observed category with fewer than five responses. Empty: for integer codes,
    # an unused integer between the smallest and largest code, i.e. the span max - min + 1
    # (in double precision, so a wide integer span cannot overflow) is larger than the
    # number of observed codes.
    sparse <- any(counts < 5L)
    empty <- n_cat >= 2L && all(lv == floor(lv)) &&
      (as.double(max(lv)) - min(lv) + 1) > n_cat
    flags[j] <- paste(c("sparse", "empty")[c(sparse, empty)], collapse = ", ")
  }

  names(variance) <- names(missing) <- names(flags) <- names(categories) <- nms
  list(variance = variance, missing = missing, flags = flags,
       categories = categories)
}

# Mardia's tests from the complete cases of the coded raw matrix `xr`; NULL when the
# complete-case covariance is singular. The tests use the maximum-likelihood (divisor-n)
# covariance S = U'U (Cholesky). The whitened data Y = Xc U^-1 give x_i' S^-1 x_j = y_i'y_j.
.screen_normality <- function(xr) {
  X <- xr[stats::complete.cases(xr), , drop = FALSE]
  n <- nrow(X)
  p <- ncol(X)
  Xc <- scale(X, scale = FALSE)
  # S is singular with n <= p complete cases (centring costs one degree of freedom). The
  # explicit test stops a rounding-positive pivot at n = p from passing chol(), which then
  # catches the remaining rank-deficient cases.
  ch <- if (n <= p) NULL else tryCatch(chol(crossprod(Xc) / n), error = function(e) NULL)
  if (is.null(ch)) return(NULL)
  Y <- Xc %*% backsolve(ch, diag(p))
  list(mardia = .mardia(Y), n_complete = n)
}

# Mardia's (1970) multivariate skewness b1p = n^-2 sum_ij (y_i'y_j)^3 and kurtosis
# b2p = n^-1 sum_i (y_i'y_i)^2 from the whitened complete cases Y, with Mardia's (1974)
# small-sample corrections.
.mardia <- function(Y) {

  n <- nrow(Y)
  p <- ncol(Y)

  # b1p costs O(n^2 p) through the Gram matrix Y Y', or O(n p^3) through the identity
  # sum_ij (y_i'y_j)^3 = sum_abc (sum_i y_ia y_ib y_ic)^2 (the squared norm of the
  # third-moment tensor). The tensor form is cheaper when n > p^2; its symmetry in a and b
  # halves the work. The Gram form runs in row blocks, which keeps memory at about 4e6
  # doubles.
  s3 <- 0
  if (n > p^2) {
    for (a in seq_len(p)) {
      f <- colSums(crossprod(Y * Y[, a], Y[, a:p, drop = FALSE])^2)   # b = a, ..., p
      s3 <- s3 + f[1L] + 2 * sum(f[-1L])
    }
  } else {
    tY <- t(Y)
    block <- max(1L, 4194304L %/% n)
    for (start in seq.int(1L, n, by = block)) {
      idx <- start:min(start + block - 1L, n)
      s3 <- s3 + sum((Y[idx, , drop = FALSE] %*% tY)^3)
    }
  }
  b1p <- s3 / n^2
  b2p <- sum(rowSums(Y^2)^2) / n

  # The skewness statistic k n b1p / 6 is chi-square with p(p + 1)(p + 2)/6 df. Mardia's
  # (1974) correction k, eq. (5.5), makes its mean equal the df for all n (eq. 5.1); without
  # it the statistic is biased low when p/n is large. Its variance can still be about 20%
  # too low or 30% too high, so the rejection rate is about .04 to .08 at a nominal .05.
  chi_df <- p * (p + 1) * (p + 2) / 6
  k <- (p + 1) * (n + 1) * (n + 3) / (n * ((n + 1) * (p + 1) - 6))
  skew <- k * n * b1p / 6

  # b2p standardised with its exact null moments, Mardia's (1974) statistic B', eq. (5.6):
  #   E[b2p]   = p(p + 2)(n - 1) / (n + 1)                  Mardia (1970), eq. (3.16)
  #   Var[b2p] = 8p(p + 2)(n - 3)(n - p - 1)(n - p + 1) /
  #              ((n + 1)^2 (n + 3)(n + 5))                 Mardia (1974), eq. (5.3)
  # The usual asymptotic pair p(p + 2) and 8p(p + 2)/n overstates both moments when p/n is
  # large, and the test then rejects normal data far too often. The variance is 0 at
  # n = p + 1, where b2p = (n - 1)^2 exactly, so the statistic is NA there.
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

# High-breakdown scale of a vector: the MAD, or the standard deviation when the MAD is 0,
# or 1 when the vector is constant.
.robust_scale <- function(v) {
  s <- stats::mad(v)
  if (!is.finite(s) || s <= 0) s <- stats::sd(v)
  if (!is.finite(s) || s <= 0) s <- 1
  s
}

# Outlier diagnostics from the complete cases of the coded raw matrix `xr`. For continuous
# data: the reweighted MCD, found with DetMCD and corrected for consistency and small samples,
# and the cases whose squared distance is above the chi-square cutoff; classical distances
# (with `fallback_reason`) when the MCD is unusable. For `categorical` data: classical
# distances without a cutoff (method "ranking"), because the MCD needs data in general
# position and a chi-square cutoff over-flags skewed rating-scale data. A classed note when no
# covariance can be inverted.
.screen_outliers <- function(xr, mcd_alpha, outlier_cutoff, categorical = FALSE) {

  # `rows` maps the complete cases back to the rows supplied.
  cc <- stats::complete.cases(xr)
  X <- xr[cc, , drop = FALSE]
  rows <- which(cc)
  n <- nrow(X)
  p <- ncol(X)

  # DetMCD works on columns standardised by their median and a robust scale (Hubert et al.,
  # 2012). The standardisation also makes the rcond() gates below test the rank of the
  # data rather than the ratio of the variables' units.
  med <- apply(X, 2L, stats::median)
  sc <- apply(X, 2L, .robust_scale)
  Z <- sweep(sweep(X, 2L, med, check.margin = FALSE), 2L, sc, "/", check.margin = FALSE)

  crit <- stats::qchisq(outlier_cutoff, p)

  # The classical mean and covariance, when mahalanobis() can invert them: its solve()
  # refuses a matrix with a small rcond(), and rcond() needs finite values.
  classical <- function(method) {
    cov_cl <- stats::cov(Z)
    if (n > p && all(is.finite(cov_cl)) && rcond(cov_cl) >= 1e-12) {
      list(center = colMeans(Z), cov = cov_cl, method = method)
    } else {
      NULL
    }
  }

  robust <- if (categorical) classical("ranking") else tryCatch(
    {
      if (n <= 2L * p) {
        cli::cli_abort("Too few complete cases for a robust covariance.",
                       class = "efa_screen_mcd_unusable")
      }
      h <- .mcd_hsize(n, p, mcd_alpha)
      fit <- .det_mcd(Z, h)
      # Consistency factor (Croux & Haesbroeck, 1999) at the realised coverage h/n and
      # small-sample factor (Pison et al., 2002) at the nominal alpha, as calibrated.
      raw_cov <- .mcd_consistency(p, h / n) *
        .mcd_cnp2(p, n, mcd_alpha, reweighted = FALSE) * fit$cov
      # Reweighting: the mean and scatter of the cases within the 0.975 cutoff of the raw
      # fit, with the matching factors. The raw fit stays when that scatter is singular.
      w <- stats::mahalanobis(Z, fit$center, raw_cov) <= stats::qchisq(0.975, p)
      rew_cov <- if (sum(w) > p) stats::cov(Z[w, , drop = FALSE]) else NULL
      if (!is.null(rew_cov) && rcond(rew_cov) >= 1e-12) {
        list(center = colMeans(Z[w, , drop = FALSE]),
             cov = .mcd_consistency(p, 0.975) *
               .mcd_cnp2(p, n, mcd_alpha, reweighted = TRUE) * rew_cov,
             method = "mcd")
      } else {
        list(center = fit$center, cov = raw_cov, method = "mcd")
      }
    },
    efa_screen_mcd_unusable = function(cnd) classical("classical")
  )

  if (is.null(robust)) {
    reason <- if (n <= p) {
      "there are too few complete cases (n <= p)"
    } else if (!all(is.finite(stats::cov(Z)))) {
      "the complete-case covariance is not finite (look for extreme values)"
    } else {
      "the variables are exactly or nearly linearly dependent"
    }
    return(.screen_note("efa_screen_no_outliers",
                        paste0("Outlier diagnostics were skipped: ", reason, ".")))
  }

  # The cause of a fallback is read off the data, not off the gate that fired: nearly
  # collinear variables and an exact fit both leave no usable h-subset, but the remedies
  # differ (drop a redundant item, against nothing to fix). The rcond() threshold is the one
  # the MCD scatter must pass in .det_mcd().
  fallback_reason <- NULL
  if (identical(robust$method, "classical")) {
    fallback_reason <- if (n <= 2L * p) {
      "there are too few complete cases (n <= 2p)"
    } else if (rcond(stats::cov(Z)) < 1e-3) {
      "the variables are nearly linearly dependent; look for redundant items"
    } else {
      "many cases lie on a hyperplane (an exact fit)"
    }
  }

  # The distances use the standardised scale, where the gates above judged the scatter. The
  # centre and scatter are returned in the variables' own units.
  d2 <- stats::mahalanobis(Z, robust$center, robust$cov)
  distances <- sqrt(d2)
  names(distances) <- as.character(rows)
  ranking <- identical(robust$method, "ranking")

  list(distances = distances,
       cutoff = if (ranking) NA_real_ else sqrt(crit),
       flagged = if (ranking) integer(0) else rows[d2 > crit],
       center = robust$center * sc + med,
       cov = outer(sc, sc) * robust$cov,
       method = robust$method,
       fallback_reason = fallback_reason,
       n_complete = n)
}

# Size of the MCD subset that covers a proportion `alpha` of the observations (Rousseeuw &
# Van Driessen, 1999): floor((n + p + 1)/2) at alpha = 0.5, rising to n at alpha = 1.
.mcd_hsize <- function(n, p, alpha) {
  n_half <- floor((n + p + 1) / 2)
  floor(2 * n_half - n + 2 * (n - n_half) * alpha)
}

# Consistency factor (Croux & Haesbroeck, 1999) that makes the raw MCD scatter consistent at
# the multivariate normal model for a subset that covers a proportion `alpha`:
# alpha / F_{chi^2_{p+2}}(qchisq(alpha, p)).
.mcd_consistency <- function(p, alpha) {
  q <- stats::qchisq(alpha, p)
  alpha / stats::pgamma(q / 2, shape = p / 2 + 1)
}

# Small-sample correction factor (Pison, Van Aelst & Willems, 2002) for the raw
# (`reweighted = FALSE`) or reweighted MCD scatter. For p > 2 the fitted p = 2 and p = 3
# coefficients are interpolated in log p; for p <= 2 the published closed forms apply. The
# factor is evaluated at alpha = 0.5 and 0.875 and interpolated linearly in alpha.
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

# One concentration step (Rousseeuw & Van Driessen, 1999): the mean and covariance of the `h`
# observations closest to (mu, S). A step cannot increase the covariance determinant.
# Returns NULL when S or the new covariance is (near-)singular.
.mcd_cstep <- function(X, mu, S, h) {
  ch <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(ch)) return(NULL)
  Xc <- sweep(X, 2L, mu, "-", check.margin = FALSE)
  d2 <- rowSums((Xc %*% chol2inv(ch)) * Xc)
  idx <- order(d2)[seq_len(h)]
  Xsub <- X[idx, , drop = FALSE]
  cov_new <- stats::cov(Xsub)
  if (rcond(cov_new) < 1e-8) return(NULL)
  list(center = colMeans(Xsub), cov = cov_new)
}

# DetMCD (Hubert, Rousseeuw & Verdonck, 2012): a deterministic search for the h-subset of the
# standardised data Z (columns centred on their median and divided by a robust scale) with
# the smallest covariance determinant. Each of six robust initial estimates gives an
# h-subset, which C-steps concentrate to a fixed point; the fixed point with the smallest
# determinant is returned. Aborts with class efa_screen_mcd_unusable on an exact fit (a
# singular h-subset has the smallest possible determinant, so the MCD is degenerate), when no
# start reaches a fixed point, or when the best scatter is near-singular.
.det_mcd <- function(Z, h, maxit = 100L) {
  n <- nrow(Z)
  p <- ncol(Z)
  unusable <- function(msg) cli::cli_abort(msg, class = "efa_screen_mcd_unusable")
  logdet <- function(S) as.numeric(determinant(S, logarithm = TRUE)$modulus)

  # The six initial scatter estimates (Hubert et al., 2012, sec. 2.2): the correlations of
  # tanh(Z), of the ranks (Spearman), and of the normal scores; the spatial-sign covariance;
  # the covariance of the half of the cases closest to the centre; and the pairwise
  # Gnanadesikan-Kettenring estimate that the OGK estimator orthogonalises.
  zn <- sqrt(rowSums(Z^2))
  nz <- zn > .Machine$double.eps
  gk <- diag(p)
  for (j in seq_len(p)[-1L]) {
    for (k in seq_len(j - 1L)) {
      gk[j, k] <- gk[k, j] <- (.robust_scale(Z[, j] + Z[, k])^2 -
                                 .robust_scale(Z[, j] - Z[, k])^2) / 4
    }
  }
  starts <- list(
    stats::cor(tanh(Z)),
    stats::cor(Z, method = "spearman"),
    stats::cor(stats::qnorm((apply(Z, 2L, rank) - 1 / 3) / (n + 1 / 3))),
    crossprod(Z[nz, , drop = FALSE] / zn[nz]),
    stats::cov(Z[order(zn)[seq_len(ceiling(n / 2))], , drop = FALSE]),
    gk
  )

  best <- NULL
  best_ld <- Inf
  for (S0 in starts) {
    if (anyNA(S0)) next
    # The start's eigenvectors P with robust scales lam of the projections B = Z P give the
    # scatter P diag(lam^2) P'. Its location is the coordinatewise median of the data
    # whitened by that scatter (B diag(1/lam) P'), mapped back. The h cases closest to it
    # form the start; (Z - loc) P = B - loc P.
    P <- eigen(S0, symmetric = TRUE)$vectors
    B <- Z %*% P
    lam <- apply(B, 2L, .robust_scale)
    loc <- drop(apply(sweep(B, 2L, lam, "/", check.margin = FALSE) %*% t(P), 2L,
                      stats::median) %*% (P %*% (lam * t(P))))
    d <- rowSums(sweep(sweep(B, 2L, drop(loc %*% P), check.margin = FALSE), 2L, lam, "/",
                       check.margin = FALSE)^2)
    idx <- order(d)[seq_len(h)]
    st <- list(center = colMeans(Z[idx, , drop = FALSE]),
               cov = stats::cov(Z[idx, , drop = FALSE]))
    if (rcond(st$cov) < 1e-8) unusable("An initial h-subset is singular (exact fit).")

    prev <- logdet(st$cov)
    stabilised <- FALSE
    for (it in seq_len(maxit)) {
      nxt <- .mcd_cstep(Z, st$center, st$cov, h)
      if (is.null(nxt)) unusable("A concentration step reached a singular subset (exact fit).")
      cur <- logdet(nxt$cov)
      st <- nxt
      if (cur >= prev - 1e-12) {
        stabilised <- TRUE
        break
      }
      prev <- cur
    }
    if (stabilised && cur < best_ld) {
      best_ld <- cur
      best <- st
    }
  }

  if (is.null(best)) unusable("No start reached a concentration fixed point.")
  # A near-singular best scatter marks a near-exact fit of the h-subset, which a
  # well-conditioned correlation matrix does not rule out (coarse items with tied answers).
  if (rcond(best$cov) < 1e-3) unusable("The MCD scatter is near-singular (exact fit).")
  best
}
