#' Hull method for determining the number of factors to retain
#'
#' Implementation of the Hull method suggested by Lorenzo-Seva, Timmerman,
#' and Kiers (2011), with an extension to principal axis factoring. See details for
#' parallelization.
#'
#' @param x matrix or data.frame. Dataframe or matrix of raw data or matrix with
#' correlations.
#' @param N numeric. Number of cases in the data. This is passed to [PARALLEL].
#'  Only has to be specified if x is a correlation matrix, otherwise it is determined
#'  based on the dimensions of x.
#' @param n_fac_theor numeric. Theoretical number of factors to retain. The maximum
#'   of this number and the number of factors suggested by [PARALLEL] plus
#'   one will be used in the Hull method.
#' @param method character. The estimation method to use. One of  `"PAF"`,
#'    `"ULS"`, or  `"ML"`, for principal axis factoring, unweighted
#'    least squares, and maximum likelihood, respectively.
#' @param gof character. The goodness of fit index to use. Either `"CAF"`,
#'   `"CFI"`, or `"RMSEA"`, or any combination of them.
#'   If `method = "PAF"` is used, only
#'   the CAF can be used as goodness of fit index. For details on the CAF, see
#'   Lorenzo-Seva, Timmerman, and Kiers (2011).
#' @param eigen_type character. On what the eigenvalues should be found in the
#'  parallel analysis. Can be one of `"SMC"`, `"PCA"`, or `"EFA"`.
#'   If using  `"SMC"` (default), the diagonal of the correlation matrices is
#'    replaced by the squared multiple correlations (SMCs) of the indicators. If
#'     using  `"PCA"`, the diagonal values of the correlation
#'  matrices are left to be 1. If using  `"EFA"`, eigenvalues are found on the
#'  correlation  matrices with the final communalities of an EFA solution as
#'  diagonal. This is passed to  [PARALLEL()].
#' @param use character. Passed to [stats::cor()] if raw data
#' is given as input. Default is `"pairwise.complete.obs"`.
#' @param cor_method character. One of `"pearson"`, `"spearman"`, or `"kendall"`,
#'   passed to [stats::cor()]. `"poly"` and `"tetra"` are not supported because
#'   `HULL` derives its factor-search bound from an internal parallel analysis
#'   against continuous reference data.
#'  Default is  `"pearson"`.
#' @param n_datasets numeric. The number of datasets to simulate. Default is 1000.
#'   This is passed to [PARALLEL()].
#' @param percent numeric. The percentile to take from the simulated eigenvalues.
#'  Default is 95. This is passed to [PARALLEL()].
#' @param decision_rule character. Which rule to use to determine the number of
#' factors to retain. Default is `"means"`, which will use the average
#' simulated eigenvalues. `"percentile"`, uses the percentiles specified
#' in percent. `"crawford"` uses the 95th percentile for the first factor
#' and the mean afterwards (based on Crawford et al, 2010). This is passed to [PARALLEL()].
#' @param n_factors numeric. Number of factors to extract if  `"EFA"` is
#' included in `eigen_type`. Default is 1. This is passed to
#' [PARALLEL()].
#' @param ... Further arguments passed to [EFA()], also in
#' [PARALLEL()].
#'
#' @details The Hull method aims to find a model with an optimal balance between
#'  model fit and number of parameters. That is, it aims to retrieve only major
#'  factors (Lorenzo-Seva, Timmerman, & Kiers, 2011). To this end, it performs
#'  the following steps (Lorenzo-Seva, Timmerman, & Kiers, 2011, p.351):
#'  \enumerate{
#'    \item It performs parallel analysis and adds one to the identified number of factors (this number is denoted *J*). *J* is taken as an upper bound of the number of factors to retain in the hull method. Alternatively, a theoretical number of factors can be entered. In this case *J* will be set to whichever of these two numbers (from parallel analysis or based on theory) is higher.
#'    \item For all 0 to *J* factors, the goodness-of-fit (one of *CAF*, *RMSEA*, or *CFI*) and the degrees of freedom (*df*) are computed.
#'    \item The solutions are ordered according to their *df*.
#'    \item Solutions that are not on the boundary of the convex hull are eliminated (see Lorenzo-Seva, Timmerman, & Kiers, 2011, for details).
#'    \item All the triplets of adjacent solutions are considered consecutively. The middle solution is excluded if its point is below or on the line connecting its neighbors in a plot of the goodness-of-fit versus the degrees of freedom.
#'    \item Step 5 is repeated until no solution can be excluded.
#'    \item The *st* values of the “hull” solutions are determined.
#'    \item The solution with the highest *st* value is selected.
#'  }
#'
#' The [PARALLEL] function and the principal axis factoring of the
#'   different number of factors can be parallelized using the future framework,
#'   by calling the [future::plan()] function. The examples
#'    provide example code on how to enable parallel processing.
#'
#'   Note that if `gof = "RMSEA"` is used, 1 - RMSEA is actually used to
#'   compare the different solutions. Thus, the threshold of .05 is then .95. This
#'   is necessary due to how the heuristic to locate the elbow of the hull works.
#'
#'   The ML estimation method uses the [psych::fa()]
#'    starting values. See also the [EFA] documentation.
#'
#'    The `HULL` function can also be called together with other factor
#'    retention criteria in the [N_FACTORS()] function.
#' @returns An object of class `efa_retention` (see [print.efa_retention()] and
#'   [plot.efa_retention()] for the print and plot methods). Its main fields are:
#' \item{n_factors}{A named numeric vector with the suggested number of factors
#'   for each requested goodness-of-fit index (`"CAF"`, `"CFI"`, and/or
#'   `"RMSEA"`).}
#' \item{results}{A list with one record per goodness-of-fit index, each holding
#'   the goodness-of-fit values, the degrees of freedom, the hull membership, and
#'   the retained solution used for printing and plotting.}
#' \item{settings}{A list of the settings used, including `n_fac_max`, the upper
#'   bound *J* of the number of factors to extract (see details).}
#'
#' @source Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011).
#' The Hull method for selecting the number of common factors. Multivariate
#' Behavioral Research, 46(2), 340-364.
#'
#' @seealso Other factor retention criteria: [CD()], [EKC()],
#' [KGC()], [PARALLEL()], [SMT()]
#'
#' [N_FACTORS()] as a wrapper function for this and all the
#' above-mentioned factor retention criteria.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # using PAF (this will throw a warning if gof is not specified manually
#' # and CAF will be used automatically)
#' HULL(test_models$baseline$cormat, N = 500, gof = "CAF")
#'
#' # using ML with all available fit indices (CAF, CFI, and RMSEA)
#' HULL(test_models$baseline$cormat, N = 500, method = "ML")
#'
#' # using ULS with only RMSEA
#' HULL(test_models$baseline$cormat, N = 500, method = "ULS", gof = "RMSEA")
#'}
#'
#'\dontrun{
#' # using parallel processing (Note: plans can be adapted, see the future
#' # package for details)
#' future::plan(future::multisession)
#' HULL(test_models$baseline$cormat, N = 500, gof = "CAF")
#' }
HULL <- function(x, N = NA, n_fac_theor = NA,
                 method = c("PAF", "ULS", "ML"), gof = c("CAF", "CFI", "RMSEA"),
                 eigen_type = c("SMC", "PCA", "EFA"),
                 use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                         "everything", "na.or.complete"),
                 cor_method = c("pearson", "spearman", "kendall", "poly", "tetra"),
                 n_datasets = 1000, percent = 95,
                 decision_rule = c("means", "percentile", "crawford"),
                 n_factors = 1, ...) {
  # Perform hull method following Lorenzo-Seva, Timmerman, and Kiers (2011)

  .assert_cor_input(x)

  method <- match.arg(method)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  # The Hull method derives its factor-search bound from an internal parallel
  # analysis, whose reference data are continuous; poly/tetra are therefore not
  # supported, consistent with PARALLEL/NEST/CD.
  .reject_poly_reference(cor_method, "HULL")
  gof <- match.arg(gof, several.ok = TRUE)
  eigen_type <- match.arg(eigen_type)
  checkmate::assert_count(n_fac_theor, na.ok = TRUE)
  checkmate::assert_count(N, na.ok = TRUE)
  decision_rule <- match.arg(decision_rule)
  checkmate::assert_count(n_factors)
  checkmate::assert_count(n_datasets)
  checkmate::assert_number(percent, lower = 0, upper = 100)

  if (ncol(x) < 6) {
    cli::cli_abort(
      c("The data has fewer than 6 indicators.",
        "i" = "The Hull method needs at least 6."),
      class = "efa_hull_min_indicators"
    )
  }

  if (method == "PAF" && !all(gof == "CAF")) {
    cli::cli_alert_info(cli::col_cyan('Only CAF can be used as gof if method "PAF" is used. Setting gof to "CAF"\n'))
    gof <- "CAF"
  }

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(
    x, N = N, use = use, cor_method = cor_method, N_policy = "required",
    singular_tail = "the Hull method cannot be executed",
    N_required_msg = "{.arg N} is not specified but is needed to compute some fit indices.")
  R <- prep$R
  N <- prep$N

  m <- ncol(R)

  # 1) perform parallel analysis to find J as n_fac_theor + 1
  par_res <- PARALLEL(R, N = N, eigen_type = eigen_type, method = method,
                      n_datasets = n_datasets, percent = percent,
                      decision_rule = decision_rule, n_factors = n_factors,
                      ...)

  n_fac_PA <- unname(par_res$n_factors[eigen_type])

  if (is.na(n_fac_PA)) {

    if (!is.na(n_fac_theor)) {
      J <- n_fac_theor + 1
    } else {
      J <- .det_max_factors(ncol(R))
    }

  } else {

    J <- max(c(n_fac_PA, n_fac_theor), na.rm = TRUE) + 1

  }

  if (J > .det_max_factors(ncol(R))) {
    J <- .det_max_factors(ncol(R))
    cli::cli_warn("Setting the maximum number of factors to {J} to ensure overidentified models.",
                  class = "efa_hull_max_factors")
  }

  if (J < 3) {
    cli::cli_warn(
      c("The suggested maximum number of factors was {J}, but must be at least 3 for the Hull method.",
        "i" = "Setting it to 3."),
      class = "efa_hull_min_factors"
    )
    J <- 3

  }

  # 2) perform factor analysis for the range of dimensions 1:J and compute f and
  #    df for every solution
  s <- matrix(0, ncol = 4, nrow = J + 1)
  s[, 1] <- 0:J

  # first for 0 factors
  if ("CAF" %in% gof) {
    s_CAF <- s
    colnames(s_CAF) <- c("nfactors", "CAF", "df", "st")
    s_CAF[1, 2] <- 1 - KMO(R)$KMO
    s_CAF[1, 3] <- (m**2 - m) / 2
  }

  if ("CFI" %in% gof) {
    s_CFI <- s
    colnames(s_CFI) <- c("nfactors", "CFI", "df", "st")
    s_CFI[1, 2] <- 0
    s_CFI[1, 3] <- (m**2 - m) / 2
  }

  if ("RMSEA" %in% gof) {
    s_RMSEA <- s
    colnames(s_RMSEA) <- c("nfactors", "RMSEA", "df", "st")
    # 0-factor (independence model) reference, on the same Bartlett-corrected
    # ML-discrepancy scale as the 1:J solutions returned by .gof()
    chi <- .null_chisq(R, N)
    df <- (m**2 - m) / 2
    # compute 1 - RMSEA
    s_RMSEA[1, 2] <- 1 - sqrt(max(0, chi - df) / (df * (N - 1)))
    s_RMSEA[1, 3] <- (m**2 - m) / 2

  }

  # Calculate loadings with EFA function
  loadings <- suppressWarnings(future.apply::future_lapply(seq_len(J), EFA,
                                                           x = R,
                                                           method = method,
                                                           N = N, ...,
                                                           future.seed = FALSE))

  # then for 1 to J factors
  for (i in seq_len(J)) {
    if (method == "PAF") {
      # compute goodness of fit "f" as CAF (common part accounted for; Eq 3)
      # compute CAF
      s_CAF[i + 1, 2] <- loadings[[i]]$fit_indices$CAF
      # compute dfs (Eq 4 provides the number of free parameters; using dfs yields
      # th same numbers, as the difference in df equals the difference in free
      # parameters)
      s_CAF[i + 1, 3] <- loadings[[i]]$fit_indices$df
    } else {
      if ("CAF" %in% gof) {
        # compute goodness of fit "f" as CAF (common part accounted for; Eq 3)
        # compute CAF
        s_CAF[i + 1, 2] <- loadings[[i]]$fit_indices$CAF
        # compute dfs (Eq 4 provides the number of free parameters; using dfs yields
        # th same numbers, as the difference in df equals the difference in free
        # parameters)
        s_CAF[i + 1, 3] <- loadings[[i]]$fit_indices$df
      }

      if ("CFI" %in% gof) {
        # compute CFI
        s_CFI[i + 1, 2] <- loadings[[i]]$fit_indices$CFI
        # compute dfs (Eq 4 provides the number of free parameters; using dfs yields
        # th same numbers, as the difference in df equals the difference in free
        # parameters)
        s_CFI[i + 1, 3] <- loadings[[i]]$fit_indices$df

      }

      if ("RMSEA" %in% gof) {
        # compute 1 - RMSEA
        s_RMSEA[i + 1, 2] <- 1 - loadings[[i]]$fit_indices$RMSEA
        # compute dfs (Eq 4 provides the number of free parameters; using dfs yields
        # th same numbers, as the difference in df equals the difference in free
        # parameters)
        s_RMSEA[i + 1, 3] <- loadings[[i]]$fit_indices$df

      }

    }

  }

  out_CAF <- list(s_complete = NA, retain = NA)
  out_CFI <- list(s_complete = NA, retain = NA)
  out_RMSEA <- list(s_complete = NA, retain = NA)

  if("CAF" %in% gof) {
    out_CAF <- .hull_calc(s = s_CAF, J = J, gof_t = "CAF")
  }
  if("CFI" %in% gof) {
    out_CFI <- .hull_calc(s = s_CFI, J = J, gof_t = "CFI")
  }
  if("RMSEA" %in% gof) {
    out_RMSEA <- .hull_calc(s = s_RMSEA, J = J, gof_t = "RMSEA")
  }

  gof_results <- list(CAF = out_CAF, CFI = out_CFI, RMSEA = out_RMSEA)

  # one record per requested goodness-of-fit index (df vs. fit, with the hull
  # membership and the retained solution used for printing and plotting)
  results <- list()
  for (g in c("CAF", "CFI", "RMSEA")) {
    if (!(g %in% gof)) next
    sol <- gof_results[[g]]$s_complete
    retain <- gof_results[[g]]$retain
    results[[g]] <- list(
      name = g,
      label = g,
      n_factors = retain,
      plot_type = "hull",
      x = unname(sol[, "df"]),
      y = unname(sol[, g]),
      reference = NULL,
      threshold = NULL,
      highlight = retain,
      point_labels = unname(sol[, "nfactors"]),
      on_hull = !is.na(sol[, "st"])
    )
  }

  out <- .new_efa_retention(
    "HULL",
    results = unname(results),
    settings = list(N = N,
                    method = method,
                    gof = gof,
                    n_fac_theor = n_fac_theor,
                    eigen_type = eigen_type,
                    use = use,
                    cor_method = cor_method,
                    n_fac_max = J),
    subtitle = paste0("Estimation method: ", method)
  )

  return(out)

}


.hull_calc <- function(s, J, gof_t){

  # 3) sort n solutions by their df values and denoted by s (already done)

  # 4) all solutions s are excluded for which a solution sj (j<i) exists such
  #    that fj > fi (eliminate solutions not on the boundary of the convex hull)

  s_complete <- s

  # A non-finite goodness-of-fit value (e.g. an undefined CFI/RMSEA for a Heywood
  # or near-singular model) cannot lie on the convex hull; drop those solutions
  # with a classed warning rather than failing in the comparisons below.
  na_rows <- !is.finite(s[, 2])
  if (any(na_rows)) {
    cli::cli_warn(
      c("{sum(na_rows)} solution{?s} had a non-finite {gof_t} value and {?was/were} excluded from the hull.",
        "i" = "Inspect the affected models, or try a different goodness-of-fit index or estimation method."),
      class = "efa_hull_na_fit"
    )
    s <- s[!na_rows, , drop = FALSE]
  }

  if (nrow(s) < 1) {
    cli::cli_abort(
      c("No solution had a finite {gof_t} value, so the Hull method cannot proceed.",
        "i" = "Try a different goodness-of-fit index or estimation method."),
      class = "efa_hull_no_fit"
    )
  }

  d_s <- diff(s[, 2])
  while (any(d_s < 0)) {
    s <- s[c(1, d_s) > 0, , drop = FALSE]
    if(nrow(s) == 1){
      break
    }
    d_s <- diff(s[, 2])
  }

  # 5) all triplets of adjacent solutions are considered consecutively.
  #    The middle solution is excluded if its point is below or on the line
  #    connecting its neighbors in GOF vs df.

  # 6) repeat 5) until no solution can be excluded

  # The `i <= nr_s - 1` bound ensures the final interior triplet is also tested;
  # the loop is a no-op while fewer than three boundary solutions remain.
  nr_s <- nrow(s)
  i <- 2

  while(i <= nr_s - 1) {

    f1 <- s[i - 1, 2]
    f2 <- s[i, 2]
    f3 <- s[i + 1, 2]
    df1 <- s[i - 1, 3]
    df2 <- s[i, 3]
    df3 <- s[i + 1, 3]

    # compute f2 if it were on the line between f1 and f3
    p_f2 <- f1 + (f3 - f1) / (df3 - df1) * (df2 - df1)

    # check if f2 is below or on the predicted line and if so, remove it
    if (f2 <= p_f2) {
      s <- s[-i, , drop = FALSE]
      nr_s <- nr_s -1
      i <- 1
    }
    i <- i + 1
  }

  if (nrow(s) < 3) {
    cli::cli_warn(
      c("Fewer than three solutions were located on the hull using {gof_t} as goodness-of-fit index.",
        "i" = "Proceeding with the maximum-{gof_t} value as a heuristic; consider additional indices or methods as a robustness check."),
      class = "efa_hull_few_solutions"
    )

    # the st values are undefined with fewer than three hull solutions
    s_complete[, 4] <- NA

    # 8) select solution with highest gof value
    retain <- s[which.max(s[, 2]), 1]


  } else {

    # 7) the st values of the hull solutions are determined (Eq 5)
    for (i in 2:(nrow(s) - 1)) {

      f_i <- s[i, 2]
      f_p <- s[i - 1, 2]
      f_n <- s[i + 1, 2]
      df_i <- s[i, 3]
      df_p <- s[i - 1, 3]
      df_n <- s[i + 1, 3]

      s[i, 4] <- ((f_i - f_p) / (df_i - df_p)) / ((f_n - f_i) / (df_n - df_i))

    }

    # combine values
    for (row_i in 0:J) {

      if (row_i %in% s[,1]) {
        s_complete[row_i + 1, 4] <- s[s[,1] == row_i, 4]
      } else {
        s_complete[row_i + 1, 4] <- NA
      }

    }

    # 8) select solution with highest st value
    retain <- s[which.max(s[, 4]), 1]

  }


  out <- list(s_complete = s_complete,
              retain = unname(retain))

  return(out)

}
