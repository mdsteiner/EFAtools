#' Various Factor Retention Criteria
#'
#' Among the most important decisions for an exploratory factor analysis (EFA) is
#' the choice of the number of factors to retain. Several factor retention
#' criteria have been developed for this. With this function, various factor
#'  retention criteria can be performed simultaneously. Additionally, the data
#'  can be checked for their suitability for factor analysis.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations. If \code{"CD"} is included as a criterion, x must be raw
#'  data.
#' @param criteria character. A vector with the factor retention methods to
#' perform. Possible inputs are: \code{"CD"}, \code{"EKC"}, \code{"HULL"},
#' \code{"KGC"}, \code{"PARALLEL"}, \code{"SCREE"}, and \code{"SMT"}
#' (see details). By default, all factor retention methods are performed.
#' @param suitability logical. Whether the data should be checked for suitability
#' for factor analysis using the Bartlett's test of sphericity and the
#' Kaiser-Guttmann criterion (see details). Default is \code{TRUE}.
#' @param N  numeric. The number of observations. Only needed if x is a
#' correlation matrix.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw
#' data is given as input. Default is \code{"pairwise.complete.obs"}.
#' @param cor_method character. Passed to \code{\link[stats:cor]{stats::cor}}
#' Default is  \code{"pearson"}.
#' @param n_factors_max numeric. Passed to \code{\link{CD}}.The maximum number
#' of factors to test against.
#' Larger numbers will increase the duration the procedure takes, but test more
#' possible solutions. Maximum possible is number of variables / 2. Default is
#' NA. If not specified, number of variables / 2 is used.
#' @param N_pop numeric. Passed to \code{\link{CD}}. Size of finite populations
#' of comparison data. Default is 10000.
#' @param N_samples numeric. Passed to \code{\link{CD}}. Number of samples drawn
#'  from each population. Default is 500.
#' @param alpha numeric. Passed to \code{\link{CD}}. The alpha level used to test
#'  the significance of the improvement added by an additional factor.
#'  Default is .30.
#' @param max_iter_CD numeric. Passed to \code{\link{CD}}. The maximum number of
#'  iterations to perform after which the iterative PAF procedure is halted.
#'   Default is 50.
#' @param n_fac_theor numeric. Passed to \code{\link{HULL}}. Theoretical number
#'  of factors to retain. The maximum of this number and the number of factors
#'  suggested by \link{PARALLEL} plus one will be used in the Hull method.
#' @param method character. Passed to \code{\link{EFA}} in \code{\link{HULL}},
#' \code{\link{KGC}}, \code{\link{SCREE}}, and \code{\link{PARALLEL}}. The
#' estimation method to use. One of  \code{"PAF"}, \code{"ULS"}, or  \code{"ML"},
#' for principal axis factoring, unweighted least squares, and maximum
#' likelihood, respectively.
#' @param gof character. Passed to \code{\link{HULL}}. The goodness of fit index
#' to use. Either \code{"CAF"}, \code{"CFI"}, or \code{"RMSEA"}, or any
#' combination of them. If \code{method = "PAF"} is used, only
#' the CAF can be used as goodness of fit index. For details on the CAF, see
#' Lorenzo-Seva, Timmerman, and Kiers (2011).
#' @param eigen_type_HULL character. Passed to  \code{\link{PARALLEL}} in
#' \code{\link{HULL}}. On what the
#' eigenvalues should be found in the parallel analysis. Can be one of
#' \code{"SMC"}, \code{"PCA"}, or \code{"EFA"}. If using  \code{"SMC"} (default),
#' the diagonal of the correlation matrices is
#' replaced by the squared multiple correlations (SMCs) of the indicators. If
#' using  \code{"PCA"}, the diagonal values of the correlation
#' matrices are left to be 1. If using  \code{"EFA"}, eigenvalues are found on the
#' correlation  matrices with the final communalities of an EFA solution as
#' diagonal.
#' @param eigen_type_other character. Passed to \code{\link{KGC}},
#' \code{\link{SCREE}}, and \code{\link{PARALLEL}}. The same as eigen_type_HULL,
#' but multiple inputs
#' are possible here. Default is to use all inputs, that is, \code{c("PCA",
#' "SMC", "EFA"})
#' @param n_factors numeric. Passed to \code{\link{PARALLEL}} (also within
#' \code{\link{HULL}}), \code{\link{KGC}}, and \code{\link{SCREE}}. Number of
#' factors to extract if \code{"EFA"} is included in \code{eigen_type_HULL} or
#'  \code{eigen_type_other}. Default is 1.
#' @param n_datasets numeric. Passed to \code{\link{PARALLEL}} (also within
#' \code{\link{HULL}}). The number of datasets to simulate. Default is 1000.
#' @param percent numeric. Passed to \code{\link{PARALLEL}} (also within
#' \code{\link{HULL}}). A vector of percentiles to take the simulated eigenvalues
#'  from. Default is 95.
#' @param decision_rule character. Passed to \code{\link{PARALLEL}} (also within
#' \code{\link{HULL}}). Which rule to use to determine the number of
#'  factors to retain. Default is \code{"means"}, which will use the average
#'  simulated eigenvalues. \code{"percentile"}, uses the percentiles specified
#'  in percent. \code{"crawford"} uses the 95th percentile for the first factor
#'  and the mean afterwards (based on Crawford et al, 2010).
#' @param show_progress logical. Whether a progress bar should be shown in the
#'   console. Default is TRUE.
#' @param ... Further arguments passed to \code{\link{EFA}} in
#' \code{\link{PARALLEL}} (also within \code{\link{HULL}}) and \code{\link{KGC}}.
#'
#' @details
#' By default, the entered data are checked for suitability for factor analysis
#' using the following methods (see respective documentations for details):
#' \itemize{
#' \item{Bartlett's test of sphericity (see \code{\link{BARTLETT}})}
#' \item{Kaiser-Meyer-Olkin criterion (see \code{\link[EFAtools]{KMO}})}}
#'
#' The available factor retention criteria are the following (see respective
#'  documentations for details):
#'  \itemize{
#' \item{Comparison data (see \code{\link{CD}})}
#' \item{Empirical Kaiser criterion (see \code{\link{EKC}})}
#' \item{Hull method (see \code{\link{HULL}})}
#' \item{Kaiser-Guttman criterion (see \code{\link{KGC}})}
#' \item{Parallel analysis (see \code{\link{PARALLEL}})}
#' \item{Scree plot (see \code{\link{SCREE}})}
#' \item{Sequential chi-square model tests, RMSEA lower bound, and AIC
#' (see \code{\link{SMT}})}
#' }
#'
#' @return A list of class N_FACTORS containing
#' \item{outputs}{A list with the outputs from \code{\link{BARTLETT}} and
#'  \code{\link[EFAtools]{KMO}} and the factor retention criteria.}
#' \item{n_factors}{A named vector containing the suggested number of factors
#' from each factor retention criterion.}
#' \item{settings}{A list of the settings used.}
#'
#' @export
#'
#' @examples
#' \donttest{
#' # All criteria, with correlation matrix and fit method "ML" (where needed)
#' # This will throw a warning for CD, as no raw data were specified
#' nfac_all <- N_FACTORS(test_models$baseline$cormat, N = 500, method = "ML")
#'
#' # The same as above, but without "CD"
#' nfac_wo_CD <- N_FACTORS(test_models$baseline$cormat, criteria = c("EKC",
#'                         "HULL", "KGC", "PARALLEL", "SCREE", "SMT"), N = 500,
#'                         method = "ML")
#'
#' # Use PAF instead of ML (this will take a lot longer). For this, gof has
#' # to be set to "CAF" for the Hull method.
#' nfac_PAF <- N_FACTORS(test_models$baseline$cormat, criteria = c("EKC",
#'                       "HULL", "KGC", "PARALLEL", "SCREE", "SMT"), N = 500,
#'                       gof = "CAF")
#'
#' # Do KGC and PARALLEL with only "PCA" type of eigenvalues
#' nfac_PCA <- N_FACTORS(test_models$baseline$cormat, criteria = c("EKC",
#'                       "HULL", "KGC", "PARALLEL", "SCREE", "SMT"), N = 500,
#'                       method = "ML", eigen_type_other = "PCA")
#'
#' # Use raw data, such that CD can also be performed
#' nfac_raw <- N_FACTORS(GRiPS_raw, method = "ML")
#'}
N_FACTORS <- function(x, criteria = c("CD", "EKC", "HULL", "KGC", "PARALLEL",
                                      "SCREE", "SMT"),
                      suitability = TRUE, N = NA,
                      use = c("pairwise.complete.obs", "all.obs",
                              "complete.obs", "everything", "na.or.complete"),
                      cor_method = c("pearson", "spearman", "kendall"),
                      n_factors_max = NA, N_pop = 10000, N_samples = 500,
                      alpha = .30, max_iter_CD = 50, n_fac_theor = NA,
                      method = c("PAF", "ULS", "ML"),
                      gof = c("CAF", "CFI", "RMSEA"),
                      eigen_type_HULL = c("SMC", "PCA", "EFA"),
                      eigen_type_other = c("PCA", "SMC", "EFA"),
                      n_factors = 1, n_datasets = 1000,
                      percent = 95,
                      decision_rule = c("means", "percentile", "crawford"),
                      show_progress = TRUE,
                      ...){

  # Perform argument checks
  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n"))

  }

  ## Perform argument checks and prepare input
  criteria <- match.arg(criteria, several.ok = TRUE)
  suitability <- checkmate::assert_flag(suitability)
  eigen_type_HULL <- match.arg(eigen_type_HULL)
  eigen_type_other <- match.arg(eigen_type_other, several.ok = TRUE)
  cor_method <- match.arg(cor_method)
  use <- match.arg(use)
  method <- match.arg(method)
  decision_rule <- match.arg(decision_rule)

  if (isTRUE(show_progress)) {
    criteria <- sort(criteria)
  }

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    R <- x

  } else {

    if (!is.na(N)) {
      warning(crayon::yellow$bold("!"), crayon::yellow(" 'N' was set and data entered. Taking N from data.\n"))
    }

    R <- stats::cor(x, use = use, method = cor_method)
    colnames(R) <- colnames(x)
    N <- nrow(x)

  }

  # Check if correlation matrix is invertable, if it is not, stop with message
  R_i <- try(solve(R), silent = TRUE)

  if (inherits(R_i, "try-error")) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Correlation matrix is singular, no further analyses are performed\n"))
  }

  # Check if correlation matrix is positive definite
  if(any(eigen(R, symmetric = TRUE, only.values = TRUE)$values <= 0)){

    R <- psych::cor.smooth(R)

  }

  # Set all outputs to NA for a start
  bart_out <- NA
  kmo_out <- NA
  cd_out <- NA
  ekc_out <- NA
  hull_out <- NA
  kgc_out <- NA
  parallel_out <- NA
  scree_out <- NA
  smt_out <- NA

  nfac_CD <- NA
  nfac_EKC <- NA
  nfac_HULL_CAF <- NA
  nfac_HULL_CFI <- NA
  nfac_HULL_RMSEA <- NA
  nfac_KGC_PCA <- NA
  nfac_KGC_SMC <- NA
  nfac_KGC_EFA <- NA
  nfac_PA_PCA <- NA
  nfac_PA_SMC <- NA
  nfac_PA_EFA <- NA
  nfac_SMT_chi <- NA
  nfac_RMSEA <- NA
  nfac_AIC <- NA

  ## Tests for suitability of factor analysis
  if(isTRUE(suitability)){

  # Bartlett's Test of Sphericity
  bart_out <- BARTLETT(R, N = N, use = use, cor_method = cor_method)

  # Kaiser-Meyer_Olkin criterion
  kmo_out <- KMO(R, use = use, cor_method = cor_method)

  }

  ## Factor retention criteria

  # Comparison data
  if("CD" %in% criteria){

    if (!.is_cormat(x)) {

      if (isTRUE(show_progress)) {
        .show_progress(criteria, "CD")
      }

      cd_out <- CD(x, n_factors_max = n_factors_max, N_pop = N_pop,
                   N_samples = N_samples, alpha = alpha, use = use,
                   cor_method = cor_method, max_iter = max_iter_CD)

      nfac_CD <- cd_out$n_factors

    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" 'x' was a correlation matrix but CD needs raw data. Skipping CD.\n"))
    }

  }

  # Empirical Kaiser Criterion
  if("EKC" %in% criteria){

    if (isTRUE(show_progress)) {
      .show_progress(criteria, "EKC")
    }

    ekc_out <- EKC(R, N = N, use = use, cor_method = cor_method)

    nfac_EKC <- ekc_out$n_factors

  }

  # HULL method
  if("HULL" %in% criteria){

    if (isTRUE(show_progress)) {
      .show_progress(criteria, "HULL")
    }

    hull_out <- HULL(R, N = N, n_fac_theor = n_fac_theor, method = method,
                     eigen_type = eigen_type_HULL, gof = gof, use = use,
                     cor_method = cor_method, n_datasets = n_datasets,
                     percent = percent, decision_rule = decision_rule,
                     n_factors = n_factors, ...)

    nfac_HULL_CAF <- hull_out$n_fac_CAF
    nfac_HULL_CFI <- hull_out$n_fac_CFI
    nfac_HULL_RMSEA <- hull_out$n_fac_RMSEA

  }

  # Kaiser-Guttman criterion
  if("KGC" %in% criteria){

    if (isTRUE(show_progress)) {
      .show_progress(criteria, "KGC")
    }

    kgc_out <- KGC(R, eigen_type = eigen_type_other, use = use,
                   cor_method = cor_method, n_factors = n_factors,
                   method = method, ...)

    nfac_KGC_PCA <- kgc_out$n_fac_PCA
    nfac_KGC_SMC <- kgc_out$n_fac_SMC
    nfac_KGC_EFA <- kgc_out$n_fac_EFA

  }

  # Parallel analysis
  if("PARALLEL" %in% criteria){

    if (isTRUE(show_progress)) {
      .show_progress(criteria, "PARALLEL")
    }

    parallel_out <- try(PARALLEL(R, N = N,
                                 n_datasets = n_datasets, percent = percent,
                                 eigen_type = eigen_type_other, use = use,
                                 cor_method = cor_method,
                                 decision_rule = decision_rule,
                                 n_factors = n_factors, method = method, ...))

    nfac_PA_PCA <- parallel_out$n_fac_PCA
    nfac_PA_SMC <- parallel_out$n_fac_SMC
    nfac_PA_EFA <- parallel_out$n_fac_EFA

  }

  # Scree plot
  if("SCREE" %in% criteria){

    if (isTRUE(show_progress)) {
      .show_progress(criteria, "SCREE")
    }

    scree_out <- SCREE(R, eigen_type = eigen_type_other, use = use,
                     cor_method = cor_method, n_factors = n_factors,
                     method = method, ...)
  }

  # Sequential chi square tests, RMSEA lower bound and AIC
  if("SMT" %in% criteria){

    if (isTRUE(show_progress)) {
      .show_progress(criteria, "SMT")
    }

    smt_out <- SMT(R, N = N, use = use, cor_method = cor_method)

    nfac_SMT_chi <- smt_out$nfac_chi
    nfac_RMSEA <- smt_out$nfac_RMSEA
    nfac_AIC <- smt_out$nfac_AIC

  }

  # Prepare settings here
  settings <- list(criteria = criteria,
                   suitability = suitability,
                   N = N,
                   use = use,
                   n_factors_max = n_factors_max,
                   N_pop = N_pop,
                   N_samples = N_samples,
                   alpha = alpha,
                   cor_method = cor_method,
                   max_iter_CD = max_iter_CD,
                   n_fac_theor = n_fac_theor,
                   method = method,
                   gof = gof,
                   eigen_type_HULL = eigen_type_HULL,
                   eigen_type_other = eigen_type_other,
                   n_factors = n_factors,
                   n_datasets = n_datasets,
                   percent = percent,
                   decision_rule = decision_rule)

  # Prepare the output
  n_factors <- c(nfac_CD = nfac_CD,
                 nfac_EKC = nfac_EKC,
                 nfac_HULL_CAF = nfac_HULL_CAF,
                 nfac_HULL_CFI = nfac_HULL_CFI,
                 nfac_HULL_RMSEA = nfac_HULL_RMSEA,
                 nfac_KGC_PCA = nfac_KGC_PCA,
                 nfac_KGC_SMC = nfac_KGC_SMC,
                 nfac_KGC_EFA = nfac_KGC_EFA,
                 nfac_PA_PCA = nfac_PA_PCA,
                 nfac_PA_SMC = nfac_PA_SMC,
                 nfac_PA_EFA = nfac_PA_EFA,
                 nfac_SMT_chi = nfac_SMT_chi,
                 nfac_RMSEA = nfac_RMSEA,
                 nfac_AIC = nfac_AIC)

  outputs <- list(bart_out = bart_out,
                  kmo_out = kmo_out,
                  cd_out = cd_out,
                  ekc_out = ekc_out,
                  hull_out = hull_out,
                  kgc_out = kgc_out,
                  parallel_out = parallel_out,
                  scree_out = scree_out,
                  smt_out = smt_out)

  output <- list(outputs = outputs,
                 n_factors = n_factors,
                 settings = settings)

  class(output) <- "N_FACTORS"

  if (isTRUE(show_progress)) {
    .show_progress(criteria, "done", TRUE)
  }

  return(output)

}
