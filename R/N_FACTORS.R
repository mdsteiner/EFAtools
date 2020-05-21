
# - say what only works with raw data (CD, PARALLEL if resampling is implemented)
# - ... further arguments passed to PARALLEL in HULL or passed to EFA in KGC, PARALLEL and HULL ->
    # for these functions, methods argument here is also passed to EFA, problem?
# - ok if same eigen_type is used for PARALELL and KGC?
# arguments set for parallel are also used in HULL parallel (because passed to
# PARALLEL there), problem?
N_FACTORS <- function(x, criteria = c("CD", "EKC", "HULL", "KGC", "MACHINE",
                                      "PARALLEL", "SMT"),
                      suitability = c("KMO", "BARTLETT"), N = NA,
                      use = c("pairwise.complete.obs", "all.obs",
                              "complete.obs", "everything", "na.or.complete"),
                      n_factors_max = NA, N_pop = 10000, N_samples = 500,
                      alpha = .30, cor_method = c("pearson", "spearman",
                                                  "kendall"),
                      max_iter = 1000, n_fac_theor = NA,
                      method = c("PAF", "ULS", "ML"),
                      gof = c("CAF", "CFI", "RMSEA"),
                      eigen_type = c("PCA", "SMC", "EFA"), n_factors = 1,
                      n_vars = NA, n_datasets = 1000, percent = 95,
                      data_type = c("sim"), # , "resample"
                      replace = TRUE, decision_rule = c("Means", "Percentile",
                                                        "Crawford"),
                      ...){

  ## Perform argument checks and prepare input
  criteria <- match.arg(criteria, several.ok = TRUE)

  # WHAT HAPPENS IF YOU ENTER NULL? SHOULD BE POSSIBLE...
  suitability <- match.arg(suitability, several.ok = TRUE)
  # OTHER CHECKS IN RESPECTIVE FUNCTIONS OK?
  # OR BETTER HERE?

  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    if(any(is.na(x))){

      stop("The correlation matrix you entered contains missing values. No
           further analyses are possible.")

    }

    R <- x

  } else {

    message("x was not a correlation matrix. Correlations are found from entered
            raw data.")

    R <- stats::cor(x, use = use)
    colnames(R) <- colnames(x)
    N <- nrow(x)

  }

  # Check if correlation matrix is invertable, if it is not, stop with message
  R_i <- try(solve(R))

  if (inherits(R_i, "try-error")) {
    stop("Correlation matrix is singular, No further analyses are possible")
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

  # Bartlett's Test of Sphericity
  if("BARTLETT" %in% suitability){

    bart_out <- BARTLETT(R, N = N, use = use)

  }

  # Kaiser-Meyer_Olkin criterion
  if("KMO" %in% suitability){

    kmo_out <- KMO(R, use = use)

  }

  ## Factor retention criteria

  # Comparison data
  if("CD" %in% criteria){

    # n_factors and method arguments needed?
    cd_out <- try(CD(x, n_factors_max = n_factors_max, N_pop = N_pop,
                     N_samples = N_samples, alpha = alpha, use = use,
                     cor_method = cor_method,
                     max_iter = max_iter, # really needed?
                     ... #(to pass to EFA?)
                     ))

    if(!inherits(cd_out, "try-error")){

      nfac_CD <- cd_out$n_factors

    }

  }

  # Empirical Kaiser Criterion
  if("EKC" %in% criteria){

    ekc_out <- EKC(x, N = N, use = use)

    nfac_EKC <- ekc_out$n_factors

  }

  # HULL method
  if("HULL" %in% criteria){

    hull_out <- HULL(R, N = N, n_fac_theor = n_fac_theor, n_factors = n_factors,
                     method = method, gof = gof, use = use, ...)

    nfac_HULL_CAF <- hull_out$n_fac_CAF
    nfac_HULL_CFI <- hull_out$n_fac_CFI
    nfac_HULL_RMSEA <- hull_out$n_fac_RMSEA

  }

  # Kaiser-Guttman criterion
  if("KGC" %in% criteria){

    kgc_out <- KGC(R, eigen_type = eigen_type, use = use, n_factors = n_factors,
                   method = method, ...)

    nfac_KGC_PCA <- kgc_out$n_fac_PCA
    nfac_KGC_SMC <- kgc_out$n_fac_SMC
    nfac_KGC_EFA <- kgc_out$n_fac_EFA

  }

  # MACHINE LEARNING ALGORITHM (TO BE DETERMINED AND ADAPTED!!!)
  if("MACHINE" %in% criteria){
  }

  # Parallel analysis
  if("PARALLEL" %in% criteria){

    # IF RESAMPLING IS NOT IMPLEMENTED, USE R INSTEAD OF X HERE AND REMOVE TRY
    parallel_out <- try(PARALLEL(x, N = N, n_vars = n_vars,
                                 n_datasets = n_datasets, percent = percent,
                                 eigen_type = eigen_type, data_type = data_type,
                                 replace = replace, use = use,
                                 decision_rule = decision_rule,
                                 n_factors = n_factors, method = method,
                                 max_iter = max_iter, ... # really needed?
                                 ))

    nfac_PA_PCA <- parallel_out$n_fac_PCA
    nfac_PA_SMC <- parallel_out$n_fac_SMC
    nfac_PA_EFA <- parallel_out$n_fac_EFA

  }

  # Sequential chi square tests, RMSEA lower bound and AIC
  if("SMT" %in% criteria){

    smt_out <- SMT(R, N = N, use = use)

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
                   max_iter = max_iter, ### REALLY NEEDED?
                   n_fac_theor = n_fac_theor,
                   method = method,
                   gof = gof,
                   eigen_type = eigen_type,
                   n_factors = n_factors,
                   n_vars = n_vars,
                   n_datasets = n_datasets,
                   percent = percent,
                   data_type = data_type, # needed if no resampling?
                   replace = replace, # needed if no resampling?
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

  ## BETTER PUT ALL OUTPUTS IN A LIST FIRST?
  output <- list(bart_out = bart_out,
                 kmo_out = kmo_out,
                 cd_out = cd_out,
                 ekc_out = ekc_out,
                 hull_out = hull_out,
                 kgc_out = kgc_out,
                 parallel_out = parallel_out,
                 smt_out = smt_out,
                 n_factors = n_factors,
                 settings = settings)

  class(output) <- "N_FACTORS"

  return(output)

}
