#' Print function for N_FACTORS objects
#'
#' @param x a list of class N_FACTORS Output from \link{N_FACTORS} function.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print N_FACTORS
#'
#' @examples
#'
print.N_FACTORS <- function(x, ...){

  suitability <- x$settings$suitability
  criteria <- x$settings$criteria
  n_fac <- x$n_factors

  if(!is.null(suitability)){

    cat("\n")
    cat(crayon::blue$bold("Tests for the suitability of the data for factor analysis"))
    cat("\n")

    if("BARTLETT" %in% suitability){
    print(x$bart_out)
    cat("\n")
    }

    if("KMO" %in% suitability){
    print(x$kmo_out)
    cat("\n")
    cat("\n")
    }

  }

  cat(crayon::blue$bold("Number of factors suggested by the different factor",
  "retention criteria"))
  cat("\n")
  cat("\n")

  if("CD" %in% criteria){

  cat("Comparison data: ", crayon::bold(n_fac["nfac_CD"]), sep = "")
  cat("\n")

  plot(x$cd_out)

  }

  if("EKC" %in% criteria){

    cat("Empirical Kaiser criterion: ", crayon::bold(n_fac["nfac_EKC"]), sep = "")
    cat("\n")

    plot(x$ekc_out)

  }

  if("HULL" %in% criteria){

    ### ADDITIONAL CHECK FOR WHICH METHOD USED? OR JUST NA IF NOT USED?
    cat("HULL method with CAF: ", crayon::bold(n_fac["nfac_HULL_CAF"]), sep = "")
    cat("\n")
    cat("HULL method with CFI: ", crayon::bold(n_fac["nfac_HULL_CFI"]), sep = "")
    cat("\n")
    cat("HULL method with RMSEA: ", crayon::bold(n_fac["nfac_HULL_RMSEA"]), sep = "")
    cat("\n")

    plot(x$hull_out)

  }

  if("KGC" %in% criteria){

    ### ADDITIONAL CHECK FOR WHICH METHOD USED? OR JUST NA IF NOT USED?
    cat("Kaiser-Guttman criterion with PCA: ",
        crayon::bold(n_fac["nfac_KGC_PCA"]), sep = "")
    cat("\n")
    cat("Kaiser-Guttman criterion with SMC: ",
        crayon::bold(n_fac["nfac_KGC_SMC"]), sep = "")
    cat("\n")
    cat("Kaiser-Guttman criterion with EFA: ",
        crayon::bold(n_fac["nfac_KGC_EFA"]), sep = "")
    cat("\n")

    plot(x$kgc_out)

  }

  if("PARALLEL" %in% criteria){

    ### ADDITIONAL CHECK FOR WHICH METHOD USED? OR JUST NA IF NOT USED?
    cat("Parallel analysis with PCA: ", crayon::bold(n_fac["nfac_PA_PCA"]), sep = "")
    cat("\n")
    cat("Parallel analysis with SMC: ", crayon::bold(n_fac["nfac_PA_SMC"]), sep = "")
    cat("\n")
    cat("Parallel analysis with EFA: ", crayon::bold(n_fac["nfac_PA_EFA"]), sep = "")
    cat("\n")

    plot(x$parallel_out)

  }

  if("SMT" %in% criteria){

    cat("Sequential \U1D712\U00B2 model tests: ", crayon::bold(n_fac["nfac_SMT_chi"]),
        sep = "")
    cat("\n")
    cat("Lower bound of RMSEA 90% confidence interval: ",
        crayon::bold(n_fac["nfac_RMSEA"]), sep = "")
    cat("\n")
    cat("Akaike Information Criterion: ", crayon::bold(n_fac["nfac_AIC"]), sep = "")
    cat("\n")

  }

}
