#' Print function for N_FACTORS objects
#'
#' @param x a list of class N_FACTORS. Output from \link{N_FACTORS} function.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print N_FACTORS
#'
#' @examples
#' # All criteria except "CD", with correlation matrix and fit method "ML"
#' (where needed)
#' N_FACTORS(test_models$baseline$cormat, criteria = c("EKC", "HULL", "KGC",
#'           "PARALLEL", "SMT"), N = 500, method = "ML")
#'
print.N_FACTORS <- function(x, ...){

  suitability <- x$settings$suitability
  criteria <- x$settings$criteria
  n_fac <- x$n_factors
  gof <- x$settings$gof
  eigen_type_KGC_PA <- x$settings$eigen_type_KGC_PA
  kmo_out <- x$outputs$kmo_out
  KMO <- kmo_out$KMO

  if(isTRUE(suitability)){

    cat("\n")
    cat(crayon::blue$bold("Tests for the suitability of the data for factor analysis"))
    cat("\n")

    print(x$output$bart_out)
    cat("\n")

    cat("\n")
    cat(crayon::bold("Kaiser-Meyer-Olkin criterion (KMO)"))
    cat("\n")

    if(!is.na(KMO) && !is.null(KMO)){

      if(KMO >= .9){
        label = crayon::green$bold("marvellous")
      } else if(KMO >= .8){
        label = crayon::green$bold("meritorious")
      } else if(KMO >= .7){
        label = crayon::green$bold("middling")
      } else if(KMO >= .6){
        label = crayon::yellow$bold("mediocre")
      } else if (KMO >= .5){
        label = crayon::red$bold("miserable")
      } else {
        label = crayon::red$bold("unacceptable")
      }

      cat("\n")
      cat("The overall KMO value for your data is ", label, " with ",
      crayon::bold(round(KMO, 3)), ".", sep = "")
      cat("\n")
      cat("\n")

      if(KMO < .5){
        cat("These data are not suitable for factor analysis.")
        cat("\n")
        cat("\n")
        cat("\n")
      } else if(KMO < .6){
        cat("These data are hardly suitable for factor analysis.")
        cat("\n")
        cat("\n")
        cat("\n")
      } else {
        cat("\n")
      }

    } else {

      cat("\n")
      cat("The overall KMO value for your data is not available.")
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

  }

  if("EKC" %in% criteria){

    cat("Empirical Kaiser criterion: ", crayon::bold(n_fac["nfac_EKC"]),
        sep = "")
    cat("\n")

  }

  if("HULL" %in% criteria){

    if("CAF" %in% gof){
    cat("HULL method with CAF: ", crayon::bold(n_fac["nfac_HULL_CAF"]),
        sep = "")
    cat("\n")
    }
    if("CFI" %in% gof){
    cat("HULL method with CFI: ", crayon::bold(n_fac["nfac_HULL_CFI"]),
        sep = "")
    cat("\n")
    }
    if("RMSEA" %in% gof){
    cat("HULL method with RMSEA: ", crayon::bold(n_fac["nfac_HULL_RMSEA"]),
        sep = "")
    cat("\n")
    }

  }

  if("KGC" %in% criteria){

    if("PCA" %in% eigen_type_KGC_PA){
    cat("Kaiser-Guttman criterion with PCA: ",
        crayon::bold(n_fac["nfac_KGC_PCA"]), sep = "")
    cat("\n")
    }
    if("SMC" %in% eigen_type_KGC_PA){
    cat("Kaiser-Guttman criterion with SMC: ",
        crayon::bold(n_fac["nfac_KGC_SMC"]), sep = "")
    cat("\n")
    }
    if("EFA" %in% eigen_type_KGC_PA){
    cat("Kaiser-Guttman criterion with EFA: ",
        crayon::bold(n_fac["nfac_KGC_EFA"]), sep = "")
    cat("\n")
    }

  }

  if("PARALLEL" %in% criteria){

    if("PCA" %in% eigen_type_KGC_PA){
    cat("Parallel analysis with PCA: ", crayon::bold(n_fac["nfac_PA_PCA"]),
        sep = "")
    cat("\n")
    }
    if("SMC" %in% eigen_type_KGC_PA){
    cat("Parallel analysis with SMC: ", crayon::bold(n_fac["nfac_PA_SMC"]),
        sep = "")
    cat("\n")
    }
    if("EFA" %in% eigen_type_KGC_PA){
    cat("Parallel analysis with EFA: ", crayon::bold(n_fac["nfac_PA_EFA"]),
        sep = "")
    cat("\n")
    }

  }

  if("SMT" %in% criteria){

    cat("Sequential \U1D712\U00B2 model tests: ",
        crayon::bold(n_fac["nfac_SMT_chi"]), sep = "")
    cat("\n")
    cat("Lower bound of RMSEA 90% confidence interval: ",
        crayon::bold(n_fac["nfac_RMSEA"]), sep = "")
    cat("\n")
    cat("Akaike Information Criterion: ", crayon::bold(n_fac["nfac_AIC"]),
        sep = "")
    cat("\n")

  }

}
