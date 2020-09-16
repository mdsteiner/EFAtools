#' Print function for N_FACTORS objects
#'
#' @param x a list of class N_FACTORS. Output from \link{N_FACTORS} function.
#' @param ... Further arguments for print.
#'
#' @export
#' @method print N_FACTORS
#'
#' @examples
#' \donttest{
#' # All criteria except "CD", with correlation matrix and fit method "ML"
#' # (where needed)
#' N_FACTORS(test_models$baseline$cormat, criteria = c("EKC", "HULL", "KGC",
#'           "PARALLEL", "SCREE", "SMT"), N = 500, method = "ML")
#'}
print.N_FACTORS <- function(x, ...){

  suitability <- x$settings$suitability
  criteria <- x$settings$criteria
  n_fac <- x$n_factors
  gof <- x$settings$gof
  eigen_type_other <- x$settings$eigen_type_other
  kmo_out <- x$outputs$kmo_out
  KMO <- kmo_out$KMO

  if(isTRUE(suitability)){

    cat("\n")
    cat(cli::rule(left = crayon::bold("Tests for the suitability of the data for factor analysis"), col = "blue"))
    cat("\n")

    cat("\n")
    cat(crayon::blue$bold("Bartlett's test of sphericity"))
    cat("\n")

    print(x$output$bart_out)
    cat("\n")

    cat(crayon::blue$bold("Kaiser-Meyer-Olkin criterion (KMO)"))
    cat("\n")

    if(!is.na(KMO) && !is.null(KMO)){

      if(KMO >= .9){
        symb <- crayon::green$bold(cli::symbol$tick)
        label <- crayon::green$bold("marvellous")
      } else if(KMO >= .8){
        symb <- crayon::green$bold(cli::symbol$tick)
        label <- crayon::green$bold("meritorious")
      } else if(KMO >= .7){
        symb <- crayon::green$bold(cli::symbol$tick)
        label <- crayon::green$bold("middling")
      } else if(KMO >= .6){
        symb <- crayon::yellow$bold("!")
        label <- crayon::yellow$bold("mediocre")
      } else if (KMO >= .5){
        symb <- crayon::red$bold(cli::symbol$cross)
        label <- crayon::red$bold("miserable")
      } else {
        symb <- crayon::red$bold(cli::symbol$cross)
        label <- crayon::red$bold("unacceptable")
      }

      cat("\n")
      cat(symb, " The overall KMO value for your data is ", label,
          " with ", crayon::bold(round(KMO, 3)), ".", sep = "")
      cat("\n")


      if(KMO < .5){
        cat(crayon::bold(" "), "These data are not suitable for factor analysis.")
        cat("\n")
      } else if(KMO < .6){
        cat(crayon::bold(" "), "These data are hardly suitable for factor analysis.")
        cat("\n")
      } else {
        cat(crayon::bold(" "), "These data are probably suitable for factor analysis.")
        cat("\n")
      }

    } else {

      cat("\n")
      cat("The overall KMO value for your data is not available.")
      cat("\n")

    }

  }

  cat("\n")
  cat(cli::rule(crayon::bold("Number of factors suggested by the different factor",
  "retention criteria"), col = "blue"))
  cat("\n")
  cat("\n")

  if("CD" %in% criteria){

  cat(crayon::blue(cli::symbol$circle_dotted, "Comparison data: "),
      crayon::bold(n_fac["nfac_CD"]), sep = "")
  cat("\n")

  }

  if("EKC" %in% criteria){

    cat(crayon::blue(cli::symbol$circle_dotted, "Empirical Kaiser criterion: "),
        crayon::bold(n_fac["nfac_EKC"]),
        sep = "")
    cat("\n")

  }

  if("HULL" %in% criteria){

    if("CAF" %in% gof){
    cat(crayon::blue(cli::symbol$circle_dotted, "Hull method with CAF: "),
        crayon::bold(n_fac["nfac_HULL_CAF"]),
        sep = "")
    cat("\n")
    }
    if("CFI" %in% gof){
    cat(crayon::blue(cli::symbol$circle_dotted, "Hull method with CFI: "),
        crayon::bold(n_fac["nfac_HULL_CFI"]),
        sep = "")
    cat("\n")
    }
    if("RMSEA" %in% gof){
    cat(crayon::blue(cli::symbol$circle_dotted, "Hull method with RMSEA: "),
        crayon::bold(n_fac["nfac_HULL_RMSEA"]),
        sep = "")
    cat("\n")
    }

  }

  if("KGC" %in% criteria){

    if("PCA" %in% eigen_type_other){
    cat(crayon::blue(cli::symbol$circle_dotted,
                     "Kaiser-Guttman criterion with PCA: "),
        crayon::bold(n_fac["nfac_KGC_PCA"]), sep = "")
    cat("\n")
    }
    if("SMC" %in% eigen_type_other){
    cat(crayon::blue(cli::symbol$circle_dotted,
                     "Kaiser-Guttman criterion with SMC: "),
        crayon::bold(n_fac["nfac_KGC_SMC"]), sep = "")
    cat("\n")
    }
    if("EFA" %in% eigen_type_other){
    cat(crayon::blue(cli::symbol$circle_dotted,
                     "Kaiser-Guttman criterion with EFA: "),
        crayon::bold(n_fac["nfac_KGC_EFA"]), sep = "")
    cat("\n")
    }

  }

  if("PARALLEL" %in% criteria){

    if("PCA" %in% eigen_type_other){
    cat(crayon::blue(cli::symbol$circle_dotted,
                     "Parallel analysis with PCA: "),
        crayon::bold(n_fac["nfac_PA_PCA"]),
        sep = "")
    cat("\n")
    }
    if("SMC" %in% eigen_type_other){
    cat(crayon::blue(cli::symbol$circle_dotted,
                     "Parallel analysis with SMC: "),
        crayon::bold(n_fac["nfac_PA_SMC"]),
        sep = "")
    cat("\n")
    }
    if("EFA" %in% eigen_type_other){
    cat(crayon::blue(cli::symbol$circle_dotted,
                     "Parallel analysis with EFA: "),
        crayon::bold(n_fac["nfac_PA_EFA"]),
        sep = "")
    cat("\n")
    }

  }

  if("SMT" %in% criteria){

    cat(crayon::blue(cli::symbol$circle_dotted,
                     "Sequential \U1D712\U00B2 model tests: "),
        crayon::bold(n_fac["nfac_SMT_chi"]), sep = "")
    cat("\n")
    cat(crayon::blue(cli::symbol$circle_dotted,
                     "Lower bound of RMSEA 90% confidence interval: "),
        crayon::bold(n_fac["nfac_RMSEA"]), sep = "")
    cat("\n")
    cat(crayon::blue(cli::symbol$circle_dotted,
                     "Akaike Information Criterion: "),
        crayon::bold(n_fac["nfac_AIC"]),
        sep = "")
    cat("\n")

  }

  if("SCREE" %in% criteria){

    graphics::plot(x$output$scree_out)

  }

}
