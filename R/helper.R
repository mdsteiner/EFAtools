#' Extract a list object by its name
#'
#' @author Andreas Soteriades
#'
#' Consider a list of named sub-lists. This function extracts, for each sub-list,
#' the sub-list element that is specified by the user. This function is useful
#' for extracting results from \code{\link{EFA}} for each permutation run in
#' \code{\link{EFA_POOLED}}.
#'
#' @param alist A list of sub-lists, typically a list of \eqn{m} objects of class
#' \code{"EFA"}, where \eqn{m} is the number of imputations passed to
#' \code{\link{EFA_POOLED}}.
#' @param object String of length 1. The name of the object to extract e.g.
#' \code{"h2"} or \code{"vars_accounted"}.
#'
#' @return A list of length \eqn{m}, with each element containing the extracted
#' \code{object} for the \eqn{k}th element (\eqn{k = 1,..., m}).
.extract_list_object <- function(alist, object) {
  lapply(
    alist,
    function(x) {
      x[[object]]
    }
  )
}

#' Calculate statistics for a list of matrices
#'
#' @author Andreas Soteriades
#'
#' Given a list of matrices, this function calculates user-supplied statistics
#' (e.g. mean, median) over the matrices. This function is useful for averaging
#' results from \code{\link{EFA}} (e.g. loadings) over all permutation runs in
#' \code{\link{EFA_POOLED}}.
#'
#' @param alist A list of sub-lists, typically a list of \eqn{m} matrices from
#' \code{"EFA"}, where \eqn{m} is the number of imputations passed to
#' \code{\link{EFA_POOLED}}.
#' @param stat A function, e.g. \code{mean} or \code{sd}.
#'
#' @return A matrix with the aggregated results
#'
.stat_over_list <- function(alist, stat) {
  # In what follows:
  # n: length(alist)
  # n_row: nrow(alist[[i]]), where i is any list element
  # n_col: ncol(alist[[i]]), where i is any list element

  # Convert list to array
  list_to_array <- simplify2array(alist)

  # Apply stat over the array. The 1:2 will apply stat on an array with n_col
  # elements (run apply(list_to_array, 1:2, function(x) x) to see said array).
  # Each element is a n x n_row matrix (run
  # dim(apply(list_to_array, 1:2, function(x) x)) to confirm said array's
  # dimensions are n x n_row x n_col, i.e. n_col matrices with dimensions
  # n x n_row). In each matrix, the i-th row is the i-th column of the matrix in
  # the j-th element of list_to_array.
  # The apply() operation takes the column sum in each element of this new array
  # that has n_col elements, and returns the means in a n_row x n_col matrix.
  apply(
    list_to_array,
    1:2,
    stat
  )
}

#' Covert a \code{"LOADINGS"} table to matrix or a matrix to \code{"LOADINGS"}
#'
#' @author Andreas Soteriades
#'
#' The loadings tables returned by \code{\link{EFA}} are of class
#' \code{"LOADINGS"}, which prevents applying functions on them. This function
#' allows to change their class to \code{"matrix"}, and to change back to
#' \code{"LOADINGS"} when done.
#'
#' @param x A table of class \code{"matrix"} or \code{"LOADINGS"}.
#' @param cl A string with the class to change the table to. Should be
#' \code{"LOADINGS"} or \code{"matrix"}.
#'
#' @return A table with the loadings, of class either \code{"LOADINGS"} or
#' \code{"matrix"}.
.change_class <- function(x, cl = 'matrix') {
  class(x) <- cl
  return(x)
}

#' Confidence intervals around mean
#'
#' @author Andreas Soteriades
#'
#' This function is used internally by \code{\link{EFA_POOLED}} to calculate
#' confidence intervals (CIs) around the pooled loadings and pooled interfactor
#' correlations.
#'
#' @details
#' The standard error (SE) to use in the creation of the CIs is calculated
#' according to Rubin's (1987) formula (eq. 11 in Hayes & Enders, 2023):
#'
#' \eqn{\sqrt{mean(SE^2) + var(\hat{\theta_{m}}) + var(\hat{\theta_{m}}) / M}}
#'
#' According to Hayes & Enders (2023) p. 42:
#'
#' \emph{[T]he first term under the radical represents the average squared
#' standard error (the within imputation sampling variance [...]), the second
#' term depends on the variance of the M parameter estimates around their
#' average (the between imputation variance [...]), and the final term
#' represents the squared standard error of the pooled estimate [...].
#' Conceptually, the first term estimates the sampling error of a complete-data
#' analysis, and the next two terms are essentially correction factors that
#' inflate the standard error to compensate for uncertainty due to the
#' imputations–that is, additional uncertainty (sampling variability) in the
#' parameter estimates caused by missing data.}
#'
#' Currently, it is not possible to calculate the first term, because
#' \code{\link{EFA}} does not calculate SEs for the loadings. Only the second
#' and third terms are used in the calculation of SE for the CIs.
#'
#' The CI is generally calculated as:
#'
#' CI = Point estimate ± Margin of error,
#'
#' where
#'
#' Margin of error = Critical value × SE of point estimate.
#'
#' To account for situations where the sample size is small, instead of using
#' z-values, the critical value is derived from the \emph{t} distribution with
#' \code{n - 1} degrees of freedom (Hazra, 2017).
#
#' @param means A matrix with the pooled loadings or pooled interfactor
#' correlations
#' @param sds A matrix with the standard deviations for the loadings or
#' interfactor correlations
#' @param p Numeric. One minus the confidence level of the CIs. Defaults to 0.05
#' for 95\% CIs.
#' @param n Numeric. Typically, the number of permutations \eqn{m}. See
#' \code{\link{EFA_POOLED}}.
#'
#' @return A list with the lower and upper CIs.
#'
#' @references
#' Hayes, T. & Enders, C. K. (2023). Maximum likelihood and multiple imputation
#' missing data handling: how they work, and how to make them work in practice.
#' In \emph{APA Handbook of Research Methods in Psychology, Second Edition} Vol.
#' 3. Data Analysis and Research Publication, H. Cooper (Editor-in-Chief).
#'
#' Hazra, A. (2017). Using the confidence interval confidently. \emph{Journal of
#' Thoracic Disease} 9(10), 4125--4130.
#'
#' Rubin, D. B. (1987). Multiple imputation for nonresponse in surveys. Wiley.
#' https://doi.org/10.1002/9780470316696
.calc_cis <- function(means, sds, p = 0.05, n) {

  # Calculate standard error
  rubin_term_1 <- 0 # Currently not calculated- see Details
  rubin_term_2 <- sds ^ 2 # The variances
  rubin_term_3 <- rubin_term_2 / n
  se <- sqrt(rubin_term_1 + rubin_term_2 + rubin_term_3)

  error <- stats::qt(1 - p / 2, n - 1) * se

  lower_ci <- means - error
  upper_ci <-  means + error

  cis <- list(
    lower = lower_ci,
    upper = upper_ci
  )

  return(cis)
}

#' Format numbers for print method
#'
#' Helper function used in the print method for class LOADINGS and SLLOADINGS.
#' Strips the 0 in front of the decimal point of a number if number < 1, only
#' keeps the first \code{digits} number of digits, and adds an empty space in
#' front of the number if the number is positive. This way all returned strings
#' (except for those > 1, which are exceptions in LOADINGS) have the same number
#' of characters.
#'
#' @param x numeric. Number to be formatted.
#' @param digits numeric. Number of digits after the comma to keep.
#' @param print_zero logical. Whether, if a number is between [-1, 1], the
#'  zero should be omitted or printed (default is FALSE, i.e. omit zeros).
#' @param pad logical. Whether, if a number starts with a 0 and the 0 is not printed
#'  a white-space should be added.
#'
#' @return A formated number
.numformat <- function(x, digits = 2, print_zero = FALSE,
                       pad = TRUE) {

  if (isFALSE(print_zero)) {

    ncode <- paste0("%.", digits, "f")
    x <- sub("^(-?)0.", "\\1.", sprintf(ncode, x))
    if (isTRUE(pad)) {
      x <- stringr::str_pad(x, digits + 2, "left")
    }


  } else {

    ncode <- paste0("%.", digits, "f")
    x <- sprintf(ncode, x)

    if (isTRUE(pad)) {
      x <- stringr::str_pad(x, digits + 3, "left")
    }

  }
  x

}


#' Compute explained variances from loadings
#'
#' From unrotated loadings compute the communalities and uniquenesses for total
#' variance. Compute explained variances per factor from rotated loadings (and
#' factor intercorrelations Phi if oblique rotation was used).
#'
#' @param L_unrot matrix. Unrotated factor loadings.
#' @param L_rot matrix. Rotated factor loadings.
#' @param Phi matrix. Factor intercorrelations. Provide only if oblique rotation
#'  is used.
#'
#' @return A matrix with sum of squared loadings, proportion explained variance
#'  from total variance per factor, same as previous but cumulative, Proportion
#'  of explained variance from total explained variance, and same as previous but
#'  cumulative.
.compute_vars <- function(L_unrot, L_rot, Phi = NULL) {

  if (is.null(Phi)) {
    # compute variance proportions
    if (ncol(L_rot) > 1) {
      vars <- colSums(L_rot^2)
    }
    else {
      vars <- sum(L_rot^2)
    }
  } else {
    # compute variance proportions
    vars <- diag(Phi %*% t(L_rot) %*% L_rot)
  }

  # Compute the explained variances. The code is based on the psych::fac() function
  # total variance (sum of communalities and uniquenesses)
  h2 <- diag(L_unrot %*% t(L_unrot))
  var_total <- sum(h2 + (1 - h2))
  vars_explained <- rbind(`SS loadings` = vars)
  vars_explained <- rbind(vars_explained, `Prop Tot Var` = vars / var_total)

  if (ncol(L_rot) > 1) {
    vars_explained <- rbind(vars_explained,
                            `Cum Prop Tot Var` = cumsum(vars / var_total))
    vars_explained <- rbind(vars_explained,
                            `Prop Comm Var` = vars / sum(vars))
    vars_explained <- rbind(vars_explained,
                            `Cum Prop Comm Var` = cumsum(vars / sum(vars)))
  }

  vars_explained
}

# varimax criterion for SPSS varimax implementation
.SV <- function(lambda) {

  n <- nrow(lambda)

  # the SPSS manual (ftp://public.dhe.ibm.com/software/analytics/spss/documentation/statistics/23.0/en/client/Manuals/IBM_SPSS_Statistics_Algorithms.pdf)
  # suggests the following formula:
  # sum(n*colSums(lambda**4) - colSums(lambda ** 2) ** 2) / n**2
  # however, the formula below produces results more in line with SPSS
  sum(n*colSums(abs(lambda)) - colSums(lambda ** 4) ** 2) / n**2

}

.get_compare_matrix <- function(x, digits = 3, r_red = .001, n_char = 10,
                                var_names = NULL, factor_names = NULL,
                                gof = FALSE) {

  # create factor names to display
  if (is.null(factor_names)) {
    if(is.null(colnames(x))){
      factor_names <- paste0("F", seq_len(ncol(x)))
    } else {
      factor_names <- colnames(x)
    }
  }

  # for equal spacing, fill the factor names such that they match the columns
  fn_nchar <- sapply(factor_names, nchar)
  factor_names[which(fn_nchar > digits + 2)] <- substr(
    factor_names[which(fn_nchar > digits + 2)] , 1, digits + 2)
  factor_names <- stringr::str_pad(factor_names, digits + 2, side = "both")

  if(gof == FALSE){

  if(is.null(var_names)) {
    if(is.null(rownames(x))){
      var_names <- paste0("V", seq_len(nrow(x)))
    } else {
    var_names <- rownames(x)
    }
  }

  max_char <- max(sapply(var_names, nchar))

  if (max_char > n_char) {
    vn_nchar <- sapply(var_names, nchar)
    var_names[which(vn_nchar > n_char)] <- substr(var_names[which(vn_nchar > n_char)],
                                              1, n_char)
    max_char <- n_char
  }

  var_names <- stringr::str_pad(var_names, max_char, side = "right")

  }

  n_col <- ncol(x)

  # create the string to paste using the crayon package
  temp <- apply(matrix(seq_len(nrow(x)), ncol = 1), 1,
                function(ind, x, cutoff, n_col, vn, digits){
                  i <- x[ind,]

                  tt <- crayon::blue(vn[ind])

                  for (kk in seq_len(n_col)) {
                    if (abs(i[kk]) <= cutoff) {
                      tt <- c(tt, .numformat(round(i[kk], digits = digits),
                                                          digits = digits,
                                             print_zero = TRUE))
                    } else {
                      tt <- c(tt,
                              crayon::red(.numformat(round(i[kk],
                                                                digits = digits),
                                                              digits = digits,
                                                              print_zero = TRUE)))
                    }
                  }
                  stringr::str_c(tt, collapse = "\t")
                }, cutoff = r_red, n_col = n_col, digits = digits, x = x,
                vn = var_names)

  factor_names <- stringr::str_c(factor_names, collapse = "\t")

  if(gof == TRUE){

    factor_names <- crayon::blue(stringr::str_c(factor_names))

  } else {

    factor_names <- crayon::blue(stringr::str_c( stringr::str_pad(" ", max_char),
                                                   "\t", factor_names))
  }


  temp <- stringr::str_c(temp, collapse = "\n")

  temp <- stringr::str_c(factor_names, "\n", temp)


  temp <- stringr::str_c(temp, "\n")

  # print the results to the console

  temp
}

.get_compare_vector <- function(x, digits = 3, r_red = .001) {

  temp_i <- NULL

  for (ii in seq_along(x)) {
    if (abs(x[ii]) > r_red) {
      temp_i <- c(temp_i, crayon::red(.numformat(round(x[ii], digits = digits),
                                                      digits = digits,
                                                      print_zero = TRUE)))
    } else {
      temp_i <- c(temp_i, .numformat(round(x[ii], digits = digits),
                                     digits = digits,
                                     print_zero = TRUE))
    }
  }

  for (ss in seq(1, length(x), 7)) {
    if (length(x) > ss + 6) {
      tt <- ss + 6
    } else {
      tt <- length(x)
    }
    if (ss == 1) {
      temp <- stringr::str_c(temp_i[ss:tt], collapse = "  ")
    } else {
      temp <- stringr::str_c(temp, "\n", stringr::str_c(temp_i[ss:tt],
                                                        collapse = "  "))
    }

  }

  temp <- stringr::str_c(temp, "\n")

  # print the results to the console
  return(temp)
}



.decimals <- function(x) {

  if ((is.null(dim(x)) && !(inherits(x, c("numeric", "integer")))) ||
      (!is.null(dim(x)) && !(inherits(x, c("matrix", "loadings", "LOADINGS",
                                          "SLLOADINGS"))))) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is of class '", class(x), "' but must be a numeric vector or matrix\n", sep = ""))
  }

  if (!is.null(dim(x))) {

    max(apply(x, 1:2, function(ll) {
      if (abs(ll - round(ll)) > .Machine$double.eps^0.5) {
        nchar(strsplit(sub('0+$', '', as.character(ll)), ".",
                       fixed = TRUE)[[1]][[2]])
      } else {
        return(0)
      }
    }))

  } else if (length(x) > 1) {

    max(sapply(x, function(ll) {
      if (abs(ll - round(ll)) > .Machine$double.eps^0.5) {
        nchar(strsplit(sub('0+$', '', as.character(ll)), ".", fixed = TRUE)[[1]][[2]])
      } else {
        return(0)
      }
    }))


  } else {

    if (abs(x - round(x)) > .Machine$double.eps^0.5) {
      nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
    } else {
      return(0)
    }

  }

}


.factor_congruence <- function(x, y, na.rm = TRUE, skip_checks = FALSE) {

  if (isFALSE(skip_checks)) {

    if (any(is.na(x) | any(is.na(y)))) {
      if (isTRUE(na.rm)) {
        warning(crayon::yellow$bold("!"), crayon::yellow(" Input contained missing values. Analysis is performed on complete cases.\n"))
        if (any(is.na(x))) {
          xc <- x[stats::complete.cases(x), ]
          y <- y[stats::complete.cases(x), ]
          x <- xc
        }
        if (any(is.na(y))) {
          yc <- y[stats::complete.cases(y), ]
          x <- x[stats::complete.cases(y), ]
          y <- yc
        }
      } else {
        warning(crayon::yellow$bold("!"), crayon::yellow(" Input contained missing values. Check your data or rerun with na.rm = TRUE.\n"))
      }
    }

  }


  nx <- dim(x)[2]
  ny <- dim(y)[2]
  cross <- t(y) %*% x
  sumsx <- sqrt(1/diag(t(x) %*% x))
  sumsy <- sqrt(1/diag(t(y) %*% y))
  result <- sumsy * (cross * rep(sumsx, each = ny))
  return(t(result))

}


.onUnload <- function (libpath) {
  library.dynam.unload("EFAtools", libpath)
}

.rmsr <- function(residuals, upper = TRUE) {
  E <- as.matrix(residuals)
  if (isTRUE(upper)) {
    vals <- E[upper.tri(E, diag = FALSE)]
  } else {
    vals <- E[row(E) != col(E)]
  }
  sqrt(mean(vals^2, na.rm = TRUE))
}

.gof <- function(L, # The loading/ pattern matrix
                 R, # The correlation matrix
                 N, # The number of cases
                 method, # The estimation method
                 Fm) { # Minimized error
  m <- nrow(L)
  q <- ncol(L)

  # dfs
  df <- ((m - q)**2 - (m + q)) / 2

  ### compute CAF
  delta_hat <- R - (L %*% t(L))
  diag(delta_hat) <- 1

  # try statement needed to fix CRAN issue MKL (R-devel on x86_64 Fedora 34 Linux
  # with alternative BLAS/LAPACK implementations) where delta_hat cannot be inverted
  # for R = test_models$baseline$cormat
  delta_hat_KMO <- try(KMO(delta_hat)$KMO, silent = TRUE)
  if (inherits(delta_hat_KMO, "try-error")) {
    CAF <- 0
    warning(crayon::yellow$bold("!"), crayon::yellow(" Problems calculating CAF, CAF set to 0 (worst value). Inspect results carefully.\n"))
  } else {
    CAF <- 1 - delta_hat_KMO
  }

  ### compute RMSR
  RMSR <- .rmsr(delta_hat)

  if (method != "PAF" && !is.na(N) && df >=0) {

    ### compute CFI

    # null model
    chi_null <- sum(R[upper.tri(R)] ^ 2) * (N - 1)
    df_null <- (m**2 - m) / 2
    delta_hat_null <- chi_null - df_null
    p_null <- stats::pchisq(chi_null, df_null, lower.tail = F)

    # current model
    # formula 12.1 from Kline 2016; Principles and practices of...
    # This should also hold for ULS solutions -> Bentler & Bonett (1980)
    # Significance Tests  and Goodness  of Fit in the Analysis of Covariance
    # Structures. Psychological Bulletin, 88(3), 588-606
    chi <- Fm * (N - 1)
    p_chi <- stats::pchisq(chi, df, lower.tail = F)
    delta_hat_m <- chi - df
    CFI <- (delta_hat_null - delta_hat_m) / delta_hat_null
    if (CFI > 1 || df == 0) CFI <- 1
    if (CFI < 0) CFI <- 0


    ### compute RMSEA, incl. 90% confidence intervals if df are not 0
    if(df != 0){

      # formula 12.6 from Kline 2016; Principles and practices of...
      RMSEA <- sqrt(max(0, chi - df) / (df * (N - 1)))

    p_chi_fun <- function(x, val, df, goal){goal - stats::pchisq(val, df, ncp = x)}

    if (stats::pchisq(chi, df = df, ncp = 0) >= .95) {
      lambda_l <- stats::uniroot(f = p_chi_fun, interval = c(1e-10, 10000), val = chi,
                          df = df, goal = .95, extendInt = "upX",
                          maxiter = 100L)$root
    } else {
      lambda_l <- 0
    }

    if (stats::pchisq(chi, df = df, ncp = 0) >= .05) {
      lambda_u <- stats::uniroot(f = p_chi_fun, interval = c(1e-10, 10000),
                          val = chi, df = df, goal = .05,
                          extendInt = "upX", maxiter = 100L)$root
    }
    else {
      lambda_u <- 0
    }

    RMSEA_LB <- sqrt(lambda_l / (df * (N - 1)))
    RMSEA_UB <- sqrt(lambda_u / (df * (N - 1)))

    if(RMSEA > 1) RMSEA <- 1
    if(RMSEA_LB > 1) RMSEA_LB <- 1
    if(RMSEA_UB > 1) RMSEA_UB <- 1

    } else {

      RMSEA <- 0
      RMSEA_LB <- 0
      RMSEA_UB <- 0

    }

    ### compute AIC and BIC based on chi square
    AIC <- chi - 2 * df
    BIC <- chi - log(N) * df

  } else {
    chi <- NA
    p_chi <- NA
    CFI <- NA
    RMSEA <- NA
    RMSEA_LB <- NA
    RMSEA_UB <- NA
    AIC <- NA
    BIC <- NA
    chi_null <- NA
    df_null <- NA
    p_null <- NA
  }

  out <- list(
    chi = chi,
    df = df,
    p_chi = p_chi,
    CAF = CAF,
    RMSR = RMSR,
    CFI = CFI,
    RMSEA = RMSEA,
    RMSEA_LB = RMSEA_LB,
    RMSEA_UB = RMSEA_UB,
    AIC = AIC,
    BIC = BIC,
    Fm = Fm,
    chi_null = chi_null,
    df_null = df_null,
    p_null = p_null
  )

}


# Checks if x is a correlation matrix
.is_cormat <- function(x){

  if(nrow(x) == ncol(x) &&
     all(x >= (-1 + .Machine$double.eps * 100), na.rm = TRUE) &&
     all(x <= (1 + .Machine$double.eps * 100), na.rm = TRUE)){

    if (round(sum(diag(x), na.rm = TRUE)) == nrow(x) && isSymmetric(unclass(unname(x)))) {

      if (any(is.na(x))) {

        stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(' "x" is likely a correlation matrix but contains missing values. Please check the entered data.\n'))

      }

      TRUE

    } else {

      FALSE

    }


  } else {

    FALSE

  }

}

# Create a string (used to print settings in print functions, to use in cat()
# with sep = "")
.settings_string <- function(x){

  n <- length(x)

if(n == 1){

  c(crayon::bold(x))

} else if (n == 2){

  c(crayon::bold(x[1]), " and ", crayon::bold(x[2]))

} else if (n > 2){

  c(paste(crayon::bold(x[seq_len(n-1)]), collapse = ", "), ", and ",
    crayon::bold(x[n]))

}

}

.det_max_factors <- function(m) {
  q <- floor((2*m + 1 - sqrt(8 * m + 9)) / 2)
  if(q < 0) q <- 0
  return(q)
}

# for progress bar in N_FACTORS
.show_progress <- function(x, what, done = FALSE) {

  cat("\r", rep(" ", ifelse(options("width") > 30, options("width"), 30)))
  to <- length(x)
  if (isFALSE(done)) {
    curr <- which(x == what)

    #cat("\r", paste0(curr, "/", to, ":"), "Running", what)
    cat("\r", rep(cli::symbol$circle_filled, curr - 1),
        "\U1F3C3", rep(cli::symbol$circle, to - (curr)),
        "Running", what)
  } else {
    cat("\r", rep(cli::symbol$circle_filled, to),
        "Done!\n")
  }

}


# for progress bar in EFA_AVERAGE
.show_av_progress <- function(emoji, what, done = FALSE) {

  cat("\r", rep(" ", ifelse(options("width") > 30, options("width"), 30)))
  if (isFALSE(done)) {
    #cat("\r", paste0(curr, "/", to, ":"), "Running", what)
    cat("\r", emoji, what)
  } else {
    cat("\r", "Done!\n")
  }

}


### extract data from efa_list
.extract_data <- function(efa_list, R, n_factors, n_efa, rotation, salience_threshold) {


  L <- array(NA_real_, c(ncol(R), n_factors, n_efa))
  L_corres <- array(NA, c(ncol(R), n_factors, n_efa))
  h2 <- matrix(NA_real_, nrow = n_efa, ncol = ncol(R))
  if (n_factors > 1) {
    vars_accounted <- array(NA_real_, c(3, n_factors, n_efa))
  } else {
    vars_accounted <- array(NA_real_, c(2, n_factors, n_efa))
  }



  if (any(rotation %in% c("promax", "oblimin", "quartimin", "simplimax",
                          "bentlerQ", "geominQ", "bifactorQ", "oblique"))) {
    extract_phi <- TRUE
    phi <- array(NA_real_, c(n_factors, n_factors, n_efa))
  } else {
    extract_phi <- FALSE
    phi <- NA
  }

  converged <- rep(NA, n_efa)
  errors <- rep(FALSE, n_efa)
  error_m <- rep(NA_character_, n_efa)
  heywood <- rep(NA, n_efa)
  admissible <- rep(NA, n_efa)
  aic <- rep(NA_real_, n_efa)
  bic <- rep(NA_real_, n_efa)
  chisq <- rep(NA_real_, n_efa)
  p_chi <- rep(NA_real_, n_efa)
  caf <- rep(NA_real_, n_efa)
  rmsea <- rep(NA_real_, n_efa)
  cfi <- rep(NA_real_, n_efa)

  if (all(rotation == "none") || n_factors == 1) {
    load_ind <- "unrot_loadings"
    var_ind <- "vars_accounted"
  } else {
    load_ind <- "rot_loadings"
    var_ind <- "vars_accounted_rot"
  }

    for (row_i in seq_len(n_efa)) {

      efa_temp <- efa_list[[row_i]]

      if (inherits(efa_temp, "try-error")) {

        errors[row_i] <- TRUE
        error_m[row_i] <- efa_temp[[1]]

      } else {
        converged[row_i] <- efa_temp$convergence

        if (efa_temp$convergence == 0) {
          has_heywood <- any(efa_temp$h2 >= 1 + .Machine$double.eps)
          heywood[row_i] <- has_heywood

          if (!has_heywood) {

            aic[row_i] <- efa_temp$fit_indices$AIC
            bic[row_i] <- efa_temp$fit_indices$BIC
            chisq[row_i] <- efa_temp$fit_indices$chi
            p_chi[row_i] <- efa_temp$fit_indices$p_chi
            caf[row_i] <- efa_temp$fit_indices$CAF
            rmsea[row_i] <- efa_temp$fit_indices$RMSEA
            cfi[row_i] <- efa_temp$fit_indices$CFI

            h2[row_i, ] <- efa_temp$h2
            L[,, row_i] <- efa_temp[[load_ind]]
            if (n_factors > 1) {
              vars_accounted[,, row_i] <- efa_temp[[var_ind]][c(1, 2, 4),]
            } else {
              vars_accounted[,, row_i] <- efa_temp[[var_ind]]
            }

            temp_corres <- abs(efa_temp[[load_ind]]) >= salience_threshold
            L_corres[,, row_i] <-temp_corres

          }


          admissible[row_i] <- ifelse(has_heywood || any(colSums(temp_corres) < 2),
                                      FALSE, TRUE)


          if (isTRUE(extract_phi)) {
            phi[,, row_i] <- efa_temp$Phi
          }
        }
      }
    }

  # remove data from nonconverged EFAs
  excl <- which(converged != 0 | errors | heywood)
  if (length(excl) > 0) {
    L <- L[,, -excl, drop = FALSE]
    L_corres <- L_corres[,, -excl, drop = FALSE]
    vars_accounted <- vars_accounted[,, -excl, drop = FALSE]
    if (isTRUE(extract_phi)) {
      phi <- phi[,, -excl, drop = FALSE]
    }
  }

  out <- list(
    L = L,
    L_corres = L_corres,
    phi = phi,
    extract_phi = extract_phi,
    h2 = h2,
    vars_accounted = vars_accounted,
    for_grid = data.frame(
      errors = errors,
      error_m = error_m,
      converged = converged,
      heywood = heywood,
      admissible = admissible,
      chisq = chisq,
      p_chi = p_chi,
      caf = caf,
      cfi = cfi,
      rmsea = rmsea,
      aic = aic,
      bic= bic
    )
  )

  return(out)

}

### average arrays
.average_values <- function(vars_accounted, L, L_corres, h2, phi, extract_phi,
                              averaging, trim, for_grid, df, ind_names) {

  if (averaging == "mean") {

    if (trim == 0) {
      # faster, but only works without trimming
      L_av <- rowMeans(L, na.rm = TRUE, dims = 2)
      h2_av <- colMeans(h2, na.rm = TRUE)
      fit_av <- colMeans(for_grid, na.rm = TRUE)
      vars_accounted_av <- rowMeans(vars_accounted, na.rm = TRUE, dims = 2)


      if (isTRUE(extract_phi)) {
        phi_av <- rowMeans(phi, na.rm = TRUE, dims = 2)
      }
    } else {
      L_av <- apply(L, 1:2, mean, na.rm = TRUE, trim = trim)
      h2_av <- apply(h2, 2, mean, na.rm = TRUE, trim = trim)
      fit_av <- apply(for_grid, 2, mean, na.rm = TRUE, trim = trim)
      vars_accounted_av <- apply(vars_accounted, 1:2, mean, na.rm = TRUE,
                                  trim = trim)

      if (isTRUE(extract_phi)) {
        phi_av <- apply(phi, 1:2, mean, na.rm = TRUE, trim = trim)
      }
    }

  } else if (averaging == "median") {
    L_av <- apply(L, 1:2, stats::median, na.rm = TRUE)
    h2_av <- apply(h2, 2, stats::median, na.rm = TRUE)
    fit_av <- apply(for_grid, 2, stats::median, na.rm = TRUE)
    vars_accounted_av <- apply(vars_accounted, 1:2, stats::median, na.rm = TRUE)
    if (isTRUE(extract_phi)) {
      phi_av <- apply(phi, 1:2, stats::median, na.rm = TRUE)
    }
  }


  nf <- ncol(L_av)
  f_names <- paste0("F", 1:nf)

  L_corres_av <- rowMeans(L_corres, na.rm = TRUE, dims = 2)
  row.names(L_corres_av) <- ind_names
  colnames(L_corres_av) <- f_names

  L_min <- apply(L, 1:2, min, na.rm = TRUE)
  L_max <- apply(L, 1:2, max, na.rm = TRUE)
  L_range <- L_max - L_min
  L_sd <- apply(L, 1:2, stats::sd, na.rm = TRUE)
  rownames(L_av) <- ind_names
  colnames(L_av) <- f_names
  class(L_av) <- "LOADINGS"
  rownames(L_min) <- ind_names
  colnames(L_min) <- f_names
  class(L_min) <- "LOADINGS"
  rownames(L_max) <- ind_names
  colnames(L_max) <- f_names
  class(L_max) <- "LOADINGS"
  rownames(L_range) <- ind_names
  colnames(L_range) <- f_names
  rownames(L_sd) <- ind_names
  colnames(L_sd) <- f_names



  vars_accounted_min <- apply(vars_accounted, 1:2, min, na.rm = TRUE)
  vars_accounted_max <- apply(vars_accounted, 1:2, max, na.rm = TRUE)
  vars_accounted_range <- vars_accounted_max - vars_accounted_min
  vars_accounted_sd <- apply(vars_accounted, 1:2, stats::sd, na.rm = TRUE)


  if (nrow(vars_accounted_av) == 2) {
    var_names <- c("SS loadings", "Prop Tot Var")
  } else {
    var_names <- c("SS loadings", "Prop Tot Var", "Prop Comm Var")
  }
  rownames(vars_accounted_av) <- var_names
  colnames(vars_accounted_av) <- f_names
  rownames(vars_accounted_min) <- var_names
  colnames(vars_accounted_min) <- f_names
  rownames(vars_accounted_max) <- var_names
  colnames(vars_accounted_max) <- f_names
  rownames(vars_accounted_range) <- var_names
  colnames(vars_accounted_range) <- f_names
  rownames(vars_accounted_sd) <- var_names
  colnames(vars_accounted_sd) <- f_names


  h2_min <- apply(h2, 2, min, na.rm = TRUE)
  h2_max <- apply(h2, 2, max, na.rm = TRUE)
  h2_range <- h2_max - h2_min
  h2_sd <- apply(h2, 2, stats::sd, na.rm = TRUE)
  names(h2_av) <- ind_names
  names(h2_min) <- ind_names
  names(h2_max) <- ind_names
  names(h2_range) <- ind_names
  names(h2_sd) <- ind_names


  fit_min <- apply(for_grid, 2, min, na.rm = TRUE)
  fit_max <- apply(for_grid, 2, max, na.rm = TRUE)
  fit_range <- fit_max - fit_min
  fit_sd <- apply(for_grid, 2, stats::sd, na.rm = TRUE)

  fit_av[is.nan(fit_av)] <- NA
  fit_min[is.infinite(fit_min)] <- NA
  fit_max[is.infinite(fit_max)] <- NA
  fit_range[is.infinite(fit_range)] <- NA

  fit_indices <- data.frame(
    index = c(names(fit_av), "df"),
    average = c(fit_av, df),
    sd = c(fit_sd, df),
    range = c(fit_range, df),
    min = c(fit_min, df),
    max = c(fit_max, df),
    stringsAsFactors = FALSE
  )

  if (isTRUE(extract_phi)) {
    phi_min <- apply(phi, 1:2, min, na.rm = TRUE)
    phi_max <- apply(phi, 1:2, max, na.rm = TRUE)
    phi_range <- phi_max - phi_min
    phi_sd <- apply(phi, 1:2, stats::sd, na.rm = TRUE)
    colnames(phi_av) <- paste0("F", 1:nf)
    rownames(phi_av) <- paste0("F", 1:nf)
    colnames(phi_min) <- paste0("F", 1:nf)
    rownames(phi_min) <- paste0("F", 1:nf)
    colnames(phi_max) <- paste0("F", 1:nf)
    rownames(phi_max) <- paste0("F", 1:nf)
    colnames(phi_range) <- paste0("F", 1:nf)
    rownames(phi_range) <- paste0("F", 1:nf)
    colnames(phi_sd) <- paste0("F", 1:nf)
    rownames(phi_sd) <- paste0("F", 1:nf)
  }


  if (isTRUE(extract_phi)) {
    phi_list <- list(
      average = phi_av,
      sd = phi_sd,
      min = phi_min,
      max = phi_max,
      range = phi_range
    )
  } else {
    phi_list <- NA
  }

  out <- list(
    h2 = list(
      average = h2_av,
      sd = h2_sd,
      min = h2_min,
      max = h2_max,
      range = h2_range
    ),
    loadings = list(
      average = L_av,
      sd = L_sd,
      min = L_min,
      max = L_max,
      range = L_range
    ),
    phi = phi_list,
    vars_accounted = list(
      average = vars_accounted_av,
      sd = vars_accounted_sd,
      min = vars_accounted_min,
      max = vars_accounted_max,
      range = vars_accounted_range
    ),
    ind_fac_corres = L_corres_av,
    fit_indices = fit_indices)

  return(out)
}


### reorder arrays according to factor congruence
.array_reorder <- function(vars_accounted, L, L_corres, phi, extract_phi, n_factors) {

  if (dim(L)[3] > 1) {
  	L1 <- L[,, 1]
    for (efa_i in 2:dim(L)[3]) {

      Ln <- L[,, efa_i]

      # reorder factors according to tuckers congruence coefficient
      # get Tucker's congruence coefficients
      congruence <- .factor_congruence(L1, Ln, skip_checks = TRUE)

      # factor order for Ln
      factor_order <- apply(abs(congruence), 1, which.max)

      # reorder
      Ln <- Ln[, factor_order]
      L_corres[,, efa_i] <- L_corres[,, efa_i][, factor_order]
      vars_accounted[,, efa_i] <- vars_accounted[,, efa_i][, factor_order]

      # get signs
      factor_sign <- diag(sign(diag(crossprod(L1, Ln))))

      # switch signs where necessary
      L[,, efa_i] <- Ln %*% factor_sign

      if (isTRUE(extract_phi)) {
        phi[,, efa_i] <- factor_sign %*% phi[,, efa_i][factor_order, factor_order] %*% factor_sign

      }

    }
  }


  return(list(L=L, L_corres = L_corres, phi = phi, vars_accounted = vars_accounted))

}

### create grid for oblique rotations in EFA_AVERAGE
.oblq_grid <- function(method, init_comm, criterion, criterion_type,
                       abs_eigen, start_method, rotation, k_promax, normalize, P_type,
                       precision, varimax_type, k_simplimax){

  g_list <- list()

  if ("promax" %in% rotation) {

    g_list[["prmx"]] <- expand.grid(method = method, init_comm = init_comm,
                                    criterion = criterion, criterion_type = criterion_type,
                                    abs_eigen = abs_eigen, start_method = start_method,
                                    rotation = "promax",
                                    k_promax = k_promax, normalize = normalize, P_type = P_type,
                                    precision = precision, varimax_type = varimax_type,
                                    k_simplimax = NA, stringsAsFactors = FALSE)

  }

  if ("simplimax" %in% rotation) {

    g_list[["smplmx"]] <- expand.grid(method = method, init_comm = init_comm,
                                      criterion = criterion, criterion_type = criterion_type,
                                      abs_eigen = abs_eigen, start_method = start_method,
                                      rotation = "simplimax",
                                      k_promax = NA, normalize = normalize, P_type = NA,
                                      precision = precision, varimax_type = NA,
                                      k_simplimax = k_simplimax, stringsAsFactors = FALSE)

  }

  rotation_temp <- rotation[!(rotation %in% c("promax", "simplimax"))]

  if (length(rotation_temp) > 0) {
    g_list[["oblq"]] <- expand.grid(method = method, init_comm = init_comm,
                                    criterion = criterion, criterion_type = criterion_type,
                                    abs_eigen = abs_eigen, start_method = start_method,
                                    rotation = rotation_temp,
                                    k_promax = NA, normalize = normalize, P_type = NA,
                                    precision = precision, varimax_type = NA,
                                    k_simplimax = NA, stringsAsFactors = FALSE)
  }

  return(do.call(rbind, g_list))

}

### create grid for orthogonal rotations in EFA_AVERAGE
.orth_grid <- function(method, init_comm, criterion, criterion_type,
                       abs_eigen, start_method, rotation, normalize,
                       precision, varimax_type){

  g_list <- list()

  if ("varimax" %in% rotation) {

    g_list[["vrmx"]] <- expand.grid(method = method, init_comm = init_comm,
                                    criterion = criterion, criterion_type = criterion_type,
                                    abs_eigen = abs_eigen, start_method = start_method,
                                    rotation = "varimax",
                                    k_promax = NA, normalize = normalize, P_type = NA,
                                    precision = precision, varimax_type = varimax_type,
                                    k_simplimax = NA, stringsAsFactors = FALSE)

  }

  rotation_temp <- rotation[!(rotation %in% c("varimax"))]

  if (length(rotation_temp) > 0) {
    g_list[["orth"]] <- expand.grid(method = method, init_comm = init_comm,
                                    criterion = criterion, criterion_type = criterion_type,
                                    abs_eigen = abs_eigen, start_method = start_method,
                                    rotation = rotation_temp,
                                    k_promax = NA, normalize = normalize, P_type = NA,
                                    precision = precision, varimax_type = NA,
                                    k_simplimax = NA, stringsAsFactors = FALSE)
  }

  return(do.call(rbind, g_list))

}

.type_grid <- function(method, init_comm, criterion, criterion_type,
                       abs_eigen, start_method, rotation, k_promax, normalize,
                       P_type, precision, varimax_type, k_simplimax) {

  t_grid_list <- list()
  if ("none" %in% rotation) {
    if (length(rotation) == 1) {

      t_grid_list[["nn"]] <- expand.grid(method = method, init_comm = init_comm,
                                         criterion = criterion, criterion_type = criterion_type,
                                         abs_eigen = abs_eigen, start_method = start_method,
                                         rotation = "none",
                                         k_promax = NA, normalize = NA, P_type = NA,
                                         precision = NA, varimax_type = NA,
                                         k_simplimax = NA, stringsAsFactors = FALSE)

    } else {

      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" rotation = 'none' is used but rotation is of length > 1. Can only average EFAs with rotations of the same type ('none', 'orthogonal', or 'oblique').\n"))

    }
  } else if ("oblique" %in% rotation) {
    if (length(rotation) == 1) {

      t_grid_list[["blq"]] <- .oblq_grid(method = method, init_comm = init_comm,
                                         criterion = criterion, criterion_type = criterion_type,
                                         abs_eigen = abs_eigen, start_method = start_method,
                                         rotation = c("promax", "oblimin", "quartimin",
                                                      "bentlerQ", "geominQ",
                                                      "bifactorQ", "simplimax"),
                                         k_promax = k_promax, normalize = normalize,
                                         P_type = P_type, precision = precision,
                                         varimax_type = varimax_type, k_simplimax = k_simplimax)

    } else {

      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" rotation = 'oblique' is used but rotation is of length > 1. Can only average EFAs with rotations of the same type ('none', 'orthogonal', or 'oblique').\n"))

    }

  } else if ("orthogonal" %in% rotation) {

    if (length(rotation) == 1) {

      t_grid_list[["rth"]] <- .orth_grid(method = method, init_comm = init_comm,
                                         criterion = criterion, criterion_type = criterion_type,
                                         abs_eigen = abs_eigen, start_method = start_method,
                                         rotation = c("varimax", "quartimax", "equamax",
                                                      "bentlerT", "geominT", "bifactorT"),
                                         normalize = normalize, precision = precision,
                                         varimax_type = varimax_type)

    } else {

      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" rotation = 'orthogonal' is used but rotation is of length > 1. Can only average EFAs with rotations of the same type ('none', 'orthogonal', or 'oblique').\n"))

    }

  } else if (all(rotation %in% c("promax", "oblimin", "quartimin", "simplimax",
                                 "bentlerQ", "geominQ", "bifactorQ"))) {

    t_grid_list[["blq2"]] <- .oblq_grid(method = method, init_comm = init_comm,
                                        criterion = criterion, criterion_type = criterion_type,
                                        abs_eigen = abs_eigen, start_method = start_method,
                                        rotation = rotation,
                                        k_promax = k_promax, normalize = normalize,
                                        P_type = P_type, precision = precision,
                                        varimax_type = varimax_type, k_simplimax = k_simplimax)

  } else if (all(rotation %in% c("varimax", "quartimax", "equamax",
                                 "bentlerT", "geominT", "bifactorT"))) {

    t_grid_list[["rth2"]] <- .orth_grid(method = method, init_comm = init_comm,
                                        criterion = criterion, criterion_type = criterion_type,
                                        abs_eigen = abs_eigen, start_method = start_method,
                                        rotation = rotation,
                                        normalize = normalize, precision = precision,
                                        varimax_type = varimax_type)

  } else if (any(rotation %in% c("promax", "oblimin", "quartimin", "simplimax",
                                 "bentlerQ", "geominQ", "bifactorQ")) &&
             any(rotation %in% c("varimax", "quartimax", "equamax",
                                 "bentlerT", "geominT", "bifactorT"))) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'rotation' contains both oblique rotations and orthogonal rotations, but can only average rotations of the same kind. Oblique rotations are 'promax', 'oblimin', 'quartimin', 'simplimax', 'bentlerQ', 'geominQ', and 'bifactorQ'. Orthogonal rotations are 'varimax', 'quartimax', 'equamax', 'bentlerT', 'geominT', and 'bifactorT'.\n"))
  }

  return(do.call(rbind, t_grid_list))
}


### Align solution according to factor congruence

.align_solution <- function(L_target, L, Phi = NULL) {

  L_target <- as.matrix(L_target)
  L <- as.matrix(L)
  if (!identical(dim(L_target), dim(L))) {
    stop("'L_target' and 'L' must have identical dimensions.", call. = FALSE)
  }

  m <- ncol(L_target)
  target_colnames <- colnames(L_target)
  target_rownames <- rownames(L_target)

  congruence <- .tucker_congruence(L_target, L)

  if (m > 1) {

    cost <- max(abs(congruence), na.rm = TRUE) - abs(congruence)

    factor_order <- as.integer(clue::solve_LSAP(cost))

    # Reorder columns of L so column j is matched to target column j.
    L <- L[, factor_order, drop = FALSE]

    if (!is.null(Phi)) {
      Phi <- Phi[factor_order, factor_order, drop = FALSE]
    }

    # Get signs after reordering.
    factor_sign <- sign(diag(crossprod(L, L_target)))

    # Avoid accidental zero columns if an exact dot product is zero.
    factor_sign[factor_sign == 0] <- 1

    factor_sign_matrix <- diag(factor_sign, nrow = m, ncol = m)

    # Apply signs.
    L <- L %*% factor_sign_matrix

    if (!is.null(Phi)) {
      Phi <- factor_sign_matrix %*% Phi %*% factor_sign_matrix
    }

  } else {

    factor_order <- 1L
    factor_sign <- as.numeric(sign(congruence[1L, 1L]))
    if (factor_sign == 0) factor_sign <- 1

    L <- L * factor_sign

    # Phi unchanged for one factor, because sign^2 = 1.
  }

  dimnames(L) <- list(target_rownames, target_colnames)
  if (!is.null(Phi)) {
    dimnames(Phi) <- list(target_colnames, target_colnames)
  }

  list(
    loadings = L,
    Phi = Phi,
    factor_order = factor_order,
    factor_sign = factor_sign
  )
}


.paste_gof_ci <- function(indices, index) {
  paste0(" [", .numformat(indices$lower[index],
                         pad = FALSE), ", ",
         .numformat(indices$upper[index],
                           pad = FALSE),
         "]")
}


#' Average a list of matrices elementwise
#'
#' @param x List of conformable matrices.
#'
#' @returns A matrix with the same dimensions as the inputs.
#'
.average_matrices <- function(x) {
  if (!is.list(x) || length(x) < 1L) {
    stop("'x' must be a non-empty list of conformable matrices.", call. = FALSE)
  }
  out <- Reduce(`+`, x) / length(x)
  dimnames(out) <- dimnames(x[[1L]])
  out
}


# Internal scalar checks for Procrustes routines
#
# These helpers intentionally use base R rather than checkmate. They are used
# by both the public wrappers and the internal consensus engine, so keeping the
# checks local gives clearer error messages and avoids partially entering the
# compiled optimizer with invalid control values.
.procrustes_is_scalar <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

.procrustes_check_numeric_scalar <- function(x, name, lower = -Inf, upper = Inf,
                                             lower_closed = TRUE,
                                             upper_closed = TRUE) {
  if (!.procrustes_is_scalar(x)) {
    stop("'", name, "' must be a finite numeric scalar.", call. = FALSE)
  }

  lower_ok <- if (lower_closed) x >= lower else x > lower
  upper_ok <- if (upper_closed) x <= upper else x < upper

  if (!lower_ok || !upper_ok) {
    left <- if (lower_closed) "[" else "("
    right <- if (upper_closed) "]" else ")"
    stop(
      "'", name, "' must be in ", left, lower, ", ", upper, right, ".",
      call. = FALSE
    )
  }

  x
}

.procrustes_check_integer_scalar <- function(x, name, lower = -Inf, upper = Inf) {
  if (!.procrustes_is_scalar(x) || abs(x - round(x)) > sqrt(.Machine$double.eps)) {
    stop("'", name, "' must be a finite integer-like scalar.", call. = FALSE)
  }
  x <- as.integer(round(x))
  if (x < lower || x > upper) {
    stop(
      "'", name, "' must be an integer in [", lower, ", ", upper, "].",
      call. = FALSE
    )
  }
  x
}

.procrustes_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.procrustes_as_matrix <- function(x, name) {
  if (inherits(x, "LOADINGS")) {
    x <- unclass(x)
  }

  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("'", name, "' must be a numeric matrix.", call. = FALSE)
  }
  storage.mode(x) <- "double"

  if (length(dim(x)) != 2L || any(dim(x) < 1L)) {
    stop("'", name, "' must have at least one row and one column.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("'", name, "' must contain only finite values.", call. = FALSE)
  }

  x
}

.procrustes_validate_matrix_pair <- function(A, B) {
  A <- .procrustes_as_matrix(A, "A")
  B <- .procrustes_as_matrix(B, "Target")

  if (!identical(dim(A), dim(B))) {
    stop("'A' and 'Target' must have identical dimensions.", call. = FALSE)
  }

  list(A = A, B = B)
}

.procrustes_validate_matrix_list <- function(x, name, min_length = 1L,
                                             expected_dim = NULL) {
  if (!is.list(x) || length(x) < min_length) {
    stop(
      "'", name, "' must be a list with at least ", min_length, " matrix/matrices.",
      call. = FALSE
    )
  }

  out <- lapply(seq_along(x), function(i) {
    .procrustes_as_matrix(x[[i]], paste0(name, "[[", i, "]]"))
  })

  ref_dim <- if (is.null(expected_dim)) dim(out[[1L]]) else expected_dim
  ok_dims <- vapply(out, function(z) identical(dim(z), ref_dim), logical(1L))
  if (!all(ok_dims)) {
    bad <- paste(which(!ok_dims), collapse = ", ")
    stop(
      "All matrices in '", name, "' must have dimensions ",
      paste(ref_dim, collapse = " x "), ". Offending element(s): ", bad, ".",
      call. = FALSE
    )
  }

  out
}

.procrustes_default_names <- function(n, prefix) {
  paste0(prefix, seq_len(n))
}

.procrustes_names <- function(A, B) {
  p <- nrow(A)
  k <- ncol(A)

  item_names <- rownames(B)
  if (is.null(item_names)) item_names <- rownames(A)
  if (is.null(item_names)) item_names <- .procrustes_default_names(p, "V")

  target_factor_names <- colnames(B)
  if (is.null(target_factor_names)) target_factor_names <- colnames(A)
  if (is.null(target_factor_names)) target_factor_names <- .procrustes_default_names(k, "F")

  source_factor_names <- colnames(A)
  if (is.null(source_factor_names)) source_factor_names <- target_factor_names

  list(
    items = item_names,
    source_factors = source_factor_names,
    target_factors = target_factor_names
  )
}

.procrustes_apply_dimnames <- function(out, A, B) {
  dn <- .procrustes_names(A, B)

  if (!is.null(out$loadings)) {
    dimnames(out$loadings) <- list(dn$items, dn$target_factors)
  }
  if (!is.null(out$T)) {
    dimnames(out$T) <- list(dn$source_factors, dn$target_factors)
  }
  if (!is.null(out$Phi)) {
    dimnames(out$Phi) <- list(dn$target_factors, dn$target_factors)
  }

  if (!is.null(out$all_values) && is.null(names(out$all_values))) {
    start_names <- if (!is.null(out$all_start_indices)) {
      paste0("start_", out$all_start_indices)
    } else {
      paste0("start_", seq_along(out$all_values))
    }
    names(out$all_values) <- start_names
  }
  if (!is.null(out$all_converged) && is.null(names(out$all_converged)) &&
      !is.null(names(out$all_values))) {
    names(out$all_converged) <- names(out$all_values)
  }
  if (!is.null(out$all_iterations) && is.null(names(out$all_iterations)) &&
      !is.null(names(out$all_values))) {
    names(out$all_iterations) <- names(out$all_values)
  }

  out
}

.procrustes_validate_crossprod <- function(S, k) {
  S <- .procrustes_as_matrix(S, "S")
  if (!identical(dim(S), c(k, k))) {
    stop("'S' must be a ", k, " x ", k, " matrix.", call. = FALSE)
  }
  if (max(abs(S - t(S))) > 1e-8 * max(1, max(abs(S)))) {
    warning("'S' is not symmetric; it should normally be crossprod(A).", call. = FALSE)
  }
  S
}

.procrustes_validate_t_init <- function(T_init, k) {
  T_init <- .procrustes_as_matrix(T_init, "T_init")
  if (!identical(dim(T_init), c(k, k))) {
    stop("'T_init' must be a ", k, " x ", k, " matrix.", call. = FALSE)
  }
  if (any(sqrt(colSums(T_init * T_init)) <= sqrt(.Machine$double.eps))) {
    stop("'T_init' must not contain zero or near-zero columns.", call. = FALSE)
  }
  if (qr(T_init)$rank < k) {
    stop("'T_init' must be nonsingular.", call. = FALSE)
  }
  T_init
}

.procrustes_validate_oblique_controls <- function(oblique_eps,
                                                  oblique_maxit,
                                                  oblique_max_line_search,
                                                  oblique_step0,
                                                  oblique_normalize,
                                                  oblique_random_starts,
                                                  oblique_screen_keep,
                                                  oblique_triage_maxit,
                                                  oblique_triage_improve_tol) {
  list(
    oblique_eps = .procrustes_check_numeric_scalar(
      oblique_eps, "oblique_eps", lower = 0, lower_closed = FALSE
    ),
    oblique_maxit = .procrustes_check_integer_scalar(
      oblique_maxit, "oblique_maxit", lower = 0L
    ),
    oblique_max_line_search = .procrustes_check_integer_scalar(
      oblique_max_line_search, "oblique_max_line_search", lower = 0L
    ),
    oblique_step0 = .procrustes_check_numeric_scalar(
      oblique_step0, "oblique_step0", lower = 0, lower_closed = FALSE
    ),
    oblique_normalize = .procrustes_check_flag(oblique_normalize, "oblique_normalize"),
    oblique_random_starts = .procrustes_check_integer_scalar(
      oblique_random_starts, "oblique_random_starts", lower = 0L
    ),
    oblique_screen_keep = .procrustes_check_integer_scalar(
      oblique_screen_keep, "oblique_screen_keep", lower = 0L
    ),
    oblique_triage_maxit = .procrustes_check_integer_scalar(
      oblique_triage_maxit, "oblique_triage_maxit", lower = 0L
    ),
    oblique_triage_improve_tol = .procrustes_check_numeric_scalar(
      oblique_triage_improve_tol, "oblique_triage_improve_tol", lower = 0
    )
  )
}

.procrustes_validate_start_index <- function(start, n_targets) {
  if (!.procrustes_is_scalar(start) || abs(start - round(start)) > sqrt(.Machine$double.eps)) {
    stop("'start' must be either an integer index or a target matrix.", call. = FALSE)
  }
  start <- as.integer(round(start))
  if (start < 1L || start > n_targets) {
    stop("'start' must be between 1 and length(init_targets).", call. = FALSE)
  }
  start
}

.procrustes_validate_starts <- function(starts, n_targets) {
  if (is.null(starts)) {
    return(seq_len(n_targets))
  }
  if (!is.numeric(starts) || length(starts) < 1L || anyNA(starts) ||
      any(!is.finite(starts)) || any(abs(starts - round(starts)) > sqrt(.Machine$double.eps))) {
    stop("'starts' must be a non-empty integer vector.", call. = FALSE)
  }
  starts <- as.integer(round(starts))
  if (any(starts < 1L | starts > n_targets)) {
    stop("All entries in 'starts' must be between 1 and length(init_targets).", call. = FALSE)
  }
  if (anyDuplicated(starts)) {
    warning("Duplicate entries in 'starts' were removed.", call. = FALSE)
    starts <- unique(starts)
  }
  starts
}

#' Tucker congruence between factors
#'
#' Compute the Tucker congruence matrix between the columns of two loading
#' matrices.
#'
#' @param L1 Numeric matrix.
#' @param L2 Numeric matrix with the same dimensions as `L1`.
#'
#' @returns A square matrix whose `(i, j)` entry is the Tucker congruence
#'   between column `i` of `L1` and column `j` of `L2`.
#'
#' @references
#' Lorenzo-Seva, U., and ten Berge, J. M. F. (2006). Tucker's congruence
#' coefficient as a meaningful index of factor similarity. *Methodology*, 2,
#' 57-64.
.tucker_congruence <- function(L1, L2) {
  L1 <- .procrustes_as_matrix(L1, "L1")
  L2 <- .procrustes_as_matrix(L2, "L2")

  if (!identical(dim(L1), dim(L2))) {
    stop("'L1' and 'L2' must have identical dimensions.", call. = FALSE)
  }

  # Column Euclidean norms.
  n1 <- sqrt(colSums(L1 * L1))
  n2 <- sqrt(colSums(L2 * L2))

  if (any(n1 <= sqrt(.Machine$double.eps))) {
    stop("L1 contains at least one zero or near-zero column.", call. = FALSE)
  }
  if (any(n2 <= sqrt(.Machine$double.eps))) {
    stop("L2 contains at least one zero or near-zero column.", call. = FALSE)
  }

  # Full pairwise congruence matrix.
  cong <- crossprod(L1, L2) / tcrossprod(n1, n2)

  rownames(cong) <- colnames(L1)
  colnames(cong) <- colnames(L2)
  cong
}

#' Count near-zero loadings
#'
#' Hyperplane count is the number of loadings with absolute value smaller than a
#' user-specified cutoff.
#'
#' @param L Numeric loading matrix.
#' @param cutoff Numeric scalar. Loadings with `abs(L) < cutoff` are counted as
#'   being in the hyperplane.
#'
#' @returns A list with the total hyperplane count and counts by factor and item.
.hyperplane_count <- function(L, cutoff = 0.15) {
  L <- .procrustes_as_matrix(L, "L")
  cutoff <- .procrustes_check_numeric_scalar(cutoff, "cutoff", lower = 0)

  in_hyperplane <- abs(L) < cutoff
  list(
    total = sum(in_hyperplane),
    by_factor = colSums(in_hyperplane),
    by_item = rowSums(in_hyperplane),
    cutoff = cutoff
  )
}

#' Mean squared discrepancy to a consensus target
#'
#' @param aligned_loadings List of aligned loading matrices.
#' @param target Consensus target matrix.
#'
#' @returns Mean sum of squared deviations from the target across matrices.
#'
.consensus_loss <- function(aligned_loadings, target) {
  # Inputs are validated by the public/internal Procrustes entry points. Keeping
  # this inner-loop helper allocation-light matters for consensus and bootstrap
  # workflows where the loss is evaluated many times.
  mean(vapply(aligned_loadings, function(L) {
    sum((L - target)^2)
  }, numeric(1L)))
}


#' Closed-form orthogonal Procrustes rotation
#'
#' Rotate `A` to the orthogonal target `B` by minimizing
#' `||A %*% T - B||_F^2` subject to `t(T) %*% T = I`.
#'
#' @param A Numeric matrix to be rotated.
#' @param B Numeric target matrix with the same dimensions as `A`.
#'
#' @returns A list with the rotated loadings, orthogonal transformation matrix,
#'   target criterion value, and basic diagnostics.
#'
#' @references
#' Schoenemann, P. H. (1966). A generalized solution of the orthogonal
#' Procrustes problem. *Psychometrika*, 31, 1-10.
.orthogonal_procrustes <- function(A, B) {
  mats <- .procrustes_validate_matrix_pair(A, B)
  A <- mats$A
  B <- mats$B

  M <- crossprod(A, B)
  s <- svd(M)
  Tmat <- s$u %*% t(s$v)

  Lrot <- A %*% Tmat
  value <- 0.5 * sum((Lrot - B)^2)

  out <- list(
    loadings = Lrot,
    T = Tmat,
    Phi = diag(ncol(A)),
    value = value,
    convergence = TRUE,
    iterations = 0L,
    kappa_T = 1,
    Table = cbind(iter = 0, f = value, log10_s = NA_real_, step = NA_real_),
    method = "orthogonal_procrustes",
    best_start_index = 1L,
    all_start_indices = 1L,
    all_values = value,
    all_converged = TRUE,
    all_iterations = 0L,
    line_search_failed = FALSE
  )

  .procrustes_apply_dimnames(out, A, B)
}


#' Internal single-start consensus engine
#'
#' This internal helper performs a single consensus-target run from one starting
#' target. Users should normally call [CONSENSUS_PROCRUSTES()] instead.
#'
#' @inheritParams CONSENSUS_PROCRUSTES
#'
.consensus_target_procrustes_single <- function(unrotated_list,
                                                init_targets = NULL,
                                                rotation = c("orthogonal", "oblique"),
                                                start = 1,
                                                tol = 1e-3,
                                                loss_tol = 1e-7,
                                                loss_patience = 5,
                                                convergence = c("either", "target", "loss", "both"),
                                                min_iter = 2,
                                                max_iter = 200,
                                                alpha = 1,
                                                match_target = TRUE,
                                                hyper_cutoff = 0.15,
                                                oblique_maxit = 300,
                                                oblique_eps = 1e-5,
                                                oblique_max_line_search = 10,
                                                oblique_step0 = 1,
                                                oblique_normalize = FALSE,
                                                oblique_random_starts = 0,
                                                oblique_random_starts_stage = c("final", "none", "outer", "both"),
                                                oblique_screen_keep = 2,
                                                oblique_triage_maxit = 25,
                                                oblique_triage_improve_tol = 0,
                                                verbose = FALSE) {
  rotation <- match.arg(rotation)
  convergence <- match.arg(convergence)
  oblique_random_starts_stage <- match.arg(oblique_random_starts_stage)

  alpha <- .procrustes_check_numeric_scalar(alpha, "alpha", lower = 0, upper = 1,
                                            lower_closed = FALSE)
  tol <- .procrustes_check_numeric_scalar(tol, "tol", lower = 0, lower_closed = FALSE)
  if (!is.null(loss_tol)) {
    loss_tol <- .procrustes_check_numeric_scalar(loss_tol, "loss_tol", lower = 0,
                                                 lower_closed = FALSE)
  }
  loss_patience <- .procrustes_check_integer_scalar(loss_patience, "loss_patience", lower = 1L)
  min_iter <- .procrustes_check_integer_scalar(min_iter, "min_iter", lower = 0L)
  max_iter <- .procrustes_check_integer_scalar(max_iter, "max_iter", lower = 1L)
  if (min_iter > max_iter) {
    stop("'min_iter' must not exceed 'max_iter'.", call. = FALSE)
  }
  if (is.null(loss_tol) && convergence %in% c("loss", "both")) {
    stop("'loss_tol' must not be NULL when convergence is 'loss' or 'both'.", call. = FALSE)
  }

  match_target <- .procrustes_check_flag(match_target, "match_target")
  verbose <- .procrustes_check_flag(verbose, "verbose")
  hyper_cutoff <- .procrustes_check_numeric_scalar(hyper_cutoff, "hyper_cutoff", lower = 0)

  oblique_controls <- .procrustes_validate_oblique_controls(
    oblique_eps = oblique_eps,
    oblique_maxit = oblique_maxit,
    oblique_max_line_search = oblique_max_line_search,
    oblique_step0 = oblique_step0,
    oblique_normalize = oblique_normalize,
    oblique_random_starts = oblique_random_starts,
    oblique_screen_keep = oblique_screen_keep,
    oblique_triage_maxit = oblique_triage_maxit,
    oblique_triage_improve_tol = oblique_triage_improve_tol
  )

  unrotated_list <- .procrustes_validate_matrix_list(
    unrotated_list, "unrotated_list", min_length = 2L
  )

  m <- length(unrotated_list)
  p <- nrow(unrotated_list[[1L]])
  k <- ncol(unrotated_list[[1L]])

  if (is.null(init_targets)) {
    init_targets <- unrotated_list
  } else {
    init_targets <- .procrustes_validate_matrix_list(
      init_targets, "init_targets", min_length = 1L, expected_dim = c(p, k)
    )
  }

  if (is.numeric(start) && length(start) == 1L) {
    start_index <- .procrustes_validate_start_index(start, length(init_targets))
    target <- init_targets[[start_index]]
    start_label <- start_index
  } else {
    target <- .procrustes_as_matrix(start, "start")
    if (!identical(dim(target), c(p, k))) {
      stop("If 'start' is a matrix, it must have the same dimensions as the loadings.", call. = FALSE)
    }
    start_index <- NA_integer_
    start_label <- "matrix"
  }

  if (!is.na(start_index)) {
    dimnames(target) <- dimnames(init_targets[[start_index]])
  }
  if (is.null(rownames(target)) || is.null(colnames(target))) {
    dn <- .procrustes_names(unrotated_list[[1L]], target)
    dimnames(target) <- list(dn$items, dn$target_factors)
  }

  S_list <- lapply(unrotated_list, crossprod)
  T_starts <- vector("list", m)
  if (rotation == "oblique" && k > 1L) {
    T_starts <- lapply(unrotated_list, function(A) .orthogonal_procrustes(A, target)$T)
  }

  outer_random_starts <- if (oblique_random_starts_stage %in% c("outer", "both")) {
    oblique_controls$oblique_random_starts
  } else {
    0L
  }

  final_random_starts <- if (oblique_random_starts_stage %in% c("final", "both")) {
    oblique_controls$oblique_random_starts
  } else {
    0L
  }

  hist_iter <- integer(max_iter)
  hist_rel_change <- numeric(max_iter)
  hist_loss <- numeric(max_iter)
  hist_rel_loss_change <- rep(NA_real_, max_iter)
  hist_inner_mean_value <- numeric(max_iter)
  hist_inner_failures <- integer(max_iter)
  hist_inner_max_kappa <- numeric(max_iter)
  hist_loss_stable_count <- integer(max_iter)
  hist_target_ok <- logical(max_iter)
  hist_loss_ok <- logical(max_iter)
  hist_stop_rule_met <- logical(max_iter)

  converged <- FALSE
  convergence_reason <- NA_character_
  prev_loss <- NA_real_
  loss_stable_count <- 0L
  n_history <- 0L

  for (iter in seq_len(max_iter)) {
    aligned <- vector("list", m)

    for (d in seq_len(m)) {
      aligned[[d]] <- PROCRUSTES(
        A = unrotated_list[[d]],
        Target = target,
        rotation = rotation,
        S = S_list[[d]],
        T_init = T_starts[[d]],
        oblique_eps = oblique_controls$oblique_eps,
        oblique_maxit = oblique_controls$oblique_maxit,
        oblique_max_line_search = oblique_controls$oblique_max_line_search,
        oblique_step0 = oblique_controls$oblique_step0,
        oblique_normalize = oblique_controls$oblique_normalize,
        oblique_random_starts = outer_random_starts,
        oblique_screen_keep = oblique_controls$oblique_screen_keep,
        oblique_triage_maxit = oblique_controls$oblique_triage_maxit,
        oblique_triage_improve_tol = oblique_controls$oblique_triage_improve_tol
      )
      if (rotation == "oblique" && k > 1L) {
        T_starts[[d]] <- aligned[[d]]$T
      }
    }

    aligned_loadings <- lapply(aligned, `[[`, "loadings")
    inner_values <- vapply(aligned, `[[`, numeric(1L), "value")
    inner_converged <- vapply(aligned, function(x) isTRUE(x$convergence), logical(1L))
    inner_kappa <- vapply(aligned, function(x) {
      if (is.null(x$kappa_T)) NA_real_ else x$kappa_T
    }, numeric(1L))

    centroid <- .average_matrices(aligned_loadings)
    if (isTRUE(match_target) && k > 1L) {
      centroid <- .align_solution(L = centroid, L_target = target)$loadings
      dimnames(centroid) <- dimnames(target)
    }

    new_target <- (1 - alpha) * target + alpha * centroid
    dimnames(new_target) <- dimnames(target)

    rel_change <- norm(new_target - target, type = "F") /
      (norm(target, type = "F") + 1e-12)

    current_loss <- .consensus_loss(aligned_loadings, new_target)

    rel_loss_change <- if (is.na(prev_loss)) {
      NA_real_
    } else {
      abs(prev_loss - current_loss) / (abs(prev_loss) + 1e-12)
    }

    if (!is.null(loss_tol) && !is.na(rel_loss_change) && rel_loss_change < loss_tol) {
      loss_stable_count <- loss_stable_count + 1L
    } else {
      loss_stable_count <- 0L
    }

    target_ok <- rel_change < tol
    loss_ok <- !is.null(loss_tol) && loss_stable_count >= loss_patience

    stop_rule_met <- switch(
      convergence,
      target = target_ok,
      loss = loss_ok,
      either = target_ok || loss_ok,
      both = target_ok && loss_ok
    )

    n_history <- iter
    hist_iter[iter] <- iter
    hist_rel_change[iter] <- rel_change
    hist_loss[iter] <- current_loss
    hist_rel_loss_change[iter] <- rel_loss_change
    hist_inner_mean_value[iter] <- mean(inner_values)
    hist_inner_failures[iter] <- sum(!inner_converged)
    hist_inner_max_kappa[iter] <- if (all(is.na(inner_kappa))) NA_real_ else max(inner_kappa, na.rm = TRUE)
    hist_loss_stable_count[iter] <- loss_stable_count
    hist_target_ok[iter] <- target_ok
    hist_loss_ok[iter] <- loss_ok
    hist_stop_rule_met[iter] <- stop_rule_met

    if (verbose) {
      message(sprintf(
        paste0(
          "iter %d | rel_change = %.3e | loss = %.8f | ",
          "rel_loss_change = %s | stable = %d | inner_failures = %d"
        ),
        iter,
        rel_change,
        current_loss,
        ifelse(is.na(rel_loss_change), "NA", sprintf("%.3e", rel_loss_change)),
        loss_stable_count,
        sum(!inner_converged)
      ))
    }

    target <- new_target
    prev_loss <- current_loss

    if (iter >= min_iter && isTRUE(stop_rule_met)) {
      converged <- TRUE
      convergence_reason <- if (target_ok && loss_ok) {
        "target_and_loss"
      } else if (target_ok) {
        "target"
      } else if (loss_ok) {
        "loss"
      } else {
        convergence
      }
      break
    }
  }

  keep <- seq_len(n_history)
  history <- data.frame(
    iter = hist_iter[keep],
    rel_change = hist_rel_change[keep],
    loss = hist_loss[keep],
    rel_loss_change = hist_rel_loss_change[keep],
    inner_mean_value = hist_inner_mean_value[keep],
    inner_failures = hist_inner_failures[keep],
    inner_max_kappa = hist_inner_max_kappa[keep],
    loss_stable_count = hist_loss_stable_count[keep],
    target_ok = hist_target_ok[keep],
    loss_ok = hist_loss_ok[keep],
    stop_rule_met = hist_stop_rule_met[keep]
  )

  final_aligned <- vector("list", m)
  for (d in seq_len(m)) {
    final_aligned[[d]] <- PROCRUSTES(
      A = unrotated_list[[d]],
      Target = target,
      rotation = rotation,
      S = S_list[[d]],
      T_init = T_starts[[d]],
      oblique_eps = oblique_controls$oblique_eps,
      oblique_maxit = oblique_controls$oblique_maxit,
      oblique_max_line_search = oblique_controls$oblique_max_line_search,
      oblique_step0 = oblique_controls$oblique_step0,
      oblique_normalize = oblique_controls$oblique_normalize,
      oblique_random_starts = final_random_starts,
      oblique_screen_keep = oblique_controls$oblique_screen_keep,
      oblique_triage_maxit = oblique_controls$oblique_triage_maxit,
      oblique_triage_improve_tol = oblique_controls$oblique_triage_improve_tol
    )
  }

  final_aligned_loadings <- lapply(final_aligned, `[[`, "loadings")
  final_T <- lapply(final_aligned, `[[`, "T")
  final_Phi <- lapply(final_aligned, `[[`, "Phi")
  final_values <- vapply(final_aligned, `[[`, numeric(1L), "value")
  final_converged <- vapply(final_aligned, function(x) isTRUE(x$convergence), logical(1L))
  final_kappa <- vapply(final_aligned, function(x) {
    if (is.null(x$kappa_T)) NA_real_ else x$kappa_T
  }, numeric(1L))

  pooled_loadings <- .average_matrices(final_aligned_loadings)
  dimnames(pooled_loadings) <- dimnames(target)

  pooled_phi <- if (rotation == "orthogonal" || k == 1L) {
    diag(k)
  } else {
    .average_matrices(final_Phi)
  }
  pooled_phi <- (pooled_phi + t(pooled_phi)) / 2
  diag(pooled_phi) <- 1
  dimnames(pooled_phi) <- list(colnames(target), colnames(target))

  final_loss <- .consensus_loss(final_aligned_loadings, target)

  outer_converged <- converged
  final_inner_converged <- all(final_converged)
  overall_converged <- outer_converged && final_inner_converged
  reported_convergence_reason <- if (overall_converged) {
    convergence_reason
  } else if (!outer_converged) {
    convergence_reason
  } else {
    "final_inner_alignment"
  }

  list(
    converged = overall_converged,
    outer_converged = outer_converged,
    final_inner_converged = final_inner_converged,
    convergence_reason = reported_convergence_reason,
    outer_convergence_reason = convergence_reason,
    iterations = nrow(history),
    tol = tol,
    loss_tol = loss_tol,
    loss_patience = loss_patience,
    convergence = convergence,
    min_iter = min_iter,
    alpha = alpha,
    match_target = match_target,
    rotation = rotation,
    start = start_label,
    start_index = start_index,
    target = target,
    history = history,
    aligned_loadings = final_aligned_loadings,
    aligned_phi = final_Phi,
    transformations = final_T,
    pooled_loadings = pooled_loadings,
    pooled_phi = pooled_phi,
    mean_loss = final_loss,
    final_values = final_values,
    final_converged = final_converged,
    final_failures = which(!final_converged),
    final_kappa_T = final_kappa,
    hyperplane_target = .hyperplane_count(target, cutoff = hyper_cutoff),
    hyperplane_pooled = .hyperplane_count(pooled_loadings, cutoff = hyper_cutoff),
    oblique_random_starts_stage = oblique_random_starts_stage,
    oblique_random_starts_outer = outer_random_starts,
    oblique_random_starts_final = final_random_starts
  )
}
