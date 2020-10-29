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
#'
#' @return A formated number
.numformat <- function(x, digits = 2, print_zero = FALSE) {

  if (isFALSE(print_zero)) {

    ncode <- paste0("%.", digits, "f")
    x <- sub("^(-?)0.", "\\1.", sprintf(ncode, x))
    x <- stringr::str_pad(x, digits + 2, "left")

  } else {

    ncode <- paste0("%.", digits, "f")
    x <- sprintf(ncode, x)
    x <- stringr::str_pad(x, digits + 3, "left")

  }

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
  CAF <- 1 - KMO(delta_hat)$KMO


  if (method != "PAF" && !is.na(N)) {

    ### compute CFI

    # null model
    chi_null <- sum(R[upper.tri(R)] ^ 2) * (N - 1)
    df_null <- (m**2 - m) / 2
    delta_hat_null <- chi_null - df_null
    p_null <- stats::pchisq(chi_null, df_null, lower.tail = F)

    # current model
    chi <- Fm * (N - 1)
    p_chi <- stats::pchisq(chi, df, lower.tail = F)
    delta_hat_m <- chi - df
    CFI <- (delta_hat_null - delta_hat_m) / delta_hat_null
    if (CFI > 1 || df == 0) CFI <- 1
    if (CFI < 0) CFI <- 0


    ### compute RMSEA, incl. 90% confidence intervals if df are not 0
    if(df != 0){

      # formula 12.6 from Kline 2015; Principles and practices of...
      RMSEA <- sqrt(max(0, chi - df) / (df * N - 1))

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

    RMSEA_LB <- sqrt(lambda_l / (df * N))
    RMSEA_UB <- sqrt(lambda_u / (df * N))

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

    if (round(sum(diag(x), na.rm = TRUE)) == nrow(x)) {

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
          has_heywood <- any(efa_temp$h2 >= .998) || any(efa_temp[[load_ind]] >= .998)
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

      # obtain signs to reflect signs of Ln if necessary
      factor_sign <- sapply(seq_len(n_factors),
                            function(ll, congruence, factor_order){
                              sign(congruence[ll, factor_order[ll]])
                            }, congruence = congruence,
                            factor_order = factor_order)

      factor_sign <- rep(factor_sign, each = nrow(L1))

      # reorder
      L[,, efa_i] <- Ln[, factor_order] * factor_sign
      L_corres[,, efa_i] <- L_corres[,, efa_i][, factor_order]
      vars_accounted[,, efa_i] <- vars_accounted[,, efa_i][, factor_order]
      if (isTRUE(extract_phi)) {
        phi[,, efa_i] <- phi[,, efa_i][factor_order, factor_order]
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


