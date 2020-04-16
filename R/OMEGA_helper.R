# Flexible omega function (e.g. to use with loadings obtained by MacOrtho)------
.OMEGA_FLEX <- function(model = NULL, type = c("EFAtools", "psych"), factor_corres = NULL,
                        var_names = NULL, fac_names = NULL, g_load = NULL,
                        s_load = NULL, u2 = NULL, cormat = NULL, pattern = NULL,
                        Phi = NULL, variance = c("correlaton", "sums_load")){

  if(all(class(model) == c("psych", "schmid"))){

    pattern <- model$oblique
    Phi <- model$phi

    model <-  model$sl

    var_names <- rownames(model)
    g_load <- model[, 1]
    s_load <- model[, 2:(ncol(model) - 3)]
    factor_names <- c("g", 1:ncol(s_load))
    u2 <- model[, "u2"]

    if(variance == "correlation"){

      if(is.null(cormat)){

        # Create the correlation matrix from the pattern coefficients and factor
        # intercorrelations
        cormat <- psych::factor.model(f = pattern, Phi = Phi, U2 = FALSE)

      }

      } else {

        # Check if it is a correlation matrix
        if(.is_cormat(cormat)){

          if(any(is.na(cormat))){

            stop("The correlation matrix you entered contains missing values.
                 Check the cormat input, specify the Phi and pattern arguments
                 instead, or set variance to 'sums_load'")

          }

        } else {

          stop("x was not a correlation matrix. Check the cormat input, specify
                the Phi and pattern arguments instead, or set variance to
                'sums_load'")

        }

    }

 } else if(all(class(model) == c("SL"))){

    cormat <- model$orig_R

    model <-  model$sl
    var_names <- rownames(model)
    g_load <- model[, 1]
    s_load <- model[, 2:(ncol(model) - 2)]
    factor_names <- c("g", 1:ncol(s_load))
    u2 <- model[, "u2"]

  } else {

    factor_names <- c("g", 1:ncol(s_load))

    if(variance == "correlation"){

      if(is.null(cormat)){

        if(is.null(Phi) | is.null(pattern)) {
          stop("If you leave model NULL, either specify the cormat
             argument or the Phi and pattern arguments, or set variance to
             'sums_load'")

        } else {

          # Create the correlation matrix from the pattern coefficients and factor
          # intercorrelations
          cormat <- psych::factor.model(f = pattern, Phi = Phi, U2 = FALSE)

        }

      } else {

        # Check if it is a correlation matrix
        if(.is_cormat(cormat)){

          if(any(is.na(cormat))){

            stop("The correlation matrix you entered contains missing values.
                 Check the cormat input, specify the Phi and pattern arguments
                 instead, or set variance to 'sums_load'")

          }

        } else {

          stop("x was not a correlation matrix. Check the cormat input, specify
                the Phi and pattern arguments instead, or set variance to
                'sums_load'")

        }

      }

    }

  }

  # Create an input dataframe
  input <- data.frame(g_load, s_load)
  colnames(input) <- factor_names
  rownames(input) <- var_names

  if(type == "EFAtools" & is.null(factor_corres)){

    stop("Either specify the factor_corres argument or set type = 'psych'
         to find variable-to-factor correspondences using the highest
         group factor loading per variable.")

  }

  if(type == "psych"){

    if(variance != "correlation"){

      warning("Variance is specified. Variance is used with value '", variance,
              "'. Results may differ from the specified type")
      }

    if(is.null(factor_corres)){
      factor_corres <- apply(input, 1,
                             function(x) which.max(abs(x[2:(ncol(s_load) + 1)])))

    } else {

      warning("Argument factor_corres is specified. Specified variable-to-factor
              correspondences are taken. To compute factor correspondences as done
              in psych, leave factor_corres = NULL.")
    }
  }

  input$u2 <- u2
  input <- cbind(factor_corres, input)
  names(input)[1] <- "factor"

  # Create all sums of factor loadings for each factor

  # Sums of all group factor loadings for each group factor
  sums_s <- colSums(input[, 3:(ncol(s_load) + 2)])

  # Sum of all g loadings
  sum_g <- sum(input$g)

  # Sum of all error variances
  sum_e <- sum(input$u2)

  # Compute sums of error variances and g-loadings for group factors
  sums_e_s <- NULL
  sums_g_s <- NULL
  sums_s_s <- NULL

  for (i in 1:ncol(s_load)){
    sums_e_s[i] <- sum(input[input$factor == i, "u2"])
    sums_g_s[i] <- sum(input[input$factor == i, "g"])
    sums_s_s[i] <- sum(input[input$factor == i, i + 2])
  }

  if(variance == "correlation"){

    # Compute omega total, hierarchical, and subscale for g-factor
    omega_tot_g <- (sum(cormat) - sum_e) / sum(cormat)
    omega_h_g <- sum_g^2 / sum(cormat)
    omega_sub_g <- sum(sums_s_s^2) / sum(cormat)

    # Compute omega total, hierarchical, and subscale for group factors
    omega_tot_sub <- NULL
    omega_h_sub <- NULL
    omega_sub_sub <- NULL

    for (i in 1:ncol(s_load)) {
      subf <- which(factor_corres == i)
      Vgr <- sum(cormat[subf, subf])
      omega_sub_sub[i] <- sums_s_s[i]^2 / Vgr
      omega_h_sub[i] <- sums_g_s[i]^2 / Vgr
      omega_tot_sub[i] <- (sums_s_s[i]^2 + sums_g_s[i]^2) / Vgr
    }

  } else if(variance == "sums_load") {

    # Compute omega total, hierarchical, and subscale for g-factor
    omega_tot_g <- (sum_g^2 + sum(sums_s_s^2)) / (sum_g^2 + sum(sums_s_s^2) +
                                                    sum_e)
    omega_h_g <- sum_g^2 / (sum_g^2 + sum(sums_s^2) + sum_e)
    omega_sub_g <- sum(sums_s_s^2) / (sum_g^2 + sum(sums_s^2) + sum_e)

    # Compute omega total, hierarchical, and subscale for group factors
    omega_tot_sub <- (sums_g_s^2 + sums_s_s^2) / (sums_g_s^2 + sums_s_s^2 +
                                                         sums_e_s)
    omega_h_sub <- sums_g_s^2 / (sums_g_s^2 + sums_s_s^2 + sums_e_s)
    omega_sub_sub <- sums_s_s^2 / (sums_g_s^2 + sums_s_s^2 + sums_e_s)
  }

  # Combine and display results in a table
  omega_tot <- c(omega_tot_g, omega_tot_sub)
  omega_h <- c(omega_h_g, omega_h_sub)
  omega_sub <- c(omega_sub_g, omega_sub_sub)

  omegas <- cbind(omega_tot, omega_h, omega_sub)
  colnames(omegas) <- c("tot", "hier", "sub")

  if(!is.null(fac_names)){

    rownames(omegas) <- c("g", fac_names)

  } else {

    if(is.null(model)){

      rownames(omegas) <- c("g", 1:ncol(s_load))

    } else {

      rownames(omegas) <- c("g", colnames(model)[2:(ncol(s_load) + 1)])

    }
  }

  omegas

}


# Omega function to use with lavaan bifactor output as input-------

.OMEGA_LAVAAN <- function(model = NULL, g_name = NULL){

  std_sol <- lavaan::inspect(model, what = "std")

  if(lavaan::inspect(model, what = "converged") == FALSE){
    stop("Model did not converge. No omegas are computed.")
  }

  if(any(is.na(std_sol$lambda))){
    stop("Some loadings are NA or NaN. No omegas are computed.")
  }

  if(any(std_sol$lambda > 1)){
    stop("A Heywood case was detected. No omegas are computed.")
  }

  if(any(std_sol$lambda == 1)){
    warning("Perfect relationship (loading equal to 1) was detected.
            At least one variable is redundant.")
  }

  if(ncol(std_sol$lambda) == 1){

    message("The model contains a single factor. Only omega total is computed")

    fac_names <- colnames(std_sol$lambda)
    var_names <- rownames(std_sol$lambda)

    # Extract sum of factor loadings and error variances
    sum_g <- sum(std_sol$lambda[, fac_names])
    e_load <- diag(std_sol$theta)
    sum_e <- sum(e_load)

    # Compute omega
    omegas <- (sum_g^2) / (sum_g^2 + sum_e)
    names(omegas) <- "omega"

  } else {

    bi_check <- std_sol$lambda
    bi_check[bi_check != 0] <- 1

    if(all(rowSums(bi_check) < 2)){

      warning("You did not fit a bifactor model. Omegas cannot be computed.
              Either provide a bifactor model or a model with a single factor.")

    }

    if(!all(rowSums(bi_check) > 1)){

      message("Some variables have less than two loadings. Did you really enter
              a bifactor model? Either provide a bifactor model or a model with a
              single factor.")

    }

    # Create list with factor and corresponding subtest names
    col_names <- colnames(std_sol$lambda)
    col_names <- col_names[!col_names %in% g_name]
    fac_names <- c(g_name, col_names)

    var_names <- list()
    for(i in 1:(length(fac_names)-1)){
      temp0 <- std_sol$lambda[, i] != 0
      var_names[[i]] <- names(temp0)[temp0]
    }

    # Create all sums of factor loadings for each factor
    sums_all <- NULL
    for (i in 1:length(fac_names)){
      temp <- sum(std_sol$lambda[, fac_names[i]])
      sums_all[i] <- temp
    }

    # Extract g-loadings and error variances, sum of g-loadings for g and sum of
    # respective loadings for every factor
    g_load <- std_sol$lambda[, g_name]
    e_load <- diag(std_sol$theta)
    sum_e <- sum(e_load)
    sum_g <- sums_all[1]
    sums_s <- sums_all[2:length(fac_names)]

    # Compute sums of error variances and g-loadings for group factors
    sums_e_s <- NULL
    sums_g_s <- NULL
    for (i in 1:length(var_names)){
      sums_e_s[i] <- sum(e_load[var_names[[i]]])
      sums_g_s[i] <- sum(g_load[var_names[[i]]])
    }

    # Compute omega total, hierarchical, and subscale for g-factor
    omega_tot_g <- (sum_g^2 + sum(sums_s^2)) / (sum_g^2 + sum(sums_s^2) + sum_e)
    omega_h_g <- sum_g^2 / (sum_g^2 + sum(sums_s^2) + sum_e)
    omega_sub_g <- sum(sums_s^2) / (sum_g^2 + sum(sums_s^2) + sum_e)

    # Compute omega total, hierarchical, and subscale for group factors
    omega_tot_sub <- (sums_s^2 + sums_g_s^2) / (sums_g_s^2 + sums_s^2 + sums_e_s)
    omega_h_sub <- sums_g_s^2 / (sums_g_s^2 + sums_s^2 + sums_e_s)
    omega_sub_sub <- sums_s^2 / (sums_g_s^2 + sums_s^2 + sums_e_s)

    # Combine and display results in a table
    omega_tot <- c(omega_tot_g, omega_tot_sub)
    omega_h <- c(omega_h_g, omega_h_sub)
    omega_sub <- c(omega_sub_g, omega_sub_sub)

    omegas <- cbind(omega_tot, omega_h, omega_sub)
    colnames(omegas) <- c("tot", "hier", "sub")

    rownames(omegas) <- fac_names

    }

  omegas

}

