# Flexible omega function (e.g. to use with loadings obtained by MacOrtho)-------
OMEGA_FLEX <- function(model = NULL, var_names, fac_names = NULL, factor_corres = NULL,
                       g_load, s_load, u2 = NULL, Phi = NULL, pattern = NULL,
                       type = "EFAdiff"){

  if(all(class(model) == c("psych", "schmid"))){

    model <-  model$sl
    var_names <- rownames(model)
    g_load <- model$sl[, 1]
    s_load <- model$sl[, 2:(ncol(model) - 3)]
    factor_names <- c("g", 1:ncol(s_load))

    if(type != "Watkins"){
      u2 <- model[, ncol(model) - 1]
    }

  } else if(all(class(model) == c("SL"))){

    model <-  model$sl
    var_names <- rownames(model)
    g_load <- model$sl[, 1]
    s_load <- model$sl[, 2:(ncol(model) - 2)]
    factor_names <- c("g", 1:ncol(s_load))

    if(type != "Watkins"){
      u2 <- model[, ncol(model)]
    }

  } else {

    factor_names <- c("g", 1:ncol(s_load))

    }

  # Create an input dataframe
  input <- data.frame(g_load, s_load)
  colnames(input) <- factor_names
  rownames(input) <- var_names

  if(type != "psych" & is.null(factor_corres)){
    stop("Either specify the factor_corres argument or set type = 'psych'")
  }

  if(type != "Watkins" & is.null(u2)){
    stop("Either specify the u2 argument or set type = 'Watkins'")
  }

  if(type == "psych"){
    if(is.null(Phi) | is.null(pattern)){
      stop("Please specify the Phi and pattern arguments")
    } else {

      # Create the correlation matrix from the pattern coefficients and factor
      # intercorrelations
      cormat <- psych::factor.model(f = pattern, Phi = Phi, U2 = FALSE)
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

  input <- cbind(factor_corres, input)
  names(input)[1] <- "factor"

  if(type == "Watkins"){

    if(is.null(u2)) {

      # Set all "non-relevant" group factor loadings to 0
      Watkins_data <- matrix(ncol = ncol(s_load), nrow = length(var_names))

      for(i in 1:ncol(s_load)){
        temp <- ifelse(input$factor == i,
                       input[input$factor == i, colnames(input)[i + 2]], 0)
        Watkins_data[, i] <- temp
      }

      input[, 1:ncol(s_load) + 2] <- Watkins_data

      # Calculate uniquenesses based only on the relevant group factor loadings
      u2_Wat <- 1 - rowSums(input[, c(2, 1:ncol(s_load) + 2)]^2)
      u2 <- u2_Wat

    } else {

      warning("Argument u2 is specified. Specified uniquenesses are taken.
              To compute uniquenesses as done in Watkins' program, leave u2
              = NULL")
    }
    }

  input$u2 <- u2

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

  if(type == "psych"){

    # Compute omega total, hierarchical, and subscale for g-factor
    omega_tot <- (sum(cormat) - sum(input$u2)) / sum(cormat)
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

  } else {

    # Compute omega total, hierarchical, and subscale for g-factor
    omega_tot_g <- (sum_g^2 + sum(sums_s^2)) / (sum_g^2 + sum(sums_s^2) + sum_e)
    omega_h_g <- sum_g^2 / (sum_g^2 + sum(sums_s^2) + sum_e)
    omega_sub_g <- sum(sums_s_s^2) / (sum_g^2 + sum(sums_s^2) + sum_e)

    # Compute omega total, hierarchical, and subscale for group factors
    omega_tot_sub <- NULL
    omega_h_sub <- NULL
    omega_sub_sub <- NULL

    omega_tot_sub <- (sums_s_s^2 + sums_g_s^2) / (sums_g_s^2 + sums_s_s^2 + sums_e_s)
    omega_h_sub <- sums_g_s^2 / (sums_g_s^2 + sums_s_s^2 + sums_e_s)
    omega_sub_sub <- sums_s_s^2 / (sums_g_s^2 + sums_s_s^2 + sums_e_s)
  }

  # Combine and display results in a table
  omega_tot <- c(omega_tot_g, omega_tot_sub)
  omega_h <- c(omega_h_g, omega_h_sub)
  omega_sub <- c(omega_sub_g, omega_sub_sub)

  omegas <- cbind(omega_tot, omega_h, omega_sub)

  if(!is.null(fac_names)){
    rownames(omegas) <- c("g", fac_names)
  } else {
    rownames(omegas) <- c("g", 1:ncol(s_load))
  }

  omegas

    }

# Omega function to use with lavaan bifactor output as input-------
OMEGA_LAVAAN <- function(model){

  if(lavaan::inspect(model, what = "converged") == FALSE){
    stop("Model did not converge. No omegas are computed.")
  }

  if(any(is.na(lavaan::inspect(model, what = "std")$lambda))){
    stop("Some loadings are NA or NaN. No omegas are computed.")
  }

  if(any(lavaan::inspect(model, what = "std")$lambda > 1)){
    stop("A Heywood case was detected. No omegas are computed.")
  }

  if(any(lavaan::inspect(model, what = "std")$lambda == 1)){
    warning("Perfect relationship (loading equal to 1) was detected.
            At least one variable is redundant.")
  }

  # Create list with factor and subtest names
  col_names <- colnames(lavaan::inspect(model, what = "std")$lambda)
  col_names <- utils::head(col_names, -1)
  fac_names <- c("g", col_names)

  var_names <- list()
  for(i in 1:(length(fac_names)-1)){
    temp0 <- lavaan::inspect(model, what = "std")$lambda[, i] != 0
    var_names[[i]] <- names(temp0)[temp0]
  }

  # Create all sums of factor loadings for each factor
  sums_all <- NULL
  for (i in 1:length(fac_names)){
    temp <- sum(lavaan::inspect(model, what = "std")$lambda[, fac_names[i]])
    sums_all[i] <- temp
  }

  # Extract g-loadings and error variances, sum of g-loadings for g and sum of
  # respective loadings for every factor
  g_load <- lavaan::inspect(model, what = "std")$lambda[,"g"]
  e_load <- diag(lavaan::inspect(model, what = "std")$theta)
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

  omegas <- matrix(omega_tot, omega_h, omega_sub)
  rownames(omegas) <- fac_names

  omegas

}

