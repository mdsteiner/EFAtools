# Flexible omega function (e.g. to use with loadings obtained by MacOrtho)------
.OMEGA_FLEX <- function(model = NULL, type = c("EFAtools", "psych"),
                        factor_corres = NULL,
                        var_names = NULL, fac_names = NULL, g_load = NULL,
                        s_load = NULL, u2 = NULL, cormat = NULL, pattern = NULL,
                        Phi = NULL, variance = c("correlation", "sums_load"),
                        add_ind = TRUE){

  if(inherits(model, "schmid")){

    pattern <- model$oblique
    Phi <- model$phi

    model <-  model$sl
    s_load_names <- setdiff(colnames(model[, -1]),
                            c("h2", "u2", "p2", "com"))

    var_names <- rownames(model)
    g_load <- model[, 1]
    s_load <- model[, s_load_names]
    factor_names <- c("g", seq_len(ncol(s_load)))
    u2 <- model[, "u2"]

 } else if(inherits(model, "SL")){

    # Honour a user-supplied cormat; only fall back to the SL object's stored
    # correlation matrix when none was given. Flexible-input SL objects store
    # orig_R = NA, in which case cormat stays NULL and the checks below apply.
    if(is.null(cormat) && is.matrix(model$orig_R)){
      cormat <- model$orig_R
    }

    model <-  model$sl
    var_names <- rownames(model)
    g_load <- model[, 1]
    s_load <- model[, 2:(ncol(model) - 2)]
    factor_names <- c("g", seq_len(ncol(s_load)))
    u2 <- model[, "u2"]

  } else {

    factor_names <- c("g", seq_len(ncol(s_load)))

  }

    if(variance == "correlation"){

      if(is.null(cormat)){

        if(is.null(Phi) | is.null(pattern)) {
          cli::cli_abort(
            c("Specify either {.arg cormat}, or {.arg Phi} and {.arg pattern}.",
              "i" = "Alternatively, set {.code variance = \"sums_load\"}."),
            class = "efa_omega_need_cormat"
          )

        } else {

          # Create the correlation matrix from the pattern coefficients and factor
          # intercorrelations
          cormat <- psych::factor.model(f = pattern, Phi = Phi, U2 = FALSE)

        }

      } else {

        # Check if it is a correlation matrix
        if(!.is_cormat(cormat)) {

          cli::cli_abort(
            c("{.arg cormat} is not a correlation matrix.",
              "i" = "Check the {.arg cormat} input, supply {.arg Phi} and {.arg pattern} instead, or set {.code variance = \"sums_load\"}."),
            class = "efa_omega_not_cormat"
          )

        }

      }

    }

  # Check if input to factor_corres is correct (g_load is a vector here, so the
  # item count is its length, equivalently nrow(s_load)).
  checkmate::assert_matrix(factor_corres, null.ok = TRUE, nrows = length(g_load),
                           ncols = ncol(s_load))

  # Create an input dataframe
  input <- data.frame(g_load, s_load)
  colnames(input) <- factor_names
  rownames(input) <- var_names

  if(type == "EFAtools" & is.null(factor_corres)){

    cli::cli_abort("Specify {.arg factor_corres}, or set {.code type = \"psych\"} to derive variable-to-factor correspondences from the highest group-factor loading per variable.",
                   class = "efa_omega_need_corres")

  }

  if(type == "psych"){

    if(variance != "correlation"){

      cli::cli_warn(
        c("{.arg variance} is specified; the value {.val {variance}} is used.",
          "i" = "Results may differ from the specified {.arg type}."),
        class = "efa_omega_variance_override"
      )
      }

    if(is.null(factor_corres)){

      factor_corres <- matrix(0, nrow = nrow(s_load), ncol= ncol(s_load))

      for(i in seq_len(nrow(s_load))){

        factor_corres[i, which.max(abs(input[i, 2:(ncol(s_load) + 1)]))] <- 1

      }



    } else {

      cli::cli_warn(
        c("{.arg factor_corres} is specified; the supplied variable-to-factor correspondences are used.",
          "i" = "To compute correspondences as in psych, leave {.code factor_corres = NULL}."),
        class = "efa_omega_corres_override"
      )
    }
  }

  input$u2 <- u2

  omega_mat <- input[, 2:(ncol(s_load) + 1)] * factor_corres

  # Sum of all g loadings
  sum_g <- sum(input$g)

  # Sum of all error variances
  sum_e <- sum(input$u2)

  # Compute sums of error variances and g-loadings for group factors
  sums_e_s <- NULL
  sums_g_s <- NULL

  for (i in seq_len(ncol(s_load))){
    sums_e_s[i] <- sum(input$u2 * factor_corres[, i])
    sums_g_s[i] <- sum(input$g * factor_corres[, i])
  }

  sums_s_s <- colSums(omega_mat)

  if(variance == "correlation"){

    # Compute omega total, hierarchical, and subscale for g-factor
    omega_tot_g <- (sum(cormat) - sum_e) / sum(cormat)
    omega_h_g <- sum_g^2 / sum(cormat)
    omega_sub_g <- sum(sums_s_s^2) / sum(cormat)

    # Compute omega total, hierarchical, and subscale for group factors
    omega_tot_sub <- NULL
    omega_h_sub <- NULL
    omega_sub_sub <- NULL

    for (i in seq_len(ncol(s_load))) {
      subf <- which(factor_corres[, i] == 1)
      Vgr <- sum(cormat[subf, subf])
      omega_sub_sub[i] <- sums_s_s[i]^2 / Vgr
      omega_h_sub[i] <- sums_g_s[i]^2 / Vgr
      omega_tot_sub[i] <- (sums_s_s[i]^2 + sums_g_s[i]^2) / Vgr
    }

  } else if(variance == "sums_load") {

    # Sums of all group factor loadings for all group factors
    sums_s <- colSums(input[, 2:(ncol(s_load) + 1)])

    # Compute omega total, hierarchical, and subscale for g-factor. All three
    # share the same (unzeroed) total-variance denominator, so that
    # tot = hier + sub for the g row.
    omega_tot_g <- (sum_g^2 + sum(sums_s_s^2)) / (sum_g^2 + sum(sums_s^2) +
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

  if(isTRUE(add_ind)){

    # Compute H index, ECV, and PUC
    sum_s_load_sq <- vector("double", ncol(omega_mat))
    h_s_load <- vector("double", ncol(omega_mat))

    # Calculate denominators for H for group factors and for ECV. The H index
    # and ECV are based on the (unzeroed) group-factor loadings, as in
    # psych::omega and Rodriguez et al. (2016); the correspondence map
    # (omega_mat) only defines item-to-factor membership, which enters the PUC.
    for (j in seq_along(omega_mat)){
      s_j <- input[[j + 1]]
      sum_s_load_sq[j] <- sum(s_j^2 / (1 - s_j^2))
      h_s_load[j] <- sum(s_j^2)
    }

    h_g <- 1 / (1 + (1 / sum(input$g^2 / (1 - input$g^2))))
    h_s <- 1 / (1 + (1 / sum_s_load_sq))

    ECV <- sum(input$g^2) / sum(sum(input$g^2), sum(h_s_load))

    # Calculate number of items per factor
    item_on_fac <- colSums(omega_mat != 0)

    # Calculate number of contaminated correlations
    cont_corrs <- sum(sapply(item_on_fac, function (x) {x*(x - 1) / 2}))

    # Calculate proportion of uncontaminated correlations (PUC)
    PUC <- 1 - cont_corrs/(nrow(omega_mat)*(nrow(omega_mat) - 1) / 2)

    # Create output
    h <- c(h_g, h_s)

    omegas <- cbind(omega_tot, omega_h, omega_sub, h, NA, NA)
    omegas[1, 5] <- ECV
    omegas[1, 6] <- PUC
    colnames(omegas) <- c("tot", "hier", "sub", "H", "ECV", "PUC")

  } else {

    omegas <- cbind(omega_tot, omega_h, omega_sub)
    colnames(omegas) <- c("tot", "hier", "sub")

  }


  if(!is.null(fac_names)){

    rownames(omegas) <- c("g", fac_names)

  } else {

    if(is.null(model)){

      rownames(omegas) <- c("g", seq_len(ncol(s_load)))

    } else {

      rownames(omegas) <- c("g", colnames(model)[2:(ncol(s_load) + 1)])

    }
  }

  class(omegas) <- "OMEGA"

  return(omegas)

}


# Omega function to use with lavaan bifactor output as input-------

.OMEGA_LAVAAN <- function(model = NULL, g_name = "g", group_names = NULL,
                          add_ind = TRUE){

  higherorder <- FALSE

  if(lavaan::lavInspect(model, what = "converged") == FALSE){
    cli::cli_abort("The model did not converge; no omegas are computed.",
                   class = "efa_omega_no_converge")
  }

  std_sol <- suppressWarnings(lavaan::lavInspect(model, what = "std",
                                drop.list.single.group = FALSE))

  ## Create empty list for output
  omegas <- list()

  # Get group names from lavaan output, if necessary
  if(is.null(group_names)){

    group_names <- names(std_sol)

  }

  # Some checks
    for(i in seq_along(std_sol)){

    if(any(is.na(std_sol[[i]][["lambda"]]))){
      cli::cli_abort("Some loadings are {.val NA} or {.val NaN}; no omegas are computed.",
                     class = "efa_omega_na_loadings")
    }

    if(any(diag(std_sol[[i]][["theta"]]) <= 0) || any(diag(std_sol[[i]][["psi"]]) <= 0)){
      cli::cli_abort("A Heywood case was detected (a variance of 0 or negative); no omegas are computed.",
                     class = "efa_omega_heywood")
    }

    if(ncol(std_sol[[i]][["lambda"]]) == 1){

      fac_names <- colnames(std_sol[[1]][["lambda"]])
      var_names <- rownames(std_sol[[1]][["lambda"]])

      # Extract sum of factor loadings and error variances
      sum_g <- sum(std_sol[[i]][["lambda"]][, fac_names])
      sum_e <- sum(diag(std_sol[[i]][["theta"]]))

      if(isTRUE(add_ind)){

        if(i == 1){

        cli::cli_inform(
          c("i" = "The model contained a single factor; only omega total and the H index are returned."),
          class = "efa_omega_single_factor"
        )

      }

      # Compute omega and H index
        omega_g <- (sum_g^2) / (sum_g^2 + sum_e)
        h_g <- 1 / (1 + (1 / sum(std_sol[[i]][["lambda"]][, fac_names]^2 /
                                   (1 - std_sol[[i]][["lambda"]][, fac_names]^2))))
        omegas[[i]] <- c(omega_g, h_g)
        names(omegas[[i]]) <- c("Omega", "H")

      } else {

        if(i == 1){

          cli::cli_inform(
            c("i" = "The model contained a single factor; only omega total is returned."),
            class = "efa_omega_single_factor"
          )

        }

        # Compute omega
        omegas[[i]] <- (sum_g^2) / (sum_g^2 + sum_e)

      }

    } else {

      # Extract colnames to see if g_name is in there
      col_names <- colnames(std_sol[[1]][["lambda"]])

      if(!any(col_names %in% g_name)){
        cli::cli_abort(
          c("Could not find the specified general-factor name in the lavaan solution.",
            "i" = "Please check the spelling."),
          class = "efa_omega_g_name"
        )
      }

      if(i == 1){

        if(all(std_sol[[1]][["lambda"]][, g_name] == 0)){

          higherorder <- TRUE

          if(sum(colSums(std_sol[[1]][["beta"]]) > 0) > 1){
            cli::cli_abort("The higher-order model has more than two latent strata or more than one second-order factor; only second-order models with one second-order factor are supported.",
                           class = "efa_omega_higher_order")
          }

          cli::cli_inform(
            c("i" = "The specified general factor is a second-order factor; omegas are computed on the Schmid-Leiman transformed second-order solution."),
            class = "efa_omega_g_second_order"
          )

        } else {

      bi_check <- std_sol[[1]][["lambda"]]
      bi_check[abs(bi_check) > 0 + .Machine$double.eps * 100] <- 1

      if(all(rowSums(bi_check) < 2)){

        cli::cli_abort(
          c("The lavaan input is invalid; no omegas are computed.",
            "i" = "Provide a bifactor model, a second-order model, or a single-factor model."),
          class = "efa_omega_invalid_lavaan"
        )

      }

      if(!all(rowSums(bi_check) > 1)){

        cli::cli_inform(
          c("i" = "Some variables have fewer than two loadings; did you enter a bifactor model? Provide a bifactor model, a second-order model, or a single-factor model."),
          class = "efa_omega_few_loadings"
        )

      }

    }

  }

      # Create list with factor and corresponding subtest names
      col_names <- col_names[!col_names %in% g_name]
      fac_names <- c(g_name, col_names)

      var_names <- list()
      for(j in seq_along(col_names)){
        temp0 <- abs(std_sol[[1]][["lambda"]][, col_names[j]]) > 0 + .Machine$double.eps * 100
        var_names[[j]] <- names(temp0)[temp0]
      }

      # Do SL-transformation if the model is a higher-order model
      if(isTRUE(higherorder)){

        # Calculate direct g loadings
        std_sol[[i]][["lambda"]][, g_name] <- std_sol[[i]][["lambda"]][, col_names] %*% std_sol[[i]][["beta"]][col_names, g_name]

        # Calculate direct group factor loadings (see .sl_group_loadings).
        std_sol[[i]][["lambda"]][, col_names] <- .sl_group_loadings(
          std_sol[[i]][["lambda"]][, col_names], std_sol[[i]][["psi"]], col_names)

      }

      # Create all sums of factor loadings for each factor
      sums_all <- vector("double", length(fac_names))
      for (j in seq_along(fac_names)){
        sums_all[j] <- sum(std_sol[[i]][["lambda"]][, fac_names[j]])
      }

      # Extract g-loadings and error variances, sum of g-loadings for g and sum of
      # respective loadings for every factor
      g_load <- std_sol[[i]][["lambda"]][, g_name]
      e_load <- diag(std_sol[[i]][["theta"]])
      sum_e <- sum(e_load)
      sum_g <- sums_all[1]
      sums_s <- sums_all[2:length(fac_names)]

      # Compute sums of error variances and g-loadings for group factors
      sums_e_s <- vector("double", length(var_names))
      sums_g_s <- vector("double", length(var_names))

      for (j in seq_along(var_names)){
        sums_e_s[j] <- sum(e_load[var_names[[j]]])
        sums_g_s[j] <- sum(g_load[var_names[[j]]])
      }

      # Compute omega total, hierarchical, and subscale for g-factor
      omega_tot_g <- (sum_g^2 + sum(sums_s^2)) / (sum_g^2 + sum(sums_s^2) + sum_e)
      omega_h_g <- sum_g^2 / (sum_g^2 + sum(sums_s^2) +
                                          sum_e)
      omega_sub_g <- sum(sums_s^2) / (sum_g^2 + sum(sums_s^2) + sum_e)

      # Compute omega total, hierarchical, and subscale for group factors
      omega_tot_sub <- (sums_s^2 + sums_g_s^2) / (sums_g_s^2 + sums_s^2 +
                                                         sums_e_s)
      omega_h_sub <- sums_g_s^2 / (sums_g_s^2 + sums_s^2 + sums_e_s)
      omega_sub_sub <- sums_s^2 / (sums_g_s^2 + sums_s^2 + sums_e_s)

      # Combine and display results in a table
      omega_tot <- c(omega_tot_g, omega_tot_sub)
      omega_h <- c(omega_h_g, omega_h_sub)
      omega_sub <- c(omega_sub_g, omega_sub_sub)

      if(isTRUE(add_ind)){

        # Compute H index, ECV, and PUC
        sum_s_load_sq <- vector("double", length(var_names))
        h_s_load <- vector("double", length(var_names))

        # Calculate denominators for H for group factors and for ECV
        for (j in seq_along(var_names)){
          sum_s_load_sq[j] <- sum(std_sol[[i]][["lambda"]][var_names[[j]],
                                                                fac_names[j+1]]^2 /
                                         (1 - std_sol[[i]][["lambda"]][var_names[[j]],
                                                                       fac_names[j+1]]^2))
          h_s_load[j] <- sum(std_sol[[i]][["lambda"]][var_names[[j]],
                                                           fac_names[j+1]]^2)
        }

        h_g <- 1 / (1 + (1 / sum(g_load^2 / (1 - g_load^2))))
        h_s <- 1 / (1 + (1 / sum_s_load_sq))

        ECV <- sum(g_load^2) / sum(sum(g_load^2), sum(h_s_load))

        # Calculate number of items per factor
        item_on_fac <- colSums(std_sol[[i]][["lambda"]] != 0)

        # Calculate number of contaminated correlations
        cont_corrs <- sum(sapply(item_on_fac, function (x) {x*(x - 1) / 2})) -
          (nrow(std_sol[[i]][["lambda"]])*(nrow(std_sol[[i]][["lambda"]]) - 1) / 2)

        # Calculate proportion of uncontaminated correlations (PUC)
        PUC <- 1 - cont_corrs /
          (nrow(std_sol[[i]][["lambda"]])*(nrow(std_sol[[i]][["lambda"]]) - 1) / 2)

        # Create output
        h <- c(h_g, h_s)

        omegas[[i]] <- cbind(omega_tot, omega_h, omega_sub, h, NA, NA)
        omegas[[i]][1, 5] <- ECV
        omegas[[i]][1, 6] <- PUC
        colnames(omegas[[i]]) <- c("tot", "hier", "sub", "H", "ECV", "PUC")

      } else {

        # Create output
        omegas[[i]] <- cbind(omega_tot, omega_h, omega_sub)
        colnames(omegas[[i]]) <- c("tot", "hier", "sub")

      }

      rownames(omegas[[i]]) <- fac_names

  }

  }

  if(length(std_sol) > 1){

    names(omegas) <- group_names

  } else {

    omegas <- omegas[[1]]

  }

  class(omegas) <- "OMEGA"

  return(omegas)

}
