# Helpers for averaging many EFA solutions: build the per-rotation parameter grid, extract
# and align factors across runs by congruence, and aggregate loadings, communalities,
# variances, and fit indices.

# for progress bar in EFA_AVERAGE
.show_av_progress <- function(emoji, what, done = FALSE) {

  cat("\r", rep(" ", ifelse(getOption("width") > 30, getOption("width"), 30)))
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



  if (any(rotation %in% c(.oblq_rotations, "oblique"))) {
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
  srmr <- rep(NA_real_, n_efa)
  tli <- rep(NA_real_, n_efa)
  ecvi <- rep(NA_real_, n_efa)
  rmsr <- rep(NA_real_, n_efa)

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
          # Use the fit's Heywood flag so improper solutions are excluded
          # consistently across estimators (a communality >= 1, or an ML/ULS
          # uniqueness pinned at the estimation boundary).
          has_heywood <- length(efa_temp$heywood) > 0
          heywood[row_i] <- has_heywood

          if (!has_heywood) {

            aic[row_i] <- efa_temp$fit_indices$AIC
            bic[row_i] <- efa_temp$fit_indices$BIC
            chisq[row_i] <- efa_temp$fit_indices$chi
            p_chi[row_i] <- efa_temp$fit_indices$p_chi
            caf[row_i] <- efa_temp$fit_indices$CAF
            rmsea[row_i] <- efa_temp$fit_indices$RMSEA
            cfi[row_i] <- efa_temp$fit_indices$CFI
            srmr[row_i] <- efa_temp$fit_indices$SRMR
            tli[row_i] <- efa_temp$fit_indices$TLI
            ecvi[row_i] <- efa_temp$fit_indices$ECVI
            rmsr[row_i] <- efa_temp$fit_indices$RMSR

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
      bic= bic,
      srmr = srmr,
      tli = tli,
      ecvi = ecvi,
      rmsr = rmsr
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

      # Match each later solution's factors to the first via the optimal
      # (linear-sum-assignment) congruence alignment, which guarantees a true
      # permutation and applies sign flips consistently to loadings and Phi.
      aligned <- .align_solution(
        L_target = L1, L = L[,, efa_i],
        Phi = if (isTRUE(extract_phi)) phi[,, efa_i] else NULL)

      factor_order <- aligned$factor_order

      # The correspondence indicators and variances are sign-invariant, so they
      # only need the same column reordering.
      L[,, efa_i] <- aligned$loadings
      L_corres[,, efa_i] <- L_corres[,, efa_i][, factor_order]
      vars_accounted[,, efa_i] <- vars_accounted[,, efa_i][, factor_order]

      if (isTRUE(extract_phi)) {
        phi[,, efa_i] <- aligned$Phi
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

      cli::cli_abort(
        c("{.code rotation = \"none\"} was used but {.arg rotation} has length > 1.",
          "i" = "Can only average EFAs with rotations of the same type: {.val none}, {.val orthogonal}, or {.val oblique}."),
        class = "efa_rotation_length"
      )

    }
  } else if ("oblique" %in% rotation) {
    if (length(rotation) == 1) {

      t_grid_list[["blq"]] <- .oblq_grid(method = method, init_comm = init_comm,
                                         criterion = criterion, criterion_type = criterion_type,
                                         abs_eigen = abs_eigen, start_method = start_method,
                                         rotation = .oblq_rotations,
                                         k_promax = k_promax, normalize = normalize,
                                         P_type = P_type, precision = precision,
                                         varimax_type = varimax_type, k_simplimax = k_simplimax)

    } else {

      cli::cli_abort(
        c("{.code rotation = \"oblique\"} was used but {.arg rotation} has length > 1.",
          "i" = "Can only average EFAs with rotations of the same type: {.val none}, {.val orthogonal}, or {.val oblique}."),
        class = "efa_rotation_length"
      )

    }

  } else if ("orthogonal" %in% rotation) {

    if (length(rotation) == 1) {

      t_grid_list[["rth"]] <- .orth_grid(method = method, init_comm = init_comm,
                                         criterion = criterion, criterion_type = criterion_type,
                                         abs_eigen = abs_eigen, start_method = start_method,
                                         rotation = .orth_rotations,
                                         normalize = normalize, precision = precision,
                                         varimax_type = varimax_type)

    } else {

      cli::cli_abort(
        c("{.code rotation = \"orthogonal\"} was used but {.arg rotation} has length > 1.",
          "i" = "Can only average EFAs with rotations of the same type: {.val none}, {.val orthogonal}, or {.val oblique}."),
        class = "efa_rotation_length"
      )

    }

  } else if (all(rotation %in% .oblq_rotations)) {

    t_grid_list[["blq2"]] <- .oblq_grid(method = method, init_comm = init_comm,
                                        criterion = criterion, criterion_type = criterion_type,
                                        abs_eigen = abs_eigen, start_method = start_method,
                                        rotation = rotation,
                                        k_promax = k_promax, normalize = normalize,
                                        P_type = P_type, precision = precision,
                                        varimax_type = varimax_type, k_simplimax = k_simplimax)

  } else if (all(rotation %in% .orth_rotations)) {

    t_grid_list[["rth2"]] <- .orth_grid(method = method, init_comm = init_comm,
                                        criterion = criterion, criterion_type = criterion_type,
                                        abs_eigen = abs_eigen, start_method = start_method,
                                        rotation = rotation,
                                        normalize = normalize, precision = precision,
                                        varimax_type = varimax_type)

  } else if (any(rotation %in% .oblq_rotations) &&
             any(rotation %in% .orth_rotations)) {
    cli::cli_abort(
      c("{.arg rotation} mixes oblique and orthogonal rotations, but only rotations of the same kind can be averaged.",
        "*" = "Oblique rotations: {.val {(.oblq_rotations)}}.",
        "*" = "Orthogonal rotations: {.val {(.orth_rotations)}}."),
      class = "efa_rotation_mismatch"
    )
  }

  return(do.call(rbind, t_grid_list))
}

### assemble the full implementation grid for EFA_AVERAGE
# Stack the per-(method, type) parameter grids into one grid of EFA() argument
# combinations. The three named types (EFAtools/psych/SPSS) take their tuning
# arguments from the shared `.efa_presets` table (the same source EFA() resolves),
# so an averaged sub-analysis matches EFA(type = <that type>); type "none" forwards
# the user-supplied arguments. Each estimator keeps only the arguments it uses
# (PAF: starting communalities and convergence criterion; ML: starting values),
# mirroring EFA() itself. Methods are processed in the order PAF, ML, ULS and types
# in the order EFAtools, psych, SPSS, none, so the row order and names of the
# assembled grid are stable.
.build_avg_grid <- function(method, type, rotation, init_comm, criterion,
                            criterion_type, abs_eigen, start_method, k_promax,
                            normalize, P_type, precision, varimax_type,
                            k_simplimax) {

  # Tags used to build the grid-list names (e.g. "ftls_pf", "nn_ml"), kept
  # identical to the names the grid is assembled under so row names stay stable.
  method_tag <- c(PAF = "pf", ML = "ml", ULS = "ls")
  type_tag   <- c(EFAtools = "ftls", psych = "psch", SPSS = "spss", none = "nn")

  # The communality/start arguments each estimator uses; the rest stay NA. For
  # named types PAF's values come from the `.efa_presets` PAF block, for type
  # "none" from the user.
  method_args <- function(m, from_user, paf = NULL) {
    switch(m,
      PAF = if (from_user) {
        list(init_comm = init_comm, criterion = criterion,
             criterion_type = criterion_type, abs_eigen = abs_eigen,
             start_method = NA)
      } else {
        list(init_comm = paf$init_comm, criterion = paf$criterion,
             criterion_type = paf$criterion_type,
             abs_eigen = paf$abs_eigen, start_method = NA)
      },
      ML  = list(init_comm = NA, criterion = NA, criterion_type = NA,
                 abs_eigen = NA,
                 start_method = if (from_user) start_method else "psych"),
      ULS = list(init_comm = NA, criterion = NA, criterion_type = NA,
                 abs_eigen = NA, start_method = NA))
  }

  grid_list <- list()

  for (m in c("PAF", "ML", "ULS")) {
    if (!(m %in% method)) next

    for (t in c("EFAtools", "psych", "SPSS", "none")) {
      if (!(t %in% type)) next

      if (t == "none") {
        ma   <- method_args(m, from_user = TRUE)
        rest <- list(k_promax = k_promax, normalize = normalize, P_type = P_type,
                     precision = precision, varimax_type = varimax_type)
      } else {
        # Resolve named-type tuning arguments through .resolve_settings() -- the
        # same resolver EFA() uses -- with EFA()'s default user arguments
        # (normalize on, everything else unset), so an averaged sub-analysis
        # cannot drift from EFA(type = t): PAF's convergence settings come from
        # the PAF block, the promax/varimax settings from the PROMAX block.
        # `precision` is not a preset, so it keeps EFA()'s default.
        paf <- .resolve_settings(
          type = t,
          user = list(init_comm = NA, criterion = NA, criterion_type = NA,
                      max_iter = NA, abs_eigen = NA),
          preset = .efa_presets$PAF)
        pro <- .resolve_settings(
          type = t,
          user = list(normalize = TRUE, P_type = NA, order_type = NA,
                      varimax_type = NA, k = NA),
          preset = .efa_presets$PROMAX)
        ma   <- method_args(m, from_user = FALSE, paf = paf)
        rest <- list(k_promax = pro$k, normalize = pro$normalize,
                     P_type = pro$P_type, precision = 1e-5,
                     varimax_type = pro$varimax_type)
      }

      grid_list[[paste0(type_tag[[t]], "_", method_tag[[m]])]] <- .type_grid(
        method = m, init_comm = ma$init_comm, criterion = ma$criterion,
        criterion_type = ma$criterion_type, abs_eigen = ma$abs_eigen,
        start_method = ma$start_method, rotation = rotation,
        k_promax = rest$k_promax, normalize = rest$normalize,
        P_type = rest$P_type, precision = rest$precision,
        varimax_type = rest$varimax_type, k_simplimax = k_simplimax)
    }
  }

  unique(do.call(rbind, grid_list))
}
