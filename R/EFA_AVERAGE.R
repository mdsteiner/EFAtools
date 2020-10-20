EFA_AVERAGE <- function(x, n_factors, N = NA, method = "PAF", rotation = "promax",
                        type = "none", aggregation = c("mean", "median"), trim = 0,
                        salience_threshold = .3,
                        max_iter = 1e4,
                        init_comm = c("smc", "mac", "unity"),
                        criterion = c(1e-3),
                        criterion_type = c("sum", "max_individual"),
                        abs_eigen = c(TRUE),
                        varimax_type = c("svd", "kaiser"),
                        normalize = TRUE,
                        k_promax = 2:4, k_simplimax = ncol(x),
                        P_type = c("norm", "unnorm"), precision = 1e-5,
                        start_method = c("psych", "factanal"),
                        use = c("pairwise.complete.obs", "all.obs",
                                "complete.obs", "everything", "na.or.complete"),
                        cor_method = c("pearson", "spearman", "kendall"),
                        show_progress = TRUE) {


  # x= IDS2_R
  # n_factors = 6
  # N = 1991
  # method = "PAF"
  # rotation = "promax"
  # type = "none"
  # aggregation = "mean"
  # trim = 0
  # salience_threshold = .3
  # max_iter = 1e4
  # init_comm = c("smc", "mac", "unity")
  # criterion = c(1e-3, 1e-6)
  # criterion_type = c("sum", "max_individual")
  # abs_eigen = c(TRUE, FALSE)
  # varimax_type = c("svd", "kaiser")
  # normalize = TRUE
  # k_promax = 2:4
  # k_simplimax = ncol(x)
  # P_type = c("norm", "unnorm")
  # precision = 1e-5
  # start_method = c("psych", "factanal")
  # use = "pairwise.complete.obs"
  # cor_method = "pearson"
  # show_progress = TRUE

  # Perform argument checks
  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n"))

  }

  checkmate::assert_count(n_factors)
  checkmate::assert_count(N, na.ok = TRUE)
  checkmate::assert_subset(method, c("PAF", "ML", "ULS"),
                           empty.ok = FALSE)
  checkmate::assert_subset(rotation, c("none", "orthogonal", "oblique", "varimax",
                                       "equamax", "quartimax", "geominT", "bentlerT",
                                       "bifactorT", "promax", "oblimin", "quartimin",
                                       "simplimax", "bentlerQ", "geominQ", "bifactorQ"),
                           empty.ok = FALSE)
  checkmate::assert_subset(type, c("none", "EFAtools", "psych", "SPSS"),
                           empty.ok = FALSE)
  aggregation <- match.arg(aggregation)
  checkmate::assert_number(trim, lower = 0, upper = 0.5)
  checkmate::assert_number(salience_threshold, lower = 0, upper = 1)
  checkmate::assert_count(max_iter)
  checkmate::assert_subset(init_comm, c("smc", "mac", "unity"),
                           empty.ok = FALSE)
  checkmate::assert_vector(criterion, strict = TRUE, any.missing = FALSE,
                           min.len = 1)
  checkmate::assert_true(all(criterion > 0 & criterion < 1))
  checkmate::assert_subset(criterion_type, c("max_individual", "sum"),
                           empty.ok = FALSE)
  checkmate::assert_subset(abs_eigen, c(TRUE, FALSE),
                           empty.ok = FALSE)
  checkmate::assert_subset(varimax_type, c("svd", "kaiser"),
                           empty.ok = FALSE)
  checkmate::assert_subset(normalize, c(TRUE, FALSE),
                           empty.ok = FALSE)
  checkmate::assert_vector(k_promax, strict = TRUE, any.missing = FALSE,
                           min.len = 1)
  checkmate::assert_vector(k_simplimax, strict = TRUE, any.missing = FALSE,
                           min.len = 1)
  checkmate::assert_subset(P_type, c("unnorm", "norm"),
                           empty.ok = FALSE)
  checkmate::assert_vector(precision, strict = TRUE, any.missing = FALSE,
                           min.len = 1)
  checkmate::assert_true(all(precision > 0 & precision < 1))
  checkmate::assert_subset(start_method, c("psych", "factanal"),
                           empty.ok = FALSE)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  checkmate::assert_flag(show_progress)



  # Check if it is a correlation matrix
  if(.is_cormat(x)){

    R <- x

  } else {

    message(cli::col_cyan(cli::symbol$info, " 'x' was not a correlation matrix. Correlations are found from entered raw data.\n"))

    if (!is.na(N)) {
      warning(crayon::yellow$bold("!"), crayon::yellow(" 'N' was set and data entered. Taking N from data.\n"))
    }

    R <- stats::cor(x, use = use, method = cor_method)
    colnames(R) <- colnames(x)
    N <- nrow(x)

  }

  # Check if correlation matrix is invertable, if it is not, stop with message
  R_i <- try(solve(R), silent = TRUE)

  if (inherits(R_i, "try-error")) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Correlation matrix is singular, no further analyses are performed\n"))
  }

  # Check if correlation matrix is positive definite, if it is not,
  # smooth the matrix (cor.smooth throws a warning)
  if (any(eigen(R, symmetric = TRUE, only.values = TRUE)$values <= 0)) {

    R <- psych::cor.smooth(R)

  }

  # Check if model is identified

  # calculate degrees of freedom
  m <- ncol(R)
  df <- ((m - n_factors)**2 - (m + n_factors)) / 2

  if(df < 0){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" The model is underidentified. Please enter a lower number of factors or use a larger number of indicators and try again.\n"))

  } else if (df == 0){

    warning(crayon::yellow$bold("!"), crayon::yellow(" The model is just identified (df = 0). We suggest to try again with a lower number of factors or a larger number of indicators.\n"))

  }

  ### create the grid with all combinations of the input arguments

  grid_list <- list()


  if ("PAF" %in% method) {

    if ("EFAtools" %in% type) {
      grid_list[["ftls_pf"]] <- .type_grid(method = "PAF", init_comm = "smc",
                                           criterion = 1e-3, criterion_type = "sum",
                                           abs_eigen = TRUE, start_method = NA,
                                           rotation = rotation, k_promax = 4,
                                           normalize = TRUE, P_type = "norm",
                                           precision = 1e-5,
                                           varimax_type = "svd",
                                           k_simplimax = k_simplimax)
    }

    if ("psych" %in% type) {
      grid_list[["psch_pf"]] <- .type_grid(method = "PAF", init_comm = "smc",
                                           criterion = 1e-3, criterion_type = "sum",
                                           abs_eigen = FALSE, start_method = NA,
                                           rotation = rotation, k_promax = 4,
                                           normalize = TRUE, P_type = "unnorm",
                                           precision = 1e-5,
                                           varimax_type = "svd",
                                           k_simplimax = k_simplimax)
    }

    if ("SPSS" %in% type) {
      grid_list[["spss_pf"]] <- .type_grid(method = "PAF", init_comm = "smc",
                                           criterion = 1e-3, criterion_type = "max_individual",
                                           abs_eigen = TRUE, start_method = NA,
                                           rotation = rotation, k_promax = 4,
                                           normalize = TRUE, P_type = "norm",
                                           precision = 1e-5,
                                           varimax_type = "kaiser",
                                           k_simplimax = k_simplimax)
    }

    if ("none" %in% type) {
      if(any(is.na(init_comm), is.na(criterion), is.na(criterion_type),
             is.na(abs_eigen))) {

        stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" type = 'none' and method = 'PAF' is used but at least one of init_comm, criterion, criterion_type, and abs_eigen is not specified. Please specify all these arguments or use 'EFAtools', 'psych', or 'SPSS' as type.\n"))

      } else if (any(c("promax", "oblique") %in% rotation) &&
                 any(is.na(k_promax), is.na(normalize), is.na(P_type),
                     is.na(precision), is.na(varimax_type))) {

        stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" type = 'none' and rotation = 'promax' or rotation = 'oblique' is used but at least one of k_promax, normalize, P_type, precision, and varimax_type is not specified. Please specify all these arguments or use 'EFAtools', 'psych', or 'SPSS' as type.\n"))


      } else if (any(c("simplimax", "oblique") %in% rotation) &&
                  any(is.na(k_simplimax), is.na(normalize),
                      is.na(precision))) {

        stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" type = 'none' and rotation = 'simplimax' or rotation = 'oblique' is used but at least one of k_simplimax, normalize, and precision is not specified. Please specify all these arguments or use 'EFAtools', 'psych', or 'SPSS' as type.\n"))


      } else if (any(c("varimax", "orthogonal") %in% rotation) &&
                  any(is.na(varimax_type), is.na(normalize),
                      is.na(precision))) {

        stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" type = 'none' and rotation = 'varimax' or rotation = 'orthogonal' is used but at least one of varimax_type, normalize, and precision is not specified. Please specify all these arguments or use 'EFAtools', 'psych', or 'SPSS' as type.\n"))


      } else {
        grid_list[["nn_pf"]] <- .type_grid(method = "PAF", init_comm = init_comm,
                                           criterion = criterion, criterion_type = criterion_type,
                                           abs_eigen = abs_eigen, start_method = NA,
                                           rotation = rotation, k_promax = k_promax,
                                           normalize = normalize, P_type = P_type,
                                           precision = precision,
                                           varimax_type = varimax_type,
                                           k_simplimax = k_simplimax)
      }

    }
  }

    if ("ML" %in% method) {

      if ("EFAtools" %in% type) {
        grid_list[["ftls_ml"]] <- .type_grid(method = "ML", init_comm = NA,
                                             criterion = NA, criterion_type = NA,
                                             abs_eigen = NA, start_method = "psych",
                                             rotation = rotation, k_promax = 4,
                                             normalize = TRUE, P_type = "norm",
                                             precision = 1e-5,
                                             varimax_type = "svd",
                                             k_simplimax = k_simplimax)
      }

      if ("psych" %in% type) {
        grid_list[["psch_ml"]] <- .type_grid(method = "ML", init_comm = NA,
                                             criterion = NA, criterion_type = NA,
                                             abs_eigen = NA, start_method = "psych",
                                             rotation = rotation, k_promax = 4,
                                             normalize = TRUE, P_type = "unnorm",
                                             precision = 1e-5,
                                             varimax_type = "svd",
                                             k_simplimax = k_simplimax)
      }

      if ("SPSS" %in% type) {
        grid_list[["spss_ml"]] <- .type_grid(method = "ML", init_comm = NA,
                                             criterion = NA, criterion_type = NA,
                                             abs_eigen = NA, start_method = "psych",
                                             rotation = rotation, k_promax = 4,
                                             normalize = TRUE, P_type = "norm",
                                             precision = 1e-5,
                                             varimax_type = "kaiser",
                                             k_simplimax = k_simplimax)
      }

      if ("none" %in% type) {
        if(any(is.na(start_method))) {

          stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" type = 'none' and method = 'ML' is used, but  start_method is not specified. Please specify all these arguments or use 'EFAtools', 'psych', or 'SPSS' as type.\n"))

        } else if (any(c("promax", "oblique") %in% rotation) &&
                   any(is.na(k_promax), is.na(normalize), is.na(P_type),
                       is.na(precision), is.na(varimax_type))) {

          stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" type = 'none' and rotation = 'promax' or rotation = 'oblique' is used but at least one of k_promax, normalize, P_type, precision, and varimax_type is not specified. Please specify all these arguments or use 'EFAtools', 'psych', or 'SPSS' as type.\n"))


        } else if (any(c("simplimax", "oblique") %in% rotation) &&
                   any(is.na(k_simplimax), is.na(normalize),
                       is.na(precision))) {

          stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" type = 'none' and rotation = 'simplimax' or rotation = 'oblique' is used but at least one of k_simplimax, normalize, and precision is not specified. Please specify all these arguments or use 'EFAtools', 'psych', or 'SPSS' as type.\n"))


        } else if (any(c("varimax", "orthogonal") %in% rotation) &&
                   any(is.na(varimax_type), is.na(normalize),
                       is.na(precision))) {

          stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" type = 'none' and rotation = 'varimax' or rotation = 'orthogonal' is used but at least one of varimax_type, normalize, and precision is not specified. Please specify all these arguments or use 'EFAtools', 'psych', or 'SPSS' as type.\n"))


        } else {
          grid_list[["nn_ml"]] <- .type_grid(method = "ML", init_comm = NA,
                                             criterion = NA, criterion_type = NA,
                                             abs_eigen = NA, start_method = start_method,
                                             rotation = rotation, k_promax = k_promax,
                                             normalize = normalize, P_type = P_type,
                                             precision = precision,
                                             varimax_type = varimax_type,
                                             k_simplimax = k_simplimax)
        }

      }
    }

      if ("ULS" %in% method) {

        if ("EFAtools" %in% type) {
          grid_list[["ftls_ls"]] <- .type_grid(method = "ULS", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, start_method = NA,
                                               rotation = rotation, k_promax = 4,
                                               normalize = TRUE, P_type = "norm",
                                               precision = 1e-5,
                                               varimax_type = "svd",
                                               k_simplimax = k_simplimax)
        }

        if ("psych" %in% type) {
          grid_list[["psch_ls"]] <- .type_grid(method = "ML", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, start_method = NA,
                                               rotation = rotation, k_promax = 4,
                                               normalize = TRUE, P_type = "unnorm",
                                               precision = 1e-5,
                                               varimax_type = "svd",
                                               k_simplimax = k_simplimax)
        }

        if ("SPSS" %in% type) {
          grid_list[["spss_ls"]] <- .type_grid(method = "ML", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, start_method = NA,
                                               rotation = rotation, k_promax = 4,
                                               normalize = TRUE, P_type = "norm",
                                               precision = 1e-5,
                                               varimax_type = "kaiser",
                                               k_simplimax = k_simplimax)
        }

        if ("none" %in% type) {
          if (any(c("promax", "oblique") %in% rotation) &&
                     any(is.na(k_promax), is.na(normalize), is.na(P_type),
                         is.na(precision), is.na(varimax_type))) {

            stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" type = 'none' and rotation = 'promax' or rotation = 'oblique' is used but at least one of k_promax, normalize, P_type, precision, and varimax_type is not specified. Please specify all these arguments or use 'EFAtools', 'psych', or 'SPSS' as type.\n"))


          } else if (any(c("simplimax", "oblique") %in% rotation) &&
                     any(is.na(k_simplimax), is.na(normalize),
                         is.na(precision))) {

            stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" type = 'none' and rotation = 'simplimax' or rotation = 'oblique' is used but at least one of k_simplimax, normalize, and precision is not specified. Please specify all these arguments or use 'EFAtools', 'psych', or 'SPSS' as type.\n"))


          } else if (any(c("varimax", "orthogonal") %in% rotation) &&
                     any(is.na(varimax_type), is.na(normalize),
                         is.na(precision))) {

            stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" type = 'none' and rotation = 'varimax' or rotation = 'orthogonal' is used but at least one of varimax_type, normalize, and precision is not specified. Please specify all these arguments or use 'EFAtools', 'psych', or 'SPSS' as type.\n"))


          } else {
            grid_list[["nn_ls"]] <- .type_grid(method = "ML", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, start_method = NA,
                                               rotation = rotation, k_promax = k_promax,
                                               normalize = normalize, P_type = P_type,
                                               precision = precision,
                                               varimax_type = varimax_type,
                                               k_simplimax = k_simplimax)
          }

        }
      }

  arg_grid <- unique(do.call(rbind, grid_list))

  ### Run all efas

  if (nrow(arg_grid) == 1) {

    warning(crayon::yellow$bold("!"), crayon::yellow(" There was only one combination of arguments, returning normal EFA() output.\n"))

    return(EFA(R, n_factors, N = N, method = arg_grid$method, rotation = arg_grid$rotation,
        type = "none", max_iter = max_iter, init_comm = arg_grid$init_comm,
        criterion = arg_grid$criterion, criterion_type = arg_grid$criterion_type,
        abs_eigen = arg_grid$abs_eigen, varimax_type = arg_grid$varimax_type,
        k = ifelse(arg_grid$rotation == "promax", arg_grid$k_promax, arg_grid$k_simplimax),
        normalize = arg_grid$normalize, P_type = arg_grid$P_type, precision = precision,
        order_type = "eigen", start_method = arg_grid$start_method))

  }

  n_cores <- future::nbrOfWorkers()

  if (isTRUE(show_progress)) {
    progressr::handlers("progress")
  }

  progressr::with_progress({
    n_efa <- nrow(arg_grid)
    stepsize <- round(n_efa / 100 * 10)
    efa_progress_bar <- progressr::progressor(steps = n_efa / stepsize)
    efa_progress_bar("Running EFAs:", class = "sticky", amount = 0)
    efa_list <- future.apply::future_lapply(1:n_efa,
                                            function(i, methods, rotations,
                                                     init_comms, criteria,
                                                     criterion_types, abs_eigens,
                                                     varimax_types, k_ps, k_ss,
                                                     normalizes, P_types, start_methods) {
      if (i %% stepsize == 0){
        efa_progress_bar(message = sprintf("Running EFA %g of %g", i, n_efa))
      }

      try(EFA(R, n_factors, N = N, method = methods[i], rotation = rotations[i],
              type = "none", max_iter = max_iter, init_comm = init_comms[i],
              criterion = criteria[i], criterion_type = criterion_types[i],
              abs_eigen = abs_eigens[i], varimax_type = varimax_types[i],
              k = ifelse(rotations[i] == "promax",k_ps[i], k_ss[i]),
              normalize = normalizes[i], P_type = P_types[i], precision = precision,
              order_type = "eigen", start_method = start_methods[i]),
          silent = TRUE)
    }, methods = arg_grid$method, rotations = arg_grid$rotation,
    init_comms = arg_grid$init_comm, criteria = arg_grid$criterion,
    criterion_types = arg_grid$criterion_type, abs_eigens = arg_grid$abs_eigen,
    varimax_types = arg_grid$varimax_type, k_ps = arg_grid$k_promax,
    k_ss = arg_grid$k_simplimax, normalizes = arg_grid$normalize,
    P_types = arg_grid$P_type, start_methods = arg_grid$start_method)

    if (n_efa %% 10 != 0){
      efa_progress_bar(message = "Done Running EFAs ")
    }
  })

  ### Extract relevant information from EFA outputs
  if (isTRUE(show_progress)) {
    .show_agg_progress("\U0001f3c3", "Extracting data...")
  }
  ext_list <- .extract_data(efa_list, R, n_factors, n_efa, rotation, salience_threshold)

  if (isTRUE(show_progress)) {
    .show_agg_progress("\U0001f6b6", "Reordering factors...")
  }

  re_list <- .array_reorder(ext_list$vars_accounted, ext_list$L, ext_list$L_corres,
                            ext_list$phi, ext_list$extract_phi, n_factors)

  if (isTRUE(show_progress)) {
    .show_agg_progress("\U0001f3c3", "Aggregating data...")
  }
  agg_list <- suppressWarnings(
    .aggregate_values(re_list$vars_accounted, re_list$L, re_list$L_corres, ext_list$h2, re_list$phi,
                      ext_list$extract_phi, aggregation, trim,
                      ext_list$for_grid[,c("chisq", "p_chi", "caf", "cfi",
                                           "rmsea", "aic", "bic")], df, colnames(R)))

  arg_grid <- cbind(arg_grid, ext_list$for_grid)


  settings <- list(
    method = method,
    rotation = rotation,
    type = type,
    n_factors = n_factors,
    N = N,
    init_comm = init_comm,
    criterion = criterion,
    criterion_type = criterion_type,
    abs_eigen = abs_eigen,
    varimax_type = varimax_type,
    normalize = normalize,
    k_promax = k_promax,
    k_simplimax = k_simplimax,
    P_type = P_type,
    precision = precision,
    start_method = start_method,
    use = use,
    cor_method = cor_method,
    max_iter = max_iter,
    aggregation = aggregation,
    trim = trim,
    salience_threshold = salience_threshold
  )

  # Create output
  output <- list(
    orig_R = R,
    h2 = agg_list$h2,
    loadings = agg_list$loadings,
    Phi = agg_list$phi,
    ind_fac_corres = agg_list$ind_fac_corres,
    vars_accounted = agg_list$vars_accounted,
    fit_indices = agg_list$fit_indices,
    implementations_grid = arg_grid,
    efa_list = efa_list,
    settings = settings
  )

  class(output) <- "EFA_AVERAGE"

  if (isTRUE(show_progress)) {
    .show_agg_progress("", "", done = TRUE)
  }

  return(output)

}

