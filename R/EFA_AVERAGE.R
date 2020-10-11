EFA_AVERAGE <- function(x, n_factors, N = NA, method = c("PAF", "ML", "ULS"),
                          rotation = c("oblique"),
                          type = c("EFAtools", "psych", "SPSS"),
                          max_iter = 1e4, init_comm = NA, criterion = NA,
                          criterion_type = NA, abs_eigen = NA,
                          use = c("pairwise.complete.obs", "all.obs",
                                  "complete.obs", "everything", "na.or.complete"),
                          varimax_type = NA, k_promax = NA, k_simplimax = ncol(x),
                          normalize = TRUE, P_type = NA, precision = 1e-5,
                          order_type = NA, start_method = NA,
                          cor_method = c("pearson", "spearman", "kendall"),
                          aggregation = c("mean", "median"), trim = 0,
                          salience_threshold = .3, ...) {

  # for testing
  x <- EFAtools::DOSPERT$cormat
  n_factors = 5
  N = EFAtools::DOSPERT$N
  method = c("PAF", "ML", "ULS")
  rotation = c("oblique")
  type = c("EFAtools", "psych", "SPSS", "none")
  max_iter = 1e4
  init_comm = c("smc", "mac")
  criterion = c(.001, .0001)
  criterion_type = c("sum", "max_individual")
  abs_eigen = c(TRUE, FALSE)
  use = c("pairwise.complete.obs")
  varimax_type = c("svd", "kaiser")
  k_promax = c(4,3)
  k_simplimax = ncol(x)
  normalize = TRUE
  P_type = c("norm", "unnorm")
  precision = 1e-5
  order_type = "eigen"
  start_method = c("psych", "factanal")
  cor_method = c("pearson")
  aggregation = "mean"
  salience_threshold = .3



  # TODO
  ### Input
  # Same as for EFA
  # Extraction method: possibility for multiple (means that it is averaged about these)
  # PAF implementations (like type, but expanded): EFAtools, R, SPSS, all, evtl. best XX or so
  # Rotation method: all oblique, all orthogonal, or spec. orthogonal or oblique
  # trimmed (passed to mean function)
  # Heywood: true or false (include solutions with heywood cases for averaging or not)
  # averaging method: mean, median



  # Perform argument checks
  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n"))

  }

  #ADAPT INPUT CHECKS

  method <- match.arg(method)
  rotation <- match.arg(rotation)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  type <- match.arg(type)
  start_method <- match.arg(start_method)

  c("none", "orthogonal", "oblique", "varimax", "equamax", "quartimax", "geominT",
  "bentlerT", "bifactorT", "promax", "oblimin",
  "quartimin", "simplimax", "bentlerQ", "geominQ",
  "bifactorQ")

  checkmate::assert_subset(rotation, c("none", "orthogonal", "oblique", "varimax",
                                       "equamax", "quartimax", "geominT", "bentlerT",
                                       "bifactorT", "promax", "oblimin", "quartimin",
                                       "simplimax", "bentlerQ", "geominQ", "bifactorQ"),
                           empty.ok = FALSE)

  checkmate::assert_count(n_factors)
  checkmate::assert_count(N, na.ok = TRUE)
  checkmate::assert_count(max_iter, null.ok = TRUE)
  checkmate::assert_choice(init_comm, c("smc", "mac", "unity"), null.ok = TRUE)
  checkmate::assert_number(criterion, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_choice(criterion_type, c("max_individual", "sums", "sum"),
                           null.ok = TRUE)
  checkmate::assert_flag(abs_eigen, null.ok = TRUE)
  checkmate::assert_number(k, null.ok = TRUE)
  checkmate::assert_choice(varimax_type, c("svd", "kaiser"), null.ok = TRUE)
  checkmate::assert_flag(normalize, null.ok = TRUE)
  checkmate::assert_choice(P_type, c("unnorm", "norm"), null.ok = TRUE)
  checkmate::assert_number(precision, lower = 0, upper = 1)
  checkmate::assert_choice(order_type, c("eigen", "ss_factors"), null.ok = TRUE)

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
                                           abs_eigen = criterion_type, start_method = NA,
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
              order_type = "eigen", start_method = start_methods[i], ...),
          silent = TRUE)
    }, methods = arg_grid$method, rotations = arg_grid$rotation,
    init_comms = arg_grid$init_comm, criteria = arg_grid$criterion,
    criterion_types = arg_grid$criterion_type, abs_eigens = arg_grid$abs_eigen,
    varimax_types = arg_grid$varimax_type, k_ps = arg_grid$k_promax,
    k_ss = arg_grid$k_simplimax, normalizes = arg_grid$normalize,
    P_types = arg_grid$P_type, start_methods = arg_grid$start_method)
  })

  ### Extract relevant information from EFA outputs
  if (isTRUE(show_progress)) {
    .show_agg_progress("Extracting data...")
  }
  ext_list <- .extract_data(efa_list, R, n_factors, n_efa, rotation)

  if (isTRUE(show_progress)) {
    .show_agg_progress("Reordering factors...")
  }

  re_list <- .array_reorder(ext_list$L, ext_list$L_corres, ext_list$phi,
                            ext_list$extract_phi, ext_list$n_factors)

  if (isTRUE(show_progress)) {
    .show_agg_progress("Aggregating data...")
  }
  agg_list <- .aggregate_values(re_list$L, re_list$L_corres, ext_list$h2, re_list$phi,
                                ext_list$extract_phi, aggregation, trim,
                                ext_list$for_grid[[c("chisq", "p_chi", "caf", "cfi",
                                                     "rmsea", "aic", "bic")]])

  arg_grid <- cbind(arg_grid, ext_list$for_grid)


  settings <- list(
    method = method,
    rotation = rotation,
    type = type,
    n_factors = n_factors,
    N = N,
    use = use,
    cor_method = cor_method,
    max_iter = max_iter,
    precision = precision,
    aggregation = aggregation,
    trim = trim,
    salience_threshold = salience_threshold
  )

  # Create output
  output <- list(
    orig_R = R,
    h2 = agg_list$h2,
    loadings = agg_list$loadings,
    phi = agg_list$phi_list,
    ind_fac_corres = agg_list$ind_fac_corres,
    fit_indices = agg_list$fit_indices,
    implementations_grid = arg_grid,
    efa_list = efa_list,
    settings = settings
  )

  class(output) <- "EFA_AVERAGE"

  if (isTRUE(show_progress)) {
    .show_agg_progress("done", done = TRUE)
  }

  return(output)

}

