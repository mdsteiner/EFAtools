#' Model averaging across different EFA methods and types
#'
#' Not all EFA procedures always arrive at the same solution. This function allows
#' you perform a number of EFAs from different methods (e.g., Maximum Likelihood
#' and Principal Axis Factoring), with different implementations (e.g., the SPSS
#' and psych implementations of Principal Axis Factoring), and across different
#' rotations of the same type (e.g., multiple oblique rotations, like promax and
#' oblimin). EFA_AVERAGE will then run all these EFAs and provide a summary across
#' the different solutions.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations. If raw data is entered, the correlation matrix is found from the
#' data.
#' @param n_factors numeric. Number of factors to extract.
#' @param N numeric. The number of observations. Needs only be specified if a
#' correlation matrix is used. If input is a correlation matrix and \code{N} = NA
#' (default), not all fit indices can be computed.
#' @param method character vector. Any combination of  "PAF", "ML", and "ULS",
#' to use principal axis factoring, maximum likelihood, or unweighted least
#' squares (also called minres), respectively, to fit the EFAs. Default is "PAF".
#' @param rotation character vector. Either perform no rotation ("none"; default),
#' any combination of orthogonal rotations ("varimax", "equamax", "quartimax", "geominT",
#' "bentlerT", and "bifactorT"; using "orthogonal" runs all of these), or of
#' oblique rotations ("promax", "oblimin", "quartimin", "simplimax", "bentlerQ",
#' "geominQ", and "bifactorQ"; using "oblique" runs all of these). Rotation types
#' (no rotation, orthogonal rotations, and oblique rotations) cannot be mixed.
#' Default is "promax".
#' @param type character vector. Any combination of "none" (default), "EFAtools",
#' "psych", and "SPSS" can be entered. "none" allows the specification of of various
#' combinations of the arguments controlling both factor extraction methods and
#' the rotations. The others ("EFAtools", "psych", and "SPSS"), control the execution
#' of the respective factor extraction method and rotation to be in line with how
#' it is executed in this package (i.e., the respective default procedure), in the
#' psych package, and in SPSS. A specific psych implementation exists for PAF, ML, varimax,
#' and promax. The SPSS implementation exists for PAF, varimax, and promax. For
#' details, see \code{\link{EFA}}.
#' @param aggregation character. One of "mean" (default), and "median". Controls
#' whether the different results should be aggregated using the (trimmed) mean,
#' or the median.
#' @param trim numeric. If aggregation is set to "mean", this argument controls
#' the trimming of extremes (for details see \code{\link[base:mean]{base::mean}}).
#' By default no trimming is done (i.e., trim = 0).
#' @param salience_threshold numeric. The threshold to use to classify a pattern
#' coefficient or loading as salient (i.e., substantial enough to assign it to
#' a factor). Default is 0.3. Indicator-to-factor correspondences will be inferred
#' based on this threshold. Note that this may not be meaningful if rotation = "none"
#' and n_factors > 1 are used, as no simple structure is present there.
#' @param max_iter numeric. The maximum number of iterations to perform after which
#' the iterative PAF procedure is halted with a warning. Default is 10,000. Note
#' that non-converged procedures are excluded from the aggregation procedure.
#' @param init_comm character vector. Any combination of "smc", "mac", and "unity".
#' Controls the methods to estimate the initial communalities in \code{PAF} if
#' "none" is among the specified types. "smc" will use squared multiple
#' correlations, "mac" will use maximum absolute correlations, "unity" will use
#' 1s (for details see \code{\link{EFA}}). Default is \code{c("smc", "mac", "unity")}.
#' @param criterion numeric vector. The convergence criterion used for PAF if
#' "none" is among the specified types.
#' If the change in communalities from one iteration to the next is smaller than
#' this criterion the solution is accepted and the procedure ends.
#' Default is \code{0.001}.
#' @param criterion_type character vector. Any combination of "max_individual" and
#' "sum". Type of convergence criterion used for PAF if "none" is among the
#' specified types. "max_individual" selects the maximum change in any of the
#' communalities from one iteration to the next and tests it against the
#' specified criterion. "sum" takes the difference of
#' the sum of all communalities in one iteration and the sum of all communalities
#' in the next iteration and tests this against the criterion
#' (for details see \code{\link{EFA}}). Default is \code{c("sum", "max_individual")}.
#' @param abs_eigen logical vector. Any combination of TRUE and FALSE.
#' Which algorithm to use in the PAF iterations if "none" is among the specified
#' types. If FALSE, the loadings are computed from the eigenvalues. This is also
#' used by the \code{\link[psych:fa]{psych::fa}} function. If TRUE the
#' loadings are computed with the absolute eigenvalues as done by SPSS
#' (for details see \code{\link{EFA}}). Default is \code{TRUE}.
#' @param varimax_type character vector. Any combination of "svd" and "kaiser".
#' The type of the varimax rotation performed if "none" is among the specified
#' types and "varimax", "promax", "orthogonal", or "oblique" is among the specified
#' rotations. "svd" uses singular value decomposition, as
#' \link[stats:varimax]{stats::varimax} does, and "kaiser" uses the varimax
#' procedure performed in SPSS. This is the original procedure from Kaiser (1958),
#' but with slight alterations in the varimax criterion (for details, see
#' \code{\link{EFA} and Grieder & Steiner, 2020).
#' Default is \code{c("svd", "kaiser")}.
#' @param normalize logical vector. Any combination of TRUE and FALSE.
#' \code{TRUE} performs a kaiser normalization before the specified rotation(s).
#' Default is \code{TRUE}.
#' @param k_promax numeric vector. The power used for computing the target matrix
#' P in the promax rotation if "none" is among the specified types and "promax"
#' or "oblique" is among the specified rotations. Default is \code{2:4}.
#' @param k_simplimax numeric. The number of 'close to zero loadings' for the
#' simplimax rotation (see \code{\link[GPArotation:GPA]{GPArotation::GPFoblq}})
#' if "simplimax" or "oblique" is among the specified rotations. Default
#' is \code{nrow(L)}, where L is the matrix of unrotated loadings.
#' @param P_type character vector. Any combination of "norm" and "unnorm".
#' This specifies how the target matrix P is computed in promax rotation if
#' "none" is among the specified types and "promax" or "oblique" is among the
#' specified rotations. "unnorm" will use the unnormalized target matrix as
#' originally done in Hendrickson and White (1964). "norm" will use a
#' normalized target matrix (for details see \code{\link{EFA}}).
#' Default is \code{c("norm", "unnorm")}.
#' @param precision numeric vector. The tolerance for stopping in the rotation
#' procedure(s). Default is 10^-5.
#' @param start_method character vector. Any combination of "psych" and "factanal".
#' How to specify the starting values for the optimization procedure for ML.
#' "psych" takes the starting values specified in \link[psych:fa]{psych::fa}.
#' "factanal" takes the starting values specified in the
#' \link[stats:factanal]{stats::factanal} function. Default is
#' \code{c("psych", "factanal")}.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
#' is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. Passed to \code{\link[stats:cor]{stats::cor}}.
#' Default is "pearson".
#' @param show_progress logical. Whether a progress bar should be shown in the
#' console. Default is TRUE.
#'
#' @details
#' TBD
#'
#' @return A list of class EFA_AVERAGE containing
#' \item{orig_R}{Original correlation matrix.}
#' \item{h2}{A list with the average, standard deviation, minimum, and maximum
#' final communality estimates.}
#'
#' HERE
#'
#' \item{orig_eigen}{Eigen values of the original correlation matrix.}
#' \item{init_eigen}{Initial eigenvalues, obtained from the correlation matrix
#'  with the initial communality estimates as diagonal in PAF.}
#' \item{final_eigen}{Eigenvalues obtained from the correlation matrix
#'  with the final communality estimates as diagonal.}
#' \item{iter}{The number of iterations needed for convergence.}
#' \item{convergence}{Integer code for convergence as returned by
#' \code{\link[stats:optim]{stats:optim}} (only for ML and ULS).
#' 0 indicates successful completion.}
#' \item{unrot_loadings}{Loading matrix containing the final unrotated loadings.}
#' \item{vars_accounted}{Matrix of explained variances and sums of squared loadings. Based on the unrotated loadings.}
#' \item{fit_indices}{For ML and ULS: Fit indices derived from the unrotated
#' factor loadings: Chi Square, including significance level, degrees of freedom
#' (df), Comparative Fit Index (CFI), Root Mean Square Error of Approximation
#' (RMSEA), including its 90\% confidence interval, and the common part accounted
#' for (CAF) index as proposed by Lorenzo-Seva, Timmerman, & Kiers (2011).
#' For PAF, only the CAF and dfs are returned.}
#' \item{rot_loadings}{Loading matrix containing the final rotated loadings
#' (pattern matrix).}
#' \item{Phi}{The factor intercorrelations (only for oblique rotations).}
#' \item{Structure}{The structure matrix (only for oblique rotations).}
#' \item{rotmat}{The rotation matrix.}
#' \item{vars_accounted_rot}{Matrix of explained variances and sums of squared
#' loadings. Based on rotated loadings and, for oblique rotations, the factor
#' intercorrelations.}
#' \item{settings}{A list of the settings used.}
#'
#' @source Grieder, S., & Steiner, M.D. (2020). Algorithmic Jingle Jungle:
#' A Comparison of Implementations of Principal Axis Factoring and Promax Rotation
#'  in R and SPSS. Manuscript in Preparation.
#' @source Hendrickson, A. E., & White, P. O. (1964). Promax: A quick method for
#' rotation to oblique simple structure. British Journal of Statistical Psychology,
#' 17 , 65–70. doi: 10.1111/j.2044-8317.1964.tb00244.x
#' @source Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. L. (2011). The
#' Hull Method for Selecting the Number of Common Factors, Multivariate Behavioral
#' Research, 46, 340-364, doi: 10.1080/00273171.2011.564527
#' @source Kaiser, H. F. (1958). The varimax criterion for analytic rotation in
#' factor analysis. Psychometrika, 23, 187–200. doi: 10.1007/BF02289233
#'
#' @export
#'
#' @examples
#' TBD
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

  if (n_factors == 1 && !all(rotation == "none")) {
    message(cli::col_cyan(cli::symbol$info, " 'n_factors' is 1, but rotation != 'none'. Setting rotation to 'none' to avoid many warnings, as 1-factor solutions cannot be rotated.\n"))
    rotation <- "none"
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
          grid_list[["psch_ls"]] <- .type_grid(method = "ULS", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, start_method = NA,
                                               rotation = rotation, k_promax = 4,
                                               normalize = TRUE, P_type = "unnorm",
                                               precision = 1e-5,
                                               varimax_type = "svd",
                                               k_simplimax = k_simplimax)
        }

        if ("SPSS" %in% type) {
          grid_list[["spss_ls"]] <- .type_grid(method = "ULS", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, start_method = NA,
                                               rotation = rotation, k_promax = 4,
                                               normalize = TRUE, P_type = "norm",
                                               precision = 1e-5,
                                               varimax_type = "kaiser",
                                               k_simplimax = k_simplimax)
        }

        if ("none" %in% type) {

            grid_list[["nn_ls"]] <- .type_grid(method = "ULS", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, start_method = NA,
                                               rotation = rotation, k_promax = k_promax,
                                               normalize = normalize, P_type = P_type,
                                               precision = precision,
                                               varimax_type = varimax_type,
                                               k_simplimax = k_simplimax)

        }
      }

  arg_grid <- unique(do.call(rbind, grid_list))

  ### Run all efas

  if (nrow(arg_grid) == 1) {

    warning(crayon::yellow$bold("!"), crayon::yellow(" There was only one combination of arguments, returning normal EFA output.\n"))

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
  } else {
    progressr::handlers("void")
  }

  progressr::with_progress({
    n_efa <- nrow(arg_grid)
    if (n_efa <= 10) {
      stepsize <- 1
    } else {
      stepsize <- round(n_efa / 100 * 10)
    }

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
              order_type = "eigen", start_method = start_methods[i], maxit = 5e4),
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

  if (n_factors > 1) {
    if (isTRUE(show_progress)) {
      .show_agg_progress("\U0001f6b6", "Reordering factors...")
    }

    re_list <- .array_reorder(ext_list$vars_accounted, ext_list$L, ext_list$L_corres,
                              ext_list$phi, ext_list$extract_phi, n_factors)
  } else {
    re_list <- ext_list
  }


  if (isTRUE(show_progress)) {
    .show_agg_progress("\U0001f3c3", "Aggregating data...")
  }
  agg_list <- suppressWarnings(
    .aggregate_values(re_list$vars_accounted, re_list$L, re_list$L_corres, ext_list$h2, re_list$phi,
                      ext_list$extract_phi, aggregation, trim,
                      ext_list$for_grid[, c("chisq", "p_chi", "caf", "cfi",
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

