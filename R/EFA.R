#' Exploratory factor analysis (EFA)
#'
#' This function does an EFA with either \code{PAF}, \code{ML},
#' or \code{ULS} with or without subsequent rotation.
#' All arguments with default value \code{NA} can be left to default if \code{type}
#' is set to one of "EFAtools", "SPSS", or "psych". The respective specifications are
#' then handled according to the specified type (see details). For all rotations
#' except varimax and promax, the \code{GPArotation} package is needed.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations. If raw data is entered, the correlation matrix is found from the
#' data.
#' @param n_factors numeric. Number of factors to extract.
#' @param N numeric. The number of observations. Needs only be specified if a
#' correlation matrix is used. If input is a correlation matrix and \code{N} = NA
#' (default), not all fit indices can be computed.
#' @param method character. One of "PAF", "ML", or "ULS" to use principal axis
#' factoring, maximum likelihood, or unweighted least squares (also called minres),
#' respectively, to fit the EFA.
#' @param rotation character. Either perform no rotation ("none"; default),
#' an orthogonal rotation ("varimax", "equamax", "quartimax", "geominT",
#' "bentlerT", or "bifactorT"), or an oblique rotation ("promax", "oblimin",
#' "quartimin", "simplimax", "bentlerQ", "geominQ", or "bifactorQ").
#' @param type character. If one of "EFAtools" (default), "psych", or "SPSS" is
#'  used, and the following arguments with default NA are left with
#'  NA, these implementations are executed according to the respective program
#'  ("psych" and "SPSS") or according to the best solution found in Grieder &
#'  Steiner (2020; "EFAtools"). Individual properties can be adapted using one of
#'  the three types and specifying some of the following arguments. If set to
#'  "none" additional arguments must be specified depending on the \code{method}
#'  and \code{rotation} used (see details).
#' @param max_iter numeric. The maximum number of iterations to perform after which
#' the iterative PAF procedure is halted with a warning. If \code{type} is one of
#' "EFAtools", "SPSS", or "psych", this is automatically specified if \code{max_iter} is
#' left to be \code{NA}, but can be overridden by entering a number. Default is
#' \code{NA}.
#' @param init_comm character. The method to estimate the initial communalities
#' in \code{PAF}. "smc" will use squared multiple correlations, "mac" will use
#' maximum absolute correlations, "unity" will use 1s (see details).
#' Default is \code{NA}.
#' @param criterion numeric. The convergence criterion used for PAF.
#' If the change in communalities from one iteration to the next is smaller than
#' this criterion the solution is accepted and the procedure ends.
#' Default is \code{NA}.
#' @param criterion_type character. Type of convergence criterion used for
#' PAF. "max_individual" selects the maximum change in any of the
#' communalities from one iteration to the next and tests it against the
#' specified criterion. This is also used by SPSS. "sum" takes the difference of
#' the sum of all communalities in one iteration and the sum of all communalities
#' in the next iteration and tests this against the criterion. This procedure is
#' used by the \code{\link[psych:fa]{psych::fa}} function. Default is \code{NA}.
#' @param abs_eigen logical. Which algorithm to use in the PAF
#' iterations. If FALSE, the loadings are computed from the eigenvalues. This is
#' also used by the \code{\link[psych:fa]{psych::fa}} function. If TRUE the
#' loadings are computed with the absolute eigenvalues as done by SPSS.
#' Default is \code{NA}.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
#' is given as input. Default is "pairwise.complete.obs".
#' @param cor_method character. Passed to \code{\link[stats:cor]{stats::cor}}.
#' Default is "pearson".
#' @param k numeric. Either the power used for computing the target matrix P in
#' the promax rotation or the number of 'close to zero loadings' for the simplimax
#' rotation (see \code{\link[GPArotation:GPA]{GPArotation::GPFoblq}}). If left to
#' \code{NA} (default), the value for promax depends on the specified type.
#' For simplimax, \code{nrow(L)}, where L is the matrix of unrotated loadings,
#' is used by default.
#' @param normalize logical. If \code{TRUE}, a kaiser normalization is
#' performed before the specified rotation. Default is \code{TRUE}.
#' @param P_type character. This specifies how the target
#' matrix P is computed in promax rotation. If "unnorm" it will use the
#' unnormalized target matrix as originally done in Hendrickson and White (1964).
#' This is also used in the psych and stats packages. If "norm" it will use the
#' normalized target matrix as used in SPSS. Default is \code{NA}.
#' @param precision numeric. The tolerance for stopping in the rotation
#' procedure. Default is 10^-5 for all rotation methods.
#' @param varimax_type character. The type of the varimax rotation performed.
#' If "svd", singular value decomposition is used, as \link[stats:varimax]{stats::varimax} does. If "kaiser", the varimax procedure performed in SPSS is used.
#' This is the original procedure from Kaiser (1958), but with slight alterations
#' in the varimax criterion (see details, and Grieder & Steiner, 2020). Default is \code{NA}.
#' @param order_type character. How to order the factors. "eigen" will reorder
#' the factors according to the largest to lowest eigenvalues of the matrix of
#' rotated loadings. "ss_factors" will reorder the factors according to descending
#' sum of squared factor loadings per factor. Default is \code{NA}.
#' @param start_method character. How to specify the starting values for the
#' optimization procedure for ML. Default is "psych" which takes the
#' starting values specified in \link[psych:fa]{psych::fa}. "factanal" takes the
#' starting values specified in the \link[stats:factanal]{stats::factanal} function.
#' Solutions are very similar.
#' @param ... Additional arguments passed to rotation functions from the \code{GPArotation} package (e.g., \code{maxit} for maximum number of iterations).
#'
#' @details There are two main ways to use this function. The easiest way is to
#' use it with a specified \code{type} (see above), which sets most of the other
#' arguments accordingly. Another way is to use it more flexibly by explicitly
#' specifying all arguments used and set \code{type} to "none" (see examples).
#' A mix of the two can also be done by specifying a \code{type} as well as
#' additional arguments. However, this will throw warnings to avoid unintentional
#' deviations from the implementations according to the specified \code{type}.
#'
#' The \code{type} argument is evaluated for PAF and for all rotations (mainly
#' important for the varimax and promax rotations). The type-specific settings
#' for these functions are detailed below.
#'
#' For PAF, the values of \code{init_comm}, \code{criterion}, \code{criterion_type},
#' and \code{abs_eigen} depend on the \code{type} argument.
#'
#' \code{type = "EFAtools"} will use the following argument specification:
#' \code{init_comm = "smc", criterion = .001, criterion_type = "sum",
#' abs_eigen = TRUE}.
#'
#' \code{type = "psych"} will use the following argument specification:
#' \code{init_comm = "smc", criterion = .001, criterion_type = "sum",
#' abs_eigen = FALSE}.
#'
#' \code{type = "SPSS"} will use the following argument specification:
#' \code{init_comm = "smc", criterion = .001, criterion_type = "max_individual",
#' abs_eigen = TRUE}.
#'
#' If SMCs fail, SPSS takes "mac". However, as SPSS takes absolute eigenvalues,
#' this is hardly ever the case. Psych, on the other hand, takes "unity" if SMCs
#' fail. The EFAtools type setting combination was the best in terms of accuracy
#' and number of Heywood cases compared to all the
#' other setting combinations tested in simulation studies in Grieder & Steiner
#' (2020), which is why this type is used as a default here.
#'
#' For varimax, the values of \code{varimax_type} and \code{order_type} depend on
#' the \code{type} argument.
#'
#' \code{type = "EFAtools"} will use the following argument specification:
#' \code{varimax_type = "svd", order_type = "eigen"}.
#'
#' \code{type = "psych"} will use the following argument specification:
#' \code{varimax_type = "svd", order_type = "eigen"}.
#'
#' \code{type = "SPSS"} will use the following argument specification:
#' \code{varimax_type = "kaiser", order_type = "ss_factors"}.
#'
#' For promax, the values of \code{P_type},
#' \code{order_type}, and \code{k} depend on the \code{type} argument.
#'
#' \code{type = "EFAtools"} will use the following argument specification:
#' \code{P_type = "norm", order_type = "eigen", k = 4}.
#'
#' \code{type = "psych"} will use the following argument specification:
#' \code{P_type = "unnorm", order_type = "eigen", k = 4}.
#'
#' \code{type = "SPSS"} will use the following argument specification:
#' \code{P_type = "norm", order_type = "ss_factors", k = 4}.
#'
#' The \code{P_type} argument can take two values, "unnorm" and "norm". It controls
#' which formula is used to compute the target matrix P in the promax rotation.
#' "unnorm" uses the formula from Hendrickson and White (1964), specifically:
#' \code{P = abs(A^(k + 1)) / A},
#' where A is the unnormalized matrix containing varimax rotated loadings.
#' "SPSS" uses the normalized varimax rotated loadings. Specifically it used the
#' following formula, which can be found in the SPSS 23 Algorithms manual:
#' \code{P = abs(A / sqrt(rowSums(A^2))) ^(k + 1) * (sqrt(rowSums(A^2)) / A)}.
#' As for PAF, the EFAtools type setting combination for promax was the best
#' compared to the other setting combinations tested in simulation studies in
#' Grieder & Steiner (2020).
#'
#' The \code{varimax_type} argument can take two values, "svd", and "kaiser". "svd" uses
#' singular value decomposition, by calling \link[stats:varimax]{stats::varimax}. "kaiser"
#' performs the varimax procedure as described in the SPSS 23 Algorithms manual and as described
#' by Kaiser (1958). However, there is a slight alteration in computing the varimax criterion, which
#' we found to better align with the results obtain from SPSS. Specifically, the original varimax
#' criterion as described in the SPSS 23 Algorithms manual is
#' \code{sum(n*colSums(lambda ^ 4) - colSums(lambda ^ 2) ^ 2) / n ^ 2}, where n is the
#' number of indicators, and lambda is the rotated loadings matrix. However, we found the following
#' to produce results more similar to those of SPSS:
#' \code{sum(n*colSums(abs(lambda)) - colSums(lambda ^ 4) ^ 2) / n^2}.
#'
#' For all other rotations except varimax and promax, the \code{type} argument
#' only controls the \code{order_type} argument with the same values as stated
#' above for the varimax and promax rotations. For these other rotations, the
#' \code{GPArotation} package is needed. Additional arguments can also be
#' specified and will be passed to the respective \code{GPArotation} function
#' (e.g., maxit to change the maximum number of iterations for the rotation procedure).
#'
#' The \code{type} argument has no effect on ULS and ML. For ULS, no additional
#' arguments are needed. For ML, an additional argument
#' \code{start_method} is needed to determine the starting values for the
#' optimization procedure. Default for this argument is "factanal" which takes
#' the starting values specified in the \link[stats:factanal]{stats::factanal} function.
#'
#'
#' @return A list of class EFA containing (a subset of) the following:
#'
#' \item{orig_R}{Original correlation matrix.}
#' \item{h2_init}{Initial communality estimates from PAF.}
#' \item{h2}{Final communality estimates from the unrotated solution.}
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
#' (RMSEA), including its 90\% confidence interval, Akaike Information Criterion
#' (AIC), Bayesian Information Criterion (BIC), and the common part accounted
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
#' # A type EFAtools (as presented in Steiner and Grieder, 2020) EFA
#' EFAtools_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                     type = "EFAtools", method = "PAF", rotation = "none")
#'
#' # A type SPSS EFA to mimick the SPSS implementation (this will throw a warning,
#' # see below)
#' SPSS_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                 type = "SPSS", method = "PAF", rotation = "none")
#'
#' # A type psych EFA to mimick the psych::fa() implementation
#' psych_PAF <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                  type = "psych", method = "PAF", rotation = "none")
#'
#' # Use ML instead of PAF with type EFAtools
#' EFAtools_ML <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                    type = "EFAtools", method = "ML", rotation = "none")
#'
#' # Use oblimin rotation instead of no rotation with type EFAtools
#' EFAtools_oblim <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                       type = "EFAtools", method = "PAF", rotation = "oblimin")
#'
#' # Do a PAF without rotation without specifying a type, so the arguments
#' # can be flexibly specified (this is only recommended if you know what your
#' # doing)
#' PAF_none <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                 type = "none", method = "PAF", rotation = "none",
#'                 max_iter = 500, init_comm = "mac", criterion = 1e-4,
#'                 criterion_type = "sum", abs_eigen = FALSE)
#'
#' # Add a promax rotation
#' PAF_pro <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
#'                type = "none", method = "PAF", rotation = "promax",
#'                max_iter = 500, init_comm = "mac", criterion = 1e-4,
#'                criterion_type = "sum", abs_eigen = FALSE, k = 3,
#'                P_type = "unnorm", precision= 1e-5, order_type = "eigen",
#'                varimax_type = "svd")
#'
EFA <- function(x, n_factors, N = NA, method = c("PAF", "ML", "ULS"),
                rotation = c("none", "varimax", "equamax", "quartimax", "geominT",
                             "bentlerT", "bifactorT", "promax", "oblimin",
                             "quartimin", "simplimax", "bentlerQ", "geominQ",
                             "bifactorQ"),
                type = c("EFAtools", "psych", "SPSS", "none"), max_iter = NA,
                init_comm = NA, criterion = NA, criterion_type = NA,
                abs_eigen = NA, use = c("pairwise.complete.obs", "all.obs",
                                          "complete.obs", "everything",
                                          "na.or.complete"),
                varimax_type = NA,
                k = NA, normalize = TRUE, P_type = NA, precision = 1e-5,
                order_type = NA, start_method = "psych",
                cor_method = c("pearson", "spearman", "kendall"),
                ...) {

  # Perform argument checks
  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n"))

  }

  method <- match.arg(method)
  rotation <- match.arg(rotation)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  type <- match.arg(type)
  start_method <- checkmate::matchArg(start_method, c("psych", "factanal", NA))

  if (is.na(start_method) && method == "ML") {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" method is 'ML' but no start method is defined. Please set start_method to 'psych' or 'factanal'.\n"))
  }

  checkmate::assert_count(n_factors)
  checkmate::assert_count(N, na.ok = TRUE)
  checkmate::assert_count(max_iter, na.ok = TRUE)
  checkmate::assert_choice(init_comm, c("smc", "mac", "unity", NA))
  checkmate::assert_number(criterion, lower = 0, upper = 1, na.ok = TRUE)
  checkmate::assert_choice(criterion_type, c("max_individual", "sums", "sum", NA))
  checkmate::assert_flag(abs_eigen, na.ok = TRUE)
  checkmate::assert_number(k, na.ok = TRUE)
  checkmate::assert_choice(varimax_type, c("svd", "kaiser", NA))
  checkmate::assert_flag(normalize, na.ok = TRUE)
  checkmate::assert_choice(P_type, c("unnorm", "norm", NA))
  checkmate::assert_number(precision, lower = 0, upper = 1)
  checkmate::assert_choice(order_type, c("eigen", "ss_factors", NA))

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

  # run factor analysis with respective fit method

  if (method == "PAF") {

  fit_out <- .PAF(R, n_factors = n_factors, N = N, type = type,
                 max_iter = max_iter, init_comm = init_comm,
                 criterion = criterion, criterion_type = criterion_type,
                 abs_eigen = abs_eigen)

  } else if (method == "ML") {

    if (is.na(N)) {

      warning(crayon::yellow$bold("!"), crayon::yellow(" Argument 'N' was NA, not all fit indices could be computed. To get all fit indices, either provide N or raw data.\n"))

    }

    fit_out <- .ML(R, n_factors = n_factors, N = N, start_method = start_method)

  } else if (method == "ULS") {

    if (is.na(N)) {

      warning(crayon::yellow$bold("!"), crayon::yellow(" Argument 'N' was NA, not all fit indices could be computed. To get all fit indices, either provide N or raw data.\n"))

    }

    fit_out <- .ULS(R, n_factors = n_factors, N = N)
  }

  # rotate factor analysis results
  if (rotation == "promax") {

    rot_out <- .PROMAX(fit_out, type = type, normalize = normalize, P_type = P_type,
                      precision = precision, order_type = order_type,
                      varimax_type = varimax_type, k = k)

  } else if (rotation == "varimax") {

    rot_out <- .VARIMAX(fit_out, type = type, normalize = normalize,
                       precision = precision, varimax_type = varimax_type,
                       order_type = order_type)

  } else if (rotation == "quartimax" || rotation == "equamax" ||
             rotation == "bentlerT" || rotation == "geominT" ||
             rotation == "bifactorT") {

    rot_out <- .ROTATE_ORTH(fit_out, type = type, rotation = rotation,
                           normalize = normalize, precision = precision,
                           order_type = order_type, ...)

  } else if (rotation == "oblimin" || rotation == "quartimin" ||
             rotation == "simplimax" || rotation == "bentlerQ" ||
             rotation == "geominQ" || rotation == "bifactorQ") {

    rot_out <- .ROTATE_OBLQ(fit_out, type = type, rotation = rotation,
                           normalize = normalize, precision = precision,
                           order_type = order_type, k = k, ...)

  } else {

    output <- fit_out

  }

  if (rotation != "none"){

    if(method == "ULS"){

      settings <- rot_out$settings
      output <- c(fit_out, within(rot_out, rm(settings)),
                  settings = list(settings))

    } else {

      settings <- c(fit_out$settings, rot_out$settings)
      output <- c(within(fit_out, rm(settings)), within(rot_out, rm(settings)),
                  settings = list(settings))

    }



  }

  # Add settings used to output
  settings_EFA <- list(
    method = method,
    rotation = rotation,
    type = type,
    n_factors = n_factors,
    N = N,
    use = use,
    cor_method = cor_method
  )

  if(method == "ULS" & rotation == "none"){

    output <- c(output, settings = list(settings_EFA))

  } else {

    settings <- c(settings_EFA, output$settings)

    output <- c(within(output, rm(settings)),
                settings = list(settings))

  }

  class(output) <- "EFA"

  return(output)

}
