#' Function to run PAF and obtain varimax or promax rotated loadings
#'
#' This is a function does an EFA with either \code{\link{PAF}}, \code{\link{ML}},
#' or \code{\link{ULS}} with or without subsequent rotation.
#' All arguments with default value \code{NULL} can be left to default if \code{type}
#' is set to one of "EFAtools", "SPSS", or "psych". The respective specifications are
#' then handled according to the specified type (see details). For all rotations except varimax and promax, the \code{\link{GPArotation} package is needed.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#' correlations.
#' @param n_factors numeric. Number of factors to extract.
#' @param cors logical. If \code{TRUE} (default) a correlation matrix is expected
#' in x, otherwise raw data is expected.
#' @param N numeric. The number of observations. Needs only be specified if a
#' correlation matrix is used. If input is a correlation matrix and \code{N} = NA
#' (default), not all fit indices can be computed. See
#'  \code{\link[psych:factor.stats]{psych::factor.stats}} for details.
#' @param method character. One of "PAF", "ML", or "ULS" to use principal axis
#' factoring, maximum likelihood, or unweighted least squares, respectively,
#' to fit the EFA.
#' @param rotation character. Either perform no rotation ("none"; default),
#' an orthogonal rotation ("varimax", "equamax", "quartimax", "geominT",
#' "bentlerT", or "bifactorT"), or an oblique rotation ("promax", "oblimin",
#' "quartimin", "simplimax", "bentlerQ", "geominQ", or "bifactorQ").
#' @param type character. If one of "EFAtools" (default), "psych", or "SPSS" is
#'  used, and the following arguments with default NULL are left with
#'  NULL, these implementations are executed according to the respective program
#'  ("psych" and "SPSS") or according to the best solution found in Grieder &
#'  Steiner (2020; "EFAtools). Individual properties can be adapted using one of
#'  the three types and specifying some of the following arguments. If set to
#'  "none" additional arguments must be specified depending on the \code{method}
#'  and \code{rotation} used (see details).
#' @param max_iter numeric. The maximum number of iterations to perform after which
#' the iterative PAF procedure is halted with a warning. If \code{type} is one of
#' "EFAtools", "SPSS", or "psych", this is automatically specified if \code{max_iter} is
#' left to be \code{NULL}, but can be overridden by entering a number. Default is
#' \code{NULL}.
#' @param init_comm character. The method to estimate the initial communalities
#' in \code{\link{PAF}}. "smc" will use squared multiple correlations, "mac" will use
#' maximum absolute correlations, "unity" will use 1s (see details).
#' Default is \code{NULL}.
#' @param criterion numeric. The convergence criterion used in \code{\link{PAF}}.
#' If the change in communalities from one iteration to the next is smaller than
#' this criterion the solution is accepted and the procedure ends. Details
#' depend on criterion_type (see \code{\link{PAF}} documentation).
#' Default is \code{NULL}.
#' @param criterion_type character. Type of convergence criterion used in
#' \code{\link{PAF}}. "max_individual" selects the maximum change in any of the
#' communalities from one iteration to the next and tests it against the
#' specified criterion. This is also used by SPSS. "sums" takes the difference of
#' the sum of all communalities in one iteration and the sum of all communalities
#' in the next iteration and tests this against the criterion. This procedure is
#' used by the \code{\link[psych:fa]{psych::fa}} function. Default is \code{NULL}.
#' @param abs_eigen logical. Which algorithm to use in the \code{\link{PAF}}
#' iterations. If FALSE, the loadings are computed from the eigenvalues. This is
#' also used by the \code{\link[psych:fa]{psych::fa}} function. If TRUE the
#' loadings are computed with the absolute eigenvalues as done by SPSS.
#' Default is \code{NULL}.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
#' is given as input. Note that in this case \code{cors} must be set to
#' \code{FALSE}. Default is "pairwise.complete.obs".
#' @param k numeric. Either the power used for computing the target matrix P in
#' the promax rotation or the number of 'close to zero loadings' for the simplimax
#' rotation (see \code{\link[GPArotation:GPFobl]{GPArotation:GPFobl}. If left to
#' \code{NULL} (default), 3 is used for promax and \code{nrow(L)}, where
#' L is the matrix of unrotated loadings, is used for simplimax.
#' @param kaiser logical. If \code{TRUE}, kaiser normalization is
#' performed before the specified rotation. Default is \code{TRUE}.
#' @param P_type character. This specifies how the target
#' matrix P is computed in promax rotation. If "unnorm" it will use the
#' unnormalized target matrix as originally done in Hendrickson and White (1964).
#' This is also used in the psych and stats packages. If "norm" it will use the
#' normalized target matrix as used in SPSS (see \code{\link{PROMAX}}
#' for details). Default is \code{NULL}.
#' @param precision numeric. The tolerance for stopping in the rotation
#' procecdure. This is passed to the "eps" argument of the
#' \code{\link[stats:varimax]{stats::varimax}} and the GPArotation functions.
#' If left \code{NULL} (default), the precision is set according to the specified
#' \code{type} for varimax and promax rotation or is set to 10^-5 for all other
#' rotations.
#' @param order_type character. How to order the factors. "eigen" will reorder
#' the factors according to the largest to lowest eigenvalues of the matrix of
#' rotated loadings. "ss_factors" will reorder the factors according to descending
#' sum of squared factor loadings per factor. Default is \code{NULL}.
#' @param start_method character. How to specify the starting values for the
#' optimization prodedure for ML. Default is "factanal" which takes the starting
#' values specified in the \link{stats}{factanal} function. "psych" takes the
#' starting values specified in \link{psych}{fa}. Solutions are very similar.
#' @param ... Additional arguments passed to rotation functions from GPArotation
#' package (e.g., \code{maxit} for maximum number of iterations).
#'
#' HIER STEHENGEBLIEBEN, ARGUMENTE OBEN IN ANDERE FUNKTIONEN REINKOPIEREN,
#' DETAILS FERTIG SCHREIBEN
#'
#' @details There are two main ways to use this function. The easiest way is to
#' use it with a specified \code{type} (see above), which sets most of the other
#' arguments accordingly. Another way is to use it more flexibly by explicitly
#' specifying all arguments used and set \code{type} to "none". A mix of the two
#' can also be done by specifing a \code{type} as well as additional arguments.
#' However, this will throw warnings to avoid unintentional deviations from the
#' implementations according to the specified \code{type}.
#'
#' type no effect on uls and ml
#'
#'what is done with what type
#'
#' init_Comm by type
#'
#' @return A list of class PAF if no rotation is used (see \code{\link{PAF}}
#'  documentation for a description), and of class PROMAX, or
#'  VARIMAX, if the respective rotation is used. In the latter cases the output
#'  list is of the following structure:
#'
#' \item{orig_R}{Original correlation matrix.}
#' \item{h2_init}{Initial communality estimates from PAF.}
#' \item{h2}{Final communality estimates from the unrotated solution}
#' \item{iter}{The number of iterations needed for convergence in PAF.}
#' \item{orig_eigen}{Eigen values of the original correlation matrix.}
#' \item{init_eigen}{Initial eigenvalues, obtained from the correlation matrix
#'  with the initial communality estimates as diagonal in PAF.}
#' \item{final_eigen}{Eigenvalues of the final iteration in PAF.}
#' \item{unrot_loadings}{Loading matrix containing the unrotated final loadings.}
#' \item{rot_loadings}{The pattern matrix of the rotated solution).}
#' \item{rotmat}{The rotation matrix.}
#' \item{Phi}{The factor intercorrelations. Only returned if promax rotation is used.}
#' \item{Structure}{The structure matrix. Only returned if promax rotation is used.}
#' \item{vars_accounted}{Matrix of explained variances and sums of squared loadings. Based on the unrotated loadings.}
#' \item{vars_accounted_rot}{Matrix of explained variances and sums of squared loadings. Based on rotated loadings and, if applicable, the factor intercorrelations.}
#' \item{fit_indices}{Fit indices derived from the unrotated factor loadings as
#'  returned by \code{\link[psych:factor.stats]{psych::factor.stats}}}
#' \item{settings}{list. The settings (arguments) used in the EFA.}
#'
#' @source Hendrickson, A. E., & White, P. O. (1964). Promax: A quick method for rotation to oblique simple structure. British Journal of Statistical Psychology, 17 , 65–70. doi: 10.1111/j.2044-8317.1964.tb00244.x
#' @source Grieder, S., & Steiner, M.D.(2019). Algorithmic Jingle Jungle: Comparison of Implementations of an EFA Procedure in R psych Versus SPSS, MacOrtho, and Omega. Submitted Manuscript.
#'
#' @export
#'
#' @examples
#' ADD EXAMPLE FOR type = "none"
#'
#' # A type EFAtools (as presented in Steiner and Grieder, 2019) EFA
#' EFA_EFAtools_5 <- EFA(IDS2_R, n_factors = 5, type = "EFAtools", method = "PAF",
#'                       rotation = "none")
#'
#' # A type SPSS EFA to mimick the SPSS implementation
#' EFA_SPSS_5 <- EFA(IDS2_R, n_factors = 5, type = "SPSS", method = "PAF",
#'                   rotation = "none")
#'
#' # A type psych EFA to mimick the psych::fa() implementation
#' EFA_psych_5 <- EFA(IDS2_R, n_factors = 5, type = "psych", method = "PAF",
#'                    rotation = "none")
EFA <- function(x, n_factors, cors = TRUE, N = NA, method = c("PAF", "ML", "ULS"),
                rotation = c("none", "varimax", "equamax", "quartimax", "geominT",
                             "bentlerT", "bifactorT", "promax", "oblimin",
                             "quartimin", "simplimax", "bentlerQ", "geominQ",
                             "bifactorQ"),
                type = c("EFAtools", "psych", "SPSS", "none"), max_iter = NULL,
                init_comm = NULL, criterion = NULL, criterion_type = NULL,
                abs_eigen = NULL, use = c("pairwise.complete.obs", "all.obs",
                                          "complete.obs", "everything",
                                          "na.or.complete"),
                k = NULL, kaiser = TRUE, P_type = NULL, precision = NULL,
                order_type = NULL, start_method = c("factanal", "psych"),
                ...) {

  # # for testing
  # x <- IDS2_R
  # N <- 1991
  # n_factors <- 5
  # method <- "PAF"
  # rotation <- "promax"
  # cors <- TRUE
  # type = "EFAtools"
  # use = "pairwise.complete.obs"
  # kaiser = TRUE
  # max_iter = NULL
  # init_comm = NULL
  # criterion = NULL
  # criterion_type = NULL
  # abs_eigen = NULL
  # k = NULL
  # precision = NULL
  # P_type = NULL
  # order_type = NULL

  method <- match.arg(method)
  rotation <- match.arg(rotation)
  use <- match.arg(use)
  type <- match.arg(type)
  start_method <- match.arg(start_method)

  # run factor analysis with respective fit method

  if (method == "PAF") {

  fit_out <- PAF(x, n_factors = n_factors, cors = cors, N = N, type = type,
                 max_iter = max_iter, init_comm = init_comm, criterion = criterion,
                 criterion_type = criterion_type, abs_eigen = abs_eigen,
                 use = use)

  } else if (method == "ML") {

    fit_out <- ML(x, n_factors = n_factors, cors = cors, N = N,
                  start_method = start_method, use = use)

  } else if (method == "ULS") {

    fit_out <- ULS(x, n_factors = n_factors, cors = cors, N = N, use = use)
  }

  # rotate factor analysis results
  if (rotation == "promax") {

    rot_out <- PROMAX(fit_out, type = type, kaiser = kaiser, P_type = P_type,
                      precision = precision, order_type = order_type, k = k)

  } else if (rotation == "varimax") {

    rot_out <- VARIMAX(fit_out, type = type, kaiser = kaiser,
                       precision = precision, order_type = order_type)

  } else if (rotation == "quartimax" || rotation == "equamax" ||
             rotation == "bentlerT" || rotation == "geominT" ||
             rotation == "bifactorT") {

    rot_out <- ROTATE_ORTH(fit_out, type = type, rotation = rotation,
                           kaiser = kaiser, precision = precision,
                           order_type = order_type, ...)

  } else if (rotation == "oblimin" || rotation == "quartimin" ||
             rotation == "simplimax" || rotation == "bentlerQ" ||
             rotation == "geominQ" || rotation == "bifactorQ") {

    rot_out <- ROTATE_OBLQ(fit_out, type = type, rotation = rotation,
                           kaiser = kaiser, precision = precision,
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
    cors = cors,
    N = N,
    use = use
  )

  if(method == "ULS" & rotation == "none"){

    output <- c(output, settings = list(settings_EFA))

  } else {

    settings <- c(settings_EFA, output$settings)

    output <- c(within(output, rm(settings)),
                settings = list(settings))

  }

  class(output) <- "EFA"

  output

}
