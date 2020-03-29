#' Function to run PAF and obtain varimax or promax rotated loadings
#'
#' This is a wrapper function to call \code{\link{PAF}} and then one of
#' \code{\link{PROMAX}} or \code{\link{VARIMAX}} or no rotation is performed.
#' All arguments with default value \code{NULL} can be left to default if \code{type}
#' is set to one of "EFAtools", "SPSS", or "psych". The respective specifications are
#' then handled by the function.
#'
#' @param x data.frame or matrix. Dataframe or matrix of raw data or matrix with
#'  correlations.
#' @param n_factors numeric. Number of factors to extract.
#' @param cors logical. If \code{TRUE} (default) a correlation matrix is expected in x,
#'  otherwise raw data is expected.
#' @param N numeric. The number of observations. Needs only be specified if a
#'  correlation matrix is used. If input is a correlation matrix and N = NA
#'  (default), not all fit indices can be computed. See
#'  \code{\link[psych:factor.stats]{psych::factor.stats}} for details.
#' @param method -...
#' @param rotation character. One of "none" (default), "promax", or "varimax".
#'  Specifies the type of rotation to perform.
#' @param type character. If one of "EFAtools" (default), "psych", or "SPSS" is
#'  used, and the following arguments with default NULL are left with
#'  NULL, these implementations are executed as reported in Grieder and Steiner
#'  (2019). Individual properties can be adapted using one of the three types and
#'  specifying some of the following arguments. If set to "none" additional
#'  arguments must be specified depending on the \code{method} and \code{rotation}
#'  used.
#' @param max_iter numeric. The maximum number of iterations to perform after which
#' the iterative PAF procedure is halted with a warning. If \code{type} is one of
#' "EFAtools", "SPSS", or "psych", this is automatically specified if \code{max_iter} is
#' left to be \code{NULL}, but can be overridden by entering a number. Default is \code{NULL}.
#' @param init_comm character. The method to estimate the initial communalities.
#'  "smc" will use squared multiple correlations. "mac" will use
#'   maximum absolute correlations. "unity" will use 1s. Default is \code{NULL}.
#' @param criterion numeric. The convergence criterion used in \code{\link{PAF}}.
#'  If the change in communalities from one iteration to the next is smaller than
#'  this criterion the solution is accepted and the procedure ends. Details
#'  depend on criterion_type (see \code{\link{PAF}} documentation). Default is \code{NULL}.
#' @param criterion_type character. Type of convergence criterion used in
#'  \code{\link{PAF}}. "max_individual" selects the maximum change in any of the
#'  communalities from one iteration to the next and tests it against the
#'  specified criterion. This is also used by SPSS. "sums" takes difference of
#'  the sum of all communalities in one iteration and the sum of all communalities
#'  in the next iteration and tests it against the criterion. This procedure is
#'  used by the \code{\link[psych:fa]{psych::fa}} function. Default is \code{NULL}.
#' @param abs_eigen logical. Which algorithm to use in the \code{\link{PAF}}
#'  iterations. If FALSE, the loadings are computed from the eigenvalues. This is
#'  also used by the \code{\link[psych:fa]{psych::fa}} function. If TRUE the
#'  loadings are computed with the absolute eigenvalues as done by SPSS. All
#'  textbooks we know use \code{abs_eigen = FALSE}. Default is \code{NULL}.
#' @param signed_loadings  logical. If \code{TRUE} (default), the sign of
#' factors after PAF with negative sum of loadings is reflected. This is done by
#' both SPSS and \code{\link[psych:fa]{psych::fa}}.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
#'  is given as input. Note that in this case \code{cors} must be set to
#'  \code{FALSE}. Default is "pairwise.complete.obs".
#' @param k numeric. The power used for computing the target matrix P in the
#'  promax rotation. Default is \code{NULL}.
#' @param kaiser logical. If \code{TRUE}, kaiser normalization is
#' performed before the varimax rotation. Default is \code{NULL}.
#' @param P_type character. This specifies how the target
#'  matrix P is computed in a promax rotation. If "unnorm" it will use the
#'  unnormalized target matrix as originally done in Hendrickson and White (1964).
#'  This is also used done the psych and stats packages. If "norm" it will use the
#'  normalized target matrix as used in SPSS (see \code{\link{PROMAX}} or
#'  Grieder and Steiner, 2019 for details). Default is \code{NULL}.
#' @param precision numeric. The tolerance for stopping in the varimax procecdure.
#'  This is passed to the "eps" argument of the
#'  \code{\link[stats:varimax]{stats::varimax}} function. Default is \code{NULL}.
#' @param order_type character. How to order the factors. "eigen" will reorder
#'  the factors according to the largest to lowest eigenvalues. "ss_factors" will
#'  reorder the factors according to descending sum of squared factor loadings
#'  per factor. See \code{\link{PROMAX}} or \code{\link{VARIMAX}} for details.
#'  @param start_method character. How to specify the starting values for the
#'  optimization prodedure. Default is "factanal" which takes the starting values
#'  specified in the \link{stats}{factanal} function. "psych" takes the starting
#'  values specified in \link{psych}{fa}. Solutions are very similar.
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
                rotation = c("none", "varimax", "promax"),
                type = c("EFAtools", "psych", "SPSS", "none"), max_iter = NULL,
                init_comm = NULL, criterion = NULL, criterion_type = NULL,
                abs_eigen = NULL, signed_loadings = TRUE,
                use = c("all.obs", "complete.obs", "pairwise.complete.obs",
                "everything", "na.or.complete"), k = NULL,
                kaiser = TRUE, P_type = NULL, precision = NULL,
                order_type = NULL, start_method = c("factanal", "psych")) {

  # # for testing
  # x <- IDS2_R
  # N <- 1991
  # n_factors <- 5
  # method <- "PAF"
  # rotation <- "promax"
  # cors <- TRUE
  # type = "EFAtools"
  # signed_loadings = TRUE
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
                 signed_loadings = signed_loadings, use = use)

  } else if (method == "ML") {

    fit_out <- ML(x, n_factors = n_factors, cors = cors, N = N,
                  signed_loadings = signed_loadings, start_method = start_method,
                  use = use)

  } else if (method == "ULS") {

    fit_out <- ULS(x, n_factors = n_factors, cors = cors, N = N,
                  signed_loadings = signed_loadings, use = use)
  }



  # rotate factor analysis results
  if (rotation == "promax") {

    rot_out <- PROMAX(fit_out, type = type, kaiser = kaiser, P_type = P_type,
                      precision = precision, order_type = order_type, k = k)

  } else if (rotation == "varimax") {

    rot_out <- VARIMAX(fit_out, type = type, kaiser = kaiser, precision = precision,
                       order_type = order_type)

  } else {

    output <- fit_out

  }

  if (rotation != "none"){

    settings <- c(fit_out$settings, rot_out$settings)
    output <- c(within(fit_out, rm(settings)), within(rot_out, rm(settings)),
                      settings = list(settings))

  }

  # Add settings used to output
  settings_EFA <- list(
    method = method,
    rotation = rotation,
    type = type
  )

  settings <- c(settings_EFA, output$settings)

  output <- c(within(output, rm(settings)),
                 settings = list(settings))

  class(output) <- "EFA"

  output

}
