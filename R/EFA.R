#' Function to run PAF and obtain varimax or promax rotated loadings
#'
#' This is just a wrapper function to call \code{\link{PAF}} and then one of
#' \code{\link{PROMAX}} or \code{\link{VARIMAX}} or no rotation is performed.
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
#' @param rotation character. One of "promax" (default), "varimax", or "none".
#'  Specifies the type of rotation to perform.
#' @param type character. If one of "EFAdiff" (default), "psych", or "SPSS" is
#'  used, and the following arguments with default NULL are left with
#'  NULL, these implementations are executed as reported in Steiner and Grieder
#'  (2019). Individual properties can be adapted using one of the three types and
#'  specifying some of the following arguments. If set to another value than one
#'  of the three specified above, all following arguments must be specified. For
#'  details see the helppages for \code{\link{PAF}}, \code{\link{PROMAX}}, or
#'  \code{\link{VARIMAX}}
#' @param max_iter numeric. The maximum number of iterations (default is 300) to
#'  perform after which the iterative PAF procedure is halted with a warning.
#' @param init_comm character. The method to estimate the initial communalities.
#'  "smc" will use squared multiple correlations. "mac" will use
#'   maximum absolute correlations. "unity" will use 1s.
#' @param criterion numeric. The convergence criterion used in \code{\link{PAF}}.
#'  If the change in communalities from one iteration to the next is smaller than
#'  this criterion the solution is accepted and the procedure ends. Details
#'  depend on criterion_type (see \code{\link{PAF}} documentation).
#' @param criterion_type character. Type of convergence criterion used in
#'  \code{\link{PAF}}. "max_individual" selects the maximum change in any of the
#'  communalities from one iteration to the next and tests it against the
#'  specified criterion. This is also used by SPSS. "sums" takes difference of
#'  the sum of all communalities in one iteration and the sum of all communalities
#'  in the next iteration and tests it against the criterion. This procedure is
#'  used by the \code{\link[psych:fa]{psych::fa}} function.
#' @param abs_eigen logical. Which algorithm to use in the \code{\link{PAF}}
#'  iterations. If FALSE, the loadings are computed from the eigenvalues. This is
#'  also used by the \code{\link[psych:fa]{psych::fa}} function. If TRUE the
#'  loadings are computed with the absolute eigenvalues as done by SPSS. All
#'  textbooks we know use \code{abs_eigen = FALSE}.
#' @param signed_loadings  logical. If \code{TRUE} (default), the sign of
#' factors after PAF with negative sum of loadings is reflected. This is done by
#' both SPSS and \code{\link[psych:fa]{psych::fa}}.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
#'  is given as input. Note that in this case \code{cors} must be set to
#'  \code{FALSE}.
#' @param k numeric. The power used for computing the target matrix P in the
#'  promax rotation.
#' @param kaiser logical. If \code{TRUE} (default), kaiser normalization is
#' performed in the varimax rotation.
#' @param P_type character. Default is \code{NULL}. This specifies how the target
#'  matrix P is computed in a promax rotation. If "HW" it will use the
#'  implementation of Hendrickson and White (1964), which is also used by
#'  the psych and stats packages. If "SPSS" it will use the SPSS version
#'  (see \code{\link{PROMAX}} or Steiner and Grieder, 2019 for details).
#' @param precision numeric. The tolerance for stopping in the varimax procecdure.
#'  This is passed to the "eps" argument of the
#'  \code{\link[stats:varimax]{stats::varimax}} function. Default is \code{NULL}.
#' @param order_type character. How to order the factors and when to reflect
#'  their signs. Default is \code{NULL}. "psych" will use the psych method, "SPSS" the
#'  SPSS method. See \code{\link{PROMAX}} or \code{\link{VARIMAX}} for details.
#'
#' @return A list of type PAF if no rotation is used (see \code{\link{PAF}}
#'  documentation for a description), and of type PROMAX, or
#'  VARIMAX, if the respective rotation is used. In the latter cases the output
#'  list is of the following structure:
#'
#' \item{orig_R}{Original correlation matrix.}
#' \item{h2_init}{Initial communality estimates from PAF.}
#' \item{h2}{Final communality estimates from unrotated loadings.}
#' \item{iter}{The number of iterations needed for convergence in PAF.}
#' \item{orig_eigen}{Eigen values of the original correlation matrix.}
#' \item{init_eigen}{Initial eigenvalues, obtained from the correlation matrix
#'  with the initial communality estimates as diagonal in PAF.}
#' \item{final_eigen}{Eigenvalues of the final iteration in PAF.}
#' \item{unrot_loadings}{Loading matrix containing the unrotated final loadings from PAF.}
#' \item{rot_loadings}{The promax or varimax rotated loadings (the pattern matrix).}
#' \item{rotmat}{The rotation matrix.}
#' \item{Phi}{The factor intercorrelations. Only returned if promax rotation is used.}
#' \item{Structure}{The structure matrix. Only returned if promax rotation is used.}
#' \item{vars_accounted}{Matrix of explained variances and sums of squared loadings. Based on rotated loadings and, if applicable, the factor intercorrelations.}
#' \item{fit_indices}{Fit indices derived from the rotated factor loadings as
#'  returned by \code{\link[psych:factor.stats]{psych::factor.stats}}}
#'
#' @export
EFA <- function(x, n_factors, cors = TRUE, N = NA, rotation = "promax",
                type = "EFAdiff", max_iter = 300, init_comm = NULL,
                criterion = NULL, criterion_type = NULL, abs_eigen = NULL,
                signed_loadings = TRUE, use = "pairwise.complete.obs",
                k = 4, kaiser = TRUE, P_type = NULL,
                precision = NULL, order_type = NULL) {

  # run principal axis factoring
  paf_out <- PAF(x, n_factors = n_factors, cors = cors, N = N, type = type,
                 max_iter = max_iter, init_comm = init_comm, criterion = criterion,
                 criterion_type = criterion_type, abs_eigen = abs_eigen,
                 signed_loadings = signed_loadings, use = use)

  # rotate factor analysis results
  if (rotation == "promax") {

    paf_rot <- PROMAX(paf_out, type = type, kaiser = kaiser, P_type = P_type,
                      precision = precision, order_type = order_type, N = N)

  } else if (rotation == "varimax") {
    paf_rot <- VARIMAX(paf_out, type = type, kaiser = kaiser, precision = precision,
                       order_type = order_type, N = N)
  } else {
    if (rotation != "none") {
      warning("You entered rotation = '", rotation, "', which is no valid input.",
              " Valid inputs are 'promax', 'varimax', or 'none'. No rotation was done.")
    }

    return(paf_out)
  }

  # prepare output
  output <- list(
    orig_R = paf_out$orig_R,
    h2_init = paf_out$h2_init,
    h2 = paf_out$h2,
    iter = paf_out$iter,
    orig_eigen = paf_out$orig_eigen,
    init_eigen = paf_out$init_eigen,
    final_eigen = paf_out$final_eigen,
    unrot_loadings = paf_out$loadings,
    rot_loadings = paf_rot$rot_loadings,
    rotmat = paf_rot$rotmat,
    vars_accounted = paf_rot$vars_accounted,
    fit_indices = paf_rot$fit_indices
  )

  if (rotation == "promax" || rotation != "varimax") {
    output$Phi = paf_rot$fit_indices
    output$Structure = paf_rot$Structure
    class(output) <- "PROMAX"
  } else if (rotation == "varimax") {
    class(output) <- "VARIMAX"
  }

  output

}
