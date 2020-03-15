#' Principal Axis Factoring
#'
#' This function implements the principal axis factoring procedure. It can
#' reproduce the results from \code{\link[psych:fa]{psych::fa}} and the SPSS
#' FACTOR algorithm. To reproduce psych or SPSS PAF, only the type argument has
#' to be specified in addition to the data and number of factors. The other
#' arguments can be used to control the procedure more flexibly by overriding
#' the default settings. If type = "none" is specified, all arguments with default
#' \code{NULL} have to be specified.
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
#' @param max_iter numeric. The maximum number of iterations to
#'  perform after which the iterative PAF procedure is halted with a warning.
#'  Default is \code{NULL}.
#' @param type character. If one of "EFAtools" (default), "psych", or "SPSS" is
#'  used, and the following arguments (except \code{signed_loadings}) are left with
#'  NULL, these implementations
#'  are executed as reported in Gieder and Steiner (2019; see details).
#'  Individual properties can be adapted using one of the three types and
#'  specifying some of the following
#'  arguments. If set to another value than one of the three specified above, all
#'  arguments with default \code{NULL} must be specified.
#' @param init_comm character. The method to estimate the initial communalities.
#'  "smc" will use squared multiple correlations. "mac" will use
#'   maximum absolute correlations. "unity" will use 1s. Default is \code{NULL}.
#' @param criterion numeric. The convergence criterion.
#'  If the change in communalities from one iteration to the next is smaller than
#'  this criterion the solution is accepted and the procedure ends. Details
#'  depend on criterion_type. Default is \code{NULL}.
#' @param criterion_type character. "max_individual" selects the
#'  maximum change in any of the communalities from one iteration to the next
#'  and tests it against the specified criterion. This is also used by SPSS.
#'  "sums" takes difference of the sum of all communalities in one iteration and
#'  the sum of all communalities in the next iteration and tests it against the
#'  criterion. This procedure is used by the \code{\link[psych:fa]{psych::fa}} function.
#'  Default is \code{NULL}.
#' @param abs_eigen logical. Which algorithm to use in the PAF iterations. If
#'  FALSE, the loadings are computed from the eigenvalues. This is
#'  also used by the \code{\link[psych:fa]{psych::fa}} function. If TRUE the
#'  loadings are computed
#'  with the absolute eigenvalues as done by SPSS. Default is \code{NULL}.
#' @param signed_loadings logical. If \code{TRUE} (default), the sign of
#' factors with negative sum of loadings is reflected. This is done by both
#' SPSS and \code{\link[psych:fa]{psych::fa}}.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
#'  is given as input. Note that in this case \code{cors} must be set to
#'  \code{FALSE}. Default is "pairwise.complete.obs".
#' @param use_cpp logical. If \code{TRUE}, the iterative PAF procedure to find the
#'  factor solution is performed using Rcpp. This is faster, but can lead to some
#'  very small differences in the output. Default is \code{NULL}.
#'
#' @details Values of \code{init_comm}, \code{criterion}, \code{criterion_type},
#' \code{abs_eigen}, and \code{use_cpp} depend on the \code{type} argument.
#'\code{type = "EFAtools"} will use the following argument specification:
#' \code{init_comm = "smc", criterion = 1e-9, criterion_type = "max_individual",
#' abs_eigen = FALSE, use_cpp = TRUE}.
#' \code{type = "psych"} will use the following argument specification:
#' \code{init_comm = "smc", criterion = .001, criterion_type = "sums",
#' abs_eigen = FALSE, use_cpp = TRUE}.
#' \code{type = "SPSS"} will use the following argument specification:
#' \code{init_comm = "smc", criterion = .001, criterion_type = "max_individual",
#' abs_eigen = TRUE, use_cpp = TRUE}.
#'
#' @return A list of class PAF containing the following
#' \item{orig_R}{Original correlation matrix.}
#' \item{h2_init}{Initial communality estimates.}
#' \item{h2}{Final communality estimates.}
#' \item{iter}{The number of iterations needed for convergence.}
#' \item{orig_eigen}{Eigen values of the original correlation matrix.}
#' \item{init_eigen}{Initial eigenvalues, obtained from the correlation matrix
#'  with the initial communality estimates as diagonal.}
#' \item{final_eigen}{Eigenvalues of the final iteration.}
#' \item{unrot_loadings}{Loading matrix containing the final loadings.}
#' \item{vars_accounted}{Matrix of explained variances and sums of squared loadings}
#' \item{fit_indices}{Common part accounted for (CAF) index as proposed by Lorenzo-Seva,
#' Timmerman, & Kiers (2011) and degrees of freedom}}
#' \item{settings}{list. The settings (arguments) used in the PAF.}
#'
#' @source Grieder, S., & Steiner, M.D.(2019). Algorithmic Jingle Jungle: Comparison of Implementations of an EFA Procedure in R psych Versus SPSS, MacOrtho, and Omega. Submitted Manuscript.
#'
#' @export
#' @examples
#' # call within EFA function:
#' EFA(IDS2_R, n_factors = 5, type = "EFAtools")
#'
#' # call as single function
#' PAF(IDS2_R, n_factors = 5, type = "EFAtools")
PAF <- function(x, n_factors, cors = TRUE, N = NA, max_iter = NULL,
                type = "EFAtools", init_comm = NULL, criterion = NULL,
                criterion_type = NULL, abs_eigen = NULL,
                signed_loadings = TRUE, use = "pairwise.complete.obs",
                use_cpp = NULL) {

  # create R correlation matrix object, if from data, using
  # pairwise binary correlations
  if (isTRUE(cors)) {
    R <- x

    # test whether a real correlation matrix is used
    if (nrow(R) != ncol(R)) {
      stop("Entered data is no correlation matrix but cors = TRUE. Either set ",
           "cors = FALSE if you entered raw data, or enter a correlation matrix.")
    }

    if (is.null(N)) {
      stop("Argument 'N' is NULL. Either provide N, N = NA, or raw data.")
    }

  } else {
    R <- stats::cor(x, use = use)
    colnames(R) <- colnames(x)
    N <- nrow(x)
  }


  if (is.null(type) || !(type %in% c("EFAtools", "psych", "SPSS"))) {

    # if type is not one of the three valid inputs, throw an error if not
    # all the other necessary arguments are specified.

    if (is.null(init_comm) || is.null(criterion) || is.null(criterion_type) ||
        is.null(abs_eigen) || is.null(signed_loadings) || is.null(max_iter) ||
        is.null(use_cpp)) {
      stop('One of "init_comm", "criterion", "criterion_type", "abs_eigen",
           "max_iter", "signed_loadings", "use_cpp", was NULL and no valid
           "type" was specified. Either use one of "EFAtools", "psych", or "SPSS"
           for type, or specify all other arguments')
    }
  } else if (type == "EFAtools") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (is.null(init_comm)) {
      init_comm <- "smc"
    } else {
      warning("Type and init_comm is specified. init_comm is used with value '",
              init_comm, "'. Results may differ from the specified type")
    }

    if (is.null(criterion)) {
      criterion <- 1e-9
    } else {
      warning("Type and criterion is specified. criterion is used with value '",
              criterion, "'. Results may differ from the specified type")
    }

    if (is.null(criterion_type)) {
      criterion_type <- "max_individual"
    } else {
      warning("Type and criterion_type is specified. criterion_type is used with value '",
              criterion_type, "'. Results may differ from the specified type")
    }

    if (is.null(max_iter)) {
      max_iter <- 1e5
    } else {
      warning("Type and max_iter is specified. max_iter is used with value '",
              max_iter, "'. Results may differ from the specified type")
    }

    if (is.null(abs_eigen)) {
      abs_eigen <- FALSE
    } else {
      warning("Type and abs_eigen is specified. abs_eigen is used with value '",
              abs_eigen, "'. Results may differ from the specified type")
    }


    if (is.null(use_cpp)) {
      use_cpp <- TRUE
    } else {
      warning("Type and use_cpp is specified. use_cpp is used with value '",
              use_cpp, "'. Results may differ from the specified type")
    }


  } else if (type == "psych") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (is.null(init_comm)) {
      init_comm <- "smc"
    } else {
      warning("Type and init_comm is specified. init_comm is used with value '",
              init_comm, "'. Results may differ from the specified type")
    }

    if (is.null(criterion)) {
      criterion <- .001
    } else {
      warning("Type and criterion is specified. criterion is used with value '",
              criterion, "'. Results may differ from the specified type")
    }

    if (is.null(criterion_type)) {
      criterion_type <- "sums"
    } else {
      warning("Type and criterion_type is specified. criterion_type is used with value '",
              criterion_type, "'. Results may differ from the specified type")
    }

    if (is.null(max_iter)) {
      max_iter <- 50
    } else {
      warning("Type and max_iter is specified. max_iter is used with value '",
              max_iter, "'. Results may differ from the specified type")
    }

    if (is.null(abs_eigen)) {
      abs_eigen <- FALSE
    } else {
      warning("Type and abs_eigen is specified. abs_eigen is used with value '",
              abs_eigen, "'. Results may differ from the specified type")
    }

    if (is.null(use_cpp)) {
      use_cpp <- TRUE
    } else {
      warning("Type and use_cpp is specified. use_cpp is used with value '",
              use_cpp, "'. Results may differ from the specified type")
    }


  } else if (type == "SPSS") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (is.null(init_comm)) {
      init_comm <- "smc"
    } else {
      warning("Type and init_comm is specified. init_comm is used with value '",
              init_comm, "'. Results may differ from the specified type")
    }

    if (is.null(criterion)) {
      criterion <- .001
    } else {
      warning("Type and criterion is specified. criterion is used with value '",
              criterion, "'. Results may differ from the specified type")
    }

    if (is.null(criterion_type)) {
      criterion_type <- "max_individual"
    } else {
      warning("Type and criterion_type is specified. criterion_type is used with value '",
              criterion_type, "'. Results may differ from the specified type")
    }

    if (is.null(max_iter)) {
      max_iter <- 25
    } else {
      warning("Type and max_iter is specified. max_iter is used with value '",
              max_iter, "'. Results may differ from the specified type")
    }

    if (is.null(abs_eigen)) {
      abs_eigen <- TRUE
    } else {
      warning("Type and abs_eigen is specified. abs_eigen is used with value '",
              abs_eigen, "'. Results may differ from the specified type")
    }

    if (is.null(use_cpp)) {
      use_cpp <- TRUE
    } else {
      warning("Type and use_cpp is specified. use_cpp is used with value '",
              use_cpp, "'. Results may differ from the specified type")
    }

  }

  # set initial communality estimates. This can be done in three ways:
  #  init_comm == "smc": uses the Squared Multiple Correlation
  #  init_comm == "unity": uses unity as inital estimates
  #  init_comm == "mac": uses Maximum Absolute Correlations
  if (init_comm == "smc") {

    # compute the inverse of R
    inverse_R <- solve(R)

    # compute and print the initial communality estimates
    h2_init <- 1 - 1 / diag(inverse_R)

    # save original correlation matrix
    orig_R <- R

    # set diagonal of R to the initial communality estimates
    diag(R) <- h2_init

  } else if (init_comm == "unity") {
    # create h2_init object with a vector of 1s
    h2_init <- diag(R)

    # save original correlation matrix
    orig_R <- R

  } else if (init_comm == "mac") {

    # save original correlation matrix
    orig_R <- R

    # avoid using the diagonal as maximum correlations
    diag(R) <- 0

    # get maximum absolute correlations
    h2_init <- apply(R, 1, function(x){
      max(abs(x))
    })

    # set diagonal of R to the initial communality estimates
    diag(R) <- h2_init

  }

  # save initial eigenvalues
  init_eigen <- eigen(R, symmetric = TRUE)$values

  # define the number of factors m
  m <- n_factors

  ### now we set up a loop to iterate through the procedure until the criterion
  #   is reached

  # counter for the number of iterations
  iter <- 1

  # set communalities to init_communalities for comparison in first iteration
  h2 <- h2_init

  # set initial delta such that at least one iteration is performed
  delta <- 1

  if (isFALSE(use_cpp)) {
    # iterative PAF
    while (delta > criterion && iter <= max_iter){

      # compute the eigenvalues and eigenvectors using the R eigen function
      eigen_list <- eigen(R, symmetric = TRUE)

      # for clarity, store the m eigenvalues and
      # m eigenvectors separately
      Lambda <- eigen_list$values[1:m]
      V <- eigen_list$vectors[, 1:m]

      if (isFALSE(abs_eigen)) {

        if(any(Lambda < 0)){
          stop("Negative Eigenvalues detected; cannot compute communality estimates.
               Try again with init_comm = 'unity' or 'mac'")
        }

        # compute the loadings from the eigenvector matrix and diagonal
        # eigenvalue matrix
        if (m > 1) {
          L <- V %*% diag(sqrt(Lambda))
        } else {
          L <- V * sqrt(Lambda)
        }

        # get the new communality estimates from the loadings
        new_h2 <- diag(L %*% t(L))

      } else if (isTRUE(abs_eigen)) {
        # this is the implementation according to SPSS, which uses
        # absolute eigenvalues to avoid problems when having to compute
        # the square root

        # compute the loadings from the eigenvector matrix and diagonal
        # eigenvalue matrix
        if (m > 1) {
          L <- V %*% diag(sqrt(abs(Lambda)))
        } else {
          L <- V * sqrt(abs(Lambda))
        }
        # get the new communality estimates from the loadings
        # in SPSS implemented as rowSums(V^2 %*% diag(abs(Lambda)))
        # which is equivalent to
        new_h2 <- diag(L %*% t(L))

      }


      if (criterion_type == "max_individual") {
        # save the maximum change in the communality estimates
        delta <- max(abs(h2 - new_h2))
      } else if (criterion_type == "sums"){
        # convergence criterion according to the psych package
        delta <- abs(sum(h2) - sum(new_h2))
      }

      # update diagonal of R with the new communality estimates
      diag(R) <- new_h2

      # update old communality estimates with new ones
      h2 <- new_h2

      # incerase iterator
      iter <- iter + 1

    }

    # break if after maximum iterations there was no convergence
    if (iter >= max_iter){
      warning("Reached maximum number of iterations without convergence.
                Results may not be interpretable.")
    }

  } else if (isTRUE(use_cpp)) {

    crit_type <- ifelse(criterion_type == "max_individual", 1, 2)

    # run the iterative PAF procedure using Rcpp
    L_list <- paf_iter(h2 = h2, criterion = criterion, R = R, n_fac = m,
                       abs_eig = abs_eigen, crit_type = crit_type,
                       max_iter = max_iter)

    h2 <- as.vector(L_list$h2)
    R <- L_list$R
    iter <- L_list$iter
    L <- L_list$L

  }

  if (signed_loadings) {
    # reverse the sign of loadings as done in the psych package,
    # and spss
    if (m > 1) {
      signs <- sign(colSums(L))
      signs[signs == 0] <- 1
      L <- L %*% diag(signs)
    } else {
      if (sum(L) < 0) {
        L <- -as.matrix(L)
      } else {
        L <- as.matrix(L)
      }

    }

  }

  if (!is.null(colnames(orig_R))) {
    # name the loading matrix so the variables can be identified
    rownames(L) <- colnames(orig_R)
  }

  colnames(L) <- paste0("F", 1:m)

  vars_accounted <- .compute_vars(L_unrot = L, L_rot = L)

  colnames(vars_accounted) <- colnames(L)

  # compute CAF and degrees of freedom
  delta_hat <- orig_R - (L %*% t(L))
  diag(delta_hat) <- 1
  CAF <- 1 - KMO(delta_hat)$KMO

  df <- ((m - nrow(L))**2 - (m + nrow(L))) / 2

  fit_ind <- list(CAF, df)

  # create the output object
  class(L) <- "LOADINGS"

  # store the settings used:

  settings <- list(
    N = N,
    max_iter = max_iter,
    type = type,
    init_comm = init_comm,
    criterion = criterion,
    criterion_type = criterion_type,
    abs_eigen = abs_eigen,
    signed_loadings = signed_loadings,
    use = use
  )


  output <- list(
    orig_R = orig_R,
    h2_init = h2_init,
    h2 = h2,
    iter = iter,
    orig_eigen = eigen(orig_R, symmetric = TRUE)$values,
    init_eigen = init_eigen,
    final_eigen = eigen(R, symmetric = TRUE)$values,
    unrot_loadings = L,
    vars_accounted = vars_accounted,
    fit_indices = fit_ind,
    settings = settings
  )

  class(output$h2_init) <- "COMMUNALITIES"
  class(output$h2) <- "COMMUNALITIES"
  class(output$orig_eigen) <- "EIGEN"
  class(output$init_eigen) <- "EIGEN"
  class(output$final_eigen) <- "EIGEN"

  class(output) <- "PAF"

  output

}
