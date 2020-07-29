#' Schmid-Leiman Transformation
#'
#' This function implements the Schmid-Leiman (SL) transformation
#' (Schmid & Leiman, 1957). It takes the pattern coefficients and factor
#' intercorrelations from an oblique factor solution as
#' input and can reproduce the results from \code{\link[psych:schmid]{psych::schmid}}
#' and from the SPSS implementation from Wolff & Preising (2005). Other arguments
#' from \code{\link{EFA}} can be used to control the procedure to find the
#' second-order loadings more flexibly. The function can also be used on a
#' second-order confirmatory factor analysis (CFA) solution from lavaan.
#'
#' @param x object of class \code{\link{EFA}}, class \code{\link[psych:fa]{psych::fa}},
#' class \code{\link[lavaan]{lavaan}} or matrix. If class \code{\link{EFA}} or
#' class \code{\link[psych:fa]{psych::fa}}, pattern coefficients and factor
#' intercorrelations are taken from this object. If class \code{\link[lavaan]{lavaan}},
#' it must be a second-order CFA solution. In this case first-order and second-order
#'  factor loadings are taken from this object and the \code{g_name} argument has
#'  to be specified.
#' x can also be a pattern matrix from an oblique factor solution (see \code{Phi})
#' or a matrix of first-order factor loadings from a higher-order confirmatory factor
#' analysis (see \code{L2}).
#' @param Phi matrix. A matrix of factor intercorrelations from an oblique factor
#' solution. Only needs to be specified if a pattern matrix is entered directly
#' into \code{x}.
#' @param type character. One of "EFAtools" (default), "psych", "SPSS", or "none".
#' This is used to control the procedure of the second-order factor analysis. See
#' \code{\link{EFA}} for details.
#' @param method character. One of "PAF", "ML", or "ULS" to use
#' principal axis factoring, maximum likelihood, or unweighted least squares
#' (also called minres), respectively, used in \code{\link{EFA}} to find the second-order
#' loadings.
#' @param g_name character. The name of the general factor. This needs only be
#' specified if \code{x} is a \code{lavaan} second-order solution. Default is "g".
#' @param ... Arguments to be passed to \code{\link{EFA}}.
#'
#' @details
#' The SL transformation (also called SL orthogonalization) is a procedure with
#' which an oblique factor solution is transformed into a hierarchical,
#' orthogonalized solution. As a first step, the factor intercorrelations are
#' again factor analyzed to find second-order factor loadings. If there is only
#' one higher-order factor, this step of the procedure stops there, resulting in
#' a second-order factor structure. The first-order factor and the second-order
#' factor are then orthogonalized, resulting in an orthogonalized factor solution
#' with proportionality constraints. The procedure thus makes a suggested
#' hierarchical data structure based on factor intercorrelations explicit. One
#' major advantage of SL transformation is that it enables variance
#' partitioning between higher-order and first-order factors, including the
#' calculation of McDonald's omegas (see \code{\link{OMEGA}}).
#'
#' @return A list of class SL containing the following
#' \item{orig_R}{Original correlation matrix.}
#' \item{sl}{A matrix with general factor loadings, group factor loadings, communalities,
#' and uniquenesses.}
#' \item{L2}{Second-order factor loadings.}
#' \item{vars_accounted}{A matrix of explained variances and sums of squared loadings.}
#' \item{iter}{The number of iterations needed for convergence in EFA.}
#' \item{settings}{list. The settings (arguments) used in EFA to get the
#' second-order loadings.}
#'
#' @source Schmid, J. & Leiman, J. M. (1957). The development of hierarchical
#' factor solutions. Psychometrika, 22(1), 53–61. doi:10.1007/BF02289209
#' @source Wolff, H.-G., & Preising, K. (2005). Exploring item and higher order
#' factor structure with the Schmid-Leiman solution: Syntax codes for SPSS and
#' SAS. Behavior Research Methods, 37 , 48–58. doi:10.3758/BF03206397
#'
#' @export
#'
#' @examples
#' ## Use with an output from the EFAtools::EFA function, both with type EFAtools
#' EFA_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                type = "EFAtools", method = "PAF", rotation = "promax")
#' SL_EFAtools <- SL(EFA_mod, type = "EFAtools", method = "PAF")
#'
#' \donttest{
#' ## Use with an output from the psych::fa function with type psych in SL
#' fa_mod <- psych::fa(test_models$baseline$cormat, nfactors = 3, n.obs = 500,
#'                     fm = "pa", rotate = "Promax")
#' SL_psych <- SL(fa_mod, type = "psych", method = "PAF")
#' }
#'
#' ## Use more flexibly by entering a pattern matrix and phi directly (useful if
#' ## a factor solution found with another program should be subjected to SL
#' ## transformation)
#'
#' ## For demonstration, take pattern matrix and phi from an EFA output
#' ## This gives the same solution as the first example
#' EFA_mod <- EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
#'                type = "EFAtools", method = "PAF", rotation = "promax")
#' SL_flex <- SL(EFA_mod$rot_loadings, Phi = EFA_mod$Phi, type = "EFAtools",
#'               method = "PAF")
#'
#' \donttest{
#' ## Use with a lavaan second-order CFA output
#'
#' # Create and fit model in lavaan (assume all variables have SDs of 1)
#' mod <- 'F1 =~ V1 + V2 + V3 + V4 + V5 + V6
#'         F2 =~ V7 + V8 + V9 + V10 + V11 + V12
#'         F3 =~ V13 + V14 + V15 + V16 + V17 + V18
#'         g =~ F1 + F2 + F3'
#' fit <- lavaan::cfa(mod, sample.cov = test_models$baseline$cormat,
#'                    sample.nobs = 500, estimator = "ml")
#'
#' SL_lav <- SL(fit, g_name = "g")
#'
#' }
SL <- function(x, Phi = NULL, type = c("EFAtools", "psych", "SPSS", "none"),
               method = c("PAF", "ML", "ULS"), g_name = "g", ...) {

  # Perform argument checks
  checkmate::assert_matrix(Phi, null.ok = TRUE)
  type <- match.arg(type)
  method <- match.arg(method)
  checkmate::assert_string(g_name)

  if(!inherits(x, c("EFA", "fa", "lavaan", "matrix", "LOADINGS", "loadings"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither an object of class EFA, fa, or lavaan, nor a matrix, nor of class LOADINGS or loadings.\n"))

  }

  if(inherits(x, "EFA")) {

    if("Phi" %in% names(x)){

      L1 <- x$rot_loadings
      n_first_fac <- ncol(x$rot_loadings)
      orig_R <- x$orig_R

      if(!is.null(Phi)){
        warning(crayon::yellow$bold("!"), crayon::yellow(" Phi argument is specified. Specified factor intercorrelations are taken. To take factor intercorrelations from the EFA output, leave Phi = NULL\n"))

      } else {

        Phi <- x$Phi

      }

    } else {

      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is either a non-rotated or orthogonal factor solution. SL needs an oblique factor solution\n"))

    }

    n_order <- order(as.numeric(gsub("F", "", colnames(L1))))
    L1 <- L1[, n_order]
    Phi <- Phi[n_order, n_order]

  } else if(inherits(x, "fa")) {

    if("Phi" %in% names(x)){

      L1 <- unclass(x$loadings)
      n_first_fac <- ncol(x$loadings)
      orig_R <- unclass(x$r)

      if(!is.null(Phi)){
        warning(crayon::yellow$bold("!"), crayon::yellow(" Phi argument is specified. Specified factor intercorrelations are taken. To take factor intercorrelations from the psych fa output, leave Phi = NULL\n"))

      } else {

        Phi <- x$Phi

      }

    } else {

      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is either a non-rotated or orthogonal factor solution. SL needs an oblique factor solution\n"))

    }

    n_order <- order(colnames(L1))
    L1 <- L1[, n_order]
    Phi <- Phi[n_order, n_order]

  } else if(inherits(x, "lavaan")){

    if(lavaan::lavInspect(x, what = "converged") == FALSE){
      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Model did not converge. No omegas are computed.\n"))
    }

    std_sol <- suppressWarnings(lavaan::lavInspect(x, what = "std"))

    if(any(is.na(std_sol$lambda))){
      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Some loadings are NA or NaN. No omegas are computed.\n"))
    }

    if(any(std_sol$lambda >= 1)){
      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" A Heywood case was detected (loading equal to or larger than 1). No omegas are computed.\n"))
    }

    # Create list with factor and corresponding subtest names
    col_names <- colnames(std_sol$lambda)

    if(!any(col_names %in% g_name)){
      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Could not find the specified name of the general factor in the entered lavaan solution. Please check the spelling.\n"))
    }

    if(!all(std_sol$lambda[, g_name] == 0)){

      warning(crayon::yellow$bold("!"), crayon::yellow(" The second-order factor you specified contains first-order loadings. Did you really enter a second-order CFA solution? Or did you enter the wrong factor name in g_name?\n"))

    }

    col_names <- col_names[!col_names %in% g_name]
    fac_names <- c(g_name, col_names)

    n_first_fac <- length(col_names)

  } else {

    if(is.null(Phi)){

      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" Phi not provided. Either enter an oblique factor solution from EFAtools::EFA or from psych::fa, or a second-order CFA solution from lavaan, or provide Phi\n"))

    }

    if (!is.null(colnames(x))) {
      n_order <- order(colnames(x))
      x <- x[, n_order]

      if(!is.null(Phi)){
      Phi <- Phi[n_order, n_order]
      }

    }

    L1 <- x
    n_first_fac <- ncol(x)
    orig_R <- NA

  }

  if(inherits(x, "lavaan")){

    # Calculate direct g loadings
    L1 <- std_sol$lambda[, col_names]
    L2 <- std_sol$beta[col_names, g_name]

    L_sls_2 <- L1 %*% L2

    # Calculate direct group factor loadings
    L_sls_1 <- L1 %*% sqrt(std_sol$psi[col_names, col_names])

    orig_R <- NA
    iter <- NA
    settings <- NA

  } else {

    # perform a factor analysis on the intercorrelation matrix of the first order
    # factors (N is only specified to avoid a warning)
    EFA_phi <- suppressWarnings(EFA(Phi, n_factors = 1, N = 100, type = type,
                                    method = method, rotation = "none", ...))

    iter <- EFA_phi$iter
    settings <- EFA_phi$settings

    # extract second order loadings
    L2 <- EFA_phi$unrot_loadings

    # Schmid-Leiman solution, direct loadings of second order factor
    L_sls_2 <- L1 %*% L2

    # compute uniqueness of higher order factor
    u2_h <- sqrt(1 - diag(L2 %*% t(L2)))

    # Schmid-Leiman solution, residualized first order factor loadings
    L_sls_1 <- L1 %*% diag(u2_h)

  }

  # Combine the Schmid-Leiman loadings in a data frame
  sl_load <- cbind(L_sls_2, L_sls_1)

  # Compute communalities and uniquenesses of the Schmid-Leiman solution
  h2_sl <- rowSums(sl_load^2)
  u2_sl <- 1 - h2_sl

  vars_accounted <- .compute_vars(L_unrot = sl_load, L_rot = sl_load)

  colnames(vars_accounted) <-c("g", paste0("F", seq_len(n_first_fac)))

  # Finalize output object
  sl <- cbind(sl_load, h2_sl, u2_sl)
  colnames(sl) <- c("g", paste0("F", seq_len(n_first_fac)), "h2", "u2")
  class(sl) <- "SLLOADINGS"

  output <- list(
    orig_R = orig_R,
    sl = sl,
    L2 = L2,
    vars_accounted = vars_accounted,
    iter = iter,
    settings = settings
    )

  class(output) <- "SL"

  output

}
