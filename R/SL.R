#' Schmid-Leiman Transformation
#'
#' This function implements the Schmid-Leiman (SL) transformation
#' (Schmid & Leiman, 1957). It takes the pattern coefficients and factor
#' intercorrelations from an oblique factor solution as
#' input and can reproduce the results from [psych::schmid()]
#' and from the SPSS implementation from Wolff & Preising (2005). Other arguments
#' from [EFA()] can be used to control the procedure to find the
#' second-order loadings more flexibly. The function can also be used on a
#' second-order confirmatory factor analysis (CFA) solution from lavaan.
#'
#' @param x object of class [EFA()], class [psych::fa()],
#' class `lavaan::lavaan()` or matrix. If class [EFA()] or
#' class [psych::fa()], pattern coefficients and factor
#' intercorrelations are taken from this object. If class `lavaan::lavaan()`,
#' it must be a second-order CFA solution. In this case first-order and second-order
#'  factor loadings are taken from this object and the `g_name` argument has
#'  to be specified.
#' x can also be a pattern matrix from an oblique factor solution (see `Phi`)
#' or a matrix of first-order factor loadings from a higher-order confirmatory factor
#' analysis (see `L2`).
#' @param Phi matrix. A matrix of factor intercorrelations from an oblique factor
#' solution. Only needs to be specified if a pattern matrix is entered directly
#' into `x`.
#' @param type character. One of "EFAtools" (default), "psych", "SPSS", or "none".
#' This is used to control the procedure of the second-order factor analysis. See
#' [EFA()] for details.
#' @param method character. One of "PAF", "ML", or "ULS" to use
#' principal axis factoring, maximum likelihood, or unweighted least squares,
#' respectively, used in [EFA()] to find the second-order loadings. "MINRES" is
#' accepted as a synonym for "ULS" (the same estimator).
#' @param g_name character. The name of the general factor. This needs only be
#' specified if `x` is a `lavaan` second-order solution. Default is "g".
#' @param ... Arguments to be passed to [EFA()].
#'
#' @details
#' The SL transformation (also called SL orthogonalization) is a procedure with
#' which an oblique factor solution is transformed into a hierarchical,
#' orthogonalized solution. As a first step, the factor intercorrelations are
#' factor analyzed to extract a single second-order (general) factor, yielding a
#' two-level hierarchical structure. The first-order factor and the second-order
#' factor are then orthogonalized, resulting in an orthogonalized factor solution
#' with proportionality constraints. The procedure thus makes a suggested
#' hierarchical data structure based on factor intercorrelations explicit. One
#' major advantage of SL transformation is that it enables variance
#' partitioning between higher-order and first-order factors, including the
#' calculation of McDonald's omegas (see [OMEGA()]).
#'
#' @return A list of class SL containing the following
#' \item{orig_R}{Original correlation matrix.}
#' \item{sl}{A matrix with general factor loadings, group factor loadings, communalities,
#' and uniquenesses.}
#' \item{L2}{Second-order factor loadings.}
#' \item{vars_accounted}{A matrix of explained variances and sums of squared loadings.}
#' \item{iter}{The number of iterations needed for convergence in EFA.}
#' \item{convergence}{Integer convergence code of the second-order EFA (0 =
#' converged); `NA` for a lavaan input. See [EFA()].}
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
#' if (requireNamespace("lavaan", quietly = TRUE)) {
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
#' }
SL <- function(x, Phi = NULL, type = c("EFAtools", "psych", "SPSS", "none"),
               method = c("PAF", "ML", "ULS", "MINRES"), g_name = "g", ...) {

  # Perform argument checks
  checkmate::assert_matrix(Phi, null.ok = TRUE)
  type <- match.arg(type)
  method <- match.arg(method)
  # "MINRES" is a synonym for "ULS" (same estimator); resolve to the canonical name.
  if (method == "MINRES") method <- "ULS"
  checkmate::assert_string(g_name)

  if(!inherits(x, c("EFA", "fa", "lavaan", "matrix", "LOADINGS", "loadings"))){

    cli::cli_abort("{.arg x} must be an {.cls EFA}, {.cls fa}, or {.cls lavaan} object, a matrix, or a {.cls LOADINGS}/{.cls loadings} object.",
                   class = "efa_sl_bad_input")

  }

  if(inherits(x, "EFA")) {

    if("Phi" %in% names(x)){

      L1 <- x$rot_loadings
      n_first_fac <- ncol(x$rot_loadings)
      orig_R <- x$orig_R

      if(!is.null(Phi)){
        cli::cli_warn(
          c("{.arg Phi} is specified; the supplied factor intercorrelations are used.",
            "i" = "To use the intercorrelations from the EFA output, leave {.code Phi = NULL}."),
          class = "efa_sl_phi_specified"
        )

      } else {

        Phi <- x$Phi

      }

    } else {

      cli::cli_abort("{.arg x} is a non-rotated or orthogonal factor solution, but SL needs an oblique solution.",
                     class = "efa_sl_not_oblique")

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
        cli::cli_warn(
          c("{.arg Phi} is specified; the supplied factor intercorrelations are used.",
            "i" = "To use the intercorrelations from the {.fn psych::fa} output, leave {.code Phi = NULL}."),
          class = "efa_sl_phi_specified"
        )

      } else {

        Phi <- x$Phi

      }

    } else {

      cli::cli_abort("{.arg x} is a non-rotated or orthogonal factor solution, but SL needs an oblique solution.",
                     class = "efa_sl_not_oblique")

    }

    n_order <- suppressWarnings(order(as.numeric(gsub("[^0-9]", "", colnames(L1)))))
    L1 <- L1[, n_order]
    Phi <- Phi[n_order, n_order]

  } else if(inherits(x, "lavaan")){

    .require_lavaan()

    if(lavaan::lavInspect(x, what = "converged") == FALSE){
      cli::cli_abort("The model did not converge; the Schmid-Leiman transformation is not performed.",
                     class = "efa_sl_no_converge")
    }

    std_sol <- suppressWarnings(lavaan::lavInspect(x, what = "std"))

    if(any(is.na(std_sol$lambda))){
      cli::cli_abort("Some loadings are {.val NA} or {.val NaN}; the Schmid-Leiman transformation is not performed.",
                     class = "efa_sl_na_loadings")
    }

    if(any(std_sol$lambda >= 1)){
      cli::cli_abort("A Heywood case was detected (a loading of 1 or larger); the Schmid-Leiman transformation is not performed.",
                     class = "efa_sl_heywood")
    }

    if(any(diag(std_sol$theta) <= 0) || any(diag(std_sol$psi) <= 0)){
      cli::cli_abort("A Heywood case was detected (a variance of 0 or negative); the Schmid-Leiman transformation is not performed.",
                     class = "efa_sl_heywood")
    }

    # Create list with factor and corresponding subtest names
    col_names <- colnames(std_sol$lambda)

    if(!any(col_names %in% g_name)){
      cli::cli_abort(
        c("Could not find the specified general-factor name in the lavaan solution.",
          "i" = "Please check the spelling."),
        class = "efa_sl_g_name"
      )
    }

    # SL needs a second-order CFA: the general factor must load on the first-order
    # factors, which lavaan stores in the `beta` (latent regression) matrix. A
    # bifactor solution has no such structure (`beta` is empty), so direct the
    # user to OMEGA() rather than failing later in the loadings algebra.
    if(is.null(std_sol$beta) || !(g_name %in% colnames(std_sol$beta)) ||
       all(std_sol$beta[, g_name] == 0)){
      cli::cli_abort(
        c("{.arg x} does not appear to be a second-order CFA solution.",
          "i" = "{.fn SL} needs a second-order model; pass a bifactor solution directly to {.fn OMEGA}."),
        class = "efa_sl_not_second_order"
      )
    }

    if(!all(std_sol$lambda[, g_name] == 0)){

      cli::cli_warn(
        c("The specified second-order factor contains first-order loadings.",
          "i" = "Did you enter a second-order CFA solution, or the wrong factor name in {.arg g_name}?"),
        class = "efa_sl_second_order_loadings"
      )

    }

    col_names <- col_names[!col_names %in% g_name]
    fac_names <- c(g_name, col_names)

    n_first_fac <- length(col_names)

  } else {

    if(is.null(Phi)){

      cli::cli_abort(
        c("{.arg Phi} was not provided.",
          "i" = "Enter an oblique solution from {.fn EFAtools::EFA} or {.fn psych::fa}, a second-order CFA from lavaan, or provide {.arg Phi}."),
        class = "efa_sl_phi_missing"
      )

    }

    if (!is.null(colnames(x))) {
      n_order <- suppressWarnings(order(as.numeric(gsub("[^0-9]", "", colnames(x)))))
      x <- x[, n_order]
      # Phi is guaranteed non-NULL here (the abort above fires otherwise).
      Phi <- Phi[n_order, n_order]
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

    # Calculate direct group factor loadings (see .sl_group_loadings).
    L_sls_1 <- .sl_group_loadings(L1, std_sol$psi, col_names)

    orig_R <- NA
    iter <- NA
    convergence <- NA
    settings <- NA

  } else {

    # perform a factor analysis on the intercorrelation matrix of the first order
    # factors (N is only specified to avoid a warning)
    EFA_phi <- suppressWarnings(EFA(Phi, n_factors = 1, N = 100, type = type,
                                    method = method, rotation = "none", ...))

    if (ncol(Phi) <= 2) {
      cli::cli_warn("The second-order EFA is underidentified.",
                    class = "efa_sl_underidentified")
    }

    iter <- EFA_phi$iter
    convergence <- EFA_phi$convergence
    settings <- EFA_phi$settings

    # extract second order loadings
    L2 <- EFA_phi$unrot_loadings

    # Schmid-Leiman solution, direct loadings of second order factor
    L_sls_2 <- L1 %*% L2

    # Communalities of the second-order factor. A value at or above 1 is a
    # Heywood case: the residualized first-order loadings would be undefined
    # (the square root below would be taken of a negative number).
    comm_h <- rowSums(L2^2)

    if(any(comm_h >= 1 + .Machine$double.eps)){
      cli::cli_abort(
        c("A Heywood case was detected in the second-order factor analysis; no Schmid-Leiman solution is computed.",
          "i" = "A second-order communality is 1 or larger, so the residualized first-order loadings are undefined."),
        class = "efa_sl_heywood"
      )
    }

    # compute uniqueness of higher order factor. A communality a hair above 1
    # (within floating-point noise, below the Heywood abort threshold above) is
    # clamped to 1 so the square root never sees a tiny negative value.
    u2_h <- sqrt(pmax(0, 1 - comm_h))

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
    convergence = convergence,
    settings = settings
    )

  class(output) <- "SL"

  output

}
