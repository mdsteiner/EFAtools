#' Schmid-Leiman Transformation
#'
#' This function implements the Schmid-Leiman transformation. It takes the pattern
#' coefficients and factor intercorrelations from an oblique factor solution as
#' input and can reproduce the results from \code{\link[psych:schmid]{psych::schmid}},
#' from the SPSS
#' implementation from Wolff & Preising (2005), and from MacOrtho (Watkins, 2004).
#' To reproduce
#' psych or SPSS, only the type argument has to be specified additional to the
#' loadings and factor intercorrelations. Other arguments from \code{\link{PAF}} can
#' be used to control the procedure to find the second order loadings more
#' flexibly.
#'
#' @param x object of class \code{\link{PROMAX}} or class \code{\link{fa}} or
#' matrix. If class \code{\link{PROMAX}} or class \code{\link{fa}},
#' pattern coefficients and factor intercorrelations are taken from this object.
#' x can also be a pattern matrix from an oblique factor solution (see \code{Phi}).
#' @param Phi matrix. A matrix of factor intercorrelations from an oblique factor
#' solution. Only needs to be specified if x is not of class \code{\link{PROMAX}}
#' or class \code{\link{fa}}.
#' @param type character. One of "EFAdiff" (default), "psych", or "SPSS". This
#' is used to control the procedure of the second order factor analysis. See
#' \code{\link{PAF}} for details.
#' @param ... Arguments to be passed to \code{PAF}.
#'
#' @return A list containing the following
#' \item{sl}{A matrix with g loadings, group factor loadings, communalities,
#' and uniquenesses.}
#' \item{L2}{Second-order factor loadings.}
#'
#' @export
SL <- function(x, Phi = NULL, type = "EFAdiff", ...) {

  if(all(class(x) == "PROMAX")) {

    L1 <- x$rot_loadings
    n_first_fac <- ncol(x$rot_loadings)

    if(!is.null(Phi)){
      warning("Phi argument is specified. Specified factor intercorrelations are
              taken. To take factor intercorrelations from the PROMAX output,
              leave Phi = NULL")

    } else {

      Phi <- x$Phi

    }

  } else if(all(class(x) == c("psych", "fa"))) {

    L1 <- unclass(x$loadings)
    n_first_fac <- ncol(x$loadings)

    if(!is.null(Phi)){
      warning("Phi argument is specified. Specified factor intercorrelations are
              taken. To take factor intercorrelations from the psych fa output,
              leave Phi = NULL")

    } else {

      Phi <- x$Phi

    }

  } else {

    L1 <- x
    n_first_fac <- ncol(x)

  }

  # perform a factor analysis on the intercorrelation matrix of the first order
  # factors
  paf_phi <- PAF(Phi, n_factors = 1, type = type, ...)

  # extract second order loadings
  L2 <- paf_phi$loadings

  # Schmid-Leiman solution, direct loadings of second order factor
  L_sls_2 <- L1 %*% L2

  # compute uniqueness of higher order factor
  u2_h <- sqrt(1 - diag(L2 %*% t(L2)))

  # Schmid-Leiman solution, residualized first order factor loadings
  L_sls_1 <- L1 %*% diag(u2_h)

  # Combine the Schmid-Leiman loadings in a data frame
  sl_load <- cbind(L_sls_2, L_sls_1)

  # Compute communalities and uniquenesses of the Schmid-Leiman solution
  h2_sl <- rowSums(sl_load^2)
  u2_sl <- 1 - h2_sl

  # Finalize output object
  sl <- cbind(sl_load, h2_sl, u2_sl)
  colnames(sl) <- c("g", 1:n_first_fac, "h2", "u2")

  output <- list(
    sl = sl,
    L2 = L2)

  class(output) <- "SL"

  output

}
