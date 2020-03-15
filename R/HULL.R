#' Hull method for determining the number of factors to retain
#'
#' Implementation of the Hull method suggested by Lorenzo-Seva, Timmereman,
#' and Kiers (2011), but using principal axis factoring. See details for
#' parallelization.
#'
#' @param x matrix or data.frame. Correlation matrix or raw data.
#' @param cors logical. Whether x is a correlation matrix. Default is set to TRUE.
#' @param n_cases numeric. Number of cases in the data. This is passed to \link{PARALLEL}.
#'  Only has to be specified if x is a correlation matrix.
#' @param n_factors numeric. Theoretical number of factors to retain. The maximum
#'   of this number and the number of factors suggested by \link{PARALLEL} olus
#'   onewill be used in the Hull method.
#' @param method character. The estimation method to use. One of  \code{"PAF"},
#'    \code{"ULS"}, or  \code{"ML"}, for principal axis factoring, unweighted
#'    least squares, and maximum likelihood, respectively.
#' @param gof character. The goodness of fit index to use. One of \code{"CAF"},
#'   \code{"CFI"}, or \code{"RMSEA"}. If \code{method = "PAF"} is used, only
#'   the CAF can be used as goodness of fit index. For details on the CAF, see
#'   Lorenzo-Seva, Timmerman, and Kiers (2011).
#' @param ... Further arguments passed to \link{PARALLEL}.
#'
#' @details The \link{PARALLEL} function and the principal axis factoring of the
#'   different number of factors can be parallelized using the future framework,
#'   by calling the \link{future}{plan} function. The examples provide example code
#'   on how to enable parallel processing.
#'
#'   Note that if \code{gof = "RMSEA"} is used, 1 - RMSEA is actually used to
#'   compare the different solutions. Thus, the threshold of .05 is now .95.
#'
#'   The ML estimation method uses the \link{stats}{factanal} starting values. See
#'   also the \link{ML} documentation.
#'
#' @return A list of class HULL containing the following objects
#' \item{n_factors}{The number of factors to retain according to the Hull method.}
#' \item{solutions}{A matrix containing the goodness of fit indices (CAF), degrees of freedom, and for the factors lying on the hull, the st values of the hull solution (see Lorenzo-Seva, Timmereman, and Kiers 2011 for details).}
#' \item{ctrl}{A list of control settings used in the print and plot methods.}
#'
#' @source Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011). The Hull method for selecting the number of common factors. Multivariate behavioral research, 46(2), 340-364.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # using PAF
#' HULL(IDS2_R, n_cases = 2000)
#'
#' # using ML with CFI
#' HULL(IDS2_R, n_cases = 2000, method = "ML", gof = "CFI")
#'
#' # using ULS with RAMSEA
#' HULL(IDS2_R, n_cases = 2000, method = "ULS", gof = "RAMSEA")
#'
#' # using parallel processing (Note: plans can be adapted, see the future
#' # package for details)
#' future::plan(future::multisession)
#' HULL(IDS2_R, n_cases = 2000)
#' }
HULL <- function(x, cors = TRUE, n_cases = NA, n_factors = NA,
                 method = c("PAF", "ULS", "ML"), gof = c("CAF", "CFI", "RMSEA"),
                 ...) {
  # Perform hull method following Lorenzo-Seva, Timmereman, and Kiers (2011)

  # # for testing
  # x <- IDS2_R
  # cors <- TRUE
  # n_cases <- 2000
  # n_factors <- 7
  # method <- "ML"
  # gof <- "CFI"

  method <- match.arg(method)
  gof <- match.arg(gof)

  if (method == "PAF" && gof != "CAF") {
    message('gof is set to "', gof, '" but must be CAF if method "PAF" is used.',
            ' Setting gof to "CAF"')
    gof <- "CAF"
  }

  if (!isTRUE(cors)) {
    n_cases <- nrow(x)
    R <- stats::cor(x, use = "pairwise")
  } else {
    R <- as.matrix(x)
  }

  if (gof != "CAF" && is.na(n_cases)) {
    stop('n_cases is not specified but is needed for computation of ', gof,
         ' fit index.')
  }

  m <- ncol(R)

  # 1) perform parallel analysis to find J as n_factors + 1
  par_res <- PARALLEL(R, n_cases = n_cases, ...)

  if (is.na(par_res$n_factor)) {

    if (!is.na(n_factors)) {
      J <- n_factors + 1
    } else {
      J <- floor(ncol(R) / 2)
    }

  } else {

    J <- max(c(par_res$n_factor, n_factors), na.rm = TRUE) + 1

  }


  # 2) perform factor analysis for the range of dimensions 1:J and compute f and
  #    df for every solution
  s <- matrix(0, ncol = 4, nrow = J + 1)
  s[, 1] <- 0:J
  colnames(s) <- c("n factors", gof, "df", "st")

  # first for 0 factors
  if (gof == "CAF") {
    s[1, 2] <- 1 - KMO(R)$KMO
  } else if (gof == "CFI") {
    s[1, 2] <- 0

    # for later use in loop
    chi_null <- sum(R[upper.tri(R)] ^ 2) * (n_cases - 1)
    df_null <- (m**2 - m) / 2
    delta_hat_null <- max(0, chi_null - df_null)

  } else if (gof == "RMSEA") {

    Fm <- sum(R[upper.tri(R)] ^ 2)
    chi <- Fm * (n_cases - 1)
    df <- (m**2 - m) / 2
    delta_hat_m <- max(0, chi - df)
    # compute 1 - RMSEA
    s[1, 2] <- 1 - sqrt(max(delta_hat_m / (df * (n_cases - 1)), 0))
  }


  s[1, 3] <- (m**2 - m) / 2


  if (method == "PAF") {
    loadings <- future.apply::future_lapply(1:J, hull_paf, R = R, criterion = .001,
                                           max_iter = 1e4)
  } else if (method == "ULS") {
    loadings <- future.apply::future_lapply(1:J, .hull_uls, R = R)
  } else if (method == "ML") {
    loadings <- future.apply::future_lapply(1:J, .hull_ml, R = R)
  }


  # then for 1 to J factors
  for (i in 1:J) {
    df <- ((m - i)**2 - (m + i)) / 2
    if (method == "PAF") {
      # compute goodness of fit "f" as CAF (common part accounted for; Eq 3)
      A_i <- loadings[[i]]
      delta_hat <- R - (A_i %*% t(A_i))
      diag(delta_hat) <- 1
      # compute CAF
      s[i + 1, 2] <- 1 - KMO(delta_hat)$KMO
    } else {
      if (gof == "CAF") {
        # compute goodness of fit "f" as CAF (common part accounted for; Eq 3)
        A_i <- loadings[[i]]$loadings
        delta_hat <- R - (A_i %*% t(A_i))
        diag(delta_hat) <- 1
        # compute CAF
        s[i + 1, 2] <- 1 - KMO(delta_hat)$KMO
      } else if (gof == "CFI") {
        Fm <- loadings[[i]]$fit
        chi <- Fm * (n_cases - 1)
        delta_hat_m <- max(0, chi - df)
        # compute CFI
        s[i + 1, 2] <- 1 - delta_hat_m / delta_hat_null

      } else if (gof == "RMSEA") {
        Fm <- loadings[[i]]$fit
        chi <- Fm * (n_cases -1)
        delta_hat_m <- max(0, chi - df)
        # compute 1 - RMSEA
        s[i + 1, 2] <- 1 - sqrt(max(delta_hat_m / (df * (n_cases - 1)), 0))

      }

    }

    # compute dfs (Eq 4 provides the number of free parameters; using dfs yields
    # th same numbers, as the difference in df equals the difference in free parameters)
    s[i + 1, 3] <- df

  }

  # 3) sort n solutions by their df values and denoted by s (already done)

  # 4) all solutions s are excluded for which a solution sj (j<i) exists such
  #    that fj > fi (eliminate solutions not on the boundary of the convex hull)

  s_complete <- s
  d_s <- diff(s[, 2])
  while (any(d_s < 0)) {
    s <- s[c(1, d_s) > 0,]
    d_s <- diff(s[, 2])
  }

  # 5) all triplets of adjacent solutions are considered consecutively.
  #    middle solution is excluded if its point is below or on the line
  #    connecting its neighbors in GOF vs df

  # 6) repeat 5) until no solution can be excluded

  nr_s <- nrow(s)
  i <- 2

  while(i < nr_s - 1) {

    f1 <- s[i - 1, 2]
    f2 <- s[i, 2]
    f3 <- s[i + 1, 2]
    df1 <- s[i - 1, 3]
    df2 <- s[i, 3]
    df3 <- s[i + 1, 3]

    # compute f2 if it were on the line between f1 and f3
    p_f2 <- f1 + (f3 - f1) / (df3 - df1) * (df2 - df1)

    # check if f2 is below or on the predicted line and if so, remove it
    if (f2 <= p_f2) {
      s <- s[-i, ]
      nr_s <- nr_s -1
      i <- 1
    }
    i <- i + 1
  }


  # 7) the st values of the hull solutions are determined (Eq 5)
  for (i in 2:(nrow(s) - 1)) {

    f_i <- s[i, 2]
    f_p <- s[i - 1, 2]
    f_n <- s[i + 1, 2]
    df_i <- s[i, 3]
    df_p <- s[i - 1, 3]
    df_n <- s[i + 1, 3]

    s[i, 4] <- ((f_i - f_p) / (df_i - df_p)) / ((f_n - f_i) / (df_n - df_i))

  }

  # combine values
  for (row_i in 0:J) {

    if (row_i %in% s[,1]) {
      s_complete[row_i + 1, 4] <- s[s[,1] == row_i, 4]
    } else {
      s_complete[row_i + 1, 4] <- NA
    }

  }


  # 8) select solution with highest st value
  retain <- s[which.max(s[, 4]), 1]


  out <- list(
    n_factors = retain,
    solutions = s_complete,
    ctrl = list(method = method,
                gof = gof,
                n_factors = n_factors)
  )

  class(out) <- "HULL"

  return(out)
}


.hull_uls <- function(n_factors, R) {
  uls <- .fit_uls(R, n_factors)
  return(list("loadings" = uls$loadings, "fit" = uls$res$value))
}

.hull_ml <- function(n_factors, R) {
  ml <- .fit_ml(R, n_factors, "factanal")
  return(list("loadings" = ml$loadings, "fit" = ml$res$value))
}

