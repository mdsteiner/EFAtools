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
#' @param ... Further arguments passed to \link{PARALLEL}.
#'
#' @details The \link{PARALLEL} function and the principal axis factoring of the
#'   different number of factors can be parallelized using the future framework,
#'   by calling the \link{future}{plan} function. The examples provide example code
#'   on how to enable parallel processing.
#'
#' @return A list of class HULL containing the following objects
#' \item{n_factors}{The number of factors to retain according to the Hull method.}
#' \item{solutions}{A matrix containing the goodness of fit indices (CAF), degrees of freedom, and for the factors lying on the hull, the st values of the hull solution (see Lorenzo-Seva, Timmereman, and Kiers 2011 for details).}
#'
#' @source Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011). The Hull method for selecting the number of common factors. Multivariate behavioral research, 46(2), 340-364.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' HULL(IDS2_R, n_cases = 2000)
#'
#' # using parallel processing (Note: plans can be adapted, see the future
#' # package for details)
#' future::plan(future::multisession)
#' HULL(IDS2_R, n_cases = 2000)
#' }
HULL <- function(x, cors = TRUE, n_cases = NA, n_factors = NA, ...) {
  # Perform hull method following Lorenzo-Seva, Timmereman, and Kiers (2011)

  if (!isTRUE(cors)) {
    n_cases <- nrow(x)
    R <- stats::cor(x, use = "pairwise")
  } else {
    R <- as.matrix(x)
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
  colnames(s) <- c("n factors", "f", "df", "st")

  # first for 0 factors
  s[1, 2] <- 1 - KMO(R)$KMO
  s[1, 3] <- ((m)**2 - (m)) / 2


  pafs <- future.apply::future_lapply(1:J, hull_paf, R = R, criterion = .001,
                                      max_iter = 1e4)

  # then for 1 to J factors
  for (i in 1:J) {
    # compute goodness of fit "f" as CAF (common part accounted for; Eq 3)
    A_i <- pafs[[i]]
    delta_hat <- R - (A_i %*% t(A_i))
    diag(delta_hat) <- 1
    s[i + 1, 2] <- 1 - KMO(delta_hat)$KMO

    # compute dfs (Eq 4 provides the number of free parameters; using dfs yields
    # th same numbers, as the difference in df equals the difference in free parameters)
    s[i + 1, 3] <- ((m - i)**2 - (m + i)) / 2

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
    solutions = s_complete
  )

  class(out) <- "HULL"

  return(out)
}
