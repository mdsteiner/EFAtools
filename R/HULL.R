#' Hull method for determining the number of factors to retain
#'
#' Implementation of the Hull method suggested by Lorenzo-Seva, Timmerman,
#' and Kiers (2011), with an extension to principal axis factoring. See details for
#' parallelization.
#'
#' @param x matrix or data.frame. Dataframe or matrix of raw data or matrix with
#' correlations.
#' @param N numeric. Number of cases in the data. This is passed to \link{PARALLEL}.
#'  Only has to be specified if x is a correlation matrix, otherwise it is determined
#'  based on the dimensions of x.
#' @param n_fac_theor numeric. Theoretical number of factors to retain. The maximum
#'   of this number and the number of factors suggested by \link{PARALLEL} plus
#'   one will be used in the Hull method.
#' @param method character. The estimation method to use. One of  \code{"PAF"},
#'    \code{"ULS"}, or  \code{"ML"}, for principal axis factoring, unweighted
#'    least squares, and maximum likelihood, respectively.
#' @param gof character. The goodness of fit index to use. Either \code{"CAF"},
#'   \code{"CFI"}, or \code{"RMSEA"}, or any combination of them.
#'   If \code{method = "PAF"} is used, only
#'   the CAF can be used as goodness of fit index. For details on the CAF, see
#'   Lorenzo-Seva, Timmerman, and Kiers (2011).
#' @param eigen_type character. On what the eigenvalues should be found in the
#'  parallel analysis. Can be one of \code{"SMC"}, \code{"PCA"}, or \code{"EFA"}.
#'   If using  \code{"SMC"} (default), the diagonal of the correlation matrices is
#'    replaced by the squared multiple correlations (SMCs) of the indicators. If
#'     using  \code{"PCA"}, the diagonal values of the correlation
#'  matrices are left to be 1. If using  \code{"EFA"}, eigenvalues are found on the
#'  correlation  matrices with the final communalities of an EFA solution as
#'  diagonal. This is passed to  \code{\link{PARALLEL}}.
#' @param use character. Passed to \code{\link[stats:cor]{stats::cor}} if raw data
#' is given as input. Default is \code{"pairwise.complete.obs"}.
#' @param cor_method character. Passed to \code{\link[stats:cor]{stats::cor}}.
#'  Default is  \code{"pearson"}.
#' @param n_datasets numeric. The number of datasets to simulate. Default is 1000.
#'   This is passed to \code{\link{PARALLEL}}.
#' @param percent numeric. A vector of percentiles to take the simulated eigenvalues from.
#'  Default is 95. This is passed to \code{\link{PARALLEL}}.
#' @param decision_rule character. Which rule to use to determine the number of
#' factors to retain. Default is \code{"means"}, which will use the average
#' simulated eigenvalues. \code{"percentile"}, uses the percentiles specified
#' in percent. \code{"crawford"} uses the 95th percentile for the first factor
#' and the mean afterwards (based on Crawford et al, 2010). This is passed to \code{\link{PARALLEL}}.
#' @param n_factors numeric. Number of factors to extract if  \code{"EFA"} is
#' included in \code{eigen_type}. Default is 1. This is passed to
#' \code{\link{PARALLEL}}.
#' @param ... Further arguments passed to \code{\link{EFA}}, also in
#' \code{\link{PARALLEL}}.
#'
#' @details The Hull method aims to find a model with an optimal balance between
#'  model fit and number of parameters. That is, it aims to retrieve only major
#'  factors (Lorenzo-Seva, Timmerman, & Kiers, 2011). To this end, it performs
#'  the following steps (Lorenzo-Seva, Timmerman, & Kiers, 2011, p.351):
#'  \enumerate{
#'    \item It performs parallel analysis and adds one to the identified number of factors (this number is denoted \emph{J}). \emph{J} is taken as an upper bound of the number of factors to retain in the hull method. Alternatively, a theoretical number of factors can be entered. In this case \emph{J} will be set to whichever of these two numbers (from parallel analysis or based on theory) is higher.
#'    \item For all 0 to \emph{J} factors, the goodness-of-fit (one of \emph{CAF}, \emph{RMSEA}, or \emph{CFI}) and the degrees of freedom (\emph{df}) are computed.
#'    \item The solutions are ordered according to their \emph{df}.
#'    \item Solutions that are not on the boundary of the convex hull are eliminated (see Lorenzo-Seva, Timmerman, & Kiers, 2011, for details).
#'    \item All the triplets of adjacent solutions are considered consecutively. The middle solution is excluded if its point is below or on the line connecting its neighbors in a plot of the goodness-of-fit versus the degrees of freedom.
#'    \item Step 5 is repeated until no solution can be excluded.
#'    \item The \emph{st} values of the “hull” solutions are determined.
#'    \item The solution with the highest \emph{st} value is selected.
#'  }
#'
#' The \link{PARALLEL} function and the principal axis factoring of the
#'   different number of factors can be parallelized using the future framework,
#'   by calling the \link[future:plan]{future::plan} function. The examples
#'    provide example code on how to enable parallel processing.
#'
#'   Note that if \code{gof = "RMSEA"} is used, 1 - RMSEA is actually used to
#'   compare the different solutions. Thus, the threshold of .05 is then .95. This
#'   is necessary due to how the heuristic to locate the elbow of the hull works.
#'
#'   The ML estimation method uses the \link[stats:factanal]{stats::factanal}
#'    starting values. See also the \link{EFA} documentation.
#'
#'    The \code{HULL} function can also be called together with other factor
#'    retention criteria in the \code{\link{N_FACTORS}} function.
#' @return A list of class HULL containing the following objects
#' \item{n_fac_CAF}{The number of factors to retain according to the Hull method
#' with the CAF.}
#' \item{n_fac_CFI}{The number of factors to retain according to the Hull method
#' with the CFI.}
#' \item{n_fac_RMSEA}{The number of factors to retain according to the Hull method
#' with the RMSEA.}
#' \item{solutions_CAF}{A matrix containing the CAFs, degrees of freedom, and for the factors lying on the hull, the st values of the hull solution (see Lorenzo-Seva, Timmerman, and Kiers 2011 for details).}
#' \item{solutions_CFI}{A matrix containing the CFIs, degrees of freedom, and for the factors lying on the hull, the st values of the hull solution (see Lorenzo-Seva, Timmerman, and Kiers 2011 for details).}
#' \item{solutions_RMSEA}{A matrix containing the RMSEAs, degrees of freedom, and for the factors lying on the hull, the st values of the hull solution (see Lorenzo-Seva, Timmerman, and Kiers 2011 for details).}
#' \item{n_fac_max}{The upper bound \emph{J} of the number of factors to extract (see details).}
#' \item{settings}{A list of the settings used.}
#'
#' @source Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011).
#' The Hull method for selecting the number of common factors. Multivariate
#' Behavioral Research, 46(2), 340-364.
#'
#' @seealso Other factor retention criteria: \code{\link{CD}}, \code{\link{EKC}},
#' \code{\link{KGC}}, \code{\link{PARALLEL}}, \code{\link{SMT}}
#'
#' \code{\link{N_FACTORS}} as a wrapper function for this and all the
#' above-mentioned factor retention criteria.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # using PAF (this will throw a warning if gof is not specified manually
#' # and CAF will be used automatically)
#' HULL(test_models$baseline$cormat, N = 500, gof = "CAF")
#'
#' # using ML with all available fit indices (CAF, CFI, and RMSEA)
#' HULL(test_models$baseline$cormat, N = 500, method = "ML")
#'
#' # using ULS with only RMSEA
#' HULL(test_models$baseline$cormat, N = 500, method = "ULS", gof = "RMSEA")
#'}
#'
#'\dontrun{
#' # using parallel processing (Note: plans can be adapted, see the future
#' # package for details)
#' future::plan(future::multisession)
#' HULL(test_models$baseline$cormat, N = 500, gof = "CAF")
#' }
HULL <- function(x, N = NA, n_fac_theor = NA,
                 method = c("PAF", "ULS", "ML"), gof = c("CAF", "CFI", "RMSEA"),
                 eigen_type = c("SMC", "PCA", "EFA"),
                 use = c("pairwise.complete.obs", "all.obs", "complete.obs",
                         "everything", "na.or.complete"),
                 cor_method = c("pearson", "spearman", "kendall"),
                 n_datasets = 1000, percent = 95,
                 decision_rule = c("means", "percentile", "crawford"),
                 n_factors = 1, ...) {
  # Perform hull method following Lorenzo-Seva, Timmerman, and Kiers (2011)

  if(!inherits(x, c("matrix", "data.frame"))){

    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(" 'x' is neither a matrix nor a dataframe. Either provide a correlation matrix or a dataframe or matrix with raw data.\n"))

  }

  method <- match.arg(method)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  gof <- match.arg(gof, several.ok = TRUE)
  eigen_type <- match.arg(eigen_type)
  checkmate::assert_count(n_fac_theor, na.ok = TRUE)
  checkmate::assert_count(N, na.ok = TRUE)
  decision_rule <- match.arg(decision_rule)
  checkmate::assert_count(n_factors)
  checkmate::assert_count(n_datasets)
  checkmate::assert_number(percent, lower = 0, upper = 100)

  if (ncol(x) < 6) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red( "Data has fewer than 6 indicators. Hull method needs at least 6.\n"))
  }

  if (method == "PAF" && !all(gof == "CAF")) {
    cli::cli_alert_info(cli::col_cyan('Only CAF can be used as gof if method "PAF" is used. Setting gof to "CAF"\n'))
    gof <- "CAF"
  }

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

  if (is.na(N)) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(' "N" is not specified but is needed for computation of some of the fit indices.\n'))
  }

  # Check if correlation matrix is invertable, if it is not, stop with message
  R_i <- try(solve(R), silent = TRUE)

  if (inherits(R_i, "try-error")) {
    stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(' Correlation matrix is singular, the HULL method cannot be exectued.\n'))
  }

  # Check if correlation matrix is positive definite
  if(any(eigen(R)$values <= 0)){

    R <- psych::cor.smooth(R)

  }

  m <- ncol(R)

  # 1) perform parallel analysis to find J as n_fac_theor + 1
  par_res <- PARALLEL(R, N = N, eigen_type = eigen_type, method = method,
                      n_datasets = n_datasets, percent = percent,
                      decision_rule = decision_rule, n_factors = n_factors,
                      ...)

  if(eigen_type == "SMC"){
    n_fac_PA <- par_res$n_fac_SMC
  } else if(eigen_type == "PCA"){
    n_fac_PA <- par_res$n_fac_PCA
  } else if(eigen_type == "EFA"){
    n_fac_PA <- par_res$n_fac_EFA
  }

  if (is.na(n_fac_PA)) {

    if (!is.na(n_fac_theor)) {
      J <- n_fac_theor + 1
    } else {
      J <- .det_max_factors(ncol(R))
    }

  } else {

    J <- max(c(n_fac_PA, n_fac_theor), na.rm = TRUE) + 1

  }

  if (J > .det_max_factors(ncol(R))) {
    J <- .det_max_factors(ncol(R))
    warning(crayon::yellow$bold("!"), crayon::yellow(' Setting maximum number of factors to',
                                                     J, 'to ensure overidentified models.\n'))
  }

  if (J < 3) {
    warning(crayon::yellow$bold("!"),
            crayon::yellow(" Suggested maximum number of factors was", J,
                           "but must be at least 3 for hull method to work.",
                           "Setting it to 3.\n"))
    J <- 3

  }

  # 2) perform factor analysis for the range of dimensions 1:J and compute f and
  #    df for every solution
  s <- matrix(0, ncol = 4, nrow = J + 1)
  s[, 1] <- 0:J

  # first for 0 factors
  if ("CAF" %in% gof) {
    s_CAF <- s
    colnames(s_CAF) <- c("nfactors", "CAF", "df", "st")
    s_CAF[1, 2] <- 1 - KMO(R)$KMO
    s_CAF[1, 3] <- (m**2 - m) / 2
  }

  if ("CFI" %in% gof) {
    s_CFI <- s
    colnames(s_CFI) <- c("nfactors", "CFI", "df", "st")
    s_CFI[1, 2] <- 0
    s_CFI[1, 3] <- (m**2 - m) / 2

    # for later use in loop
    chi_null <- sum(R[upper.tri(R)] ^ 2) * (N - 1)
    df_null <- (m**2 - m) / 2
    delta_hat_null <- max(0, chi_null - df_null)

  }

  if ("RMSEA" %in% gof) {
    s_RMSEA <- s
    colnames(s_RMSEA) <- c("nfactors", "RMSEA", "df", "st")
    Fm <- sum(R[upper.tri(R)] ^ 2)
    chi <- Fm * (N - 1)
    df <- (m**2 - m) / 2
    # compute 1 - RMSEA
    s_RMSEA[1, 2] <- 1 - sqrt(max(0, chi - df) / (df * N - 1))
    s_RMSEA[1, 3] <- (m**2 - m) / 2

  }

  # Calculate loadings with EFA function
  loadings <- suppressWarnings(future.apply::future_lapply(seq_len(J), EFA,
                                                           x = R,
                                                           method = method,
                                                           N = N, ...,
                                                           future.seed = FALSE))

  # then for 1 to J factors
  for (i in seq_len(J)) {
    if (method == "PAF") {
      # compute goodness of fit "f" as CAF (common part accounted for; Eq 3)
      # compute CAF
      s_CAF[i + 1, 2] <- loadings[[i]]$fit_indices$CAF
      # compute dfs (Eq 4 provides the number of free parameters; using dfs yields
      # th same numbers, as the difference in df equals the difference in free
      # parameters)
      s_CAF[i + 1, 3] <- loadings[[i]]$fit_indices$df
    } else {
      if ("CAF" %in% gof) {
        # compute goodness of fit "f" as CAF (common part accounted for; Eq 3)
        # compute CAF
        s_CAF[i + 1, 2] <- loadings[[i]]$fit_indices$CAF
        # compute dfs (Eq 4 provides the number of free parameters; using dfs yields
        # th same numbers, as the difference in df equals the difference in free
        # parameters)
        s_CAF[i + 1, 3] <- loadings[[i]]$fit_indices$df
      }

      if ("CFI" %in% gof) {
        # compute CFI
        s_CFI[i + 1, 2] <- loadings[[i]]$fit_indices$CFI
        # compute dfs (Eq 4 provides the number of free parameters; using dfs yields
        # th same numbers, as the difference in df equals the difference in free
        # parameters)
        s_CFI[i + 1, 3] <- loadings[[i]]$fit_indices$df

      }

      if ("RMSEA" %in% gof) {
        # compute 1 - RMSEA
        s_RMSEA[i + 1, 2] <- 1 - loadings[[i]]$fit_indices$RMSEA
        # compute dfs (Eq 4 provides the number of free parameters; using dfs yields
        # th same numbers, as the difference in df equals the difference in free
        # parameters)
        s_RMSEA[i + 1, 3] <- loadings[[i]]$fit_indices$df

      }

    }

  }

  out_CAF <- list(s_complete = NA, retain = NA)
  out_CFI <- list(s_complete = NA, retain = NA)
  out_RMSEA <- list(s_complete = NA, retain = NA)

  if("CAF" %in% gof) {
    out_CAF <- .hull_calc(s = s_CAF, J = J, gof_t = "CAF")
  }
  if("CFI" %in% gof) {
    out_CFI <- .hull_calc(s = s_CFI, J = J, gof_t = "CFI")
  }
  if("RMSEA" %in% gof) {
    out_RMSEA <- .hull_calc(s = s_RMSEA, J = J, gof_t = "RMSEA")
  }

  out <- list(
    n_fac_CAF = out_CAF$retain,
    n_fac_CFI = out_CFI$retain,
    n_fac_RMSEA = out_RMSEA$retain,
    solutions_CAF = out_CAF$s_complete,
    solutions_CFI = out_CFI$s_complete,
    solutions_RMSEA = out_RMSEA$s_complete,
    n_fac_max = J,
    settings = list(N = N,
                    method = method,
                    gof = gof,
                    n_fac_theor = n_fac_theor,
                    eigen_type = eigen_type,
                    use = use,
                    cor_method = cor_method)
  )

  class(out) <- "HULL"

  return(out)

}


.hull_calc <- function(s, J, gof_t){

  # 3) sort n solutions by their df values and denoted by s (already done)

  # 4) all solutions s are excluded for which a solution sj (j<i) exists such
  #    that fj > fi (eliminate solutions not on the boundary of the convex hull)

  s_complete <- s
  d_s <- diff(s[, 2])
  while (any(d_s < 0)) {
    s <- s[c(1, d_s) > 0, , drop = FALSE]
    if(nrow(s) == 1){
      break
    }
    d_s <- diff(s[, 2])
  }

  if (nrow(s) < 3) {
    warning(crayon::yellow$bold("!"),
            crayon::yellow(" Less than three solutions located on the hull have",
            "been identified when using", gof_t, "as goodness of fit index.",
                           "Proceeding by taking the value with the maximum",
                           gof_t, "as heuristic. You may want to consider",
                           "additional indices or methods as robustness check.\n"))

    # combine values
    for (row_i in 0:J) {

      if (row_i %in% s[,1]) {
        s_complete[row_i + 1, 4] <- s[s[,1] == row_i, 4]
      } else {
        s_complete[row_i + 1, 4] <- NA
      }

    }
    s_complete[, 4] <- rep(NA, nrow(s_complete))

    # 8) select solution with highest gof value
    retain <- s[which.max(s[, 2]), 1]


  } else {


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

  }


  out <- list(s_complete = s_complete,
              retain = unname(retain))

  return(out)

}
