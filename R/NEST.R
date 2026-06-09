#' Next eigenvalue sufficiency test (NEST)
#'
#' NEST uses many synthetic datasets to generate reference eigenvalues against
#' which to compare the empirical eigenvalues. This is similar to parallel
#' analysis, but other than parallel analysis, NEST does not just rely on
#' synthetic eigenvalues based on an identity matrix as null model.
#' It was introduced by Achim (2017), see also Brandenburg and Papenberg (2024) and
#' Caron (2025) for further simulation studies including NEST.
#'
#' @param x data.frame or matrix. data.frame or matrix of raw data or matrix with
#' correlations.
#' @param N numeric. The number of observations. Only needed if x is a correlation
#'  matrix.
#' @param alpha numeric. The alpha level to use (i.e., 1-alpha percentile of eigenvalues is used for reference values).
#' @param use character. Passed to [stats::cor()] if raw
#'  data is given as input. Default is  `"pairwise.complete.obs"`.
#' @param cor_method character. Passed to [stats::cor()].
#'  Default is  `"pearson"`.
#' @param n_datasets numeric. The number of datasets to simulate. Default is 1000.
#' @param ... Additional arguments passed to [EFA()]. For example,
#' the extraction method can be changed here (default is "PAF"). PAF is more
#' robust, but it will take longer compared to the other estimation methods
#' available ("ML" and "ULS").
#'
#' @details NEST compares the first empirical eigenvalue against the first eigenvalues
#' of `n_dataset` synthetic datasets based on a null model  (i.e.,
#'  with uncorrelated variables; same as in parallel analysis, see [PARALLEL()]).
#'  The following eigenvalues are compared against synthetic datasets based on an EFA-model with one fewer factors
#'  than the position of the respective empirical eigenvalue. E.g, the second
#'  empirical eigenvalue is compared against synthetic data based on a one-factor
#'  model. The `alpha`-level defines against which percentile of the synthetic
#'  eigenvalue distribution to compare the empirical eigenvalues against, i.e., an
#'  alpha of .05 (the default) uses the 95th percentile as reference value.
#'
#'  For details on the method, including simulation studies, see Achim (2017),
#'  Brandenburg and Papenberg (2024), and Caron (2025).
#'
#'  The `NEST` function can also be called together with other factor
#'   retention criteria in the [N_FACTORS()] function.
#'
#'
#' @return A list of class NEST containing the following objects
#' \item{eigenvalues}{A vector containing the empirical eigenvalues of the entered data.}
#' \item{n_factors}{The number of factors to retain according to the NEST procedure.}
#' \item{references}{A vector containing the reference eigenvalues.}
#' \item{prob}{For the first n_factors + 1 empirical eigenvalues, the proportion <= the set of n_datasets  reference eigenvalues.}
#' \item{settings}{A list of control settings used in the print function.}
#'
#' @source Achim, A. (2017). Testing the number of required dimensions in exploratory factor analysis. The Quantitative Methods for Psychology, 13(1), 64–74. https://doi.org/10.20982/tqmp.13.1.p064
#' @source Brandenburg, N., & Papenberg, M. (2024). Reassessment of innovative methods to determine the number of factors: A simulation-based comparison of exploratory graph analysis and Next Eigenvalue Sufficiency Test. Psychological Methods, 29(1), 21–47. https://doi.org/10.1037/met0000527
#' @source Caron, P.-O. (2025). A Comparison of the Next Eigenvalue Sufficiency Test to Other Stopping Rules for the Number of Factors in Factor Analysis.
#' Educational and Psychological Measurement, Online-first publication. https://doi.org/10.1177/00131644241308528
#'
#' @export
#'
#' @examples
#'
#' # with correlation matrix
#' NEST(test_models$baseline$cormat, N = 500)
#'
#' # with raw data
#' NEST(GRiPS_raw)
NEST <- function(x, N = NA,
                 alpha = .05,
                 use = c("pairwise.complete.obs", "all.obs",
                         "complete.obs", "everything",
                         "na.or.complete"),
                 cor_method = c("pearson", "spearman", "kendall"),
                 n_datasets = 1000,
                 ...) {


  # Perform argument checks
  .assert_cor_input(x)

  checkmate::assert_count(N, na.ok = TRUE)
  checkmate::assert_number(alpha, lower = 0, upper = 1)
  use <- match.arg(use)
  cor_method <- match.arg(cor_method)
  checkmate::assert_count(n_datasets, na.ok = FALSE,
                          positive = TRUE)

  # Detect or compute the correlation matrix, check it, and smooth it if needed
  prep <- .prepare_cor_input(x, N = N, use = use, cor_method = cor_method,
                             N_policy = "required")
  R <- prep$R
  N <- prep$N

  emp_eigen <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  nvar <- ncol(R)

  # rule by Achim (2017)
  max_fac <- floor(.8 * nvar)

  references <- vector(mode = "double", length = max_fac)
  prob <- vector(mode = "double", length = max_fac)

  for (nf in 1:max_fac) {

    if (nf == 1) {

      # for the first factor, use identity matrix as reference data
      M <- diag(nvar)

    } else {


      # For further factors, use model with nf - 1 as reference data
      mm <- EFA(R, N = N, n_factors = nf - 1, ...)
      L <- mm$unrot_loadings
      u2 <- diag(sqrt(1-mm$h2))

      # bind loadings and uniquenesses together to create data according to
      # factor score rule X = F*L + e, which is faster than using the model
      # implied R to draw from a multivariate normal distribution.
      M <- cbind(L, u2)


    }

    # use C++ to simulate datasets
    ref_values <- .nest_sym(nf, N, t(M), n_datasets)

    references[nf] <- stats::quantile(ref_values, probs = 1 - alpha)
    prob[nf] <- mean(emp_eigen[nf] < ref_values)
    if (emp_eigen[nf] <= references[nf]) {
      break
    }


  }

  n_factors <- nf - 1
  references <- references[1:nf]
  prob <- prob[1:nf]

  out <- list(
    eigenvalues = emp_eigen,
    n_factors = n_factors,
    references = references,
    prob = prob,
    settings = list(
      alpha = alpha,
      N = N,
      n_datasets = n_datasets,
      use = use,
      cor_method = cor_method
    )
  )

  class(out) <- "NEST"

  out

}
