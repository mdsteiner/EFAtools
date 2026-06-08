## Principal Axis Factoring
.PAF <- function(x, n_factors, N = NA, max_iter = NA,
                type = c("EFAtools", "psych", "SPSS", "none"),
                init_comm = NA, criterion = NA,
                criterion_type = NA, abs_eigen = NA) {

  # Get correlation matrix entered or created in EFA
  R <- x

  # Fill the type preset defaults and warn about any pinned tuning arguments.
  resolved <- .resolve_settings(
    type = type,
    user = list(init_comm = init_comm, criterion = criterion,
                criterion_type = criterion_type, max_iter = max_iter,
                abs_eigen = abs_eigen),
    preset = .efa_presets$PAF
  )
  init_comm <- resolved$init_comm
  criterion <- resolved$criterion
  criterion_type <- resolved$criterion_type
  max_iter <- resolved$max_iter
  abs_eigen <- resolved$abs_eigen

  # set initial communality estimates. This can be done in three ways:
  #  init_comm == "smc": uses the Squared Multiple Correlation
  #  init_comm == "unity": uses unity as inital estimates
  #  init_comm == "mac": uses Maximum Absolute Correlations
  if (init_comm == "smc") {

    # save original correlation matrix
    orig_R <- R

    # compute smcs
    if (type == "psych") {

      h2_init <- psych::smc(R)

    } else {

      h2_init <- 1 - 1 / diag(solve(R))

    }

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
  init_eigen <- eigen(R, symmetric = TRUE, only.values = TRUE)$values

  # define the number of factors m
  m <- n_factors

  # set communalities to init_communalities for comparison in first iteration
  h2 <- h2_init

  crit_type <- ifelse(criterion_type == "max_individual", 1, 2)

  # run the iterative PAF procedure using Rcpp
  L_list <- .paf_iter(h2 = h2, criterion = criterion, R = R, n_fac = m,
                     abs_eig = abs_eigen, crit_type = crit_type,
                     max_iter = max_iter)

  h2 <- as.vector(L_list$h2)
  R <- L_list$R
  iter <- L_list$iter
  L <- L_list$L

  # save convergence status (0 = converged, 1 = not converged (maximum number of
  # iterations reached))
  if (iter >= max_iter){
    convergence <- 1
  } else {
    convergence <- 0
  }

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

  if (!is.null(colnames(orig_R))) {
    # name the loading matrix so the variables can be identified
    rownames(L) <- colnames(orig_R)
  } else {
    varnames <- paste0("V", seq_len(ncol(orig_R)))
    colnames(orig_R) <- varnames
    rownames(orig_R) <- varnames
    rownames(L) <- varnames
  }

  colnames(L) <- paste0("F", seq_len(m))

  vars_accounted <- .compute_vars(L_unrot = L, L_rot = L)

  colnames(vars_accounted) <- colnames(L)

  fit_ind <- .gof(L, orig_R, N, "PAF", NA)

  # calculate model implied R
  model_implied_R <- L %*% t(L) + diag(1 - h2)

  # create the output object
  class(L) <- "LOADINGS"

  # store the settings used:

  settings <- list(
    max_iter = max_iter,
    init_comm = init_comm,
    criterion = criterion,
    criterion_type = criterion_type,
    abs_eigen = abs_eigen
  )

  # Name communalities
  names(h2) <- colnames(orig_R)
  names(h2_init) <- colnames(orig_R)

  # Create output
  output <- list(
    orig_R = orig_R,
    h2_init = h2_init,
    h2 = h2,
    orig_eigen = eigen(orig_R, symmetric = TRUE)$values,
    init_eigen = init_eigen,
    final_eigen = eigen(R, symmetric = TRUE)$values,
    iter = iter,
    convergence = convergence,
    unrot_loadings = L,
    vars_accounted = vars_accounted,
    fit_indices = fit_ind,
    model_implied_R = model_implied_R,
    residuals = orig_R - model_implied_R,
    settings = settings
  )

  output

}
