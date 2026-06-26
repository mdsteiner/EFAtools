## Principal Axis Factoring
## Thin fitter: returns the raw PAF results; shared post-processing happens in
## .finalize_fit() / .estimate_model().
.PAF <- function(x, n_factors, max_iter = NA,
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

      h2_init <- .smc_start(R)

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

  # Negative eigenvalues make the communality estimates (square roots of the
  # eigenvalues) undefined. The kernel flags this rather than aborting so the
  # condition can be raised here as a classed R condition.
  if (isTRUE(L_list$neg_eigen)) {
    cli::cli_abort(
      c("Negative eigenvalues detected; cannot compute communality estimates.",
        "i" = 'Try {.code init_comm = "unity"} or {.code init_comm = "mac"}.'),
      class = "efa_paf_negative_eigen"
    )
  }

  h2 <- as.vector(L_list$h2)
  iter <- L_list$iter
  L <- L_list$L

  # save convergence status (0 = converged, 1 = not converged): the procedure has
  # not converged if the final change in the communality estimates still exceeds
  # the convergence criterion.
  if (isTRUE(L_list$delta > criterion)) {
    convergence <- 1
    cli::cli_warn(
      "Reached the maximum number of iterations without convergence; results may not be interpretable.",
      class = "efa_paf_nonconvergence"
    )
  } else {
    convergence <- 0
  }

  # store the settings used:
  settings <- list(
    max_iter = max_iter,
    init_comm = init_comm,
    criterion = criterion,
    criterion_type = criterion_type,
    abs_eigen = abs_eigen
  )

  # raw fit, finalized by .estimate_model()
  list(
    L = L,
    h2 = h2,
    Fm = NA_real_,
    iter = iter,
    convergence = convergence,
    orig_R = orig_R,
    R_final = L_list$R,
    h2_init = h2_init,
    init_eigen = init_eigen,
    settings = settings
  )

}
