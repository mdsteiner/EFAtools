## Principal Axis Factoring
.PAF <- function(x, n_factors, N = NA, max_iter = NA,
                type = c("EFAtools", "psych", "SPSS", "none"),
                init_comm = NA, criterion = NA,
                criterion_type = NA, abs_eigen = NA) {

  # Get correlation matrix entered or created in EFA
  R <- x

  if (type == "none") {

    # if type is none, throw an error if not
    # all the other necessary arguments are specified.

    if (is.na(init_comm) || is.na(criterion) || is.na(criterion_type) ||
        is.na(abs_eigen) || is.na(max_iter)) {
      stop(crayon::red$bold(cli::symbol$circle_cross), crayon::red(' One of "init_comm", "criterion", "criterion_type", "abs_eigen", "max_iter" was NA and no valid "type" was specified. Either use one of "EFAtools", "psych", or "SPSS" for type, or specify all other arguments\n'))
    }

  } else if (type == "EFAtools") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (is.na(init_comm)) {
      init_comm <- "smc"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and init_comm is specified. init_comm is used with value '", init_comm, "'. Results may differ from the specified type\n"))
    }

    if (is.na(criterion)) {
      criterion <- 1e-3
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and criterion is specified. criterion is used with value '", criterion, "'. Results may differ from the specified type\n"))
    }

    if (is.na(criterion_type)) {
      criterion_type <- "sum"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and criterion_type is specified. criterion_type is used with value '", criterion_type, "'. Results may differ from the specified type\n"))
    }

    if (is.na(max_iter)) {
      max_iter <- 300
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and max_iter is specified. max_iter is used with value '", max_iter, "'. Results may differ from the specified type\n"))
    }

    if (is.na(abs_eigen)) {
      abs_eigen <- TRUE
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and abs_eigen is specified. abs_eigen is used with value '", abs_eigen, "'. Results may differ from the specified type\n"))
    }


  } else if (type == "psych") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (is.na(init_comm)) {
      init_comm <- "smc"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and init_comm is specified. init_comm is used with value '", init_comm, "'. Results may differ from the specified type\n"))
    }

    if (is.na(criterion)) {
      criterion <- .001
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and criterion is specified. criterion is used with value '", criterion, "'. Results may differ from the specified type\n"))
    }

    if (is.na(criterion_type)) {
      criterion_type <- "sum"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and criterion_type is specified. criterion_type is used with value '", criterion_type, "'. Results may differ from the specified type\n"))
    }

    if (is.na(max_iter)) {
      max_iter <- 50
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and max_iter is specified. max_iter is used with value '", max_iter, "'. Results may differ from the specified type\n"))
    }

    if (is.na(abs_eigen)) {
      abs_eigen <- FALSE
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and abs_eigen is specified. abs_eigen is used with value '", abs_eigen, "'. Results may differ from the specified type\n"))
    }


  } else if (type == "SPSS") {

    # if not specified, set PAF properties. If specified, throw warning that
    # results may not exactly match the specified type

    if (is.na(init_comm)) {
      init_comm <- "smc"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and init_comm is specified. init_comm is used with value '", init_comm, "'. Results may differ from the specified type\n"))
    }

    if (is.na(criterion)) {
      criterion <- .001
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and criterion is specified. criterion is used with value '", criterion, "'. Results may differ from the specified type\n"))
    }

    if (is.na(criterion_type)) {
      criterion_type <- "max_individual"
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and criterion_type is specified. criterion_type is used with value '", criterion_type, "'. Results may differ from the specified type\n"))
    }

    if (is.na(max_iter)) {
      max_iter <- 25
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and max_iter is specified. max_iter is used with value '", max_iter, "'. Results may differ from the specified type\n"))
    }

    if (is.na(abs_eigen)) {
      abs_eigen <- TRUE
    } else {
      warning(crayon::yellow$bold("!"), crayon::yellow(" Type and abs_eigen is specified. abs_eigen is used with value '", abs_eigen, "'. Results may differ from the specified type\n"))
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
    settings = settings
  )

  output

}
