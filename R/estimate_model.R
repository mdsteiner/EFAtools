# Single estimation core: dispatch to a method fitter, then post-process once.
#
# `.PAF()`, `.ML()`, and `.ULS()` are thin fitters that return only their raw
# results (`L`, `h2`, objective `Fm`, `iter`, `convergence`, the original
# correlation matrix `orig_R`, and the matrix `R_final` whose eigenvalues are the
# final eigenvalues), plus any method-specific extras. `.finalize_fit()` performs
# the shared post-processing common to all methods, and `.estimate_model()`
# assembles the method-specific output object.

# Lower bound the ML and ULS optimisers impose on the uniquenesses. A solution
# whose uniqueness is pinned at this floor is an improper (boundary) solution:
# the unconstrained optimum would place it at or below zero. Shared by the ML/ULS
# fitters and the Heywood detector in .finalize_fit() so the value has one source.
.uniqueness_floor <- 0.005

# Squared multiple correlations (1 - 1/diag(R^-1)) used as starting communalities by the
# estimators. Shared by the PAF, ML, and ULS fitters so the start has one source.
.smc_start <- function(R) {
  1 - 1 / diag(solve(R))
}

# Shared post-processing for an unrotated solution. Reflects the loadings to a
# consistent sign, names them (with a V-fallback when the input is unnamed),
# computes the explained variances, fit indices, communalities, model-implied
# correlation matrix, residuals, and the original/final eigenvalues.
.finalize_fit <- function(fit, N, method, lean = FALSE) {

  L <- fit$L
  orig_R <- fit$orig_R
  h2 <- fit$h2
  n_factors <- ncol(L)

  # reverse the sign of loadings as done in the psych package and SPSS
  if (n_factors > 1) {
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

  # Bootstrap replicate path: .boot_se_ci() aggregates only each replicate's
  # unrotated loadings, fit indices, and residuals, so only those are computed.
  # The two eigendecompositions, variable naming, explained variances, and
  # Heywood detection are skipped (none are aggregated per replicate), and the
  # analytic RMSEA confidence bounds are not solved (.gof(ci = FALSE)). The sign
  # reflection above is kept so the loadings are identical to the full path,
  # leaving the downstream target alignment unchanged.
  if (lean) {
    model_implied_R <- L %*% t(L) + diag(1 - h2)
    return(list(
      unrot_loadings = L,
      fit_indices = .gof(L, orig_R, N, method, fit$Fm, ci = FALSE),
      residuals = orig_R - model_implied_R,
      convergence = fit$convergence
    ))
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

  colnames(L) <- paste0("F", seq_len(n_factors))

  vars_accounted <- .compute_vars(L_unrot = L, L_rot = L)
  colnames(vars_accounted) <- colnames(L)

  fit_ind <- .gof(L, orig_R, N, method, fit$Fm)

  # calculate model implied R
  model_implied_R <- L %*% t(L) + diag(1 - h2)

  # create the output object
  class(L) <- "LOADINGS"

  # Name communalities
  names(h2) <- colnames(orig_R)

  # Detect Heywood (improper) cases. Named integer vector of the affected
  # variables (empty if none); surfaced to the user by EFA() and shown in
  # summary(). Under PAF the communality can reach or exceed 1 directly. Under
  # ML/ULS the optimiser constrains the uniquenesses to [floor, 1], so an improper
  # solution instead shows up as a uniqueness pinned at the lower floor (the
  # boundary case); flag those too so detection is consistent across estimators.
  heywood_comm <- h2 >= 1
  heywood_boundary <- if (!is.null(fit$psi)) {
    fit$psi <= .uniqueness_floor + sqrt(.Machine$double.eps)
  } else {
    rep(FALSE, length(h2))
  }
  heywood <- which(heywood_comm | heywood_boundary)

  list(
    orig_R = orig_R,
    h2 = h2,
    orig_eigen = eigen(orig_R, symmetric = TRUE)$values,
    final_eigen = eigen(fit$R_final, symmetric = TRUE)$values,
    iter = fit$iter,
    convergence = fit$convergence,
    heywood = heywood,
    unrot_loadings = L,
    vars_accounted = vars_accounted,
    fit_indices = fit_ind,
    model_implied_R = model_implied_R,
    residuals = orig_R - model_implied_R
  )
}

# Dispatch to the requested fitter, run the shared post-processor, and assemble
# the method-specific output object (field set and order differ per method).
.estimate_model <- function(R, method, n_factors, N = NA,
                            type = "none", max_iter = NA, init_comm = NA,
                            criterion = NA, criterion_type = NA, abs_eigen = NA,
                            start_method = NA, lean = FALSE) {

  fit <- switch(
    method,
    PAF = .PAF(R, n_factors = n_factors, type = type, max_iter = max_iter,
               init_comm = init_comm, criterion = criterion,
               criterion_type = criterion_type, abs_eigen = abs_eigen),
    ML = .ML(R, n_factors = n_factors, start_method = start_method),
    ULS = .ULS(R, n_factors = n_factors)
  )

  common <- .finalize_fit(fit, N = N, method = method, lean = lean)

  # The bootstrap replicate path needs only the post-processed common fields;
  # skip the method-specific output assembly, which the aggregation never reads.
  if (lean) {
    return(common)
  }

  if (method == "PAF") {

    h2_init <- fit$h2_init
    names(h2_init) <- colnames(common$orig_R)

    output <- list(
      orig_R = common$orig_R,
      h2_init = h2_init,
      h2 = common$h2,
      orig_eigen = common$orig_eigen,
      init_eigen = fit$init_eigen,
      final_eigen = common$final_eigen,
      iter = common$iter,
      convergence = common$convergence,
      heywood = common$heywood,
      unrot_loadings = common$unrot_loadings,
      vars_accounted = common$vars_accounted,
      fit_indices = common$fit_indices,
      model_implied_R = common$model_implied_R,
      residuals = common$residuals,
      settings = fit$settings
    )

  } else if (method == "ML") {

    output <- c(common, list(settings = fit$settings))

  } else {

    output <- common

  }

  output
}
