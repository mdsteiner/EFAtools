#' Exploratory factor analysis on multiple data imputations
#'
#' @author Andreas Soteriades, Markus Steiner
#'
#' @description
#' Fits \code{\link{EFA}} to each of several imputed datasets, aligns the
#' factor solutions to a common factor space, and pools the resulting estimates
#' and selected fit quantities across imputations.
#'
#' @details
#' The function first fits the same \code{\link{EFA}} model to each imputed
#' dataset. Unrotated loading matrices can optionally be put into a common
#' signed/permuted factor order before averaging. Rotated loading matrices are
#' aligned with either a consensus Procrustes target or with the first
#' imputation's rotated solution as a fixed target. For oblique solutions,
#' factor intercorrelations are transformed/aligned together with the loading
#' matrices so that the factor model remains internally consistent.
#'
#' Point estimates are pooled by arithmetic averaging after alignment. For
#' oblique rotations, the returned structure matrix is computed from the pooled
#' aligned pattern matrix and the pooled factor correlation matrix,
#' \eqn{Structure = Lambda Phi}. Communalities are always computed as the
#' diagonal of the common-factor reproduced correlation matrix,
#' \eqn{diag(Lambda Phi Lambda')} for oblique rotations and
#' \eqn{diag(Lambda Lambda')} otherwise.
#'
#' Residuals are not averaged from the per-imputation residual matrices. Instead,
#' the observed correlation matrices are averaged across imputations and residuals
#' are calculated from the pooled observed correlation matrix minus the
#' model-implied correlation matrix of the pooled solution. Consequently,
#' residual-based fit indices such as RMSR/SRMR are based on pooled residuals.
#'
#' Fit indices based on the model chi-square are not arithmetic means of the
#' per-imputation fit indices. If possible, chi-square-type fit is pooled with
#' the D2 rule. The asymptotic chi-square approximation to D2 is then used for
#' RMSEA and CFI. AIC and BIC, if returned, are chi-square-derived descriptive
#' quantities based on this D2 approximation and should not be interpreted as
#' likelihood-based MI information criteria. D3/D4 pooling is not implemented
#' here because the current \code{EFA} objects do not expose the likelihood
#' quantities needed for those methods.
#'
#' If each component \code{EFA} call was run with \code{se = "np-boot"} and
#' returned \code{boot.arrays}, pooled bootstrap SEs and Wald-type MI confidence
#' intervals are computed for loadings, communalities, residuals, and, when
#' applicable, factor correlations and structure coefficients. Importantly, the
#' rotated bootstrap loading matrices stored by the component \code{EFA} calls
#' are not reused directly, because they were aligned to imputation-specific
#' targets. Instead, the unrotated bootstrap loading matrices are re-aligned to
#' the final MI target before estimating within-imputation bootstrap covariance
#' matrices. Rubin-type MI pooling is then applied with
#' \eqn{T = Ubar + (1 + 1 / m) B}. Confidence intervals for loadings, Phi,
#' communalities, residuals, and structure coefficients are Wald-type MI
#' intervals; \code{boot.CI$fit_indices_pooled_algorithm}, when available, is a
#' percentile-style summary obtained by re-running the pooled-fit algorithm over
#' matched bootstrap replicate indices.
#'
#' @param data_list A list of length \eqn{m}, where \eqn{m} is the number of
#' imputations. Each list element is a data frame or matrix of raw data, or a
#' correlation matrix. See argument \code{x} in \code{\link{EFA}}.
#' @param p Numeric. One minus the confidence level used for pooled Wald-type
#' bootstrap/MI confidence intervals when bootstrap arrays are available.
#' @param target_method Character. \code{"consensus"} aligns all solutions to an
#' iteratively updated consensus target. \code{"first_target"} aligns all
#' solutions to the first imputation's rotated solution.
#' @param align_unrotated Character. How to align unrotated loadings before
#' pooling. \code{"signed_tucker_congruence"} preserves the unrotated axes up
#' to factor reordering and sign changes using Tucker congruence.
#' \code{"procrustes"} aligns the unrotated matrices to the first imputation by
#' orthogonal Procrustes rotation. \code{"none"} averages unrotated loadings as
#' returned by \code{EFA}.
#' @param fit_pool_method Character. Currently only \code{"D2"} is implemented
#' for chi-square-type fit. If no chi-square is available, only residual-based
#' fit and descriptive quantities are returned.
#' @param consensus_args List of additional arguments passed to
#' \code{CONSENSUS_PROCRUSTES}.
#' @param procrustes_args List of additional arguments passed to \code{PROCRUSTES}
#' for fixed-target alignment.
#' @param rmsea_ci_level Numeric. Confidence level for the RMSEA CI.
#' @param rmsr_upper Logical. If \code{TRUE}, compute RMSR from the unique
#' off-diagonal residual correlations. If \code{FALSE}, use the full off-diagonal
#' matrix.
#' @param ... Additional arguments passed to \code{\link{EFA}}.
#'
#' @return A list of class \code{"EFA_POOLED"} containing pooled estimates,
#' residuals, fit indices, the individual fits, and MI diagnostics.
#'
#' @export
EFA_POOLED <- function(data_list,
                       p = 0.05,
                       target_method = c("consensus", "first_target"),
                       align_unrotated = c("signed_tucker_congruence", "none", "procrustes"),
                       fit_pool_method = c("D2"),
                       consensus_args = list(),
                       procrustes_args = list(),
                       rmsea_ci_level = .90,
                       rmsr_upper = TRUE,
                       ...) {

  efa_args <- list(...)

  checkmate::assert_list(data_list, min.len = 2, null.ok = FALSE)
  lapply(data_list, checkmate::assert_multi_class, c("matrix", "data.frame"))
  checkmate::assert_number(p, na.ok = FALSE, lower = 0, upper = 1)
  checkmate::assert_list(consensus_args, null.ok = FALSE)
  checkmate::assert_list(procrustes_args, null.ok = FALSE)
  checkmate::assert_number(rmsea_ci_level, na.ok = FALSE, lower = 0, upper = 1)
  checkmate::assert_flag(rmsr_upper)

  target_method <- match.arg(target_method)
  align_unrotated <- match.arg(align_unrotated)
  fit_pool_method <- match.arg(fit_pool_method)

  m_imp <- length(data_list)

  ## -------------------------------------------------------------------------
  ## Fit EFA to each imputed dataset
  ## -------------------------------------------------------------------------

  fits <- lapply(data_list, function(data_list_subset) {
    do.call(EFA, c(list(x = data_list_subset), efa_args))
  })

  .efa_pooled_check_fits(fits)

  settings <- fits[[1]]$settings
  method <- settings$method
  rotation <- settings$rotation
  var_names <- rownames(fits[[1]]$orig_R)
  if (is.null(var_names)) {
    var_names <- colnames(data_list[[1]])
  }

  if (is.null(rotation) || identical(rotation, "none")) {
    rotation_type <- "none"
  } else if (rotation %in% c("varimax", "quartimax", "equamax", "bentlerT",
                             "geominT", "bifactorT")) {
    rotation_type <- "orthogonal"
  } else if (rotation %in% c("promax", "oblimin", "quartimin", "simplimax",
                             "bentlerQ", "geominQ", "bifactorQ")) {
    rotation_type <- "oblique"
  } else {
    stop("Unknown rotation type: ", rotation, call. = FALSE)
  }

  ## -------------------------------------------------------------------------
  ## Extract and align unrotated loadings
  ## -------------------------------------------------------------------------

  unrot_loadings <- lapply(.extract_list_object(fits, "unrot_loadings"),
                           .change_class, "matrix")

  unrot_loadings_aligned <- .efa_pooled_align_unrotated_list(
    unrot_loadings = unrot_loadings,
    align_unrotated = align_unrotated
  )

  mean_unrot_loadings <- .average_matrices(unrot_loadings_aligned)

  ## -------------------------------------------------------------------------
  ## Align rotated loadings and Phi, if a rotation was requested
  ## -------------------------------------------------------------------------

  phis <- NULL
  structure_loadings <- NULL

  if (rotation_type != "none") {
    rot_loadings_initial <- lapply(.extract_list_object(fits, "rot_loadings"),
                                   .change_class, "matrix")

    if (target_method == "first_target") {
      target_rotations <- vector("list", m_imp)
      rot_loadings <- vector("list", m_imp)
      phis <- vector("list", m_imp)

      rot_loadings[[1]] <- rot_loadings_initial[[1]]
      if (rotation_type == "oblique") {
        phis[[1]] <- fits[[1]]$Phi
      }

      for (d in 2:m_imp) {
        target_rotations[[d]] <- do.call(
          PROCRUSTES,
          c(list(A = unrot_loadings[[d]],
                 Target = rot_loadings_initial[[1]],
                 rotation = rotation_type),
            procrustes_args)
        )
        rot_loadings[[d]] <- target_rotations[[d]]$loadings
        if (rotation_type == "oblique") {
          phis[[d]] <- target_rotations[[d]]$Phi
        }
      }

      final_target <- rot_loadings_initial[[1]]
      alignment <- list(method = "first_target",
                        target = final_target,
                        target_rotations = target_rotations)

    } else if (target_method == "consensus") {
      consensus <- do.call(
        CONSENSUS_PROCRUSTES,
        c(list(unrotated_list = unrot_loadings,
               init_targets = rot_loadings_initial,
               rotation = rotation_type),
          consensus_args)
      )

      rot_loadings <- consensus$aligned_loadings
      phis <- consensus$aligned_phi
      final_target <- consensus$target
      alignment <- consensus
    }

    mean_rot_loadings <- .average_matrices(rot_loadings)
    rownames(mean_rot_loadings) <- var_names

    if (rotation_type == "oblique") {
      mean_phis <- .average_matrices(phis)
      structure_loadings <- Map(function(L, Phi) L %*% Phi, rot_loadings, phis)
      # Keep Structure parallel to the returned pooled pattern matrix and Phi:
      # it is the plug-in structure of the pooled solution, not the arithmetic
      # mean of the imputation-specific structure matrices.
      mean_structure_loadings <- mean_rot_loadings %*% mean_phis
      dimnames(mean_structure_loadings) <- dimnames(mean_rot_loadings)
    } else {
      mean_phis <- NULL
      mean_structure_loadings <- NULL
    }

  } else {
    rot_loadings <- NULL
    mean_rot_loadings <- NULL
    mean_phis <- NULL
    mean_structure_loadings <- NULL
    final_target <- NULL
    alignment <- NULL
  }

  ## -------------------------------------------------------------------------
  ## Pooled observed and model-implied correlation matrices; residuals and RMSR
  ## -------------------------------------------------------------------------

  orig_R_list <- .extract_list_object(fits, "orig_R")

  pooled_orig_R <- .average_matrices(orig_R_list)

  if (rotation_type == "oblique") {
    common_R <- mean_rot_loadings %*% mean_phis %*% t(mean_rot_loadings)
  } else if (rotation_type == "orthogonal") {
    common_R <- mean_rot_loadings %*% t(mean_rot_loadings)
  } else {
    common_R <- mean_unrot_loadings %*% t(mean_unrot_loadings)
  }

  h2 <- diag(common_R)
  names(h2) <- rownames(common_R)

  model_implied_R <- common_R
  diag(model_implied_R) <- 1
  dimnames(model_implied_R) <- list(var_names, var_names)

  residuals <- pooled_orig_R - model_implied_R
  diag(residuals) <- 0
  dimnames(residuals) <- list(var_names, var_names)

  RMSR <- .rmsr(residuals, upper = rmsr_upper)

  ## -------------------------------------------------------------------------
  ## Variance-accounted tables
  ## -------------------------------------------------------------------------

  mean_vars_accounted <- .compute_vars(
    L_unrot = mean_unrot_loadings,
    L_rot = mean_unrot_loadings,
    Phi = NULL
  )

  if (rotation_type == "oblique") {
    mean_vars_accounted_rot <- .compute_vars(
      L_unrot = mean_unrot_loadings,
      L_rot = mean_rot_loadings,
      Phi = mean_phis
    )
  } else if (rotation_type == "orthogonal") {
    mean_vars_accounted_rot <- .compute_vars(
      L_unrot = mean_unrot_loadings,
      L_rot = mean_rot_loadings,
      Phi = NULL
    )
  } else {
    mean_vars_accounted_rot <- NULL
  }

  ## -------------------------------------------------------------------------
  ## Fit indices: residual-based from pooled residuals, chi-square-based via D2
  ## -------------------------------------------------------------------------

  Ns <- .efa_pooled_get_Ns(data_list, fits, efa_args)
  Ns_ok <- Ns[is.finite(Ns)]
  if (length(unique(Ns_ok)) > 1L) {
    warning("The imputed datasets appear to have different N. Fit indices use the mean N across imputations.",
            call. = FALSE)
  }
  N_pool <- mean(Ns, na.rm = TRUE)

  fit_indices <- .efa_pooled_fit_indices(
    fits = fits,
    pooled_R = pooled_orig_R,
    residuals = residuals,
    RMSR = RMSR,
    N = N_pool,
    method = method,
    pool_method = fit_pool_method,
    rmsea_ci_level = rmsea_ci_level
  )

  ## -------------------------------------------------------------------------
  ## Bootstrap SEs and CIs for pooled MI estimates, if available
  ## -------------------------------------------------------------------------

  boot_pooled <- .efa_pooled_bootstrap_pool(
    fits = fits,
    orig_R_list = orig_R_list,
    unrot_loadings_aligned = unrot_loadings_aligned,
    mean_unrot_loadings = mean_unrot_loadings,
    rot_loadings = rot_loadings,
    phis = phis,
    structure_loadings = structure_loadings,
    mean_structure_loadings = mean_structure_loadings,
    final_target = final_target,
    rotation_type = rotation_type,
    align_unrotated = align_unrotated,
    procrustes_args = procrustes_args,
    h2 = h2,
    residuals = residuals,
    fit_indices = fit_indices,
    pooled_orig_R = pooled_orig_R,
    N = N_pool,
    method = method,
    pool_method = fit_pool_method,
    rmsea_ci_level = rmsea_ci_level,
    alpha = p,
    rmsr_upper = rmsr_upper
  )

  ## -------------------------------------------------------------------------
  ## Return object
  ## -------------------------------------------------------------------------

  mean_unrot_loadings <- .change_class(mean_unrot_loadings, "LOADINGS")
  if (!is.null(mean_rot_loadings)) {
    mean_rot_loadings <- .change_class(mean_rot_loadings, "LOADINGS")
  }
  if (!is.null(mean_structure_loadings)) {
    mean_structure_loadings <- .change_class(mean_structure_loadings, "LOADINGS")
  }

  results <- list(
    h2 = h2,
    unrot_loadings = mean_unrot_loadings,
    vars_accounted = mean_vars_accounted,
    fit_indices = fit_indices,
    model_implied_R = model_implied_R,
    residuals = residuals,
    orig_R = pooled_orig_R,
    settings = c(settings, list(
      pooled = TRUE,
      n_imputations = m_imp,
      target_method = target_method,
      align_unrotated = align_unrotated,
      fit_pool_method = fit_pool_method
    )),
    fits = fits,
    alignment = alignment,
    mi_diagnostics = fit_indices$mi_diagnostics
  )

  if (rotation_type != "none") {
    results$rot_loadings <- mean_rot_loadings
    results$vars_accounted_rot <- mean_vars_accounted_rot
  }

  if (rotation_type == "oblique") {
    results$Phi <- mean_phis
    results$Structure <- mean_structure_loadings
  }

  if (!is.null(boot_pooled)) {
    results$boot.SE <- boot_pooled$SE
    results$boot.CI <- boot_pooled$CI
    results$boot.arrays <- boot_pooled$arrays
    results$boot.MI <- boot_pooled$MI
    results$standardized_residuals <- results$residuals / boot_pooled$SE$residuals
  }

  class(results) <- c("EFA_POOLED", "EFA")
  results
}

## =============================================================================
## Internal helpers for EFA_POOLED
## =============================================================================

.efa_pooled_get_Ns <- function(data_list, fits, efa_args) {
  # Recover the N used in each EFA fit. Correlation-matrix input may not carry N,
  # so return NA unless N was supplied to EFA() or stored in settings.
  vapply(seq_along(data_list), function(d) {
    if (!is.null(fits[[d]]$settings$N) && !is.na(fits[[d]]$settings$N)) {
      return(as.numeric(fits[[d]]$settings$N))
    }
    if (!is.null(efa_args$N) && !is.na(efa_args$N)) {
      return(as.numeric(efa_args$N))
    }
    if (!.is_cormat(data_list[[d]])) {
      return(nrow(data_list[[d]]))
    }
    NA_real_
  }, numeric(1))
}

.efa_pooled_check_fits <- function(fits) {
  # Fail early if the fitted EFA objects are not conformable. Pooling only makes
  # sense when all imputations estimate the same model on the same variables.
  if (length(fits) < 2L) {
    stop("At least two EFA fits are required for MI pooling.", call. = FALSE)
  }

  dims <- vapply(fits, function(x) {
    paste(dim(as.matrix(x$unrot_loadings)), collapse = "x")
  }, character(1))
  if (length(unique(dims)) != 1L) {
    stop("All unrotated loading matrices must have the same dimensions.", call. = FALSE)
  }

  var_names <- lapply(fits, function(x) rownames(as.matrix(x$orig_R)))
  if (!all(vapply(var_names[-1], identical, logical(1), var_names[[1]]))) {
    stop("All imputations must contain the same variables in the same order.", call. = FALSE)
  }

  setting_value <- function(x, name) {
    val <- x$settings[[name]]
    if (is.null(val)) NA_character_ else as.character(val)
  }

  for (nm in c("method", "rotation", "n_factors")) {
    vals <- vapply(fits, setting_value, character(1), name = nm)
    vals <- vals[!is.na(vals)]
    if (length(unique(vals)) > 1L) {
      stop("All imputations must use the same ", nm, ".", call. = FALSE)
    }
  }

  invisible(TRUE)
}

.efa_pooled_align_unrotated_list <- function(unrot_loadings,
                                             align_unrotated = c("signed_tucker_congruence", "none", "procrustes")) {
  # Unrotated factor axes are arbitrary up to ordering and signs. This helper
  # puts them into a common orientation before simple arithmetic averaging.
  # Unlike the rotated solution below, this step should not seek simple
  # structure; it only removes indeterminacy in the unrotated axes.
  align_unrotated <- match.arg(align_unrotated)

  if (align_unrotated == "none") {
    return(unrot_loadings)
  }

  target <- as.matrix(unrot_loadings[[1]])
  out <- vector("list", length(unrot_loadings))
  out[[1]] <- target

  for (d in seq_along(unrot_loadings)[-1]) {
    Ld <- as.matrix(unrot_loadings[[d]])

    if (align_unrotated == "signed_tucker_congruence") {
      # Preserve the axes up to a one-to-one permutation/sign change based on
      # Tucker congruence. .align_solution() implements the assignment and sign
      # correction without applying a continuous rotation.
      out[[d]] <- .align_solution(L_target = target, L = Ld)$loadings
    } else if (align_unrotated == "procrustes") {
      # Optional continuous orthogonal alignment for users who prefer a least-
      # squares Procrustes match for the unrotated matrices.
      out[[d]] <- .change_class(PROCRUSTES(
        A = Ld,
        Target = target,
        rotation = "orthogonal"
      )$loadings, "matrix")
    }
  }

  out
}

.efa_pooled_rmsea_ci <- function(chi, df, N, level = .90) {
  # RMSEA CI from the noncentral chi-square inversion used by EFAtools.
  # The chi-square supplied here should already be the pooled test statistic.
  if (is.na(chi) || is.na(df) || is.na(N) || df <= 0 || N <= 1) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  alpha <- 1 - level
  lower_goal <- 1 - alpha / 2
  upper_goal <- alpha / 2

  pfun <- function(lambda, goal) {
    goal - stats::pchisq(chi, df = df, ncp = lambda)
  }

  lambda_l <- 0
  lambda_u <- 0

  if (stats::pchisq(chi, df = df, ncp = 0) >= lower_goal) {
    lambda_l <- tryCatch(
      stats::uniroot(pfun, interval = c(1e-10, 10000), goal = lower_goal,
                     extendInt = "upX", maxiter = 100L)$root,
      error = function(e) NA_real_
    )
  }

  if (stats::pchisq(chi, df = df, ncp = 0) >= upper_goal) {
    lambda_u <- tryCatch(
      stats::uniroot(pfun, interval = c(1e-10, 10000), goal = upper_goal,
                     extendInt = "upX", maxiter = 100L)$root,
      error = function(e) NA_real_
    )
  }

  denom <- df * (N - 1)
  c(lower = sqrt(lambda_l / denom), upper = sqrt(lambda_u / denom))
}

.efa_pooled_D2 <- function(chis, df) {
  # D2 pools complete-data chi-square statistics across imputations. The
  # returned chi value is the asymptotic chi-square approximation df * F_D2,
  # which is used downstream for RMSEA/CFI.
  chis <- chis[is.finite(chis)]
  M <- length(chis)

  if (M < 2L || !is.finite(df) || df <= 0) {
    return(NULL)
  }

  Tbar <- mean(chis)

  if (Tbar <= 0) {
    return(list(
      F = 0,
      chi = 0,
      df1 = df,
      df2 = Inf,
      p = 1,
      ARIV = 0,
      FMI = 0,
      Tbar = Tbar,
      M = M
    ))
  }

  ARIV <- (1 + 1 / M) * sum((sqrt(chis) - sqrt(Tbar))^2) / (M - 1)

  if (ARIV <= .Machine$double.eps) {
    F_D2 <- Tbar / df
    df2 <- Inf
    p_val <- stats::pchisq(Tbar, df = df, lower.tail = FALSE)
  } else {
    F_D2 <- (Tbar / df - ((M + 1) / (M - 1)) * ARIV) / (1 + ARIV)
    df2 <- df^(-3 / M) * (M - 1) * (1 + 1 / ARIV)^2
    p_val <- if (F_D2 < 0) 1 else stats::pf(F_D2, df1 = df, df2 = df2, lower.tail = FALSE)
  }

  chi_D2 <- max(0, df * F_D2)
  FMI <- ARIV / (1 + ARIV)

  list(
    F = F_D2,
    chi = chi_D2,
    df1 = df,
    df2 = df2,
    p = p_val,
    ARIV = ARIV,
    FMI = FMI,
    Tbar = Tbar,
    M = M
  )
}

.efa_pooled_fit_indices <- function(fits,
                                    pooled_R,
                                    residuals,
                                    RMSR,
                                    N,
                                    method,
                                    pool_method = "D2",
                                    rmsea_ci_level = .90) {

  fit_list <- .extract_list_object(fits, "fit_indices")
  fit_list <- fit_list[!vapply(fit_list, is.null, logical(1))]

  p_vars <- nrow(pooled_R)
  df <- NA_real_
  chis <- numeric(0)

  if (length(fit_list) > 0L) {
    dfs <- vapply(fit_list, function(x) if (!is.null(x$df)) x$df else NA_real_, numeric(1))
    finite_dfs <- unique(stats::na.omit(dfs))
    if (length(finite_dfs) > 1L) {
      stop("Cannot D2-pool chi-square fit because the imputation-specific dfs differ.",
           call. = FALSE)
    }
    df <- finite_dfs[1]
    chis <- vapply(fit_list, function(x) if (!is.null(x$chi)) x$chi else NA_real_, numeric(1))
  }

  D2 <- NULL
  if (identical(pool_method, "D2") && length(chis) > 0L && is.finite(df)) {
    D2 <- .efa_pooled_D2(chis, df)
  }

  if (!is.null(D2)) {
    chi <- D2$chi
    p_chi <- D2$p
    Fm <- if (is.finite(N) && N > 1) chi / (N - 1) else NA_real_
  } else {
    chi <- NA_real_
    p_chi <- NA_real_
    Fm <- NA_real_
  }

  df_null <- p_vars * (p_vars - 1) / 2
  chi_null <- if (is.finite(N) && N > 1) {
    sum(pooled_R[upper.tri(pooled_R)]^2) * (N - 1)
  } else {
    NA_real_
  }
  p_null <- if (is.finite(chi_null)) stats::pchisq(chi_null, df_null, lower.tail = FALSE) else NA_real_

  if (is.finite(chi) && is.finite(chi_null) && is.finite(df) && df >= 0) {
    delta_null <- chi_null - df_null
    delta_model <- chi - df
    CFI <- (delta_null - delta_model) / delta_null
    CFI <- max(0, min(1, CFI))
  } else {
    CFI <- NA_real_
  }

  if (is.finite(chi) && is.finite(df) && is.finite(N) && df > 0 && N > 1) {
    RMSEA <- sqrt(max((chi - df) / (df * (N - 1)), 0))
    rmsea_ci <- .efa_pooled_rmsea_ci(chi, df, N, level = rmsea_ci_level)
  } else {
    RMSEA <- NA_real_
    rmsea_ci <- c(lower = NA_real_, upper = NA_real_)
  }

  # These mirror the chi-square-derived quantities used elsewhere in EFAtools.
  # Under MI/D2 pooling they are descriptive only, not likelihood-based MI AIC/BIC.
  AIC <- if (is.finite(chi) && is.finite(df)) chi - 2 * df else NA_real_
  BIC <- if (is.finite(chi) && is.finite(df) && is.finite(N)) chi - log(N) * df else NA_real_

  ## CAF in EFAtools is 1 - KMO(delta_hat), with diagonal set to 1.
  delta_hat <- residuals
  diag(delta_hat) <- 1
  CAF <- tryCatch(1 - KMO(delta_hat)$KMO,
                  error = function(e) NA_real_)

  out <- list(
    chi = chi,
    df = df,
    p_chi = p_chi,
    CAF = CAF,
    CFI = CFI,
    RMSEA = RMSEA,
    RMSEA_LB = rmsea_ci[["lower"]],
    RMSEA_UB = rmsea_ci[["upper"]],
    AIC = AIC,
    BIC = BIC,
    Fm = Fm,
    RMSR = RMSR,
    SRMR = RMSR,
    chi_null = chi_null,
    df_null = df_null,
    p_null = p_null,
    pool_method = pool_method,
    mi_diagnostics = list(
      D2_F = if (!is.null(D2)) D2$F else NA_real_,
      D2_df1 = if (!is.null(D2)) D2$df1 else NA_real_,
      D2_df2 = if (!is.null(D2)) D2$df2 else NA_real_,
      D2_chi_asymptotic = if (!is.null(D2)) D2$chi else NA_real_,
      chi_bar_naive = if (!is.null(D2)) D2$Tbar else if (length(chis) > 0L) mean(chis, na.rm = TRUE) else NA_real_,
      ARIV = if (!is.null(D2)) D2$ARIV else NA_real_,
      FMI = if (!is.null(D2)) D2$FMI else NA_real_,
      m = if (!is.null(D2)) D2$M else length(fits)
    )
  )

  out
}

## -----------------------------------------------------------------------------
## Bootstrap pooling helpers
## -----------------------------------------------------------------------------

.efa_pooled_has_boot_arrays <- function(fits) {
  # Bootstrap MI pooling requires the raw unrotated bootstrap loading arrays.
  # Scalar bootstrap SEs/CIs are not enough because the replicates must be
  # re-expressed in the final MI target space.
  vapply(fits, function(x) {
    !is.null(x$boot.arrays) &&
      !is.null(x$boot.arrays$unrot_loadings)
  }, logical(1))
}

.efa_pooled_vec <- function(M) {
  # Column-major vectorization so matrix estimates can be Rubin-pooled.
  as.vector(as.matrix(M))
}

.efa_pooled_vech <- function(M) {
  # Vectorize the lower triangle of a symmetric Phi matrix without duplicates.
  M <- as.matrix(M)
  M[lower.tri(M, diag = TRUE)]
}

.efa_pooled_unvech <- function(v, k) {
  # Reconstruct a symmetric matrix from vech() output.
  M <- matrix(NA_real_, k, k)
  M[lower.tri(M, diag = TRUE)] <- v
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  M
}

.efa_pooled_communalities <- function(Lambda, Phi = NULL) {
  # Communalities are the diagonal of the common-factor reproduced matrix.
  Lambda <- as.matrix(Lambda)
  common_R <- if (is.null(Phi)) {
    Lambda %*% t(Lambda)
  } else {
    Lambda %*% as.matrix(Phi) %*% t(Lambda)
  }
  diag(common_R)
}

.efa_pooled_model_implied <- function(Lambda, Phi = NULL) {
  # Correlation-model implied matrix: common-factor part plus uniquenesses,
  # implemented by setting the diagonal of Lambda Phi Lambda' to one.
  Lambda <- as.matrix(Lambda)
  common_R <- if (is.null(Phi)) {
    Lambda %*% t(Lambda)
  } else {
    Lambda %*% as.matrix(Phi) %*% t(Lambda)
  }
  implied <- common_R
  diag(implied) <- 1
  implied
}

.efa_pooled_residual_from_solution <- function(R, Lambda, Phi = NULL) {
  # Residuals for one imputation in the same target space as its aligned solution.
  E <- as.matrix(R) - .efa_pooled_model_implied(Lambda, Phi)
  diag(E) <- 0
  E
}

.efa_pooled_make_ci <- function(est, se, df, alpha) {
  # Wald-type confidence intervals using Rubin degrees of freedom when finite.
  crit <- ifelse(
    is.finite(df),
    stats::qt(1 - alpha / 2, df = df),
    stats::qnorm(1 - alpha / 2)
  )
  list(
    lower = est - crit * se,
    upper = est + crit * se
  )
}

.efa_pooled_rubin_pool <- function(q_list, boot_mat_list, alpha = 0.05,
                                   est_override = NULL) {
  # q_list contains one aligned point-estimate vector per imputation.
  # boot_mat_list contains aligned bootstrap replicates for estimating U_d.
  m <- length(q_list)
  q_mat <- do.call(rbind, q_list)
  q_bar <- colMeans(q_mat, na.rm = TRUE)

  U_list <- lapply(boot_mat_list, function(Bmat) {
    Bmat <- as.matrix(Bmat)
    if (nrow(Bmat) < 2L) {
      return(matrix(NA_real_, ncol(Bmat), ncol(Bmat)))
    }
    stats::cov(Bmat, use = "pairwise.complete.obs")
  })

  U_bar <- Reduce("+", U_list) / m
  B_mi <- if (m > 1L) {
    stats::cov(q_mat, use = "pairwise.complete.obs")
  } else {
    matrix(0, ncol(q_mat), ncol(q_mat))
  }

  T_mi <- U_bar + (1 + 1 / m) * B_mi
  se <- sqrt(diag(T_mi))

  U_diag <- diag(U_bar)
  B_diag <- diag(B_mi)
  r <- rep(0, length(se))

  has_U <- is.finite(U_diag) & U_diag > 0
  r[has_U] <- ((1 + 1 / m) * B_diag[has_U]) / U_diag[has_U]
  r[!has_U & is.finite(B_diag) & B_diag > 0] <- Inf

  df <- rep(Inf, length(se))
  finite_pos_r <- is.finite(r) & r > 0
  df[finite_pos_r] <- (m - 1) * (1 + 1 / r[finite_pos_r])^2
  df[is.infinite(r)] <- m - 1

  fmi <- r / (1 + r)
  fmi[is.nan(fmi)] <- 0

  est <- if (is.null(est_override)) q_bar else est_override
  ci <- .efa_pooled_make_ci(est, se, df, alpha)

  list(
    est = est,
    q_bar = q_bar,
    se = se,
    ci = ci,
    U_bar = U_bar,
    B_mi = B_mi,
    T_mi = T_mi,
    RIV = r,
    FMI = fmi,
    df = df
  )
}

.efa_pooled_align_unrot_boot <- function(L_boot, target, align_unrotated) {
  # Align one bootstrap unrotated loading matrix to the pooled unrotated target.
  if (align_unrotated == "none") {
    return(as.matrix(L_boot))
  }

  if (align_unrotated == "procrustes") {
    return(.change_class(PROCRUSTES(
      A = as.matrix(L_boot),
      Target = as.matrix(target),
      rotation = "orthogonal"
    )$loadings, "matrix"))
  }

  .align_solution(
    L_target = as.matrix(target),
    L = as.matrix(L_boot)
  )$loadings
}

.efa_pooled_align_rot_boot <- function(L_boot_unrot, final_target, rotation_type,
                                       procrustes_args = list()) {
  # Re-rotate/re-align the bootstrap unrotated loadings to the final MI target.
  # This is the key step that makes within-imputation bootstrap variances refer
  # to the same estimand as the pooled rotated solution.
  if (rotation_type == "none") {
    return(NULL)
  }

  pr <- do.call(
    PROCRUSTES,
    c(list(A = as.matrix(L_boot_unrot),
           Target = as.matrix(final_target),
           rotation = rotation_type),
      procrustes_args)
  )

  if (!is.null(pr$convergence) && isFALSE(pr$convergence)) {
    return(NULL)
  }

  list(
    Lambda = .change_class(pr$loadings, "matrix"),
    Phi = if (!is.null(pr$Phi)) as.matrix(pr$Phi) else NULL
  )
}

.efa_pooled_rubin_matrix_result <- function(pool, nrow, ncol, dimnames = NULL) {
  # Convert vectorized Rubin output back into matrix-shaped estimates.
  est <- matrix(pool$est, nrow = nrow, ncol = ncol)
  se <- matrix(pool$se, nrow = nrow, ncol = ncol)
  lower <- matrix(pool$ci$lower, nrow = nrow, ncol = ncol)
  upper <- matrix(pool$ci$upper, nrow = nrow, ncol = ncol)

  if (!is.null(dimnames)) {
    dimnames(est) <- dimnames
    dimnames(se) <- dimnames
    dimnames(lower) <- dimnames
    dimnames(upper) <- dimnames
  }

  list(est = est, se = se, ci = list(lower = lower, upper = upper))
}

.efa_pooled_rubin_symmetric_result <- function(pool, k, dimnames = NULL) {
  # Convert vech-based Rubin output back into symmetric matrix form.
  est <- .efa_pooled_unvech(pool$est, k)
  se <- .efa_pooled_unvech(pool$se, k)
  lower <- .efa_pooled_unvech(pool$ci$lower, k)
  upper <- .efa_pooled_unvech(pool$ci$upper, k)

  if (!is.null(dimnames)) {
    dimnames(est) <- dimnames
    dimnames(se) <- dimnames
    dimnames(lower) <- dimnames
    dimnames(upper) <- dimnames
  }

  list(est = est, se = se, ci = list(lower = lower, upper = upper))
}

.efa_pooled_scalar_boot_summary <- function(x, alpha = 0.05) {
  # Percentile-style summaries for scalar quantities generated by the full
  # pooled bootstrap algorithm.
  x <- as.matrix(x)
  se <- apply(x, 2L, stats::sd, na.rm = TRUE)
  ci <- list(
    lower = apply(x, 2L, stats::quantile, probs = alpha / 2, na.rm = TRUE),
    upper = apply(x, 2L, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  )
  list(se = se, ci = ci)
}


.efa_pooled_bootstrap_pool <- function(fits,
                                       orig_R_list,
                                       unrot_loadings_aligned,
                                       mean_unrot_loadings,
                                       rot_loadings = NULL,
                                       phis = NULL,
                                       structure_loadings = NULL,
                                       mean_structure_loadings = NULL,
                                       final_target = NULL,
                                       rotation_type = c("none", "orthogonal", "oblique"),
                                       align_unrotated = c("signed_tucker_congruence", "none", "procrustes"),
                                       procrustes_args = list(),
                                       h2,
                                       residuals,
                                       fit_indices,
                                       pooled_orig_R,
                                       N,
                                       method,
                                       pool_method = "D2",
                                       rmsea_ci_level = .90,
                                       alpha = 0.05,
                                       rmsr_upper = TRUE) {

  rotation_type <- match.arg(rotation_type)
  align_unrotated <- match.arg(align_unrotated)

  has_boot <- .efa_pooled_has_boot_arrays(fits)
  if (!all(has_boot)) {
    if (any(has_boot)) {
      warning(
        "Bootstrap arrays were found for only some imputations. Pooled bootstrap SEs/CIs were not computed.",
        call. = FALSE
      )
    }
    return(NULL)
  }

  m <- length(fits)
  p_vars <- nrow(as.matrix(mean_unrot_loadings))
  k <- ncol(as.matrix(mean_unrot_loadings))
  loading_dimnames <- dimnames(as.matrix(mean_unrot_loadings))

  B_vec <- vapply(fits, function(x) dim(x$boot.arrays$unrot_loadings)[3], integer(1))
  if (length(unique(B_vec)) > 1L) {
    B_use <- min(B_vec)
    warning(
      "The number of bootstrap replicates differs across imputations. Using the minimum number available in all imputations.",
      call. = FALSE
    )
  } else {
    B_use <- B_vec[[1]]
  }

  q_unrot <- lapply(unrot_loadings_aligned, .efa_pooled_vec)
  boot_unrot <- vector("list", m)

  q_rot <- boot_rot <- NULL
  q_phi <- boot_phi <- NULL
  q_structure <- boot_structure <- NULL
  q_h2 <- vector("list", m)
  boot_h2 <- vector("list", m)
  q_residuals <- vector("list", m)
  boot_residuals <- vector("list", m)

  if (rotation_type != "none") {
    q_rot <- lapply(rot_loadings, .efa_pooled_vec)
    boot_rot <- vector("list", m)
  }
  if (rotation_type == "oblique") {
    q_phi <- lapply(phis, .efa_pooled_vech)
    boot_phi <- vector("list", m)
    q_structure <- lapply(structure_loadings, .efa_pooled_vec)
    boot_structure <- vector("list", m)
  }

  scalar_fit_names <- names(fit_indices)[vapply(fit_indices, function(x) {
    is.numeric(x) && length(x) == 1L
  }, logical(1))]
  can_boot_fit_algorithm <- length(scalar_fit_names) > 0L &&
    all(vapply(fits, function(x) {
      !is.null(x$boot.arrays$fit_indices) && !is.null(x$boot.arrays$residuals)
    }, logical(1)))
  pooled_fit_boot <- NULL
  if (can_boot_fit_algorithm) {
    pooled_fit_boot <- matrix(
      NA_real_, nrow = B_use, ncol = length(scalar_fit_names),
      dimnames = list(NULL, scalar_fit_names)
    )
  }

  nonconv_procrustes <- integer(m)

  for (d in seq_len(m)) {
    arr <- fits[[d]]$boot.arrays
    unrot_arr <- arr$unrot_loadings

    boot_unrot[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars * k)

    if (rotation_type != "none") {
      boot_rot[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars * k)
    }
    if (rotation_type == "oblique") {
      boot_phi[[d]] <- matrix(NA_real_, nrow = B_use, ncol = k * (k + 1) / 2)
      boot_structure[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars * k)
    }

    boot_h2[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars)
    boot_residuals[[d]] <- matrix(NA_real_, nrow = B_use, ncol = p_vars * p_vars)

    if (rotation_type == "oblique") {
      q_h2[[d]] <- .efa_pooled_communalities(rot_loadings[[d]], phis[[d]])
      q_residuals[[d]] <- .efa_pooled_vec(.efa_pooled_residual_from_solution(
        orig_R_list[[d]], rot_loadings[[d]], phis[[d]]
      ))
    } else if (rotation_type == "orthogonal") {
      q_h2[[d]] <- .efa_pooled_communalities(rot_loadings[[d]], NULL)
      q_residuals[[d]] <- .efa_pooled_vec(.efa_pooled_residual_from_solution(
        orig_R_list[[d]], rot_loadings[[d]], NULL
      ))
    } else {
      q_h2[[d]] <- .efa_pooled_communalities(unrot_loadings_aligned[[d]], NULL)
      q_residuals[[d]] <- .efa_pooled_vec(.efa_pooled_residual_from_solution(
        orig_R_list[[d]], unrot_loadings_aligned[[d]], NULL
      ))
    }

    for (b in seq_len(B_use)) {
      Lb_unrot0 <- unrot_arr[, , b]
      Lb_unrot <- .efa_pooled_align_unrot_boot(
        L_boot = Lb_unrot0,
        target = mean_unrot_loadings,
        align_unrotated = align_unrotated
      )
      boot_unrot[[d]][b, ] <- .efa_pooled_vec(Lb_unrot)

      if (rotation_type != "none") {
        # Do not use arr$rot_loadings here. Those matrices were aligned to the
        # imputation-specific EFA target. MI pooling requires every bootstrap
        # replicate to be represented in the final MI target space.
        rot_b <- .efa_pooled_align_rot_boot(
          L_boot_unrot = Lb_unrot0,
          final_target = final_target,
          rotation_type = rotation_type,
          procrustes_args = procrustes_args
        )

        if (!is.null(rot_b)) {
          boot_rot[[d]][b, ] <- .efa_pooled_vec(rot_b$Lambda)

          if (rotation_type == "oblique") {
            boot_phi[[d]][b, ] <- .efa_pooled_vech(rot_b$Phi)
            boot_structure[[d]][b, ] <- .efa_pooled_vec(rot_b$Lambda %*% rot_b$Phi)
            boot_h2[[d]][b, ] <- .efa_pooled_communalities(rot_b$Lambda, rot_b$Phi)
          } else {
            boot_h2[[d]][b, ] <- .efa_pooled_communalities(rot_b$Lambda, NULL)
          }
        } else {
          nonconv_procrustes[d] <- nonconv_procrustes[d] + 1L
        }
      } else {
        boot_h2[[d]][b, ] <- .efa_pooled_communalities(Lb_unrot, NULL)
      }

      if (!is.null(arr$residuals)) {
        boot_residuals[[d]][b, ] <- .efa_pooled_vec(arr$residuals[, , b])
      }
    }
  }

  if (any(nonconv_procrustes > 0L)) {
    warning(
      "Some bootstrap Procrustes rotations did not converge. ",
      "Pooled bootstrap SEs/CIs are based on the valid aligned replicates.",
      call. = FALSE
    )
  }

  pool_unrot <- .efa_pooled_rubin_pool(q_unrot, boot_unrot, alpha = alpha)
  pool_h2 <- .efa_pooled_rubin_pool(
    q_h2, boot_h2, alpha = alpha,
    est_override = as.vector(h2)
  )
  pool_residuals <- .efa_pooled_rubin_pool(
    q_residuals, boot_residuals, alpha = alpha,
    est_override = .efa_pooled_vec(residuals)
  )

  unrot_res <- .efa_pooled_rubin_matrix_result(
    pool_unrot, p_vars, k, loading_dimnames
  )
  residual_res <- .efa_pooled_rubin_matrix_result(
    pool_residuals, p_vars, p_vars, dimnames(as.matrix(residuals))
  )

  SE <- list(
    h2 = pool_h2$se,
    unrot_loadings = unrot_res$se,
    residuals = residual_res$se
  )

  CI <- list(
    h2 = list(lower = pool_h2$ci$lower, upper = pool_h2$ci$upper),
    unrot_loadings = unrot_res$ci,
    residuals = residual_res$ci
  )

  arrays <- list(
    unrot_loadings = boot_unrot,
    h2 = boot_h2,
    residuals = boot_residuals
  )

  MI <- list(
    unrot_loadings = list(RIV = pool_unrot$RIV, FMI = pool_unrot$FMI, df = pool_unrot$df),
    h2 = list(RIV = pool_h2$RIV, FMI = pool_h2$FMI, df = pool_h2$df),
    residuals = list(RIV = pool_residuals$RIV, FMI = pool_residuals$FMI, df = pool_residuals$df),
    bootstrap_rotation_failures = nonconv_procrustes,
    bootstrap_rotation_valid = B_use - nonconv_procrustes
  )

  if (rotation_type != "none") {
    pool_rot <- .efa_pooled_rubin_pool(q_rot, boot_rot, alpha = alpha)
    rot_res <- .efa_pooled_rubin_matrix_result(pool_rot, p_vars, k, loading_dimnames)
    SE$rot_loadings <- rot_res$se
    CI$rot_loadings <- rot_res$ci
    arrays$rot_loadings <- boot_rot
    MI$rot_loadings <- list(RIV = pool_rot$RIV, FMI = pool_rot$FMI, df = pool_rot$df)
  }

  if (rotation_type == "oblique") {
    pool_phi <- .efa_pooled_rubin_pool(q_phi, boot_phi, alpha = alpha)
    phi_res <- .efa_pooled_rubin_symmetric_result(
      pool_phi, k, dimnames(as.matrix(phis[[1]]))
    )
    SE$Phi <- phi_res$se
    CI$Phi <- phi_res$ci
    arrays$Phi <- boot_phi
    MI$Phi <- list(RIV = pool_phi$RIV, FMI = pool_phi$FMI, df = pool_phi$df)

    pool_structure <- .efa_pooled_rubin_pool(
      q_structure, boot_structure, alpha = alpha,
      est_override = .efa_pooled_vec(mean_structure_loadings)
    )
    structure_res <- .efa_pooled_rubin_matrix_result(
      pool_structure, p_vars, k, loading_dimnames
    )
    SE$Structure <- structure_res$se
    CI$Structure <- structure_res$ci
    arrays$Structure <- boot_structure
    MI$Structure <- list(
      RIV = pool_structure$RIV,
      FMI = pool_structure$FMI,
      df = pool_structure$df
    )
  }

  ## Fit-index summaries ------------------------------------------------------
  ## fit_indices_descriptive: Rubin-Wald summaries of per-imputation fit indices.
  ## fit_indices_pooled_algorithm: percentile summaries obtained by re-running
  ## the pooled residual/D2 calculations bootstrap-index by bootstrap-index.

  fit_q_names <- Reduce(intersect, lapply(fits, function(x) names(x$fit_indices)))
  fit_q_names <- fit_q_names[vapply(fits[[1]]$fit_indices[fit_q_names], function(x) {
    is.numeric(x) && length(x) == 1L
  }, logical(1))]

  if (length(fit_q_names) > 0L && all(vapply(fits, function(x) {
    !is.null(x$boot.arrays$fit_indices)
  }, logical(1)))) {
    q_fit <- lapply(fits, function(x) unlist(x$fit_indices[fit_q_names]))
    boot_fit <- vector("list", m)

    for (d in seq_len(m)) {
      fit_arr <- as.matrix(fits[[d]]$boot.arrays$fit_indices)
      if (is.null(colnames(fit_arr))) {
        colnames(fit_arr) <- names(fits[[d]]$fit_indices)
      }
      boot_fit[[d]] <- fit_arr[seq_len(B_use), fit_q_names, drop = FALSE]
    }

    pool_fit_desc <- .efa_pooled_rubin_pool(q_fit, boot_fit, alpha = alpha)
    SE$fit_indices_descriptive <- stats::setNames(pool_fit_desc$se, fit_q_names)
    CI$fit_indices_descriptive <- list(
      lower = stats::setNames(pool_fit_desc$ci$lower, fit_q_names),
      upper = stats::setNames(pool_fit_desc$ci$upper, fit_q_names)
    )
    arrays$fit_indices_descriptive <- boot_fit
    MI$fit_indices_descriptive <- list(
      RIV = stats::setNames(pool_fit_desc$RIV, fit_q_names),
      FMI = stats::setNames(pool_fit_desc$FMI, fit_q_names),
      df = stats::setNames(pool_fit_desc$df, fit_q_names)
    )
  }

  if (!is.null(pooled_fit_boot)) {
    # Use bootstrap index b as a Monte Carlo pairing across imputations. The
    # bootstrap samples are independent across imputations, so the pairing has
    # no substantive meaning; it simply creates B replicate runs of the complete
    # pooled-fit algorithm.
    for (b in seq_len(B_use)) {
      Lb_unrot_list <- vector("list", m)
      Lb_rot_list <- if (rotation_type != "none") vector("list", m) else NULL
      Phib_list <- if (rotation_type == "oblique") vector("list", m) else NULL
      Rb_list <- vector("list", m)
      boot_fit_objects <- vector("list", m)

      for (d in seq_len(m)) {
        arr <- fits[[d]]$boot.arrays
        Lb_unrot <- matrix(boot_unrot[[d]][b, ], nrow = p_vars, ncol = k)
        Lb_unrot_list[[d]] <- Lb_unrot

        if (rotation_type != "none") {
          Lrow <- boot_rot[[d]][b, ]
          if (anyNA(Lrow)) {
            next
          }
          Lb_rot <- matrix(Lrow, nrow = p_vars, ncol = k)
          Lb_rot_list[[d]] <- Lb_rot

          if (rotation_type == "oblique") {
            Phirow <- boot_phi[[d]][b, ]
            if (anyNA(Phirow)) {
              next
            }
            Phib <- .efa_pooled_unvech(Phirow, k)
            Phib_list[[d]] <- Phib
            model_b <- .efa_pooled_model_implied(Lb_rot, Phib)
          } else {
            model_b <- .efa_pooled_model_implied(Lb_rot, NULL)
          }
        } else {
          model_b <- .efa_pooled_model_implied(Lb_unrot, NULL)
        }

        # EFA() stores bootstrap residuals. Adding the model-implied matrix for
        # the aligned bootstrap solution reconstructs the bootstrap observed R;
        # because rotations preserve the reproduced correlation matrix, this is
        # consistent with the final target-space representation.
        Rb <- arr$residuals[, , b] + model_b
        diag(Rb) <- 1
        Rb_list[[d]] <- Rb

        fit_arr <- as.matrix(arr$fit_indices)
        if (is.null(colnames(fit_arr))) {
          colnames(fit_arr) <- names(fits[[d]]$fit_indices)
        }
        boot_fit_objects[[d]] <- list(fit_indices = as.list(fit_arr[b, , drop = TRUE]))
      }

      if (any(vapply(Rb_list, is.null, logical(1)))) {
        next
      }

      pooled_R_b <- .average_matrices(Rb_list)

      if (rotation_type == "oblique") {
        if (any(vapply(Lb_rot_list, is.null, logical(1))) ||
            any(vapply(Phib_list, is.null, logical(1)))) {
          next
        }
        L_pool_b <- .average_matrices(Lb_rot_list)
        Phi_pool_b <- .average_matrices(Phib_list)
        model_pool_b <- .efa_pooled_model_implied(L_pool_b, Phi_pool_b)
      } else if (rotation_type == "orthogonal") {
        if (any(vapply(Lb_rot_list, is.null, logical(1)))) {
          next
        }
        L_pool_b <- .average_matrices(Lb_rot_list)
        model_pool_b <- .efa_pooled_model_implied(L_pool_b, NULL)
      } else {
        L_pool_b <- .average_matrices(Lb_unrot_list)
        model_pool_b <- .efa_pooled_model_implied(L_pool_b, NULL)
      }

      residual_b <- pooled_R_b - model_pool_b
      diag(residual_b) <- 0
      RMSR_b <- .rmsr(residual_b, upper = rmsr_upper)

      if (any(vapply(boot_fit_objects, is.null, logical(1)))) {
        next
      }

      fit_b <- .efa_pooled_fit_indices(
        fits = boot_fit_objects,
        pooled_R = pooled_R_b,
        residuals = residual_b,
        RMSR = RMSR_b,
        N = N,
        method = method,
        pool_method = pool_method,
        rmsea_ci_level = rmsea_ci_level
      )

      pooled_fit_boot[b, ] <- unlist(fit_b[scalar_fit_names])
    }

    fit_boot_summary <- .efa_pooled_scalar_boot_summary(pooled_fit_boot, alpha = alpha)
    SE$fit_indices_pooled_algorithm <- fit_boot_summary$se
    CI$fit_indices_pooled_algorithm <- fit_boot_summary$ci
    arrays$fit_indices_pooled_algorithm <- pooled_fit_boot
  }

  list(
    SE = SE,
    CI = CI,
    arrays = arrays,
    MI = MI
  )
}
