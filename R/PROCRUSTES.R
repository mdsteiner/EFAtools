#' Rotate a loading matrix to a target using Procrustes alignment
#'
#' `PROCRUSTES()` aligns one loading matrix to a target loading matrix with the
#' same dimensions. It is used internally by `EFA_POOLED()`, but can also be used
#' directly when factor columns must be brought into a common orientation before
#' averaging or comparing solutions.
#'
#' For `rotation = "orthogonal"`, the function solves the closed-form orthogonal
#' Procrustes problem
#'
#' \deqn{\min_T \frac{1}{2}\|A T - B\|_F^2 \quad \textrm{subject to}\quad T'T = I,}
#'
#' where `A` is the loading matrix and `B` is `Target`.
#'
#' For `rotation = "oblique"`, the function calls the compiled
#' `.oblique_procrustes()` optimizer. The oblique convention is the same as in
#' `GPArotation::targetQ()`:
#'
#' \deqn{L = A T^{-T}, \qquad \Phi = T'T, \qquad diag(\Phi) = 1.}
#'
#' Random starts are only used for oblique alignment. For one-factor models,
#' oblique and orthogonal alignment are equivalent, so the function uses the
#' stable one-factor orthogonal solution instead of calling the oblique optimizer.
#'
#' @param A Numeric loading matrix to be aligned.
#' @param Target Numeric target matrix with the same dimensions as `A`.
#' @param rotation Character string, either `"orthogonal"` or `"oblique"`.
#' @param S Optional `k x k` cross-product matrix `crossprod(A)`. Supplying this
#'   is useful when the same `A` is rotated repeatedly. `S` is used only when
#'   `oblique_normalize = FALSE`; if Kaiser normalization is requested, the
#'   cross-product must be recomputed on the normalized matrix.
#' @param T_init Optional `k x k` nonsingular starting transformation matrix for
#'   the oblique solver. Its columns are normalized internally.
#' @param oblique_eps Positive convergence tolerance for the projected-gradient
#'   norm in the oblique solver.
#' @param oblique_maxit Non-negative integer. Maximum number of projected-gradient
#'   updates in the full oblique solver.
#' @param oblique_max_line_search Non-negative integer. Maximum number of
#'   step-halving attempts after the initial line-search step.
#' @param oblique_step0 Positive initial step size for the oblique solver.
#' @param oblique_normalize Logical; if `TRUE`, apply Kaiser row normalization in
#'   the oblique solver and back-transform the aligned loadings afterwards.
#' @param oblique_random_starts Non-negative integer. Number of additional random
#'   starts used by the oblique solver.
#' @param oblique_screen_keep Non-negative integer. Number of random starts
#'   retained after cheap objective screening and sent to triage optimization.
#' @param oblique_triage_maxit Non-negative integer. Number of short optimization
#'   iterations used in the triage stage.
#' @param oblique_triage_improve_tol Non-negative scalar. Relative improvement
#'   required for a triaged start to be promoted to full optimization.
#'
#' @returns A list containing aligned `loadings`, transformation matrix `T`,
#'   factor intercorrelation matrix `Phi`, target criterion `value`, convergence
#'   diagnostics, line-search diagnostics, and multi-start summaries. Row and
#'   column names are preserved where possible.
#'
#' @export
#'
PROCRUSTES <- function(A,
                       Target,
                       rotation = c("orthogonal", "oblique"),
                       S = NULL,
                       T_init = NULL,
                       oblique_eps = 1e-5,
                       oblique_maxit = 1000,
                       oblique_max_line_search = 10,
                       oblique_step0 = 1,
                       oblique_normalize = FALSE,
                       oblique_random_starts = 0,
                       oblique_screen_keep = 2,
                       oblique_triage_maxit = 25,
                       oblique_triage_improve_tol = 0) {
  rotation <- match.arg(rotation)

  mats <- .procrustes_validate_matrix_pair(A, Target)
  A <- mats$A
  Target <- mats$B
  k <- ncol(A)

  if (rotation == "orthogonal" || k == 1L) {
    out <- .orthogonal_procrustes(A, Target)
    if (rotation == "oblique" && k == 1L) {
      out$method <- "single_factor_procrustes"
    }
    return(out)
  }

  controls <- .procrustes_validate_oblique_controls(
    oblique_eps = oblique_eps,
    oblique_maxit = oblique_maxit,
    oblique_max_line_search = oblique_max_line_search,
    oblique_step0 = oblique_step0,
    oblique_normalize = oblique_normalize,
    oblique_random_starts = oblique_random_starts,
    oblique_screen_keep = oblique_screen_keep,
    oblique_triage_maxit = oblique_triage_maxit,
    oblique_triage_improve_tol = oblique_triage_improve_tol
  )

  S_arg <- if (is.null(S)) NULL else .procrustes_validate_crossprod(S, k)
  T_arg <- if (is.null(T_init)) NULL else .procrustes_validate_t_init(T_init, k)

  out <- .oblique_procrustes(
    A = A,
    B = Target,
    S_r = S_arg,
    T_init_r = T_arg,
    eps = controls$oblique_eps,
    maxit = controls$oblique_maxit,
    max_line_search = controls$oblique_max_line_search,
    step0 = controls$oblique_step0,
    normalize = controls$oblique_normalize,
    random_starts = controls$oblique_random_starts,
    screen_keep = controls$oblique_screen_keep,
    triage_maxit = controls$oblique_triage_maxit,
    triage_improve_tol = controls$oblique_triage_improve_tol
  )

  out$method <- "oblique_procrustes"
  .procrustes_apply_dimnames(out, A, Target)
}


#' Consensus Procrustes alignment across multiple loading matrices
#'
#' Align several loading matrices to a common Procrustes consensus target.
#'
#' The function iterates between two steps:
#'
#' 1. each loading matrix is independently aligned to the current target with
#'    `PROCRUSTES()`; and
#' 2. the target is updated to the centroid, i.e., the elementwise average of the
#'    aligned matrices.
#'
#' This makes the target symmetric across imputations or samples: no single
#' solution is permanently privileged as the reference. The outer loop can stop
#' when the target stabilizes, when the consensus loss stabilizes, or when both
#' criteria are satisfied.
#'
#' If `multi_start = FALSE`, one consensus run is performed. If
#' `multi_start = TRUE`, the same single-start engine is repeated for the
#' selected starting targets, and the run with the smallest final mean loss is
#' returned as the main result. All runs and a between-run congruence summary are
#' retained in the `multi_start` component.
#'
#' @param unrotated_list List of unrotated loading matrices to be aligned. All
#'   matrices must be numeric, finite, and have identical dimensions.
#' @param init_targets Optional list of starting target matrices. These are
#'   typically rotated loading matrices from the corresponding analyses. If
#'   `NULL`, `unrotated_list` is used.
#' @param rotation Character string, either `"orthogonal"` or `"oblique"`.
#' @param start Either a single integer selecting an element of `init_targets`,
#'   or an explicit target matrix. Used when `multi_start = FALSE`.
#' @param multi_start Logical. If `FALSE`, perform one consensus-target run. If
#'   `TRUE`, repeat the single-start algorithm for each element of `starts`.
#' @param starts Integer vector selecting elements of `init_targets` used as
#'   starting targets when `multi_start = TRUE`. If `NULL`, all elements of
#'   `init_targets` are used. Duplicate entries are removed.
#' @param tol Positive relative Frobenius-norm convergence tolerance for the
#'   outer target update.
#' @param loss_tol Positive tolerance for the relative change in the outer
#'   consensus loss. If `NULL`, loss-based convergence is disabled. It cannot be
#'   `NULL` when `convergence` is `"loss"` or `"both"`.
#' @param loss_patience Positive integer. Number of consecutive iterations with
#'   relative loss change below `loss_tol` required for loss-based convergence.
#' @param convergence Character string controlling the stopping rule. `"either"`
#'   stops when either target or loss convergence is satisfied; `"target"` uses
#'   only target change; `"loss"` uses only loss change; `"both"` requires both.
#' @param min_iter Non-negative integer. Minimum number of outer iterations
#'   before convergence can be declared.
#' @param max_iter Positive integer. Maximum number of outer consensus
#'   iterations.
#' @param alpha Damping factor for the target update. `alpha = 1` uses the full
#'   centroid update. Smaller values, such as `0.5`, can reduce oscillation.
#' @param match_target Logical. If `TRUE`, the updated centroid is signed and
#'   column-matched to the previous target before convergence is evaluated.
#' @param hyper_cutoff Non-negative cutoff used by `.hyperplane_count()` for
#'   summary output.
#' @param oblique_maxit,oblique_eps,oblique_max_line_search,oblique_step0,oblique_normalize Parameters
#'   passed to `PROCRUSTES()` when `rotation = "oblique"`.
#' @param oblique_random_starts Number of random starts used by the inner oblique
#'   solver. See `oblique_random_starts_stage` for when these are used.
#' @param oblique_random_starts_stage Character string controlling whether random
#'   starts are used during the outer consensus loop, the final alignment pass,
#'   both, or neither. The default `"final"` keeps the outer loop smooth by using
#'   warm starts during iteration and applies random-start protection only in the
#'   final pass.
#' @param oblique_screen_keep,oblique_triage_maxit,oblique_triage_improve_tol
#'   Screening and triage parameters passed to the compiled oblique solver when
#'   random starts are used.
#' @param verbose Logical; if `TRUE`, print convergence messages for the outer
#'   loop.
#'
#' @returns A list with the converged target, aligned matrices, pooled loadings,
#'   pooled `Phi`, convergence history, inner-alignment diagnostics, and
#'   hyperplane summaries. If `multi_start = TRUE`, the `multi_start` element also
#'   contains the per-start losses, convergence indicators, run summaries, all
#'   run objects, and between-run Tucker congruence matrices.
#'
#' @export
#'
#' @references
#' Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika*, 40,
#' 33-51.
CONSENSUS_PROCRUSTES <- function(unrotated_list,
                                 init_targets = NULL,
                                 rotation = c("orthogonal", "oblique"),
                                 start = 1,
                                 multi_start = FALSE,
                                 starts = NULL,
                                 tol = 1e-3,
                                 loss_tol = 1e-6,
                                 loss_patience = 5,
                                 convergence = c("either", "target", "loss", "both"),
                                 min_iter = 2,
                                 max_iter = 200,
                                 alpha = 1,
                                 match_target = TRUE,
                                 hyper_cutoff = 0.15,
                                 oblique_maxit = 500,
                                 oblique_eps = 1e-5,
                                 oblique_max_line_search = 10,
                                 oblique_step0 = 1,
                                 oblique_normalize = FALSE,
                                 oblique_random_starts = 0,
                                 oblique_random_starts_stage = c("final", "none", "outer", "both"),
                                 oblique_screen_keep = 2,
                                 oblique_triage_maxit = 25,
                                 oblique_triage_improve_tol = 0,
                                 verbose = FALSE) {
  rotation <- match.arg(rotation)
  convergence <- match.arg(convergence)
  oblique_random_starts_stage <- match.arg(oblique_random_starts_stage)
  multi_start <- .procrustes_check_flag(multi_start, "multi_start")

  single_args <- list(
    unrotated_list = unrotated_list,
    init_targets = init_targets,
    rotation = rotation,
    tol = tol,
    loss_tol = loss_tol,
    loss_patience = loss_patience,
    convergence = convergence,
    min_iter = min_iter,
    max_iter = max_iter,
    alpha = alpha,
    match_target = match_target,
    hyper_cutoff = hyper_cutoff,
    oblique_maxit = oblique_maxit,
    oblique_eps = oblique_eps,
    oblique_max_line_search = oblique_max_line_search,
    oblique_step0 = oblique_step0,
    oblique_normalize = oblique_normalize,
    oblique_random_starts = oblique_random_starts,
    oblique_random_starts_stage = oblique_random_starts_stage,
    oblique_screen_keep = oblique_screen_keep,
    oblique_triage_maxit = oblique_triage_maxit,
    oblique_triage_improve_tol = oblique_triage_improve_tol,
    verbose = verbose
  )

  if (!multi_start) {
    out <- do.call(
      .consensus_target_procrustes_single,
      c(single_args, list(start = start))
    )

    out$multi_start <- list(
      enabled = FALSE,
      starts = out$start,
      best_start = out$start,
      best_index = 1L,
      losses = out$mean_loss,
      converged = out$converged,
      iterations = out$iterations,
      summary = data.frame(
        start = out$start,
        loss = out$mean_loss,
        converged = out$converged,
        iterations = out$iterations
      ),
      runs = NULL,
      congruence = NULL
    )

    return(out)
  }

  if (is.null(init_targets)) {
    single_args$init_targets <- unrotated_list
    n_targets <- length(unrotated_list)
  } else {
    n_targets <- length(init_targets)
  }
  starts <- .procrustes_validate_starts(starts, n_targets)

  runs <- lapply(starts, function(s) {
    do.call(
      .consensus_target_procrustes_single,
      c(single_args, list(start = s))
    )
  })

  losses <- vapply(runs, function(x) x$mean_loss, numeric(1L))
  names(losses) <- paste0("start_", starts)

  finite_losses <- is.finite(losses)
  if (!any(finite_losses)) {
    stop("No multi-start consensus run produced a finite final loss.", call. = FALSE)
  }
  best_idx <- which.min(ifelse(finite_losses, losses, Inf))
  best <- runs[[best_idx]]

  n_runs <- length(runs)
  run_names <- paste0("start_", starts)
  congruence <- vector("list", n_runs)
  names(congruence) <- run_names

  for (i in seq_len(n_runs)) {
    congruence[[i]] <- vector("list", n_runs)
    names(congruence[[i]]) <- run_names
  }

  for (i in seq_len(n_runs)) {
    for (j in i:n_runs) {
      cij <- .tucker_congruence(
        runs[[i]]$pooled_loadings,
        runs[[j]]$pooled_loadings
      )
      congruence[[i]][[j]] <- cij
      congruence[[j]][[i]] <- if (i == j) cij else t(cij)
    }
  }

  run_summary <- data.frame(
    start = starts,
    loss = losses,
    converged = vapply(runs, function(x) isTRUE(x$converged), logical(1L)),
    iterations = vapply(runs, function(x) x$iterations, integer(1L)),
    final_failures = vapply(runs, function(x) length(x$final_failures), integer(1L)),
    row.names = run_names
  )

  best$multi_start <- list(
    enabled = TRUE,
    starts = starts,
    best_start = starts[best_idx],
    best_index = best_idx,
    losses = losses,
    converged = run_summary$converged,
    iterations = run_summary$iterations,
    summary = run_summary,
    runs = runs,
    congruence = congruence
  )

  best
}
