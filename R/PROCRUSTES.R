#' Rotate a loading matrix to a target using Procrustes alignment
#'
#' Matrix alignment function. For orthogonal rotation it
#' uses the closed-form orthogonal Procrustes solution. For oblique rotation it
#' calls \link{.oblique_procrustes}, which uses the same oblique
#' transformation convention as \code{\link[GPArotation:rotations]{GPArotation::targetQ}}:
#' \code{L = A \%*\% solve(t(T))} and \code{Phi = t(T) \%*\% T}.
#'
#' Random starts are only used in the oblique case.
#'
#' @param A Numeric loading matrix to be aligned.
#' @param Target Numeric target matrix with the same dimensions as \code{A}.
#' @param rotation Character string, either \code{"orthogonal"} or \code{"oblique"}.
#' @param S Optional \code{k x k} cross-product matrix \code{crossprod(A)}. Supplying this
#'   is useful when the same \code{A} is rotated repeatedly.
#' @param T_init Optional \code{k x k} starting transformation matrix for the oblique
#'   solver.
#' @param oblique_eps Convergence tolerance for the oblique solver.
#' @param oblique_maxit Maximum number of full projected-gradient iterations for
#'   the oblique solver.
#' @param oblique_max_line_search Maximum number of step-halving attempts in each
#'   oblique line-search step.
#' @param oblique_step0 Initial step size for the oblique solver.
#' @param oblique_normalize Logical; if \code{TRUE}, apply Kaiser normalization in the
#'   oblique solver.
#' @param oblique_random_starts Number of additional random starts used by the
#'   oblique solver.
#' @param oblique_screen_keep Number of random starts retained after cheap
#'   screening and sent to triage optimization.
#' @param oblique_triage_maxit Number of short optimization iterations used in
#'   the triage stage.
#' @param oblique_triage_improve_tol Relative improvement required for a triaged
#'   start to be promoted to full optimization.
#'
#' @returns A list containing the rotated loadings, transformation matrix,
#'   \code{Phi}, target criterion value, convergence diagnostics, and multi-start
#'   summaries.
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

  if (rotation == "orthogonal") {
    return(.orthogonal_procrustes(A, Target))
  }

  S_arg <- if (is.null(S)) NULL else as.matrix(S)
  T_arg <- if (is.null(T_init)) NULL else as.matrix(T_init)

  .oblique_procrustes(
    A = as.matrix(A),
    B = as.matrix(Target),
    S_r = S_arg,
    T_init_r = T_arg,
    eps = oblique_eps,
    maxit = as.integer(oblique_maxit),
    max_line_search = as.integer(oblique_max_line_search),
    step0 = oblique_step0,
    normalize = oblique_normalize,
    random_starts = as.integer(oblique_random_starts),
    screen_keep = as.integer(oblique_screen_keep),
    triage_maxit = as.integer(oblique_triage_maxit),
    triage_improve_tol = oblique_triage_improve_tol
  )
}


#' Consensus Procrustes alignment across multiple loading matrices
#'
#' Align several loading matrices to a common Procrustes consensus target.
#'
#' For a fixed current target, each loading matrix is aligned independently using
#' \link{PROCRUSTES}. The target is then updated as the average of
#' the aligned loading matrices. This is repeated until the outer consensus loop
#' converges.
#'
#' Compared with a simple fixed-target approach, the consensus update is
#' symmetric across imputations or samples: no single loading matrix is treated
#' as the permanent reference. Compared with earlier versions of this function,
#' convergence can be assessed not only by the relative Frobenius change in the
#' target, but also by stabilization of the consensus loss. This is useful for
#' oblique rotations, where very small elementwise target changes can be harder
#' to achieve than a stable consensus objective.
#'
#' The function can be used in two modes. If \code{multi_start = FALSE}, a single
#' consensus run is performed using \code{start} as the initial target. If
#' \code{multi_start = TRUE}, the single-start engine is repeated for each element of
#' \code{starts}, and the function returns all runs together with the best run,
#' chosen by the smallest final mean loss.
#'
#' @param unrotated_list List of unrotated loading matrices to be aligned.
#' @param init_targets Optional list of starting target matrices. These are
#'   typically rotated loading matrices from the corresponding analyses. If
#'   \code{NULL}, \code{unrotated_list} is used.
#' @param rotation Character string, either \code{"orthogonal"} or \code{"oblique"}.
#' @param start Either a single integer selecting an element of \code{init_targets},
#'   or an explicit target matrix. Used when \code{multi_start = FALSE}.
#' @param multi_start Logical. If \code{FALSE}, perform one consensus-target run. If
#'   \code{TRUE}, repeat the single-start algorithm for each element of \code{starts}.
#' @param starts Integer vector selecting the elements of \code{init_targets} used as
#'   starting targets when \code{multi_start = TRUE}. If \code{NULL}, all elements of
#'   \code{init_targets} are used. Ignored when \code{multi_start = FALSE}.
#' @param tol Relative Frobenius-norm convergence tolerance for the outer target
#'   update. The default \code{1e-3} is intended as a practical outer-loop tolerance;
#'   the inner oblique optimizer has its own tolerance \code{oblique_eps}.
#' @param loss_tol Optional tolerance for the relative change in the outer
#'   consensus loss. If \code{NULL}, loss-based convergence is disabled.
#' @param loss_patience Number of consecutive iterations with relative loss
#'   change below \code{loss_tol} required for loss-based convergence.
#' @param convergence Character string controlling the stopping rule. \code{"either"}
#'   stops when either \code{tol} or \code{loss_tol} is satisfied; \code{"target"} uses only
#'   the target-change criterion; \code{"loss"} uses only the loss criterion; \code{"both"}
#'   requires both.
#' @param min_iter Minimum number of outer iterations before convergence can be
#'   declared.
#' @param max_iter Maximum number of outer consensus iterations.
#' @param alpha Damping factor for the target update. \code{alpha = 1} uses the full
#'   centroid update. Smaller values, such as \code{0.5}, can reduce oscillation.
#' @param match_target Logical. If \code{TRUE}, the updated centroid is signed and
#'   column-matched to the previous target before convergence is evaluated.
#' @param hyper_cutoff Cutoff used by \link{.hyperplane_count} for summary output.
#' @param oblique_maxit,oblique_eps,oblique_max_line_search,oblique_step0,oblique_normalize Parameters
#' passed to \link{PROCRUSTES} when
#'   \code{rotation = "oblique"}.
#' @param oblique_random_starts Number of random starts used by the inner oblique
#'   solver. See \code{oblique_random_starts_stage} for when these are used.
#' @param oblique_random_starts_stage Character string controlling whether random
#'   starts are used during the outer consensus loop, the final alignment pass,
#'   both, or neither. The default \code{"final"} keeps the outer loop smooth by using
#'   warm starts during iteration and applies random-start protection only in the
#'   final pass.
#' @param oblique_screen_keep,oblique_triage_maxit,oblique_triage_improve_tol
#'   Screening and triage parameters passed to the compiled oblique solver when
#'   random starts are used.
#' @param verbose Logical; if \code{TRUE}, print simple convergence messages for the
#'   outer loop.
#'
#' @returns
#' A list with the converged target, aligned matrices,
#' pooled loadings, pooled \code{Phi}, convergence history, and summary diagnostics.
#' If \code{multi_start = TRUE}, additional infos in the multi_start list with infos on
#' the different starts.
#'
#' @export
#'
#' @references
#' Gower, J. C. (1975). Generalized Procrustes analysis. Psychometrika, 40,
#' 33-51.
CONSENSUS_PROCRUSTES <- function(unrotated_list,
                                 init_targets = NULL,
                                 rotation = c("orthogonal", "oblique"),
                                 start = 1,
                                 multi_start = FALSE,
                                 starts = NULL,
                                 tol = 1e-3,
                                 loss_tol = 1e-7,
                                 loss_patience = 5,
                                 convergence = c("either", "target", "loss", "both"),
                                 min_iter = 2,
                                 max_iter = 200,
                                 alpha = 1,
                                 match_target = TRUE,
                                 hyper_cutoff = 0.15,
                                 oblique_maxit = 300,
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

  if (!isTRUE(multi_start)) {
    out <- do.call(
      .consensus_target_procrustes_single,
      c(single_args, list(start = start))
    )

    out$multi_start <- list(
      enabled = FALSE,
      starts = start,
      best_start = start,
      best_index = 1L,
      losses = out$mean_loss,
      runs = NULL,
      congruence = NULL
    )

    return(out)
  }

  if (is.null(init_targets)) {
    single_args$init_targets <- unrotated_list
  }
  if (is.null(starts)) {
    starts <- seq_along(single_args$init_targets)
  }

  runs <- lapply(starts, function(s) {
    do.call(
      .consensus_target_procrustes_single,
      c(single_args, list(start = s))
    )
  })

  losses <- vapply(runs, function(x) x$mean_loss, numeric(1))
  names(losses) <- paste0("start_", starts)

  best_idx <- which.min(losses)
  best <- runs[[best_idx]]

  n_runs <- length(runs)
  congruence <- vector("list", n_runs)
  names(congruence) <- paste0("start_", starts)

  for (i in seq_len(n_runs)) {
    congruence[[i]] <- vector("list", n_runs)
    names(congruence[[i]]) <- paste0("start_", starts)

    for (j in seq_len(n_runs)) {
      congruence[[i]][[j]] <- .tucker_congruence(
        runs[[i]]$pooled_loadings,
        runs[[j]]$pooled_loadings
      )
    }
  }

  best$multi_start <- list(
    enabled = TRUE,
    starts = starts,
    best_start = starts[best_idx],
    best_index = best_idx,
    losses = losses,
    runs = runs,
    congruence = congruence
  )

  best
}
