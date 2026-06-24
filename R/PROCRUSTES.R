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
#' By default the oblique solver is warm-started from the closed-form orthogonal
#' Procrustes solution, which resolves the factor permutation and sign
#' indeterminacy and avoids the poor local minima an identity start can fall
#' into. Supply `T_init` to override this start. Random starts are only used for
#' oblique alignment. For one-factor models, oblique and orthogonal alignment are
#' equivalent, so the function uses the stable one-factor orthogonal solution
#' instead of calling the oblique optimizer.
#'
#' @param A Numeric loading matrix to be aligned.
#' @param Target Numeric target matrix with the same dimensions as `A`.
#' @param rotation Character string, either `"orthogonal"` or `"oblique"`.
#' @param S Optional `k x k` cross-product matrix `crossprod(A)`. Supplying this
#'   is useful when the same `A` is rotated repeatedly. `S` is used only when
#'   `oblique_normalize = FALSE`; if Kaiser normalization is requested, the
#'   cross-product must be recomputed on the normalized matrix.
#' @param T_init Optional `k x k` nonsingular starting transformation matrix for
#'   the oblique solver. Its columns are normalized internally. If `NULL` (the
#'   default), the oblique solver is warm-started from the closed-form orthogonal
#'   Procrustes solution.
#' @param oblique_eps Positive convergence tolerance for the projected-gradient
#'   norm in the oblique solver.
#' @param oblique_maxit Non-negative integer. Maximum number of projected-gradient
#'   updates in the full oblique solver.
#' @param oblique_max_line_search Non-negative integer. Maximum number of
#'   step-halving attempts after the initial line-search step.
#' @param oblique_step0 Positive initial step size for the oblique solver.
#' @param oblique_normalize Logical; if `TRUE`, apply Kaiser row normalization to
#'   the loadings (only) in the oblique solver and back-transform the aligned
#'   loadings afterwards, leaving `Target` unnormalized (as in
#'   `GPArotation::targetQ(normalize = TRUE)`).
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
#'   column names are preserved where possible. When `oblique_normalize = TRUE`
#'   the returned `loadings` are back-transformed to the original scale, but
#'   `value` is the criterion on the Kaiser-normalized loadings, so it is not
#'   `0.5 * sum((loadings - Target)^2)`.
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

  # Warm-start the oblique solver from the closed-form orthogonal Procrustes
  # solution when the user supplies no start. The identity start can be trapped
  # in poor local minima on well-conditioned problems (and trapped solutions
  # still report convergence), whereas the orthogonal solution resolves the
  # factor permutation and sign structure and lands in the basin of the global
  # oblique optimum. The consensus engine warm-starts its inner alignments the
  # same way. An explicit T_init overrides this.
  if (is.null(T_arg)) {
    if (controls$oblique_normalize) {
      # Compute the start in the same Kaiser row-normalized space the compiled
      # solver optimizes in: the loadings are row-weighted, the target is left raw
      # (mirrors the row weighting in .oblique_procrustes()).
      w <- sqrt(rowSums(A^2))
      w[!is.finite(w) | w < 1e-15] <- 1
      T_arg <- .procrustes_orthogonal_T(A / w, Target)
    } else {
      T_arg <- .procrustes_orthogonal_T(A, Target)
    }
  }

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


#' Generalized Procrustes Analysis consensus target across loading matrices
#'
#' Internal helper that constructs a Generalized Procrustes Analysis (GPA)
#' consensus target across a list of loading matrices and returns the aligned
#' loadings, the centroid target, and convergence diagnostics. Used by
#' [EFA_POOLED()] under `target_method = "consensus"` to build a common
#' rotation target across imputations. Oblique rotations are not supported
#' here: the iteration is degenerate for oblique transforms with more than
#' one factor (cf. Lorenzo-Seva & Van Ginkel 2016, who use a Promin step on
#' top of the centroid rather than iterated oblique Procrustes); callers
#' should pass the unrotated solutions of an orthogonal rotation, or use
#' `target_method = "first_target"`.
#'
#' The iteration alternates two steps:
#'
#' 1. each loading matrix is aligned to the current target with `PROCRUSTES()`;
#' 2. the target is updated to the elementwise centroid of the aligned matrices.
#'
#' The outer loop stops when the target stabilises, when the consensus loss
#' stabilises, or when both criteria are satisfied.
#'
#' If `multi_start = FALSE`, one consensus run is performed. If
#' `multi_start = TRUE`, the same engine is repeated for the selected starting
#' targets and the run with the smallest final mean loss is returned as the
#' main result; all runs and a between-run congruence summary are retained in
#' the `multi_start` component.
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
#' @param verbose Logical; if `TRUE`, print convergence messages for the outer
#'   loop.
#'
#' @returns A list with the converged target, aligned matrices, pooled loadings,
#'   pooled `Phi`, convergence history, inner-alignment diagnostics, and
#'   hyperplane summaries. If `multi_start = TRUE`, the `multi_start` element also
#'   contains the per-start losses, convergence indicators, run summaries, all
#'   run objects, and between-run Tucker congruence matrices.
#'
#' @references
#' Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika*, 40,
#' 33-51.
#'
#' Van Ginkel, J. R., & Kroonenberg, P. M. (2014). Using Generalized
#' Procrustes Analysis for Multiple Imputation in Principal Component
#' Analysis. *Journal of Classification*, 31, 242-269.
#'
#' Lorenzo-Seva, U., & Van Ginkel, J. R. (2016). Multiple Imputation of
#' missing values in exploratory factor analysis of multidimensional scales:
#' estimating latent trait scores. *Anales de Psicologia*, 32, 596-608.
#'
.gpa_consensus_target <- function(unrotated_list,
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
                                 verbose = FALSE) {
  rotation <- match.arg(rotation)
  convergence <- match.arg(convergence)
  multi_start <- .procrustes_check_flag(multi_start, "multi_start")

  # The iteration is well-behaved for orthogonal Procrustes (closed-form per
  # step; classical Gower 1975 GPA) but the extra column-skew freedom of an
  # oblique transform lets the centroid loss collapse to degenerate targets
  # for k > 1; Lorenzo-Seva & Van Ginkel (2016) sidestep this by running a
  # Promin target rotation on top of the centroid, which this engine does
  # not implement. NCOL() treats NULL/vectors/1-D arrays as a single column
  # so the guard fires only on genuinely multi-column inputs.
  if (rotation == "oblique" && is.list(unrotated_list) &&
      length(unrotated_list) >= 1L &&
      isTRUE(NCOL(unrotated_list[[1L]]) >= 2L)) {
    cli::cli_abort(
      c("GPA-consensus alignment does not support oblique rotations with more than one factor.",
        "i" = "Use {.code target_method = \"first_target\"} in {.fn EFA_POOLED}, or pass orthogonal unrotated loadings."),
      class = "efa_consensus_oblique_unsupported"
    )
  }

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
    cli::cli_abort("No multi-start consensus run produced a finite final loss.",
                   class = "efa_consensus_no_finite_loss")
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
