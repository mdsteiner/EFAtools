# Procrustes rotation and the consensus-target engine: input validation, the closed-form
# orthogonal Procrustes solution, congruence/hyperplane diagnostics, and the single-start
# consensus iterator behind CONSENSUS_PROCRUSTES().

#' Average a list of matrices elementwise
#'
#' @param x List of conformable matrices.
#'
#' @returns A matrix with the same dimensions as the inputs.
#'
.average_matrices <- function(x) {
  if (!is.list(x) || length(x) < 1L) {
    cli::cli_abort("{.arg x} must be a non-empty list of conformable matrices.",
                   class = "efa_not_matrix_list")
  }
  out <- Reduce(`+`, x) / length(x)
  dimnames(out) <- dimnames(x[[1L]])
  out
}


# Internal scalar checks for Procrustes routines
#
# These helpers intentionally use base R rather than checkmate. They are used
# by both the public wrappers and the internal consensus engine, so keeping the
# checks local gives clearer error messages and avoids partially entering the
# compiled optimizer with invalid control values.
.procrustes_is_scalar <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

.procrustes_check_numeric_scalar <- function(x, name, lower = -Inf, upper = Inf,
                                             lower_closed = TRUE,
                                             upper_closed = TRUE) {
  if (!.procrustes_is_scalar(x)) {
    cli::cli_abort("{.arg {name}} must be a finite numeric scalar.",
                   class = "efa_not_scalar")
  }

  lower_ok <- if (lower_closed) x >= lower else x > lower
  upper_ok <- if (upper_closed) x <= upper else x < upper

  if (!lower_ok || !upper_ok) {
    left <- if (lower_closed) "[" else "("
    right <- if (upper_closed) "]" else ")"
    cli::cli_abort("{.arg {name}} must be in {left}{lower}, {upper}{right}.",
                   class = "efa_out_of_range")
  }

  x
}

.procrustes_check_integer_scalar <- function(x, name, lower = -Inf, upper = Inf) {
  if (!.procrustes_is_scalar(x) || abs(x - round(x)) > sqrt(.Machine$double.eps)) {
    cli::cli_abort("{.arg {name}} must be a finite integer-like scalar.",
                   class = "efa_not_integer")
  }
  x <- as.integer(round(x))
  if (x < lower || x > upper) {
    cli::cli_abort("{.arg {name}} must be an integer in [{lower}, {upper}].",
                   class = "efa_out_of_range")
  }
  x
}

.procrustes_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    cli::cli_abort("{.arg {name}} must be {.code TRUE} or {.code FALSE}.",
                   class = "efa_not_flag")
  }
  x
}

.procrustes_as_matrix <- function(x, name) {
  if (inherits(x, "LOADINGS")) {
    x <- unclass(x)
  }

  x <- as.matrix(x)
  if (!is.numeric(x)) {
    cli::cli_abort("{.arg {name}} must be a numeric matrix.", class = "efa_not_matrix")
  }
  storage.mode(x) <- "double"

  if (length(dim(x)) != 2L || any(dim(x) < 1L)) {
    cli::cli_abort("{.arg {name}} must have at least one row and one column.",
                   class = "efa_empty_matrix")
  }
  if (any(!is.finite(x))) {
    cli::cli_abort("{.arg {name}} must contain only finite values.",
                   class = "efa_nonfinite")
  }

  x
}

.procrustes_validate_matrix_pair <- function(A, B) {
  A <- .procrustes_as_matrix(A, "A")
  B <- .procrustes_as_matrix(B, "Target")

  if (!identical(dim(A), dim(B))) {
    cli::cli_abort("{.arg A} and {.arg Target} must have identical dimensions.",
                   class = "efa_dim_mismatch")
  }

  list(A = A, B = B)
}

.procrustes_validate_matrix_list <- function(x, name, min_length = 1L,
                                             expected_dim = NULL) {
  if (!is.list(x) || length(x) < min_length) {
    cli::cli_abort("{.arg {name}} must be a list with at least {min_length} matri{?x/ces}.",
                   class = "efa_not_matrix_list")
  }

  out <- lapply(seq_along(x), function(i) {
    .procrustes_as_matrix(x[[i]], paste0(name, "[[", i, "]]"))
  })

  ref_dim <- if (is.null(expected_dim)) dim(out[[1L]]) else expected_dim
  ok_dims <- vapply(out, function(z) identical(dim(z), ref_dim), logical(1L))
  if (!all(ok_dims)) {
    bad <- paste(which(!ok_dims), collapse = ", ")
    cli::cli_abort(
      c("All matrices in {.arg {name}} must have dimensions {paste(ref_dim, collapse = ' x ')}.",
        "x" = "Offending element(s): {bad}."),
      class = "efa_dim_mismatch"
    )
  }

  out
}

.procrustes_default_names <- function(n, prefix) {
  paste0(prefix, seq_len(n))
}

.procrustes_names <- function(A, B) {
  p <- nrow(A)
  k <- ncol(A)

  item_names <- rownames(B)
  if (is.null(item_names)) item_names <- rownames(A)
  if (is.null(item_names)) item_names <- .procrustes_default_names(p, "V")

  target_factor_names <- colnames(B)
  if (is.null(target_factor_names)) target_factor_names <- colnames(A)
  if (is.null(target_factor_names)) target_factor_names <- .procrustes_default_names(k, "F")

  source_factor_names <- colnames(A)
  if (is.null(source_factor_names)) source_factor_names <- target_factor_names

  list(
    items = item_names,
    source_factors = source_factor_names,
    target_factors = target_factor_names
  )
}

.procrustes_apply_dimnames <- function(out, A, B) {
  dn <- .procrustes_names(A, B)

  if (!is.null(out$loadings)) {
    dimnames(out$loadings) <- list(dn$items, dn$target_factors)
  }
  if (!is.null(out$T)) {
    dimnames(out$T) <- list(dn$source_factors, dn$target_factors)
  }
  if (!is.null(out$Phi)) {
    dimnames(out$Phi) <- list(dn$target_factors, dn$target_factors)
  }

  if (!is.null(out$all_values) && is.null(names(out$all_values))) {
    start_names <- if (!is.null(out$all_start_indices)) {
      paste0("start_", out$all_start_indices)
    } else {
      paste0("start_", seq_along(out$all_values))
    }
    names(out$all_values) <- start_names
  }
  if (!is.null(out$all_converged) && is.null(names(out$all_converged)) &&
      !is.null(names(out$all_values))) {
    names(out$all_converged) <- names(out$all_values)
  }
  if (!is.null(out$all_iterations) && is.null(names(out$all_iterations)) &&
      !is.null(names(out$all_values))) {
    names(out$all_iterations) <- names(out$all_values)
  }

  out
}

.procrustes_validate_crossprod <- function(S, k) {
  S <- .procrustes_as_matrix(S, "S")
  if (!identical(dim(S), c(k, k))) {
    cli::cli_abort("{.arg S} must be a {k} x {k} matrix.", class = "efa_dim_mismatch")
  }
  if (max(abs(S - t(S))) > 1e-8 * max(1, max(abs(S)))) {
    cli::cli_warn("{.arg S} is not symmetric; it should normally be {.code crossprod(A)}.",
                  class = "efa_not_symmetric")
  }
  S
}

.procrustes_validate_t_init <- function(T_init, k) {
  T_init <- .procrustes_as_matrix(T_init, "T_init")
  if (!identical(dim(T_init), c(k, k))) {
    cli::cli_abort("{.arg T_init} must be a {k} x {k} matrix.", class = "efa_dim_mismatch")
  }
  if (any(sqrt(colSums(T_init * T_init)) <= sqrt(.Machine$double.eps))) {
    cli::cli_abort("{.arg T_init} must not contain zero or near-zero columns.",
                   class = "efa_zero_column")
  }
  if (qr(T_init)$rank < k) {
    cli::cli_abort("{.arg T_init} must be nonsingular.", class = "efa_singular")
  }
  T_init
}

.procrustes_validate_oblique_controls <- function(oblique_eps,
                                                  oblique_maxit,
                                                  oblique_max_line_search,
                                                  oblique_step0,
                                                  oblique_normalize,
                                                  oblique_random_starts,
                                                  oblique_screen_keep,
                                                  oblique_triage_maxit,
                                                  oblique_triage_improve_tol) {
  list(
    oblique_eps = .procrustes_check_numeric_scalar(
      oblique_eps, "oblique_eps", lower = 0, lower_closed = FALSE
    ),
    oblique_maxit = .procrustes_check_integer_scalar(
      oblique_maxit, "oblique_maxit", lower = 0L
    ),
    oblique_max_line_search = .procrustes_check_integer_scalar(
      oblique_max_line_search, "oblique_max_line_search", lower = 0L
    ),
    oblique_step0 = .procrustes_check_numeric_scalar(
      oblique_step0, "oblique_step0", lower = 0, lower_closed = FALSE
    ),
    oblique_normalize = .procrustes_check_flag(oblique_normalize, "oblique_normalize"),
    oblique_random_starts = .procrustes_check_integer_scalar(
      oblique_random_starts, "oblique_random_starts", lower = 0L
    ),
    oblique_screen_keep = .procrustes_check_integer_scalar(
      oblique_screen_keep, "oblique_screen_keep", lower = 0L
    ),
    oblique_triage_maxit = .procrustes_check_integer_scalar(
      oblique_triage_maxit, "oblique_triage_maxit", lower = 0L
    ),
    oblique_triage_improve_tol = .procrustes_check_numeric_scalar(
      oblique_triage_improve_tol, "oblique_triage_improve_tol", lower = 0
    )
  )
}

.procrustes_validate_start_index <- function(start, n_targets) {
  if (!.procrustes_is_scalar(start) || abs(start - round(start)) > sqrt(.Machine$double.eps)) {
    cli::cli_abort("{.arg start} must be either an integer index or a target matrix.",
                   class = "efa_bad_start")
  }
  start <- as.integer(round(start))
  if (start < 1L || start > n_targets) {
    cli::cli_abort("{.arg start} must be between 1 and {n_targets}.",
                   class = "efa_out_of_range")
  }
  start
}

.procrustes_validate_starts <- function(starts, n_targets) {
  if (is.null(starts)) {
    return(seq_len(n_targets))
  }
  if (!is.numeric(starts) || length(starts) < 1L || anyNA(starts) ||
      any(!is.finite(starts)) || any(abs(starts - round(starts)) > sqrt(.Machine$double.eps))) {
    cli::cli_abort("{.arg starts} must be a non-empty integer vector.",
                   class = "efa_bad_start")
  }
  starts <- as.integer(round(starts))
  if (any(starts < 1L | starts > n_targets)) {
    cli::cli_abort("All entries in {.arg starts} must be between 1 and {n_targets}.",
                   class = "efa_out_of_range")
  }
  if (anyDuplicated(starts)) {
    cli::cli_warn("Duplicate entries in {.arg starts} were removed.",
                  class = "efa_duplicate_starts")
    starts <- unique(starts)
  }
  starts
}

#' Count near-zero loadings
#'
#' Hyperplane count is the number of loadings with absolute value smaller than a
#' user-specified cutoff.
#'
#' @param L Numeric loading matrix.
#' @param cutoff Numeric scalar. Loadings with `abs(L) < cutoff` are counted as
#'   being in the hyperplane.
#'
#' @returns A list with the total hyperplane count and counts by factor and item.
.hyperplane_count <- function(L, cutoff = 0.15) {
  L <- .procrustes_as_matrix(L, "L")
  cutoff <- .procrustes_check_numeric_scalar(cutoff, "cutoff", lower = 0)

  in_hyperplane <- abs(L) < cutoff
  list(
    total = sum(in_hyperplane),
    by_factor = colSums(in_hyperplane),
    by_item = rowSums(in_hyperplane),
    cutoff = cutoff
  )
}

#' Mean squared discrepancy to a consensus target
#'
#' @param aligned_loadings List of aligned loading matrices.
#' @param target Consensus target matrix.
#'
#' @returns Mean sum of squared deviations from the target across matrices.
#'
.consensus_loss <- function(aligned_loadings, target) {
  # Inputs are validated by the public/internal Procrustes entry points. Keeping
  # this inner-loop helper allocation-light matters for consensus and bootstrap
  # workflows where the loss is evaluated many times.
  mean(vapply(aligned_loadings, function(L) {
    sum((L - target)^2)
  }, numeric(1L)))
}


#' Closed-form orthogonal Procrustes rotation
#'
#' Rotate `A` to the orthogonal target `B` by minimizing
#' `||A %*% T - B||_F^2` subject to `t(T) %*% T = I`.
#'
#' @param A Numeric matrix to be rotated.
#' @param B Numeric target matrix with the same dimensions as `A`.
#'
#' @returns A list with the rotated loadings, orthogonal transformation matrix,
#'   target criterion value, and basic diagnostics.
#'
#' @references
#' Schoenemann, P. H. (1966). A generalized solution of the orthogonal
#' Procrustes problem. *Psychometrika*, 31, 1-10.
.orthogonal_procrustes <- function(A, B) {
  mats <- .procrustes_validate_matrix_pair(A, B)
  A <- mats$A
  B <- mats$B

  M <- crossprod(A, B)
  s <- svd(M)
  Tmat <- s$u %*% t(s$v)

  Lrot <- A %*% Tmat
  value <- 0.5 * sum((Lrot - B)^2)

  out <- list(
    loadings = Lrot,
    T = Tmat,
    Phi = diag(ncol(A)),
    value = value,
    convergence = TRUE,
    iterations = 0L,
    kappa_T = 1,
    Table = cbind(iter = 0, f = value, log10_s = NA_real_, step = NA_real_),
    method = "orthogonal_procrustes",
    best_start_index = 1L,
    all_start_indices = 1L,
    all_values = value,
    all_converged = TRUE,
    all_iterations = 0L,
    line_search_failed = FALSE
  )

  .procrustes_apply_dimnames(out, A, B)
}


#' Internal single-start consensus engine
#'
#' This internal helper performs a single consensus-target run from one starting
#' target. Users should normally call [CONSENSUS_PROCRUSTES()] instead.
#'
#' @inheritParams CONSENSUS_PROCRUSTES
#'
.consensus_target_procrustes_single <- function(unrotated_list,
                                                init_targets = NULL,
                                                rotation = c("orthogonal", "oblique"),
                                                start = 1,
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

  alpha <- .procrustes_check_numeric_scalar(alpha, "alpha", lower = 0, upper = 1,
                                            lower_closed = FALSE)
  tol <- .procrustes_check_numeric_scalar(tol, "tol", lower = 0, lower_closed = FALSE)
  if (!is.null(loss_tol)) {
    loss_tol <- .procrustes_check_numeric_scalar(loss_tol, "loss_tol", lower = 0,
                                                 lower_closed = FALSE)
  }
  loss_patience <- .procrustes_check_integer_scalar(loss_patience, "loss_patience", lower = 1L)
  min_iter <- .procrustes_check_integer_scalar(min_iter, "min_iter", lower = 0L)
  max_iter <- .procrustes_check_integer_scalar(max_iter, "max_iter", lower = 1L)
  if (min_iter > max_iter) {
    cli::cli_abort("{.arg min_iter} must not exceed {.arg max_iter}.",
                   class = "efa_out_of_range")
  }
  if (is.null(loss_tol) && convergence %in% c("loss", "both")) {
    cli::cli_abort("{.arg loss_tol} must not be {.code NULL} when {.arg convergence} is {.val loss} or {.val both}.",
                   class = "efa_bad_loss_tol")
  }

  match_target <- .procrustes_check_flag(match_target, "match_target")
  verbose <- .procrustes_check_flag(verbose, "verbose")
  hyper_cutoff <- .procrustes_check_numeric_scalar(hyper_cutoff, "hyper_cutoff", lower = 0)

  oblique_controls <- .procrustes_validate_oblique_controls(
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

  unrotated_list <- .procrustes_validate_matrix_list(
    unrotated_list, "unrotated_list", min_length = 2L
  )

  m <- length(unrotated_list)
  p <- nrow(unrotated_list[[1L]])
  k <- ncol(unrotated_list[[1L]])

  if (is.null(init_targets)) {
    init_targets <- unrotated_list
  } else {
    init_targets <- .procrustes_validate_matrix_list(
      init_targets, "init_targets", min_length = 1L, expected_dim = c(p, k)
    )
  }

  if (is.numeric(start) && length(start) == 1L) {
    start_index <- .procrustes_validate_start_index(start, length(init_targets))
    target <- init_targets[[start_index]]
    start_label <- start_index
  } else {
    target <- .procrustes_as_matrix(start, "start")
    if (!identical(dim(target), c(p, k))) {
      cli::cli_abort("If {.arg start} is a matrix, it must have the same dimensions as the loadings.",
                     class = "efa_dim_mismatch")
    }
    start_index <- NA_integer_
    start_label <- "matrix"
  }

  if (!is.na(start_index)) {
    dimnames(target) <- dimnames(init_targets[[start_index]])
  }
  if (is.null(rownames(target)) || is.null(colnames(target))) {
    dn <- .procrustes_names(unrotated_list[[1L]], target)
    dimnames(target) <- list(dn$items, dn$target_factors)
  }

  S_list <- lapply(unrotated_list, crossprod)
  T_starts <- vector("list", m)
  if (rotation == "oblique" && k > 1L) {
    T_starts <- lapply(unrotated_list, function(A) .orthogonal_procrustes(A, target)$T)
  }

  outer_random_starts <- if (oblique_random_starts_stage %in% c("outer", "both")) {
    oblique_controls$oblique_random_starts
  } else {
    0L
  }

  final_random_starts <- if (oblique_random_starts_stage %in% c("final", "both")) {
    oblique_controls$oblique_random_starts
  } else {
    0L
  }

  hist_iter <- integer(max_iter)
  hist_rel_change <- numeric(max_iter)
  hist_loss <- numeric(max_iter)
  hist_rel_loss_change <- rep(NA_real_, max_iter)
  hist_inner_mean_value <- numeric(max_iter)
  hist_inner_failures <- integer(max_iter)
  hist_inner_max_kappa <- numeric(max_iter)
  hist_loss_stable_count <- integer(max_iter)
  hist_target_ok <- logical(max_iter)
  hist_loss_ok <- logical(max_iter)
  hist_stop_rule_met <- logical(max_iter)

  converged <- FALSE
  convergence_reason <- NA_character_
  prev_loss <- NA_real_
  loss_stable_count <- 0L
  n_history <- 0L

  for (iter in seq_len(max_iter)) {
    aligned <- vector("list", m)

    for (d in seq_len(m)) {
      aligned[[d]] <- PROCRUSTES(
        A = unrotated_list[[d]],
        Target = target,
        rotation = rotation,
        S = S_list[[d]],
        T_init = T_starts[[d]],
        oblique_eps = oblique_controls$oblique_eps,
        oblique_maxit = oblique_controls$oblique_maxit,
        oblique_max_line_search = oblique_controls$oblique_max_line_search,
        oblique_step0 = oblique_controls$oblique_step0,
        oblique_normalize = oblique_controls$oblique_normalize,
        oblique_random_starts = outer_random_starts,
        oblique_screen_keep = oblique_controls$oblique_screen_keep,
        oblique_triage_maxit = oblique_controls$oblique_triage_maxit,
        oblique_triage_improve_tol = oblique_controls$oblique_triage_improve_tol
      )
      if (rotation == "oblique" && k > 1L) {
        T_starts[[d]] <- aligned[[d]]$T
      }
    }

    aligned_loadings <- lapply(aligned, `[[`, "loadings")
    inner_values <- vapply(aligned, `[[`, numeric(1L), "value")
    inner_converged <- vapply(aligned, function(x) isTRUE(x$convergence), logical(1L))
    inner_kappa <- vapply(aligned, function(x) {
      if (is.null(x$kappa_T)) NA_real_ else x$kappa_T
    }, numeric(1L))

    centroid <- .average_matrices(aligned_loadings)
    if (isTRUE(match_target) && k > 1L) {
      centroid <- .align_solution(L = centroid, L_target = target)$loadings
      dimnames(centroid) <- dimnames(target)
    }

    new_target <- (1 - alpha) * target + alpha * centroid
    dimnames(new_target) <- dimnames(target)

    rel_change <- norm(new_target - target, type = "F") /
      (norm(target, type = "F") + 1e-12)

    current_loss <- .consensus_loss(aligned_loadings, new_target)

    rel_loss_change <- if (is.na(prev_loss)) {
      NA_real_
    } else {
      abs(prev_loss - current_loss) / (abs(prev_loss) + 1e-12)
    }

    if (!is.null(loss_tol) && !is.na(rel_loss_change) && rel_loss_change < loss_tol) {
      loss_stable_count <- loss_stable_count + 1L
    } else {
      loss_stable_count <- 0L
    }

    target_ok <- rel_change < tol
    loss_ok <- !is.null(loss_tol) && loss_stable_count >= loss_patience

    stop_rule_met <- switch(
      convergence,
      target = target_ok,
      loss = loss_ok,
      either = target_ok || loss_ok,
      both = target_ok && loss_ok
    )

    n_history <- iter
    hist_iter[iter] <- iter
    hist_rel_change[iter] <- rel_change
    hist_loss[iter] <- current_loss
    hist_rel_loss_change[iter] <- rel_loss_change
    hist_inner_mean_value[iter] <- mean(inner_values)
    hist_inner_failures[iter] <- sum(!inner_converged)
    hist_inner_max_kappa[iter] <- if (all(is.na(inner_kappa))) NA_real_ else max(inner_kappa, na.rm = TRUE)
    hist_loss_stable_count[iter] <- loss_stable_count
    hist_target_ok[iter] <- target_ok
    hist_loss_ok[iter] <- loss_ok
    hist_stop_rule_met[iter] <- stop_rule_met

    if (verbose) {
      cli::cli_inform(sprintf(
        paste0(
          "iter %d | rel_change = %.3e | loss = %.8f | ",
          "rel_loss_change = %s | stable = %d | inner_failures = %d"
        ),
        iter,
        rel_change,
        current_loss,
        ifelse(is.na(rel_loss_change), "NA", sprintf("%.3e", rel_loss_change)),
        loss_stable_count,
        sum(!inner_converged)
      ))
    }

    target <- new_target
    prev_loss <- current_loss

    if (iter >= min_iter && isTRUE(stop_rule_met)) {
      converged <- TRUE
      convergence_reason <- if (target_ok && loss_ok) {
        "target_and_loss"
      } else if (target_ok) {
        "target"
      } else if (loss_ok) {
        "loss"
      } else {
        convergence
      }
      break
    }
  }

  keep <- seq_len(n_history)
  history <- data.frame(
    iter = hist_iter[keep],
    rel_change = hist_rel_change[keep],
    loss = hist_loss[keep],
    rel_loss_change = hist_rel_loss_change[keep],
    inner_mean_value = hist_inner_mean_value[keep],
    inner_failures = hist_inner_failures[keep],
    inner_max_kappa = hist_inner_max_kappa[keep],
    loss_stable_count = hist_loss_stable_count[keep],
    target_ok = hist_target_ok[keep],
    loss_ok = hist_loss_ok[keep],
    stop_rule_met = hist_stop_rule_met[keep]
  )

  final_aligned <- vector("list", m)
  for (d in seq_len(m)) {
    final_aligned[[d]] <- PROCRUSTES(
      A = unrotated_list[[d]],
      Target = target,
      rotation = rotation,
      S = S_list[[d]],
      T_init = T_starts[[d]],
      oblique_eps = oblique_controls$oblique_eps,
      oblique_maxit = oblique_controls$oblique_maxit,
      oblique_max_line_search = oblique_controls$oblique_max_line_search,
      oblique_step0 = oblique_controls$oblique_step0,
      oblique_normalize = oblique_controls$oblique_normalize,
      oblique_random_starts = final_random_starts,
      oblique_screen_keep = oblique_controls$oblique_screen_keep,
      oblique_triage_maxit = oblique_controls$oblique_triage_maxit,
      oblique_triage_improve_tol = oblique_controls$oblique_triage_improve_tol
    )
  }

  final_aligned_loadings <- lapply(final_aligned, `[[`, "loadings")
  final_T <- lapply(final_aligned, `[[`, "T")
  final_Phi <- lapply(final_aligned, `[[`, "Phi")
  final_values <- vapply(final_aligned, `[[`, numeric(1L), "value")
  final_converged <- vapply(final_aligned, function(x) isTRUE(x$convergence), logical(1L))
  final_kappa <- vapply(final_aligned, function(x) {
    if (is.null(x$kappa_T)) NA_real_ else x$kappa_T
  }, numeric(1L))

  pooled_loadings <- .average_matrices(final_aligned_loadings)
  dimnames(pooled_loadings) <- dimnames(target)

  pooled_phi <- if (rotation == "orthogonal" || k == 1L) {
    diag(k)
  } else {
    .average_matrices(final_Phi)
  }
  pooled_phi <- (pooled_phi + t(pooled_phi)) / 2
  diag(pooled_phi) <- 1
  dimnames(pooled_phi) <- list(colnames(target), colnames(target))

  final_loss <- .consensus_loss(final_aligned_loadings, target)

  outer_converged <- converged
  final_inner_converged <- all(final_converged)
  overall_converged <- outer_converged && final_inner_converged
  reported_convergence_reason <- if (overall_converged) {
    convergence_reason
  } else if (!outer_converged) {
    convergence_reason
  } else {
    "final_inner_alignment"
  }

  list(
    converged = overall_converged,
    outer_converged = outer_converged,
    final_inner_converged = final_inner_converged,
    convergence_reason = reported_convergence_reason,
    outer_convergence_reason = convergence_reason,
    iterations = nrow(history),
    tol = tol,
    loss_tol = loss_tol,
    loss_patience = loss_patience,
    convergence = convergence,
    min_iter = min_iter,
    alpha = alpha,
    match_target = match_target,
    rotation = rotation,
    start = start_label,
    start_index = start_index,
    target = target,
    history = history,
    aligned_loadings = final_aligned_loadings,
    aligned_phi = final_Phi,
    transformations = final_T,
    pooled_loadings = pooled_loadings,
    pooled_phi = pooled_phi,
    mean_loss = final_loss,
    final_values = final_values,
    final_converged = final_converged,
    final_failures = which(!final_converged),
    final_kappa_T = final_kappa,
    hyperplane_target = .hyperplane_count(target, cutoff = hyper_cutoff),
    hyperplane_pooled = .hyperplane_count(pooled_loadings, cutoff = hyper_cutoff),
    oblique_random_starts_stage = oblique_random_starts_stage,
    oblique_random_starts_outer = outer_random_starts,
    oblique_random_starts_final = final_random_starts
  )
}
