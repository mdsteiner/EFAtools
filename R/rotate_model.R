# Single rotation core: pick the engine by name, then post-process once.
#
# `.rotate_model()` resolves the `type` preset for the requested rotation, guards
# the single-factor case, runs the rotation engine that the rotation name maps to
# in a lookup table, and (for every engine except promax) hands the raw solution
# to `.reflect_and_order()` for the shared sign reflection, factor ordering,
# structure matrix, and explained variances. Selecting the engine from a name
# table keeps the rotation name and the criterion it runs in lockstep. Promax
# reflects and reorders its varimax base before the oblique fit, so it finalizes
# its own solution.

# Orthogonal rotation engines, keyed by rotation name. Each resolves its criterion-specific
# argument and hands the matching compiled entry to `.gpf_native()`, which returns the raw
# rotated solution ($loadings, $Th); `eps`, `normalize`, `randomStarts`, and any extra
# arguments are forwarded through `...`. The Crawford-Ferguson criteria (quartimax with
# kappa = 0, equamax with kappa = k / (2 p)), geomin, Bentler, and bifactor all use the native
# engine. The geomin offset `delta` defaults to GPArotation's 0.01 and may be overridden by an
# exact-named `delta` in `...`.
.orth_engines <- list(
  equamax   = function(L, ...) .gpf_native(.rotate_cf_orth, L,
                                           list(kappa = ncol(L) / (2 * nrow(L))), ...),
  quartimax = function(L, ...) .gpf_native(.rotate_cf_orth, L, list(kappa = 0), ...),
  bentlerT  = function(L, ...) .gpf_native(.rotate_bentler_orth, L, list(), ...),
  geominT   = function(L, ...) .gpf_native(.rotate_geomin_orth, L,
                                           list(delta = .gpf_crit(list(...), "delta", 0.01)),
                                           ...),
  bifactorT = function(L, ...) .gpf_native(.rotate_bifactor_orth, L, list(), ...)
)

# Warn when a native gradient-projection rotation did not reach the convergence tolerance
# within `maxit` iterations, mirroring the warning the GPArotation engines emit so that the
# native and GPArotation-backed rotations behave the same way on a non-convergent fit.
.warn_rotation_no_convergence <- function(converged, maxit) {
  if (!isTRUE(converged)) {
    cli::cli_warn("The rotation did not converge within {maxit} iterations.",
                  class = "efa_rotation_no_convergence")
  }
}

# Read an overridable criterion parameter (the geomin `delta`, the oblimin `gam`) from a
# rotation engine's `...` by EXACT name, falling back to `default`. Looking the value up by
# exact name -- rather than declaring the parameter as a named formal before `...` -- keeps
# R's partial matching from silently binding an abbreviated or misspelled dot argument (e.g.
# `EFA(rotation = "geominT", del = 1)`) to the criterion parameter; an unrecognised argument is
# ignored, as it already is for the criteria that take no tuning parameter.
.gpf_crit <- function(dots, name, default) {
  if (name %in% names(dots)) dots[[name]] else default
}

# Shared driver for the native gradient-projection rotations. `entry` is the compiled rotation
# (e.g. `.rotate_cf_orth`); `crit` is a named list of the criterion-specific argument(s) for
# that entry (empty for the criteria that take none). The GPArotation-style controls (`eps`,
# `normalize`, `randomStarts` -> `random_starts`, `maxit`) are mapped onto the compiled entry, a
# non-convergent fit warns, and the raw rotated solution is returned in the convention
# `.reflect_and_order()` expects: the rotated `$loadings`, the rotation matrix `$Th`, and -- for
# oblique rotations -- the factor correlations `$Phi` (orthogonal entries return no `Phi`).
# Passing the criterion arguments explicitly via `crit` keeps them out of the partial-matching
# path; the trailing `...` simply absorbs the controls the engine-table closures forward
# through. The screen-and-triage multistart controls (`screen_keep`, `triage_maxit`,
# `triage_improve_tol`) are forwarded only when a criterion overrides the compiled entry's default
# (a `NULL` leaves that default in place); they are set per criterion by `.gpf_multistart_defaults`
# below. simplimax's entry takes no such controls (it always runs full multistart), and its closure
# passes none, so they are never forwarded to it.
.gpf_native <- function(entry, L, crit = list(), eps = 1e-5, normalize = TRUE,
                        randomStarts = 100, maxit = 1000L,
                        screen_keep = NULL, triage_maxit = NULL,
                        triage_improve_tol = NULL, ...) {
  args <- list(eps = eps, normalize = normalize,
               random_starts = randomStarts, maxit = maxit)
  if (!is.null(screen_keep)) args$screen_keep <- screen_keep
  if (!is.null(triage_maxit)) args$triage_maxit <- triage_maxit
  if (!is.null(triage_improve_tol)) args$triage_improve_tol <- triage_improve_tol
  res <- do.call(entry, c(list(unclass(L)), crit, args))
  .warn_rotation_no_convergence(res$convergence, maxit)
  out <- list(loadings = res$loadings, Th = res$Th,
              all_values = res$all_values, all_converged = res$all_converged)
  if (!is.null(res$Phi)) out$Phi <- res$Phi
  out
}

# Oblique rotation engines, keyed by rotation name. Each resolves its criterion-specific
# argument and hands the matching compiled entry to `.gpf_native()`, which returns the raw
# rotated solution ($loadings, $Th, $Phi); `k` is the simplimax target count (consumed here and
# ignored by the others). The oblimin and quartimin criteria (quartimin is oblimin pinned at
# `gam = 0`, not overridable; oblimin's `gam` defaults to 0 and may be overridden by an
# exact-named `gam` in `...`), geomin, Bentler, bifactor, and simplimax all use the native
# engine.
.oblq_engines <- list(
  oblimin   = function(L, k, ...) .gpf_native(.rotate_oblimin, L,
                                              list(gam = .gpf_crit(list(...), "gam", 0)), ...),
  quartimin = function(L, k, ...) .gpf_native(.rotate_oblimin, L, list(gam = 0), ...),
  simplimax = function(L, k, ...) .gpf_native(.rotate_simplimax_oblq, L, list(k = k), ...),
  bentlerQ  = function(L, k, ...) .gpf_native(.rotate_bentler_oblq, L, list(), ...),
  geominQ   = function(L, k, ...) .gpf_native(.rotate_geomin_oblq, L,
                                              list(delta = .gpf_crit(list(...), "delta", 0.01)),
                                              ...),
  bifactorQ = function(L, k, ...) .gpf_native(.rotate_bifactor_oblq, L, list(), ...)
)

# Per-criterion screen-and-triage multistart overrides for the native gradient-projection
# rotations. The compiled entries triage the five best-screened random starts by default
# (screen_keep = 5, triage_maxit = 25); two already reach the global criterion optimum for every
# criterion except geominQ across a panel of real correlation matrices, and five is a cheap safety
# margin (a small, single-digit-millisecond cost
# per rotation) against harder, unseen data. The oblique geomin criterion (geominQ) is multimodal
# enough to need more: its best-screened starts can all miss the global basin, so it triages more
# starts for longer (screen_keep = 10, triage_maxit = 50), recovering the global optimum on every
# tested matrix while staying far cheaper than fully optimizing all starts. A criterion absent from
# this list uses the compiled default. (geominT, the orthogonal geomin, reaches the global optimum
# at the default; simplimax is omitted because it always runs full multistart, not
# screen-and-triage.)
.gpf_multistart_defaults <- list(
  geominQ = list(screen_keep = 10, triage_maxit = 50)
)

# Canonical rotation names by family. The engine tables above key the GPArotation
# engines and exclude the special-cased `varimax`/`promax`; these vectors are the
# *full* family memberships used to classify a rotation name. `.rotation_family()`
# maps a single rotation name to "none", "orthogonal", or "oblique".
.orth_rotations <- c("varimax", "quartimax", "equamax", "bentlerT", "geominT", "bifactorT")
.oblq_rotations <- c("promax", "oblimin", "quartimin", "simplimax",
                     "bentlerQ", "geominQ", "bifactorQ")

.rotation_family <- function(rotation) {
  if (identical(rotation, "none")) return("none")
  if (rotation %in% .orth_rotations) return("orthogonal")
  if (rotation %in% .oblq_rotations) return("oblique")
  cli::cli_abort("Unknown rotation: {.val {rotation}}.", class = "efa_unknown_rotation")
}

# Shared post-processing for a rotated solution. Reflects each factor to a
# consistent (positive) orientation, orders the factors by their sum of squared
# loadings, and assembles the rotated loadings, factor correlations, structure
# matrix, and explained variances. The sign reflection and the factor ordering are
# applied to the loadings, the rotation matrix, and the factor intercorrelations
# alike, so the structure matrix and reported correlations stay consistent with the
# loadings. `name_factors` controls whether the factor columns are labelled
# (varimax labels them; the GPArotation engines leave them unnamed). Used by every
# engine except promax.
.reflect_and_order <- function(loadings, Phi = NULL, rotmat, L_unrot,
                               name_factors) {

  oblique <- !is.null(Phi)
  var_names <- rownames(L_unrot)

  # reflect factors with negative sums
  signs <- .reflect_signs(loadings)
  loadings <- loadings %*% diag(signs, nrow = length(signs))

  # order the factors by their sum of squared loadings (largest first)
  ord <- order(colSums(loadings^2), decreasing = TRUE)
  loadings <- loadings[, ord, drop = FALSE]

  fac_names <- if (isTRUE(name_factors)) colnames(L_unrot)[ord] else NULL
  dimnames(loadings) <- list(var_names, fac_names)

  # the rotation matrix follows the same sign reflection and factor ordering as
  # the loadings, so the documented reproduction identity holds (orthogonal:
  # L_unrot %*% rotmat; oblique: L_unrot %*% t(solve(rotmat)))
  rotmat <- (rotmat %*% diag(signs))[, ord, drop = FALSE]

  if (oblique) {
    # reflect and reorder the factor intercorrelations the same way as the
    # loadings so the structure matrix and reported correlations stay in sync
    Phi <- diag(signs) %*% Phi %*% diag(signs)
    Phi <- Phi[ord, ord]
    dimnames(Phi) <- NULL
  }

  vars_accounted_rot <- .compute_vars(L_unrot = L_unrot, L_rot = loadings, Phi = Phi)
  colnames(vars_accounted_rot) <- fac_names

  if (!oblique) {
    class(loadings) <- "LOADINGS"
    return(list(rot_loadings = loadings,
                rotmat = rotmat,
                vars_accounted_rot = vars_accounted_rot))
  }

  Structure <- loadings %*% Phi
  dimnames(Structure) <- list(var_names, fac_names)

  class(loadings) <- "LOADINGS"
  class(Structure) <- "LOADINGS"

  list(rot_loadings = loadings,
       Phi = Phi,
       Structure = Structure,
       rotmat = rotmat,
       vars_accounted_rot = vars_accounted_rot)
}

# Single-factor solutions cannot be rotated; warn and return the unrotated
# loadings in the family's output shape. The rotation-only outputs (factor
# intercorrelations, structure matrix, rotated variances) are returned as NULL so
# the print/summary sections that guard on `is.null` skip them.
.rotate_single_factor <- function(L, settings, oblique) {

  cli::cli_warn("A single factor cannot be rotated; returning the unrotated loadings.",
                class = "efa_single_factor")

  if (isTRUE(oblique)) {
    return(list(rot_loadings = L, Phi = NULL, Structure = NULL, rotmat = NA,
                vars_accounted_rot = NULL, settings = settings))
  }

  list(rot_loadings = L, rotmat = NA, vars_accounted_rot = NULL, settings = settings)
}

# Count distinct local minima among the criterion values reached by the rotation's
# random starts. Values within a small relative tolerance of each other are treated
# as the same minimum (different starts converging to one optimum differ only at the
# convergence tolerance), so the count reflects how many genuinely different optima
# the starts found.
.count_distinct_minima <- function(values, tol = 1e-6) {
  v <- sort(values[is.finite(values)])
  if (length(v) == 0L) return(0L)
  n <- 1L
  for (i in seq_along(v)[-1]) {
    if (v[i] - v[i - 1L] > tol * (1 + abs(v[i]))) n <- n + 1L
  }
  n
}

# Summarise a multistart rotation run for the output: the number of random starts
# requested, how many converged to a genuine local optimum, how many distinct optima they
# found, and the spread and best value of the attained criterion. `all_values` and
# `all_converged` are the per-start criterion values and convergence flags returned by the
# native rotation engines. The distinct-optima count and spread are taken over the
# CONVERGED starts only: the screen-and-triage engine leaves the starts it does not promote
# at a short, unconverged triage iterate, which is not a local optimum and must not be
# counted as one (those starts have `all_converged == 0`). The best criterion value is the
# lowest attained over all finite starts (the selected solution).
.rotation_diagnostics <- function(all_values, all_converged, randomStarts) {
  finite <- is.finite(all_values)
  converged <- all_values[finite & all_converged == 1L]
  list(n_starts = randomStarts,
       n_converged = length(converged),
       n_distinct_minima = .count_distinct_minima(converged),
       criterion_spread = if (length(converged) > 1L) max(converged) - min(converged) else 0,
       criterion_best = if (any(finite)) min(all_values[finite]) else NA_real_)
}

# Resolve the preset, run the named engine, and post-process. The field set and
# order of the returned list differ per rotation family.
.rotate_model <- function(x, rotation, type = c("EFAtools", "psych", "SPSS", "none"),
                          normalize = TRUE, precision = 1e-5, order_type = NA,
                          varimax_type = NA, P_type = NA, k = NA,
                          randomStarts = 100, ...) {

  L <- x$unrot_loadings

  if (rotation == "varimax") {

    resolved <- .resolve_settings(
      type = type,
      user = list(normalize = normalize, order_type = order_type,
                  varimax_type = varimax_type),
      preset = .efa_presets$VARIMAX
    )
    settings <- list(normalize = resolved$normalize, precision = precision,
                     order_type = resolved$order_type,
                     varimax_type = resolved$varimax_type)

    if (ncol(L) < 2) return(.rotate_single_factor(L, settings, oblique = FALSE))

    AV <- if (resolved$varimax_type == "svd") {
      stats::varimax(L, normalize = resolved$normalize, eps = precision)
    } else {
      .VARIMAX_SPSS(L, normalize = resolved$normalize, precision = precision)
    }
    out <- .reflect_and_order(AV$loadings, rotmat = AV$rotmat, L_unrot = L,
                              name_factors = TRUE)
    return(c(out, list(settings = settings)))

  }

  if (rotation == "promax") {

    resolved <- .resolve_settings(
      type = type,
      user = list(normalize = normalize, P_type = P_type, order_type = order_type,
                  varimax_type = varimax_type, k = k),
      preset = .efa_presets$PROMAX
    )
    settings <- list(normalize = resolved$normalize, P_type = resolved$P_type,
                     precision = precision, order_type = resolved$order_type,
                     varimax_type = resolved$varimax_type, k = resolved$k)

    if (ncol(L) < 2) return(.rotate_single_factor(L, settings, oblique = TRUE))

    out <- .rotate_promax(L, normalize = resolved$normalize, P_type = resolved$P_type,
                          precision = precision, order_type = resolved$order_type,
                          varimax_type = resolved$varimax_type, k = resolved$k)
    return(c(out, list(settings = settings)))

  }

  if (rotation %in% names(.orth_engines)) {

    resolved <- .resolve_settings(
      type = type,
      user = list(normalize = normalize, order_type = order_type),
      preset = .efa_presets$ROTATE_ORTH
    )
    settings <- list(normalize = resolved$normalize, precision = precision,
                     order_type = resolved$order_type, randomStarts = randomStarts)

    if (ncol(L) < 2) return(.rotate_single_factor(L, settings, oblique = FALSE))

    ms <- .gpf_multistart_defaults[[rotation]]
    AV <- .orth_engines[[rotation]](L, eps = precision,
                                    normalize = resolved$normalize,
                                    randomStarts = randomStarts,
                                    screen_keep = ms$screen_keep,
                                    triage_maxit = ms$triage_maxit,
                                    triage_improve_tol = ms$triage_improve_tol, ...)
    settings$rotation_diagnostics <- .rotation_diagnostics(AV$all_values,
                                                           AV$all_converged, randomStarts)
    out <- .reflect_and_order(AV$loadings, rotmat = AV$Th, L_unrot = L,
                              name_factors = FALSE)
    return(c(out, list(settings = settings)))

  }

  if (rotation %in% names(.oblq_engines)) {

    resolved <- .resolve_settings(
      type = type,
      user = list(normalize = normalize, order_type = order_type),
      preset = .efa_presets$ROTATE_OBLQ
    )
    settings <- list(normalize = resolved$normalize, precision = precision,
                     order_type = resolved$order_type, k = k,
                     randomStarts = randomStarts)

    if (ncol(L) < 2) return(.rotate_single_factor(L, settings, oblique = TRUE))

    if (rotation == "simplimax" && is.na(k)) {
      k <- nrow(L)
      settings$k <- k
    }

    ms <- .gpf_multistart_defaults[[rotation]]
    AV <- .oblq_engines[[rotation]](L, k = k, eps = precision,
                                    normalize = resolved$normalize,
                                    randomStarts = randomStarts,
                                    screen_keep = ms$screen_keep,
                                    triage_maxit = ms$triage_maxit,
                                    triage_improve_tol = ms$triage_improve_tol, ...)
    settings$rotation_diagnostics <- .rotation_diagnostics(AV$all_values,
                                                           AV$all_converged, randomStarts)
    out <- .reflect_and_order(AV$loadings, Phi = AV$Phi, rotmat = AV$Th,
                              L_unrot = L, name_factors = FALSE)
    return(c(out, list(settings = settings)))

  }
}
