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

# Warn when an oblique rotation returns near-collinear factors. A factor correlation
# above `tol` in absolute value means two rotated factors are barely distinguishable,
# which is the usual signature of extracting more factors than the data support (and,
# for oblimin, of a large `gam`). The solution is still returned -- an extreme
# correlation is a property of the fit rather than a failure of the rotation -- but it
# should not be interpreted without inspecting the factor intercorrelations. The 0.9
# cut-off follows the extreme-correlation diagnostic GPArotation reports from its own
# oblique rotations (Bernaards & Jennrich, 2005).
#
# The reported value is TRUNCATED, not rounded: a correlation of 0.9999 must not print as
# "1.00", which is the perfectly collinear state a rejected singular transformation would
# have produced, and a correlation of 0.905 must not print as "0.90", which reads as
# sitting at the cut-off rather than above it.
.warn_extreme_phi <- function(Phi, tol = 0.9) {

  if (is.null(Phi) || !is.matrix(Phi) || nrow(Phi) < 2L) return(invisible(NULL))

  off <- abs(Phi[upper.tri(Phi)])
  off <- off[is.finite(off)]
  if (length(off) == 0L || max(off) <= tol) return(invisible(NULL))

  max_r <- formatC(trunc(max(off) * 1000) / 1000, format = "f", digits = 3)
  cli::cli_warn(
    c("Extreme factor correlations (max |r| = {max_r}).",
      "i" = "This usually indicates factor over-extraction; inspect {.field Phi}."),
    class = "efa_rotation_extreme_phi"
  )
  invisible(NULL)
}

# Warn when a bifactor criterion is requested for a two-factor solution. The
# Jennrich-Bentler bifactor criterion sums the products of squared loadings over pairs of
# distinct GROUP factors (Jennrich & Bentler, 2011, 2012), so with a general factor and a
# single group factor there is no such pair: the criterion is identically zero, every
# rotation attains it, and the identity start wins the tie. The returned "rotated" loadings
# are then the unrotated ones, while the multistart diagnostics -- one distinct optimum, zero
# criterion spread -- read exactly like a well-identified solution. The solution is still
# returned, so the situation is signalled rather than turned into an error.
.warn_bifactor_trivial <- function(rotation, n_factors) {
  if (!rotation %in% c("bifactorT", "bifactorQ") || n_factors != 2L) {
    return(invisible(NULL))
  }
  cli::cli_warn(
    c("The bifactor criterion is identically zero with only one group factor.",
      "i" = "With two factors every rotation is optimal, so the returned loadings are the
             unrotated ones.",
      "i" = "Extract at least three factors for a bifactor rotation."),
    class = "efa_bifactor_trivial"
  )
  invisible(NULL)
}

# Validate the arguments the compiled rotation entries take, in R, before the entry is
# called. The compiled entries check the same quantities and stay in place as backstops for
# direct internal callers, but they raise unclassed exceptions whose messages name neither
# the function nor the argument, and a non-numeric value never reaches them at all (the
# numeric coercion at the boundary fails first, with a message naming only the C type).
# Checking here gives every one of them the `efa_rotation_arg` class and the argument name
# the user typed. `crit` carries only the criterion parameter of the entry actually being
# run, so a criterion is never held to another criterion's contract; `maxit` is checked for
# every criterion because every entry takes it.
.assert_rotation_args <- function(crit, maxit, L) {

  bad_arg <- function(arg, requirement) {
    cli::cli_abort("{.arg {arg}} must be {requirement}.", class = "efa_rotation_arg")
  }

  # A fractional `maxit` is truncated at the C boundary, so the engine would run a different
  # budget from the one the non-convergence warning reports.
  if (!isTRUE(checkmate::test_int(maxit, lower = 0))) {
    bad_arg("maxit", "a single whole number of at least 0")
  }

  if ("gam" %in% names(crit) && !isTRUE(checkmate::test_number(crit$gam, finite = TRUE))) {
    bad_arg("gam", "a single finite number")
  }

  if ("delta" %in% names(crit) &&
      !(isTRUE(checkmate::test_number(crit$delta, finite = TRUE)) && crit$delta > 0)) {
    bad_arg("delta", "a single positive finite number")
  }

  # simplimax counts the loadings it drives toward zero, so its upper bound is the number of
  # loadings in the solution; the message quotes that bound rather than leaving the user to
  # work it out.
  if ("k" %in% names(crit)) {
    n_loadings <- nrow(L) * ncol(L)
    if (!isTRUE(checkmate::test_int(crit$k, lower = 1, upper = n_loadings))) {
      bad_arg("k", paste0("a single whole number between 1 and ", n_loadings,
                          ", the number of loadings"))
    }
  }

  invisible(NULL)
}

# Read an overridable criterion parameter (the geomin `delta`, the oblimin `gam`) from a
# rotation engine's `...` by EXACT name, falling back to `default`. Looking the value up by
# exact name -- rather than declaring the parameter as a named formal before `...` -- keeps
# R's partial matching from silently binding an abbreviated or misspelled dot argument (e.g.
# `EFA(rotation = "geominT", del = 1)`) to the criterion parameter. efa_fit() and
# rotate_control() reject a name their engines cannot consume up front (see
# `.rotation_dot_extras()`); the exact-name lookup remains the guard against partial matching
# for direct internal callers, where an unrecognised argument is simply ignored.
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
  .assert_rotation_args(crit, maxit, L)
  args <- list(eps = eps, normalize = normalize,
               random_starts = randomStarts, maxit = maxit)
  if (!is.null(screen_keep)) args$screen_keep <- screen_keep
  if (!is.null(triage_maxit)) args$triage_maxit <- triage_maxit
  if (!is.null(triage_improve_tol)) args$triage_improve_tol <- triage_improve_tol
  res <- do.call(entry, c(list(unclass(L)), crit, args))
  .warn_rotation_no_convergence(res$convergence, maxit)
  # `converged` is the flag of the start whose solution is RETURNED, which the aggregate
  # per-start counts below cannot reconstruct: the winner is chosen by lowest criterion
  # value, so a non-converged start with a lower value beats a converged one.
  out <- list(loadings = res$loadings, Th = res$Th, converged = isTRUE(res$convergence),
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

# Names a rotation's engine can consume from efa_fit()'s `...`, kept in lockstep with the
# engine tables above: `maxit` is the one `.gpf_native()` control `.rotate_model()` does not
# already pass explicitly, and the criterion extras are read by exact name via `.gpf_crit()`
# (`gam` for oblimin, `delta` for the geomin pair). varimax and promax run outside the native
# engines and take no extras, as does rotation = "none", where no engine runs at all.
# efa_fit() rejects any other dot argument up front; `.rotation_extra_union` is the
# across-rotation union rotate_control() validates its stored extras against.
.rotation_dot_extras <- function(rotation) {
  switch(rotation,
         none = ,
         varimax = ,
         promax = character(0),
         oblimin = c("maxit", "gam"),
         geominT = ,
         geominQ = c("maxit", "delta"),
         "maxit")
}

.rotation_extra_union <- c("maxit", "gam", "delta")

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
# consistent (positive) orientation, orders the factors as requested by
# `order_type`, and assembles the rotated loadings, factor correlations, structure
# matrix, and explained variances. The sign reflection and the factor ordering are
# applied to the loadings, the rotation matrix, and the factor intercorrelations
# alike, so the structure matrix and reported correlations stay consistent with the
# loadings. The factor columns are labelled "F1".."Fk" by their position in the
# ordered solution, so every rotation returns labelled loadings, factor
# intercorrelations, and structure coefficients. The labels deliberately do not
# follow the unrotated columns: a rotation mixes all k factors, so a rotated column
# corresponds to no single unrotated factor, and the engines return their columns in
# an arbitrary order (for the multistart engines it depends on which random start
# wins). Carrying the unrotated labels through the sort would therefore attach a
# run-dependent permutation of "F1".."Fk" to an otherwise stable solution.
# Used by every engine except promax.
.reflect_and_order <- function(loadings, Phi = NULL, rotmat, L_unrot, order_type) {

  oblique <- !is.null(Phi)
  var_names <- rownames(L_unrot)

  # reflect factors with negative sums
  signs <- .reflect_signs(loadings)
  loadings <- loadings %*% diag(signs, nrow = length(signs))

  if (oblique) {
    # reflect the factor intercorrelations the same way as the loadings so the
    # structure matrix and reported correlations stay in sync
    Phi <- diag(signs, nrow = length(signs)) %*% Phi %*% diag(signs, nrow = length(signs))
  }

  # order the factors by their explained variance (largest first). "eigen" orders
  # by the reported "SS loadings": colSums(L^2) for orthogonal solutions, the
  # Phi-weighted sum of squares diag(Phi L'L) for oblique ones, so the reported
  # variances decrease monotonically (as in psych). "ss_factors" orders oblique
  # factors by the raw pattern sum of squares colSums(L^2) instead. For orthogonal
  # solutions Phi = I, so the two keys coincide and order_type has no effect.
  ss <- if (oblique && order_type == "eigen") {
    diag(Phi %*% crossprod(loadings))
  } else {
    colSums(loadings^2)
  }
  ord <- order(ss, decreasing = TRUE)
  loadings <- loadings[, ord, drop = FALSE]

  fac_names <- colnames(L_unrot)
  dimnames(loadings) <- list(var_names, fac_names)

  # the rotation matrix follows the same sign reflection and factor ordering as
  # the loadings, so the documented reproduction identity holds (orthogonal:
  # L_unrot %*% rotmat; oblique: L_unrot %*% t(solve(rotmat)))
  rotmat <- (rotmat %*% diag(signs))[, ord, drop = FALSE]

  if (oblique) {
    Phi <- Phi[ord, ord, drop = FALSE]
    dimnames(Phi) <- list(fac_names, fac_names)
  }

  vars_accounted_rot <- .compute_vars(L_unrot = L_unrot, L_rot = loadings, Phi = Phi)
  colnames(vars_accounted_rot) <- fac_names

  if (!oblique) {
    class(loadings) <- c("efa_loadings", "LOADINGS")
    return(list(rot_loadings = loadings,
                rotmat = rotmat,
                vars_accounted_rot = vars_accounted_rot))
  }

  Structure <- loadings %*% Phi
  dimnames(Structure) <- list(var_names, fac_names)

  class(loadings) <- c("efa_loadings", "LOADINGS")
  class(Structure) <- c("efa_loadings", "LOADINGS")

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

# Summarise a multistart rotation run for the output: how many starts were available in
# total, how many of them were actually optimized, how many converged to a genuine local
# optimum, how many distinct optima they found, and the spread and best value of the
# attained criterion. `all_values` and `all_converged` are the per-start criterion values
# and convergence flags returned by the native rotation engines, one entry per OPTIMIZED
# start, so `length(all_values)` is the number of optimizations that actually ran. Under the
# screen-and-triage strategy that is far fewer than the number of random starts: the starts
# that are not among the best-screened have their objective evaluated at the start point and
# are then discarded, so they can never yield a local optimum and are not counted here.
# `n_starts_total` includes the rational (identity) start alongside the `randomStarts` random
# ones, so it bounds `n_converged` from above. The distinct-optima count and spread are taken
# over the CONVERGED starts only: the screen-and-triage engine leaves the starts it does not
# promote at a short, unconverged triage iterate, which is not a local optimum and must not
# be counted as one (those starts have `all_converged == 0`). The best criterion value is the
# lowest attained over all finite starts (the selected solution).
# `converged` is the convergence flag of the start that WON, i.e. of the solution actually
# returned. It is reported separately from `n_converged` because the two can disagree: the
# winner is the lowest-criterion start, and a start that stopped short of the tolerance at a
# lower criterion value beats a converged start at a higher one. Collapsing the two would
# present such a solution as converged because its rivals were.
.rotation_diagnostics <- function(all_values, all_converged, randomStarts, converged) {
  finite <- is.finite(all_values)
  converged_values <- all_values[finite & all_converged == 1L]
  list(converged = isTRUE(converged),
       n_starts_total = as.integer(randomStarts) + 1L,
       n_optimized = length(all_values),
       n_converged = length(converged_values),
       n_distinct_minima = .count_distinct_minima(converged_values),
       criterion_spread = if (length(converged_values) > 1L) diff(range(converged_values)) else 0,
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
      .varimax_svd(L, normalize = resolved$normalize, precision = precision)
    } else {
      .VARIMAX_SPSS(L, normalize = resolved$normalize, precision = precision)
    }
    out <- .reflect_and_order(AV$loadings, rotmat = AV$rotmat, L_unrot = L,
                              order_type = resolved$order_type)
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
    .warn_extreme_phi(out$Phi)
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

    .warn_bifactor_trivial(rotation, ncol(L))

    ms <- .gpf_multistart_defaults[[rotation]]
    AV <- .orth_engines[[rotation]](L, eps = precision,
                                    normalize = resolved$normalize,
                                    randomStarts = randomStarts,
                                    screen_keep = ms$screen_keep,
                                    triage_maxit = ms$triage_maxit,
                                    triage_improve_tol = ms$triage_improve_tol, ...)
    settings$rotation_diagnostics <- .rotation_diagnostics(AV$all_values,
                                                           AV$all_converged, randomStarts,
                                                           AV$converged)
    out <- .reflect_and_order(AV$loadings, rotmat = AV$Th, L_unrot = L,
                              order_type = resolved$order_type)
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

    .warn_bifactor_trivial(rotation, ncol(L))

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
                                                           AV$all_converged, randomStarts,
                                                           AV$converged)
    out <- .reflect_and_order(AV$loadings, Phi = AV$Phi, rotmat = AV$Th,
                              L_unrot = L, order_type = resolved$order_type)
    .warn_extreme_phi(out$Phi)
    return(c(out, list(settings = settings)))

  }

  cli::cli_abort("Unknown rotation: {.val {rotation}}.", class = "efa_unknown_rotation")
}
