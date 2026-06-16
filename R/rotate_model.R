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

# Orthogonal rotation engines, keyed by rotation name. Each returns the raw rotated
# solution ($loadings, $Th); `eps`, `normalize`, `randomStarts`, and any extra
# arguments are forwarded through `...`. The Crawford-Ferguson criteria (quartimax,
# equamax) use the native engine; the remaining criteria use GPArotation.
.orth_engines <- list(
  equamax   = function(L, ...) .gpf_cf_orth(L, kappa = ncol(L) / (2 * nrow(L)), ...),
  quartimax = function(L, ...) .gpf_cf_orth(L, kappa = 0, ...),
  bentlerT  = function(L, ...) GPArotation::bentlerT(L, ...),
  geominT   = function(L, ...) GPArotation::geominT(L, ...),
  bifactorT = function(L, ...) GPArotation::bifactorT(L, ...)
)

# Thin wrapper over the native Crawford-Ferguson orthogonal rotation. Maps the
# GPArotation-style engine signature (`eps`, `normalize`, `randomStarts`, `maxit`) onto
# the compiled entry and returns the raw rotated solution ($loadings, $Th) in the
# convention `.reflect_and_order()` expects (the rotation matrix reproduces the rotated
# loadings via `L_unrot %*% Th`). `maxit` is forwarded so the documented EFA `...`
# rotation control reaches the optimizer as it does for the GPArotation engines; the
# criterion itself takes no further tuning arguments, so any other `...` are ignored.
.gpf_cf_orth <- function(L, kappa, eps = 1e-5, normalize = TRUE,
                         randomStarts = 100, maxit = 1000L, ...) {
  res <- .rotate_cf_orth(unclass(L), kappa = kappa, eps = eps,
                         normalize = normalize, random_starts = randomStarts,
                         maxit = maxit)
  list(loadings = res$loadings, Th = res$Th)
}

# Oblique GPArotation engines, keyed by rotation name. Each returns the raw
# rotated solution ($loadings, $Th, $Phi). `k` is consumed by simplimax and
# ignored by the others.
.oblq_engines <- list(
  oblimin   = function(L, k, ...) GPArotation::oblimin(L, ...),
  quartimin = function(L, k, ...) GPArotation::quartimin(L, ...),
  simplimax = function(L, k, ...) GPArotation::simplimax(L, k = k, ...),
  bentlerQ  = function(L, k, ...) GPArotation::bentlerQ(L, ...),
  geominQ   = function(L, k, ...) GPArotation::geominQ(L, ...),
  bifactorQ = function(L, k, ...) GPArotation::bifactorQ(L, ...)
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
  signs <- sign(colSums(loadings))
  signs[signs == 0] <- 1
  loadings <- loadings %*% diag(signs)

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

    AV <- .orth_engines[[rotation]](L, eps = precision,
                                    normalize = resolved$normalize,
                                    randomStarts = randomStarts, ...)
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

    AV <- .oblq_engines[[rotation]](L, k = k, eps = precision,
                                    normalize = resolved$normalize,
                                    randomStarts = randomStarts, ...)
    out <- .reflect_and_order(AV$loadings, Phi = AV$Phi, rotmat = AV$Th,
                              L_unrot = L, name_factors = FALSE)
    return(c(out, list(settings = settings)))

  }
}
