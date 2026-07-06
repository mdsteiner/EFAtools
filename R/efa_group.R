#' Multigroup exploratory factor analysis
#'
#' Fit an exploratory factor analysis in each of several groups at a common
#' number of factors and bring the per-group solutions into one shared
#' orientation so their loadings can be compared. Each group is fitted with
#' [EFA()]; the solutions are then aligned either to a symmetric consensus target
#' or to a chosen reference group (see *Alignment*).
#'
#' @details
#' ## Input
#' Groups can be supplied in two ways: raw data together with a grouping vector
#' (`x` a data frame or matrix, `groups` one value per row), or a named list of
#' per-group data sets in `x` (with `groups` left `NULL`). The list may hold raw
#' data frames or correlation matrices (supply `N`), but not a mix of the two.
#' All groups must contain the same items in the same order; a different item set
#' or order is an error rather than being silently reordered.
#'
#' Every group is fitted at the same `n_factors` (a common-\eqn{k} multigroup
#' model). Extra arguments in `...` (for example `method`, `rotation`,
#' `cor_method`, or `type`) are forwarded unchanged to each [EFA()] call, so the
#' estimator and rotation are common to all groups.
#'
#' The \eqn{k}-factor model must be identified for the shared item set: its
#' degrees of freedom \eqn{((p - k)^2 - (p + k)) / 2} must be non-negative. Unlike
#' a single [EFA()] fit -- which only warns on an under-identified model -- a
#' multigroup fit aborts, because a shared alignment target across an
#' under-identified group is not interpretable.
#'
#' ## Alignment
#' A factor solution is identified only up to a rotation of its factors, so the
#' per-group solutions must be brought into a common orientation before their
#' loadings can be compared. Two strategies are available and are chosen
#' automatically:
#'
#' - **Consensus** (the default for orthogonal rotations and for unrotated
#'   solutions): a symmetric Generalized Procrustes Analysis target is built
#'   across all groups (Gower, 1975), and every group's loadings are rotated to
#'   it. No group is privileged.
#' - **Reference**: every group's loadings are aligned by Procrustes rotation to
#'   one reference group's loadings, which are kept fixed. This path is used when
#'   `reference_group` is given, and is used automatically for oblique rotations
#'   because the consensus iteration is not defined for oblique transforms with
#'   more than one factor. When an oblique rotation triggers the reference path
#'   without an explicit `reference_group`, the first group is used and a message
#'   reports this; the requested rotation is never silently changed.
#'
#' In both cases the returned per-group loadings share the column order and sign
#' of the returned `target`.
#'
#' ## Comparing the aligned loadings
#' Because the per-group loadings share one orientation, they can be compared cell by cell.
#' `efa_group()` reports a per-pair summary of their differences (`diffs`) and a per-item,
#' per-factor flag table (`flags`) marking cells whose absolute difference reaches `delta`.
#' `delta` is a *descriptive salience heuristic*, not a significance test; a bootstrap
#' (`b_boot > 0`) additionally reports whether each difference's confidence interval excludes
#' zero. With `invariance = TRUE`, a per-factor verdict grades the matched Tucker congruence
#' by the Lorenzo-Seva and ten Berge (2006) similarity bands -- "equal" (`>= .95`) and "fair"
#' (`[.85, .95)`) -- labelling weaker congruences "incongruent"; when a bootstrap is available
#' the verdict uses the congruence CI lower bound, so a factor is judged "equal" only if even
#' the lower bound clears the band.
#'
#' @param x A data frame or matrix of raw data (with `groups`), or a named list
#'   of per-group data sets -- either raw data frames/matrices or correlation
#'   matrices (all of one kind).
#' @param groups A vector with one value per row of `x`, giving each row's group.
#'   Only used when `x` is a single raw data set; leave `NULL` when `x` is a list.
#'   Rows with a missing group value are dropped with a warning.
#' @param n_factors numeric. The common number of factors extracted in every
#'   group.
#' @param N numeric. The number of observations per group, used only for
#'   correlation-matrix input: either a single value applied to all groups or one
#'   value per group. Ignored for raw data, where `N` is taken from each group's
#'   data. Default is `NA`.
#' @param reference_group The group to align the others to (a group name or an
#'   integer index). If `NULL` (default), orthogonal and unrotated solutions use
#'   the symmetric consensus target; oblique solutions fall back to the first
#'   group as reference. Supplying a value forces the reference alignment.
#' @param b_boot numeric. The number of non-parametric bootstrap replicates used to
#'   form percentile confidence intervals for the between-group Tucker congruences.
#'   `0` (the default) skips the bootstrap and returns the congruence point estimates
#'   only. Bootstrapping requires raw data; it is skipped with a warning for
#'   correlation-matrix input.
#' @param ci numeric. The confidence level for the bootstrap congruence intervals, a
#'   single value in `(0, 1)`. Default is `0.95`.
#' @param seed numeric or `NULL`. An optional seed making the bootstrap reproducible
#'   and independent of the number of parallel workers (the replicate fits run with
#'   [future_lapply()][future.apply::future_lapply], for which a parallel plan can be
#'   set via [future::plan()]). When supplied, the caller's random-number stream is
#'   restored afterwards, leaving no side effect. Default is `NULL`.
#' @param delta numeric. The salience threshold for the per-item loading-difference flag
#'   table: an item's loading on a factor is flagged for a group pair when the groups'
#'   aligned loadings differ by at least `delta` in absolute value. This is a descriptive
#'   salience heuristic, not a significance test; common alternatives are `0.15` and `0.20`.
#'   The threshold applies to whatever loading metric the chosen rotation produces (pattern
#'   coefficients for an oblique rotation). `0` flags every cell. Default is `0.1`.
#' @param invariance logical. Whether to add an approximate-invariance verdict per factor and
#'   group pair from the Lorenzo-Seva and ten Berge (2006) congruence bands (see *Value*).
#'   Default is `FALSE`.
#' @param ... Additional arguments passed to [EFA()] for every group (for
#'   example `method`, `rotation`, `cor_method`, `type`).
#'
#' @returns An object of class `efa_group`, a list containing:
#' \item{loadings}{A named list of the aligned per-group loading matrices. Their
#'   columns match the columns of `target` in order and sign.}
#' \item{target}{The alignment target: the symmetric consensus target, or the
#'   reference group's own loadings.}
#' \item{Phi}{A named list of the aligned per-group factor intercorrelations for
#'   an oblique rotation; `NULL` otherwise.}
#' \item{congruence}{Tucker congruence between the aligned group loadings, a list
#'   with: `matrices`, a nested list whose `[[g]][[h]]` element is the
#'   factor-by-factor congruence matrix between the aligned loadings of groups
#'   `g` and `h`; `matched`, a groups-by-groups-by-factors array of the
#'   matched-factor congruences (the diagonal of each pairwise matrix); and
#'   `degenerate`, a groups-by-groups logical matrix flagging pairs whose
#'   congruence is undefined (for example, a near-zero factor), for which the
#'   corresponding entries are `NA`. When `b_boot > 0` (raw data), three further
#'   elements are added: `matched_se`, the bootstrap standard error of each
#'   matched congruence; `matched_ci`, a list of `lower` and `upper` percentile
#'   confidence limits (each a groups-by-groups-by-factors array); and `n_boot`,
#'   the number of bootstrap replicates the intervals are based on.}
#' \item{diffs}{A data frame with one row per group pair summarising the differences
#'   between their aligned loadings: the mean, median, minimum, and maximum absolute
#'   difference, the root-mean-square difference (`rmse`), and `n_flagged`, the number of
#'   loading cells whose absolute difference reaches `delta`.}
#' \item{flags}{A data frame with one row per group pair, item, and factor giving the signed
#'   loading difference (`diff`), its absolute value (`abs_diff`), and `flagged`, whether it
#'   reaches `delta`. When a bootstrap was run (`b_boot > 0`, raw data), `ci_lower`,
#'   `ci_upper`, and `ci_excludes_0` add the percentile confidence interval for the
#'   difference and whether it excludes zero; otherwise these are `NA`.}
#' \item{invariance}{When `invariance = TRUE`, a data frame with one row per group pair and
#'   factor giving the matched Tucker congruence (`phi`), its bootstrap CI lower bound
#'   (`phi_lower`, `NA` without a bootstrap), and an approximate-invariance `verdict` based on
#'   the Lorenzo-Seva and ten Berge (2006) similarity bands: `phi >= 0.95` is "equal" and
#'   `[0.85, 0.95)` is "fair"; congruences `< 0.85`, below their bands, are labelled
#'   "incongruent". The verdict is read from `phi_lower` when a bootstrap is available
#'   (conservative) and from `phi` otherwise. `NULL` when `invariance = FALSE`.}
#' \item{efa}{The named list of per-group [EFA()] objects (each retains its own
#'   diagnostics, e.g. `heywood`).}
#' \item{alignment}{The alignment result: the consensus object (see
#'   [PROCRUSTES()]), or a list with the reference group, the target, and the
#'   per-group Procrustes results.}
#' \item{settings}{A list of the settings used, including the per-group `N`, the
#'   alignment method, the rotation, the estimator, the input type, and whether a
#'   bootstrap is available (`can_bootstrap`, `FALSE` for correlation-matrix
#'   input).}
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the Bootstrap*.
#' Chapman & Hall.
#'
#' Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika*, 40,
#' 33-51. doi: 10.1007/BF02291478
#'
#' Lorenzo-Seva, U., and ten Berge, J. M. F. (2006). Tucker's congruence coefficient as a
#' meaningful index of factor similarity. *Methodology*, 2, 57-64.
#' doi: 10.1027/1614-2241.2.2.57
#'
#' @family factor analysis
#'
#' @export
#'
#' @examples
#' # Raw data split by a grouping vector (unrotated, consensus alignment)
#' g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
#' mg <- efa_group(GRiPS_raw, groups = g, n_factors = 1)
#' mg$loadings
#'
#' # Per-pair difference summary and the per-item salience-flag table
#' mg$diffs
#' mg$flags
#'
#' \donttest{
#' # Percentile bootstrap confidence intervals for the between-group congruences, with an
#' # approximate-invariance verdict read conservatively off the congruence CI lower bound
#' mg_ci <- efa_group(GRiPS_raw, groups = g, n_factors = 1, b_boot = 100, seed = 42,
#'                    invariance = TRUE)
#' mg_ci$congruence$matched_ci
#' mg_ci$invariance
#'
#' # A named list of correlation matrices sharing the same items, common
#' # three-factor model, orthogonal rotation -> symmetric consensus target
#' bands <- list(age_6_8 = WJIV_ages_6_8$cormat, age_14_19 = WJIV_ages_14_19$cormat)
#' Ns <- c(WJIV_ages_6_8$N, WJIV_ages_14_19$N)
#' efa_group(bands, n_factors = 3, N = Ns, rotation = "varimax")
#'
#' # An oblique rotation aligns to a reference group (reported via a message)
#' efa_group(bands, n_factors = 3, N = Ns, rotation = "promax")
#' }
#'
efa_group <- function(x, groups = NULL, n_factors, N = NA,
                      reference_group = NULL, b_boot = 0L, ci = 0.95,
                      seed = NULL, delta = 0.1, invariance = FALSE, ...) {

  efa_args <- list(...)
  checkmate::assert_count(n_factors, positive = TRUE)
  checkmate::assert_count(b_boot)
  b_boot <- as.integer(b_boot)
  if (!checkmate::test_number(ci) || ci <= 0 || ci >= 1) {
    cli::cli_abort("{.arg ci} must be a single number strictly between 0 and 1.",
                   class = "efa_group_bad_ci")
  }
  if (!is.null(seed)) checkmate::assert_number(seed, finite = TRUE)
  checkmate::assert_number(delta, lower = 0, finite = TRUE)
  checkmate::assert_flag(invariance)

  # Resolve the input into a named list of per-group data sets, tagging whether
  # they are raw data or correlation matrices (which cannot be bootstrapped).
  resolved <- .efa_group_resolve_input(x, groups)
  group_data <- resolved$group_data
  group_names <- names(group_data)
  input_type <- resolved$input_type
  m <- length(group_data)

  # Bootstrap confidence intervals for the congruences need raw data to resample;
  # correlation-matrix input carries no cases, so a requested bootstrap is skipped
  # with a note and only the point estimates are returned.
  do_boot <- b_boot > 0L && input_type == "raw"
  if (b_boot > 0L && input_type != "raw") {
    cli::cli_warn(
      c("Bootstrap confidence intervals for the congruences need raw data.",
        "i" = "{.arg x} holds correlation matrices, which carry no cases to resample; the point estimates are returned without intervals."),
      class = "efa_group_boot_unavailable"
    )
  }

  if (m < 2L) {
    cli::cli_abort(
      c("{.fn efa_group} needs at least two groups.",
        "x" = "Found {m} group{?s}."),
      class = "efa_group_too_few_groups"
    )
  }

  # Enforce an identical item set and order across groups.
  .efa_group_check_items(group_data)
  p <- ncol(group_data[[1L]])

  # Validate an explicit reference group now, before the per-group fits, so a bad
  # name/index fails fast. The oblique fall-back reference is resolved later, once
  # the (common) rotation is known.
  ref_requested <- !is.null(reference_group)
  ref_idx_requested <- if (ref_requested) {
    .efa_group_resolve_reference(reference_group, group_names)
  }

  # The common-k model must be identified for the shared items: non-negative EFA
  # degrees of freedom (the Ledermann bound). This also rules out k >= p.
  df <- .efa_df(p, n_factors)
  if (df < 0) {
    cli::cli_abort(
      c("The {n_factors}-factor model is under-identified for {p} variable{?s}.",
        "x" = "It has {df} degree{?s} of freedom; a non-negative value is required.",
        "i" = "Extract fewer factors or include more variables."),
      class = "efa_group_under_identified"
    )
  }

  # N to pass to EFA: derived from the data for raw groups, supplied for
  # correlation matrices (scalar recycled or one value per group).
  Ns_in <- .efa_group_resolve_N(N, m, input_type)

  # A supplied `seed` makes the whole congruence bootstrap reproducible and
  # independent of the number of parallel workers, and leaves the caller's RNG stream
  # untouched: it is saved and restored on exit (or, if none existed, the state
  # set.seed() creates is removed again). The per-group EFA() fits are called with
  # seed = NULL so they advance this one seeded stream in sequence rather than each
  # resetting it; EFA()'s own future.seed = TRUE keeps the replicate fits
  # worker-count-independent. Mirrors the seed handling in EFA().
  if (do_boot) .set_local_seed(seed)

  # Fit every group at the common number of factors. A degenerate group (too few
  # cases, a constant item, a non-computable matrix) makes EFA abort; re-raise it
  # with the group's name so the failure is attributable. EFA warnings (e.g. a
  # Heywood case) are left to surface; the flagged variables stay in the fit's
  # `heywood` element.
  #
  # When a bootstrap is requested, each group is additionally fitted with
  # se = "np-boot" so its replicate unrotated-loading cube is produced; the cube is
  # captured before the fit is stripped back to its point-estimate form (out$efa stays
  # lean and identical in shape whether or not a bootstrap was run). The loop runs
  # strictly serially in a fixed group order: each group advances the shared (seeded)
  # RNG stream, so group g+1's resampling depends on group g -- do not parallelise it.
  #
  # efa_group controls the SE method itself, so an `se` passed through `...` is dropped
  # (with a note): otherwise it would trigger an unrequested, unseeded per-group
  # bootstrap or leave SE/replicate payload in out$efa. Congruence intervals are
  # requested with `b_boot`, not `se`.
  boot_efa_args <- efa_args
  if ("se" %in% names(boot_efa_args)) {
    cli::cli_warn(
      "{.arg se} is ignored by {.fn efa_group}; use {.arg b_boot} for congruence confidence intervals.",
      class = "efa_group_se_ignored"
    )
    boot_efa_args$se <- NULL
  }
  boot_cubes <- if (do_boot) stats::setNames(vector("list", m), group_names)

  fits <- vector("list", m)
  names(fits) <- group_names
  for (g in seq_len(m)) {
    fit_args <- c(list(x = group_data[[g]], n_factors = n_factors, N = Ns_in[[g]]),
                  boot_efa_args)
    if (do_boot) {
      fit_args <- c(fit_args,
                    list(se = "np-boot", b_boot = b_boot, ci = ci, seed = NULL))
    }
    fit_g <- tryCatch(
      do.call(EFA, fit_args),
      error = function(e) {
        cli::cli_abort(
          c("The {.fn EFA} fit failed for group {.val {group_names[[g]]}}.",
            "x" = conditionMessage(e)),
          class = "efa_group_fit_failed", parent = e
        )
      }
    )
    if (do_boot) {
      boot_cubes[[g]] <- fit_g$replicates$unrot_loadings
      fit_g <- .efa_strip_boot(fit_g)
    }
    fits[[g]] <- fit_g
  }

  # The estimator and rotation are common (one `...` applied to every group), so
  # read them from the first fit. `settings$rotation` is always populated.
  settings1 <- fits[[1L]]$settings
  rotation <- settings1$rotation
  rotation_type <- if (is.null(rotation)) "none" else .rotation_family(rotation)
  oblique <- rotation_type == "oblique" && n_factors >= 2L

  # Loadings as plain matrices (EFA returns them classed as LOADINGS).
  unrot_loadings <- lapply(fits, function(f) .change_class(f$unrot_loadings, "matrix"))
  rot_loadings <- if (rotation_type != "none") {
    lapply(fits, function(f) .change_class(f$rot_loadings, "matrix"))
  }

  # The reference path is used when a reference group is named, and always for an
  # oblique rotation (the consensus iteration is undefined for oblique k > 1).
  use_reference <- ref_requested || oblique

  if (use_reference) {
    ref_idx <- if (ref_requested) ref_idx_requested else 1L

    # Surface the automatic routing when an oblique rotation forced it without an
    # explicit choice; the rotation itself is kept, not switched to orthogonal.
    if (oblique && !ref_requested) {
      cli::cli_inform(
        c("Oblique rotations are aligned to a reference group, not a symmetric consensus target.",
          "i" = "The consensus target is undefined for oblique rotations with more than one factor.",
          "i" = "Group {.val {group_names[[ref_idx]]}} is used as the reference; set {.arg reference_group} to choose another."),
        class = "efa_group_oblique_reference"
      )
    }

    ref_loadings <- if (rotation_type == "none") {
      unrot_loadings[[ref_idx]]
    } else {
      rot_loadings[[ref_idx]]
    }
    ref_phi <- if (oblique) .change_class(fits[[ref_idx]]$Phi, "matrix")

    aligned <- .efa_group_align_reference(
      unrot_loadings = unrot_loadings, ref_loadings = ref_loadings,
      ref_idx = ref_idx, ref_phi = ref_phi, oblique = oblique,
      group_names = group_names
    )
    alignment_method <- "reference"
  } else {
    aligned <- .efa_group_align_consensus(
      unrot_loadings = unrot_loadings, rot_loadings = rot_loadings,
      rotation_type = rotation_type, group_names = group_names
    )
    alignment_method <- "consensus"
  }

  # Tucker congruence between the aligned group loadings (full pairwise matrices,
  # the matched-factor diagonals, and a flag for any degenerate pair).
  congruence <- .efa_group_congruence(aligned$loadings)

  # Percentile-bootstrap confidence intervals for the matched congruences (raw data
  # and b_boot > 0). Each group's replicate unrotated loadings are re-aligned to the
  # frozen point-estimate target with the same Procrustes mode as the point estimate,
  # and the per-factor congruence is recomputed per replicate. This runs in the main
  # process (serial) under the seed umbrella above, so it is reproducible and
  # worker-count-independent.
  if (do_boot) {
    proc_rotation <- if (oblique) "oblique" else "orthogonal"
    boot <- .efa_group_boot_congruence(boot_cubes, aligned$target, proc_rotation, ci)
    congruence$matched_se <- boot$matched_se
    congruence$matched_ci <- boot$matched_ci
    congruence$n_boot <- boot$n_boot
  }

  # Cross-group loading differences on the aligned solutions: a per-pair magnitude summary
  # and a per-item, per-factor salience-flag table, paired with the bootstrap difference
  # intervals when a bootstrap was run.
  diffs <- .efa_group_diffs(aligned$loadings, delta,
                            diff_ci = if (do_boot) boot$diff_ci else NULL)

  # Optional approximate-invariance verdict per factor and pair (Lorenzo-Seva & ten Berge,
  # 2006), read conservatively off the congruence CI lower bound when bootstrapped.
  invariance_tbl <- if (isTRUE(invariance)) .efa_group_invariance(congruence)

  # The N actually used per group (derived for raw data, supplied otherwise).
  Ns_used <- vapply(fits, function(f) {
    n <- f$settings$N
    if (is.null(n)) NA_real_ else as.numeric(n)
  }, numeric(1L))
  names(Ns_used) <- group_names

  settings <- list(
    n_factors = n_factors,
    N = Ns_used,
    reference_group = if (use_reference) group_names[[ref_idx]] else NULL,
    alignment = alignment_method,
    rotation = rotation,
    rotation_family = rotation_type,
    method = settings1$method,
    cor_method = settings1$cor_method,
    input_type = input_type,
    can_bootstrap = input_type == "raw",
    b_boot = if (do_boot) b_boot else 0L,
    ci = ci,
    delta = delta,
    invariance = invariance,
    groups = group_names,
    efa_args = efa_args
  )

  .new_efa_group(
    loadings = aligned$loadings,
    target = aligned$target,
    Phi = aligned$Phi,
    congruence = congruence,
    diffs = diffs$diffs,
    flags = diffs$flags,
    invariance = invariance_tbl,
    efa = fits,
    alignment = aligned$alignment,
    settings = settings
  )
}


# Assemble an efa_group object with its stable field set and class. The single place the
# object's shape is defined, so every field is listed here in a fixed order.
.new_efa_group <- function(loadings, target, Phi, congruence, diffs, flags,
                           invariance, efa, alignment, settings) {
  structure(
    list(
      loadings = loadings,
      target = target,
      Phi = Phi,
      congruence = congruence,
      diffs = diffs,
      flags = flags,
      invariance = invariance,
      efa = efa,
      alignment = alignment,
      settings = settings
    ),
    class = "efa_group"
  )
}


# Resolve `x`/`groups` into a named list of per-group data plus the input type
# ("raw" or "cormat"). Single raw data are split by `groups`; a list is taken as
# already-split groups (all raw or all correlation matrices, never a mix).
.efa_group_resolve_input <- function(x, groups) {

  if (is.data.frame(x) || is.matrix(x)) {
    if (.is_cormat(x)) {
      cli::cli_abort(
        c("A correlation matrix cannot be split into groups.",
          "x" = "{.arg x} looks like a correlation matrix, which carries no cases to split.",
          "i" = "Supply raw data with a {.arg groups} vector, or a named list of per-group correlation matrices."),
        class = "efa_group_cormat_needs_list"
      )
    }
    if (is.null(groups)) {
      cli::cli_abort(
        c("A single data set needs a {.arg groups} vector.",
          "i" = "Supply {.arg groups} (one value per row of {.arg x}), or pass a named list of per-group data sets."),
        class = "efa_group_needs_groups"
      )
    }
    if (length(groups) != nrow(x)) {
      cli::cli_abort(
        c("{.arg groups} must have one value per row of {.arg x}.",
          "x" = "{.arg x} has {nrow(x)} row{?s} but {.arg groups} has length {length(groups)}."),
        class = "efa_group_groups_length"
      )
    }

    g <- as.factor(groups)
    na_g <- is.na(g)
    if (any(na_g)) {
      cli::cli_warn(
        "{sum(na_g)} row{?s} with a missing group value {?was/were} dropped.",
        class = "efa_group_na_group"
      )
      x <- x[!na_g, , drop = FALSE]
      g <- g[!na_g]
    }
    # Drop empty levels -- both those left by the NA removal above and any unused
    # levels a factor `groups` already carried -- so they do not become phantom
    # zero-row groups.
    g <- droplevels(g)

    lv <- levels(g)
    group_data <- lapply(lv, function(l) x[g == l, , drop = FALSE])
    names(group_data) <- lv
    return(list(group_data = group_data, input_type = "raw"))
  }

  if (is.list(x)) {
    if (!is.null(groups)) {
      cli::cli_abort(
        c("{.arg groups} is only used with a single data set.",
          "i" = "When {.arg x} is a list of groups, drop {.arg groups}."),
        class = "efa_group_groups_with_list"
      )
    }

    for (i in seq_along(x)) {
      if (!inherits(x[[i]], c("matrix", "data.frame"))) {
        cli::cli_abort(
          c("Every element of the group list must be a matrix or data frame.",
            "x" = "Element {i} is {.obj_type_friendly {x[[i]]}}."),
          class = "efa_group_input"
        )
      }
    }

    is_cm <- vapply(x, .is_cormat, logical(1L))
    if (length(unique(is_cm)) > 1L) {
      cli::cli_abort(
        c("The group list mixes correlation matrices and raw data.",
          "i" = "Supply either all raw data sets or all correlation matrices."),
        class = "efa_group_mixed_input"
      )
    }

    nm <- names(x)
    gen <- paste0("group", seq_along(x))
    if (is.null(nm)) {
      nm <- gen
    } else {
      nm[!nzchar(nm) | is.na(nm)] <- gen[!nzchar(nm) | is.na(nm)]
    }

    # Duplicated names would leave the second group unaddressable by name and make
    # a name-based `reference_group` ambiguous, so require them to be unique.
    if (anyDuplicated(nm)) {
      dup <- unique(nm[duplicated(nm)])
      cli::cli_abort(
        c("The groups must have unique names.",
          "x" = "Duplicated group name{?s}: {.val {dup}}."),
        class = "efa_group_duplicate_groups"
      )
    }

    group_data <- x
    names(group_data) <- nm
    return(list(group_data = group_data,
                input_type = if (isTRUE(is_cm[[1L]])) "cormat" else "raw"))
  }

  cli::cli_abort(
    c("{.arg x} must be raw data (with {.arg groups}) or a named list of per-group data sets.",
      "x" = "You supplied {.obj_type_friendly {x}}."),
    class = "efa_group_input"
  )
}


# Abort unless every group has the same items in the same order. Item identity is
# compared by column name when names are present, otherwise only the count.
.efa_group_check_items <- function(group_data) {
  ref_names <- colnames(group_data[[1L]])
  ref_p <- ncol(group_data[[1L]])

  for (i in seq_along(group_data)[-1L]) {
    gi <- group_data[[i]]
    if (ncol(gi) != ref_p) {
      cli::cli_abort(
        c("All groups must have the same items.",
          "x" = "Group {i} has {ncol(gi)} item{?s} but group 1 has {ref_p}."),
        class = "efa_group_unequal_items"
      )
    }
    ni <- colnames(gi)
    if (!identical(ni, ref_names)) {
      msg <- if (!is.null(ni) && !is.null(ref_names) && setequal(ni, ref_names)) {
        c("All groups must have their items in the same order.",
          "x" = "Group {i} has the same items as group 1 but in a different order.",
          "i" = "Reorder the columns so every group matches.")
      } else {
        c("All groups must have the same items.",
          "x" = "The item names in group {i} do not match group 1.")
      }
      cli::cli_abort(msg, class = "efa_group_unequal_items")
    }
  }
  invisible(TRUE)
}


# The N to pass to each EFA call: NA for raw data (EFA derives it from the data),
# a single value recycled or one value per group for correlation matrices.
.efa_group_resolve_N <- function(N, m, input_type) {
  if (input_type == "raw") {
    return(rep(list(NA), m))
  }
  if (length(N) == 1L) {
    return(rep(list(N), m))
  }
  if (length(N) != m) {
    cli::cli_abort(
      c("{.arg N} must be a single value or one value per group.",
        "x" = "You supplied {length(N)} value{?s} for {m} group{?s}."),
      class = "efa_group_bad_n"
    )
  }
  as.list(N)
}


# Resolve `reference_group` (a group name or an integer index) to a group index.
.efa_group_resolve_reference <- function(reference_group, group_names) {
  m <- length(group_names)
  if (length(reference_group) != 1L) {
    cli::cli_abort("{.arg reference_group} must be a single group name or index.",
                   class = "efa_group_bad_reference")
  }
  if (is.character(reference_group)) {
    idx <- match(reference_group, group_names)
    if (is.na(idx)) {
      cli::cli_abort(
        c("{.arg reference_group} must name one of the groups.",
          "x" = "{.val {reference_group}} is not one of {.val {group_names}}."),
        class = "efa_group_bad_reference"
      )
    }
    return(idx)
  }
  if (is.numeric(reference_group) && is.finite(reference_group) &&
      reference_group == round(reference_group) &&
      reference_group >= 1L && reference_group <= m) {
    return(as.integer(reference_group))
  }
  cli::cli_abort(
    c("{.arg reference_group} must be a group name or an integer between 1 and {m}.",
      "x" = "You supplied {reference_group}."),
    class = "efa_group_bad_reference"
  )
}


# Consensus alignment: rotate every group's unrotated loadings to a symmetric
# GPA-consensus target. Used for orthogonal and unrotated solutions only, so the
# factor correlations are always the identity and are reported as NULL.
.efa_group_align_consensus <- function(unrot_loadings, rot_loadings,
                                       rotation_type, group_names) {
  consensus <- .gpa_consensus_target(
    unrotated_list = unrot_loadings,
    init_targets = if (rotation_type == "none") NULL else rot_loadings,
    rotation = "orthogonal"
  )

  if (!isTRUE(consensus$converged)) {
    cli::cli_warn(
      c("The consensus alignment did not meet its convergence criterion.",
        "i" = "Inspect {.code alignment$history}; consider a reference group."),
      class = "efa_group_align_failed"
    )
  }

  loadings <- consensus$aligned_loadings
  names(loadings) <- group_names

  list(loadings = loadings, target = consensus$target, Phi = NULL,
       alignment = consensus)
}


# Reference alignment: keep the reference group's loadings fixed and Procrustes-
# rotate every other group's unrotated loadings onto them (orthogonally, or
# obliquely for an oblique rotation with more than one factor).
.efa_group_align_reference <- function(unrot_loadings, ref_loadings, ref_idx,
                                       ref_phi, oblique, group_names) {
  m <- length(unrot_loadings)
  proc_rotation <- if (oblique) "oblique" else "orthogonal"

  loadings <- vector("list", m)
  names(loadings) <- group_names
  phis <- if (oblique) stats::setNames(vector("list", m), group_names)
  procrustes <- stats::setNames(vector("list", m), group_names)
  valid <- rep(TRUE, m)

  loadings[[ref_idx]] <- ref_loadings
  if (oblique) phis[[ref_idx]] <- ref_phi

  for (g in seq_len(m)[-ref_idx]) {
    pr <- PROCRUSTES(A = unrot_loadings[[g]], Target = ref_loadings,
                     rotation = proc_rotation)
    loadings[[g]] <- pr$loadings
    if (oblique) phis[[g]] <- pr$Phi
    procrustes[[g]] <- pr
    valid[g] <- isTRUE(pr$valid)
  }

  if (any(!valid)) {
    cli::cli_warn(
      c("At least one group could not be aligned to a valid rotation.",
        "i" = "Inspect {.code alignment$procrustes}; the best available alignment is used."),
      class = "efa_group_align_failed"
    )
  }

  alignment <- list(
    method = "reference",
    reference = group_names[[ref_idx]],
    reference_index = ref_idx,
    target = ref_loadings,
    procrustes = procrustes,
    valid = stats::setNames(valid, group_names)
  )

  list(loadings = loadings, target = ref_loadings,
       Phi = if (oblique) phis, alignment = alignment)
}


# Pairwise Tucker congruence between the aligned group loadings: the full
# factor-by-factor matrix for every group pair, the matched-factor congruence
# (each matrix's diagonal) as a groups-by-groups-by-factors array, and a flag for
# pairs whose congruence is undefined. A near-zero (or otherwise non-finite) factor
# makes `.tucker_congruence()` abort; that condition is caught and turned into an
# NA-flagged pair rather than failing the run, while a structural error (a dimension
# or finiteness violation of the loadings) is left to surface.
.efa_group_congruence <- function(loadings) {
  m <- length(loadings)
  group_names <- names(loadings)
  k <- ncol(loadings[[1L]])

  fac_names <- colnames(loadings[[1L]])
  if (is.null(fac_names)) fac_names <- paste0("F", seq_len(k))

  na_mat <- matrix(NA_real_, k, k, dimnames = list(fac_names, fac_names))

  matrices <- stats::setNames(vector("list", m), group_names)
  for (g in seq_len(m)) {
    matrices[[g]] <- stats::setNames(vector("list", m), group_names)
  }

  matched <- array(
    NA_real_, dim = c(m, m, k),
    dimnames = list(group_names, group_names, fac_names)
  )
  degenerate <- matrix(FALSE, m, m, dimnames = list(group_names, group_names))

  # Upper triangle and diagonal; the reverse pair is the transpose. A near-zero
  # factor makes the congruence undefined -- `.tucker_congruence()` aborts, which is
  # caught (only the undefined-congruence conditions, not a structural error) and
  # turned into an NA-filled, flagged pair rather than propagating the abort.
  for (i in seq_len(m)) {
    for (j in i:m) {
      cij <- tryCatch(
        .tucker_congruence(loadings[[i]], loadings[[j]]),
        efa_zero_column = function(e) NULL,
        efa_undefined_congruence = function(e) NULL
      )

      if (is.null(cij)) {
        matrices[[i]][[j]] <- na_mat
        degenerate[i, j] <- TRUE
        if (i != j) {
          matrices[[j]][[i]] <- na_mat
          degenerate[j, i] <- TRUE
        }
        next  # matched stays NA (the array is initialised to NA).
      }

      dimnames(cij) <- list(fac_names, fac_names)
      d <- diag(cij)
      matrices[[i]][[j]] <- cij
      matched[i, j, ] <- d
      if (i != j) {
        matrices[[j]][[i]] <- t(cij)
        matched[j, i, ] <- d
      }
    }
  }

  list(matrices = matrices, matched = matched, degenerate = degenerate)
}


# Drop the bootstrap SE/CI/replicate payload from a group fit so the stored EFA
# object reads as a point-estimate fit (out$efa is identical in shape whether or not
# a bootstrap was run; the replicate cube is captured separately before stripping).
.efa_strip_boot <- function(fit) {
  fit$SE <- NULL
  fit$CI <- NULL
  fit$replicates <- NULL
  fit$standardized_residuals <- NULL
  fit$settings$se <- "none"
  fit
}


# Percentile-bootstrap confidence intervals for the matched Tucker congruences and the
# per-item, per-factor loading differences between groups.
# `boot_cubes` holds one p x k x b cube of replicate unrotated loadings per group
# (from EFA(se = "np-boot")). Each replicate is re-aligned to the frozen point-estimate
# `target` with the same Procrustes mode as the point estimate ("orthogonal" for
# consensus / reference-orthogonal, "oblique" for reference-oblique), and the
# per-factor congruence and the pairwise loading differences are recomputed from the same
# re-aligned replicates for every group pair. A replicate contributes only if every group
# aligned in it; incomplete replicates are dropped with a single classed warning. Aggregation
# matches .boot_se_ci()/.array_se_ci(): the percentile interval at level `ci`, and the
# replicate SD as the standard error.
.efa_group_boot_congruence <- function(boot_cubes, target, proc_rotation, ci) {
  m <- length(boot_cubes)
  group_names <- names(boot_cubes)
  b <- dim(boot_cubes[[1L]])[3L]

  target <- as.matrix(target)
  p <- nrow(target)
  k <- ncol(target)
  fac_names <- colnames(target)
  if (is.null(fac_names)) fac_names <- paste0("F", seq_len(k))
  item_names <- rownames(target)
  if (is.null(item_names)) item_names <- paste0("V", seq_len(p))

  # Re-align every replicate of every group to the frozen target. The oblique solver
  # is left at its default random_starts = 0: warm-started from the closed-form
  # orthogonal Procrustes solution it is deterministic (draws no RNG), which keeps the
  # bootstrap worker-count-independent and makes the alignment invariant to the sign /
  # column permutation the replicate cubes already carry from their own alignment.
  aligned <- vector("list", m)
  valid <- matrix(FALSE, nrow = m, ncol = b)
  for (g in seq_len(m)) {
    cube <- boot_cubes[[g]]
    ag <- array(NA_real_, dim = c(p, k, b))
    for (i in seq_len(b)) {
      Li <- cube[, , i]
      # A replicate EFA() could not fit is NA-filled by EFA; skip it.
      if (!all(is.finite(Li))) next
      pr <- tryCatch(
        PROCRUSTES(A = Li, Target = target, rotation = proc_rotation),
        error = function(e) NULL
      )
      if (is.null(pr) || !isTRUE(pr$valid)) next
      ag[, , i] <- pr$loadings
      valid[g, i] <- TRUE
    }
    aligned[[g]] <- ag
  }

  # A replicate is usable only when every group aligned in it.
  complete <- apply(valid, 2L, all)
  n_valid <- sum(complete)

  if (n_valid == 0L) {
    cli::cli_abort(
      c("All {b} bootstrap replicate{?s} failed; no congruence confidence intervals could be computed.",
        "i" = "The resampled correlation matrices may be degenerate; try more observations or fewer factors."),
      class = "efa_group_boot_all_failed"
    )
  }
  if (n_valid < b) {
    n_failed <- b - n_valid
    cli::cli_warn(
      c("{n_failed} bootstrap replicate{?s} failed and {?was/were} excluded.",
        "i" = "The congruence confidence intervals are based on {n_valid} replicate{?s}."),
      class = "efa_group_boot_failed"
    )
  }

  # Matched (per-factor) congruence between every group pair, per usable replicate.
  # Both members of a pair are aligned to the same target, so factor f of one matches
  # factor f of the other and the matched congruence is the diagonal of the pairwise
  # matrix. A near-zero factor makes `.tucker_congruence()` abort; only the degeneracy
  # conditions are caught (as NA), so a structural error still surfaces.
  matched_boot <- array(
    NA_real_, dim = c(m, m, k, b),
    dimnames = list(group_names, group_names, fac_names, NULL)
  )
  for (i in which(complete)) {
    for (a in seq_len(m)) {
      for (h in a:m) {
        d <- tryCatch(
          diag(.tucker_congruence(aligned[[a]][, , i], aligned[[h]][, , i])),
          efa_zero_column = function(e) rep(NA_real_, k),
          efa_undefined_congruence = function(e) rep(NA_real_, k)
        )
        matched_boot[a, h, , i] <- d
        if (a != h) matched_boot[h, a, , i] <- d
      }
    }
  }

  # se = SD of the replicates; ci = percentile interval at level `ci` (matches
  # .boot_se_ci()). Aggregation is over the replicate margin; a pair-factor cell that
  # is all-NA (every usable replicate degenerate there) yields a silent NA.
  l_ci <- (1 - ci) / 2
  probs <- c(l_ci, ci + l_ci)
  agg <- .array_se_ci(matched_boot, probs, M = c(1L, 2L, 3L))

  # Percentile CIs for the per-item, per-factor loading differences L_g - L_h between every
  # group pair, from the same re-aligned replicates. Both members of a pair share the frozen
  # target, so element (v, f) of one lines up with element (v, f) of the other and the
  # cellwise difference is well defined. Only the unordered pairs (g < h) are used downstream,
  # so the diagonal (a zero self-difference) and the lower triangle are left NA. Every complete
  # replicate has finite aligned loadings, so no all-NA quantile arises here. Percentile limits
  # only are needed (no SE), so aggregate the replicate margin directly rather than via
  # .array_se_ci(), which would also build a discarded SD array per pair.
  comp <- which(complete)
  diff_lower <- array(
    NA_real_, dim = c(m, m, p, k),
    dimnames = list(group_names, group_names, item_names, fac_names)
  )
  diff_upper <- diff_lower
  for (a in seq_len(m - 1L)) {
    for (h in (a + 1L):m) {
      dab <- aligned[[a]][, , comp, drop = FALSE] - aligned[[h]][, , comp, drop = FALSE]
      qs <- apply(dab, c(1L, 2L), stats::quantile, probs = probs, na.rm = TRUE)
      diff_lower[a, h, , ] <- qs[1L, , ]
      diff_upper[a, h, , ] <- qs[2L, , ]
    }
  }

  list(matched_se = agg$se, matched_ci = agg$ci, n_boot = n_valid,
       diff_ci = list(lower = diff_lower, upper = diff_upper))
}


# Cross-group loading differences on the already-aligned per-group loadings: a per-pair
# magnitude summary and a per-item, per-factor flag table. `loadings` is the named list of
# aligned matrices (all p x k, sharing the target's orientation), so the pairwise comparison
# is a pure difference -- the printless `.compare_loadings()` core does no reordering. A cell
# is flagged when its absolute loading difference reaches `delta` (a descriptive salience
# heuristic, not a test). When `diff_ci` (bootstrap percentile limits from
# .efa_group_boot_congruence()) is supplied, each flag is paired with whether that cell's
# difference interval excludes zero.
.efa_group_diffs <- function(loadings, delta, diff_ci = NULL) {
  m <- length(loadings)
  group_names <- names(loadings)
  p <- nrow(loadings[[1L]])
  k <- ncol(loadings[[1L]])

  item_names <- rownames(loadings[[1L]])
  if (is.null(item_names)) item_names <- paste0("V", seq_len(p))
  fac_names <- colnames(loadings[[1L]])
  if (is.null(fac_names)) fac_names <- paste0("F", seq_len(k))

  # Column-major labels (item within factor), matching as.vector() of a p x k matrix.
  indicator <- rep(item_names, times = k)
  factor <- rep(fac_names, each = p)

  n_pairs <- m * (m - 1L) / 2L
  diffs_rows <- vector("list", n_pairs)
  flags_rows <- vector("list", n_pairs)
  r <- 0L
  for (g in seq_len(m - 1L)) {
    for (h in (g + 1L):m) {
      r <- r + 1L
      cmp <- .compare_loadings(loadings[[g]], loadings[[h]], corres = FALSE)
      d <- cmp$diff
      ad <- abs(d)
      flagged <- ad >= delta

      if (is.null(diff_ci)) {
        lo <- rep(NA_real_, p * k)
        up <- rep(NA_real_, p * k)
        exc <- rep(NA, p * k)
      } else {
        lo <- as.vector(diff_ci$lower[g, h, , ])
        up <- as.vector(diff_ci$upper[g, h, , ])
        exc <- (lo > 0) | (up < 0)
      }

      flags_rows[[r]] <- data.frame(
        group_1 = group_names[g],
        group_2 = group_names[h],
        indicator = indicator,
        factor = factor,
        diff = as.vector(d),
        abs_diff = as.vector(ad),
        flagged = as.vector(flagged),
        ci_lower = lo,
        ci_upper = up,
        ci_excludes_0 = exc,
        stringsAsFactors = FALSE
      )

      diffs_rows[[r]] <- data.frame(
        group_1 = group_names[g],
        group_2 = group_names[h],
        mean_abs_diff = cmp$mean_abs_diff,
        median_abs_diff = cmp$median_abs_diff,
        min_abs_diff = cmp$min_abs_diff,
        max_abs_diff = cmp$max_abs_diff,
        rmse = cmp$g,
        n_flagged = sum(flagged),
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    diffs = do.call(rbind, diffs_rows),
    flags = do.call(rbind, flags_rows)
  )
}


# Approximate-invariance bands for a Tucker congruence, following Lorenzo-Seva & ten Berge
# (2006): phi >= .95 the factors are "equal", [.85, .95) "fair" similarity; congruences below
# .85 (below their named bands) are labelled "incongruent" here. Vectorised and NA-safe -- an
# undefined (e.g. degenerate) congruence stays NA rather than being classified.
.invariance_band <- function(phi) {
  out <- rep(NA_character_, length(phi))
  ok <- !is.na(phi)
  out[ok & phi >= 0.95] <- "equal"
  out[ok & phi >= 0.85 & phi < 0.95] <- "fair"
  out[ok & phi < 0.85] <- "incongruent"
  out
}


# Approximate-invariance verdict per group pair and factor from the matched Tucker
# congruences. Reads the point congruences (`congruence$matched`) and, when a bootstrap was
# run, their percentile CI lower bounds (`congruence$matched_ci$lower`). The verdict is read
# conservatively off the CI lower bound when available -- a factor is judged "equal" only if
# even the lower bound clears the .95 band -- otherwise off the point estimate, with the
# `phi_lower` column left NA to signal a non-conservative, point-based reading. The basis is
# chosen once by whether a bootstrap ran, so the whole table is uniformly conservative or
# uniformly point-based. A degenerate pair (NA congruence) yields an NA verdict.
.efa_group_invariance <- function(congruence) {
  matched <- congruence$matched
  dn <- dimnames(matched)
  group_names <- dn[[1L]]
  fac_names <- dn[[3L]]
  m <- dim(matched)[1L]
  k <- dim(matched)[3L]

  has_ci <- !is.null(congruence$matched_ci)
  lower <- if (has_ci) congruence$matched_ci$lower

  rows <- vector("list", m * (m - 1L) / 2L * k)
  r <- 0L
  for (g in seq_len(m - 1L)) {
    for (h in (g + 1L):m) {
      for (f in seq_len(k)) {
        r <- r + 1L
        phi <- matched[g, h, f]
        phi_lower <- if (has_ci) lower[g, h, f] else NA_real_
        basis <- if (has_ci) phi_lower else phi
        rows[[r]] <- data.frame(
          group_1 = group_names[g],
          group_2 = group_names[h],
          factor = fac_names[f],
          phi = phi,
          phi_lower = phi_lower,
          verdict = .invariance_band(basis),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}
