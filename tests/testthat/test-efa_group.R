# Tests for efa_group(): per-group fitting at a common number of factors and the
# shared consensus/reference alignment. Deterministic rotations (varimax, promax)
# and identical-group invariants are used so no random-start behaviour is relied
# on. Conditions are asserted by class, not message text.

cmat <- test_models$baseline$cormat
p <- ncol(cmat)


test_that("efa_group fits every group at a common k and returns an efa_group object", {
  bands <- list(a = cmat, b = cmat, c = cmat)
  mg <- efa_group(bands, n_factors = 3, N = 500, rotation = "varimax", estimator = "PAF")

  expect_s3_class(mg, "efa_group")
  expect_length(mg$loadings, 3)
  expect_named(mg$loadings, c("a", "b", "c"))
  expect_true(all(vapply(mg$efa, inherits, logical(1), "EFA")))
  expect_true(all(vapply(mg$efa, function(f) as.numeric(f$settings$n_factors),
                         numeric(1)) == 3))
  expect_identical(dim(mg$loadings[[1]]), c(p, 3L))
  expect_identical(dim(mg$target), c(p, 3L))
  expect_identical(mg$settings$alignment, "consensus")
  expect_null(mg$Phi)
})


test_that("identical groups align to a shared consensus target", {
  bands <- list(a = cmat, b = cmat, c = cmat)
  mg <- efa_group(bands, n_factors = 3, N = 500, rotation = "varimax")

  # With identical inputs every aligned solution collapses onto the target.
  expect_equal(mg$loadings$a, mg$loadings$b, tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(mg$loadings$a, mg$loadings$c, tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(mg$loadings$a, mg$target, tolerance = 1e-6, ignore_attr = TRUE)
})


test_that("the alignment method follows the rotation and reference_group", {
  bands <- list(a = cmat, b = cmat)
  expect_identical(
    efa_group(bands, n_factors = 3, N = 500, rotation = "varimax")$settings$alignment,
    "consensus")
  expect_identical(
    efa_group(bands, n_factors = 3, N = 500, rotation = "varimax",
              reference_group = 1)$settings$alignment,
    "reference")
})


test_that("reference alignment keeps the reference group fixed", {
  bands <- list(a = cmat, b = cmat)
  mg <- efa_group(bands, n_factors = 3, N = 500, rotation = "varimax",
                  reference_group = "b")

  expect_identical(mg$settings$alignment, "reference")
  expect_identical(mg$settings$reference_group, "b")
  # The reference group's aligned loadings are its own rotated loadings, and are
  # the shared target.
  expect_equal(mg$loadings$b, mg$efa$b$rot_loadings, ignore_attr = TRUE)
  expect_equal(mg$target, mg$efa$b$rot_loadings, ignore_attr = TRUE)
  # Identical data -> the other group aligns exactly onto the reference.
  expect_equal(mg$loadings$a, mg$loadings$b, tolerance = 1e-6, ignore_attr = TRUE)
})


test_that("an oblique rotation routes to reference-Procrustes with a note", {
  bands <- list(a = cmat, b = cmat)

  # The note fires and the requested rotation is not silently changed.
  expect_message(
    efa_group(bands, n_factors = 3, N = 500, rotation = "promax"),
    class = "efa_group_oblique_reference"
  )

  mg <- suppressMessages(
    efa_group(bands, n_factors = 3, N = 500, rotation = "promax"))
  expect_identical(mg$settings$alignment, "reference")
  expect_identical(mg$settings$rotation, "promax")
  expect_false(is.null(mg$Phi))
  expect_length(mg$Phi, 2)
  expect_identical(dim(mg$Phi[[1]]), c(3L, 3L))
})


test_that("a list of correlation matrices disables the bootstrap and threads N", {
  b1 <- WJIV_ages_14_19$cormat[1:12, 1:12]
  b2 <- WJIV_ages_20_39$cormat[1:12, 1:12]
  mg <- suppressWarnings(
    efa_group(list(young = b1, adult = b2), n_factors = 2,
              N = c(WJIV_ages_14_19$N, WJIV_ages_20_39$N), rotation = "varimax"))

  expect_identical(mg$settings$input_type, "cormat")
  expect_false(mg$settings$can_bootstrap)
  expect_equal(unname(mg$settings$N),
               c(WJIV_ages_14_19$N, WJIV_ages_20_39$N))
  expect_length(mg$loadings, 2)
  expect_identical(dim(mg$loadings$young), c(12L, 2L))
})


test_that("a single factor aligns by sign only", {
  g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
  mg <- suppressMessages(efa_group(GRiPS_raw, groups = g, n_factors = 1))

  expect_identical(mg$settings$input_type, "raw")
  expect_true(mg$settings$can_bootstrap)
  expect_identical(ncol(mg$target), 1L)
  # k = 1: the aligned loadings differ from the per-group unrotated loadings only
  # by a global sign.
  expect_equal(abs(mg$loadings$g1), abs(unclass(mg$efa$g1$unrot_loadings)),
               tolerance = 1e-8, ignore_attr = TRUE)
})


test_that("unused factor levels in groups do not create phantom groups", {
  g <- factor(rep(c("a", "b"), length.out = nrow(GRiPS_raw)),
              levels = c("a", "b", "c"))
  mg <- suppressMessages(efa_group(GRiPS_raw, groups = g, n_factors = 1))

  expect_length(mg$loadings, 2)
  expect_named(mg$loadings, c("a", "b"))
  expect_named(mg$efa, c("a", "b"))
})


test_that("Tucker congruence matches psych::factor.congruence and is symmetric", {
  skip_on_cran()
  skip_if_not_installed("psych")

  b1 <- WJIV_ages_14_19$cormat[1:12, 1:12]
  b2 <- WJIV_ages_20_39$cormat[1:12, 1:12]
  mg <- efa_group(list(young = b1, adult = b2), n_factors = 2,
                  N = c(WJIV_ages_14_19$N, WJIV_ages_20_39$N),
                  rotation = "varimax")

  cong <- mg$congruence

  # Full pairwise matrices, a groups x groups x factors matched diagonal, and no
  # degeneracy on a regular fit.
  expect_named(cong, c("matrices", "matched", "degenerate"))
  expect_identical(dim(cong$matched), c(2L, 2L, 2L))
  expect_false(any(cong$degenerate))

  # The full factor-by-factor matrix equals psych's Tucker congruence to ~1e-12.
  # psych rounds its return value to `digits`, so ask for full precision.
  ref <- psych::factor.congruence(mg$loadings$young, mg$loadings$adult,
                                  digits = 15)
  expect_equal(cong$matrices$young$adult, ref, tolerance = 1e-12,
               ignore_attr = TRUE)
  # The matched diagonal is that matrix's diagonal.
  expect_equal(cong$matched["young", "adult", ], diag(ref), tolerance = 1e-12,
               ignore_attr = TRUE)

  # The reverse pair is checked against an independent psych call with the groups
  # swapped (not merely the transpose of the forward pair), so a dropped or wrong
  # transpose would be caught.
  ref_rev <- psych::factor.congruence(mg$loadings$adult, mg$loadings$young,
                                      digits = 15)
  expect_equal(cong$matrices$adult$young, ref_rev, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(cong$matched["adult", "young", ], diag(ref_rev), tolerance = 1e-12,
               ignore_attr = TRUE)
})


test_that("identical groups have matched congruence of one", {
  bands <- list(a = cmat, b = cmat, c = cmat)
  mg <- efa_group(bands, n_factors = 3, N = 500, rotation = "varimax")

  # Matched-factor congruence between identical groups is 1 across every pair.
  expect_equal(as.numeric(mg$congruence$matched), rep(1, 3 * 3 * 3),
               tolerance = 1e-6)
  # A group with itself is exactly congruent, factor by factor.
  expect_equal(diag(mg$congruence$matrices$a$a), rep(1, 3), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_false(any(mg$congruence$degenerate))
})


test_that("a degenerate factor yields NA congruence and a flag, not an error", {
  set.seed(123)
  clean <- matrix(stats::rnorm(12), nrow = 6, ncol = 2)
  degen <- clean
  degen[, 2] <- 1e-14  # an all-near-zero factor

  # The near-zero factor must not abort the run.
  cong <- .efa_group_congruence(list(g1 = clean, g2 = degen))

  # Every pair touching the degenerate group is NA and flagged in both directions.
  expect_true(all(is.na(cong$matched["g2", , ])))
  expect_true(all(is.na(cong$matched[, "g2", ])))
  expect_true(all(is.na(cong$matrices$g1$g2)))
  expect_true(all(is.na(cong$matrices$g2$g1)))
  expect_true(cong$degenerate["g1", "g2"])
  expect_true(cong$degenerate["g2", "g1"])
  expect_true(cong$degenerate["g2", "g2"])

  # The clean group's self-pair is unaffected and finite.
  expect_false(cong$degenerate["g1", "g1"])
  expect_equal(diag(cong$matrices$g1$g1), rep(1, 2), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_true(all(is.finite(cong$matched["g1", "g1", ])))
})


test_that("efa_group guards its inputs with classed conditions", {
  # Same item count, different names -> unequal items.
  m1 <- cmat[1:5, 1:5]
  m2 <- m1
  dimnames(m2) <- list(paste0("z", 1:5), paste0("z", 1:5))
  expect_error(efa_group(list(m1, m2), n_factors = 1, N = 500),
               class = "efa_group_unequal_items")

  # Fewer than two groups.
  expect_error(efa_group(list(only = cmat), n_factors = 3, N = 500),
               class = "efa_group_too_few_groups")

  # Under-identified model (p = 4, k = 2 -> df = -1).
  m4 <- cmat[1:4, 1:4]
  expect_error(efa_group(list(m4, m4), n_factors = 2, N = 500),
               class = "efa_group_under_identified")

  # A correlation matrix cannot be split into groups.
  expect_error(efa_group(cmat, groups = rep(1:2, length.out = ncol(cmat)),
                         n_factors = 2, N = 500),
               class = "efa_group_cormat_needs_list")

  # Single raw data set without a grouping vector.
  expect_error(efa_group(GRiPS_raw, n_factors = 1),
               class = "efa_group_needs_groups")

  # A grouping vector supplied together with a list.
  expect_error(efa_group(list(cmat, cmat), groups = c(1, 2), n_factors = 3, N = 500),
               class = "efa_group_groups_with_list")

  # A list mixing raw data and a correlation matrix.
  raw <- as.data.frame(GRiPS_raw[1:50, ])
  cor_raw <- stats::cor(GRiPS_raw, use = "pairwise")
  expect_error(efa_group(list(raw, cor_raw), n_factors = 1),
               class = "efa_group_mixed_input")

  # An unknown reference group.
  expect_error(efa_group(list(a = cmat, b = cmat), n_factors = 3, N = 500,
                         reference_group = "nope"),
               class = "efa_group_bad_reference")

  # Duplicated group names in a list.
  expect_error(efa_group(list(a = cmat, a = cmat), n_factors = 3, N = 500),
               class = "efa_group_duplicate_groups")

  # A confidence level outside (0, 1).
  expect_error(efa_group(list(a = cmat, b = cmat), n_factors = 3, N = 500, ci = 1.5),
               class = "efa_group_bad_ci")

  # N length neither 1 nor one-per-group (correlation-matrix input).
  expect_error(efa_group(list(a = cmat, b = cmat), n_factors = 3, N = c(100, 200, 300)),
               class = "efa_group_bad_n")
})


# Bootstrap percentile CIs for the between-group Tucker congruences. GRiPS_raw is well
# conditioned, so every replicate should fit; a genuine (not identical-group) split
# keeps the congruences comfortably below the 1.0 ceiling so the point estimate sits
# inside the percentile interval.

test_that("a bootstrap adds percentile CIs to the matched congruences", {
  g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
  mg <- suppressMessages(suppressWarnings(
    efa_group(GRiPS_raw, groups = g, n_factors = 2, rotation = "varimax",
              estimator = "PAF", b_boot = 60, seed = 2024)))

  cong <- mg$congruence
  expect_named(cong, c("matrices", "matched", "degenerate",
                       "matched_se", "matched_ci", "n_boot"))
  expect_identical(dim(cong$matched_ci$lower), c(2L, 2L, 2L))
  expect_identical(dim(cong$matched_ci$upper), c(2L, 2L, 2L))
  expect_identical(dim(cong$matched_se), c(2L, 2L, 2L))
  expect_identical(dimnames(cong$matched_ci$lower),
                   list(c("g1", "g2"), c("g1", "g2"), c("F1", "F2")))

  # Every replicate fits on well-conditioned data.
  expect_identical(cong$n_boot, 60L)
  expect_identical(mg$settings$b_boot, 60L)
  expect_identical(mg$settings$ci, 0.95)

  # A valid interval everywhere; a self-pair is exactly congruent with itself.
  expect_true(all(cong$matched_ci$lower <= cong$matched_ci$upper))
  expect_equal(diag(cong$matched_ci$lower[, , 1]), c(1, 1), ignore_attr = TRUE)
  expect_equal(diag(cong$matched_ci$upper[, , 1]), c(1, 1), ignore_attr = TRUE)
})


test_that("the bootstrap CIs bracket the point congruence", {
  skip_on_cran()
  g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
  mg <- suppressMessages(suppressWarnings(
    efa_group(GRiPS_raw, groups = g, n_factors = 2, rotation = "varimax",
              estimator = "PAF", b_boot = 300, seed = 2024)))

  pt <- mg$congruence$matched
  lo <- mg$congruence$matched_ci$lower
  up <- mg$congruence$matched_ci$upper

  # The point estimate lies inside the percentile interval (a small tolerance guards
  # the Monte Carlo edge; the congruences here are ~0.97, comfortably below 1).
  expect_true(all(lo - 0.01 <= pt & pt <= up + 0.01, na.rm = TRUE))
  # Off-diagonal intervals are non-degenerate; the self-pairs are the [1, 1] ceiling.
  expect_true(up["g1", "g2", 1] > lo["g1", "g2", 1])
  expect_equal(as.numeric(pt[cbind(1:2, 1:2, 1)]), c(1, 1), tolerance = 1e-8)
})


test_that("np-boot congruence CIs are reproducible at 1 vs 2 workers", {
  skip_on_cran()
  # The replicate fits are parallelised across workers inside each per-group EFA() with
  # future.apply (future.seed = TRUE), while efa_group's re-alignment is deterministic,
  # so a fixed `seed` must return the same intervals regardless of the number of
  # workers. The multisession workers are fresh R processes that load the installed
  # package, so run this under devtools::check() / after devtools::install() for the
  # worker fits to match the main process. (multicore is unavailable on Windows.)
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
  run <- function() suppressMessages(suppressWarnings(
    efa_group(GRiPS_raw, groups = g, n_factors = 2, rotation = "varimax",
              estimator = "PAF", b_boot = 12, seed = 2024)))

  future::plan(future::sequential)
  one <- run()

  future::plan(future::multisession, workers = 2)
  two <- run()

  expect_equal(one$congruence$matched_ci, two$congruence$matched_ci, tolerance = 1e-10)
  expect_equal(one$congruence$matched_se, two$congruence$matched_se, tolerance = 1e-10)
  expect_identical(one$congruence$n_boot, two$congruence$n_boot)
})


test_that("a supplied seed makes the CIs reproducible and restores the RNG", {
  g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
  boot <- function() suppressMessages(suppressWarnings(
    efa_group(GRiPS_raw, groups = g, n_factors = 2, rotation = "varimax",
              estimator = "PAF", b_boot = 30, seed = 99)))

  set.seed(1)
  state_before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  one <- boot()
  # The seeded bootstrap leaves the caller's stream untouched.
  expect_identical(state_before,
                   get(".Random.seed", envir = globalenv(), inherits = FALSE))

  two <- boot()
  expect_equal(one$congruence$matched_ci, two$congruence$matched_ci)
})


test_that("bootstrap replicate failures are dropped with a classed warning", {
  set.seed(1)
  p <- 6L; k <- 2L; b <- 10L
  target <- matrix(stats::rnorm(p * k), p, k)
  jitter_cube <- function() {
    a <- array(NA_real_, c(p, k, b))
    for (i in seq_len(b)) a[, , i] <- target + matrix(stats::rnorm(p * k, sd = 0.05), p, k)
    a
  }
  c1 <- jitter_cube()
  c2 <- jitter_cube()
  c2[, , c(2, 5, 8)] <- NA_real_  # three unfittable replicates in group 2
  cubes <- list(g1 = c1, g2 = c2)

  expect_warning(
    res <- .efa_group_boot_congruence(cubes, target, "orthogonal", 0.95),
    class = "efa_group_boot_failed"
  )
  # Only the replicates complete in both groups are used.
  expect_identical(res$n_boot, 7L)
  expect_identical(dim(res$matched_ci$lower), c(2L, 2L, 2L))

  # Every replicate failing is a hard error, not an empty interval.
  all_na <- list(g1 = array(NA_real_, c(p, k, b)), g2 = array(NA_real_, c(p, k, b)))
  expect_error(
    .efa_group_boot_congruence(all_na, target, "orthogonal", 0.95),
    class = "efa_group_boot_all_failed"
  )
})


test_that("a correlation-matrix input skips the bootstrap with a warning", {
  b1 <- WJIV_ages_14_19$cormat[1:12, 1:12]
  b2 <- WJIV_ages_20_39$cormat[1:12, 1:12]

  expect_warning(
    mg <- efa_group(list(young = b1, adult = b2), n_factors = 2,
                    N = c(WJIV_ages_14_19$N, WJIV_ages_20_39$N),
                    rotation = "varimax", b_boot = 50),
    class = "efa_group_boot_unavailable"
  )
  # No interval fields are added on the point-estimate-only path.
  expect_named(mg$congruence, c("matrices", "matched", "degenerate"))
  expect_identical(mg$settings$b_boot, 0L)
})


test_that("an se passed through ... is dropped with a classed warning", {
  # efa_group manages the SE method itself (bootstrap via b_boot); a stray se must not
  # trigger an unrequested, unseeded per-group bootstrap or leave payload in out$efa.
  # `seed` is passed explicitly so `se` is not partial-matched to the `seed` formal and
  # actually reaches `...`.
  g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
  expect_warning(
    mg <- suppressMessages(
      efa_group(GRiPS_raw, groups = g, n_factors = 1, se = "np-boot", seed = 1)),
    class = "efa_group_se_ignored"
  )
  expect_null(mg$efa$g1$SE)
  expect_null(mg$efa$g1$replicates)
  expect_named(mg$congruence, c("matrices", "matched", "degenerate"))
})


# Cross-group loading differences, per-item salience flags, and the approximate-invariance
# verdict. The table builders are exercised directly on hand-built aligned loadings (so an
# injected non-invariant item is fully controlled), and the classifier against known phi.

test_that(".invariance_band maps congruences to the Lorenzo-Seva & ten Berge bands", {
  expect_identical(
    .invariance_band(c(0.98, 0.95, 0.90, 0.85, 0.8499, 0.50, NA)),
    c("equal", "equal", "fair", "fair", "incongruent", "incongruent", NA)
  )
})


test_that(".efa_group_diffs flags an injected non-invariant loading and summarises the pair", {
  L1 <- matrix(c(0.70, 0.65, 0.10, 0.15,
                 0.12, 0.10, 0.72, 0.68), nrow = 4, ncol = 2,
               dimnames = list(paste0("i", 1:4), c("F1", "F2")))
  L2 <- L1
  L2["i2", "F1"] <- L1["i2", "F1"] + 0.25  # one salient non-invariant cell
  loadings <- list(g1 = L1, g2 = L2)

  dd <- .efa_group_diffs(loadings, delta = 0.10)

  # Flag table: one row per item x factor, exactly the injected cell flagged.
  expect_named(dd$flags, c("group_1", "group_2", "indicator", "factor", "diff",
                           "abs_diff", "flagged", "ci_lower", "ci_upper", "ci_excludes_0"))
  expect_identical(nrow(dd$flags), 8L)
  expect_equal(sum(dd$flags$flagged), 1)
  row <- dd$flags[dd$flags$flagged, ]
  expect_identical(row$indicator, "i2")
  expect_identical(row$factor, "F1")
  expect_equal(row$abs_diff, 0.25)
  expect_equal(row$diff, -0.25)  # g1 - g2
  # No bootstrap supplied -> the interval columns are NA.
  expect_true(all(is.na(dd$flags$ci_lower)))
  expect_true(all(is.na(dd$flags$ci_upper)))
  expect_true(all(is.na(dd$flags$ci_excludes_0)))

  # Diffs summary: one row for the single pair.
  expect_named(dd$diffs, c("group_1", "group_2", "mean_abs_diff", "median_abs_diff",
                           "min_abs_diff", "max_abs_diff", "rmse", "n_flagged"))
  expect_identical(nrow(dd$diffs), 1L)
  expect_identical(dd$diffs$group_1, "g1")
  expect_identical(dd$diffs$group_2, "g2")
  expect_equal(dd$diffs$n_flagged, 1)
  expect_equal(dd$diffs$max_abs_diff, 0.25)
  expect_equal(dd$diffs$mean_abs_diff, 0.25 / 8)

  # delta is a movable salience threshold: a larger one clears the flag, a smaller one keeps it.
  expect_equal(sum(.efa_group_diffs(loadings, delta = 0.30)$flags$flagged), 0)
  expect_equal(sum(.efa_group_diffs(loadings, delta = 0.20)$flags$flagged), 1)
})


test_that(".efa_group_diffs pairs each flag with whether the bootstrap CI excludes zero", {
  L1 <- matrix(c(0.70, 0.65, 0.10, 0.15,
                 0.12, 0.10, 0.72, 0.68), nrow = 4, ncol = 2,
               dimnames = list(paste0("i", 1:4), c("F1", "F2")))
  loadings <- list(g1 = L1, g2 = L1)

  # A synthetic difference-CI cube (groups x groups x items x factors) for the g1-g2 pair:
  # most cells bracket zero; (i2, F1) is a strictly positive interval and (i3, F2) a strictly
  # negative one, so both exclude zero.
  gn <- c("g1", "g2"); it <- paste0("i", 1:4); fn <- c("F1", "F2")
  lower <- array(NA_real_, dim = c(2, 2, 4, 2), dimnames = list(gn, gn, it, fn))
  upper <- lower
  lo_slice <- matrix(-0.05, 4, 2); up_slice <- matrix(0.05, 4, 2)
  lo_slice[2, 1] <- 0.10; up_slice[2, 1] <- 0.30    # (i2, F1): excludes 0 (positive)
  lo_slice[3, 2] <- -0.30; up_slice[3, 2] <- -0.10  # (i3, F2): excludes 0 (negative)
  lower[1, 2, , ] <- lo_slice
  upper[1, 2, , ] <- up_slice
  diff_ci <- list(lower = lower, upper = upper)

  fl <- .efa_group_diffs(loadings, delta = 0.10, diff_ci = diff_ci)$flags

  expect_equal(sum(fl$ci_excludes_0), 2)
  i2F1 <- fl$indicator == "i2" & fl$factor == "F1"
  i3F2 <- fl$indicator == "i3" & fl$factor == "F2"
  expect_true(fl$ci_excludes_0[i2F1])
  expect_true(fl$ci_excludes_0[i3F2])
  expect_equal(fl$ci_lower[i2F1], 0.10)
  expect_equal(fl$ci_upper[i2F1], 0.30)
  # A cell whose interval brackets zero is not marked.
  expect_false(fl$ci_excludes_0[fl$indicator == "i1" & fl$factor == "F1"])
})


test_that(".efa_group_invariance reads the verdict conservatively off the CI lower bound", {
  gn <- c("g1", "g2"); fn <- c("F1", "F2")
  matched <- array(NA_real_, dim = c(2, 2, 2), dimnames = list(gn, gn, fn))
  matched[1, 2, ] <- c(0.99, 0.97); matched[2, 1, ] <- c(0.99, 0.97)
  matched[1, 1, ] <- 1; matched[2, 2, ] <- 1
  lower <- array(NA_real_, dim = c(2, 2, 2), dimnames = list(gn, gn, fn))
  # F1's lower bound stays in the "equal" band; F2's drops into "fair" though its point is "equal".
  lower[1, 2, ] <- c(0.96, 0.88); lower[2, 1, ] <- c(0.96, 0.88)
  lower[1, 1, ] <- 1; lower[2, 2, ] <- 1

  inv <- .efa_group_invariance(list(matched = matched,
                                    matched_ci = list(lower = lower)))

  expect_named(inv, c("group_1", "group_2", "factor", "phi", "phi_lower", "verdict"))
  expect_identical(nrow(inv), 2L)  # one pair x two factors
  expect_identical(inv$verdict[inv$factor == "F1"], "equal")
  expect_identical(inv$verdict[inv$factor == "F2"], "fair")  # conservative: point 0.97 would be "equal"
  expect_equal(inv$phi[inv$factor == "F2"], 0.97)
  expect_equal(inv$phi_lower[inv$factor == "F2"], 0.88)

  # Without a bootstrap the verdict is point-based and phi_lower is left NA.
  inv_pt <- .efa_group_invariance(list(matched = matched))
  expect_identical(inv_pt$verdict, c("equal", "equal"))
  expect_true(all(is.na(inv_pt$phi_lower)))

  # A degenerate (NA) congruence yields an NA verdict rather than "incongruent".
  matched_deg <- matched
  matched_deg[1, 2, 2] <- NA; matched_deg[2, 1, 2] <- NA
  inv_deg <- .efa_group_invariance(list(matched = matched_deg))
  expect_true(is.na(inv_deg$verdict[inv_deg$factor == "F2"]))
})


test_that("efa_group always returns diffs/flags and gates the invariance verdict", {
  bands <- list(a = cmat, b = cmat, c = cmat)

  mg <- efa_group(bands, n_factors = 3, N = 500, rotation = "varimax", invariance = TRUE)

  expect_s3_class(mg$diffs, "data.frame")
  expect_s3_class(mg$flags, "data.frame")
  expect_s3_class(mg$invariance, "data.frame")
  # choose(3, 2) = 3 pairs; the invariance table adds a row per factor.
  expect_identical(nrow(mg$diffs), 3L)
  expect_identical(nrow(mg$invariance), 9L)
  expect_identical(nrow(mg$flags), 3L * p * 3L)

  # Identical groups: nothing salient, congruence 1, every factor "equal".
  expect_equal(sum(mg$flags$flagged), 0)
  expect_true(all(mg$diffs$n_flagged == 0))
  expect_true(all(mg$invariance$verdict == "equal"))
  expect_equal(mg$invariance$phi, rep(1, 9), tolerance = 1e-6)
  # No bootstrap -> the interval-based columns are NA.
  expect_true(all(is.na(mg$invariance$phi_lower)))
  expect_true(all(is.na(mg$flags$ci_excludes_0)))

  expect_equal(mg$settings$delta, 0.1)
  expect_true(mg$settings$invariance)

  # The verdict is gated: off by default, and diffs/flags are still present.
  mg0 <- efa_group(bands, n_factors = 3, N = 500, rotation = "varimax")
  expect_null(mg0$invariance)
  expect_false(mg0$settings$invariance)
  expect_s3_class(mg0$diffs, "data.frame")
  expect_s3_class(mg0$flags, "data.frame")
})


test_that("a bootstrap fills the flag CIs and the point difference sits inside them", {
  skip_on_cran()
  # An oblique rotation routes to the reference-Procrustes path. The displayed diff (from the
  # aligned point loadings) and the bootstrap difference CI (from the re-aligned replicate
  # cubes) must share one coordinate frame, so every point diff must lie inside its interval.
  g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
  mg <- suppressMessages(suppressWarnings(
    efa_group(GRiPS_raw, groups = g, n_factors = 2, rotation = "promax",
              estimator = "PAF", b_boot = 60, seed = 2024)))

  expect_identical(mg$settings$alignment, "reference")

  fl <- mg$flags
  # The bootstrap populates the interval columns.
  expect_false(all(is.na(fl$ci_lower)))
  expect_false(all(is.na(fl$ci_excludes_0)))
  # Every point-estimate difference lies within its own percentile interval (a small
  # tolerance guards the Monte Carlo edge; a genuine frame offset would push it far outside).
  expect_true(all(fl$ci_lower - 0.02 <= fl$diff & fl$diff <= fl$ci_upper + 0.02,
                  na.rm = TRUE))
})


# ---- print / format / plot -------------------------------------------------

# scrub_num, plus a mask for the salience-flag count `n_flagged / (p*k)`: n_flagged is a hard
# threshold count over BLAS-sensitive aligned loadings, so it can drift by one across
# platforms, and scrub_num leaves it verbatim (it carries no decimal point). Pin the layout
# and wording, not the fragile count.
scrub_group <- function(lines) sub("[0-9]+/[0-9]+", "<flagged>", scrub_num(lines))

test_that("print and format render the efa_group report", {
  local_reproducible_output()

  bands <- list(age_6_8 = WJIV_ages_6_8$cormat, age_14_19 = WJIV_ages_14_19$cormat)
  Ns <- c(WJIV_ages_6_8$N, WJIV_ages_14_19$N)

  # consensus alignment, orthogonal rotation, invariance verdicts, no bootstrap
  mg <- efa_group(bands, n_factors = 3, N = Ns, rotation = "varimax",
                  invariance = TRUE)
  expect_snapshot(print(mg), transform = scrub_group)

  # print() is cat(format()); format() returns a plain character vector
  expect_type(format(mg), "character")

  # three groups -> one row per pair; identical groups -> perfect congruence
  cmat <- test_models$baseline$cormat
  mg_ident <- efa_group(list(a = cmat, b = cmat, c = cmat), n_factors = 2, N = 500,
                        rotation = "varimax", invariance = TRUE)
  expect_snapshot(print(mg_ident), transform = scrub_group)

  # oblique rotation -> reference-alignment header, no invariance section
  mg_ob <- suppressMessages(suppressWarnings(
    efa_group(bands, n_factors = 3, N = Ns, rotation = "promax")))
  expect_snapshot(print(mg_ob), transform = scrub_group)
})


test_that("print reports the bootstrap congruence intervals and verdicts", {
  skip_on_cran()
  local_reproducible_output()

  g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
  mg <- suppressMessages(suppressWarnings(
    efa_group(GRiPS_raw, groups = g, n_factors = 1, b_boot = 50, seed = 42,
              invariance = TRUE)))
  # the congruence section notes the CIs, the differences section notes the
  # CI-excluding-0 count, and the invariance heading reads off the CI lower bound
  expect_snapshot(print(mg), transform = scrub_group)
})


test_that("plot methods return ggplot objects", {
  bands <- list(age_6_8 = WJIV_ages_6_8$cormat, age_14_19 = WJIV_ages_14_19$cormat)
  Ns <- c(WJIV_ages_6_8$N, WJIV_ages_14_19$N)
  mg <- efa_group(bands, n_factors = 3, N = Ns, rotation = "varimax")

  expect_s3_class(plot(mg, type = "congruence"), "ggplot")
  expect_s3_class(plot(mg, type = "differences"), "ggplot")
  expect_s3_class(plot(mg), "ggplot")            # default is "congruence"
  expect_error(plot(mg, type = "nope"))          # match.arg rejects an unknown type
})


test_that("the bootstrap congruence plot builds", {
  skip_on_cran()
  g <- rep(c("g1", "g2"), length.out = nrow(GRiPS_raw))
  mg <- suppressMessages(suppressWarnings(
    efa_group(GRiPS_raw, groups = g, n_factors = 2, rotation = "promax",
              estimator = "PAF", b_boot = 30, seed = 7)))
  # with a bootstrap the congruence plot carries a point-range (CI) layer; smoke only
  expect_s3_class(plot(mg, type = "congruence"), "ggplot")
  expect_s3_class(plot(mg, type = "differences"), "ggplot")
})


test_that("the deterministic plots match their vdiffr baselines", {
  skip_if_not_installed("vdiffr")
  bands <- list(age_6_8 = WJIV_ages_6_8$cormat, age_14_19 = WJIV_ages_14_19$cormat)
  Ns <- c(WJIV_ages_6_8$N, WJIV_ages_14_19$N)
  mg <- efa_group(bands, n_factors = 3, N = Ns, rotation = "varimax")

  vdiffr::expect_doppelganger("efa_group congruence plot",
                              plot(mg, type = "congruence"))
  vdiffr::expect_doppelganger("efa_group differences plot",
                              plot(mg, type = "differences"))
})
