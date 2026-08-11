# Regression net against the bundled SPSS FACTOR reference output (`SPSS_23`, `SPSS_27`).
# Reproducing SPSS's principal axis factoring, varimax, and promax is what the `"SPSS"`
# presets exist for, and these reference objects are the only oracles in the suite that are
# SPSS output rather than another R implementation, so the presets are otherwise pinned only
# against their own literals.
#
# Two things are checked, because they can fail independently:
#
#  (a) the ROTATIONS on their own, by rotating SPSS's own `paf_load` under `type = "SPSS"`.
#      The extraction is then out of the comparison, so any deviation belongs to the Kaiser
#      varimax, the "norm" promax target, or the ss_factors ordering.
#  (b) the FULL preset path, by extracting from the correlation matrix with
#      `estimate_control(type = "SPSS")` and rotating with `rotate_control(type = "SPSS")`,
#      which adds the PAF preset itself to the comparison.
#
# Columns and signs are matched with the package's own congruence alignment before
# comparing, and agreement is reported as the largest absolute deviation of any single
# coefficient, which is the quantity the tolerances below were set from; a relative
# difference averaged over the matrix would hide a single badly placed loading.

spss_refs <- list(SPSS_23 = SPSS_23, SPSS_27 = SPSS_27)

# The simulated cases whose PAF iteration needs more than SPSS's own default of 25
# iterations to reach the stored reference (see `?SPSS_23`).
spss_slow_cases <- c("case_1a", "case_11b")

max_abs_diff <- function(x, y) max(abs(unclass(x) - unclass(y)))

# largest absolute deviation after matching columns and signs to the reference
aligned_diff <- function(reference, loadings) {
  max_abs_diff(.align_solution(L_target = reference, L = unclass(loadings))$loadings,
               reference)
}

test_that("the SPSS rotation presets reproduce the bundled SPSS FACTOR rotations", {
  skip_on_cran()

  # Tolerance: the largest deviation observed over all 17 reference solutions is 1.2e-2, on
  # the 29-variable seven-factor WJIV_3_5 promax, and it is dominated by convergence slack
  # rather than by an algorithmic difference -- tightening `precision` to 1e-8 takes it to
  # 1.1e-3. Three of the solutions reproduce to machine precision. 2e-2 therefore sits above
  # the worst honest residual and far below any genuine regression, which moves loadings in
  # the second decimal.
  for (ver in names(spss_refs)) {
    for (nm in names(spss_refs[[ver]])) {

      ent <- spss_refs[[ver]][[nm]]
      lab <- paste(ver, nm)
      unrot <- list(unrot_loadings = ent$paf_load)

      vm <- .rotate_model(unrot, rotation = "varimax", type = "SPSS")
      expect_lt(aligned_diff(ent$var_load, vm$rot_loadings), 2e-2,
                label = paste("varimax deviation,", lab))

      pm <- .rotate_model(unrot, rotation = "promax", type = "SPSS")
      al <- .align_solution(L_target = ent$pro_load, L = unclass(pm$rot_loadings),
                            Phi = unclass(pm$Phi))
      expect_lt(max_abs_diff(al$loadings, ent$pro_load), 2e-2,
                label = paste("promax deviation,", lab))
      expect_lt(max_abs_diff(al$Phi, ent$pro_phi), 2e-2,
                label = paste("promax Phi deviation,", lab))
    }
  }
})

test_that("the SPSS presets reproduce SPSS FACTOR from the correlation matrix", {
  skip_on_cran()

  # The four simulated cases are the entries whose correlation matrices ship as well
  # (`test_models`), so they are the only ones where the extraction can be checked rather
  # than taken from the reference. The stored solutions were produced with the iteration
  # limit raised above SPSS's default of 25 (see `?SPSS_23`), so the extraction is run the
  # same way; at that budget the loadings agree to about 1e-14 and the pin is set at 1e-8,
  # six orders below the 2e-2 to 3e-2 disagreement the capped iteration produces.
  #
  # The fits depend only on the case, so each is run once and then compared against both
  # stored reference versions.
  for (cs in names(test_models)) {

    tm <- test_models[[cs]]
    fit <- function(rotation, max_iter) {
      suppressWarnings(efa_fit(
        tm$cormat, n_factors = tm$n_factors, N = tm$N, estimator = "PAF",
        rotation = rotation,
        estimate_control = estimate_control(type = "SPSS", max_iter = max_iter),
        rotate_control = rotate_control(type = "SPSS")))
    }

    unrot <- fit("none", 500)
    vm <- fit("varimax", 500)
    pm <- fit("promax", 500)
    # Only the slow cases need the capped run: where the extraction converges within
    # SPSS's own limit, `unrot$iter <= 25` already says the cap never bit and the capped
    # fit would be the identical solution.
    slow <- cs %in% spss_slow_cases
    capped <- if (slow) fit("none", 25) else NULL

    for (ver in names(spss_refs)) {

      ent <- spss_refs[[ver]][[cs]]
      lab <- paste(ver, cs)

      expect_lt(aligned_diff(ent$paf_load, unrot$unrot_loadings), 1e-8,
                label = paste("PAF deviation,", lab))

      # Rotating the preset's own extraction rather than SPSS's costs nothing beyond the
      # deviation already measured in the rotation-only block above, so the same class of
      # tolerance applies; the worst case here is 3.4e-4.
      expect_lt(aligned_diff(ent$var_load, vm$rot_loadings), 1e-3,
                label = paste("varimax deviation,", lab))

      al <- .align_solution(L_target = ent$pro_load, L = unclass(pm$rot_loadings),
                            Phi = unclass(pm$Phi))
      expect_lt(max_abs_diff(al$loadings, ent$pro_load), 1e-3,
                label = paste("promax deviation,", lab))
      expect_lt(max_abs_diff(al$Phi, ent$pro_phi), 1e-3,
                label = paste("promax Phi deviation,", lab))

      # At SPSS's own iteration limit the two slow cases miss the stored reference by two
      # orders of magnitude, which is why the documentation records that the references
      # need a raised limit.
      if (slow) {
        expect_gt(aligned_diff(ent$paf_load, capped$unrot_loadings), 1e-3,
                  label = paste("capped PAF deviation,", lab))
      }
    }

    # Convergence is a property of the fit, not of the reference version it is compared
    # against, so it is asserted once per case.
    expect_equal(unrot$convergence, 0, label = paste("PAF convergence,", cs))
    if (slow) {
      expect_equal(capped$convergence, 1, label = paste("capped convergence,", cs))
      expect_gt(unrot$iter, 25, label = paste("iterations used,", cs))
    } else {
      expect_lte(unrot$iter, 25, label = paste("iterations used,", cs))
    }
  }
})

rm(spss_refs, spss_slow_cases, max_abs_diff, aligned_diff)
