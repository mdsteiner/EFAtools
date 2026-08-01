unrot <- EFA(test_models$baseline$cormat, 3, N = 500)
obli <- .rotate_model(unrot, rotation = "oblimin", type = "EFAtools")

unrot_1 <- EFA(test_models$baseline$cormat, 1, N = 500)
obli_1 <- suppressWarnings(.rotate_model(unrot_1, rotation = "oblimin",
                                        type = "EFAtools"))

quarti <- suppressWarnings(.rotate_model(unrot, rotation = "quartimin", type = "psych"))
simpli <- suppressWarnings(.rotate_model(unrot, rotation = "simplimax", type = "SPSS",
                       maxit = 2000, randomStarts = 1))
bentQ <- suppressWarnings(.rotate_model(unrot, rotation = "bentlerQ", type = "none",
                       order_type = "eigen"))
geoQ <- suppressWarnings(.rotate_model(unrot, rotation = "geominQ", type = "EFAtools"))
bifacQ <- suppressWarnings(.rotate_model(unrot, rotation = "bifactorQ", type = "EFAtools"))

test_that("output class and dimensions are correct", {
  expect_s3_class(obli$rot_loadings, "LOADINGS")
  expect_s3_class(obli_1$rot_loadings, "LOADINGS") # The unrotated loadings here!
  expect_s3_class(quarti$rot_loadings, "LOADINGS")
  expect_s3_class(simpli$rot_loadings, "LOADINGS")
  expect_s3_class(bentQ$rot_loadings, "LOADINGS")
  expect_s3_class(geoQ$rot_loadings, "LOADINGS")
  expect_s3_class(bifacQ$rot_loadings, "LOADINGS")

  expect_output(str(obli), "List of 6")
  expect_output(str(obli_1), "List of 6")
  expect_output(str(quarti), "List of 6")
  expect_output(str(simpli), "List of 6")
  expect_output(str(bentQ), "List of 6")
  expect_output(str(geoQ), "List of 6")
  expect_output(str(bifacQ), "List of 6")

  checkmate::expect_matrix(obli$Phi)
  expect_s3_class(obli$Structure, "LOADINGS")
  checkmate::expect_matrix(obli$rotmat)
  checkmate::expect_matrix(obli$vars_accounted_rot)

  checkmate::expect_matrix(quarti$Phi)
  expect_s3_class(quarti$Structure, "LOADINGS")
  checkmate::expect_matrix(quarti$rotmat)
  checkmate::expect_matrix(quarti$vars_accounted_rot)

  checkmate::expect_matrix(simpli$Phi)
  expect_s3_class(simpli$Structure, "LOADINGS")
  checkmate::expect_matrix(simpli$rotmat)
  checkmate::expect_matrix(simpli$vars_accounted_rot)

  checkmate::expect_matrix(bentQ$Phi)
  expect_s3_class(bentQ$Structure, "LOADINGS")
  checkmate::expect_matrix(bentQ$rotmat)
  checkmate::expect_matrix(bentQ$vars_accounted_rot)

  checkmate::expect_matrix(geoQ$Phi)
  expect_s3_class(geoQ$Structure, "LOADINGS")
  checkmate::expect_matrix(geoQ$rotmat)
  checkmate::expect_matrix(geoQ$vars_accounted_rot)

  checkmate::expect_matrix(bifacQ$Phi)
  expect_s3_class(bifacQ$Structure, "LOADINGS")
  checkmate::expect_matrix(bifacQ$rotmat)
  checkmate::expect_matrix(bifacQ$vars_accounted_rot)

  expect_null(obli_1$Phi)
  expect_null(obli_1$Structure)
  expect_equal(obli_1$rotmat, NA)
  expect_null(obli_1$vars_accounted_rot)
})

test_that("every oblique rotation labels its factor columns", {
  # Downstream consumers match factors by column name (reliability coefficients and the
  # Schmid-Leiman transformation both reorder by the "F<j>" labels), so every rotation --
  # not only the varimax-based ones -- must label the rotated pattern, the factor
  # intercorrelations, and the structure matrix. The labels are the unrotated ones
  # permuted into the new factor order, so they are a permutation of F1..Fk.
  for (rot in .oblq_rotations) {
    set.seed(42)
    fit <- suppressWarnings(
      efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
              estimator = "PAF", rotation = rot)
    )
    fac_names <- colnames(fit$rot_loadings)
    expect_equal(sort(fac_names), paste0("F", 1:3),
                 label = paste0("sorted rot_loadings colnames for ", rot))
    expect_equal(dimnames(fit$Phi), list(fac_names, fac_names),
                 label = paste0("Phi dimnames for ", rot))
    expect_equal(colnames(fit$Structure), fac_names,
                 label = paste0("Structure colnames for ", rot))
  }
})

test_that("maxit bounds every optimization the multistart solver runs", {
  skip_on_cran()

  # The screen-and-triage stage has its own short iteration budget, but it must still
  # respect a user-supplied `maxit`: with maxit = 1 no start can converge, so the fit warns
  # whether or not random starts are in play. Before the cap was applied to the triage
  # stage, the random-start fit silently returned a fully converged triage solution and the
  # two calls below disagreed about whether the same maxit had converged.
  expect_warning(
    efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500, estimator = "PAF",
            rotation = "geominQ", maxit = 1),
    class = "efa_rotation_no_convergence"
  )
  expect_warning(
    efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500, estimator = "PAF",
            rotation = "geominQ", maxit = 1,
            rotate_control = rotate_control(random_starts = 0)),
    class = "efa_rotation_no_convergence"
  )
})

test_that("settings are returned correctly", {
  expect_named(obli$settings, c("normalize", "precision", "order_type", "k", "randomStarts", "rotation_diagnostics"))
  expect_named(obli_1$settings, c("normalize", "precision", "order_type", "k", "randomStarts"))
  expect_named(quarti$settings, c("normalize", "precision", "order_type", "k", "randomStarts", "rotation_diagnostics"))
  expect_named(simpli$settings, c("normalize", "precision", "order_type", "k", "randomStarts", "rotation_diagnostics"))
  expect_named(bentQ$settings, c("normalize", "precision", "order_type", "k", "randomStarts", "rotation_diagnostics"))
  expect_named(geoQ$settings, c("normalize", "precision", "order_type", "k", "randomStarts", "rotation_diagnostics"))
  expect_named(bifacQ$settings, c("normalize", "precision", "order_type", "k", "randomStarts", "rotation_diagnostics"))

  # the diagnostic carries forwarded per-start criterion values for every native oblique
  # rotation (a dropped all_values would leave criterion_best at NA)
  for (obj in list(obli, quarti, simpli, bentQ, geoQ, bifacQ)) {
    expect_true(is.finite(obj$settings$rotation_diagnostics$criterion_best))
  }

  expect_equal(obli$settings$normalize, TRUE)
  expect_equal(obli_1$settings$normalize, TRUE)
  expect_equal(quarti$settings$normalize, TRUE)
  expect_equal(simpli$settings$normalize, TRUE)
  expect_equal(bentQ$settings$normalize, TRUE)
  expect_equal(geoQ$settings$normalize, TRUE)
  expect_equal(bifacQ$settings$normalize, TRUE)

  expect_equal(obli$settings$precision, 1e-5)
  expect_equal(obli_1$settings$precision, 1e-5)
  expect_equal(quarti$settings$precision, 1e-5)
  expect_equal(simpli$settings$precision, 1e-5)
  expect_equal(bentQ$settings$precision, 1e-5)
  expect_equal(geoQ$settings$precision, 1e-5)
  expect_equal(bifacQ$settings$precision, 1e-5)

  expect_equal(obli$settings$order_type, "eigen")
  expect_equal(obli_1$settings$order_type, "eigen")
  expect_equal(quarti$settings$order_type, "eigen")
  expect_equal(simpli$settings$order_type, "ss_factors")
  expect_equal(bentQ$settings$order_type, "eigen")
  expect_equal(geoQ$settings$order_type, "eigen")
  expect_equal(bifacQ$settings$order_type, "eigen")

  expect_equal(obli$settings$k, NA)
  expect_equal(obli_1$settings$k, NA)
  expect_equal(quarti$settings$k, NA)
  expect_equal(simpli$settings$k, 18)
  expect_equal(bentQ$settings$k, NA)
  expect_equal(geoQ$settings$k, NA)
  expect_equal(bifacQ$settings$k, NA)
})

test_that("errors etc. are thrown correctly", {
  expect_error(.rotate_model(unrot, rotation = "oblimin", type = "none"), class = "efa_type_none")

  expect_warning(.rotate_model(unrot, rotation = "oblimin", type = "EFAtools",
                              normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "oblimin", type = "EFAtools",
                              order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot, rotation = "oblimin", type = "psych",
                              normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "oblimin", type = "psych",
                              order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot, rotation = "oblimin", type = "SPSS",
                              normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "oblimin", type = "SPSS",
                              order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot_1, rotation = "oblimin", type = "EFAtools"), class = "efa_single_factor")

  # A rotation name no engine table knows aborts rather than falling off the end of the
  # dispatch and returning NULL, which would surface far downstream as a missing solution.
  expect_error(.rotate_model(unrot, rotation = "obimin", type = "EFAtools"),
               class = "efa_unknown_rotation")
})

test_that("a singular oblique transformation is rejected rather than pseudo-inverted", {
  # The oblique manifold evaluates a candidate transformation T through its inverse, so a
  # rank-deficient candidate has to be rejected: substituting an approximate
  # (pseudo-inverse) solution returns something that is not a rotation at all, with a
  # rank-deficient Phi and pattern coefficients many orders of magnitude too large. The
  # oblique Procrustes entry exposes the same checked inverse directly through its T_init
  # nonsingularity check, which is the deterministic way to exercise the guard. (The
  # rejection is raised by the compiled entry, so it carries no efa_* class.)
  set.seed(1)
  A <- matrix(stats::rnorm(12), 6, 2)
  B <- matrix(stats::rnorm(12), 6, 2)

  # rank 1: column normalization maps both columns to (1, 0)
  expect_error(.oblique_procrustes(A, B, T_init = cbind(c(1, 0), c(2, 0))),
               class = "Rcpp::exception")
  expect_type(.oblique_procrustes(A, B, T_init = diag(2)), "list")
})

test_that("oblimin at a degenerate gam returns a valid rotation and reports the collapse", {
  # gam > 0 increasingly rewards correlated factors, and by gam = 1 the oblimin solution
  # on this fixture collapses toward a single factor: the transformation goes singular.
  # Rejecting the singular candidate keeps the returned solution a genuine rotation --
  # Phi positive definite with off-diagonals strictly inside (-1, 1) -- instead of the
  # pseudo-inverse artefact (Phi off-diagonals of exactly 1, pattern coefficients of order
  # 1e16). The solution is still degenerate, so it must not be returned silently.
  set.seed(42)
  res <- suppressWarnings(.rotate_model(unrot, rotation = "oblimin", type = "EFAtools",
                                        gam = 1))
  Phi <- res$Phi

  expect_true(all(is.finite(unclass(res$rot_loadings))))
  expect_lt(max(abs(unclass(res$rot_loadings))), 1e3)
  expect_lt(max(abs(Phi[upper.tri(Phi)])), 1)
  expect_gt(min(eigen(Phi, symmetric = TRUE, only.values = TRUE)$values), 0)

  set.seed(42)
  expect_warning(
    withCallingHandlers(
      .rotate_model(unrot, rotation = "oblimin", type = "EFAtools", gam = 1),
      efa_rotation_no_convergence = function(w) invokeRestart("muffleWarning")),
    class = "efa_rotation_extreme_phi"
  )
})

test_that("an oblique rotation warns when the factors are near-collinear", {
  # A factor correlation above .9 in absolute value means two rotated factors are barely
  # distinguishable, which is the usual signature of extracting more factors than the data
  # support. The fit is still returned -- an extreme correlation is a property of the
  # solution, not a failure of the rotation -- so only the interpretation is flagged.
  expect_warning(.warn_extreme_phi(matrix(c(1, -.95, -.95, 1), 2, 2)),
                 class = "efa_rotation_extreme_phi")
  expect_silent(.warn_extreme_phi(matrix(c(1, .6, .6, 1), 2, 2)))

  # the unit diagonal never triggers the check, and a single-factor solution carries no
  # factor correlations at all
  expect_silent(.warn_extreme_phi(matrix(1, 1, 1)))
  expect_silent(.warn_extreme_phi(NULL))

  # The reported correlation is truncated, not rounded. A near-collapse must not print as
  # "1.000" -- that is the perfectly collinear Phi a rejected singular transformation would
  # have produced, and reporting it would undo the distinction. A value just above the
  # cut-off must likewise not print as "0.900", which reads as sitting at the cut-off.
  # (The condition is identified by class; only the number it carries is read back.)
  reported <- function(r) {
    w <- tryCatch(.warn_extreme_phi(matrix(c(1, r, r, 1), 2, 2)),
                  efa_rotation_extreme_phi = function(w) conditionMessage(w))
    sub(".*max \\|r\\| = ([0-9.]+).*", "\\1", w)
  }
  expect_equal(reported(0.9999846784), "0.999")
  expect_equal(reported(0.905), "0.905")

  # a well-conditioned oblique fit is not flagged. Asserted on this class alone rather
  # than with expect_silent(): the multistart rotation's own convergence flag is
  # BLAS-sensitive, and that warning is not what this test is about.
  set.seed(42)
  expect_no_warning(.rotate_model(unrot, rotation = "oblimin", type = "EFAtools"),
                    class = "efa_rotation_extreme_phi")

  # the check sits on both oblique return paths: promax finalizes its own solution and is
  # flagged the same way when its two factors are all but the same factor
  L <- cbind(c(.80, .75, .70, .65, .60, .55),
             c(.78, .74, .69, .66, .58, .56))
  expect_warning(.rotate_model(list(unrot_loadings = L), rotation = "promax",
                               type = "EFAtools"),
                 class = "efa_rotation_extreme_phi")
})

test_that("order_type orders oblique factors by the requested key", {
  # "eigen" orders by the reported SS loadings (the factor-intercorrelation-weighted
  # sum of squares diag(Phi L'L)); "ss_factors" orders by the unweighted pattern sum
  # of squares colSums(L^2). For an oblique solution these keys can rank the factors
  # differently, so the two order_type values can yield different factor orders.

  # On a real oblimin solution whose two keys rank the factors the same way, both
  # order_type values produce the identical, internally consistent solution sorted
  # (descending) by the relevant key.
  set.seed(42)
  eig <- suppressWarnings(.rotate_model(unrot, rotation = "oblimin", type = "EFAtools"))
  set.seed(42)
  ss <- suppressWarnings(.rotate_model(unrot, rotation = "oblimin", type = "SPSS"))
  expect_equal(unclass(ss$rot_loadings), unclass(eig$rot_loadings), ignore_attr = TRUE)
  expect_equal(unclass(eig$Structure), unclass(eig$rot_loadings) %*% eig$Phi,
               ignore_attr = TRUE)
  expect_equal(unclass(ss$Structure), unclass(ss$rot_loadings) %*% ss$Phi,
               ignore_attr = TRUE)
  expect_false(is.unsorted(rev(diag(eig$Phi %*% crossprod(unclass(eig$rot_loadings))))))
  expect_false(is.unsorted(rev(colSums(unclass(ss$rot_loadings)^2))))

  # Feed .reflect_and_order a fixed oblique solution whose two keys rank the factors
  # differently, and check that each order_type sorts by its own key, producing
  # different orders while keeping the solution consistent (structure = pattern %*% Phi).
  L <- rbind(
    c(0.85, 0.10, 0.05),
    c(0.80, 0.05, 0.00),
    c(0.05, 0.75, 0.55),
    c(0.10, 0.55, 0.75),
    c(0.00, 0.50, 0.70),
    c(0.05, 0.45, 0.65)
  )
  Phi <- matrix(c(1.00, 0.20, 0.20,
                  0.20, 1.00, 0.65,
                  0.20, 0.65, 1.00), 3, 3)
  expect_false(identical(order(colSums(L^2), decreasing = TRUE),
                         order(diag(Phi %*% crossprod(L)), decreasing = TRUE)))

  eig2 <- .reflect_and_order(L, Phi = Phi, rotmat = diag(3), L_unrot = L,
                             order_type = "eigen")
  ss2 <- .reflect_and_order(L, Phi = Phi, rotmat = diag(3), L_unrot = L,
                            order_type = "ss_factors")
  expect_false(is.unsorted(rev(diag(eig2$Phi %*% crossprod(unclass(eig2$rot_loadings))))))
  expect_false(is.unsorted(rev(colSums(unclass(ss2$rot_loadings)^2))))
  expect_false(isTRUE(all.equal(unclass(eig2$rot_loadings), unclass(ss2$rot_loadings),
                                check.attributes = TRUE)))
  expect_equal(unclass(eig2$Structure), unclass(eig2$rot_loadings) %*% eig2$Phi,
               ignore_attr = TRUE)
  expect_equal(unclass(ss2$Structure), unclass(ss2$rot_loadings) %*% ss2$Phi,
               ignore_attr = TRUE)
})

test_that("oblique Phi, structure, and rotmat are reflected/reordered with the loadings", {
  skip_on_cran()

  # Feed .reflect_and_order a raw oblique solution whose first factor has a negative column
  # sum, and check that the documented sign reflection is applied consistently to the
  # loadings, the factor intercorrelations, and the rotation matrix (a
  # structure == pattern %*% Phi check alone misses the Phi/rotmat reflection). This
  # exercises the package's own reflection branch directly and deterministically, rather
  # than relying on a rotation engine to return a negative-column-sum factor (an orientation
  # that varies across GPArotation versions). The raw input is a real oblimin solution with
  # its first factor negated -- still a valid oblique solution (loadings = L %*% t(solve(Th)),
  # Phi = t(Th) %*% Th) -- so reflecting it back must recover the canonical solution.
  L <- unrot$unrot_loadings
  o <- suppressWarnings(.rotate_model(unrot, rotation = "oblimin", type = "EFAtools"))
  k <- ncol(o$rot_loadings)

  neg <- replace(rep(1, k), 1L, -1)
  raw_loadings <- unclass(o$rot_loadings) %*% diag(neg)
  raw_Phi      <- diag(neg) %*% o$Phi %*% diag(neg)
  raw_rotmat   <- o$rotmat %*% diag(neg)

  res <- .reflect_and_order(raw_loadings, Phi = raw_Phi, rotmat = raw_rotmat,
                            L_unrot = L, order_type = "eigen")

  # the negative factor is reflected back to a non-negative column sum
  expect_true(all(colSums(unclass(res$rot_loadings)) >= 0))
  # loadings, Phi, and the rotation matrix are reflected consistently, recovering the
  # canonical solution
  expect_equal(unclass(res$rot_loadings), unclass(o$rot_loadings), ignore_attr = TRUE)
  expect_equal(res$Phi, o$Phi, ignore_attr = TRUE)
  expect_equal(res$rotmat, o$rotmat, ignore_attr = TRUE)
  # the structure matrix stays pattern %*% Phi, and the rotation matrix reproduces the
  # rotated pattern via L_unrot %*% t(solve(rotmat))
  expect_equal(unclass(res$Structure), unclass(res$rot_loadings) %*% res$Phi,
               ignore_attr = TRUE)
  expect_equal(unclass(L) %*% t(solve(res$rotmat)), unclass(res$rot_loadings),
               ignore_attr = TRUE, tolerance = 1e-6)
})

test_that("the bentlerQ solution satisfies the oblique structure invariants", {
  skip_on_cran()

  # bentlerQ is computed by the native gradient-projection engine. These invariants hold for any
  # valid oblique solution and need no reference package (so coverage survives the criterion
  # moving off GPArotation): the structure matrix equals pattern %*% Phi, the rotation matrix
  # reproduces the rotated pattern via the documented identity L_unrot %*% t(solve(Th)), and the
  # factor correlation matrix has a unit diagonal.
  L <- unrot$unrot_loadings

  expect_equal(unclass(bentQ$Structure), unclass(bentQ$rot_loadings) %*% bentQ$Phi,
               ignore_attr = TRUE)
  expect_equal(unclass(L) %*% t(solve(bentQ$rotmat)), unclass(bentQ$rot_loadings),
               ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(diag(bentQ$Phi), rep(1, ncol(bentQ$rot_loadings)), ignore_attr = TRUE)
})

test_that("the bifactorQ solution satisfies the oblique structure invariants", {
  skip_on_cran()

  # bifactorQ is computed by the native gradient-projection engine. These invariants hold for any
  # valid oblique solution and need no reference package (so coverage survives the criterion moving
  # off GPArotation): the structure matrix equals pattern %*% Phi, the rotation matrix reproduces
  # the rotated pattern via the documented identity L_unrot %*% t(solve(Th)), and the factor
  # correlation matrix has a unit diagonal.
  L <- unrot$unrot_loadings

  expect_equal(unclass(bifacQ$Structure), unclass(bifacQ$rot_loadings) %*% bifacQ$Phi,
               ignore_attr = TRUE)
  expect_equal(unclass(L) %*% t(solve(bifacQ$rotmat)), unclass(bifacQ$rot_loadings),
               ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(diag(bifacQ$Phi), rep(1, ncol(bifacQ$rot_loadings)), ignore_attr = TRUE)
})

test_that("the simplimax solution satisfies the oblique structure invariants", {
  skip_on_cran()

  # simplimax is computed by the native gradient-projection engine. These invariants hold for any
  # valid oblique solution and need no reference package (so coverage survives the criterion moving
  # off GPArotation): the structure matrix equals pattern %*% Phi, the rotation matrix reproduces
  # the rotated pattern via the documented identity L_unrot %*% t(solve(Th)), and the factor
  # correlation matrix has a unit diagonal.
  L <- unrot$unrot_loadings

  expect_equal(unclass(simpli$Structure), unclass(simpli$rot_loadings) %*% simpli$Phi,
               ignore_attr = TRUE)
  expect_equal(unclass(L) %*% t(solve(simpli$rotmat)), unclass(simpli$rot_loadings),
               ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(diag(simpli$Phi), rep(1, ncol(simpli$rot_loadings)), ignore_attr = TRUE)
})

rm(unrot, obli, unrot_1, obli_1, quarti, simpli, bentQ, geoQ, bifacQ)
