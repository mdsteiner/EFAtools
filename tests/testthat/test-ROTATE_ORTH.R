unrot <- EFA(test_models$baseline$cormat, 3, N = 500)
equa <- .rotate_model(unrot, rotation = "equamax", type = "EFAtools")

unrot_1 <- EFA(test_models$baseline$cormat, 1, N = 500)
equa_1 <- suppressWarnings(.rotate_model(unrot_1, rotation = "equamax",
                                        type = "EFAtools"))

quarti <- .rotate_model(unrot, rotation = "quartimax", type = "psych")
bentT <- .rotate_model(unrot, rotation = "bentlerT", type = "none",
                      order_type = "eigen")
geoT <- .rotate_model(unrot, rotation = "geominT", type = "SPSS")
bifacT <- .rotate_model(unrot, rotation = "bifactorT", type = "EFAtools")

test_that("output class and dimensions are correct", {
  expect_s3_class(equa$rot_loadings, "LOADINGS")
  expect_s3_class(equa_1$rot_loadings, "LOADINGS") # The unrotated loadings here!
  expect_s3_class(quarti$rot_loadings, "LOADINGS")
  expect_s3_class(bentT$rot_loadings, "LOADINGS")
  expect_s3_class(geoT$rot_loadings, "LOADINGS")
  expect_s3_class(bifacT$rot_loadings, "LOADINGS")

  expect_output(str(equa), "List of 4")
  expect_output(str(equa_1), "List of 4")
  expect_output(str(quarti), "List of 4")
  expect_output(str(bentT), "List of 4")
  expect_output(str(geoT), "List of 4")
  expect_output(str(bifacT), "List of 4")

  checkmate::expect_matrix(equa$rotmat)
  checkmate::expect_matrix(equa$vars_accounted_rot)

  checkmate::expect_matrix(quarti$rotmat)
  checkmate::expect_matrix(quarti$vars_accounted_rot)

  checkmate::expect_matrix(bentT$rotmat)
  checkmate::expect_matrix(bentT$vars_accounted_rot)

  checkmate::expect_matrix(geoT$rotmat)
  checkmate::expect_matrix(geoT$vars_accounted_rot)

  checkmate::expect_matrix(bifacT$rotmat)
  checkmate::expect_matrix(bifacT$vars_accounted_rot)

  expect_equal(equa_1$rotmat, NA)
  expect_null(equa_1$vars_accounted_rot)
})

test_that("every orthogonal rotation labels its factor columns", {
  # The oblique counterpart of this test lives in test-ROTATE_OBLQ.R; together they cover
  # every rotation efa_fit() offers. Downstream consumers match factors by column name, so
  # the labels must be present for the criterion rotations too, not just for varimax.
  for (rot in .orth_rotations) {
    set.seed(42)
    fit <- suppressWarnings(
      efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
              estimator = "PAF", rotation = rot)
    )
    expect_equal(sort(colnames(fit$rot_loadings)), paste0("F", 1:3),
                 label = paste0("sorted rot_loadings colnames for ", rot))
  }
})

test_that("settings are returned correctly", {
  expect_named(equa$settings, c("normalize", "precision", "order_type", "randomStarts", "rotation_diagnostics"))
  expect_named(equa_1$settings,c("normalize", "precision", "order_type", "randomStarts"))
  expect_named(quarti$settings, c("normalize", "precision", "order_type", "randomStarts", "rotation_diagnostics"))
  expect_named(bentT$settings, c("normalize", "precision", "order_type", "randomStarts", "rotation_diagnostics"))
  expect_named(geoT$settings, c("normalize", "precision", "order_type", "randomStarts", "rotation_diagnostics"))
  expect_named(bifacT$settings, c("normalize", "precision", "order_type", "randomStarts", "rotation_diagnostics"))

  # the diagnostic carries forwarded per-start criterion values for every native
  # orthogonal rotation (a dropped all_values would leave criterion_best at NA)
  for (obj in list(equa, quarti, bentT, geoT, bifacT)) {
    expect_true(is.finite(obj$settings$rotation_diagnostics$criterion_best))
  }

  expect_equal(equa$settings$normalize, TRUE)
  expect_equal(equa_1$settings$normalize, TRUE)
  expect_equal(quarti$settings$normalize, TRUE)
  expect_equal(bentT$settings$normalize, TRUE)
  expect_equal(geoT$settings$normalize, TRUE)
  expect_equal(bifacT$settings$normalize, TRUE)

  expect_equal(equa$settings$precision, 1e-5)
  expect_equal(equa_1$settings$precision, 1e-5)
  expect_equal(quarti$settings$precision, 1e-5)
  expect_equal(bentT$settings$precision, 1e-5)
  expect_equal(geoT$settings$precision, 1e-5)
  expect_equal(bifacT$settings$precision, 1e-5)

  expect_equal(equa$settings$order_type, "eigen")
  expect_equal(equa_1$settings$order_type, "eigen")
  expect_equal(quarti$settings$order_type, "eigen")
  expect_equal(bentT$settings$order_type, "eigen")
  expect_equal(geoT$settings$order_type, "ss_factors")
  expect_equal(bifacT$settings$order_type, "eigen")
})

test_that("errors etc. are thrown correctly", {
  expect_error(.rotate_model(unrot, rotation = "equamax", type = "none"), class = "efa_type_none")

  expect_warning(.rotate_model(unrot, rotation = "equamax", type = "EFAtools",
                              normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "equamax", type = "EFAtools",
                              order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot, rotation = "equamax", type = "psych",
                              normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "equamax", type = "psych",
                              order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot, rotation = "equamax", type = "SPSS",
                              normalize = FALSE), class = "efa_type_override")
  expect_warning(.rotate_model(unrot, rotation = "equamax", type = "SPSS",
                              order_type = "ss_factors"), class = "efa_type_override")

  expect_warning(.rotate_model(unrot_1, rotation = "equamax", type = "EFAtools"), class = "efa_single_factor")
})

test_that("a two-factor bifactorT request warns that the criterion is degenerate", {
  # The Jennrich-Bentler criterion sums the products of squared loadings over pairs of
  # distinct group factors, so with a general factor and a single group factor there is no
  # such pair: the criterion is identically zero, every rotation attains it, and the
  # identity start wins the tie. The returned solution is then the unrotated one, while the
  # multistart diagnostics (one distinct optimum, zero spread) read like a perfectly
  # identified rotation -- so the degeneracy has to be signalled. The oblique counterpart of
  # this test lives in test-ROTATE_OBLQ.R.
  unrot_2 <- EFA(test_models$baseline$cormat, 2, N = 500)

  expect_warning(.rotate_model(unrot_2, rotation = "bifactorT", type = "EFAtools"),
                 class = "efa_bifactor_trivial")

  # the warning reaches the public surface too, not just the rotation core
  expect_warning(
    efa_fit(test_models$baseline$cormat, n_factors = 2, N = 500, estimator = "PAF",
            rotation = "bifactorT"),
    class = "efa_bifactor_trivial"
  )

  res <- suppressWarnings(.rotate_model(unrot_2, rotation = "bifactorT",
                                        type = "EFAtools"))
  expect_equal(res$settings$rotation_diagnostics$criterion_best, 0)
  expect_equal(res$settings$rotation_diagnostics$criterion_spread, 0)

  # no rotation was performed: the returned loadings are the unrotated ones put through the
  # shared sign/order post-processing, i.e. exactly what an identity rotation produces
  ident <- .reflect_and_order(unclass(unrot_2$unrot_loadings), rotmat = diag(2),
                              L_unrot = unrot_2$unrot_loadings, order_type = "eigen")
  expect_equal(unclass(res$rot_loadings), unclass(ident$rot_loadings), ignore_attr = TRUE)

  # with three factors there are two group factors, so the criterion is non-degenerate
  set.seed(42)
  expect_no_warning(.rotate_model(unrot, rotation = "bifactorT", type = "EFAtools"),
                    class = "efa_bifactor_trivial")
})

test_that("the orthogonal engines validate their arguments before the compiled entry", {
  # The compiled entries check the same quantities, but they raise unclassed exceptions
  # whose messages name neither the function nor the argument, and a non-numeric value never
  # reaches them at all (the numeric coercion at the C boundary fails first). Validating in
  # R gives all of them one class that user code can catch. The oblique criterion parameters
  # are covered in test-ROTATE_OBLQ.R.
  cm <- test_models$baseline$cormat
  fit <- function(...) efa_fit(cm, n_factors = 3, N = 500, estimator = "PAF", ...)

  expect_error(fit(rotation = "geominT", delta = 0), class = "efa_rotation_arg")
  expect_error(fit(rotation = "geominT", delta = -0.01), class = "efa_rotation_arg")
  expect_error(fit(rotation = "geominT", delta = "a"), class = "efa_rotation_arg")

  expect_error(fit(rotation = "quartimax", maxit = -1), class = "efa_rotation_arg")
  # a fractional maxit was truncated at the C boundary, so the engine ran a different
  # budget from the one the non-convergence warning quoted
  expect_error(fit(rotation = "bentlerT", maxit = 2.7), class = "efa_rotation_arg")
})

test_that("the native orthogonal rotation matrices are orthogonal and reproduce the loadings", {
  skip_on_cran()

  # bentlerT and bifactorT are computed by the native gradient-projection engine. These
  # invariants hold for any valid orthogonal solution and need no reference package (so
  # coverage survives a criterion moving off GPArotation): the rotation matrix is orthogonal
  # (t(Th) %*% Th == I), and it reproduces the rotated loadings via the documented identity
  # L_unrot %*% Th == rot_loadings. Both solutions are checked in one loop, labelled by
  # criterion so a failure still names the engine that produced it. The
  # sign-reflection/reordering of the rotation matrix is exercised deterministically by the
  # orthogonal reflect-and-order test below.
  L <- unrot$unrot_loadings
  k <- ncol(L)
  solutions <- list(bentlerT = bentT, bifactorT = bifacT)

  for (nm in names(solutions)) {
    res <- solutions[[nm]]
    expect_equal(crossprod(res$rotmat), diag(k), ignore_attr = TRUE, tolerance = 1e-6,
                 info = nm)
    expect_equal(unclass(L) %*% res$rotmat, unclass(res$rot_loadings),
                 ignore_attr = TRUE, tolerance = 1e-6, info = nm)
  }
})

test_that("orthogonal loadings and rotmat are reflected/reordered consistently", {
  skip_on_cran()

  # Feed .reflect_and_order a raw orthogonal solution whose first factor has a negative column
  # sum, and check that the documented sign reflection is applied consistently to the loadings
  # and the rotation matrix (the orthogonal counterpart of the oblique reflect-and-order test).
  # This exercises the package's orthogonal reflection branch directly and deterministically,
  # rather than relying on the bentlerT engine to return a negative-column-sum factor (an
  # orientation that varies with the random starts). The raw input is a real bentlerT solution
  # with its first factor negated -- still a valid orthogonal solution (loadings = L %*% Th) --
  # so reflecting it back must recover the canonical solution.
  L <- unrot$unrot_loadings
  k <- ncol(bentT$rot_loadings)

  neg <- replace(rep(1, k), 1L, -1)
  raw_loadings <- unclass(bentT$rot_loadings) %*% diag(neg)
  raw_rotmat   <- bentT$rotmat %*% diag(neg)

  res <- .reflect_and_order(raw_loadings, rotmat = raw_rotmat, L_unrot = L,
                            order_type = "eigen")

  # the negative factor is reflected back to a non-negative column sum
  expect_true(all(colSums(unclass(res$rot_loadings)) >= 0))
  # loadings and the rotation matrix are reflected consistently, recovering the canonical solution
  expect_equal(unclass(res$rot_loadings), unclass(bentT$rot_loadings), ignore_attr = TRUE)
  expect_equal(res$rotmat, bentT$rotmat, ignore_attr = TRUE)
  # the rotation matrix stays orthogonal and reproduces the rotated loadings via L_unrot %*% rotmat
  expect_equal(crossprod(res$rotmat), diag(k), ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(unclass(L) %*% res$rotmat, unclass(res$rot_loadings),
               ignore_attr = TRUE, tolerance = 1e-6)
})

rm(unrot, equa, unrot_1, equa_1, quarti, bentT, geoT, bifacT)
