# Fixtures: oblique and orthogonal EFA fits on a correlation matrix, plus a
# raw-data single-factor fit for the scores path.
efa_ob <- suppressMessages(suppressWarnings(
  EFA(test_models$baseline$cormat, n_factors = 3, N = 500, type = "EFAtools",
      method = "PAF", rotation = "oblimin")))
efa_or <- suppressMessages(suppressWarnings(
  EFA(test_models$baseline$cormat, n_factors = 3, N = 500, type = "EFAtools",
      method = "PAF", rotation = "none")))
efa_grips <- suppressMessages(EFA(GRiPS_raw, n_factors = 1, type = "EFAtools",
                                  method = "PAF"))

fs_ob_cor  <- suppressMessages(efa_scores(test_models$baseline$cormat, f = efa_ob))
fs_grips   <- suppressMessages(efa_scores(GRiPS_raw, f = efa_grips, method = "Bartlett"))


test_that("output has the expected class, names, and shapes", {
  expect_s3_class(fs_ob_cor, "efa_scores")
  expect_s3_class(fs_grips, "efa_scores")

  expect_named(fs_ob_cor, c("weights", "scores", "r.scores", "score_cor",
                            "determinacy", "settings"))
  expect_named(fs_grips, c("weights", "scores", "r.scores", "score_cor",
                           "determinacy", "settings"))

  # Weights keep the variable rownames and carry factor labels; correlation input
  # returns no scores.
  L <- unclass(efa_ob$rot_loadings)
  expect_equal(dim(fs_ob_cor$weights), dim(L))
  expect_identical(rownames(fs_ob_cor$weights), rownames(L))
  expect_identical(colnames(fs_ob_cor$weights), c("F1", "F2", "F3"))
  expect_null(fs_ob_cor$scores)

  # Raw input returns one score row per observation, one column per factor.
  expect_equal(dim(fs_grips$scores), c(nrow(GRiPS_raw), 1L))
  expect_false(anyNA(fs_grips$scores))

  # Diagnostics matrices are m x m and share the factor names.
  m <- ncol(L)
  expect_equal(dim(fs_ob_cor$r.scores), c(m, m))
  expect_equal(dim(fs_ob_cor$score_cor), c(m, m))
  expect_identical(dimnames(fs_ob_cor$r.scores), dimnames(fs_ob_cor$score_cor))
})

test_that("the determinacy table is consistent with score_cor", {
  expect_s3_class(fs_ob_cor$determinacy, "data.frame")
  expect_named(fs_ob_cor$determinacy, c("rho", "rho2", "guttman"))
  expect_identical(rownames(fs_ob_cor$determinacy), colnames(fs_ob_cor$weights))

  # The determinacy is the score_cor diagonal; rho2 and guttman are its transforms.
  expect_equal(fs_ob_cor$determinacy$rho, diag(fs_ob_cor$score_cor),
               ignore_attr = TRUE)
  expect_equal(fs_ob_cor$determinacy$rho2, fs_ob_cor$determinacy$rho^2)
  expect_equal(fs_ob_cor$determinacy$guttman, 2 * fs_ob_cor$determinacy$rho^2 - 1)
})

test_that("settings record the resolved method and dimensions", {
  expect_equal(fs_grips$settings$method, "Bartlett")
  expect_equal(fs_ob_cor$settings$method, "regression")   # the default
  expect_equal(fs_ob_cor$settings$n_factors, 3L)
  expect_true(is.na(fs_ob_cor$settings$n_obs))            # correlation input
  expect_true(is.na(fs_ob_cor$settings$n_scored))
  # n_obs counts the supplied rows, n_scored the rows that could be scored; with
  # complete data the two agree.
  expect_equal(fs_grips$settings$n_obs, nrow(GRiPS_raw))
  expect_equal(fs_grips$settings$n_scored, nrow(GRiPS_raw))
})

test_that("a constant model variable in the scoring data is an error", {
  x_const <- as.matrix(GRiPS_raw)
  x_const[, 2] <- 3
  expect_error(
    suppressMessages(efa_scores(x_const, f = efa_grips, method = "Bartlett")),
    class = "efa_scores_constant_var")

  # A data.frame takes the same route, and two constant variables render the
  # plural branch of the message.
  df_const <- as.data.frame(GRiPS_raw)
  df_const[[1]] <- 2
  df_const[[3]] <- 5
  expect_error(
    suppressMessages(efa_scores(df_const, f = efa_grips, method = "Bartlett")),
    class = "efa_scores_constant_var")

  # A model variable that is entirely missing, or that carries an infinite value, has
  # no usable spread either; each renders its own cause in the same condition.
  x_all_na <- as.matrix(GRiPS_raw)
  x_all_na[, 3] <- NA_real_
  expect_error(
    suppressMessages(efa_scores(x_all_na, f = efa_grips, method = "Bartlett")),
    class = "efa_scores_constant_var")

  x_inf <- as.matrix(GRiPS_raw)
  x_inf[1, 4] <- Inf
  expect_error(
    suppressMessages(efa_scores(x_inf, f = efa_grips, method = "Bartlett")),
    class = "efa_scores_constant_var")

  # Unnamed input reports the column position instead of a name.
  x_unnamed <- unname(as.matrix(GRiPS_raw))
  x_unnamed[, 2] <- 3
  expect_error(
    suppressMessages(efa_scores(x_unnamed, f = efa_grips, method = "Bartlett")),
    class = "efa_scores_constant_var")

  # A constant column that is not a model variable is dropped by name alignment
  # and must not trip the guard.
  x_extra <- cbind(as.matrix(GRiPS_raw), constant_extra = 1)
  expect_no_error(
    suppressMessages(efa_scores(x_extra, f = efa_grips, method = "Bartlett")))
})

test_that("cases with a missing model variable are counted and reported", {
  x_na <- as.matrix(GRiPS_raw)
  x_na[1, 1] <- NA

  expect_warning(
    fs_one <- suppressMessages(efa_scores(x_na, f = efa_grips, method = "Bartlett")),
    class = "efa_scores_incomplete_cases")

  # n_obs stays the number of supplied rows; only n_scored drops.
  expect_equal(fs_one$settings$n_obs, nrow(GRiPS_raw))
  expect_equal(fs_one$settings$n_scored, nrow(GRiPS_raw) - 1L)
  expect_true(is.na(fs_one$scores[1, 1]))
  expect_false(anyNA(fs_one$scores[-1, , drop = FALSE]))

  # The plural branch of the message is a separate rendering.
  x_na[5, 2] <- NA
  x_na[9, 3] <- NA
  expect_warning(
    fs_many <- suppressMessages(efa_scores(x_na, f = efa_grips, method = "Bartlett")),
    class = "efa_scores_incomplete_cases")
  expect_equal(fs_many$settings$n_scored, nrow(GRiPS_raw) - 3L)
})

test_that("Phi supplied alongside an efa object is ignored with a warning", {
  expect_warning(
    fs <- suppressMessages(efa_scores(test_models$baseline$cormat, f = efa_ob,
                                      Phi = diag(3))),
    class = "efa_scores_phi_ignored")
  expect_equal(fs$weights, fs_ob_cor$weights)
})

test_that("classed conditions fire as expected", {
  expect_error(efa_scores(1:5), class = "efa_input_not_matrix")
  expect_error(efa_scores(GRiPS_raw, f = 1:5), class = "efa_scores_bad_f")
  expect_message(efa_scores(test_models$baseline$cormat, f = efa_ob),
                 class = "efa_scores_needs_raw")
  # Raw input isolates the phi_null message (a correlation matrix would also emit
  # needs_raw).
  expect_message(efa_scores(GRiPS_raw, f = unclass(efa_grips$unrot_loadings)),
                 class = "efa_scores_phi_null")

  # Anderson-Rubin is undefined for a single factor.
  expect_error(efa_scores(GRiPS_raw, f = efa_grips, method = "Anderson"),
               class = "efa_scores_anderson_single")

  # Anderson-Rubin forces uncorrelated scores, so it is inappropriate for an oblique
  # solution: warn there, but not for an orthogonal one (whose Phi is the identity).
  expect_warning(
    suppressMessages(efa_scores(test_models$baseline$cormat, f = efa_ob,
                                method = "Anderson")),
    class = "efa_scores_anderson_oblique")
  expect_no_warning(
    suppressMessages(efa_scores(test_models$baseline$cormat, f = efa_or,
                                method = "Anderson")))
})

test_that("a correlation matrix in `x` is checked against the model it is scored with", {
  cm <- test_models$baseline$cormat

  # A fitted solution carries the correlations it was estimated from, so `x` is not read
  # for the weights. Unchecked, a matrix of any size and of entirely different variables
  # returned the fitted solution's weights as though it had produced them.
  five <- matrix(.2, 5, 5)
  diag(five) <- 1
  dimnames(five) <- list(paste0("X", 1:5), paste0("X", 1:5))
  expect_error(suppressMessages(efa_scores(five, f = efa_ob)),
               class = "efa_scores_matrix_dim")

  # The right size, but the variables of some other model.
  foreign <- cm
  dimnames(foreign) <- list(paste0("Q", 1:18), paste0("Q", 1:18))
  expect_error(suppressMessages(efa_scores(foreign, f = efa_ob)),
               class = "efa_scores_matrix_names")

  # Named on one axis only, which identifies neither by name nor by position: refused
  # here as it is for `rho`, rather than passed over as it was when `x` went unchecked.
  half <- cm
  colnames(half) <- NULL
  expect_error(suppressMessages(efa_scores(half, f = efa_ob)),
               class = "efa_scores_matrix_names")

  # Where `f` is a loading matrix the solution carries no correlations, so `x` is the
  # scoring matrix. Named for the model and aligned to it there, so a permuted matrix is
  # reordered rather than scored against the wrong variables. (With a fitted `f` the same
  # pair agrees whatever the order, `x` not being read, so the model axis is what this
  # can be asserted on.)
  L <- unclass(efa_ob$rot_loadings)
  p <- c(4:18, 1:3)
  expect_equal(suppressMessages(efa_scores(cm[p, p], f = L, Phi = efa_ob$Phi)),
               suppressMessages(efa_scores(cm, f = L, Phi = efa_ob$Phi)))

  # The inputs the check must not turn away: a data frame, which `x` is documented to
  # take, and a wholly unnamed matrix, which is matched by position.
  expect_equal(suppressMessages(efa_scores(as.data.frame(cm), f = efa_ob)),
               suppressMessages(efa_scores(cm, f = efa_ob)))
  unnamed <- cm
  dimnames(unnamed) <- NULL
  expect_s3_class(suppressMessages(efa_scores(unnamed, f = efa_ob)), "efa_scores")

  # The wrapper forwards to the same check, so it is refused there too.
  expect_error(suppressMessages(FACTOR_SCORES(five, f = efa_ob)),
               class = "efa_scores_matrix_dim")
})

test_that("loadings with missing or duplicate factor names are relabelled F1..Fm", {
  L <- unclass(efa_or$unrot_loadings)          # NULL colnames
  colnames(L) <- c("F1", "F1", "F2")           # duplicate factor labels

  # Must not abort with the opaque base-R "duplicate row.names" error.
  fs <- suppressMessages(efa_scores(test_models$baseline$cormat, f = L, Phi = diag(3)))
  expect_identical(colnames(fs$weights), c("F1", "F2", "F3"))
  expect_identical(rownames(fs$determinacy), c("F1", "F2", "F3"))
})

test_that("a non-Pearson EFA with rho = NULL on raw data warns, rho silences it", {
  skip_on_cran()

  efa_sp <- suppressMessages(EFA(GRiPS_raw, n_factors = 1, type = "EFAtools",
                                 method = "PAF", cor_method = "spearman"))
  expect_warning(efa_scores(GRiPS_raw, f = efa_sp, method = "regression"),
                 class = "efa_scores_cor_method")

  rho_sp <- stats::cor(GRiPS_raw, method = "spearman")
  expect_no_warning(efa_scores(GRiPS_raw, f = efa_sp, rho = rho_sp,
                               method = "regression"))
})

test_that("raw-data columns are aligned to the model by name, not position", {
  # Reversing the column order must not change the scores or weights: efa_scores
  # aligns x to the variable names carried by the model, not by column position.
  rev_idx <- rev(seq_len(ncol(GRiPS_raw)))

  fs_ordered  <- suppressMessages(
    efa_scores(GRiPS_raw, f = efa_grips, method = "Bartlett"))
  fs_shuffled <- suppressMessages(
    efa_scores(GRiPS_raw[, rev_idx], f = efa_grips, method = "Bartlett"))

  expect_equal(fs_shuffled$scores, fs_ordered$scores)
  expect_equal(fs_shuffled$weights, fs_ordered$weights)

  # Guard against a positional regression: scoring the shuffled columns by position
  # (the pre-alignment behaviour) gives materially different scores.
  positional <- scale(as.matrix(GRiPS_raw[, rev_idx])) %*% fs_ordered$weights
  expect_gt(max(abs(positional - fs_shuffled$scores)), 0.1)
})

test_that("column name alignment handles missing, extra, and unnamed columns", {
  # A model variable missing from x cannot be scored.
  expect_error(
    suppressMessages(efa_scores(GRiPS_raw[, -1], f = efa_grips, method = "Bartlett")),
    class = "efa_scores_missing_vars")

  fs_plain <- suppressMessages(
    efa_scores(GRiPS_raw, f = efa_grips, method = "Bartlett"))

  # Extra, non-model columns are ignored (dropped by name alignment).
  x_extra <- cbind(as.matrix(GRiPS_raw), extra_col = 0)
  fs_extra <- suppressMessages(
    efa_scores(x_extra, f = efa_grips, method = "Bartlett"))
  expect_equal(fs_extra$scores, fs_plain$scores)
  expect_identical(rownames(fs_extra$weights), rownames(fs_plain$weights))

  # Unnamed data fall back to positional matching (order already matches the model).
  x_unnamed <- unname(as.matrix(GRiPS_raw))
  fs_unnamed <- suppressMessages(
    efa_scores(x_unnamed, f = efa_grips, method = "Bartlett"))
  expect_equal(fs_unnamed$scores, fs_plain$scores, ignore_attr = TRUE)
})

test_that("alignment holds for the loadings-matrix path and multiple factors", {
  # A bare loading matrix with more than one factor exercises the path where the
  # scoring correlation is cor(x), recomputed from x: the reorder must precede it
  # so the weights stay consistent. Two factors also check that the score columns
  # stay matched under a column reorder. Both the standardized (Bartlett) and the
  # uncentered (components) score multiplications are covered.
  efa2 <- suppressMessages(suppressWarnings(
    EFA(GRiPS_raw, n_factors = 2, type = "EFAtools", method = "PAF",
        rotation = "oblimin")))
  L   <- unclass(efa2$rot_loadings)
  Phi <- efa2$Phi
  rev_idx <- rev(seq_len(ncol(GRiPS_raw)))

  for (m in c("Bartlett", "components")) {
    fs_ordered <- suppressMessages(suppressWarnings(
      efa_scores(GRiPS_raw, f = L, Phi = Phi, method = m)))
    fs_shuffled <- suppressMessages(suppressWarnings(
      efa_scores(GRiPS_raw[, rev_idx], f = L, Phi = Phi, method = m)))

    expect_equal(fs_shuffled$scores, fs_ordered$scores, info = m)
    expect_equal(fs_shuffled$weights, fs_ordered$weights, info = m)
    expect_equal(dim(fs_ordered$scores), c(nrow(GRiPS_raw), 2L), info = m)
  }
})

test_that("named scoring and factor correlations are aligned to the loadings", {
  L <- unclass(efa_ob$rot_loadings)
  R <- test_models$baseline$cormat
  Phi <- efa_ob$Phi

  base <- suppressMessages(efa_scores(R, f = L, Phi = Phi))

  var_perm <- rev(seq_len(nrow(R)))
  R_perm <- R[var_perm, var_perm]
  from_permuted_R <- suppressMessages(efa_scores(R_perm, f = L, Phi = Phi))
  expect_equal(from_permuted_R$weights, base$weights)
  expect_equal(from_permuted_R$score_cor, base$score_cor)

  # A correlation matrix must use the same variable order on its two labelled axes.
  col_perm <- c(seq.int(2L, ncol(R)), 1L)
  R_axis_perm <- R[var_perm, col_perm]
  expect_error(
    suppressMessages(efa_scores(R, f = L, Phi = Phi, rho = R_axis_perm)),
    class = "efa_scores_matrix_names"
  )

  fac_perm <- c(3L, 1L, 2L)
  Phi_perm <- Phi[fac_perm, fac_perm]
  from_permuted_Phi <- suppressMessages(efa_scores(R, f = L, Phi = Phi_perm))
  expect_equal(from_permuted_Phi$weights, base$weights)
  expect_equal(from_permuted_Phi$score_cor, base$score_cor)

  # rho follows the same name alignment when raw observations are scored.
  R_grips <- stats::cor(GRiPS_raw)
  rho_base <- suppressMessages(
    efa_scores(GRiPS_raw, f = efa_grips, rho = R_grips, method = "regression"))
  grips_perm <- rev(seq_len(nrow(R_grips)))
  rho_perm <- suppressMessages(
    efa_scores(GRiPS_raw, f = efa_grips,
               rho = R_grips[grips_perm, grips_perm],
               method = "regression"))
  expect_equal(rho_perm$weights, rho_base$weights)
  expect_equal(rho_perm$scores, rho_base$scores)
})

test_that("invalid scoring-correlation matrices fail before matrix algebra", {
  L <- unclass(efa_ob$rot_loadings)
  R <- test_models$baseline$cormat

  expect_error(
    suppressMessages(efa_scores(R, f = L, Phi = efa_ob$Phi,
                                rho = diag(nrow(R) - 1L))),
    class = "efa_scores_matrix_dim")

  R_asym <- R
  R_asym[1, 2] <- R_asym[1, 2] + 0.1
  expect_error(
    suppressMessages(efa_scores(R, f = L, Phi = efa_ob$Phi, rho = R_asym)),
    class = "efa_scores_matrix_symmetric")

  R_bad_names <- R
  dimnames(R_bad_names) <- list(paste0("bad", seq_len(nrow(R))),
                                paste0("bad", seq_len(ncol(R))))
  expect_error(
    suppressMessages(efa_scores(R_bad_names, f = L, Phi = efa_ob$Phi)),
    class = "efa_scores_matrix_names")

  # Matching name sets do not excuse inconsistent labels on the two axes.
  R_bad_axes <- R
  colnames(R_bad_axes) <- rev(colnames(R_bad_axes))
  expect_error(
    suppressMessages(efa_scores(R_bad_axes, f = L, Phi = efa_ob$Phi)),
    class = "efa_scores_matrix_names")
})

test_that("score labels and uniqueness inputs fail or fall back cleanly", {
  R2 <- diag(2)
  dimnames(R2) <- list(c("v1", "v2"), c("v1", "v2"))

  # Missing factor labels are not valid data-frame row names; replace the whole
  # factor axis with deterministic labels before diagnostics are constructed.
  L_bad_factor_names <- matrix(c(.6, .2, .1, .5), 2, 2,
                               dimnames = list(c("v1", "v2"), c("F1", NA)))
  relabelled <- suppressMessages(efa_scores(R2, f = L_bad_factor_names))
  expect_identical(colnames(relabelled$weights), c("F1", "F2"))

  # Named raw data can only be aligned when its labels identify columns uniquely.
  x_bad_names <- GRiPS_raw
  colnames(x_bad_names)[2] <- colnames(x_bad_names)[1]
  expect_error(
    suppressMessages(efa_scores(x_bad_names, f = efa_grips)),
    class = "efa_scores_x_names"
  )

  # Bartlett/Anderson uniqueness weighting is undefined at exactly zero; do not
  # continue into 1 / 0 and return non-finite weights or an eigensolver error.
  L_zero_u2 <- matrix(c(1, .5), 2, 1,
                      dimnames = list(c("v1", "v2"), "F1"))
  expect_error(
    suppressMessages(efa_scores(R2, f = L_zero_u2, method = "Bartlett")),
    class = "efa_scores_nonpositive_uniqueness"
  )

  expect_error(
    suppressMessages(efa_scores(R2, f = matrix("x", 2, 1))),
    class = "efa_scores_bad_loadings"
  )
  L_nonfinite <- matrix(c(.5, Inf), 2, 1)
  expect_error(
    suppressMessages(efa_scores(R2, f = L_nonfinite)),
    class = "efa_scores_bad_loadings"
  )
  expect_error(
    suppressMessages(efa_scores(matrix("x", 4, 2), f = L_zero_u2)),
    class = "efa_scores_bad_x"
  )
})

test_that("accepted correlation round-off is canonicalized", {
  R_near <- diag(3)
  R_near[1, 2] <- .2 + 1e-10
  R_near[2, 1] <- .2
  dimnames(R_near) <- list(letters[1:3], letters[1:3])
  out <- .align_correlation_axis(R_near, 3L, letters[1:3], "rho")
  expect_identical(out, t(out))
  expect_identical(unname(diag(out)), rep(1, 3))

  R_outside <- R_near
  R_outside[1, 2] <- R_outside[2, 1] <- 1 + 1e-12
  expect_error(
    .align_correlation_axis(R_outside, 3L, letters[1:3], "rho"),
    class = "efa_scores_matrix_correlation"
  )
})

test_that("component scores stay raw-scale while their diagnostics are standardized", {
  # method = "components" multiplies the raw, uncentered data by the loadings, but
  # r.scores / score_cor / determinacy are computed from W' R W, i.e. they describe
  # scale(x) %*% W. With unequal item variances the two diverge, which is what the
  # help page states; every other method returns standardized scores, where they agree.
  withr::local_seed(42)
  R <- test_models$baseline$cormat
  L <- unclass(efa_ob$rot_loadings)
  p <- ncol(R)
  x <- matrix(stats::rnorm(300 * p), 300, p) %*% chol(R)
  x <- sweep(x, 2L, seq(1, 10, length.out = p), "*")  # unequal variances
  x <- x + 10                                         # and non-zero means
  colnames(x) <- rownames(L)

  fs <- suppressMessages(suppressWarnings(
    efa_scores(x, f = L, Phi = efa_ob$Phi, method = "components")))

  expect_equal(fs$scores, x %*% fs$weights, ignore_attr = TRUE)
  # Uncentered. Taken on the absolute mean because a rotated factor carries no fixed
  # sign, so a column of the solution may be negated on another platform.
  expect_gt(min(abs(colMeans(fs$scores))), 1)
  expect_equal(unname(fs$r.scores),
               unname(stats::cor(scale(x) %*% fs$weights)), tolerance = 1e-8)
  expect_gt(max(abs(fs$r.scores - stats::cor(fs$scores))), 0.01)

  fs_b <- suppressMessages(suppressWarnings(
    efa_scores(x, f = L, Phi = efa_ob$Phi, method = "Bartlett")))
  expect_equal(unname(fs_b$r.scores), unname(stats::cor(fs_b$scores)),
               tolerance = 1e-8)
})

test_that("print and summary output are stable", {
  skip_on_cran()
  testthat::local_reproducible_output()

  # The printed values are eigen- and inverse-derived, so their last digits can drift
  # across BLAS implementations; scrub_num pins the layout, headers, and wording instead.

  # Brief print, raw-data scores.
  expect_snapshot(print(fs_grips), transform = scrub_num)

  # Brief print, correlation input (weights only).
  expect_snapshot(print(fs_ob_cor), transform = scrub_num)

  # Full summary: weight matrix, validity/univocality, and score intercorrelations.
  expect_snapshot(print(summary(fs_ob_cor)), transform = scrub_num)
})


# --- FACTOR_SCORES delegation: byte-identical regression against psych --------
#
# The rewritten FACTOR_SCORES() delegates to efa_scores(). Its weights/scores/
# r.scores must reproduce the original psych::factor.scores() pass-through
# bit-for-bit on the non-Heywood fixtures; only the R2 slot changes, now the
# corrected native determinacy^2 (which for oblique fits diverges from psych's
# rowSums(pattern^2) approximation but equals the textbook diag(S' R^-1 S)).
test_that("FACTOR_SCORES reproduces the original psych engine bit-for-bit", {
  skip_on_cran()
  skip_if_not_installed("psych")

  efa_raw <- suppressMessages(suppressWarnings(
    EFA(DOSPERT_raw, n_factors = 10, type = "EFAtools", method = "PAF",
        rotation = "oblimin")))

  fixtures <- list(
    list(x = DOSPERT_raw, f = efa_raw, method = "Bartlett"),
    list(x = test_models$baseline$cormat, f = efa_ob, method = "Thurstone"),
    list(x = test_models$baseline$cormat, f = efa_or, method = "Thurstone")
  )

  for (fx in fixtures) {
    lbl <- paste0(fx$method, "/", if (.is_cormat(fx$x)) "cor" else "raw")
    new <- suppressMessages(suppressWarnings(
      FACTOR_SCORES(fx$x, f = fx$f, method = fx$method)))

    rot <- fx$f$settings$rotation
    L <- unclass(if (!identical(rot, "none")) fx$f$rot_loadings else fx$f$unrot_loadings)
    ps <- psych::factor.scores(x = as.matrix(fx$x), f = L, Phi = fx$f$Phi,
                               method = fx$method)

    expect_equal(new$weights, ps$weights, tolerance = 1e-12,
                 ignore_attr = TRUE, info = lbl)
    expect_equal(new$r.scores, ps$r.scores, tolerance = 1e-12,
                 ignore_attr = TRUE, info = lbl)
    if (is.null(ps$scores)) {
      expect_null(new$scores, info = lbl)
    } else {
      expect_equal(new$scores, ps$scores, tolerance = 1e-12,
                   ignore_attr = TRUE, info = lbl)
    }

    # R2 is the corrected native determinacy^2; for regression it is diag(S' R^-1 S).
    es <- suppressMessages(suppressWarnings(
      efa_scores(fx$x, f = fx$f,
                 method = if (fx$method == "Thurstone") "regression" else fx$method)))
    expect_equal(unname(new$R2), unname(es$determinacy$rho2),
                 tolerance = 1e-12, info = lbl)
    if (fx$method == "Thurstone") {
      Phi <- if (is.null(fx$f$Phi)) diag(ncol(L)) else unclass(fx$f$Phi)
      S <- L %*% Phi
      Rr <- unclass(fx$f$orig_R)
      expect_equal(unname(new$R2), unname(diag(t(S) %*% solve(Rr, S))),
                   tolerance = 1e-10, info = lbl)
    }
  }
})

test_that("FACTOR_SCORES keeps its legacy output shape and Thurstone alias", {
  # Legacy names: raw carries `missing`, correlation input does not.
  fs_raw <- suppressMessages(suppressWarnings(
    FACTOR_SCORES(GRiPS_raw, f = efa_grips, method = "Bartlett")))
  fs_cor <- suppressMessages(FACTOR_SCORES(test_models$baseline$cormat, f = efa_ob))

  expect_named(fs_raw, c("scores", "weights", "r.scores", "missing", "R2", "settings"))
  expect_named(fs_cor, c("scores", "weights", "r.scores", "R2", "settings"))
  expect_s3_class(fs_raw, "FACTOR_SCORES")

  # "Thurstone" maps to efa_scores' "regression": identical weights.
  fs_thu <- suppressMessages(FACTOR_SCORES(test_models$baseline$cormat, f = efa_ob,
                                           method = "Thurstone"))
  es_reg <- suppressMessages(efa_scores(test_models$baseline$cormat, f = efa_ob,
                                        method = "regression"))
  expect_equal(fs_thu$weights, es_reg$weights, ignore_attr = TRUE)
})

test_that("the f error names the averaged matrix an efa_average object carries", {
  # An efa_average object is rejected while "a matrix" is listed as acceptable, and its
  # averaged loadings are exactly that. The message now says where they are.
  avg <- suppressWarnings(suppressMessages(
    efa_average(test_models$baseline$cormat, n_factors = 3, N = 500,
                estimator = c("PAF", "ULS"), rotation = "promax", type = "EFAtools",
                start_method = "psych", show_progress = FALSE)))

  e <- tryCatch(efa_scores(GRiPS_raw, f = avg), efa_scores_bad_f = function(e) e)
  expect_s3_class(e, "efa_scores_bad_f")
  expect_snapshot(cat(conditionMessage(e)))

  # The component the message names is accepted. Supplied on its own it carries no factor
  # correlations, which efa_scores() reports itself (efa_scores_phi_null), so the message
  # does not have to.
  expect_s3_class(
    suppressMessages(efa_scores(test_models$baseline$cormat, f = avg$loadings$average)),
    "efa_scores")
})

rm(efa_ob, efa_or, efa_grips, fs_ob_cor, fs_grips)
