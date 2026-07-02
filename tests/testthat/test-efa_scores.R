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
  expect_equal(fs_grips$settings$n_obs, nrow(GRiPS_raw))
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

test_that("print and summary output are stable", {
  skip_on_cran()
  testthat::local_reproducible_output()

  # Brief print, raw-data scores.
  expect_snapshot(print(fs_grips))

  # Brief print, correlation input (weights only).
  expect_snapshot(print(fs_ob_cor))

  # Full summary: weight matrix, validity/univocality, and score intercorrelations.
  expect_snapshot(print(summary(fs_ob_cor)))
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

rm(efa_ob, efa_or, efa_grips, fs_ob_cor, fs_grips)
