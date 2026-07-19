# Case-insensitive choice matching for the efa_* interface: .match_arg_ci() resolves
# match.arg()-style choice arguments accepting any capitalization and unambiguous
# abbreviations, .map_subset_ci() maps vector-valued subset arguments. Both always
# return the canonical spelling from the choices, so stored settings, registry keys,
# and printed output are unaffected by how a value was capitalized.

test_that(".match_arg_ci returns the default when the argument is left at the formals", {
  f <- function(estimator = c("PAF", "ML", "ULS")) .match_arg_ci(estimator)
  expect_identical(f(), "PAF")

  # with several.ok, the full default resolves to all choices, as in match.arg()
  g <- function(gof = c("CAF", "CFI", "RMSEA")) .match_arg_ci(gof, several.ok = TRUE)
  expect_identical(g(), c("CAF", "CFI", "RMSEA"))
})

test_that(".match_arg_ci matches case-insensitively and returns the canonical spelling", {
  expect_identical(.match_arg_ci("ml", c("PAF", "ML", "ULS")), "ML")
  expect_identical(.match_arg_ci("Paf", c("PAF", "ML", "ULS")), "PAF")
  expect_identical(.match_arg_ci("GEOMINT", c("geominT", "geominQ")), "geominT")
  # exact spellings keep working unchanged
  expect_identical(.match_arg_ci("ULS", c("PAF", "ML", "ULS")), "ULS")
})

test_that(".match_arg_ci resolves abbreviations, with an exact match beating a prefix", {
  expect_identical(.match_arg_ci("pa", c("PAF", "ML", "ULS")), "PAF")
  expect_identical(.match_arg_ci("minr", c("PAF", "ML", "ULS", "MINRES")), "MINRES")
  # "efa" is both a full choice and a prefix of another one: the exact match wins
  expect_identical(.match_arg_ci("efa", c("EFAtools", "EFA")), "EFA")
})

test_that(".match_arg_ci with several.ok matches every element", {
  crits <- c("CD", "EKC", "HULL", "KGC", "PARALLEL", "SCREE", "SMT", "NEST", "MAP")
  expect_identical(.match_arg_ci(c("hull", "parallel"), crits, several.ok = TRUE),
                   c("HULL", "PARALLEL"))
})

test_that(".match_arg_ci rejects ambiguous and unmatched values with a classed error", {
  expect_error(.match_arg_ci("geomin", c("geominT", "geominQ")),
               class = "efa_bad_choice")
  expect_error(.match_arg_ci("bogus", c("PAF", "ML", "ULS")),
               class = "efa_bad_choice")
  # a bad element is an error even when the others match (stricter than match.arg(),
  # which silently drops it)
  expect_error(.match_arg_ci(c("HULL", "bogus"), c("HULL", "EKC"), several.ok = TRUE),
               class = "efa_bad_choice")
  # several values without several.ok, an empty vector, and a non-character input
  expect_error(.match_arg_ci(c("PAF", "ML"), c("PAF", "ML", "ULS")),
               class = "efa_bad_choice")
  expect_error(.match_arg_ci(character(), c("PAF", "ML"), several.ok = TRUE),
               class = "efa_bad_choice")
  expect_error(.match_arg_ci(1L, c("PAF", "ML")), class = "efa_bad_choice")
})

test_that(".match_arg_ci case-folds locale-independently and lets the caller set the class", {
  # ASCII folding, not toupper(): a value with an "i" must resolve regardless of the
  # locale's upper-casing rules (toupper("i") is the dotted capital under a Turkish locale).
  expect_identical(.match_arg_ci("minres", c("PAF", "ML", "ULS", "MINRES")), "MINRES")
  expect_identical(.match_arg_ci("oblimin", c("oblimin", "quartimin")), "oblimin")
  # the caller can raise its own condition class (the control constructors do)
  expect_error(.match_arg_ci("bogus", c("PAF", "ML"), class = "efa_control_input"),
               class = "efa_control_input")
})

test_that(".map_subset_ci maps onto the canonical spellings and leaves non-matches alone", {
  expect_identical(.map_subset_ci(c("paf", "MINRES"), c("PAF", "ML", "ULS", "MINRES")),
                   c("PAF", "MINRES"))
  # exact case-folded matches only: no abbreviations, and an unknown value passes
  # through unchanged for the caller's subset assertion to reject
  expect_identical(.map_subset_ci(c("pa", "bogus"), c("PAF", "ML")), c("pa", "bogus"))
  expect_identical(.map_subset_ci(TRUE, c("PAF", "ML")), TRUE)
})

test_that("efa_fit stores the canonical estimator and rotation from lowercase input", {
  fit <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
                 estimator = "ml", rotation = "geomint")
  expect_identical(fit$settings$estimator, "ML")
  # the settings keep a `method` entry duplicating `estimator`; it must stay in sync
  expect_identical(fit$settings$method, "ML")
  expect_identical(fit$settings$rotation, "geominT")
})

test_that("the retention criteria canonicalize case-insensitive choices (fast path)", {
  cmat <- test_models$baseline$cormat

  # efa_kgc's eigen_type (several.ok) and efa_ekc's type (several.ok) resolve without an
  # EFA fit, so these exercise the .match_arg_ci swap on a default run
  kgc <- efa_kgc(cmat, eigen_type = c("pca", "smc"))
  expect_identical(kgc$settings$eigen_type, c("PCA", "SMC"))

  ekc <- efa_ekc(cmat, N = 500, type = c("bva2017", "am2019"))
  expect_identical(ekc$settings$type, c("BvA2017", "AM2019"))
})

test_that("efa_retain matches criteria case-insensitively and keeps canonical names", {
  skip_if_not_slow()

  set.seed(500)
  out <- efa_retain(test_models$baseline$cormat, N = 500,
                    criteria = c("hull", "parallel"), gof = "caf",
                    eigen_type_other = "smc", suitability = FALSE,
                    show_progress = FALSE)
  expect_identical(out$settings$criteria, c("HULL", "PARALLEL"))
  expect_identical(out$settings$gof, "CAF")
  expect_identical(out$settings$eigen_type_other, "SMC")
  # the per-criterion results are keyed by the canonical uppercase names
  expect_true(any(grepl("^HULL", names(out$n_factors))))
  expect_true(any(grepl("^PARALLEL", names(out$n_factors))))
})

test_that("the control constructors match the type preset case-insensitively", {
  expect_identical(estimate_control(type = "spss")$type, "SPSS")
  expect_identical(rotate_control(type = "efatools")$type, "EFAtools")
})

test_that("efa_scores matches the factor-score method case-insensitively", {
  fit <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500)
  expect_message(
    fs <- efa_scores(test_models$baseline$cormat, f = fit, method = "bartlett"),
    class = "efa_scores_needs_raw"
  )
  expect_identical(fs$settings$method, "Bartlett")
})
