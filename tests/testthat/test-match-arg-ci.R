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

# The simulation-side choice arguments. `return_pop = TRUE` stops after the population is
# built, which is all these need: the value is canonicalized before any data are drawn.
sim_Lambda <- population_models$loadings$baseline
sim_Phi <- population_models$phis_3$moderate

test_that("efa_simulate matches marginals case-insensitively", {
  sim <- efa_simulate(Lambda = sim_Lambda, Phi = sim_Phi,
                      marginals = "vm", return_pop = TRUE)
  expect_identical(sim$settings$marginals, "VM")
  expect_error(efa_simulate(Lambda = sim_Lambda, marginals = "bogus", return_pop = TRUE),
               class = "efa_bad_choice")
})

test_that("efa_simulate matches the model-error method case-insensitively", {
  sim <- efa_simulate(Lambda = sim_Lambda, Phi = sim_Phi,
                      model_error = "tkl", target_rmsea = 0.05, return_pop = TRUE)
  expect_identical(sim$model_error$method, "TKL")
  expect_error(efa_simulate(Lambda = sim_Lambda, model_error = "bogus",
                            target_rmsea = 0.05, return_pop = TRUE),
               class = "efa_bad_choice")
})

test_that("efa_simulate matches the missing-data mechanism case-insensitively", {
  sim <- efa_simulate(Lambda = sim_Lambda, Phi = sim_Phi,
                      missing = "mcar", missing_prop = 0.1, return_pop = TRUE)
  expect_identical(sim$settings$missing, "MCAR")
  expect_error(efa_simulate(Lambda = sim_Lambda, missing = "bogus",
                            missing_prop = 0.1, return_pop = TRUE),
               class = "efa_bad_choice")
})

test_that("efa_simulate matches the threshold-matching rule case-insensitively", {
  sim <- efa_simulate(Lambda = sim_Lambda, Phi = sim_Phi, categories = 4,
                      match = "Polychoric", return_pop = TRUE)
  expect_identical(sim$settings$match, "polychoric")
  expect_error(efa_simulate(Lambda = sim_Lambda, categories = 4, match = "bogus",
                            return_pop = TRUE),
               class = "efa_bad_choice")
})

test_that("efa_power matches mode and type case-insensitively", {
  # the analytic path resolves both without any simulation
  pw <- efa_power(mode = "RMSEA", type = "NotClose", p = 18, k = 3, N = 200)
  expect_identical(pw$settings$mode, "rmsea")
  expect_identical(pw$settings$type, "notclose")
  expect_error(efa_power(mode = "bogus"), class = "efa_bad_choice")
  expect_error(efa_power(type = "bogus", p = 18, k = 3, N = 200),
               class = "efa_bad_choice")
})

test_that("efa_power matches the model-error method case-insensitively", {
  pw <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 150,
                  n_datasets = 3, criteria = "EKC", model_error = "tkl",
                  target_rmsea = 0.05, seed = 42)
  expect_identical(pw$settings$model_error, "TKL")
  expect_error(efa_power("simulation", Lambda = sim_Lambda, N = 150, n_datasets = 3,
                         criteria = "EKC", model_error = "bogus", target_rmsea = 0.05),
               class = "efa_bad_choice")
})

test_that("efa_mi matches the fit-pooling method case-insensitively", {
  cmats <- list(test_models$baseline$cormat, test_models$baseline$cormat)
  pooled <- efa_mi(cmats, n_factors = 3, N = 500, estimator = "PAF",
                   rotation = "promax", fit_pool_method = "d2")
  expect_identical(pooled$settings$fit_pool_method, "D2")
  expect_error(efa_mi(cmats, n_factors = 3, N = 500, fit_pool_method = "bogus"),
               class = "efa_bad_choice")
})

test_that("efa_scores matches the factor-score method case-insensitively", {
  fit <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500)
  expect_message(
    fs <- efa_scores(test_models$baseline$cormat, f = fit, method = "bartlett"),
    class = "efa_scores_needs_raw"
  )
  expect_identical(fs$settings$method, "Bartlett")
})

# --------------------------------------------------------------------------------------
# The rule, enumerated
# --------------------------------------------------------------------------------------
# One table listing every choice-valued argument a user can supply on the efa_* interface.
# It is the rule in executable form: a choice argument that is missing from it, or that is
# present and does not reject an unmatched value with the interface's condition class, is
# what this block exists to catch. Rejection happens in each function's argument-check
# block, before anything is computed, so the whole table runs in well under a second.

.choice_arg_cases <- function() {
  cmat <- test_models$baseline$cormat
  fit <- efa_fit(cmat, n_factors = 3, N = 500, rotation = "promax")
  loads <- unclass(fit$rot_loadings)
  case <- function(fn, args, ..., class = "efa_bad_choice") {
    lapply(c(...), function(a) list(fn = fn, args = args, arg = a, class = class))
  }
  c(
    case("efa_fit", list(cmat, n_factors = 3, N = 500), "se", "use", "cor_method"),
    case("efa_bartlett", list(cmat, N = 500), "use", "cor_method"),
    case("efa_kmo", list(cmat), "use", "cor_method"),
    case("efa_map", list(cmat), "use", "cor_method"),
    case("efa_smt", list(cmat, N = 500), "use", "cor_method"),
    case("efa_ekc", list(cmat, N = 500), "use", "cor_method", "type"),
    case("efa_kgc", list(cmat), "use", "cor_method", "eigen_type"),
    case("efa_scree", list(cmat), "use", "cor_method", "eigen_type"),
    case("efa_hull", list(cmat, N = 500),
         "use", "cor_method", "decision_rule", "estimator", "gof", "eigen_type"),
    case("efa_parallel", list(cmat, N = 500, n_datasets = 5),
         "use", "cor_method", "decision_rule", "eigen_type"),
    case("efa_nest", list(cmat, N = 500, n_datasets = 5), "use", "cor_method"),
    case("efa_cd", list(GRiPS_raw[1:50, ]), "cor_method"),
    case("efa_screen", list(cmat, N = 500), "use", "cor_method"),
    case("efa_retain", list(cmat, N = 500, suitability = FALSE),
         "use", "cor_method", "decision_rule", "criteria", "estimator", "gof",
         "eigen_type_HULL", "eigen_type_other", "ekc_type"),
    case("efa_average", list(cmat, n_factors = 3, N = 500, show_progress = FALSE),
         "use", "cor_method", "averaging"),
    case("efa_compare", list(loads, loads), "reorder"),
    case("efa_reliability", list(fit), "variance"),
    case("efa_procrustes", list(loads, loads), "rotation"),
    case("efa_mi", list(list(cmat, cmat), n_factors = 3, N = 500),
         "target_method", "align_unrotated", "fit_pool_method"),
    case("efa_scores", list(cmat, f = fit), "method"),
    case("efa_power", list(p = 18, k = 3, N = 200), "mode", "type"),
    # the simulation branch validates its own choices before it draws anything
    case("efa_power", list("simulation"), "criteria", "estimator", "model_error"),
    case("estimate_control", list(),
         "type", "init_comm", "criterion_type", "start_method",
         class = "efa_control_input"),
    case("rotate_control", list(),
         "type", "order_type", "varimax_type", "p_type",
         class = "efa_control_input"),
    # the print and plot layer, which takes its own choice arguments
    case("format", list(fit), "sort_loadings"),
    case("summary", list(fit), "ci", "ci_filter", "sort_loadings"),
    # `view` is not listed here: it is set internally ("brief" by format.efa(), "full" by
    # format.summary.efa()) and never reaches the renderer from a user's call
    case("print", list(summary(fit)), "ci", "ci_filter", "sort_loadings"),
    case("format", list(fit$rot_loadings), "name_style", "sort_loadings"),
    case("residuals", list(fit), "type")
  )
}

test_that("every choice-valued argument rejects an unmatched value with a classed error", {
  for (case in .choice_arg_cases()) {
    expect_error(
      do.call(case$fn, c(case$args, stats::setNames(list("bogus"), case$arg))),
      class = case$class,
      info = paste0(case$fn, "(", case$arg, " = \"bogus\")")
    )
  }
})

test_that("efa_average's grid arguments are case-insensitive but keep the subset check", {
  # These take a whole vector of values to build the averaging grid from, so they are
  # validated by checkmate::assert_subset() rather than by .match_arg_ci(): the value is
  # mapped onto the canonical spellings first, and an element matching none of them falls
  # through to the subset assertion, whose message already names the argument and lists the
  # admissible values. The condition is therefore deliberately unclassed here.
  cmat <- test_models$baseline$cormat
  avg <- efa_average(cmat, n_factors = 3, N = 500, estimator = c("paf", "PAF"),
                     rotation = c("Promax", "Oblimin"), type = "psych", init_comm = "SMC",
                     criterion_type = "SUM", varimax_type = "SVD", p_type = "NORM",
                     start_method = "Psych", show_progress = FALSE)
  # "paf" and "PAF" fold onto one grid entry, so the grid runs each requested value once
  expect_identical(avg$settings$estimator, "PAF")
  expect_identical(avg$settings$rotation, c("promax", "oblimin"))

  for (arg in c("estimator", "rotation", "type", "init_comm", "criterion_type",
                "varimax_type", "p_type", "start_method")) {
    expect_error(
      do.call(efa_average, c(list(cmat, n_factors = 3, N = 500, show_progress = FALSE),
                             stats::setNames(list("bogus"), arg))),
      info = arg
    )
  }
})

test_that("plot.efa_group rejects an unmatched type with a classed error", {
  mg <- efa_group(GRiPS_raw[1:100, ],
                  groups = rep(c("g1", "g2"), length.out = 100), n_factors = 1)
  expect_error(plot(mg, type = "bogus"), class = "efa_bad_choice")
  expect_s3_class(plot(mg, type = "Congruence"), "ggplot")
})

test_that("the correlation and inference arguments canonicalize a flipped spelling", {
  cmat <- test_models$baseline$cormat

  fit <- efa_fit(cmat, n_factors = 3, N = 500, estimator = "ML",
                 cor_method = "Pearson", use = "Complete.Obs", se = "Information")
  expect_identical(fit$settings$cor_method, "pearson")
  expect_identical(fit$settings$use, "complete.obs")
  expect_identical(fit$settings$se, "information")

  expect_identical(efa_kmo(cmat, cor_method = "Kendall")$settings$cor_method, "kendall")
  # efa_hull() does not record `decision_rule`, so the flipped spelling is pinned by the
  # result it produces being the one the canonical spelling produces
  expect_identical(efa_parallel(cmat, N = 500, n_datasets = 20, decision_rule = "Crawford",
                                eigen_type = "PCA")$settings$decision_rule, "crawford")
  expect_identical(efa_average(cmat, n_factors = 3, N = 500, rotation = "none",
                               type = "none", averaging = "Median",
                               show_progress = FALSE)$settings$averaging, "median")

  # the tuning knobs, which reach the constructors from the flat interface as well
  expect_identical(estimate_control(init_comm = "SMC")$init_comm, "smc")
  expect_identical(estimate_control(criterion_type = "SUM")$criterion_type, "sum")
  expect_identical(estimate_control(start_method = "Factanal")$start_method, "factanal")
  expect_identical(rotate_control(order_type = "Eigen")$order_type, "eigen")
  expect_identical(rotate_control(varimax_type = "SVD")$varimax_type, "svd")
  expect_identical(rotate_control(p_type = "NORM")$p_type, "norm")

  # the print layer returns text, so the flipped value is asserted by it not erroring
  expect_identical(format(fit, sort_loadings = "Primary"),
                   format(fit, sort_loadings = "primary"))
  expect_identical(residuals(fit, type = "Raw"), residuals(fit, type = "raw"))
})

test_that("the superseded uppercase functions keep matching their choices exactly", {
  # These four are frozen on the pre-rename behaviour (see R/EFAtools-superseded.R), so a
  # flipped spelling must still be refused there even though the successors accept it.
  cmat <- test_models$baseline$cormat
  efa <- EFA(cmat, n_factors = 3, N = 500, rotation = "promax")

  # A bare expect_error() would also be satisfied by the successors' own condition, which is
  # exactly what must not appear here, so assert the class the frozen path may never raise.
  refuses_exactly <- function(expr) {
    cnd <- tryCatch({expr; NULL}, error = function(e) e)
    !is.null(cnd) && !inherits(cnd, "efa_bad_choice")
  }

  expect_true(refuses_exactly(EFA(cmat, n_factors = 3, N = 500, type = "efatools")))
  expect_true(refuses_exactly(EFA(cmat, n_factors = 3, N = 500, rotation = "Promax")))
  expect_true(refuses_exactly(SL(efa, n_factors = 3, type = "efatools")))
  expect_true(refuses_exactly(OMEGA(efa, type = "efatools")))
  expect_true(refuses_exactly(FACTOR_SCORES(cmat, f = efa, method = "Regression")))
})

# --------------------------------------------------------------------------------------
# NULL selects the documented default
# --------------------------------------------------------------------------------------
# .match_arg_ci() resolves a NULL argument to the caller's formal default, not to the first
# admissible value. The two differ wherever an argument takes several values at once, or
# wherever its choices are listed in another order than its default, and taking the first
# choice there silently ran something narrower than the help page documents.

test_that(".match_arg_ci resolves NULL to the caller's formal default", {
  f <- function(estimator = c("ML", "PAF", "ULS")) .match_arg_ci(estimator)
  expect_identical(f(NULL), "ML")

  # several.ok: the whole default vector, not its first element
  g <- function(gof = c("CAF", "CFI", "RMSEA")) .match_arg_ci(gof, several.ok = TRUE)
  expect_identical(g(NULL), c("CAF", "CFI", "RMSEA"))

  # an explicit `choices` ordered differently from the default
  h <- function(eigen_type = "SMC") {
    .match_arg_ci(eigen_type, c("PCA", "SMC", "EFA"), several.ok = TRUE)
  }
  expect_identical(h(NULL), "SMC")

  # a formal that is itself NULL, or has no default at all, falls back to the first choice
  k <- function(match = NULL) .match_arg_ci(match, c("thresholds", "polychoric"))
  expect_identical(k(NULL), "thresholds")
  m <- function(match) .match_arg_ci(match, c("thresholds", "polychoric"))
  expect_identical(m(NULL), "thresholds")
})

test_that("NULL selects the documented default at every affected call site", {
  cmat <- test_models$baseline$cormat

  expect_identical(efa_kgc(cmat, eigen_type = NULL)$settings$eigen_type,
                   efa_kgc(cmat)$settings$eigen_type)
  expect_identical(efa_scree(cmat, eigen_type = NULL)$settings$eigen_type,
                   efa_scree(cmat)$settings$eigen_type)
  expect_identical(efa_parallel(cmat, N = 500, n_datasets = 5,
                                eigen_type = NULL)$settings$eigen_type,
                   c("PCA", "SMC", "EFA"))
  expect_identical(efa_hull(cmat, N = 500, estimator = "ML", gof = NULL)$settings$gof,
                   c("CAF", "CFI", "RMSEA"))
  expect_identical(
    efa_retain(cmat, N = 500, criteria = "KGC", suitability = FALSE,
               eigen_type_other = NULL)$settings$eigen_type_other,
    "SMC")
  expect_identical(
    efa_retain(cmat, N = 500, criteria = "HULL", suitability = FALSE, estimator = "ML",
               gof = NULL)$settings$gof,
    c("CAF", "CFI", "RMSEA"))

  # efa_simulate()'s `match` documents NULL as "thresholds"; that must keep working
  sim <- efa_simulate(Lambda = sim_Lambda, Phi = sim_Phi, categories = 4,
                      match = NULL, return_pop = TRUE)
  expect_identical(sim$settings$match, "thresholds")
})

test_that("efa_retain(criteria = NULL) runs the six documented default criteria", {
  skip_if_not_slow()

  set.seed(500)
  # CD is one of the six defaults and needs raw data, so a correlation matrix skips it with
  # a warning -- which is itself evidence that all six were requested
  expect_warning(
    out <- efa_retain(test_models$baseline$cormat, N = 500, criteria = NULL,
                      suitability = FALSE, n_datasets = 50, n_datasets_nest = 50,
                      show_progress = FALSE),
    class = "efa_criterion_skipped"
  )
  expect_identical(out$settings$criteria,
                   eval(formals(efa_retain)$criteria))
})

test_that("efa_power(criteria = NULL) evaluates the documented default criteria", {
  skip_if_not_slow()

  pw <- efa_power("simulation", Lambda = sim_Lambda, Phi = sim_Phi, N = 150,
                  n_datasets = 3, criteria = NULL, seed = 42)
  expect_identical(pw$settings$criteria, eval(formals(efa_power)$criteria))
})
