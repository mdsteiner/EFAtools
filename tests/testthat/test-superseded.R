# The superseded UPPERCASE names must be thin, silent bridges to their efa_*
# implementations: identical result, no added condition, and the old class string kept so
# `inherits(x, "OLD")` and old S3 dispatch still resolve.

# `EFA()` is the *translating* wrapper of record: efa_fit() collects the estimation and
# rotation tuning knobs into estimate_control() / rotate_control(), so EFA() repacks its frozen
# flat argument list into those objects rather than forwarding it verbatim. It is therefore not
# part of the static `superseded_contract` below; its contract is that the old flat call and
# the new control-object call return byte-identical objects, across several type / method /
# rotation combinations, and that EFA() adds no condition of its own. A shared seed pins the
# random-start rotations so the two calls draw from the same RNG state.
#
# The wrappers whose `...` reaches a fit (PARALLEL, KGC, HULL, SCREE, NEST, SL, N_FACTORS,
# EFA_POOLED) translate for the same reason -- the flat knobs a legacy call passes through
# their dots would otherwise match no efa_fit() formal and be dropped -- so the static contract
# pins only their frozen signature and the knob-forwarding tests at the end of this file pin
# their behaviour.

test_that("EFA() and efa_fit() return byte-identical objects", {
  cormat <- test_models$baseline$cormat
  combos <- expand.grid(
    method = c("PAF", "ML", "ULS"),
    rotation = c("none", "varimax", "promax", "oblimin"),
    type = c("EFAtools", "psych", "SPSS"),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(combos))) {
    m <- combos$method[i]
    r <- combos$rotation[i]
    ty <- combos$type[i]
    set.seed(42L)
    old <- suppressWarnings(EFA(cormat, n_factors = 3, N = 500, method = m,
                                rotation = r, type = ty))
    set.seed(42L)
    new <- suppressWarnings(efa_fit(cormat, n_factors = 3, N = 500, method = m,
                                    rotation = r,
                                    estimate_control = estimate_control(type = ty),
                                    rotate_control = rotate_control(type = ty)))
    expect_identical(old, new,
                     info = paste("method =", m, "; rotation =", r, "; type =", ty))
  }
})

test_that("efa_fit() carries the new classes and its methods dispatch", {
  mod <- efa_fit(test_models$baseline$cormat, n_factors = 3, N = 500,
                 rotation = "promax")

  expect_identical(class(mod), c("efa", "EFA"))
  expect_identical(class(mod$unrot_loadings), c("efa_loadings", "LOADINGS"))
  expect_identical(class(mod$rot_loadings), c("efa_loadings", "LOADINGS"))
  expect_true(inherits(mod, "EFA"))

  # the renamed S3 methods resolve on the new class
  expect_s3_class(summary(mod), "summary.efa")
  expect_type(residuals(mod), "double")
  expect_no_error(invisible(format(mod)))
})

test_that("efa_fit() rejects a former flat knob forwarded through `...`, keeps rotation extras", {
  cm <- test_models$baseline$cormat

  # A knob that moved into the control objects (`k`, `type`, `max_iter`, ...) is no longer an
  # efa_fit() formal. Left to reach `.efa_core()` -- which still declares those names -- it
  # would supply the formal twice and abort with "matched by multiple actual arguments"; left
  # to be dropped, it would quietly run the default preset instead of the requested one.
  # Neither is acceptable, so it is rejected up front with a classed error that names the
  # control object owning the knob.
  expect_error(efa_fit(cm, n_factors = 3, N = 500, rotation = "promax", k = 4),
               class = "efa_flat_knob_in_dots")

  # a caller that forwards its user `...` (here efa_group) surfaces the same error rather
  # than silently ignoring the knob
  expect_error(
    suppressWarnings(efa_group(list(a = cm, b = cm), n_factors = 2, type = "SPSS")),
    class = "efa_flat_knob_in_dots"
  )
  # and the control objects ride through those same dots to reach the fit
  expect_s3_class(
    suppressWarnings(efa_group(list(a = cm, b = cm), n_factors = 2,
                               estimate_control = estimate_control(type = "SPSS"),
                               rotate_control = rotate_control(type = "SPSS"))),
    "efa_group"
  )

  # Genuine rotation extras still reach the rotation engine, via `...` or via rotate_control().
  via_dots <- efa_fit(cm, n_factors = 3, N = 500, rotation = "geominQ", maxit = 3000)
  via_ctrl <- efa_fit(cm, n_factors = 3, N = 500, rotation = "geominQ",
                      rotate_control = rotate_control(maxit = 3000))
  expect_identical(via_dots$rot_loadings, via_ctrl$rot_loadings)
})

test_that("EFA() stays transparent to efa_fit()'s conditions", {
  # no condition of its own on a clean correlation-matrix path
  expect_no_condition(EFA(test_models$baseline$cormat, n_factors = 3, N = 500))

  # it neither swallows nor adds to efa_fit()'s conditions: the raw-data path still
  # signals that the correlations were computed from the data
  expect_message(EFA(GRiPS_raw, n_factors = 1), class = "efa_cor_from_data")
})

test_that("BARTLETT() forwards to efa_bartlett() identically", {
  old <- BARTLETT(test_models$baseline$cormat, N = 500)
  new <- efa_bartlett(test_models$baseline$cormat, N = 500)

  expect_identical(old, new)
  expect_identical(class(old), c("efa_bartlett", "BARTLETT"))
  expect_true(inherits(old, "BARTLETT"))
})

test_that("BARTLETT() adds no condition and stays transparent", {
  # no deprecation signal (or any other condition) on a clean path
  expect_no_condition(BARTLETT(test_models$baseline$cormat, N = 500))

  # but it neither swallows nor adds to the conditions efa_bartlett() itself
  # raises: the raw-data path still signals that the correlations were computed
  expect_message(BARTLETT(GRiPS_raw), class = "efa_cor_from_data")
})

test_that("KMO() forwards to efa_kmo() identically", {
  old <- KMO(test_models$baseline$cormat)
  new <- efa_kmo(test_models$baseline$cormat)

  expect_identical(old, new)
  expect_identical(class(old), c("efa_kmo", "KMO"))
  expect_true(inherits(old, "KMO"))
})

test_that("KMO() adds no condition and stays transparent", {
  expect_no_condition(KMO(test_models$baseline$cormat))
  expect_message(KMO(GRiPS_raw), class = "efa_cor_from_data")
})

# The retention criteria share the lowercase `efa_retention` class and gain no
# extra class string on the old name, so the wrapper must return a byte-identical
# object with that single class and raise no condition of its own. EKC and MAP
# stand in for the family: both are deterministic, so identity is exact.

test_that("EKC() forwards to efa_ekc() identically", {
  old <- EKC(test_models$baseline$cormat, N = 500)
  new <- efa_ekc(test_models$baseline$cormat, N = 500)

  expect_identical(old, new)
  expect_identical(class(old), "efa_retention")
})

test_that("EKC() adds no condition and stays transparent", {
  expect_no_condition(EKC(test_models$baseline$cormat, N = 500))
  expect_message(EKC(GRiPS_raw), class = "efa_cor_from_data")
})

test_that("MAP() forwards to efa_map() identically", {
  old <- MAP(test_models$baseline$cormat)
  new <- efa_map(test_models$baseline$cormat)

  expect_identical(old, new)
  expect_identical(class(old), "efa_retention")
})

test_that("MAP() adds no condition and stays transparent", {
  expect_no_condition(MAP(test_models$baseline$cormat))
  expect_message(MAP(GRiPS_raw), class = "efa_cor_from_data")
})

# The retention orchestrator keeps the dual class its efa_retain() successor
# carries (so both `inherits(x, "N_FACTORS")` and the new S3 dispatch resolve),
# forwards byte-identically, and raises no condition of its own. EKC and MAP
# stand in for the criteria it runs: both are deterministic, so identity is exact.

test_that("N_FACTORS() forwards to efa_retain() identically", {
  old <- N_FACTORS(test_models$baseline$cormat, N = 500, suitability = FALSE,
                   criteria = c("EKC", "MAP"))
  new <- efa_retain(test_models$baseline$cormat, N = 500, suitability = FALSE,
                    criteria = c("EKC", "MAP"))

  expect_identical(old, new)
  expect_identical(class(old), c("efa_retain", "N_FACTORS"))
  expect_true(inherits(old, "N_FACTORS"))
})

test_that("N_FACTORS() adds no condition and stays transparent", {
  expect_no_condition(N_FACTORS(test_models$baseline$cormat, N = 500,
                                suitability = FALSE, criteria = c("EKC", "MAP")))
})

# The Schmid-Leiman transformation keeps the dual class on the returned list and
# on the nested loading matrix, so `inherits(x, "SL")`, `inherits(x$sl,
# "SLLOADINGS")`, and the OMEGA() / efa_reliability() input branches keyed on them
# all still resolve. An oblique EFA is the input every Schmid-Leiman path needs.

oblique_efa <- function() {
  EFA(test_models$baseline$cormat, N = 500, n_factors = 3,
      type = "EFAtools", method = "PAF", rotation = "promax")
}

test_that("SL() forwards to efa_schmid_leiman() identically", {
  efa_mod <- oblique_efa()
  old <- SL(efa_mod, type = "EFAtools", method = "PAF")
  new <- efa_schmid_leiman(efa_mod, method = "PAF",
                           estimate_control = estimate_control(type = "EFAtools"))

  expect_identical(old, new)
  expect_identical(class(old), c("efa_schmid_leiman", "SL"))
  expect_identical(class(old$sl), c("efa_sl_loadings", "SLLOADINGS"))
})

test_that("SL() adds no condition and stays transparent", {
  efa_mod <- oblique_efa()

  expect_no_condition(SL(efa_mod, type = "EFAtools", method = "PAF"))

  # it neither swallows nor adds to the conditions efa_schmid_leiman() raises:
  # supplying Phi alongside an EFA object still warns
  expect_warning(SL(efa_mod, type = "EFAtools", method = "PAF", Phi = efa_mod$Phi),
                 class = "efa_sl_phi_specified")
})

# OMEGA() and FACTOR_SCORES() carry the superseded badge but keep their own
# implementation, output shape, and class. They must gain no runtime signal.

test_that("OMEGA() adds no condition and keeps its class", {
  sl_mod <- efa_schmid_leiman(oblique_efa(), method = "PAF",
                              estimate_control = estimate_control(type = "EFAtools"))
  fc <- sl_mod$sl[, c("F1", "F2", "F3")] >= .2

  expect_no_condition(om <- OMEGA(sl_mod, type = "EFAtools", factor_corres = fc))
  expect_s3_class(om, "OMEGA")
})

test_that("FACTOR_SCORES() adds no condition and keeps its class", {
  efa_mod <- oblique_efa()

  # the correlation-matrix path signals that only weights can be returned; the
  # wrapper neither swallows that message nor adds one of its own
  expect_message(FACTOR_SCORES(test_models$baseline$cormat, f = efa_mod),
                 class = "efa_scores_needs_raw")

  fs <- suppressMessages(FACTOR_SCORES(test_models$baseline$cormat, f = efa_mod))
  expect_s3_class(fs, "FACTOR_SCORES")
})

# The comparison helper keeps the dual class its efa_compare() successor carries, so
# `inherits(x, "COMPARE")` and the old print / format / plot dispatch all still resolve.

test_that("COMPARE() forwards to efa_compare() identically", {
  old <- COMPARE(matrix(c(1, 1, 1, 2), ncol = 2), matrix(c(1, 1, 1, 1), ncol = 2))
  new <- efa_compare(matrix(c(1, 1, 1, 2), ncol = 2), matrix(c(1, 1, 1, 1), ncol = 2))

  expect_identical(old, new)
  expect_identical(class(old), c("efa_compare", "COMPARE"))
  expect_true(inherits(old, "COMPARE"))
})

test_that("COMPARE() adds no condition and stays transparent", {
  expect_no_condition(COMPARE(1:10, 1:10))

  # it neither swallows nor adds to the conditions efa_compare() itself raises: a named
  # vector under the default congruence reordering still warns that it cannot be reordered
  expect_warning(COMPARE(c(a = 1, b = 2, c = 4), c(a = 1, b = 2, c = 4)),
                 class = "efa_compare_reorder_vectors")
})

# The Procrustes alignment returns a plain, unclassed list, so the old name has no class
# string to preserve -- only the frozen signature and the silent pass-through.

# A fixed loading matrix and a quarter-turn rotation of it. The default orthogonal solver
# recovers that rotation exactly whatever the loadings are, so the two tests below need no
# random draw (and leave the RNG stream untouched for whatever runs next).
procrustes_A <- matrix(c(0.8, 0.7, 0.6, 0.1, 0.2, 0.1,
                         0.1, 0.2, 0.1, 0.8, 0.7, 0.6), ncol = 2)
procrustes_target <- procrustes_A %*% matrix(c(0, 1, -1, 0), ncol = 2)

test_that("PROCRUSTES() forwards to efa_procrustes() identically", {
  old <- PROCRUSTES(procrustes_A, Target = procrustes_target)
  new <- efa_procrustes(procrustes_A, Target = procrustes_target)

  expect_identical(old, new)
  expect_identical(class(old), "list")
})

test_that("PROCRUSTES() adds no condition and stays transparent", {
  expect_no_condition(PROCRUSTES(procrustes_A, Target = procrustes_target))

  # the dimension check efa_procrustes() performs still reaches the caller
  expect_error(PROCRUSTES(procrustes_A, Target = procrustes_target[, 1, drop = FALSE]),
               class = "efa_dim_mismatch")
})

# Model averaging keeps the dual class its efa_average() successor carries, so
# `inherits(x, "EFA_AVERAGE")` and the old print / format / plot dispatch all still
# resolve. Two grid cells (PAF and ML, both the EFAtools implementation) are enough to
# average over, and every EFA in them is deterministic, so identity is exact.

average_grid <- function(f) {
  f(test_models$baseline$cormat, n_factors = 3, N = 500,
    method = c("PAF", "ML"), type = "EFAtools", start_method = "psych",
    show_progress = FALSE)
}

# Averaging a grid is the most expensive wrapper check in this file, so the silence
# assertion rides along with the identity one rather than fitting the grid a third time.
test_that("EFA_AVERAGE() forwards to efa_average() identically and adds no condition", {
  expect_no_condition(old <- average_grid(EFA_AVERAGE))
  new <- average_grid(efa_average)

  expect_identical(old, new)
  expect_identical(class(old), c("efa_average", "EFA_AVERAGE"))
  expect_true(inherits(old, "EFA_AVERAGE"))
})

test_that("EFA_AVERAGE() stays transparent to efa_average()'s conditions", {
  # it neither swallows nor adds to the conditions efa_average() itself raises: a grid
  # that collapses to a single EFA still warns and returns that EFA
  expect_warning(EFA_AVERAGE(test_models$baseline$cormat, n_factors = 3, N = 500,
                             method = "PAF", type = "EFAtools", show_progress = FALSE),
                 class = "efa_avg_single_combination")
})

# EFA_POOLED() keeps the dual class its efa_mi() successor carries (so both
# `inherits(x, "EFA_POOLED")` and the new S3 dispatch resolve), forwards
# byte-identically, and raises no condition of its own. A two-imputation list of
# identical correlation matrices is deterministic (se = "none"), so identity is exact.

test_that("EFA_POOLED() forwards to efa_mi() identically and adds no condition", {
  imps <- list(test_models$baseline$cormat, test_models$baseline$cormat)

  expect_no_condition(
    old <- EFA_POOLED(imps, n_factors = 3, N = 500, method = "PAF", rotation = "promax")
  )
  new <- efa_mi(imps, n_factors = 3, N = 500, method = "PAF", rotation = "promax")

  expect_identical(old, new)
  expect_identical(class(old), c("efa_mi", "EFA_POOLED", "efa", "EFA"))
  expect_true(inherits(old, "EFA_POOLED"))
})

# Each old name is frozen: its argument list may never change. That half is checked
# statically for all seventeen names below.
#
# A name marked `translating` also repacks the old flat tuning knobs its `...` may carry into
# the control objects the fit now takes (`.repack_flat_dots()`), so it is not a one-line
# forwarder and only its signature is pinned here; that its knobs still reach the fit is
# pinned behaviourally further down. Every other name forwards each of its arguments,
# unchanged, to the identically named parameter of its efa_* successor -- checked statically on
# the wrapper's own body, so a dropped, renamed, or mis-mapped argument is caught even for the
# simulation-heavy criteria that are too expensive to run twice in an identity test.
superseded_contract <- list(
  BARTLETT = list(
    callee = "efa_bartlett",
    formals = c("x", "N", "use", "cor_method")),
  KMO = list(
    callee = "efa_kmo",
    formals = c("x", "use", "cor_method")),
  PARALLEL = list(
    callee = "efa_parallel",
    translating = TRUE,
    formals = c("x", "N", "n_vars", "n_datasets", "percent", "eigen_type", "use",
                "cor_method", "decision_rule", "n_factors", "...")),
  EKC = list(
    callee = "efa_ekc",
    formals = c("x", "N", "use", "cor_method", "type")),
  KGC = list(
    callee = "efa_kgc",
    translating = TRUE,
    formals = c("x", "eigen_type", "use", "cor_method", "n_factors", "...")),
  HULL = list(
    callee = "efa_hull",
    translating = TRUE,
    formals = c("x", "N", "n_fac_theor", "method", "gof", "eigen_type", "use",
                "cor_method", "n_datasets", "percent", "decision_rule",
                "n_factors", "...")),
  SCREE = list(
    callee = "efa_scree",
    translating = TRUE,
    formals = c("x", "eigen_type", "use", "cor_method", "n_factors", "...")),
  MAP = list(
    callee = "efa_map",
    formals = c("x", "use", "cor_method")),
  NEST = list(
    callee = "efa_nest",
    translating = TRUE,
    formals = c("x", "N", "alpha", "use", "cor_method", "n_datasets", "...")),
  SMT = list(
    callee = "efa_smt",
    formals = c("x", "N", "use", "cor_method")),
  CD = list(
    callee = "efa_cd",
    formals = c("x", "n_factors_max", "N_pop", "N_samples", "alpha", "cor_method",
                "max_iter")),
  N_FACTORS = list(
    callee = "efa_retain",
    translating = TRUE,
    formals = c("x", "criteria", "suitability", "N", "use", "cor_method",
                "n_factors_max", "N_pop", "N_samples", "alpha", "max_iter_CD",
                "n_fac_theor", "method", "gof", "eigen_type_HULL",
                "eigen_type_other", "n_factors", "n_datasets", "percent",
                "decision_rule", "ekc_type", "n_datasets_nest", "alpha_nest",
                "show_progress", "...")),
  SL = list(
    callee = "efa_schmid_leiman",
    # `type` is frozen here but is no longer a bare formal of the successor: it is repacked
    # into the estimation control alongside any flat knob the dots carry.
    translating = TRUE,
    formals = c("x", "Phi", "type", "method", "g_name", "...")),
  COMPARE = list(
    callee = "efa_compare",
    formals = c("x", "y", "reorder", "corres", "thresh", "digits", "m_red",
                "range_red", "round_red", "print_diff", "na.rm", "x_labels",
                "plot", "plot_red")),
  EFA_AVERAGE = list(
    callee = "efa_average",
    # the frozen `P_type` formal forwards to the successor's renamed `p_type`.
    renamed = c(P_type = "p_type"),
    formals = c("x", "n_factors", "N", "method", "rotation", "type", "averaging",
                "trim", "salience_threshold", "max_iter", "init_comm", "criterion",
                "criterion_type", "abs_eigen", "varimax_type", "normalize",
                "k_promax", "k_simplimax", "P_type", "precision", "start_method",
                "use", "cor_method", "show_progress")),
  PROCRUSTES = list(
    callee = "efa_procrustes",
    formals = c("A", "Target", "rotation", "S", "T_init", "oblique_eps",
                "oblique_maxit", "oblique_max_line_search", "oblique_step0",
                "oblique_normalize", "oblique_random_starts",
                "oblique_screen_keep", "oblique_triage_maxit",
                "oblique_triage_improve_tol")),
  EFA_POOLED = list(
    callee = "efa_mi",
    translating = TRUE,
    formals = c("data_list", "p", "target_method", "align_unrotated",
                "fit_pool_method", "consensus_args", "procrustes_args",
                "rmsea_ci_level", "rmsr_upper", "..."))
)

test_that("every superseded wrapper keeps a frozen signature and forwards it whole", {
  for (old_name in names(superseded_contract)) {
    spec <- superseded_contract[[old_name]]
    wrapper <- get(old_name)

    expect_identical(names(formals(wrapper)), spec$formals)

    # a translating wrapper repacks the flat tuning knobs out of its `...` first, so it is not
    # a one-line forwarder; the frozen signature above is the static half of its contract, and
    # the knob-forwarding tests below are the behavioural half.
    if (isTRUE(spec$translating)) next

    # the body is `{ efa_*(<forwarded arguments>) }`: exactly one call, to the
    # successor and to nothing else
    b <- body(wrapper)
    expect_identical(as.character(b[[1]]), "{")
    expect_length(b, 2L)

    fwd <- b[[2]]
    expect_identical(as.character(fwd[[1]]), spec$callee)

    # the first formal (`x` for every name but PROCRUSTES, whose data argument is
    # `A`) passed positionally first, `...` positionally last, and every other
    # formal passed by name to the parameter of the same name
    args <- as.list(fwd)[-1]
    arg_names <- names(args)
    if (is.null(arg_names)) arg_names <- rep("", length(args))
    passed <- vapply(seq_along(args), function(i) {
      value <- paste(deparse(args[[i]]), collapse = "")
      if (nzchar(arg_names[i])) paste0(arg_names[i], " = ", value) else value
    }, character(1))

    positional <- spec$formals[[1]]
    forwarded <- setdiff(spec$formals, c(positional, "..."))
    # a renamed argument forwards `<successor param> = <frozen formal>`; every
    # other formal forwards unchanged as `<formal> = <formal>`.
    callee_param <- forwarded
    if (!is.null(spec$renamed)) {
      hit <- forwarded %in% names(spec$renamed)
      callee_param[hit] <- spec$renamed[forwarded[hit]]
    }
    expected <- c(positional, paste0(callee_param, " = ", forwarded))
    if ("..." %in% spec$formals) expected <- c(expected, "...")
    expect_identical(passed, expected)
  }
})

# The `...`-forwarding wrappers translate rather than forward verbatim: the old flat tuning
# knobs their dots may carry are repacked into the control objects the fit now takes. Forwarded
# bare, `type` (and every other moved knob) would match no `efa_fit()` formal and be dropped,
# and the fit would quietly run the default preset instead of the requested one.

test_that(".repack_flat_dots() splits the flat knobs into the two controls", {
  out <- .repack_flat_dots(list(type = "SPSS", max_iter = 42, k = 3, maxit = 900))

  expect_s3_class(out$estimate_control, "efa_estimate_control")
  expect_s3_class(out$rotate_control, "efa_rotate_control")
  expect_identical(out$estimate_control$type, "SPSS")
  expect_identical(out$rotate_control$type, "SPSS")
  expect_equal(out$estimate_control$max_iter, 42)
  expect_equal(out$rotate_control$k, 3)
  # a genuine rotation extra is not a control knob and rides on untouched
  expect_equal(out$maxit, 900)

  # an unnamed dot rides on too: the leftovers are kept by position, not by name
  pos <- .repack_flat_dots(list(42, max_iter = 5))
  expect_identical(pos[[3]], 42)
  expect_identical(names(pos), c("estimate_control", "rotate_control", ""))

  # the former argument spellings land on the current names
  ren <- .repack_flat_dots(list(P_type = "norm", randomStarts = 7))
  expect_identical(ren$rotate_control$p_type, "norm")
  expect_equal(ren$rotate_control$random_starts, 7)

  # a knob named without a `type` resolves against the preset the flat interface defaulted to
  expect_identical(.repack_flat_dots(list(k = 2))$rotate_control$type, "EFAtools")

  # nothing to translate: the dots are handed on exactly as they came in
  expect_identical(.repack_flat_dots(list(maxit = 5)), list(maxit = 5))

  # an explicit control object is the current interface and is never second-guessed
  ctl <- list(estimate_control = estimate_control(type = "psych"), type = "SPSS")
  expect_identical(.repack_flat_dots(ctl), ctl)
})

test_that("EFA_POOLED() still tunes the fit through its dots", {
  imps <- list(test_models$baseline$cormat, test_models$baseline$cormat)
  args <- list(n_factors = 3, N = 500, method = "PAF", rotation = "promax")

  old <- do.call(EFA_POOLED, c(list(imps), args, list(type = "psych")))
  new <- do.call(efa_mi, c(list(imps), args,
                           list(estimate_control = estimate_control(type = "psych"),
                                rotate_control = rotate_control(type = "psych"))))
  expect_equal(old$rot_loadings, new$rot_loadings)

  # the knob is honoured, not merely ignored in both calls: psych is not the default preset
  default <- do.call(efa_mi, c(list(imps), args))
  expect_false(isTRUE(all.equal(old$rot_loadings, default$rot_loadings)))
})

test_that("PARALLEL() still tunes the EFA eigenvalues through its dots", {
  cormat <- test_models$baseline$cormat
  args <- list(N = 500, eigen_type = "EFA", n_datasets = 2)

  set.seed(42L)
  old <- suppressWarnings(do.call(PARALLEL, c(list(cormat), args, list(type = "psych"))))
  set.seed(42L)
  new <- suppressWarnings(do.call(
    efa_parallel, c(list(cormat), args,
                    list(estimate_control = estimate_control(type = "psych"),
                         rotate_control = rotate_control(type = "psych")))))
  expect_equal(old$eigenvalues, new$eigenvalues)
})

test_that("SL() still routes a non-default `type` into the second-order fit", {
  efa_mod <- oblique_efa()

  # `efa_schmid_leiman()` has no `type` formal any more, so the frozen argument only reaches the
  # second-order fit if the wrapper repacks it into the estimation control. A default-valued
  # `type` would pass either way, so the preset pinned here has to be a non-default one.
  old <- SL(efa_mod, type = "SPSS", method = "ULS")
  new <- efa_schmid_leiman(efa_mod, method = "ULS",
                           estimate_control = estimate_control(type = "SPSS"))
  expect_identical(old, new)

  # and it is honoured, not merely accepted: SPSS is not the preset the fit defaults to
  default <- efa_schmid_leiman(efa_mod, method = "ULS")
  expect_false(isTRUE(all.equal(old$settings, default$settings)))
})

test_that("SL() can still supply the PAF knobs that type = 'none' requires", {
  efa_mod <- oblique_efa()

  # forwarded bare these knobs were dropped, leaving type = "none" with nothing to resolve
  expect_no_error(
    sl <- SL(efa_mod, type = "none", method = "PAF", init_comm = "smc",
             criterion = 1e-3, criterion_type = "sum", max_iter = 300, abs_eigen = TRUE)
  )
  expect_s3_class(sl, "efa_schmid_leiman")
})

test_that("HULL() still tunes the fit through its dots", {
  cormat <- test_models$baseline$cormat
  args <- list(N = 500, method = "PAF", gof = "CAF")

  # HULL() fits every 1..J model through the estimation engine, so a flat knob passed
  # through its dots has to reach that fit. It is repacked into the control objects, and
  # efa_hull() must be able to take the controls it is handed.
  old <- do.call(HULL, c(list(cormat), args, list(type = "psych")))
  new <- do.call(efa_hull, c(list(cormat), args,
                             list(estimate_control = estimate_control(type = "psych"),
                                  rotate_control = rotate_control(type = "psych"))))
  expect_identical(old, new)

  # the knob is honoured, not merely accepted: max_iter = 1 halts PAF after a single
  # iteration, so the goodness-of-fit values the hull is built from must change
  pinned <- suppressWarnings(
    do.call(efa_hull, c(list(cormat), args,
                        list(estimate_control = estimate_control(max_iter = 1)))))
  default <- do.call(efa_hull, c(list(cormat), args))
  expect_false(isTRUE(all.equal(pinned$results, default$results)))
})

# Silence: a superseded name must add no condition of its own. The criteria below take a
# correlation matrix, so nothing is emitted about deriving one from raw data.

test_that("the retention criteria wrappers add no condition", {
  cormat <- test_models$baseline$cormat

  expect_no_condition(PARALLEL(cormat, N = 500, n_datasets = 2))
  expect_no_condition(KGC(cormat))
  expect_no_condition(SCREE(cormat))
  expect_no_condition(SMT(cormat, N = 500))
  expect_no_condition(NEST(cormat, N = 500, n_datasets = 2))
  expect_no_condition(HULL(cormat, N = 500, method = "ML", gof = "CAF"))
})

test_that("CD() adds no condition and stays transparent", {
  set.seed(42L)
  old <- CD(GRiPS_raw, N_pop = 500, N_samples = 50)
  set.seed(42L)
  new <- efa_cd(GRiPS_raw, N_pop = 500, N_samples = 50)

  expect_identical(old, new)
  expect_identical(class(old), "efa_retention")

  set.seed(42L)
  expect_no_condition(CD(GRiPS_raw, N_pop = 500, N_samples = 50))
})

# A flat tuning knob no longer has an efa_fit() formal, so a bare copy would land in the
# rotation extras and be dropped -- the fit would silently run the default preset. It is
# rejected instead, and the superseded names (which repack it into the controls) keep working.

test_that("efa_fit() rejects a flat tuning knob passed as a bare dot", {
  cormat <- test_models$baseline$cormat

  for (knob in list(list(type = "SPSS"), list(max_iter = 500), list(init_comm = "mac"),
                    list(criterion = 1e-4), list(criterion_type = "sum"),
                    list(abs_eigen = TRUE), list(start_method = "factanal"),
                    list(normalize = FALSE), list(precision = 1e-6),
                    list(order_type = "eigen"), list(varimax_type = "kaiser"),
                    list(p_type = "norm"), list(k = 3), list(random_starts = 2),
                    list(P_type = "norm"), list(randomStarts = 2))) {
    expect_error(
      do.call(efa_fit, c(list(cormat, n_factors = 3, N = 500, method = "PAF",
                              rotation = "promax"), knob)),
      class = "efa_flat_knob_in_dots")
  }

  # the functions that forward their dots to a fit surface the same classed error, rather
  # than an opaque "criterion could not be run"
  expect_error(efa_retain(cormat, N = 500, criteria = "HULL", type = "psych"),
               class = "efa_flat_knob_in_dots")
  expect_error(efa_hull(cormat, N = 500, method = "PAF", gof = "CAF", max_iter = 500),
               class = "efa_flat_knob_in_dots")
  expect_error(efa_kgc(cormat, eigen_type = "EFA", init_comm = "mac"),
               class = "efa_flat_knob_in_dots")
  expect_error(efa_mi(list(cormat, cormat), n_factors = 2, type = "psych"),
               class = "efa_flat_knob_in_dots")

  # arguments that really are efa_fit() formals still pass through the same dots
  expect_no_error(suppressWarnings(
    efa_retain(cormat, N = 500, criteria = "HULL", method = "PAF")))

  # a genuine rotation extra is still forwarded, not rejected
  expect_no_error(efa_fit(cormat, n_factors = 3, N = 500, method = "PAF",
                          rotation = "oblimin", maxit = 750))
  # and the control objects remain the supported route (pinning a knob alongside a `type`
  # raises the usual preset-override notice, which is not what this test is about)
  expect_no_error(suppressWarnings(
    efa_fit(cormat, n_factors = 3, N = 500, method = "PAF", rotation = "promax",
            estimate_control = estimate_control(type = "SPSS", max_iter = 500))))
})

test_that("the superseded names still tune the fit despite the flat-knob guard", {
  cormat <- test_models$baseline$cormat

  # EFA() repacks its frozen flat arguments into the controls, so the guard never sees them
  expect_no_error(old <- suppressWarnings(
    EFA(cormat, n_factors = 3, N = 500, method = "PAF", rotation = "promax",
        type = "SPSS", max_iter = 500)))
  new <- suppressWarnings(
    efa_fit(cormat, n_factors = 3, N = 500, method = "PAF", rotation = "promax",
            estimate_control = estimate_control(type = "SPSS", max_iter = 500),
            rotate_control = rotate_control(type = "SPSS")))
  expect_identical(old$rot_loadings, new$rot_loadings)

  # the same holds for the wrappers whose dots reach a fit
  expect_no_error(HULL(cormat, N = 500, method = "PAF", gof = "CAF", max_iter = 500))
  expect_no_error(N_FACTORS(cormat, N = 500, criteria = "EKC", type = "psych"))
})

test_that("N_FACTORS() keeps partial-matching `max_iter` to its frozen `max_iter_CD`", {
  cormat <- test_models$baseline$cormat

  # `max_iter` is a unique prefix of the frozen `max_iter_CD` formal, so R has always matched it
  # there rather than passing it on through the dots. efa_retain() now rejects that spelling --
  # the tuning knob it looks like belongs in estimate_control() -- but the old name is frozen,
  # so the argument keeps the meaning every existing script gave it.
  nf <- suppressWarnings(suppressMessages(
    N_FACTORS(cormat, N = 500, suitability = FALSE, criteria = "EKC", max_iter = 5)))
  expect_identical(nf$settings$max_iter_CD, 5)

  expect_error(efa_retain(cormat, N = 500, criteria = "EKC", max_iter = 5),
               class = "efa_flat_knob_in_dots")
})
