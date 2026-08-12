# Declarative control objects for the estimation and rotation tuning knobs.
#
# `estimate_control()` and `rotate_control()` bundle the tuning arguments that are
# otherwise passed one by one to a fit. Each returns a validated, classed list that
# is a 1:1 declarative surface over the settings resolved internally by
# `.estimate_model()` (via `.PAF()`) and `.rotate_model()`. The constructors only
# validate and store: they carry the chosen `type` preset and the declared knobs
# (with `NA` marking "unset, resolve from the preset") and leave the actual preset
# resolution -- and any override warning -- to the fitting cores, because which
# preset applies depends on the estimator and rotation, which a control object does
# not know.

# TRUE when a knob is left unset: a single genuine NA, the "resolve from the type
# preset" sentinel. NaN is deliberately not treated as unset (is.na(NaN) is TRUE) so a
# stray NaN falls through to the type check and is rejected rather than silently stored.
.is_control_unset <- function(x) length(x) == 1L && is.na(x) && !is.nan(x)

# Index of a single knob value among its choices, case-folded and allowing unambiguous
# abbreviations. NA marks "no usable match": a non-string input, no match at all, or -- from
# charmatch()'s 0 -- an abbreviation that prefixes several choices.
.control_choice_index <- function(x, choices) {
  if (!checkmate::test_string(x)) return(NA_integer_)
  idx <- charmatch(.fold_upper(x), .fold_upper(choices))
  if (!is.na(idx) && idx == 0L) NA_integer_ else idx
}

# Resolve a choice-valued knob (or NA, when `na_ok`) against its choices, case-insensitively
# and allowing unambiguous abbreviations as everywhere else on the interface, and return the
# canonical spelling so the stored knob and the printed control always carry it. Classed abort
# otherwise.
.assert_control_choice <- function(x, choices, arg, na_ok = TRUE) {
  if (na_ok && .is_control_unset(x)) return(x)
  idx <- .control_choice_index(x, choices)
  if (is.na(idx)) {
    na_txt <- if (na_ok) " (or NA to resolve from the type preset)" else ""
    cli::cli_abort(
      "{.arg {arg}} must be one of {.val {choices}}{na_txt}.",
      class = "efa_control_input"
    )
  }
  choices[idx]
}

# Validate a numeric knob against its admissible range (optionally a whole number, optionally
# NA); classed abort otherwise. `positive` demands x > 0; `lower` bounds it inclusively and
# `upper` inclusively unless `upper_strict`, which demands x < upper. The bounds mirror the ones
# the fitting core enforces, so a control object can never carry a value the fit would go on to
# reject.
.assert_control_number <- function(x, arg, positive = TRUE, lower = -Inf, upper = Inf,
                                   na_ok = TRUE, int = FALSE, upper_strict = FALSE) {
  if (na_ok && .is_control_unset(x)) return(invisible(NULL))
  ok <- if (int) checkmate::test_int(x) else checkmate::test_number(x, finite = TRUE)
  if (isTRUE(ok) && positive) ok <- x > 0
  if (isTRUE(ok)) ok <- x >= lower && if (upper_strict) x < upper else x <= upper
  if (!isTRUE(ok)) {
    kind <- if (int) "whole number" else "number"
    bounds <- character()
    if (positive) bounds <- c(bounds, "greater than 0")
    if (!positive && is.finite(lower)) bounds <- c(bounds, paste("of at least", lower))
    if (is.finite(upper)) {
      bounds <- c(bounds, paste(if (upper_strict) "smaller than" else "of at most", upper))
    }
    bound <- if (length(bounds)) paste0(" ", paste(bounds, collapse = " and ")) else ""
    na_txt <- if (na_ok) " (or NA to resolve from the type preset)" else ""
    cli::cli_abort(
      "{.arg {arg}} must be a single {kind}{bound}{na_txt}.",
      class = "efa_control_input"
    )
  }
  invisible(NULL)
}

# Resolve a choice-valued knob that the flat interface partial-matched, keeping that
# convenience: an abbreviation resolves to the full value and NA stays unset (it means "not
# needed by this estimator"). Returns the resolved value; classed abort otherwise.
.match_control_choice <- function(x, choices, arg) {
  if (.is_control_unset(x)) return(x)
  idx <- .control_choice_index(x, choices)
  if (is.na(idx)) {
    cli::cli_abort(
      "{.arg {arg}} must be one of {.val {choices}} (or NA to leave it unset).",
      class = "efa_control_input"
    )
  }
  choices[idx]
}

# ASCII-only upper-casing for case-folded matching. toupper() is locale-dependent -- under a
# Turkish/Azeri locale toupper("i") is the dotted capital "İ", so an input like "minres"
# or "oblimin" would fold to a string the (already-uppercase) choice can never equal. The
# choice values are all ASCII, so folding both sides with a fixed a-z -> A-Z map keeps the
# matching identical on every platform.
.fold_upper <- function(x) chartr("abcdefghijklmnopqrstuvwxyz",
                                  "ABCDEFGHIJKLMNOPQRSTUVWXYZ", x)

# Case-insensitive match.arg(): resolve a choice-valued argument against its choices and
# return the canonical spelling(s) from `choices`, so downstream code, stored settings, and
# printed output always carry the canonical value. Values are compared case-folded, and an
# exact (case-folded) match beats a prefix match, so unambiguous abbreviations still resolve
# as they do under match.arg(). It is deliberately stricter than match.arg() in `several.ok`
# mode: an unmatched element aborts here, whereas match.arg() silently drops it. Mirroring
# match.arg(), an `arg` identical to `choices` selects the default -- the first choice, or
# all of them when `several.ok` -- and a missing `choices` is taken from the caller's formal
# default; a NULL `arg` selects the documented default too, see the comment on that branch.
# Ambiguous or unmatched values are a classed abort listing the choices; `class` lets a
# caller raise its own condition class (the control constructors reuse their
# `efa_control_input`).
.match_arg_ci <- function(arg, choices, several.ok = FALSE,
                          arg_name = deparse1(substitute(arg)),
                          class = "efa_bad_choice") {
  caller <- sys.parent()
  from_formals <- missing(choices)
  if (from_formals) {
    choices <- eval(formals(sys.function(caller))[[arg_name]],
                    envir = sys.frame(caller))
  }
  # NULL means "use the default", as in match.arg(), but the default is the caller's formal
  # default, NOT choices[1L]: the two part company under `several.ok` whenever the formal
  # default is a whole vector (efa_hull()'s `gof`, efa_kgc()'s `eigen_type`), and wherever an
  # explicit `choices` is a superset of the formal default (efa_retain()'s `criteria`) or
  # orders it differently (`eigen_type_other`). Taking choices[1L] there silently runs
  # something narrower than the help page documents, with no condition raised. A formal
  # without a default, or one that is itself NULL (efa_simulate()'s `match`, documented as
  # resolving to "thresholds"), falls back to choices[1L].
  if (is.null(arg)) {
    if (from_formals) {
      # `choices` already is that default -- no need to look it up a second time
      arg <- choices
    } else {
      # Keep the default inside a length-one list while testing it: a formal without a
      # default is the empty symbol, and binding that to a variable (or passing it on)
      # forces it into R's "argument is missing" error. Comparing the enclosing list never
      # evaluates it.
      default <- formals(sys.function(caller))[arg_name]
      arg <- if (identical(unname(as.list(default)), list(quote(expr = )))) {
        NULL
      } else {
        eval(default[[1L]], envir = sys.frame(caller))
      }
    }
    if (is.null(arg)) return(choices[1L])
  }
  if (identical(arg, choices)) return(if (several.ok) choices else choices[1L])

  # charmatch(): NA = no match, 0 = ambiguous (a prefix of several choices)
  idx <- if (is.character(arg)) charmatch(.fold_upper(arg), .fold_upper(choices)) else NA_integer_
  len_ok <- if (several.ok) length(arg) >= 1L else length(arg) == 1L
  if (!len_ok || anyNA(idx) || any(idx == 0L)) {
    no_match <- if (is.character(arg)) unique(arg[is.na(idx)]) else character()
    ambiguous <- if (is.character(arg)) unique(arg[!is.na(idx) & idx == 0L]) else character()
    what <- if (several.ok) "one or more" else "one"
    cli::cli_abort(
      c("{.arg {arg_name}} must be {what} of {.val {choices}}.",
        if (!is.character(arg)) c("x" = "You supplied {.obj_type_friendly {arg}}."),
        if (!several.ok && is.character(arg) && length(arg) > 1L) {
          c("x" = "You supplied {length(arg)} values.")
        },
        if (length(no_match) > 0L) {
          c("x" = "{.val {no_match}} {?does/do} not match any of the choices.")
        },
        if (length(ambiguous) > 0L) {
          c("x" = "{.val {ambiguous}} {?matches/match} more than one choice.")
        }),
      class = class
    )
  }
  choices[idx]
}

# Case-insensitively map the elements of a vector-valued choice argument onto the canonical
# spelling of `choices` (exact case-folded matches only -- no abbreviations, matching the
# subset semantics of the callers). An element matching no choice is returned unchanged, so
# the caller's subset assertion still rejects it with its usual message.
.map_subset_ci <- function(x, choices) {
  if (!is.character(x)) return(x)
  idx <- match(.fold_upper(x), .fold_upper(choices))
  found <- !is.na(idx)
  x[found] <- choices[idx[found]]
  x
}

# Validate a logical knob (or NA, when `na_ok`); classed abort otherwise.
.assert_control_flag <- function(x, arg, na_ok = TRUE) {
  if (na_ok && .is_control_unset(x)) return(invisible(NULL))
  if (!checkmate::test_flag(x)) {
    na_txt <- if (na_ok) " (or NA to resolve from the type preset)" else ""
    cli::cli_abort(
      "{.arg {arg}} must be a single logical value{na_txt}.",
      class = "efa_control_input"
    )
  }
  invisible(NULL)
}

#' Control objects for estimation and rotation settings
#'
#' `estimate_control()` and `rotate_control()` collect the estimation and rotation
#' *tuning* arguments of a factor analysis into two small, validated objects. They
#' are a declarative surface over the same settings resolved internally by the
#' package's estimation and rotation engines, so that a fit's many tuning knobs can
#' be prepared, inspected, and reused as a single value instead of being passed one
#' by one.
#'
#' Each argument that governs a `type` preset defaults to `NA`, meaning "leave this
#' knob to the preset". Setting `type` to one of `"EFAtools"`, `"psych"`, or
#' `"SPSS"` fills those knobs from the corresponding preset when the fit is run;
#' setting `type = "none"` requires the relevant knobs to be supplied explicitly.
#' The control object only records the chosen `type` and the knobs you set: the
#' preset is resolved (and any "argument set alongside `type`" warning issued) when
#' the object is used to fit a model, exactly as it is today, because which preset
#' applies depends on the estimator and rotation.
#'
#' @param type character. One of `"EFAtools"` (default), `"psych"`, `"SPSS"`, or
#'   `"none"`. Selects the preset that fills the `NA`-defaulted knobs below when the
#'   control is used to fit a model.
#' @param init_comm character. Method for the initial communalities in principal
#'   axis factoring: `"smc"` (squared multiple correlations), `"mac"` (maximum
#'   absolute correlations), or `"unity"`. `NA` (default) resolves from `type`.
#' @param criterion numeric. The convergence criterion for principal axis
#'   factoring: iteration stops once the change in communalities falls below it. A
#'   single number greater than 0 and smaller than 1; `NA` (default) resolves from `type`.
#' @param criterion_type character. The convergence criterion type for principal
#'   axis factoring: `"max_individual"` (the largest change in any communality, as
#'   in SPSS) or `"sum"` (the change in the summed communalities, as in
#'   [psych::fa()]). `NA` (default) resolves from `type`.
#' @param max_iter numeric. The maximum number of principal-axis-factoring
#'   iterations before the procedure is halted with a warning. A single whole
#'   number of at least 1; `NA` (default) resolves from `type`.
#' @param abs_eigen logical. Which algorithm the principal-axis-factoring
#'   iterations use: `FALSE` computes the loadings from the eigenvalues (as in
#'   [psych::fa()]); `TRUE` uses the absolute eigenvalues (as in SPSS). `NA`
#'   (default) resolves from `type`.
#' @param start_method character. Starting values for the maximum-likelihood
#'   optimiser: `"psych"` (default, the [psych::fa()] starts) or `"factanal"` (the
#'   [stats::factanal()] starts); abbreviations are matched. Not governed by `type`.
#'   Only maximum likelihood uses it, so `NA` leaves it unset and is rejected only by
#'   a fit that is actually run with `estimator = "ML"`.
#' @param fiml_max_iter numeric. The maximum number of EM iterations used to estimate the
#'   two-stage full-information maximum-likelihood moments from raw data with missing
#'   values (`cor_method = "fiml"`); the last iterate is returned, with a warning, if the
#'   cap is reached. A single whole number of at least 1; default `500`. Not governed by
#'   `type`, and unused by every other correlation method. The EM converges linearly and
#'   needs more iterations the larger the fraction of missing information, so raise it when
#'   a fit reports that the moments did not converge.
#' @param fiml_tol numeric. The convergence tolerance of that EM: iteration stops once the
#'   largest change in the standardized moments (the standardized means, log-variances, and
#'   correlations) falls below it, so it does not depend on the variables' measurement
#'   scale. A single number greater than 0 and smaller than 1 (at or above 1 the criterion is
#'   met immediately and the starting moments would be returned as converged); default `1e-5`.
#'   Not governed by `type`.
#'
#' @returns `estimate_control()` returns a list of class `efa_estimate_control` with
#'   the components `type`, `init_comm`, `criterion`, `criterion_type`, `max_iter`,
#'   `abs_eigen`, `start_method`, `fiml_max_iter`, and `fiml_tol`.
#'   `rotate_control()` returns a list of class
#'   `efa_rotate_control` with the components `type`, `normalize`, `precision`,
#'   `order_type`, `varimax_type`, `p_type`, `k`, `random_starts`, and `extra_args`
#'   (a named list of any additional arguments forwarded to the rotation engine).
#'
#' @seealso [efa_fit()], which takes both controls; [efa_retain()], the retention
#'   criteria, and [efa_schmid_leiman()], which take an `estimate_control` for the
#'   fits they run.
#'
#' @family Control functions
#'
#' @examples
#' # Estimation knobs taken entirely from a preset:
#' estimate_control(type = "SPSS")
#'
#' # A preset with one knob pinned to a non-preset value:
#' estimate_control(type = "EFAtools", max_iter = 500)
#'
#' # Every knob supplied explicitly (type = "none"):
#' estimate_control(type = "none", init_comm = "smc", criterion = 1e-3,
#'                  criterion_type = "sum", max_iter = 300, abs_eigen = TRUE)
#'
#' @export
estimate_control <- function(type = c("EFAtools", "psych", "SPSS", "none"),
                             init_comm = NA, criterion = NA, criterion_type = NA,
                             max_iter = NA, abs_eigen = NA, start_method = "psych",
                             fiml_max_iter = 500, fiml_tol = 1e-5) {

  # Match `type` case-insensitively, but keep the constructor's own condition class so a bad
  # `type` is caught alongside the other knobs (all `efa_control_input`).
  type <- .match_arg_ci(type, class = "efa_control_input")

  init_comm <- .assert_control_choice(init_comm, c("smc", "mac", "unity"), "init_comm")
  # The convergence tolerance is compared against a change in communalities, so a value at or
  # above 1 stops the iterations immediately; the fit rejects it (see .PAF()), and the control
  # must not be able to carry a value the fit would go on to reject.
  .assert_control_number(criterion, "criterion", upper = 1, upper_strict = TRUE)
  criterion_type <- .assert_control_choice(criterion_type, c("max_individual", "sum"),
                                           "criterion_type")
  .assert_control_number(max_iter, "max_iter", int = TRUE)
  .assert_control_flag(abs_eigen, "abs_eigen")
  # `start_method` only governs the ML optimiser, so NA ("not needed here") is admissible and
  # is rejected by the fit itself, and only when `estimator = "ML"`.
  start_method <- .match_control_choice(start_method, c("psych", "factanal"), "start_method")
  # The two FIML knobs are not preset-driven, so they carry real defaults rather than the NA
  # sentinel and NA is not admissible. The tolerance is bounded above at 1 for the same reason
  # `criterion` is: it is compared against a change in standardized moments -- standardized means,
  # log-variances, and correlations -- so a value at or above 1 is met on the first iteration and
  # the EM would return its starting moments while reporting itself converged. Unlike `criterion`
  # the fit does not reject such a value downstream (the EM only demands tol > 0), so the bound has
  # to be enforced here or not at all.
  .assert_control_number(fiml_max_iter, "fiml_max_iter", int = TRUE, na_ok = FALSE)
  .assert_control_number(fiml_tol, "fiml_tol", upper = 1, upper_strict = TRUE, na_ok = FALSE)

  structure(
    list(type = type, init_comm = init_comm, criterion = criterion,
         criterion_type = criterion_type, max_iter = max_iter,
         abs_eigen = abs_eigen, start_method = start_method,
         fiml_max_iter = fiml_max_iter, fiml_tol = fiml_tol),
    class = "efa_estimate_control"
  )
}

#' @rdname estimate_control
#'
#' @param normalize logical. If `TRUE` (default), a Kaiser normalization is
#'   performed before the rotation. The one knob that is always on unless you turn
#'   it off with `FALSE`.
#' @param precision numeric. The convergence tolerance of the rotation procedure. A
#'   single number greater than 0 and at most 1; default `1e-5`. Each rotation stage
#'   monitors its own quantity, so the same number is not the same tolerance everywhere:
#'   `varimax_type = "kaiser"` stops on the *absolute* change in the varimax simplicity
#'   criterion, which is an average over variables (and so does not scale with how many
#'   there are) but rises with the number of factors, roughly toward
#'   `1 - 1 / n_factors`, so a fixed value is a relatively weaker tolerance the more
#'   factors are extracted; `varimax_type = "svd"` stops on the *relative* change in the
#'   singular values (as in [stats::varimax()]); and the criterion rotations fitted by
#'   gradient projection stop when the *projected-gradient norm* falls below it. Promax
#'   inherits whichever of the two varimax tests its `varimax_type` selects, because it
#'   rotates a varimax base.
#' @param order_type character. How the factors are ordered: `"eigen"` (by
#'   descending sum of squared loadings, as in [psych::fa()]) or `"ss_factors"` (by
#'   descending unweighted sum of squared loadings). `NA` (default) resolves from
#'   `type`.
#' @param varimax_type character. The varimax variant used (for the varimax and
#'   promax rotations): `"svd"` (as in [stats::varimax()]) or `"kaiser"` (the SPSS
#'   / Kaiser (1958) procedure). `NA` (default) resolves from `type`.
#' @param p_type character. How the promax target matrix is computed: `"unnorm"`
#'   (the unnormalized target of Hendrickson & White (1964), also used by psych and
#'   stats) or `"norm"` (the normalized target used by SPSS). `NA` (default)
#'   resolves from `type`.
#' @param k numeric. The promax power (for the target matrix) or the number of
#'   near-zero loadings for simplimax. A single number greater than 0; `NA`
#'   (default) leaves it to the fit (the `type`-dependent promax value, or
#'   `nrow(loadings)` for simplimax). Simplimax counts loadings, so a fit using it
#'   additionally requires a whole number no larger than the number of loadings in
#'   the solution; promax's power has no such restriction.
#' @param random_starts numeric. The number of random starts used by the
#'   criterion-based rotations to guard against local minima. A single whole number
#'   of at least 0, where `0` runs the rotation from its warm start only; default
#'   `100`. The default suffices for the smooth criteria; `simplimax` remains
#'   materially start-dependent at it, so raise it there (see the *Rotations* section
#'   of [efa_fit()]).
#' @param ... Additional arguments forwarded to the rotation engine. Only the names
#'   a rotation engine can consume are accepted: `maxit` (a whole number of at least 0
#'   bounding a *single* gradient-projection optimization -- the multi-start search runs
#'   several of them and each is bounded separately, so it is not a budget for the run as a
#'   whole; varimax and promax have no such stage and take only `precision`), and the
#'   criterion parameters `gam` (oblimin; `gam = 0` is the
#'   recommended default, and larger values increasingly reward correlated factors and
#'   can drive the solution toward factor collapse, so inspect `Phi` before interpreting
#'   a fit with `gam > 0`) and `delta` (geomin; a positive number, default `0.01`);
#'   anything else is rejected as a
#'   misspelling. They are stored in `extra_args` and passed on to the rotation engine
#'   when the control is used to fit a model; an extra a given fit's rotation does not
#'   consume is ignored by that fit, so one control can serve fits with different
#'   rotations. An estimation knob (which belongs in [estimate_control()]) or one of the
#'   former spellings `P_type` and `randomStarts` is likewise rejected here, because the
#'   fit would silently drop it.
#'
#' @examples
#'
#' # Rotation knobs taken from a preset:
#' rotate_control(type = "psych")
#'
#' # A criterion-specific extra argument, forwarded to the rotation engine:
#' rotate_control(type = "EFAtools", k = 3, gam = 0.5)
#'
#' @export
rotate_control <- function(type = c("EFAtools", "psych", "SPSS", "none"),
                           normalize = TRUE, precision = 1e-5, order_type = NA,
                           varimax_type = NA, p_type = NA, k = NA,
                           random_starts = 100, ...) {

  type <- .match_arg_ci(type, class = "efa_control_input")

  .assert_control_flag(normalize, "normalize", na_ok = FALSE)
  .assert_control_number(precision, "precision", upper = 1, na_ok = FALSE)
  order_type <- .assert_control_choice(order_type, c("eigen", "ss_factors"), "order_type")
  varimax_type <- .assert_control_choice(varimax_type, c("svd", "kaiser"), "varimax_type")
  p_type <- .assert_control_choice(p_type, c("norm", "unnorm"), "p_type")
  .assert_control_number(k, "k")
  # 0 is a meaningful setting, not a missing one: it runs the rotation from its warm start
  # only, with no random restarts.
  .assert_control_number(random_starts, "random_starts", positive = FALSE, lower = 0,
                         int = TRUE, na_ok = FALSE)

  extra_args <- list(...)
  # A misplaced tuning knob is not a rotation engine extra: its name collides with a formal
  # of the fit, so the extras splice in efa_fit() would drop it and the fit would quietly run
  # the preset value instead. Refuse to carry it rather than let it be silently ignored.
  bad_est <- intersect(names(extra_args), .flat_estimate_knobs)
  bad_flat <- intersect(names(extra_args), c("P_type", "randomStarts"))
  if (length(bad_est) > 0L || length(bad_flat) > 0L) {
    cli::cli_abort(
      c("{.arg {c(bad_est, bad_flat)}} cannot be passed to {.fn rotate_control} as a rotation
         engine extra.",
        if (length(bad_est) > 0L) c("i" = "The estimation knobs live in {.fn estimate_control}."),
        if (length(bad_flat) > 0L) c("i" = "The former {.arg P_type} and {.arg randomStarts} are
                                            the formals {.arg p_type} and {.arg random_starts}.")),
      class = "efa_control_input"
    )
  }
  # Any other extra must be a name some rotation's engine can consume (the engines read their
  # extras by exact name); an unknown one is a misspelling that the fit would silently drop.
  # Deliberately validated against the across-rotation union, not a single rotation: a control
  # is reusable across fits, so it may legitimately carry an extra that some of its fits (e.g.
  # the unrotated or promax rows of an efa_average() grid) never consume.
  extra_nms <- names(extra_args)
  if (is.null(extra_nms)) extra_nms <- rep("", length(extra_args))
  bad_extra <- unique(setdiff(extra_nms[nzchar(extra_nms)], .rotation_extra_union))
  if (length(bad_extra) > 0L || any(!nzchar(extra_nms))) {
    msg <- if (length(bad_extra) > 0L) {
      "{.arg {bad_extra}} {?is/are} not {?a name/names} a rotation engine can consume."
    } else {
      "The extra engine arguments in {.arg ...} must be named."
    }
    cli::cli_abort(
      c(msg,
        "i" = "The accepted engine extras are {.arg {(.rotation_extra_union)}}.",
        "i" = "Check for a misspelled engine argument (for example {.arg gam}, not
               {.arg gamma}); the rotation tuning knobs are {.fn rotate_control}'s named
               arguments."),
      class = "efa_control_input"
    )
  }

  structure(
    list(type = type, normalize = normalize, precision = precision,
         order_type = order_type, varimax_type = varimax_type, p_type = p_type,
         k = k, random_starts = random_starts, extra_args = extra_args),
    class = "efa_rotate_control"
  )
}

# The flat estimation and rotation tuning knobs the pre-control interface exposed, split by the
# control object each one now lives in.
.flat_estimate_knobs <- c("init_comm", "criterion", "criterion_type", "max_iter",
                          "abs_eigen", "start_method")
.flat_rotate_knobs <- c("normalize", "precision", "order_type", "varimax_type", "p_type",
                        "k", "random_starts")

# The efa_fit() arguments that govern its inference machinery (standard errors and the
# random-number stream feeding them). They are efa_fit() formals, so a caller forwarding
# `...` into a fit accepts them by construction -- which is wrong wherever the fit is only
# an internal step whose standard errors are discarded (see .reject_inference_dots()).
.fit_inference_args <- c("se", "b_boot", "ci", "seed")

# Reject a flat tuning knob that reached a fit as a bare dot. No fitting function has a formal
# for any of these any more -- they live in the control objects -- so a bare copy would be taken
# for a rotation extra and dropped, and the fit would quietly run the default preset instead
# of the requested settings. Takes the argument NAMES (`...names()`), so the dots are never
# forced; `fn` names the function the caller actually reached, so the message points at the call
# that has to change.
.reject_flat_knobs <- function(nms, fn = "efa_fit") {
  # The renamed `method` argument is checked first: it is the more specific mistake, and
  # every function that must reject the flat knobs must reject the former estimator
  # spelling too, so one call site covers both.
  .reject_renamed_method(nms, fn = fn)
  bad <- intersect(nms, c("type", .flat_estimate_knobs, .flat_rotate_knobs,
                          "P_type", "randomStarts"))
  if (length(bad) == 0L) return(invisible(NULL))
  # point the example at the constructor the offending knob actually lives in
  example <- if (all(bad %in% c(.flat_rotate_knobs, "P_type", "randomStarts"))) {
    paste0(fn, "(x, ..., rotate_control = rotate_control(k = 3))")
  } else {
    paste0(fn, "(x, ..., estimate_control = estimate_control(max_iter = 500))")
  }
  cli::cli_abort(
    c("{.arg {bad}} cannot be passed to {.fn {fn}} directly.",
      "i" = "The estimation and rotation tuning knobs live in {.fn estimate_control} and
             {.fn rotate_control}.",
      "i" = "For example: {.code {example}}."),
    class = "efa_flat_knob_in_dots"
  )
}

# Reject a `...` name that no fit could ever consume. The retention criteria and their
# orchestrator forward their `...` into efa_fit(), but only when a fit actually runs (e.g.
# eigen_type includes "EFA", or a criterion that fits a model is selected). A misspelled name
# is otherwise silently dropped and the criterion runs with the default it looks like it was
# told to change. A name that is neither an efa_fit() formal nor a rotation-engine extra is
# such a misspelling and is refused here, independently of whether a fit runs. Takes the
# argument NAMES (`...names()`) so the dots are never forced; `fn` names the function reached.
#
# It must stay lenient for the superseded wrappers, which repack their flat tuning knobs into
# `estimate_control` / `rotate_control` objects and splice those into these same dots: both
# control names are efa_fit() formals, so they are in the accepted set by construction.
.reject_unknown_fit_dots <- function(nms, fn, unrotated = FALSE) {
  nms <- nms[nzchar(nms)]
  if (length(nms) == 0L) return(invisible(NULL))
  known <- c(setdiff(names(formals(efa_fit)), "..."), .rotation_extra_union)
  # A caller whose fits are always unrotated runs no rotation engine, so an engine extra
  # can never be consumed there. Refusing it here, rather than letting it reach the fit,
  # is what makes the error name the function the user called and arrive before the
  # criterion does its work (`rotation` itself is refused by .reject_rotation_dots()).
  # efa_fit()'s inference arguments are refused for the same reason, but with a message of
  # their own: they ARE efa_fit() arguments, so "check for a misspelled name" would be wrong.
  if (unrotated) {
    known <- setdiff(known, .rotation_extra_union)
    .reject_inference_dots(nms, fn = fn)
  }
  # setdiff() already returns unique values, so the offenders are named once each
  bad <- setdiff(nms, known)
  if (length(bad) == 0L) return(invisible(NULL))
  cli::cli_abort(
    c("{.arg {bad}} {?is/are} not {?an argument/arguments} of {.fn {fn}} or {.fn efa_fit}.",
      "i" = "Check for a misspelled argument name."),
    class = "efa_unused_dots"
  )
}

# Reject efa_fit()'s inference arguments in the `...` of a caller whose fits are internal
# steps. Two distinct hazards, both silent before this guard:
#   * `se`, `b_boot`, `ci`: the fit does compute the standard errors -- they are efa_fit()
#     formals, so they pass any whitelist built from them -- but the caller keeps only a few
#     quantities from the fit and never reports its standard errors, so the result is
#     identical and the work is thrown away (measurably: an information-based SE roughly
#     doubles an efa_kgc() run).
#   * `seed`: efa_fit() seeds itself with .set_local_seed(), which restores the caller's
#     stream when the fit returns, so the seed is spent on a fit that draws no random
#     numbers. Where the caller draws at all, that draw still comes from the caller's stream.
#     The name therefore looks like it pins the run and pins nothing at all; set.seed() is
#     what governs it.
# The hints must hold for every caller, so neither claims that this particular call fits, or
# simulates, or is a retention criterion: efa_retain() with a criterion such as MAP fits no
# model at all, efa_kgc() and efa_scree() fit one but simulate nothing, and
# efa_schmid_leiman() is not a criterion and keeps its internal fit's loadings rather than
# eigenvalues.
# Takes the argument NAMES (`...names()`), so the dots are never forced.
.reject_inference_dots <- function(nms, fn) {
  bad <- intersect(nms, .fit_inference_args)
  if (length(bad) == 0L) return(invisible(NULL))
  msg <- c("{.arg {bad}} cannot be passed to {.fn {fn}}.")
  if (any(c("se", "b_boot", "ci") %in% bad)) {
    msg <- c(msg, "i" = "Any fit run here is an internal step whose standard errors are not
                         reported, so they would be computed and then discarded.")
  }
  if ("seed" %in% bad) {
    msg <- c(msg, "i" = "The fits run here draw no random numbers, so {.arg seed} has
                         nothing to govern. Whatever else the call draws comes from the
                         caller's stream: call {.fn set.seed} before {.fn {fn}}.")
  }
  cli::cli_abort(msg, class = "efa_unused_dots")
}

# The retention criteria's fits are always unrotated, so a rotation setting passed through
# their `...` is meaningless: a user's `rotation = ...` makes a criterion fit spin up a
# rotation engine for nothing -- and a criterion-based rotation draws its random starts from
# the caller's stream, so it silently moves a seeded result -- while a `rotate_control` is
# simply ignored by the unrotated fit. Refuse
# both. The superseded N_FACTORS() wrapper still repacks a frozen `type` (or a flat rotation
# knob) into a rotate_control() object that legitimately rides through these dots, so a
# `rotate_control` that IS such a control object is exempt, as is a `NULL` (efa_fit()'s own
# "not supplied" default, which the repack passes through untouched) -- only some other value
# under that name, or any `rotation`, is refused. Takes the dots by value for that check.
.reject_rotation_dots <- function(dots, fn) {
  nms <- names(dots)
  if (is.null(nms)) return(invisible(NULL))
  bad <- character(0)
  if ("rotation" %in% nms) bad <- c(bad, "rotation")
  rc <- dots[nms == "rotate_control"]
  ok_rc <- vapply(rc, function(v) is.null(v) || inherits(v, "efa_rotate_control"), logical(1))
  if (!all(ok_rc)) bad <- c(bad, "rotate_control")
  if (length(bad) == 0L) return(invisible(NULL))
  cli::cli_abort(
    c("{.arg {bad}} cannot be passed to {.fn {fn}}.",
      "i" = "The fits run by the retention criteria are always unrotated."),
    class = "efa_unused_dots"
  )
}

# Reject the former `method` argument (the estimator is selected with `estimator`) when it
# lands in a fit function's `...`: it would otherwise be forwarded into the rotation extras
# or a criterion fit -- surfacing as an opaque downstream error or being silently ignored --
# instead of selecting the estimator.
.reject_renamed_method <- function(nms, fn = "efa_fit") {
  if (!"method" %in% nms) return(invisible(NULL))
  cli::cli_abort(
    c("{.arg method} is not an argument of {.fn {fn}}.",
      "i" = "The estimator is selected with {.arg estimator}, e.g.
             {.code {fn}(x, ..., estimator = \"ML\")}."),
    class = "efa_renamed_arg"
  )
}

# Validate an estimation control (NULL means "the efa_fit() default"). Shared by efa_fit() and by
# the functions that only pass the control on to a fit, so a bogus object is rejected even when
# the call ends up running no fit at all (e.g. efa_kgc() with eigen_type = "PCA").
.assert_estimate_control <- function(x) {
  if (is.null(x) || inherits(x, "efa_estimate_control")) return(invisible(NULL))
  cli::cli_abort(
    "{.arg estimate_control} must be a control object from {.fn estimate_control}.",
    class = "efa_control_input"
  )
}

# The legacy alias spellings .repack_flat_dots() translates rather than forwards. Kept as a
# named constant so the drop filter below and the repack agree by construction.
.flat_alias_names <- c("method", "type", "P_type", "randomStarts")

# The frozen wrappers that forward their `...` (N_FACTORS(), PARALLEL(), KGC(), SCREE(),
# HULL(), NEST(), EFA_POOLED()) keep the flat interface's silent-ignore contract: a name
# that nothing on the old surface could ever consume was dropped without comment, so it is
# dropped here too -- before .repack_flat_dots() translates the flat knobs and before the
# successors' guards (.reject_unknown_fit_dots(), .reject_rotation_dots(), efa_fit()'s
# per-rotation extras whitelist) would refuse it. The kept universe is everything the old
# surface consumed: the efa_fit() formals, the flat tuning knobs, and the alias spellings.
# Keeping the flat knobs matters when a control object rides alongside one (the repack then
# translates nothing), so such a mixed call still fails loudly in the successor's flat-knob
# guard instead of losing the knob here. The retention wrappers' criterion fits are always
# unrotated (the default), so a `rotation` setting and the rotation-engine extras
# (`maxit`/`gam`/`delta`) were never consumable there -- 0.8.0 ignored them without effect on
# the result -- and are dropped with the junk; EFA_POOLED() instead passes `unrotated =
# FALSE` plus the extras its selected rotation's engine reads, mirroring EFA()'s own
# rotation-aware filter. Keeping a name here is not the same as accepting it: efa_fit()'s
# inference arguments (`se`/`b_boot`/`ci`/`seed`) pass this filter as efa_fit() formals and
# are then refused by .reject_inference_dots() in the retention successors. That is a
# DELIBERATE exception to the silent-ignore rule above, not an instance of it: unlike the
# rotation extras, the flat EFA() did take all four and a pre-rename call could name them.
# It changes no legacy result, because on this path they were inert -- these wrappers hand
# EFA() a correlation matrix, which never reached its seeding branch at all, and a standard
# error the fit computed was thrown away by the criterion (`se` even aborted outright through
# some of them: "sandwich" for want of an acov, "information" through KGC()/SCREE() for want
# of an N). See .reject_inference_dots() for the two hazards themselves.
# A successor-only name carrying a malformed value (e.g.
# `rotate_control = "SPSS"`) is deliberately KEPT: no pre-rename code could have used that
# name, so the successor's loud classed validation is the right outcome, exactly as EFA()
# rejects the successor-only names outright. Unnamed dots pass through untouched, as on the
# reject side.
.drop_unknown_frozen_dots <- function(dots, extras = character(), unrotated = TRUE) {
  nms <- names(dots)
  if (is.null(nms)) return(dots)
  known <- c(setdiff(names(formals(efa_fit)), "..."),
             .flat_estimate_knobs, .flat_rotate_knobs, .flat_alias_names, extras)
  if (unrotated) known <- setdiff(known, "rotation")
  # charmatch(), not %in%: the flat interface reached its target through do.call(), where R
  # partial-matched an abbreviated name against the formals ahead of `...`, so an
  # abbreviation was a name the old surface could consume. Dropping it here would turn such
  # a call into a silently different fit; keeping it (an ambiguous abbreviation too, which
  # base R would also have refused) lets it travel on and fail loudly instead.
  dots[!nzchar(nms) | !is.na(charmatch(nms, known))]
}

# Repack the flat tuning knobs a call may still carry in its `...` into the control objects the
# fitting interface now takes. `efa_fit()` has no flat tuning formals left, so a bare knob
# forwarded through `...` would otherwise be dropped and the fit would quietly run the default
# preset instead of the requested one. Returns the arguments to splice into the call: the two
# controls -- built only when the caller actually named `type` or a knob -- followed by every
# remaining dot (the genuine rotation extras, e.g. `maxit`/`gam`) unchanged. One `type` governs
# both presets, exactly as the single flat `type` always did.
.repack_flat_dots <- function(dots, type = NULL) {

  # The flat interface selected the estimator with `method`; the fitting functions take
  # `estimator`. Translate the name first -- unlike the knob repacking below it must happen
  # even when a control object is passed alongside, because `method` never lives in a control.
  # Both names at once is a half-migrated call: silently letting either win would fit with
  # an estimator the caller did not (knowingly) choose, so it is rejected instead.
  if ("method" %in% names(dots)) {
    if ("estimator" %in% names(dots)) {
      cli::cli_abort(
        c("{.arg method} and {.arg estimator} cannot both be supplied.",
          "i" = "{.arg method} is the former name of {.arg estimator}; pass only one."),
        class = "efa_renamed_arg"
      )
    }
    dots$estimator <- dots$method
    dots$method <- NULL
  }

  # A control object passed explicitly is the current interface; never second-guess it.
  if (any(c("estimate_control", "rotate_control") %in% names(dots))) return(dots)

  # the former argument spellings reach the controls under their current names
  if ("P_type" %in% names(dots)) {
    dots$p_type <- dots$P_type
    dots$P_type <- NULL
  }
  if ("randomStarts" %in% names(dots)) {
    dots$random_starts <- dots$randomStarts
    dots$randomStarts <- NULL
  }

  if (is.null(type) && "type" %in% names(dots)) type <- dots$type
  dots$type <- NULL

  est_named <- intersect(names(dots), .flat_estimate_knobs)
  rot_named <- intersect(names(dots), .flat_rotate_knobs)

  # nothing to translate: hand the dots on exactly as they came in
  if (is.null(type) && length(est_named) == 0L && length(rot_named) == 0L) return(dots)

  # the preset the flat interface applied when a call tuned a knob without naming a `type`
  if (is.null(type)) type <- "EFAtools"

  # keep the leftovers by logical indexing: indexing by name would corrupt an unnamed dot
  # (its "" matches nothing, so it would come back as a NULL named NA)
  nms <- names(dots)
  if (is.null(nms)) nms <- character(length(dots))
  c(list(estimate_control = do.call(estimate_control, c(list(type = type), dots[est_named])),
         rotate_control = do.call(rotate_control, c(list(type = type), dots[rot_named]))),
    dots[!(nms %in% c(est_named, rot_named))])
}

# Render one "<label>: <value>" line inside a cli container. An unset (NA)
# preset-driven knob shows as a dim marker; other values are shown verbatim.
.control_knob_line <- function(label, value) {
  shown <- if (.is_control_unset(value)) {
    cli::col_grey("<from type preset>")
  } else {
    paste(format(value), collapse = ", ")
  }
  cli::cli_text("{.field {label}}: {shown}")
}

#' Print and format a control object
#'
#' `print()` shows the chosen `type` and each tuning knob, with an unset (`NA`)
#' preset-driven knob marked as resolved from the `type` preset. `format()`
#' assembles the same report and returns it as a character vector; `print()` is
#' `cat(format(x), sep = "\n")`. The lines follow the active console theme, so they
#' are plain when colours are disabled.
#'
#' @param x A control object from [estimate_control()] or [rotate_control()].
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines.
#'
#' @seealso [estimate_control()], [rotate_control()]
#'
#' @family Control functions
#'
#' @name print.efa_control
#'
#' @examples
#' est <- estimate_control(type = "SPSS")
#' est
#' writeLines(format(est))
#'
NULL

#' @rdname print.efa_control
#' @export
#' @method print efa_estimate_control
print.efa_estimate_control <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.efa_control
#' @export
#' @method format efa_estimate_control
format.efa_estimate_control <- function(x, ...) {
  cli::cli_format_method({
    cli::cli_text("{.strong Estimation control} ({.field type}: {.val {x$type}})")
    cli::cli_text("")
    .control_knob_line("init_comm", x$init_comm)
    .control_knob_line("criterion", x$criterion)
    .control_knob_line("criterion_type", x$criterion_type)
    .control_knob_line("max_iter", x$max_iter)
    .control_knob_line("abs_eigen", x$abs_eigen)
    .control_knob_line("start_method", x$start_method)
    .control_knob_line("fiml_max_iter", x$fiml_max_iter)
    .control_knob_line("fiml_tol", x$fiml_tol)
  })
}

#' @rdname print.efa_control
#' @export
#' @method print efa_rotate_control
print.efa_rotate_control <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.efa_control
#' @export
#' @method format efa_rotate_control
format.efa_rotate_control <- function(x, ...) {
  cli::cli_format_method({
    cli::cli_text("{.strong Rotation control} ({.field type}: {.val {x$type}})")
    cli::cli_text("")
    .control_knob_line("normalize", x$normalize)
    .control_knob_line("precision", x$precision)
    .control_knob_line("order_type", x$order_type)
    .control_knob_line("varimax_type", x$varimax_type)
    .control_knob_line("p_type", x$p_type)
    .control_knob_line("k", x$k)
    .control_knob_line("random_starts", x$random_starts)

    if (length(x$extra_args) > 0) {
      nms <- names(x$extra_args)
      if (is.null(nms)) nms <- rep("", length(x$extra_args))
      vals <- vapply(x$extra_args,
                     function(v) paste(format(v), collapse = ", "), character(1))
      # `pairs` is interpolated as a variable so any user-supplied value cannot be
      # parsed as cli markup.
      pairs <- paste0(nms, " = ", vals)
      cli::cli_text("{.field extra_args}: {pairs}")
    }
  })
}
