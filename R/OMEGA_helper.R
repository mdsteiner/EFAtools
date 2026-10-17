# Percent of uncontaminated correlations (PUC; Reise et al., 2013; Bonifay et al., 2015): the
# share of variable pairs that have no group factor in common. `grp_ind` is the logical
# variables-by-group-factors membership matrix.
.puc <- function(grp_ind) {
  grp_ind <- as.matrix(grp_ind)
  shares_grp <- (grp_ind %*% t(grp_ind)) > 0
  n_items <- nrow(grp_ind)
  cont_corrs <- sum(shares_grp[upper.tri(shares_grp)])
  1 - cont_corrs / (n_items * (n_items - 1) / 2)
}

# H index, 1 / (1 + 1 / sum(lambda^2 / (1 - lambda^2))) (Hancock & Mueller, 2001). NA when a
# loading has |lambda| >= 1, because its error term 1 - lambda^2 is then not positive.
.h_index <- function(lambda) {
  if (length(lambda) == 0L || any(abs(lambda) >= 1, na.rm = TRUE)) return(NA_real_)
  1 / (1 + 1 / sum(lambda^2 / (1 - lambda^2)))
}

# Reliability coefficients of a normalized spec (see the adapters below): omega total,
# hierarchical, and subscale of the whole-scale composite and of each group factor's composite;
# with `add_ind`, the H index, ECV, and PUC; with `add_rel`, standardized alpha.
#
# Every numerator is read off the model-implied common variance Lambda Psi Lambda' (Zinbarg et
# al., 2005, Eq. 8; Zinbarg et al., 2006, Eq. 6). Psi is `Phi` for a correlated-factors spec,
# whose `g_load` is zero, and the identity otherwise. A composite's variance is the sum of its
# block of `cormat` ("correlation"; McDonald, 1999, Eq. 6.2.1) or its model-implied variance
# ("sums_load"): its common variance plus its uniquenesses, or its block of `implied` where the
# spec has residual covariances. `arg` names the map argument of the calling front end.
.reliability_core <- function(spec, variance = c("correlation", "sums_load"),
                              add_ind = TRUE, add_rel = FALSE,
                              arg = "factor_corres") {

  variance <- match.arg(variance)

  g_load <- spec$g_load
  s_load <- spec$s_load
  u2 <- spec$u2
  factor_corres <- spec$map
  cormat <- spec$cormat
  var_names <- spec$var_names
  Phi <- spec$Phi
  implied <- spec$implied

  # The map selects the composites, so only 0 and 1 have a meaning. %in% gives FALSE for NA.
  if (!(is.logical(factor_corres) || is.numeric(factor_corres)) ||
      !all(factor_corres %in% c(0, 1))) {
    vals <- unique(as.vector(as.matrix(factor_corres)))
    other <- utils::head(vals[!vals %in% c(0, 1)], 3L)
    cli::cli_abort(
      c("{.arg {arg}} must hold only 0 and 1.",
        "x" = "It also holds {.val {other}}.",
        "i" = "Mark each variable's group factors with 1 and every other entry with 0."),
      class = "efa_reliability_map_values"
    )
  }

  g <- as.numeric(g_load)
  s_mat <- matrix(as.numeric(unclass(s_load)), nrow = length(g))
  all_idx <- seq_along(g)
  members <- lapply(seq_len(ncol(s_mat)), function(i) which(factor_corres[, i] == 1))

  group_common <- if (is.null(Phi)) tcrossprod(s_mat) else s_mat %*% Phi %*% t(s_mat)
  common <- tcrossprod(g) + group_common

  .rel_check_keying(common, var_names)

  composite_var <- if (variance == "correlation") {
    function(idx) sum(cormat[idx, idx])
  } else if (is.null(implied)) {
    function(idx) sum(common[idx, idx]) + sum(u2[idx])
  } else {
    function(idx) sum(implied[idx, idx])
  }

  # The whole-scale composite holds every variable, so all loadings of a group factor enter its
  # group variance 1' S Psi S' 1, and omega total = omega hierarchical + omega subscale there.
  V <- composite_var(all_idx)
  omega_tot <- sum(common) / V
  omega_h <- sum(g)^2 / V
  omega_sub <- sum(group_common) / V

  # A subscale composite holds the variables assigned to its factor. Its omega total counts
  # the common variance from every factor, including cross-loadings; its omega subscale counts
  # its own factor only.
  for (i in seq_along(members)) {
    idx <- members[[i]]
    V_i <- composite_var(idx)
    omega_tot[i + 1L] <- sum(common[idx, idx]) / V_i
    omega_h[i + 1L] <- sum(g[idx])^2 / V_i
    omega_sub[i + 1L] <- sum(s_mat[idx, i])^2 / V_i
  }

  # A group factor without assigned variables has no composite.
  empty_fac <- which(lengths(members) == 0L)
  if (length(empty_fac) > 0) {
    omega_tot[empty_fac + 1L] <- NA
    omega_h[empty_fac + 1L] <- NA
    omega_sub[empty_fac + 1L] <- NA
    cli::cli_warn(
      c("Some group factors have no assigned items in {.arg {arg}}.",
        "i" = "Their omega coefficients and H index are returned as {.val NA}."),
      class = "efa_omega_empty_factor"
    )
  }

  # A share of variance above 1, or omega hierarchical above omega total, is not admissible.
  # Without a general factor, only omega total is reported (see .rel_drop_general()).
  has_general <- !isTRUE(all(g == 0))
  tol <- .Machine$double.eps^0.5
  out_of_range <- if (has_general) {
    any(c(omega_tot, omega_h, omega_sub) > 1 + tol, na.rm = TRUE) ||
      any(omega_h > omega_tot + tol, na.rm = TRUE)
  } else {
    any(omega_tot > 1 + tol, na.rm = TRUE)
  }
  if (out_of_range) {
    cli::cli_warn(
      c(if (has_general) {
          "Some omega coefficients are out of range (above 1, or omega hierarchical above omega total)."
        } else {
          "Some omega total values are above 1."
        },
        "i" = if (variance == "correlation") {
          "Check that {.arg cormat} is consistent with the loadings."
        } else {
          "Check the uniquenesses of the affected composite for a negative value."
        }),
      class = "efa_omega_out_of_range"
    )
  }

  # Standardized alpha (Cronbach, 1951) of the same composites, from `cormat`, else from the
  # model-implied correlation matrix. NA for fewer than two variables.
  if (isTRUE(add_rel)) {
    R_rel <- if (!is.null(cormat)) {
      cormat
    } else if (!is.null(implied)) {
      implied
    } else {
      stats::cov2cor(common + diag(u2, nrow = length(u2)))
    }
    alpha <- vapply(c(list(all_idx), members), function(idx) {
      k <- length(idx)
      if (k < 2L) return(NA_real_)
      Rk <- R_rel[idx, idx, drop = FALSE]
      k / (k - 1) * (1 - sum(diag(Rk)) / sum(Rk))
    }, numeric(1))
  }

  # Heywood case: a uniqueness at or below 0. which() drops a missing uniqueness.
  heywood_idx <- which(u2 < .Machine$double.eps)
  heywood_vars <- if (is.null(var_names)) as.character(heywood_idx) else
    var_names[heywood_idx]
  if (length(heywood_vars) > 0) {
    cli::cli_warn(
      c(paste("{cli::qty(heywood_vars)}Heywood case{?s} detected for {.val {heywood_vars}}:",
              "a communality at or above 1 leaves a non-positive unique variance."),
        "i" = "{cli::qty(heywood_vars)}Coefficients involving the affected variable{?s} may not be interpretable."),
      class = "efa_reliability_heywood"
    )
  }

  omegas <- cbind(tot = omega_tot, hier = omega_h, sub = omega_sub)

  if (isTRUE(add_ind)) {

    # The H index of a group factor uses the loadings of its assigned variables. `hiload_*`
    # flag a loading |lambda| >= 1 on a variable with a positive uniqueness (an oblique pattern
    # coefficient), whose undefined H is not a Heywood case.
    h_s <- numeric(ncol(s_mat))
    hiload_s <- logical(ncol(s_mat))
    for (j in seq_along(members)) {
      idx <- members[[j]]
      h_s[j] <- .h_index(s_mat[idx, j])
      hiload_s[j] <- any(abs(s_mat[idx, j]) >= 1 & u2[idx] >= .Machine$double.eps,
                         na.rm = TRUE)
    }
    h_g <- .h_index(g)
    hiload_g <- any(abs(g) >= 1 & u2 >= .Machine$double.eps, na.rm = TRUE)

    if (isTRUE(hiload_g) || any(hiload_s)) {
      cli::cli_warn(
        c("The H index could not be computed for some constructs (a loading with an absolute value of at least 1).",
          "i" = "The H index is returned as {.val NA} for those constructs."),
        class = "efa_reliability_h_undefined"
      )
    }

    # ECV (Rodriguez et al., 2016) uses the full loading columns.
    ECV <- sum(g^2) / (sum(g^2) + sum(s_mat^2))
    PUC <- .puc(factor_corres == 1)

    omegas <- cbind(omegas, H = c(h_g, h_s), ECV = NA, PUC = NA)
    omegas[1, "ECV"] <- ECV
    omegas[1, "PUC"] <- PUC
  }

  if (isTRUE(add_rel)) omegas <- cbind(omegas, alpha = alpha)

  rownames(omegas) <- c("g", spec$fac_names)
  class(omegas) <- "OMEGA"
  omegas

}

# Front-end adapters: each normalizes one kind of input to the spec that .reliability_core()
# reads (g_load, s_load, u2, map, cormat, var_names, fac_names, and optionally Phi and
# implied).

# The variable-to-factor map: a supplied `factor_corres`, checked against the loadings, or, for
# `type = "psych"`, each variable's largest absolute group loading (as in psych::omega).
# `type = "EFAtools"` needs an explicit map.
.rel_map <- function(s_load, factor_corres = NULL, type = c("EFAtools", "psych"),
                     arg = "factor_map") {

  type <- match.arg(type)
  s_load <- as.matrix(s_load)

  if (!is.null(factor_corres)) {
    .rel_assert_map_dim(factor_corres, s_load, arg = arg)
    .rel_check_map(s_load, factor_corres, arg = arg)
    return(factor_corres)
  }

  if (type == "EFAtools") {
    cli::cli_abort(
      "Specify {.arg {arg}}, or set {.code type = \"psych\"} to derive variable-to-factor correspondences from the highest group-factor loading per variable.",
      class = "efa_reliability_need_corres"
    )
  }

  map <- matrix(0, nrow = nrow(s_load), ncol = ncol(s_load))
  for (i in seq_len(nrow(s_load))) {
    # a variable without a group loading (e.g., a direct indicator of the general factor)
    # belongs to no group factor
    if (any(abs(s_load[i, ]) > .Machine$double.eps * 100, na.rm = TRUE)) {
      map[i, which.max(abs(s_load[i, ]))] <- 1
    }
  }
  map

}

# Default map of group loadings that can carry structural zeros: the non-zero pattern if most
# variables have an exact zero group loading, as a target matrix has; else, as for an estimated
# solution (also one rounded for print, which has a few incidental zeros), each variable's
# largest loading. A single group column has no choice to make and keeps its pattern.
.rel_default_map <- function(s_load) {
  nonzero <- abs(as.matrix(s_load)) > .Machine$double.eps * 100
  if (ncol(nonzero) <= 1L || mean(rowSums(!nonzero) > 0) >= 0.5) {
    nonzero
  } else {
    .rel_map(s_load, type = "psych")
  }
}

# Abort unless a supplied map has one row per variable and one column per group factor; a
# shorter map would otherwise be recycled.
.rel_assert_map_dim <- function(map, s_load, arg = "factor_map") {

  s_load <- as.matrix(s_load)

  if (!is.matrix(map) || !identical(dim(map), dim(s_load))) {
    got <- if (is.matrix(map)) {
      "{.arg {arg}} has {nrow(map)} row{?s} and {ncol(map)} column{?s}."
    } else {
      "{.arg {arg}} is not a matrix."
    }
    cli::cli_abort(
      c("{.arg {arg}} must have one row per variable and one column per group factor.",
        "x" = got,
        "i" = "The solution has {nrow(s_load)} variable{?s} and {ncol(s_load)} group factor{?s}."),
      class = "efa_reliability_map_dim"
    )
  }

  invisible(map)

}

# Warn about a supplied map whose columns seem permuted: a column whose assigned variables
# hardly load on it (mean |loading| below `min_load`) but do load on another column.
.rel_check_map <- function(s_load, map, min_load = 0.1, arg = "factor_map") {

  s_load <- abs(as.matrix(s_load))
  map <- as.matrix(map)
  if (!identical(dim(map), dim(s_load))) return(invisible(NULL))
  assigned <- map == 1

  labs <- colnames(s_load)
  if (is.null(labs)) labs <- as.character(seq_len(ncol(s_load)))

  implausible <- vapply(seq_len(ncol(assigned)), function(j) {
    idx <- which(assigned[, j])
    if (length(idx) == 0L) return(FALSE)
    own <- mean(s_load[idx, j])
    others <- colMeans(s_load[idx, -j, drop = FALSE])
    isTRUE(own < min_load) && any(others > min_load, na.rm = TRUE)
  }, logical(1))

  if (any(implausible)) {
    bad <- labs[implausible]
    cli::cli_warn(
      c("{cli::qty(bad)}The items assigned to {.arg {arg}} column{?s} {.val {bad}} hardly load on {?that/those} group factor{?s}, but do load on another one.",
        "i" = "The columns of {.arg {arg}} are matched to the group factors by position; check that they are in the column order of the solution."),
      class = "efa_reliability_implausible_map"
    )
  }

  invisible(NULL)

}

# Warn when the variables are not all keyed in the same direction. The coefficients describe
# the raw sum of the variables, so they are not reflected first, as psych::omega does (Flora,
# 2020). The ratio sum(C) / sum(|C|) of the common variance C is 1 when no variable works
# against the sum and falls toward 0 as the keying becomes balanced. It reads C rather than
# loading signs, because a reverse-keyed factor of a correlated-factors solution shows as a
# negative factor correlation, not as negative loadings. The cutoff .75 is below every
# correctly keyed solution of the bundled data (lowest .89) and above every unrecoded one
# (highest .32). Variables with a negative row sum of C are named as the ones to check.
.rel_check_keying <- function(common, var_names = NULL, cutoff = 0.75) {

  aligned <- sum(abs(common))
  if (!is.finite(aligned) || aligned <= .Machine$double.eps^0.5) return(invisible(NULL))

  if (isTRUE(sum(common) < cutoff * aligned)) {

    against <- which(rowSums(common) < 0)
    n_bad <- length(against)

    lead <- if (n_bad > 0) {
      labs <- if (is.null(var_names)) as.character(against) else var_names[against]
      cap <- .cap_label_list(labs)
      paste("The variables are not all keyed in the same direction:",
            "{.val {cap$shown}}{cap$rest} {cli::qty(n_bad)}{?is/are} keyed against the rest.")
    } else {
      "The variables are not all keyed in the same direction."
    }

    cli::cli_warn(
      c(lead,
        "i" = "The coefficients describe the raw unit-weighted sum of the variables as supplied, in which the reversed variables cancel rather than add.",
        "i" = "Reverse-code them before fitting the solution (Flora, 2020).",
        "i" = "If the factors measure opposite constructs, use only the coefficients of each factor. The sum of all variables does not measure one construct."),
      class = "efa_reliability_mixed_keying"
    )
  }

  invisible(NULL)

}

# Warn when supplied uniquenesses do not complete the loadings to unit variance (h2 + u2 = 1),
# most often because communalities were given instead. A warning, not an error: unstandardized
# components give a valid raw-score coefficient. `tol` is above the rounding of a published
# two-decimal table (about .015) and below the error of supplying communalities (at least .40).
.rel_check_u2 <- function(h2, u2, tol = 0.05) {

  if (length(h2) == 0L || length(h2) != length(u2)) return(invisible(NULL))

  resid <- max(abs(h2 + u2 - 1))
  if (!isTRUE(resid > tol)) return(invisible(NULL))

  cli::cli_warn(
    c("The supplied {.arg u2} do not complete the loadings to unit variance.",
      "x" = "A communality plus its uniqueness differs from 1 by up to {round(resid, 3)}.",
      "i" = "A standardized solution has {.code u2 = 1 - h2}; check that {.arg u2} holds uniquenesses rather than communalities, and that they belong to these loadings.",
      "i" = "Unstandardized components are scored as given, and yield the coefficients of the raw composite."),
    class = "efa_reliability_u2_not_standardized"
  )

  invisible(NULL)

}

# Abort unless `Phi` is a k x k correlation matrix: symmetric, unit diagonal, and positive
# semi-definite (within the usual tolerance, as an estimated Phi can be singular to rounding).
# NULL means uncorrelated group factors.
.rel_check_phi <- function(Phi, k, arg = "Phi") {

  if (is.null(Phi)) return(invisible(NULL))

  bad <- function(detail) {
    cli::cli_abort(
      c("{.arg {arg}} must be the correlation matrix of the {k} group factor{?s}.",
        "x" = detail,
        "i" = "Leave it {.val NULL} for uncorrelated group factors."),
      class = "efa_reliability_bad_phi"
    )
  }

  if (!is.matrix(Phi) || !is.numeric(Phi)) {
    bad("It is {.obj_type_friendly {Phi}}, not a numeric matrix.")
  }
  if (nrow(Phi) != k || ncol(Phi) != k) {
    bad("It is {nrow(Phi)} by {ncol(Phi)}, not {k} by {k}.")
  }
  if (anyNA(Phi)) {
    bad("It contains missing values.")
  }

  tol <- .Machine$double.eps^0.5
  if (any(abs(Phi - t(Phi)) > tol)) {
    bad("It is not symmetric.")
  }
  if (any(abs(diag(Phi) - 1) > tol)) {
    bad("Its diagonal is not all ones.")
  }
  if (min(eigen(Phi, symmetric = TRUE, only.values = TRUE)$values) < -tol) {
    bad("It is not positive semi-definite, so it is not a correlation matrix.")
  }

  invisible(NULL)

}

# Check a supplied correlation matrix against the solution and return it as a matrix (a data
# frame is accepted). Named variables in another order are reordered to `var_names`; a different
# set of variables, or a different number, aborts. Unnamed matrices are matched by position.
.rel_check_cormat <- function(cormat, var_names = NULL, n_items = NULL) {

  if (!.is_cormat(cormat)) {
    cli::cli_abort(
      c("{.arg cormat} is not a correlation matrix.",
        "i" = "Check the {.arg cormat} input, or leave it {.code NULL}."),
      class = "efa_reliability_not_cormat"
    )
  }

  cormat <- as.matrix(cormat)

  if (is.null(n_items) && !is.null(var_names)) n_items <- length(var_names)

  if (!is.null(n_items) && nrow(cormat) != n_items) {
    cli::cli_abort(
      c("{.arg cormat} must have one row and column per variable in the solution.",
        "x" = "{.arg cormat} has {nrow(cormat)} variable{?s}, but the solution has {n_items}.",
        "i" = "Supply the correlation matrix of the variables in the solution."),
      class = "efa_reliability_cormat_dim"
    )
  }

  rn <- rownames(cormat)
  cn <- colnames(cormat)
  if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) {
    cli::cli_abort(
      c("The row and column names of {.arg cormat} differ.",
        "i" = "Use the same names on both axes, or remove them to match by position."),
      class = "efa_reliability_cormat_names"
    )
  }
  if (is.null(rn)) rn <- cn

  if (!is.null(var_names) && !is.null(rn) && !identical(rn, var_names)) {
    if (anyDuplicated(rn) == 0L && anyDuplicated(var_names) == 0L &&
        setequal(rn, var_names)) {
      idx <- match(var_names, rn)
      cormat <- cormat[idx, idx, drop = FALSE]
    } else {
      cli::cli_abort(
        c("The variables of {.arg cormat} are not the variables of the solution.",
          "i" = "Supply the correlation matrix of the solution's variables, in any order."),
        class = "efa_reliability_cormat_names"
      )
    }
  }

  cormat
}

# Adapter: an efa_schmid_leiman() (SL) object. The correlation matrix is a supplied `cormat`,
# else the object's `orig_R`; an SL object built from loadings has none, so only
# variance = "sums_load" applies to it.
.rel_adapt_SL <- function(model, factor_corres = NULL, type = "EFAtools",
                          cormat = NULL, fac_names = NULL) {

  sl <- model$sl
  var_names <- rownames(sl)
  g_load <- sl[, 1]
  s_load <- sl[, 2:(ncol(sl) - 2), drop = FALSE]
  u2 <- sl[, "u2"]

  if (is.null(cormat)) {
    if (is.matrix(model$orig_R)) cormat <- model$orig_R
  } else {
    cormat <- .rel_check_cormat(cormat, var_names, nrow(sl))
  }

  if (is.null(fac_names)) fac_names <- colnames(s_load)

  # Derived as for a loading matrix, so that the structural zeros of a transformed lavaan
  # solution keep its pattern, as on the lavaan route.
  map <- if (is.null(factor_corres) && type == "psych") {
    .rel_default_map(s_load)
  } else {
    .rel_map(s_load, factor_corres, type)
  }

  list(g_load = g_load, s_load = s_load, u2 = u2, map = map,
       cormat = cormat, var_names = var_names, fac_names = fac_names)

}

# Adapter: a psych::schmid() object. Without a supplied `cormat`, the correlation matrix is
# rebuilt from the oblique pattern and its factor correlations.
.rel_adapt_schmid <- function(model, factor_corres = NULL, type = "psych",
                              cormat = NULL, fac_names = NULL) {

  pattern <- model$oblique
  Phi <- model$phi
  sl <- model$sl

  s_load_names <- setdiff(colnames(sl[, -1]), c("h2", "u2", "p2", "com"))
  var_names <- rownames(sl)
  g_load <- sl[, 1]
  s_load <- sl[, s_load_names, drop = FALSE]
  u2 <- sl[, "u2"]

  if (is.null(cormat)) {
    cormat <- psych::factor.model(f = pattern, Phi = Phi, U2 = FALSE)
  } else {
    cormat <- .rel_check_cormat(cormat, var_names, nrow(sl))
  }

  if (is.null(fac_names)) fac_names <- s_load_names

  list(g_load = g_load, s_load = s_load, u2 = u2,
       map = .rel_map(s_load, factor_corres, type),
       cormat = cormat, var_names = var_names, fac_names = fac_names)

}

# Adapter: components (g_load, s_load, u2). The correlation matrix is `cormat`, or is rebuilt
# from `pattern` and `Phi` of the oblique solution that a Schmid-Leiman `s_load` came from, or
# stays NULL (only variance = "sums_load" then applies). Without `pattern`, `Phi` is the factor
# correlation matrix of `s_load` itself: the components are then a correlated-factors solution,
# which has no general factor.
.rel_adapt_manual <- function(g_load, s_load, u2, var_names, factor_corres = NULL,
                              type = "EFAtools", cormat = NULL, pattern = NULL,
                              Phi = NULL, fac_names = NULL) {

  # unclass() first: as.matrix() keeps the class of a Schmid-Leiman loading table.
  s_load <- as.matrix(unclass(s_load))

  if (!is.null(cormat) && !is.null(pattern)) {
    cli::cli_abort(
      c("{.arg cormat} and {.arg pattern} are two ways to give the same correlation matrix, so supply one or the other.",
        "i" = "Drop {.arg pattern} to score the components against {.arg cormat}.",
        "i" = "Drop {.arg cormat} to reconstruct one instead, from {.arg pattern} and the {.arg Phi} of the oblique solution it came from."),
      class = "efa_reliability_cormat_and_pattern"
    )
  }

  phi_solution <- !is.null(Phi) && is.null(pattern)

  if (!is.null(Phi)) {
    .rel_check_phi(Phi, if (phi_solution) ncol(s_load) else NCOL(pattern))
  }

  # A general factor and correlated group factors do not partition the composite variance.
  if (phi_solution && !isTRUE(all(g_load == 0))) {
    cli::cli_abort(
      c("{.arg Phi} describes correlated group factors, which cannot be combined with a general factor.",
        "i" = "Supply {.arg Phi} for a correlated-factors solution, whose {.arg g_load} is zero throughout.",
        "i" = "For a Schmid-Leiman or bifactor solution the group factors are orthogonal; pass {.arg pattern} as well if {.arg Phi} belongs to the oblique solution it came from."),
      class = "efa_reliability_phi_with_general"
    )
  }

  # Without a general factor, a `pattern` that repeats `s_load` looks like a parent oblique
  # solution, so state the reading that is taken.
  if (!is.null(Phi) && !is.null(pattern) && isTRUE(all(g_load == 0))) {
    cli::cli_warn(
      c("{.arg Phi} is read as the factor correlations of {.arg pattern}, not of {.arg s_load}.",
        "i" = "It is used to reconstruct the correlation matrix only, so the coefficients describe uncorrelated group factors.",
        "i" = "Omit {.arg pattern} to score {.arg s_load} as a correlated-factors solution instead."),
      class = "efa_reliability_phi_pattern"
    )
  }

  if (is.null(cormat)) {
    if (!is.null(Phi) && !is.null(pattern)) {
      cormat <- psych::factor.model(f = pattern, Phi = Phi, U2 = FALSE)
    }
  } else {
    # Matched to the loadings' row names; `var_names` only labels the output.
    cormat <- .rel_check_cormat(cormat, rownames(s_load), nrow(s_load))
  }

  group_h2 <- if (phi_solution) rowSums((s_load %*% Phi) * s_load) else rowSums(s_load^2)
  .rel_check_u2(g_load^2 + group_h2, u2)

  if (is.null(fac_names)) fac_names <- seq_len(ncol(s_load))

  # Derived as for a loading matrix, so that a target with exact zeros keeps its pattern.
  map <- if (is.null(factor_corres) && type == "psych") {
    .rel_default_map(s_load)
  } else {
    .rel_map(s_load, factor_corres, type)
  }

  list(g_load = g_load, s_load = s_load, u2 = u2, map = map,
       Phi = if (phi_solution) Phi else NULL,
       cormat = cormat, var_names = var_names, fac_names = fac_names)

}

# Default map of a `lavaan` loading matrix: its non-zero pattern, except that in an ESEM block
# (`efa()` in the model syntax), whose loadings are all free, each variable is assigned to the
# block's factor with its largest loading. `blocks` holds the ESEM block of each column, "" for
# none.
.rel_lavaan_map <- function(lambda, blocks) {
  map <- abs(lambda) > .Machine$double.eps * 100
  for (b in setdiff(unique(blocks), "")) {
    cols <- which(blocks == b)
    for (r in which(rowSums(map[, cols, drop = FALSE]) > 0)) {
      map[r, cols[-which.max(abs(lambda[r, cols]))]] <- FALSE
    }
  }
  map
}

# Adapter: a `lavaan` fit, one spec per group, scored with the model-implied composite variance.
# A fit is a bifactor solution when variables load on a factor named `g_name` and on another
# factor; a second-order solution (structure in `beta`) is Schmid-Leiman transformed; a
# single-factor group is flagged `single = TRUE`; every other fit is a set of correlated factors
# (`correlated = TRUE`), with or without cross-loadings. `implied` is read from lavaan only for
# a group with residual covariances, which the diagonal of Theta misses.
.rel_adapt_lavaan <- function(model, g_name = "g", group_names = NULL) {

  .require_lavaan()

  if (isFALSE(lavaan::lavInspect(model, what = "converged"))) {
    cli::cli_abort("The model did not converge; no reliability coefficients are computed.",
                   class = "efa_reliability_no_converge")
  }

  std_sol <- suppressWarnings(lavaan::lavInspect(model, what = "std",
                                drop.list.single.group = FALSE))
  est <- lavaan::lavInspect(model, what = "est", drop.list.single.group = FALSE)
  implied_cor <- NULL

  # The ESEM block of each factor; a factor outside every block has none.
  pt <- lavaan::parTable(model)
  esem <- if (!is.null(pt[["efa"]])) {
    in_block <- pt[["op"]] == "=~" & nzchar(pt[["efa"]])
    stats::setNames(pt[["efa"]][in_block], pt[["lhs"]][in_block])
  }
  blocks_of <- function(cols) {
    b <- unname(esem[cols])
    if (is.null(b)) b <- rep("", length(cols))
    b[is.na(b)] <- ""
    b
  }

  if (is.null(group_names)) {
    group_names <- names(std_sol)
  } else if (length(group_names) != length(std_sol)) {
    cli::cli_abort(
      c("{.arg group_names} does not match the number of groups in the {.cls lavaan} model.",
        "i" = "The model has {length(std_sol)} group{?s}, but {.arg group_names} has {length(group_names)}."),
      class = "efa_reliability_group_names"
    )
  }

  tol <- .Machine$double.eps * 100
  higherorder <- FALSE
  correlated <- FALSE
  few_loadings <- FALSE
  groups <- vector("list", length(std_sol))

  for (i in seq_along(std_sol)) {

    # A factor without positive (residual) variance leaves no standardized solution. A
    # non-positive residual variance of a variable is a Heywood case, which the core reports.
    if (any(diag(est[[i]][["psi"]]) <= 0)) {
      cli::cli_abort(
        "A factor has a variance or residual variance of 0 or less; no reliability coefficients are computed.",
        class = "efa_reliability_factor_heywood"
      )
    }

    lambda <- std_sol[[i]][["lambda"]]
    theta <- std_sol[[i]][["theta"]]
    psi <- std_sol[[i]][["psi"]]

    if (any(is.na(lambda))) {
      cli::cli_abort("Some loadings are {.val NA} or {.val NaN}; no reliability coefficients are computed.",
                     class = "efa_reliability_na_loadings")
    }

    var_names_i <- rownames(lambda)

    implied_i <- if (any(abs(theta[upper.tri(theta)]) > tol, na.rm = TRUE)) {
      if (is.null(implied_cor)) {
        implied_cor <- lavaan::lavInspect(model, what = "cor.ov",
                                          drop.list.single.group = FALSE)
      }
      implied_cor[[i]]
    }

    if (ncol(lambda) == 1) {
      # `fac_label` is kept apart from the spec, whose `fac_names` label group-factor rows.
      groups[[i]] <- c(
        list(single = TRUE, fac_label = colnames(lambda)),
        .rel_single_factor_spec(list(g_load = lambda[, 1],
                                     s_load = lambda[, 0, drop = FALSE],
                                     u2 = diag(theta), cormat = NULL,
                                     implied = implied_i,
                                     var_names = var_names_i)))
      next
    }

    col_names <- colnames(lambda)

    # All groups share the structure, so detect it once.
    if (i == 1) {

      beta <- std_sol[[i]][["beta"]]
      higherorder <- !is.null(beta) && any(abs(beta) > tol, na.rm = TRUE)
      on <- abs(lambda) > tol
      has_g <- g_name %in% col_names
      bifactor <- has_g && any(on[, g_name] & rowSums(on) >= 2)
      correlated <- !higherorder && !bifactor

      # A general factor under another name: outside the ESEM blocks, a factor that has
      # variables in common with every other factor, and that loads on every variable while
      # other factors do not, or (with three or more factors) is uncorrelated with them. A
      # second-order model needs `g_name` to find its general factor.
      cfa <- on[, blocks_of(col_names) == "", drop = FALSE]
      full <- colSums(cfa) == nrow(cfa)
      general_like <- ncol(cfa) > 1L && any(colSums(crossprod(cfa) > 0) == ncol(cfa))
      orthogonal <- all(abs(psi[upper.tri(psi)]) <= tol, na.rm = TRUE)
      misnamed <- general_like &&
        ((any(full) && !all(full)) || (orthogonal && ncol(cfa) >= 3L))
      if (!has_g && (higherorder || misnamed)) {
        cli::cli_abort(
          c("The {.cls lavaan} solution has no factor named {.val {g_name}}.",
            "i" = if (higherorder) {
              "The fit is a second-order model; set {.arg g_name} to the name of its second-order factor."
            } else {
              "One factor has variables in common with every other factor, so the fit looks like a bifactor model; set {.arg g_name} to the name of its general factor."
            }),
          class = "efa_reliability_g_name"
        )
      }

      if (higherorder) {
        if (sum(colSums(abs(beta) > tol, na.rm = TRUE) > 0) > 1) {
          cli::cli_abort("The higher-order model has more than two latent strata or more than one second-order factor; only second-order models with one second-order factor are supported.",
                         class = "efa_reliability_higher_order")
        }
      } else if (!correlated) {
        few_loadings <- !all(rowSums(on) > 1)
      }
    }

    if (isTRUE(correlated)) {
      # The standardized `psi` is the factor correlation matrix, indexed by name to match the
      # columns of `lambda`.
      Phi_i <- if (is.null(col_names)) psi else psi[col_names, col_names, drop = FALSE]
      groups[[i]] <- list(single = FALSE, g_load = rep(0, nrow(lambda)),
                          s_load = lambda, u2 = diag(theta), Phi = Phi_i,
                          map = .rel_lavaan_map(lambda, blocks_of(col_names)),
                          cormat = NULL, implied = implied_i,
                          var_names = var_names_i, fac_names = col_names)
      next
    }

    # The decomposition needs uncorrelated factors: a bifactor model fitted with orthogonal
    # factors, or a second-order model whose first-order covariances run through `beta`.
    if (any(abs(psi[upper.tri(psi)]) > tol, na.rm = TRUE)) {
      cli::cli_abort(
        c("The factors of the {.cls lavaan} solution are correlated; the coefficients need uncorrelated factors.",
          "i" = "A bifactor model needs {.code orthogonal = TRUE}; a second-order model needs no freed covariances between the first-order factors.",
          "i" = "If {.val {g_name}} is not a general factor, set {.arg g_name} to a name that the model does not use."),
        class = "efa_reliability_correlated_factors"
      )
    }

    col_names <- col_names[!col_names %in% g_name]

    if (isTRUE(higherorder)) {
      # Schmid-Leiman transform: general loadings are the direct loadings of the general
      # factor plus the first-order loadings times the second-order loadings, computed
      # before the group loadings are overwritten.
      lambda[, g_name] <- lambda[, g_name] +
        lambda[, col_names] %*% std_sol[[i]][["beta"]][col_names, g_name]
      lambda[, col_names] <- .sl_group_loadings(lambda[, col_names], psi, col_names)
    }

    s_load <- lambda[, col_names, drop = FALSE]

    groups[[i]] <- list(single = FALSE, g_load = lambda[, g_name], s_load = s_load,
                        u2 = diag(theta), map = .rel_lavaan_map(s_load, blocks_of(col_names)),
                        cormat = NULL, implied = implied_i,
                        var_names = var_names_i, fac_names = col_names)

  }

  list(groups = groups, group_names = group_names, variance = "sums_load",
       higher_order = higherorder, few_loadings = few_loadings,
       correlated = correlated)

}

# Blank the cells that a solution without a general factor does not define. Omega
# hierarchical, ECV, and PUC need a general factor. Omega subscale is the share in addition to
# a general factor (Rodriguez et al., 2016), and once the factors correlate it is not a share of
# the composite variance at all, because the cross terms belong to no single factor. The
# whole-scale row has no factor, so no H index.
.rel_drop_general <- function(x) {
  x[, c("hier", "sub", "ECV", "PUC")] <- NA
  x["g", "H"] <- NA
  x
}

# The spec of a solution with exactly one factor, or NULL if `spec` has more. The factor is
# moved to `g_load`, whichever slot it came in (the general column of a one-column matrix, or the
# only group factor of a one-factor efa_fit() solution), so every route gives the same
# coefficients. A 1 x 1 `Phi` carries no information and is dropped.
.rel_single_factor_spec <- function(spec) {

  s_load <- as.matrix(spec$s_load)
  n_g <- if (isTRUE(all(spec$g_load == 0))) 0L else 1L
  if (n_g + ncol(s_load) != 1L) return(NULL)

  g_load <- if (n_g == 1L) spec$g_load else s_load[, 1]
  none <- matrix(numeric(0), nrow = length(g_load), ncol = 0)

  list(g_load = g_load, s_load = none, u2 = spec$u2, map = none, Phi = NULL,
       cormat = spec$cormat, implied = spec$implied, var_names = spec$var_names,
       fac_names = character(0))

}

# The row label of a single factor: its own name, else "F1" (not "g": one factor is not the
# general factor of a hierarchy). A position, as the manual route gives, is not a name.
.rel_single_factor_label <- function(fac_names) {
  if (length(fac_names) != 1L || is.na(fac_names) || is.numeric(fac_names) ||
      !nzchar(as.character(fac_names))) {
    return("F1")
  }
  as.character(fac_names)
}

# Blank the cells that a single-factor solution does not define, and label its row. Omega
# subscale needs group factors, omega hierarchical equals omega total, and ECV and PUC are 1
# by construction.
.rel_drop_single_factor <- function(x, fac_names = NULL) {
  x[, c("hier", "sub", "ECV", "PUC")] <- NA
  rownames(x) <- .rel_single_factor_label(fac_names)
  x
}

# Adapter: an efa_fit() solution. An oblique or orthogonal solution is a correlated-factors
# spec (zero general factor, the rotated pattern, and Phi, the identity for an orthogonal
# rotation), not a Schmid-Leiman transform, which is underidentified for two factors and biased
# by cross-loadings (Flora, 2020; Mansolf & Reise, 2016). A "bifactorT" solution is a bifactor
# spec, with the factor that loads on every variable as the general factor. A single-factor
# solution is read from its unrotated loadings.
.rel_adapt_efa <- function(model, factor_corres = NULL, type = "psych",
                           cormat = NULL, fac_names = NULL) {

  L_unrot <- unclass(model$unrot_loadings)
  single <- NCOL(L_unrot) == 1L
  rotation <- model$settings$rotation %||% ""

  if (!single && rotation == "none") {
    cli::cli_abort(
      c("{.arg model} is an unrotated solution with {NCOL(L_unrot)} factors.",
        "i" = "Refit it with a rotation."),
      class = "efa_reliability_not_rotated"
    )
  }
  if (!single && rotation == "bifactorQ") {
    cli::cli_abort(
      c("{.arg model} is an oblique bifactor solution ({.code rotation = \"bifactorQ\"}), whose group factors are correlated.",
        "i" = "The coefficients need uncorrelated group factors: refit with {.code rotation = \"bifactorT\"}, or transform an oblique solution with {.fn efa_schmid_leiman}."),
      class = "efa_reliability_bifactor_oblique"
    )
  }

  orig_R <- if (is.matrix(model$orig_R)) model$orig_R

  # The rotation orders the factors by their sums of squares, so the general factor need not
  # come first. It is the factor that loads on every variable: the column whose smallest
  # absolute loading is largest, as for the analytic standard errors.
  if (!single && rotation == "bifactorT") {
    L <- unclass(model$rot_loadings)
    g_col <- which.max(apply(abs(L), 2, min))
    return(.rel_adapt_bifactor(L[, c(g_col, seq_len(ncol(L))[-g_col]), drop = FALSE],
                               factor_corres = factor_corres,
                               cormat = cormat %||% orig_R, fac_names = fac_names))
  }

  L1 <- if (single) L_unrot else unclass(model$rot_loadings)
  Phi <- model[["Phi"]]
  if (single || is.null(Phi)) Phi <- diag(ncol(L1))

  s_load <- L1

  common <- s_load %*% Phi %*% t(s_load)
  u2 <- 1 - diag(common)

  cormat <- if (!is.null(cormat)) {
    .rel_check_cormat(cormat, rownames(s_load), nrow(s_load))
  } else {
    orig_R %||% (common + diag(u2, nrow = length(u2)))
  }

  if (is.null(fac_names)) {
    fac_names <- colnames(s_load)
    if (is.null(fac_names)) fac_names <- paste0("F", seq_len(ncol(s_load)))
  }

  list(g_load = rep(0, nrow(s_load)), s_load = s_load, u2 = u2,
       map = .rel_map(s_load, factor_corres, type), Phi = Phi,
       cormat = cormat, var_names = rownames(s_load), fac_names = fac_names)

}

# Adapter: a bifactor loading matrix (general factor first). Uniquenesses default to
# 1 - rowSums(L^2); the map defaults to .rel_default_map(); without a `cormat`, the model-implied
# correlation matrix is used.
.rel_adapt_bifactor <- function(loadings, factor_corres = NULL, u2 = NULL,
                                cormat = NULL, fac_names = NULL) {

  loadings <- as.matrix(loadings)
  g_load <- loadings[, 1]
  s_load <- loadings[, -1, drop = FALSE]

  if (is.null(u2)) {
    u2 <- 1 - rowSums(loadings^2)
  } else {
    .rel_check_u2(rowSums(loadings^2), u2)
  }

  if (is.null(factor_corres)) {
    factor_corres <- .rel_default_map(s_load)
  } else {
    .rel_assert_map_dim(factor_corres, s_load, arg = "factor_map")
    .rel_check_map(s_load, factor_corres, arg = "factor_map")
  }

  if (is.null(cormat)) {
    cormat <- loadings %*% t(loadings) + diag(u2, nrow = length(u2))
  } else {
    # Matched to the loadings' own row names, so an unnamed matrix is matched by position.
    cormat <- .rel_check_cormat(cormat, rownames(loadings), nrow(loadings))
  }

  var_names <- rownames(loadings)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(nrow(loadings)))

  if (is.null(fac_names)) {
    # A one-column matrix is a single factor, named by its only column.
    fac_names <- if (ncol(s_load) == 0L) colnames(loadings) else colnames(s_load)
    if (is.null(fac_names)) fac_names <- seq_len(ncol(s_load))
  }

  list(g_load = g_load, s_load = s_load, u2 = u2, map = factor_corres,
       cormat = cormat, var_names = var_names, fac_names = fac_names)

}

# Map from a .reliability_core() column to its coefficient name, kind, and printed label in an
# efa_reliability result. The row order sets the order of the coefficients.
.reliability_registry <- function() {
  data.frame(
    core        = c("tot", "hier", "sub", "alpha", "H", "ECV", "PUC"),
    coefficient = c("omega_total", "omega_hierarchical", "omega_subscale",
                    "alpha", "H", "ECV", "PUC"),
    kind        = c("reliability", "reliability", "reliability", "reliability",
                    "reliability", "common_variance", "common_variance"),
    label       = c("tot", "hier", "sub", "alpha", "H", "ECV", "PUC"),
    stringsAsFactors = FALSE
  )
}

# Shape coefficient matrices into the long efa_reliability result. `x` is one matrix (an
# ungrouped solution) or a named list of matrices (one per group). The first row of a matrix is
# the general factor (or the single factor); NA cells are dropped.
.reliability_result <- function(x, settings = NULL) {

  if (is.list(x)) {
    groups <- x
    group_names <- names(x)
    if (is.null(group_names)) group_names <- rep("", length(x))
    blank <- !nzchar(group_names)
    group_names[blank] <- as.character(seq_along(x))[blank]
  } else {
    groups <- list(x)
    group_names <- NA_character_
  }

  reg <- .reliability_registry()

  parts <- lapply(seq_along(groups), function(g) {
    mat <- as.matrix(groups[[g]])
    factors <- rownames(mat)
    if (is.null(factors)) factors <- paste0("F", seq_len(nrow(mat)))
    keep <- reg[reg$core %in% colnames(mat), , drop = FALSE]

    # The level is read by position, so a group factor named "g" stays at level "group".
    row_level <- ifelse(seq_along(factors) == 1L, "general", "group")

    # Without a general factor, the first row describes the sum of all variables.
    if (isTRUE(settings$no_general) && length(factors) > 0L) {
      factors[1L] <- "total"
      row_level[1L] <- "total"
    }

    data.frame(
      coefficient = rep(keep$coefficient, each = length(factors)),
      level = rep(row_level, times = nrow(keep)),
      factor = rep(factors, times = nrow(keep)),
      group = rep(group_names[g], length(factors) * nrow(keep)),
      value = as.numeric(mat[, keep$core]),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, parts)
  out <- out[!is.na(out$value), , drop = FALSE]
  rownames(out) <- NULL

  kind <- stats::setNames(reg$kind, reg$coefficient)[unique(out$coefficient)]

  attr(out, "settings") <- settings
  attr(out, "kind") <- kind
  class(out) <- c("efa_reliability", "data.frame")
  out

}

# Flexible omega function (e.g. to use with loadings obtained by MacOrtho)------
.OMEGA_FLEX <- function(model = NULL, type = c("EFAtools", "psych"),
                        factor_corres = NULL,
                        var_names = NULL, fac_names = NULL, g_load = NULL,
                        s_load = NULL, u2 = NULL, cormat = NULL, pattern = NULL,
                        Phi = NULL, variance = c("correlation", "sums_load"),
                        add_ind = TRUE){

  if(inherits(model, "schmid")){

    pattern <- model$oblique
    Phi <- model$phi

    model <-  model$sl
    s_load_names <- setdiff(colnames(model[, -1]),
                            c("h2", "u2", "p2", "com"))

    var_names <- rownames(model)
    g_load <- model[, 1]
    s_load <- model[, s_load_names]
    u2 <- model[, "u2"]

 } else if(inherits(model, "SL")){

    # Honour a user-supplied cormat; only fall back to the SL object's stored
    # correlation matrix when none was given. Flexible-input SL objects store
    # orig_R = NA, in which case cormat stays NULL and the checks below apply.
    if(is.null(cormat) && is.matrix(model$orig_R)){
      cormat <- model$orig_R
    }

    model <-  model$sl
    var_names <- rownames(model)
    g_load <- model[, 1]
    s_load <- model[, 2:(ncol(model) - 2)]
    u2 <- model[, "u2"]

  }

  # Same general + group factor labels regardless of how s_load was obtained.
  factor_names <- c("g", seq_len(ncol(s_load)))

    if(variance == "correlation"){

      if(is.null(cormat)){

        if(is.null(Phi) | is.null(pattern)) {
          cli::cli_abort(
            c("Specify either {.arg cormat}, or {.arg Phi} and {.arg pattern}.",
              "i" = "Alternatively, set {.code variance = \"sums_load\"}."),
            class = "efa_omega_need_cormat"
          )

        } else {

          # Create the correlation matrix from the pattern coefficients and factor
          # intercorrelations
          cormat <- psych::factor.model(f = pattern, Phi = Phi, U2 = FALSE)

        }

      } else {

        # Check if it is a correlation matrix
        if(!.is_cormat(cormat)) {

          cli::cli_abort(
            c("{.arg cormat} is not a correlation matrix.",
              "i" = "Check the {.arg cormat} input, supply {.arg Phi} and {.arg pattern} instead, or set {.code variance = \"sums_load\"}."),
            class = "efa_omega_not_cormat"
          )

        }

      }

    }

  # Check if input to factor_corres is correct (g_load is a vector here, so the
  # item count is its length, equivalently nrow(s_load)).
  checkmate::assert_matrix(factor_corres, null.ok = TRUE, nrows = length(g_load),
                           ncols = ncol(s_load))

  # Create an input dataframe
  input <- data.frame(g_load, s_load)
  colnames(input) <- factor_names
  rownames(input) <- var_names

  if(type == "EFAtools" & is.null(factor_corres)){

    cli::cli_abort("Specify {.arg factor_corres}, or set {.code type = \"psych\"} to derive variable-to-factor correspondences from the highest group-factor loading per variable.",
                   class = "efa_omega_need_corres")

  }

  if(type == "psych"){

    if(variance != "correlation"){

      cli::cli_warn(
        c("{.arg variance} is specified; the value {.val {variance}} is used.",
          "i" = "Results may differ from the specified {.arg type}."),
        class = "efa_omega_variance_override"
      )
      }

    if(is.null(factor_corres)){

      factor_corres <- .rel_map(s_load, type = "psych")


    } else {

      cli::cli_warn(
        c("{.arg factor_corres} is specified; the supplied variable-to-factor correspondences are used.",
          "i" = "To compute correspondences as in psych, leave {.code factor_corres = NULL}."),
        class = "efa_omega_corres_override"
      )
    }
  }

  # Group-factor row labels: user-supplied names, else the factor names carried by a
  # model object, else integer positions for manually specified components.
  fac_names_out <- if (!is.null(fac_names)) {
    fac_names
  } else if (is.null(model)) {
    seq_len(ncol(s_load))
  } else {
    colnames(model)[seq_len(ncol(s_load)) + 1L]
  }

  # Hand the normalized components to the reliability engine.
  spec <- list(g_load = g_load, s_load = s_load, u2 = u2, map = factor_corres,
               cormat = cormat, var_names = var_names, fac_names = fac_names_out)

  .reliability_core(spec, variance = variance, add_ind = add_ind,
                    arg = "factor_corres")

}


# Reshape a lavaan reliability solution into OMEGA's legacy output ------------
#
# OMEGA's lavaan path is a thin front-end over the shared reliability machinery:
# the lavaan adapter normalizes the fitted model to a per-group spec (Schmid-
# Leiman transforming a second-order model) and the reliability core scores every
# group, single-factor ones included. This function only reassembles those results
# into OMEGA's historical shapes -- a coefficient matrix per group (unwrapped for a
# single group), a named list for several groups, and a named c(Omega, H) vector
# (or a bare omega) for a single factor -- and emits the user-facing notes about the
# model structure. All reliability math lives in .reliability_core.
#
# A single factor defines standardized alpha as well, which the core computes and
# efa_reliability() surfaces. OMEGA's single-factor output is a named c(Omega, H) vector
# rather than a coefficient matrix, so reporting alpha would change the shape of a
# superseded function's return value and not only its coefficient menu; it keeps the two
# coefficients it has always returned.
.OMEGA_LAVAAN <- function(model = NULL, g_name = "g", group_names = NULL,
                          add_ind = TRUE){

  adapt <- .rel_adapt_lavaan(model, g_name = g_name, group_names = group_names)

  # OMEGA's input is a Schmid-Leiman, bifactor, second-order, or single-factor solution,
  # each of which has a general factor, and its wide per-factor output reports that
  # factor's omega hierarchical, ECV, and PUC. A correlated-factors fit defines none of
  # them, so it is refused here rather than scored into a table whose general-factor
  # columns would be zeros. efa_reliability() scores it and omits those coefficients.
  if (isTRUE(adapt$correlated)) {
    cli::cli_abort(
      c("The lavaan input is invalid; no reliability coefficients are computed.",
        "i" = "Provide a bifactor model, a second-order model, or a single-factor model.",
        "i" = "To score a correlated-factors solution, use {.fn efa_reliability}."),
      class = "efa_reliability_invalid_lavaan"
    )
  }

  # A second-order general factor was Schmid-Leiman transformed by the adapter.
  if (isTRUE(adapt$higher_order)) {
    cli::cli_inform(
      c("i" = "The specified general factor is a second-order factor; omegas are computed on the Schmid-Leiman transformed second-order solution."),
      class = "efa_omega_g_second_order"
    )
  }

  # A supplied "bifactor" in which some item loads on fewer than two factors.
  if (isTRUE(adapt$few_loadings)) {
    cli::cli_inform(
      c("i" = "Some variables have fewer than two loadings; did you enter a bifactor model? Provide a bifactor model, a second-order model, or a single-factor model."),
      class = "efa_omega_few_loadings"
    )
  }

  omegas <- vector("list", length(adapt$groups))
  informed_single <- FALSE

  for (i in seq_along(adapt$groups)) {

    grp <- adapt$groups[[i]]

    if (isTRUE(grp$single)) {

      # A single factor: omega total (and, with add_ind, the H index) only.
      if (!informed_single) {
        msg <- if (isTRUE(add_ind)) {
          "The model contained a single factor; only omega total and the H index are returned."
        } else {
          "The model contained a single factor; only omega total is returned."
        }
        cli::cli_inform(c("i" = msg), class = "efa_omega_single_factor")
        informed_single <- TRUE
      }

      # Scored through the core on the spec the adapter normalized it to, as every other
      # group is; only the two coefficients OMEGA reports for a single factor are kept.
      sf <- unclass(.reliability_core(grp, "sums_load", add_ind = add_ind,
                                      arg = "factor_corres"))
      omegas[[i]] <- if (isTRUE(add_ind)) {
        stats::setNames(c(sf["g", "tot"], sf["g", "H"]), c("Omega", "H"))
      } else {
        sf["g", "tot"]
      }

    } else {

      # A multi-factor (bifactor / Schmid-Leiman) group: score through the core,
      # storing the unclassed matrix -- the OMEGA class is attached to the whole
      # output below, matching the historical list-of-matrices shape. The core
      # labels the general-factor row "g"; restore the user's general-factor name.
      mat <- unclass(.reliability_core(grp, "sums_load", add_ind = add_ind,
                                       arg = "factor_corres"))
      rownames(mat)[1] <- g_name
      omegas[[i]] <- mat

    }
  }

  if (length(omegas) > 1) {
    names(omegas) <- adapt$group_names
  } else {
    omegas <- omegas[[1]]
  }

  class(omegas) <- "OMEGA"
  omegas

}
