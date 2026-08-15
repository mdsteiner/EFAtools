# Percent of uncontaminated correlations (PUC): the proportion of item pairs that
# do not share a group factor (Reise et al., 2013; Bonifay et al., 2015). `grp_ind`
# is a logical items-by-group-factors membership matrix (the general factor is
# excluded, as it contaminates no correlation). Each item pair sharing at least one
# group factor is a single contaminated correlation, regardless of how many group
# factors it shares.
.puc <- function(grp_ind) {
  grp_ind <- as.matrix(grp_ind)
  shares_grp <- (grp_ind %*% t(grp_ind)) > 0
  n_items <- nrow(grp_ind)
  cont_corrs <- sum(shares_grp[upper.tri(shares_grp)])
  1 - cont_corrs / (n_items * (n_items - 1) / 2)
}

# H index over a construct's own loadings: the correlation between an optimally-weighted
# composite and the factor, 1 / (1 + 1 / sum(lambda^2 / (1 - lambda^2))) (Hancock & Mueller,
# 2001). A loading with |lambda| >= 1 makes the error term 1 - lambda^2 <= 0 and the summand
# undefined, so H is not interpretable for that construct; return NA -- as for a construct
# with no indicators -- rather than a negative or infinite value. For standardized loadings
# this only happens in an improper (Heywood) solution, but an oblique pattern coefficient can
# exceed 1 in a proper solution; callers flag genuine Heywood cases separately, via the unique
# variances. A missing loading (NA) is excluded from the >= 1 test and flows through the sum
# to yield NA.
.h_index <- function(lambda) {
  if (length(lambda) == 0L || any(abs(lambda) >= 1, na.rm = TRUE)) return(NA_real_)
  1 / (1 + 1 / sum(lambda^2 / (1 - lambda^2)))
}

# Reliability coefficients over a normalized model spec. `spec` carries the general-factor
# loadings (`g_load`), the group-factor loading matrix (`s_load`), the uniquenesses (`u2`),
# the item-to-factor correspondence map (`map`), the group-factor row labels (`fac_names`),
# the item names (`var_names`), and -- for `variance = "correlation"` -- a correlation
# matrix (`cormat`) plus, for a correlated-factors solution, the factor
# intercorrelations (`Phi`, in the column order of `s_load`). `Phi` belongs to a spec
# with no general factor: it marks the group loadings as an oblique pattern whose
# factors covary, and enters the model-implied common variance in either variance mode.
# A spec that carries both a general factor and a `Phi` is not supported -- the callers that
# build a `Phi` set `g_load` to zero, and reject the combination where a user could supply
# one, because the general and group parts would no longer partition the composite.
# Given these, it computes
# McDonald's omega total, hierarchical, and subscale for the general factor and each group
# factor, and, when `add_ind = TRUE`, the H index (Hancock & Mueller, 2001) together with the
# ECV and PUC bifactor indices (Rodriguez, Reise & Haviland, 2016). Every omega numerator is
# the true score variance the model attributes to a composite, read off the model-implied
# common variance Lambda Psi Lambda'. Two total-variance denominators are supported:
# `"correlation"` takes a composite's variance from the correlation matrix, as McDonald's
# omega total does (1999, Test Theory, Eq. 6.2.1; Zinbarg, Yovel, Revelle & McDonald, 2006,
# Eqs. 4 and 6); `"sums_load"` uses the model-implied composite variance instead, the common
# variance plus the uniquenesses, and so needs no correlation matrix. When
# `add_rel = TRUE`, three further columns are appended: standardized Cronbach's alpha
# (Cronbach, 1951) for the whole scale and each subscale, and -- per group factor, over its
# assigned (simple-structure) composite -- the composite reliability (congeneric omega;
# Joreskog, 1971; Raykov, 2001) and the average variance extracted (AVE, a convergent-validity
# index rather than a reliability; Fornell & Larcker, 1981). The function is purely
# computational: `spec` is assumed already normalized by the calling adapter, which owns all
# front-end input handling. `arg` names the user-facing map argument of the calling
# front-end, so the empty-factor warning points at the right one: `factor_map` for
# efa_reliability(), `factor_corres` on the frozen OMEGA() surface. Both front-ends pass
# it explicitly; the default is the frozen spelling so that a caller who forgets leaves
# OMEGA()'s wording as it has always been rather than silently changing it.
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

  # The map states membership, so only 0 and 1 (or the logicals they stand for) mean
  # anything. Any other value is read two ways below and cannot be read one way: the
  # composites are the entries equal to 1, while the map is also multiplied into the
  # loadings and the uniquenesses, where it weights them. A map of 2s therefore names the
  # same composites but empties every one of them -- the coefficients come back NA under a
  # whole-scale row that still looks right -- so it is refused rather than resolved.
  # Checked here rather than in the front ends because this is where every route meets: the
  # derived maps are 0/1 or logical by construction, but a user's arrives through several
  # callers, not all of which validate it.
  # %in% is match-based, so a missing entry answers FALSE rather than propagating.
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

  # Same general + group factor labels regardless of how s_load was obtained.
  factor_names <- c("g", seq_len(ncol(s_load)))

  # Create an input dataframe
  input <- data.frame(g_load, s_load)
  colnames(input) <- factor_names
  rownames(input) <- var_names
  input$u2 <- u2

  # seq_len(k) + 1 rather than 2:(k + 1): the latter counts down to c(2, 1) when a spec
  # arrives with no group factors, selecting the general column twice over.
  omega_mat <- input[, seq_len(ncol(s_load)) + 1L, drop = FALSE] * factor_corres

  # Sum of all g loadings
  sum_g <- sum(input$g)

  # Sum of all error variances
  sum_e <- sum(input$u2)

  # Sums of error variances and of loadings over each group factor's assigned items.
  sums_e_s <- colSums(input$u2 * factor_corres)
  sums_s_s <- colSums(omega_mat)

  # The variables a group factor's composite is made of. The omega numerators and
  # denominators are formed from this one index vector, so each coefficient's two halves
  # always describe the same composite. The map is 0/1 by the check above, so these are the
  # same items the masked sums above are taken over.
  members <- lapply(seq_len(ncol(s_load)), function(i) which(factor_corres[, i] == 1))

  # Model-implied common variance of the variables, Lambda Psi Lambda'. The general factor
  # is an extra column orthogonal to the group factors, so the quadratic form splits into a
  # general and a group part; Psi is the identity for a bifactor or Schmid-Leiman solution
  # and the factor intercorrelations Phi for a correlated-factors one, whose general column
  # is zero throughout (see above). Summing the block belonging to a set of variables gives
  # that composite's model-implied true score variance, 1' Lambda Psi Lambda' 1, which is the
  # numerator of every omega total below (Zinbarg, Revelle, Yovel & Li, 2005, Eq. 8; Zinbarg,
  # Yovel, Revelle & McDonald, 2006, Eq. 6). Its diagonal holds the communalities.
  s_mat <- as.matrix(input[, seq_len(ncol(s_load)) + 1L, drop = FALSE])
  Phi <- spec$Phi
  group_common <- if (is.null(Phi)) tcrossprod(s_mat) else s_mat %*% Phi %*% t(s_mat)
  common <- tcrossprod(input$g) + group_common

  # The coefficients below describe the unit-weighted sum of the variables as supplied, so
  # flag a solution whose variables are not all keyed in the same direction before any of
  # them is reported.
  .rel_check_keying(common, var_names)

  if(variance == "correlation"){

    # Compute omega total, hierarchical, and subscale for g-factor
    omega_tot_g <- sum(common) / sum(cormat)
    omega_h_g <- sum_g^2 / sum(cormat)
    omega_sub_g <- sum(sums_s_s^2) / sum(cormat)

    # Compute omega total, hierarchical, and subscale for group factors
    omega_tot_sub <- NULL
    omega_h_sub <- NULL
    omega_sub_sub <- NULL

    # A composite's omega total is the true score variance the model attributes to it over
    # its observed variance. Every factor the composite's variables load on contributes:
    # with correlated factors through the cross-loadings and the factor correlations, and in
    # an estimated bifactor or Schmid-Leiman solution through the cross-loadings alone,
    # which such a solution always has -- only a confirmatory bifactor model carries
    # structural zeros off a composite's own column. Counting the general factor and the
    # composite's own column only would report that factor's congeneric omega, which the
    # omega subscale column already reports, under the name of a total. Reading the
    # composite's block of `common` counts every factor, whichever kind of solution the spec
    # holds.
    for (i in seq_len(ncol(s_load))) {
      subf <- members[[i]]
      Vgr <- sum(cormat[subf, subf])
      omega_tot_sub[i] <- sum(common[subf, subf]) / Vgr
      omega_h_sub[i] <- sum(input$g[subf])^2 / Vgr
      omega_sub_sub[i] <- sum(s_mat[subf, i])^2 / Vgr
    }

  } else if(variance == "sums_load") {

    # Sums of all group factor loadings for all group factors
    sums_s <- colSums(s_mat)

    # What the group factors contribute to the whole-scale composite, 1' S Psi S' 1. For
    # uncorrelated group factors that is the sum of the squared loading-column sums; for a
    # correlated-factors spec the factor correlations enter it as the quadratic form
    # (S'1)' Phi (S'1), which the squared column sums alone would drop. The subscale rows
    # below read the Phi-aware `common`, so scoring this row from the column sums would
    # describe an orthogonal model and a correlated one in the same table. The
    # uncorrelated branch keeps the literal expression, so nothing about the arithmetic of
    # a bifactor or Schmid-Leiman solution changes.
    group_var <- if (is.null(Phi)) sum(sums_s^2) else sum(sums_s * (Phi %*% sums_s))

    # Compute omega total, hierarchical, and subscale for the whole scale. The
    # composite here is all variables, so every loading on every group factor
    # contributes to its variance: the general and group terms use the full
    # loading-column sums, which add up to the whole of the model-implied common
    # variance. omega total is then 1 - sum(u^2) / V, with V
    # the model-implied composite variance, and the general and group variances
    # partition it exactly (tot = hier + sub for the g row).
    omega_tot_g <- (sum_g^2 + group_var) / (sum_g^2 + group_var + sum_e)
    omega_h_g <- sum_g^2 / (sum_g^2 + group_var + sum_e)
    omega_sub_g <- group_var / (sum_g^2 + group_var + sum_e)

    # Compute omega total, hierarchical, and subscale for group factors. A subscale
    # composite's model-implied variance is the common variance it receives from every
    # factor plus its unique variances, so the denominator reads that composite's block
    # of `common` exactly as the whole-scale row above reads the whole matrix. Counting
    # only the general and own-column terms would leave the variance the composite
    # receives through cross-loadings on the other group factors out of both the
    # numerator and the denominator, which understates omega total and inflates omega
    # hierarchical and omega subscale. The three add up on the whole-scale row but not
    # here, where omega subscale stays the factor's own contribution while omega total
    # counts them all.
    omega_tot_sub <- NULL
    omega_h_sub <- NULL
    omega_sub_sub <- NULL

    for (i in seq_len(ncol(s_load))) {
      subf <- members[[i]]
      common_sub <- sum(common[subf, subf])
      Vsub <- common_sub + sum(u2[subf])
      omega_tot_sub[i] <- common_sub / Vsub
      omega_h_sub[i] <- sum(input$g[subf])^2 / Vsub
      omega_sub_sub[i] <- sum(s_mat[subf, i])^2 / Vsub
    }
  }

  # A group factor with no assigned items (an all-zero map column, e.g. a
  # type = "psych" factor that is never the largest loading for any item) has a
  # zero subscale-variance denominator and hence undefined coefficients. Report
  # them as NA -- consistent with the H index below -- instead of a silent NaN.
  empty_fac <- which(colSums(factor_corres == 1) == 0)
  if (length(empty_fac) > 0) {
    omega_tot_sub[empty_fac] <- NA
    omega_h_sub[empty_fac] <- NA
    omega_sub_sub[empty_fac] <- NA
    cli::cli_warn(
      c("Some group factors have no assigned items in {.arg {arg}}.",
        "i" = "Their omega coefficients and H index are returned as {.val NA}."),
      class = "efa_omega_empty_factor"
    )
  }

  # Combine and display results in a table
  omega_tot <- c(omega_tot_g, omega_tot_sub)
  omega_h <- c(omega_h_g, omega_h_sub)
  omega_sub <- c(omega_sub_g, omega_sub_sub)

  # A share of variance above 1, or a general-factor share above the share of all factors
  # together, is not an admissible reliability whichever total variance it was divided by.
  # The two conventions get there by different routes, so the remedy differs with them.
  #
  # In "correlation" mode the denominator (the correlation matrix) is independent of the
  # loading-based numerators, so a correlation matrix inconsistent with the supplied loadings
  # can push a coefficient above 1. It takes a composite whose loadings over-predict its
  # observed variance by more than the whole of its unique variance, which cannot happen when
  # the loadings were estimated from that same matrix, but can when the two come from
  # different sources.
  #
  # In "sums_load" mode the denominator is that composite's own common variance plus its own
  # unique variances, so an omega total above 1 takes those uniquenesses to sum to something
  # negative -- an improper solution, which the Heywood check below reports from the
  # uniquenesses themselves, except that that check is gated on `add_ind` while this one is
  # not. Omega subscale reaches it by a second route that needs no improper solution at all:
  # it divides a group factor's own squared loading sum by its composite's whole variance,
  # and once the group factors are correlated that sum is not a part of it, because the
  # composite's cross-loadings on the other factors can contribute negative common variance
  # and take the total below what its own factor puts in.
  #
  # Either way, a set of factor correlations that is not positive semi-definite makes the
  # group factors contribute negative common variance, which can put omega hierarchical above
  # omega total. Warn rather than return any of these silently.
  tol <- .Machine$double.eps^0.5
  if (any(c(omega_tot, omega_h, omega_sub) > 1 + tol, na.rm = TRUE) ||
      any(omega_h > omega_tot + tol, na.rm = TRUE)) {

    # Without a general factor omega hierarchical is zero throughout, so the second clause
    # there says an omega total came out negative -- and the front-end drops the hierarchical
    # column for such a solution, so naming it would point at a coefficient the result does
    # not carry.
    lead <- if (isTRUE(all(input$g == 0))) {
      "Some omega coefficients are out of range (above 1, or a negative omega total)."
    } else {
      "Some omega coefficients are out of range (above 1, or omega hierarchical above omega total)."
    }

    # The uniquenesses are named rather than pointed at as {.arg u2}: they are only an
    # argument on the manual and bare-matrix routes, and come from the model on the others.
    hints <- if (variance == "correlation") {
      "Check that {.arg cormat} is consistent with the loadings."
    } else {
      "Check the uniquenesses of the affected composite -- its own subtotal, not their sum over all variables -- for a negative value."
    }
    if (!is.null(Phi)) {
      hints <- c(hints, paste("Omega subscale can also exceed 1 in a proper solution, when a",
                              "composite's cross-loadings on the other group factors reduce its",
                              "variance below what its own factor contributes; inspect those and",
                              "{.field Phi}."))
    }
    names(hints) <- rep("i", length(hints))

    cli::cli_warn(c(lead, hints), class = "efa_omega_out_of_range")
  }

  # Optional reliability / convergent-validity coefficients, appended as columns when
  # requested. Kept behind add_rel so the OMEGA coefficient menu is unchanged otherwise.
  if (isTRUE(add_rel)) {

    # Standardized Cronbach's alpha (Cronbach, 1951), k / (k - 1) *
    # (1 - sum(diag(R_sub)) / sum(R_sub)), for the whole scale (all items) and each
    # subscale (the map's assigned items). It needs a correlation matrix: the supplied
    # cormat, or -- when none was given -- the model-implied common variance plus the
    # uniquenesses (which assumes the model holds), standardized to a unit diagonal with
    # cov2cor so the formula returns the standardized coefficient even when the loadings
    # and uniquenesses do not complete to unit item variance. Standardized alpha uses the
    # correlation matrix; raw, covariance-based alpha would need item standard deviations
    # the spec does not carry. Alpha is undefined for fewer than two items (returned NA).
    R_rel <- if (!is.null(cormat)) {
      cormat
    } else {
      stats::cov2cor(common + diag(u2, nrow = length(u2)))
    }

    composites <- c(list(seq_len(nrow(s_mat))), members)
    alpha <- vapply(composites, function(idx) {
      k <- length(idx)
      if (k < 2L) return(NA_real_)
      Rk <- R_rel[idx, idx, drop = FALSE]
      k / (k - 1) * (1 - sum(diag(Rk)) / sum(Rk))
    }, numeric(1))

    # Composite reliability (congeneric omega; Joreskog, 1971; Raykov, 2001) and average
    # variance extracted (AVE; Fornell & Larcker, 1981) per group factor, over that
    # factor's assigned composite (simple structure): CR = (sum L)^2 / ((sum L)^2 +
    # sum u2), AVE = sum L^2 / (sum L^2 + sum u2), summed over the assigned items only --
    # a foreign item's cross-loading never enters the subscale score, so it must not enter
    # its reliability. AVE is a convergent-validity index, not a reliability. Both describe
    # a single group factor, so the general-factor (whole-scale) row is NA. sums_s_s and
    # sums_e_s are the assigned loading and uniqueness sums already formed for the omegas,
    # and omega_mat is the map-masked group loadings, so colSums(omega_mat^2) is the
    # assigned sum of squared loadings.
    sq_assigned <- colSums(omega_mat^2)
    CR  <- c(NA_real_, sums_s_s^2 / (sums_s_s^2 + sums_e_s))
    AVE <- c(NA_real_, sq_assigned / (sq_assigned + sums_e_s))

    # A group factor with no assigned items has no composite: NA, as for the omegas.
    # CR and AVE would otherwise be 0 / 0 = NaN; alpha is already NA for an empty
    # composite via its k < 2 guard, so only CR and AVE need the explicit override.
    if (length(empty_fac) > 0) {
      CR[empty_fac + 1] <- NA
      AVE[empty_fac + 1] <- NA
    }
  }

  if(isTRUE(add_ind)){

    # Compute H index, ECV, and PUC
    h_s <- vector("double", ncol(s_load))
    h_s_load <- vector("double", ncol(s_load))
    hiload_s <- logical(ncol(s_load))

    # The H index (Hancock & Mueller, 2001) of a group factor is defined over that factor's
    # own indicators, so its loadings are restricted to the items assigned to the factor by
    # the map; .h_index guards its denominator against a loading of absolute value at least 1.
    # `hiload_s` flags, per group factor, such a loading on an item with a positive unique
    # variance -- i.e. an oblique pattern coefficient >= 1 in a proper solution, not a Heywood
    # case -- so its undefined H is reported separately from (and independently of) any Heywood
    # case below. ECV (Rodriguez et al., 2016) is a ratio of common variances and uses the
    # full group-factor loading columns. The map defines item-to-factor membership, which also
    # enters the PUC.
    for (j in seq_len(ncol(s_load))){
      s_j <- input[[j + 1]]
      mem <- factor_corres[, j] == 1
      s_mem <- s_j[mem]
      h_s[j] <- .h_index(s_mem)
      hiload_s[j] <- length(s_mem) > 0 &&
        any(abs(s_mem) >= 1 & u2[mem] >= .Machine$double.eps, na.rm = TRUE)
      h_s_load[j] <- sum(s_j^2)
    }

    h_g <- .h_index(input$g)
    hiload_g <- any(abs(input$g) >= 1 & u2 >= .Machine$double.eps, na.rm = TRUE)

    # Heywood case: an improper solution, i.e. a unique variance at or below zero
    # (equivalently a communality at or above 1). This is the standard definition, as in
    # EFA() and the lavaan adapter, and is a property of the solution rather than of the H
    # index, so it is flagged from the uniquenesses (not the loadings). Coefficients
    # involving the affected variables may not be interpretable. which() drops any NA
    # uniqueness (a malformed component, not a Heywood case) so it neither triggers a
    # false Heywood warning nor suppresses the h-undefined note below.
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

    # Separately, the H index is undefined for any construct with a loading whose absolute
    # value is at least 1 (its 1 - lambda^2 error term is non-positive), so .h_index has
    # returned NA for it. For standardized (orthogonal) loadings a loading >= 1 only occurs
    # in a Heywood case, already flagged above from the uniquenesses; hiload_* count only the
    # constructs whose offending loading has a positive unique variance (an oblique pattern
    # coefficient in a proper solution), so this note is per-construct -- it fires for such a
    # construct even when an unrelated Heywood case exists elsewhere, and stays silent when
    # the loading >= 1 is itself the Heywood case.
    if (isTRUE(hiload_g) || any(hiload_s)) {
      cli::cli_warn(
        c("The H index could not be computed for some constructs (a loading with an absolute value of at least 1).",
          "i" = "The H index is returned as {.val NA} for those constructs."),
        class = "efa_reliability_h_undefined"
      )
    }

    ECV <- sum(input$g^2) / sum(sum(input$g^2), sum(h_s_load))

    # Proportion of uncontaminated correlations from the map membership (matching
    # the H index; the general factor is not a contaminant; see .puc).
    PUC <- .puc(factor_corres == 1)

    # Create output
    h <- c(h_g, h_s)

    omegas <- cbind(omega_tot, omega_h, omega_sub, h, NA, NA)
    omegas[1, 5] <- ECV
    omegas[1, 6] <- PUC
    colnames(omegas) <- c("tot", "hier", "sub", "H", "ECV", "PUC")

  } else {

    omegas <- cbind(omega_tot, omega_h, omega_sub)
    colnames(omegas) <- c("tot", "hier", "sub")

  }

  # Append the reliability / convergent-validity columns when computed.
  if (isTRUE(add_rel)) {
    omegas <- cbind(omegas, alpha = alpha, CR = CR, AVE = AVE)
  }

  rownames(omegas) <- c("g", spec$fac_names)

  class(omegas) <- "OMEGA"

  omegas

}

# Front-end adapters: normalize each reliability input source to the spec consumed
# by .reliability_core(). Every adapter returns the same list
# (g_load, s_load, u2, map, cormat, var_names, fac_names), so the core stays
# input-agnostic and the front-end only has to pick the right adapter.

# Resolve the item-to-factor correspondence map for adapters that take a separate
# `factor_corres`/`type` (SL, schmid, manual, efa): with a supplied `factor_corres`,
# validate its dimensions and use it; otherwise, for `type = "psych"`, assign each
# item to its highest-|loading| group factor (as psych::omega does); `type =
# "EFAtools"` requires an explicit `factor_corres`. `arg` names the user-facing
# argument the map came from, as in .rel_check_map().
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

  # type == "psych": each variable's correspondence is its largest |group loading|.
  map <- matrix(0, nrow = nrow(s_load), ncol = ncol(s_load))
  for (i in seq_len(nrow(s_load))) {
    map[i, which.max(abs(s_load[i, ]))] <- 1
  }
  map

}

# Abort when a user-supplied item-to-factor map does not conform to the group loadings
# it will be applied to: one row per variable and one column per group factor. Without
# this check a map with too few rows is recycled against the loadings and returns
# coefficients above 1 rather than an error. `arg` names the user-facing argument the map
# came from, as in .rel_check_map().
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

# Sanity-check a user-supplied item-to-factor map against the loadings it will be
# applied to. The map's columns are matched to the group factors by position, so a map
# written in a different factor order than the solution yields well-formed but
# meaningless subscale coefficients rather than an error. Flag each column whose
# assigned items barely load on it (mean |loading| below `min_load`) while the same
# items do load on some other column -- the signature of a permuted or mistranscribed
# map. Requiring a better home for the items keeps the check to that signature: a map
# whose items load on nothing is a weak solution, not a misaligned map, and is left
# alone. `min_load` is deliberately far below the salience conventions used elsewhere
# in the package (`salience_threshold`, .3 by default): a deliberate map may well
# disagree with the loadings, so only assignments that are essentially unsupported are
# flagged. Silent by construction for a map derived by .rel_map(): there each item is
# assigned to its own largest |loading|, so a column's assigned mean can never fall
# below another column's mean for those items. `arg` names the user-facing argument the
# map came from, so the message points at the right one for each front-end. A map that
# does not conform to the loadings is left to the caller, which validates the dimensions.
.rel_check_map <- function(s_load, map, min_load = 0.1, arg = "factor_map") {

  s_load <- abs(as.matrix(s_load))
  map <- as.matrix(map)
  if (!identical(dim(map), dim(s_load))) return(invisible(NULL))
  assigned <- map == 1

  labs <- colnames(s_load)
  if (is.null(labs)) labs <- as.character(seq_len(ncol(s_load)))

  implausible <- vapply(seq_len(ncol(assigned)), function(j) {
    idx <- which(assigned[, j])
    # An empty column carries no items to judge; the core reports it separately.
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

# Warn when the variables of a solution are not all keyed in the same direction -- most
# often a scale whose reverse-worded items were never reverse-coded. Every coefficient
# computed here describes the raw unit-weighted sum of the variables as supplied, so a
# variable keyed against the rest subtracts from that sum's true score variance instead of
# adding to it and the coefficients collapse: correct for the sum that was scored, but
# reading as a poor scale rather than as a missing step. The variables are deliberately not
# reflected first, as psych::omega does, because reflection describes a different composite
# than the one the coefficients are defined on (Flora, 2020).
#
# `common` is the model-implied common variance Lambda Psi Lambda' of the variables, whose
# total sum(common) = 1' Lambda Psi Lambda' 1 is the composite's true score variance. Under
# any re-keying D = diag(+-1) of the variables that variance becomes sum_ij d_i d_j
# common_ij, which is at most sum(abs(common)); the ratio of the two is therefore 1 exactly
# when no variable works against the sum and falls towards 0 as the keying becomes balanced.
# The rule reads the common variance rather than the loading signs on purpose: `g_load` is
# zero throughout for a correlated-factors spec, so a sign test on it never fires there, and
# the group loadings carry no sign information either, because the factor columns are
# reflected -- a reverse-keyed block comes back with positive loadings and a negative factor
# correlation, which the common variance shows and the loadings do not.
#
# `cutoff` is deliberately far below any well-keyed solution, as `min_load` is in
# .rel_check_map(): the ratio is (P - N) / (P + N) over the positive and negative common
# variance, so firing takes negative common variance exceeding about a seventh of the
# positive. On the bundled data no correctly keyed solution reaches that: the lowest is .89
# (DOSPERT_raw at two factors, whose domains are only weakly related), whereas a scale whose
# reverse-worded items were never reverse-coded stays below .32 (UPPS_raw, two to six
# factors) and reverse-coding one factor's block of a well-keyed solution takes it to .23.
#
# The variables named are those whose own common variance with the composite, the row sum of
# `common`, is negative -- the ones subtracting from the sum, and so the ones to reverse-code.
# They are named for guidance and not as the trigger: a variable can drag the total down
# without its own row sum turning negative, so the ratio decides whether to warn and the row
# sums only say where to look. The list is therefore allowed to come out empty, in which case
# the message states the finding without it. `var_names` may be NULL, as for the Heywood
# report above, in which case the positions stand in for the names.
.rel_check_keying <- function(common, var_names = NULL, cutoff = 0.75) {

  aligned <- sum(abs(common))

  # Nothing to judge: no loadings at all, or a malformed component (NA) the callers
  # report separately.
  if (!is.finite(aligned) || aligned <= .Machine$double.eps^0.5) return(invisible(NULL))

  if (isTRUE(sum(common) < cutoff * aligned)) {

    against <- which(rowSums(common) < 0)
    n_bad <- length(against)

    lead <- if (n_bad > 0) {
      labs <- if (is.null(var_names)) as.character(against) else var_names[against]
      # Capped as elsewhere: a long scale can have dozens of reverse-worded items, and a
      # message that prints them all is unreadable.
      cap <- .cap_label_list(labs)
      paste("The variables are not all keyed in the same direction:",
            "{.val {cap$shown}}{cap$rest} {cli::qty(n_bad)}{?is/are} keyed against the rest.")
    } else {
      "The variables are not all keyed in the same direction."
    }

    cli::cli_warn(
      c(lead,
        "i" = "The coefficients describe the raw unit-weighted sum of the variables as supplied, in which the reversed variables cancel rather than add.",
        "i" = "Reverse-code them before fitting the solution (Flora, 2020)."),
      class = "efa_reliability_mixed_keying"
    )
  }

  invisible(NULL)

}

# Warn (classed) when supplied uniquenesses do not complete the supplied loadings to unit
# variance. `h2` is the diagonal of the model-implied common variance Lambda Psi Lambda',
# so a standardized solution has h2 + u2 = 1 for every variable. Only the routes whose
# uniquenesses come from the user need this: every fitted route derives them as 1 - h2, and
# so satisfies it by construction.
#
# What it catches is a slip on adjacent columns. The printed Schmid-Leiman table ends in
# `h2` and `u2`, and many published tables report only `h2`, leaving the reader to form
# `1 - h2` by hand -- so passing communalities where uniquenesses belong is a one-character
# mistake rather than a contrived one. It is silent otherwise: `u2` does not enter the
# coefficients at all under `variance = "correlation"` with a correlation matrix, and that
# is the mode a user without one is turned away from, so the mistake is inert exactly where
# it would be harmless and load-bearing exactly where it is not.
#
# A warning rather than an error, because unstandardized (covariance-metric) components are
# a defensible input that yields a legitimate raw-score coefficient; the user needs to know
# the components are not standardized, not to be stopped. `tol` sits well clear of both
# sides: a published table rounded to two decimals leaves a residual of about .015, while
# supplying communalities leaves at least .40.
.rel_check_u2 <- function(h2, u2, tol = 0.05) {

  # Nothing to compare: no variables at all, or a `u2` that does not line up with the
  # loadings, which is the callers' error to report and not a statement about
  # standardization. Either would otherwise reach `max()` over an empty vector or compare a
  # recycled one.
  if (length(h2) == 0L || length(h2) != length(u2)) return(invisible(NULL))

  resid <- max(abs(h2 + u2 - 1))

  # A non-finite residual is a malformed component, which the callers report separately.
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

# Abort (classed) unless `Phi` can be the factor intercorrelation matrix of the `k` group
# factors it will be scored with. `NULL` means the identity -- uncorrelated group factors, as a
# Schmid-Leiman or bifactor solution has by construction -- and is always admissible, so this
# only runs where a matrix was actually supplied.
#
# The checks matter because `Phi` enters the model-implied common variance as the quadratic form
# S Phi S', where a matrix that is not a correlation matrix does not merely mislabel the solution:
# an asymmetric one is silently symmetrised by the quadratic form, a non-unit diagonal rescales
# the factors, and one that is not positive semi-definite attributes negative common variance to a
# composite and returns coefficients outside [0, 1] with nothing else to catch them.
#
# The shape tests are written out rather than delegated to .is_cormat(): that helper classifies
# raw data against a correlation matrix, so it rejects by design the 1 x 1 matrix a single group
# factor takes, and it aborts on a missing value instead of reporting one. Positive
# semi-definiteness is judged on the smallest eigenvalue against the package's usual tolerance --
# an estimated factor correlation matrix is often singular to rounding, which is admissible, while
# a genuinely indefinite one is not.
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
  # Symmetric by the test above, so the eigenvalues are real.
  if (min(eigen(Phi, symmetric = TRUE, only.values = TRUE)$values) < -tol) {
    bad("It is not positive semi-definite, so it is not a correlation matrix.")
  }

  invisible(NULL)

}

# Abort (classed) when a supplied correlation matrix fails the .is_cormat() check;
# shared guard for the adapters that accept a user cormat. A correlation matrix supplied as a
# data frame -- what read.csv() returns for a published correlation table -- is accepted and
# handed back as a matrix, as .prepare_cor_input() does on the other input route. The coercion
# is load-bearing for one operation only: the standardized-alpha branch of .reliability_core()
# takes diag() of a subset of it, which a data frame does not support. The omega quantities
# only sum and subset it, which a data frame does support -- so without the coercion they
# would return correct numbers and only add_rel = TRUE would fail. Keep the coercion here
# rather than at the one operation, so the spec carries a single type.
#
# The matrix is also checked against the solution it will be used with, given as the
# number of variables (`n_items`) and, where the solution carries them, their names
# (`var_names`). The subscale coefficients divide loading sums by the correlations of the
# same items, so a correlation matrix of other variables, or of the same variables in
# another order, yields well-formed but wrong subscale coefficients while leaving the
# whole-scale row -- which only sums the whole matrix -- untouched, and one of the wrong
# size fails much later as an out-of-bounds subscript. Named variables in a different
# order are an unambiguous permutation and are reordered to the solution; a different set
# of names cannot be resolved and aborts.
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

  # Only labelled variables can be matched. An unlabelled correlation matrix, or an
  # unlabelled solution, is matched by position as before.
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

# Adapter: normalize an SL() object to the reliability spec. Mirrors the SL branch of
# .OMEGA_FLEX. Honours a user cormat, else the SL object's stored correlation matrix
# (orig_R = NA for flexible-input SL objects, in which case cormat stays NULL and
# only variance = "sums_load" applies).
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

  list(g_load = g_load, s_load = s_load, u2 = u2,
       map = .rel_map(s_load, factor_corres, type),
       cormat = cormat, var_names = var_names, fac_names = fac_names)

}

# Adapter: normalize a psych::schmid() object to the reliability spec. Mirrors the
# schmid branch of .OMEGA_FLEX. When no cormat is supplied, it is reconstructed from
# the schmid pattern matrix and factor intercorrelations via psych::factor.model.
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

# Adapter: normalize manually supplied SL/bifactor components (g_load, s_load, u2) to
# the reliability spec. Mirrors the model = NULL branch of .OMEGA_FLEX. A cormat may be
# given directly or reconstructed from an oblique pattern matrix and factor
# intercorrelations -- one or the other, not both -- or left NULL, which only
# variance = "sums_load" can then score.
#
# `Phi` describes the loadings it is paired with. Supplied together with `pattern` it belongs
# to that separate oblique solution and is used only to rebuild the correlation matrix, as it
# always has been -- the pairing is what a Schmid-Leiman input needs, where `s_load` holds the
# orthogonalized group loadings and `pattern` the parent oblique ones, two different matrices
# whose model-implied correlations differ once more than three factors make the second-order
# model overidentified. Supplied without a `pattern` there is no other loading matrix for it
# to describe, so it is the factor intercorrelation matrix of `s_load` itself: the components
# then specify a correlated-factors solution, and it enters the model-implied common variance
# as the quadratic form S Phi S'. That is the one way this route can express correlated group
# factors; without it the coefficients would describe an orthogonal model instead, and fall
# short of the ones the same solution gets through `efa_fit()`.
.rel_adapt_manual <- function(g_load, s_load, u2, var_names, factor_corres = NULL,
                              type = "EFAtools", cormat = NULL, pattern = NULL,
                              Phi = NULL, fac_names = NULL) {

  # unclass() before as.matrix(): a Schmid-Leiman loading table is a matrix already, so
  # as.matrix() returns it with its class intact, and the data.frame() the core builds from
  # it then keeps it as a single column instead of expanding it -- which fails only later,
  # on a names/length mismatch that says nothing about the input.
  s_load <- as.matrix(unclass(s_load))

  # `cormat` and `pattern` are two ways to give one thing, so taking both is an ambiguous
  # request rather than a redundant one. It also cannot be resolved quietly: with a `cormat`
  # in hand the reconstruction never runs and `pattern` is read by nothing, yet its presence
  # still decides whether `Phi` describes this solution or that one -- an argument with no
  # effect of its own silently choosing between two different sets of coefficients. Refuse the
  # combination so the caller says which correlation matrix they mean.
  if (!is.null(cormat) && !is.null(pattern)) {
    cli::cli_abort(
      c("{.arg cormat} and {.arg pattern} are two ways to give the same correlation matrix, so supply one or the other.",
        "i" = "Drop {.arg pattern} to score the components against {.arg cormat}.",
        "i" = "Drop {.arg cormat} to reconstruct one instead, from {.arg pattern} and the {.arg Phi} of the oblique solution it came from."),
      class = "efa_reliability_cormat_and_pattern"
    )
  }

  # Whether `Phi` describes this solution's own group factors, rather than a `pattern`'s.
  phi_solution <- !is.null(Phi) && is.null(pattern)

  # Checked on either reading: paired with a `pattern` it is what the reconstructed
  # correlation matrix is built from, and so ends up in every coefficient's denominator,
  # which makes it exactly as load-bearing as it is on its own. It has to be the correlation
  # matrix of whichever loadings it describes, hence the two column counts.
  if (!is.null(Phi)) {
    .rel_check_phi(Phi, if (phi_solution) ncol(s_load) else NCOL(pattern))
  }

  if (phi_solution) {
    # A general factor and correlated group factors together do not give the variance
    # decomposition these coefficients report: omega hierarchical and omega subscale would
    # no longer partition the composite, and the PUC presupposes that two variables of
    # different group factors share only the general factor. The core makes the same
    # statement in its own terms; refuse the combination here, where it can be named.
    if (!isTRUE(all(g_load == 0))) {
      cli::cli_abort(
        c("{.arg Phi} describes correlated group factors, which cannot be combined with a general factor.",
          "i" = "Supply {.arg Phi} for a correlated-factors solution, whose {.arg g_load} is zero throughout.",
          "i" = "For a Schmid-Leiman or bifactor solution the group factors are orthogonal; pass {.arg pattern} as well if {.arg Phi} belongs to the oblique solution it came from."),
        class = "efa_reliability_phi_with_general"
      )
    }
  }

  # A `pattern` alongside `Phi` says the two belong together and `s_load` is a separate,
  # orthogonalized matrix -- which is exactly right for a Schmid-Leiman input, and exactly
  # wrong for components that carry no general factor, where `s_load` is itself the oblique
  # pattern and the factor correlations belong in its common variance. Here `pattern` is
  # load-bearing (it is reconstructing the correlation matrix, or the check above would have
  # aborted), so the reading cannot be settled from the arguments: a genuine parent oblique
  # solution and a `pattern` that merely repeats `s_load` look the same. Take the documented
  # one and say what was assumed, since the two give materially different coefficients.
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
    # Matched against the loadings' own row names, not `var_names`: here `var_names` are
    # the labels the caller wants on the output rows, which need not be the names the
    # correlation matrix carries. Unlabelled components are matched by position, as the
    # manual contract -- every component in the row order of the loadings -- implies.
    cormat <- .rel_check_cormat(cormat, rownames(s_load), nrow(s_load))
  }

  # The communalities the supplied loadings imply -- the diagonal of Lambda Psi Lambda',
  # under the factor correlations where this spec carries them. The uniquenesses are the
  # caller's own here, so they are checked against these.
  group_h2 <- if (phi_solution) rowSums((s_load %*% Phi) * s_load) else rowSums(s_load^2)
  .rel_check_u2(g_load^2 + group_h2, u2)

  if (is.null(fac_names)) fac_names <- seq_len(ncol(s_load))

  list(g_load = g_load, s_load = s_load, u2 = u2,
       map = .rel_map(s_load, factor_corres, type),
       Phi = if (phi_solution) Phi else NULL,
       cormat = cormat, var_names = var_names, fac_names = fac_names)

}

# Adapter: normalize a lavaan single-factor, correlated-factors, bifactor, or
# second-order solution to the reliability spec, one entry per fitted group. Mirrors the
# structural detection and Schmid-Leiman transform of .OMEGA_LAVAAN (second-order group
# loadings via .sl_group_loadings), and drives the core with variance = "sums_load" (the
# composite variances are model-implied). Returns a list with, per group, a full spec
# (correlated factors / bifactor / second-order / single factor); a single-factor group is
# additionally flagged `single = TRUE`, since a solution with one factor defines fewer
# coefficients than the core computes for it. A correlated-factors spec is flagged by
# `correlated = TRUE` on the returned list: it carries no general factor, so a front-end
# that reports the general-factor coefficients has to drop them (efa_reliability) or
# refuse the input (OMEGA).
.rel_adapt_lavaan <- function(model, g_name = "g", group_names = NULL) {

  .require_lavaan()

  if (isFALSE(lavaan::lavInspect(model, what = "converged"))) {
    cli::cli_abort("The model did not converge; no reliability coefficients are computed.",
                   class = "efa_reliability_no_converge")
  }

  std_sol <- suppressWarnings(lavaan::lavInspect(model, what = "std",
                                drop.list.single.group = FALSE))

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

    lambda <- std_sol[[i]][["lambda"]]
    theta <- std_sol[[i]][["theta"]]
    psi <- std_sol[[i]][["psi"]]

    if (any(is.na(lambda))) {
      cli::cli_abort("Some loadings are {.val NA} or {.val NaN}; no reliability coefficients are computed.",
                     class = "efa_reliability_na_loadings")
    }

    if (any(diag(theta) <= 0) || any(diag(psi) <= 0)) {
      cli::cli_abort("A Heywood case was detected (a variance of 0 or negative); no reliability coefficients are computed.",
                     class = "efa_reliability_heywood")
    }

    var_names_i <- rownames(lambda)

    # A single factor has no group factors: normalize it to the spec such a solution takes on
    # every other input route, which the core scores like any of them. Marked `single = TRUE`
    # so the front-ends can report the coefficients it defines.
    if (ncol(lambda) == 1) {
      groups[[i]] <- c(
        # Carried beside the spec rather than in its `fac_names`, which the core reads to
        # build the group-factor rows this spec has none of.
        list(single = TRUE, fac_label = colnames(lambda)),
        .rel_single_factor_spec(list(g_load = lambda[, 1],
                                     s_load = lambda[, 0, drop = FALSE],
                                     u2 = diag(theta), cormat = NULL,
                                     var_names = var_names_i)))
      next
    }

    col_names <- colnames(lambda)

    # Detect the model type once (all groups share the fixed-zero structure).
    if (i == 1) {

      # A second-order model routes the covariances of its first-order factors through
      # `beta`, which none of the other supported structures has. Read from `beta` rather
      # than from an all-zero general-factor column, because that column cannot be found
      # at all when the general factor is misnamed -- and such a fit would then pass the
      # simple-structure test below and be scored as a set of correlated factors, whose
      # `psi` is here the residual covariance matrix of the first-order factors rather
      # than their correlation matrix.
      beta <- std_sol[[i]][["beta"]]
      higherorder <- !is.null(beta) && any(abs(beta) > tol, na.rm = TRUE)

      # A variable loading on two or more factors is what a general factor over and above
      # the group factors looks like in a loading matrix, and is the structure a bifactor
      # model is recognized by. Without one -- and without the second-order structure
      # above -- there is no general factor whose variance a composite could be
      # decomposed into, and the fit is an ordinary set of correlated factors. Decided
      # from the structure rather than from `g_name`, which says nothing on its own: a
      # factor may be named "g" without being a general factor, and a general factor need
      # not be named "g".
      bi_check <- lambda
      bi_check[abs(bi_check) > tol] <- 1
      correlated <- !higherorder && all(rowSums(bi_check) < 2)

      if (isTRUE(higherorder)) {
        if (sum(colSums(beta) > 0) > 1) {
          cli::cli_abort("The higher-order model has more than two latent strata or more than one second-order factor; only second-order models with one second-order factor are supported.",
                         class = "efa_reliability_higher_order")
        }
      } else if (!correlated) {
        # A genuine bifactor has each item loading on the general and a group
        # factor; flag the borderline case where some item loads on fewer than
        # two factors so the front-end can warn the user (it is still scored).
        few_loadings <- !all(rowSums(bi_check) > 1)
      }
    }

    # Both structures with a general factor locate it by name; a correlated-factors
    # solution has none to find, so it needs no name and reads none.
    if (!correlated && !any(col_names %in% g_name)) {
      cli::cli_abort(
        c("Could not find the specified general-factor name in the lavaan solution.",
          "i" = "Please check the spelling.",
          "i" = if (!higherorder) {
            "Some variables load on two or more factors, so the fit is read as a bifactor solution; a correlated-factors solution, which needs no general factor, has each variable loading on one factor only."
          }),
        class = "efa_reliability_g_name"
      )
    }

    if (isTRUE(correlated)) {
      # An ordinary set of correlated factors: no general factor, the loadings are the
      # oblique pattern, and the standardized `psi` is the factor correlation matrix
      # belonging to it. The core reads the pair as the model-implied common variance
      # Lambda Psi Lambda', exactly as it does for an oblique `efa_fit()` solution; the
      # front-end drops the coefficients such a solution does not define.
      # `psi` is indexed by name rather than taken whole, as the second-order branch below
      # does: the core reads it as the correlations of the columns of `s_load` in their
      # order, which the two matrices' shared latent ordering supplies but nothing checks.
      # Unlabelled columns, which a fitted model does not have, keep the positional
      # reading rather than subscripting `psi` down to nothing.
      Phi_i <- if (is.null(col_names)) psi else psi[col_names, col_names, drop = FALSE]
      groups[[i]] <- list(single = FALSE, g_load = rep(0, nrow(lambda)),
                          s_load = lambda, u2 = diag(theta), Phi = Phi_i,
                          map = abs(lambda) > tol, cormat = NULL,
                          var_names = var_names_i, fac_names = col_names)
      next
    }

    # The coefficients decompose a composite's variance into a general part and one part
    # per group factor, which requires the latent variables to be uncorrelated: a
    # bifactor model is fitted with orthogonal factors, and a second-order model routes
    # the first-order covariances through `beta`, leaving `psi` diagonal either way.
    # `lavaan::cfa()` does not impose that by default, so a "bifactor" model fitted
    # without `orthogonal = TRUE` returns correlated factors whose covariances the
    # coefficients below would silently drop, as does a second-order model with a freed
    # residual covariance between first-order factors. Reject either rather than report
    # a decomposition of a model that does not admit one. An NA here is a malformed
    # solution rather than a correlation, and is left to the checks above.
    if (any(abs(psi[upper.tri(psi)]) > tol, na.rm = TRUE)) {
      cli::cli_abort(
        c("The factors of the {.cls lavaan} solution are correlated; the coefficients need uncorrelated factors.",
          "i" = "A bifactor model needs {.code orthogonal = TRUE}; a second-order model needs no freed covariances between the first-order factors."),
        class = "efa_reliability_correlated_factors"
      )
    }

    col_names <- col_names[!col_names %in% g_name]

    if (isTRUE(higherorder)) {
      # Schmid-Leiman the second-order solution: direct general-factor loadings from
      # the first-order loadings times the second-order (beta) loadings, and direct
      # group-factor loadings via .sl_group_loadings. The general column is computed
      # from the original first-order loadings before they are overwritten.
      lambda[, g_name] <- lambda[, col_names] %*% std_sol[[i]][["beta"]][col_names, g_name]
      lambda[, col_names] <- .sl_group_loadings(lambda[, col_names], psi, col_names)
    }

    s_load <- lambda[, col_names, drop = FALSE]

    groups[[i]] <- list(single = FALSE, g_load = lambda[, g_name], s_load = s_load,
                        u2 = diag(theta), map = abs(s_load) > tol, cormat = NULL,
                        var_names = var_names_i, fac_names = col_names)

  }

  list(groups = groups, group_names = group_names, variance = "sums_load",
       higher_order = higherorder, few_loadings = few_loadings,
       correlated = correlated)

}

# Blank the cells of a computed coefficient matrix that a correlated-factors solution
# (one with no general factor) does not define, so the result builder omits them as it
# omits any other undefined coefficient.
#
# The general-factor decomposition is not identified without a general factor: omega
# hierarchical and ECV are structurally zero, and PUC ("percent of uncontaminated
# correlations") presupposes a general factor for the cross-factor correlations to be
# uncontaminated of. On the whole-scale (g) row the omega subscale and H index are
# further artifacts of the synthetic all-zero general-factor column rather than
# coefficients (the H of an all-zero loading vector is 0, and the whole-scale subscale
# omega does not partition the composite without a general factor). What remains is what
# such a solution does define: whole-scale omega total and alpha, and each group factor's
# congeneric omega, H, and alpha, which stay on their own rows.
.rel_drop_general <- function(x) {
  x[, c("hier", "ECV", "PUC")] <- NA
  x["g", c("sub", "H")] <- NA
  x
}

# The canonical spec of a solution with exactly one factor, or NULL when `spec` does not
# describe one.
#
# One factor is one factor whichever slot it arrived in. An input that names it the general
# factor (a single-factor `lavaan` fit, a one-column bifactor matrix) and one that names it
# the only group factor (a one-factor `efa_fit()` solution, manual components with a zero
# `g_load`) are the same model, so both are rewritten the same way -- the loadings as
# `g_load`, no group-factor columns -- and every route reaches the same coefficients rather
# than one reading of the same solution per input format. Any `Phi` goes with them: the
# correlation matrix of a single factor can only be the 1 x 1 identity, which .rel_check_phi()
# has already established and which says nothing about the solution.
#
# `fac_names` is emptied here only because the core builds the row labels as
# c("g", fac_names) and this spec has one row; the name the input gave the factor is applied
# to that row afterwards, by .rel_drop_single_factor(), from .rel_single_factor_label().
#
# The core scores the result as it scores every other spec -- its `seq_len(k) + 1` indexing
# and its empty-map guards already cover a spec with no group factors -- so a single factor
# needs no arithmetic of its own. What such a solution does not define is dropped afterwards,
# by .rel_drop_single_factor().
.rel_single_factor_spec <- function(spec) {

  s_load <- as.matrix(spec$s_load)
  n_g <- if (isTRUE(all(spec$g_load == 0))) 0L else 1L
  if (n_g + ncol(s_load) != 1L) return(NULL)

  g_load <- if (n_g == 1L) spec$g_load else s_load[, 1]
  # No group factors, in the p x 0 shape the core reads as none.
  none <- matrix(numeric(0), nrow = length(g_load), ncol = 0)

  list(g_load = g_load, s_load = none, u2 = spec$u2, map = none, Phi = NULL,
       cormat = spec$cormat, var_names = spec$var_names, fac_names = character(0))

}

# The row label of a single factor: the name the input gives it, or "F1" where the input gives
# none. A factor supplied as the only group factor is named in `fac_names` -- from the
# solution's loading columns, or from the user's own argument -- and keeps that name, as does
# a `lavaan` factor, which the model syntax always names.
#
# The fallback is the default first-factor label rather than the general-factor "g" that the
# multi-factor solutions use on this row. A single factor is the whole model, not the general
# factor of a bifactor or hierarchical one, which is what "g" would state; a factor supplied
# as the general factor of a one-column matrix carries no name saying otherwise, so it takes
# the neutral label too. That also leaves every route agreeing on the label, and not only on
# the coefficients, wherever the input names the factor nothing.
#
# A name is a name: an empty or missing one is not, and neither is the bare column position
# the manual route falls back to when the caller names no factors (`seq_len(ncol(s_load))`,
# which is 1 for a single factor).
.rel_single_factor_label <- function(fac_names) {
  if (length(fac_names) != 1L || is.na(fac_names) || is.numeric(fac_names) ||
      !nzchar(as.character(fac_names))) {
    return("F1")
  }
  as.character(fac_names)
}

# Blank the cells of a computed coefficient matrix that a single-factor solution does not
# define, so the result builder omits them as it omits any other undefined coefficient, and
# label its one row with the name the factor came in under (see .rel_single_factor_label).
#
# What one factor leaves is omega total, standardized alpha, and the H index. Alpha assumes
# essentially tau-equivalent items, which is nested in a one-factor model, so this is the one
# solution for which reporting it is defensible rather than merely computable. The others are
# not withheld for want of an estimate but because each would state something the solution
# cannot: omega subscale is the variance due to the group factors, of which there are none;
# omega hierarchical is the same quantity as omega total here, the one factor accounting for
# all of the common variance, and two columns holding the same number invite a reader to
# compare a number with itself; and the ECV and the PUC are 1 by construction, which reads as
# evidence of unidimensionality rather than as the arithmetic of a model with one factor.
.rel_drop_single_factor <- function(x, fac_names = NULL) {
  x[, c("hier", "sub", "ECV", "PUC")] <- NA
  rownames(x) <- .rel_single_factor_label(fac_names)
  x
}

# The single-factor note, stated once for both efa_reliability() paths that reach it (the
# lavaan front-end and the shared spec route), so the two cannot drift apart.
.rel_inform_single_factor <- function() {
  cli::cli_inform(
    c("i" = "The solution has a single factor; omega total, alpha, and the H index are returned.",
      "i" = "Omega subscale needs group factors, omega hierarchical is the same quantity as omega total here, and the ECV and the PUC are 1 by construction, so they are omitted."),
    class = "efa_reliability_single_factor"
  )
}

# Adapter: normalize an oblique EFA() object to the reliability spec, treating it as
# the correlated-factors model it is. A correlated-factors solution identifies
# whole-scale omega total and each factor's congeneric omega/H, but not the
# hierarchical/bifactor indices (omega hierarchical, ECV, PUC), which require a
# general factor. Rather than manufacture one by a Schmid-Leiman transformation --
# which is underidentified for fewer than three factors, imposes proportionality
# constraints, and biases the general loadings under the cross-loadings EFA solutions
# almost always have (Flora, 2020, Adv. Methods Pract. Psychol. Sci.; Mansolf &
# Reise, 2016, Multivariate Behav. Res.) -- the spec carries a zero general factor,
# the oblique pattern as the group loadings, and Phi. The core then returns the
# whole-scale omega total as 1' L Phi L' 1 / 1'R1, each factor's omega total as the
# Phi-aware common variance of its composite over that composite's observed variance,
# each factor's congeneric omega subscale and H, and omega hierarchical / ECV / PUC 0.
#
# A one-factor solution is the exception to the oblique requirement below, having nothing to
# rotate: `efa_fit()` returns it unrotated, and where a rotation was asked for anyway leaves
# `rot_loadings` a copy of the unrotated loadings and no factor intercorrelations. It is read
# from `unrot_loadings`, which such an object always carries, under the 1 x 1 identity its
# single factor takes -- so the same expressions serve it -- and the front-end then rewrites
# the one-column spec into the single-factor one (.rel_single_factor_spec), as it does for
# every other route that can carry one factor. Refusing it for want of an oblique rotation
# would be the one piece of advice a one-factor solution cannot follow.
.rel_adapt_efa <- function(model, factor_corres = NULL, type = "psych",
                           cormat = NULL, fac_names = NULL) {

  L_unrot <- unclass(model$unrot_loadings)
  single <- NCOL(L_unrot) == 1L

  if (!single && !("Phi" %in% names(model))) {
    cli::cli_abort(
      c("{.arg model} is not an oblique EFA solution.",
        "i" = "Reliability from an EFA needs correlated factors; refit with an oblique rotation."),
      class = "efa_reliability_not_oblique"
    )
  }

  L1 <- if (single) L_unrot else unclass(model$rot_loadings)
  Phi <- if (single) diag(1) else model$Phi

  # Order the factor columns by number, as SL() does, for stable labels. Only some
  # rotations label their loading columns, so the order falls back to the columns'
  # own order when the labels are absent or carry no factor number.
  n_order <- .sl_factor_order(colnames(L1), ncol(L1))
  s_load <- L1[, n_order, drop = FALSE]
  Phi <- Phi[n_order, n_order, drop = FALSE]

  # Model-implied common variance L Phi L'; its diagonal gives the communalities.
  common <- s_load %*% Phi %*% t(s_load)
  u2 <- 1 - diag(common)

  if (is.null(cormat)) {
    cormat <- if (is.matrix(model$orig_R)) model$orig_R
              else common + diag(u2, nrow = length(u2))
  } else {
    cormat <- .rel_check_cormat(cormat, rownames(s_load), nrow(s_load))
  }

  if (is.null(fac_names)) {
    fac_names <- colnames(s_load)
    # Unlabelled loading columns still need one row label per group factor.
    if (is.null(fac_names)) fac_names <- paste0("F", seq_len(ncol(s_load)))
  }

  list(g_load = rep(0, nrow(s_load)), s_load = s_load, u2 = u2,
       map = .rel_map(s_load, factor_corres, type), Phi = Phi,
       cormat = cormat, var_names = rownames(s_load), fac_names = fac_names)

}

# Adapter: normalize a raw bifactor loading matrix (general factor in the first
# column, group factors in the rest) to the reliability spec. Uniquenesses default to
# the orthogonal-bifactor communalities (1 - rowSums(L^2)); the correspondence map
# defaults to the nonzero group-loading pattern (overridable via factor_corres); and,
# when no observed cormat is supplied, the model-implied L L' + diag(u2) is used so
# both variance conventions apply.
.rel_adapt_bifactor <- function(loadings, factor_corres = NULL, u2 = NULL,
                                cormat = NULL, fac_names = NULL) {

  loadings <- as.matrix(loadings)
  g_load <- loadings[, 1]
  s_load <- loadings[, -1, drop = FALSE]

  # Derived uniquenesses complete the loadings to unit variance by construction; only ones
  # the caller supplied are worth checking against them.
  if (is.null(u2)) {
    u2 <- 1 - rowSums(loadings^2)
  } else {
    .rel_check_u2(rowSums(loadings^2), u2)
  }

  if (is.null(factor_corres)) {
    factor_corres <- abs(s_load) > .Machine$double.eps * 100
  } else {
    # This adapter does not route through .rel_map(), so a supplied map gets the same
    # dimension check and plausibility check here. Only efa_reliability() reaches this
    # adapter, so the map is always its `factor_map`.
    .rel_assert_map_dim(factor_corres, s_load, arg = "factor_map")
    .rel_check_map(s_load, factor_corres, arg = "factor_map")
  }

  if (is.null(cormat)) {
    cormat <- loadings %*% t(loadings) + diag(u2, nrow = length(u2))
  } else {
    # The loadings' own row names, so an unlabelled matrix is matched by position
    # rather than against the V1, V2, ... fallback labels below.
    cormat <- .rel_check_cormat(cormat, rownames(loadings), nrow(loadings))
  }

  var_names <- rownames(loadings)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(nrow(loadings)))

  if (is.null(fac_names)) {
    # A matrix with no group-factor columns holds a single factor, and the only place the
    # matrix names it is the general column. Take the name from there, so the label the
    # front-end puts on that solution's one row is the one the input gave it -- the same
    # label the same solution gets through the components.
    fac_names <- if (ncol(s_load) == 0L) colnames(loadings) else colnames(s_load)
    if (is.null(fac_names)) fac_names <- seq_len(ncol(s_load))
  }

  list(g_load = g_load, s_load = s_load, u2 = u2, map = factor_corres,
       cormat = cormat, var_names = var_names, fac_names = fac_names)

}

# Canonical map from a `.reliability_core()` output column to its public identity in an
# `efa_reliability` result: the long-format coefficient name, its kind (a reliability
# coefficient versus a common-variance / bifactor index), and the compact label used in
# the printed tables. Columns the core can also produce but that a reliability result does
# not surface (CR, AVE) are intentionally absent, so they never enter the output. The row
# order fixes the coefficient order in the tidy result and the column order in the print.
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

# Shape one or more computed coefficient matrices into the tidy `efa_reliability` result:
# a long data.frame `{coefficient, level, factor, group, value}` with a `settings`
# attribute, a `kind` attribute (each surfaced coefficient tagged reliability vs
# common-variance), and class `efa_reliability`. `x` is a single coefficient matrix (one
# unnamed group) or a named list of them (one entry per group); each matrix has factor
# rows (the general factor `"g"` first, then the group factors) and columns from the
# `.reliability_core()` menu. Only registry columns are carried, so the same builder serves
# a full multi-factor matrix and a single-factor `g`-only matrix alike; NA cells (a
# structurally undefined index such as ECV on a group row, or a Heywood / empty-factor
# coefficient the core has already warned about) are dropped, so the result holds only
# realized coefficients. Purely a reshaping helper -- it computes nothing.
.reliability_result <- function(x, settings = NULL) {

  if (is.list(x)) {
    # A list is always a multigroup result: every group gets a distinct label so the
    # print renders one block per group. Unnamed or blank entries fall back to their
    # position -- NA is reserved for the single, ungrouped matrix case below.
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

    # A regular (coefficient x factor) grid: each kept coefficient contributes one value
    # per factor. Built in one shot -- the columns of `mat[, keep$core]` unroll
    # column-major, matching the coefficient-major row order -- so a matrix with no
    # surfaced columns yields a 0-row frame rather than a NULL that would break the
    # assembly below.
    # The general factor is the first row of every matrix the core builds, whatever that row
    # is labelled: it is "g" for a solution with group factors, and the factor's own name for
    # a single-factor one, which reports that factor and no group factors at all. Read by
    # position rather than by that label, which also keeps a group factor a user happens to
    # name "g" at level "group". Vectorized over the row index so a matrix with no rows
    # yields no levels, matching the empty coefficient and factor columns below.
    row_level <- ifelse(seq_along(factors) == 1L, "general", "group")

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

      factor_corres <- matrix(0, nrow = nrow(s_load), ncol= ncol(s_load))

      for(i in seq_len(nrow(s_load))){

        # seq_len(k) + 1, not 2:(k + 1): the latter counts down to c(2, 1) with no group
        # factors and selects the general column, as .reliability_core() documents.
        factor_corres[i, which.max(abs(input[i, seq_len(ncol(s_load)) + 1L]))] <- 1

      }



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
