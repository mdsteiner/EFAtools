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
# matrix (`cormat`). Given these, it computes
# McDonald's omega total, hierarchical, and subscale for the general factor and each group
# factor, and, when `add_ind = TRUE`, the H index (Hancock & Mueller, 2001) together with the
# ECV and PUC bifactor indices (Rodriguez, Reise & Haviland, 2016). Two total-variance
# conventions are supported: `"correlation"` uses the correlation matrix; `"sums_load"` uses
# the model-implied composite variance from the squared loading-column sums and the
# uniquenesses, McDonald's (1999, Test Theory, Eq. 6.20a) model-implied total. When
# `add_rel = TRUE`, three further columns are appended: standardized Cronbach's alpha
# (Cronbach, 1951) for the whole scale and each subscale, and -- per group factor, over its
# assigned (simple-structure) composite -- the composite reliability (congeneric omega;
# Joreskog, 1971; Raykov, 2001) and the average variance extracted (AVE, a convergent-validity
# index rather than a reliability; Fornell & Larcker, 1981). The function is purely
# computational: `spec` is assumed already normalized by the calling adapter, which owns all
# front-end input handling.
.reliability_core <- function(spec, variance = c("correlation", "sums_load"),
                              add_ind = TRUE, add_rel = FALSE) {

  variance <- match.arg(variance)

  g_load <- spec$g_load
  s_load <- spec$s_load
  u2 <- spec$u2
  factor_corres <- spec$map
  cormat <- spec$cormat
  var_names <- spec$var_names

  # Same general + group factor labels regardless of how s_load was obtained.
  factor_names <- c("g", seq_len(ncol(s_load)))

  # Create an input dataframe
  input <- data.frame(g_load, s_load)
  colnames(input) <- factor_names
  rownames(input) <- var_names
  input$u2 <- u2

  omega_mat <- input[, 2:(ncol(s_load) + 1), drop = FALSE] * factor_corres

  # Sum of all g loadings
  sum_g <- sum(input$g)

  # Sum of all error variances
  sum_e <- sum(input$u2)

  # Compute sums of error variances and g-loadings for group factors
  sums_e_s <- NULL
  sums_g_s <- NULL

  for (i in seq_len(ncol(s_load))){
    sums_e_s[i] <- sum(input$u2 * factor_corres[, i])
    sums_g_s[i] <- sum(input$g * factor_corres[, i])
  }

  sums_s_s <- colSums(omega_mat)

  if(variance == "correlation"){

    # Compute omega total, hierarchical, and subscale for g-factor
    omega_tot_g <- (sum(cormat) - sum_e) / sum(cormat)
    omega_h_g <- sum_g^2 / sum(cormat)
    omega_sub_g <- sum(sums_s_s^2) / sum(cormat)

    # Compute omega total, hierarchical, and subscale for group factors
    omega_tot_sub <- NULL
    omega_h_sub <- NULL
    omega_sub_sub <- NULL

    for (i in seq_len(ncol(s_load))) {
      subf <- which(factor_corres[, i] == 1)
      Vgr <- sum(cormat[subf, subf])
      omega_sub_sub[i] <- sums_s_s[i]^2 / Vgr
      omega_h_sub[i] <- sums_g_s[i]^2 / Vgr
      omega_tot_sub[i] <- (sums_s_s[i]^2 + sums_g_s[i]^2) / Vgr
    }

  } else if(variance == "sums_load") {

    # Sums of all group factor loadings for all group factors
    sums_s <- colSums(input[, 2:(ncol(s_load) + 1), drop = FALSE])

    # Compute omega total, hierarchical, and subscale for the whole scale. The
    # composite here is all items, so every item's loading on every group factor
    # contributes to its variance: the general and group terms use the full
    # loading-column sums. omega total is then McDonald's model-implied total,
    # 1 - sum(u^2) / V, with V the model-implied composite variance, and the
    # general and group variances partition it exactly (tot = hier + sub for the
    # g row). McDonald (1999, Test Theory, Eq. 6.20a).
    omega_tot_g <- (sum_g^2 + sum(sums_s^2)) / (sum_g^2 + sum(sums_s^2) + sum_e)
    omega_h_g <- sum_g^2 / (sum_g^2 + sum(sums_s^2) + sum_e)
    omega_sub_g <- sum(sums_s^2) / (sum_g^2 + sum(sums_s^2) + sum_e)

    # Compute omega total, hierarchical, and subscale for group factors
    omega_tot_sub <- (sums_g_s^2 + sums_s_s^2) / (sums_g_s^2 + sums_s_s^2 +
                                                         sums_e_s)
    omega_h_sub <- sums_g_s^2 / (sums_g_s^2 + sums_s_s^2 + sums_e_s)
    omega_sub_sub <- sums_s_s^2 / (sums_g_s^2 + sums_s_s^2 + sums_e_s)
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
      c("Some group factors have no assigned items in {.arg factor_corres}.",
        "i" = "Their omega coefficients and H index are returned as {.val NA}."),
      class = "efa_omega_empty_factor"
    )
  }

  # Combine and display results in a table
  omega_tot <- c(omega_tot_g, omega_tot_sub)
  omega_h <- c(omega_h_g, omega_h_sub)
  omega_sub <- c(omega_sub_g, omega_sub_sub)

  # In "correlation" mode the total-variance denominator (the correlation matrix)
  # is independent of the loading-based numerators, so a correlation matrix that
  # is inconsistent with the supplied loadings and uniquenesses can push a
  # coefficient above 1 or make omega hierarchical exceed omega total. These are
  # not admissible reliabilities; warn rather than return them silently.
  tol <- .Machine$double.eps^0.5
  if (variance == "correlation" &&
      (any(c(omega_tot, omega_h, omega_sub) > 1 + tol, na.rm = TRUE) ||
       any(omega_h > omega_tot + tol, na.rm = TRUE))) {
    cli::cli_warn(
      c("Some omega coefficients are out of range (above 1, or omega hierarchical above omega total).",
        "i" = "Check that {.arg cormat} is consistent with the loadings and uniquenesses, or use {.code variance = \"sums_load\"}."),
      class = "efa_omega_out_of_range"
    )
  }

  # Optional reliability / convergent-validity coefficients, appended as columns when
  # requested. Kept behind add_rel so the OMEGA coefficient menu is unchanged otherwise.
  if (isTRUE(add_rel)) {

    s_mat <- as.matrix(input[, 2:(ncol(s_load) + 1), drop = FALSE])

    # Standardized Cronbach's alpha (Cronbach, 1951), k / (k - 1) *
    # (1 - sum(diag(R_sub)) / sum(R_sub)), for the whole scale (all items) and each
    # subscale (the map's assigned items). It needs a correlation matrix: the supplied
    # cormat, or -- when none was given -- the model-implied Lambda Lambda' + diag(u2)
    # built from the (orthogonal general + group) loadings and uniquenesses (which
    # assumes the model holds), standardized to a unit diagonal with cov2cor so the
    # formula returns the standardized coefficient even when the loadings and
    # uniquenesses do not complete to unit item variance. Standardized alpha uses the
    # correlation matrix; raw, covariance-based alpha would need item standard deviations
    # the spec does not carry. Alpha is undefined for fewer than two items (returned NA).
    R_rel <- if (!is.null(cormat)) {
      cormat
    } else {
      Lambda <- cbind(input$g, s_mat)
      stats::cov2cor(Lambda %*% t(Lambda) + diag(u2, nrow = length(u2)))
    }

    composites <- c(list(seq_len(nrow(s_mat))),
                    lapply(seq_len(ncol(s_load)),
                           function(j) which(factor_corres[, j] == 1)))
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
# "EFAtools"` requires an explicit `factor_corres`.
.rel_map <- function(s_load, factor_corres = NULL, type = c("EFAtools", "psych")) {

  type <- match.arg(type)
  s_load <- as.matrix(s_load)

  if (!is.null(factor_corres)) {
    checkmate::assert_matrix(factor_corres, nrows = nrow(s_load),
                             ncols = ncol(s_load))
    return(factor_corres)
  }

  if (type == "EFAtools") {
    cli::cli_abort(
      "Specify {.arg factor_corres}, or set {.code type = \"psych\"} to derive variable-to-factor correspondences from the highest group-factor loading per variable.",
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

# Abort (classed) when a supplied correlation matrix fails the .is_cormat() check;
# shared guard for the adapters that accept a user cormat.
.rel_check_cormat <- function(cormat) {
  if (!.is_cormat(cormat)) {
    cli::cli_abort(
      c("{.arg cormat} is not a correlation matrix.",
        "i" = "Check the {.arg cormat} input, or leave it {.code NULL}."),
      class = "efa_reliability_not_cormat"
    )
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
    cormat <- .rel_check_cormat(cormat)
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
    cormat <- .rel_check_cormat(cormat)
  }

  if (is.null(fac_names)) fac_names <- s_load_names

  list(g_load = g_load, s_load = s_load, u2 = u2,
       map = .rel_map(s_load, factor_corres, type),
       cormat = cormat, var_names = var_names, fac_names = fac_names)

}

# Adapter: normalize manually supplied SL/bifactor components (g_load, s_load, u2) to
# the reliability spec. Mirrors the model = NULL branch of .OMEGA_FLEX. A cormat may
# be given directly, reconstructed from an oblique pattern matrix and factor
# intercorrelations, or left NULL for variance = "sums_load".
.rel_adapt_manual <- function(g_load, s_load, u2, var_names, factor_corres = NULL,
                              type = "EFAtools", cormat = NULL, pattern = NULL,
                              Phi = NULL, fac_names = NULL) {

  s_load <- as.matrix(s_load)

  if (is.null(cormat)) {
    if (!is.null(Phi) && !is.null(pattern)) {
      cormat <- psych::factor.model(f = pattern, Phi = Phi, U2 = FALSE)
    }
  } else {
    cormat <- .rel_check_cormat(cormat)
  }

  if (is.null(fac_names)) fac_names <- seq_len(ncol(s_load))

  list(g_load = g_load, s_load = s_load, u2 = u2,
       map = .rel_map(s_load, factor_corres, type),
       cormat = cormat, var_names = var_names, fac_names = fac_names)

}

# Adapter: normalize a lavaan single-factor, bifactor, or second-order solution to
# the reliability spec, one entry per fitted group. Mirrors the structural detection
# and Schmid-Leiman transform of .OMEGA_LAVAAN (second-order group loadings via
# .sl_group_loadings), and drives the core with variance = "sums_load" (the composite
# variances are model-implied). Returns a list with, per group, either a full spec
# (bifactor / second-order) or a single-factor marker (`single = TRUE`) carrying the
# loadings and uniquenesses -- a single factor has no group factors and so is scored
# directly (omega = sum(g)^2 / (sum(g)^2 + sum(u2)), H via .h_index), not through the
# multi-factor core.
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

    # A single factor has no group factors: mark it for direct scoring.
    if (ncol(lambda) == 1) {
      groups[[i]] <- list(single = TRUE, g_load = lambda[, 1], u2 = diag(theta),
                          var_names = var_names_i, fac_names = colnames(lambda))
      next
    }

    col_names <- colnames(lambda)

    if (!any(col_names %in% g_name)) {
      cli::cli_abort(
        c("Could not find the specified general-factor name in the lavaan solution.",
          "i" = "Please check the spelling."),
        class = "efa_reliability_g_name"
      )
    }

    # Detect the model type once (all groups share the fixed-zero structure): a
    # general-factor column of zeros marks a second-order model (the general factor
    # loads the first-order factors via `beta`, not the items); otherwise the solution
    # must be a genuine bifactor (at least one item loading on two factors).
    if (i == 1) {
      if (all(lambda[, g_name] == 0)) {
        higherorder <- TRUE
        if (sum(colSums(std_sol[[1]][["beta"]]) > 0) > 1) {
          cli::cli_abort("The higher-order model has more than two latent strata or more than one second-order factor; only second-order models with one second-order factor are supported.",
                         class = "efa_reliability_higher_order")
        }
      } else {
        bi_check <- lambda
        bi_check[abs(bi_check) > tol] <- 1
        if (all(rowSums(bi_check) < 2)) {
          cli::cli_abort(
            c("The lavaan input is invalid; no reliability coefficients are computed.",
              "i" = "Provide a bifactor model, a second-order model, or a single-factor model."),
            class = "efa_reliability_invalid_lavaan"
          )
        }
        # A genuine bifactor has each item loading on the general and a group
        # factor; flag the borderline case where some item loads on fewer than
        # two factors so the front-end can warn the user (it is still scored).
        few_loadings <- !all(rowSums(bi_check) > 1)
      }
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
       higher_order = higherorder, few_loadings = few_loadings)

}

# Direct-scoring omega total and H index for a single-factor solution: a single
# factor has no group factors, so it is scored directly rather than through the
# multi-factor core. omega total is sum(loadings)^2 / (sum(loadings)^2 + sum(u2))
# and H is the Heywood-guarded H index over the loadings. Used by the OMEGA and
# efa_reliability front-ends for the lavaan single-factor path.
.rel_single_factor <- function(g_load, u2) {
  sum_g <- sum(g_load)
  c(omega = sum_g^2 / (sum_g^2 + sum(u2)), H = .h_index(g_load))
}

# Adapter: normalize an oblique EFA() object to the reliability spec, treating it as
# the correlated-factors model it is. A correlated-factors solution identifies
# whole-scale omega total and each factor's congeneric omega/H, but not the
# hierarchical/bifactor indices (omega hierarchical, ECV, PUC), which require a
# general factor. Rather than manufacture one by a Schmid-Leiman transformation --
# which is underidentified for fewer than three factors, imposes proportionality
# constraints, and biases the general loadings under the cross-loadings EFA solutions
# almost always have (Flora, 2020, Adv. Methods Pract. Psychol. Sci.; Mansolf &
# Reise, 2016, Multivariate Behav. Res.) -- the spec carries a zero general factor
# and the oblique pattern as the group loadings. The core then returns omega total
# (1' L Phi L' 1 / 1' R 1) and per-factor congeneric omega/H, with omega
# hierarchical / ECV / PUC 0.
.rel_adapt_efa <- function(model, factor_corres = NULL, type = "psych",
                           cormat = NULL, fac_names = NULL) {

  if (!("Phi" %in% names(model))) {
    cli::cli_abort(
      c("{.arg model} is not an oblique EFA solution.",
        "i" = "Reliability from an EFA needs correlated factors; refit with an oblique rotation."),
      class = "efa_reliability_not_oblique"
    )
  }

  L1 <- unclass(model$rot_loadings)
  Phi <- model$Phi

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
    cormat <- .rel_check_cormat(cormat)
  }

  if (is.null(fac_names)) {
    fac_names <- colnames(s_load)
    # Unlabelled loading columns still need one row label per group factor.
    if (is.null(fac_names)) fac_names <- paste0("F", seq_len(ncol(s_load)))
  }

  list(g_load = rep(0, nrow(s_load)), s_load = s_load, u2 = u2,
       map = .rel_map(s_load, factor_corres, type),
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

  if (is.null(u2)) u2 <- 1 - rowSums(loadings^2)

  if (is.null(factor_corres)) {
    factor_corres <- abs(s_load) > .Machine$double.eps * 100
  }

  if (is.null(cormat)) {
    cormat <- loadings %*% t(loadings) + diag(u2, nrow = length(u2))
  } else {
    cormat <- .rel_check_cormat(cormat)
  }

  var_names <- rownames(loadings)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(nrow(loadings)))

  if (is.null(fac_names)) {
    fac_names <- colnames(s_load)
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
    data.frame(
      coefficient = rep(keep$coefficient, each = length(factors)),
      level = rep(ifelse(factors == "g", "general", "group"), times = nrow(keep)),
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

        factor_corres[i, which.max(abs(input[i, 2:(ncol(s_load) + 1)]))] <- 1

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
    colnames(model)[2:(ncol(s_load) + 1)]
  }

  # Hand the normalized components to the reliability engine.
  spec <- list(g_load = g_load, s_load = s_load, u2 = u2, map = factor_corres,
               cormat = cormat, var_names = var_names, fac_names = fac_names_out)

  .reliability_core(spec, variance = variance, add_ind = add_ind)

}


# Reshape a lavaan reliability solution into OMEGA's legacy output ------------
#
# OMEGA's lavaan path is a thin front-end over the shared reliability machinery:
# the lavaan adapter normalizes the fitted model to a per-group spec (Schmid-
# Leiman transforming a second-order model), the reliability core scores each
# multi-factor group, and single-factor groups are scored directly. This function
# only reassembles those results into OMEGA's historical shapes -- a coefficient
# matrix per group (unwrapped for a single group), a named list for several
# groups, and a named c(Omega, H) vector (or a bare omega) for a single factor --
# and emits the user-facing notes about the model structure. All reliability math
# lives in .reliability_core / .rel_single_factor.
.OMEGA_LAVAAN <- function(model = NULL, g_name = "g", group_names = NULL,
                          add_ind = TRUE){

  adapt <- .rel_adapt_lavaan(model, g_name = g_name, group_names = group_names)

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

      sf <- .rel_single_factor(grp$g_load, grp$u2)
      omegas[[i]] <- if (isTRUE(add_ind)) {
        stats::setNames(c(sf[["omega"]], sf[["H"]]), c("Omega", "H"))
      } else {
        sf[["omega"]]
      }

    } else {

      # A multi-factor (bifactor / Schmid-Leiman) group: score through the core,
      # storing the unclassed matrix -- the OMEGA class is attached to the whole
      # output below, matching the historical list-of-matrices shape. The core
      # labels the general-factor row "g"; restore the user's general-factor name.
      mat <- unclass(.reliability_core(grp, "sums_load", add_ind = add_ind))
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
