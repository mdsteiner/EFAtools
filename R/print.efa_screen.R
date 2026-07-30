#' Print and format an efa_screen object
#'
#' `print()` turns the factor-analysis screening diagnostics computed by
#' [efa_screen()] into a sectioned report with banded, colour-coded verdicts:
#' sampling adequacy and sphericity (the Kaiser-Meyer-Olkin measure and Bartlett's
#' test of sphericity), multicollinearity (the determinant and condition number of
#' the correlation matrix), the per-variable diagnostics, and, when raw data were
#' supplied, multivariate normality and multivariate outliers. It closes with a
#' consolidated list of actionable recommendations (for example, which items to
#' consider dropping, whether to prefer an ordinal or a robust estimator, and a
#' caveat that keeps an over-powered Bartlett's test from being over-trusted).
#' `format()` assembles the same report and returns it as a character vector;
#' `print()` is `cat(format(x), sep = "\n")`. The lines follow the active console
#' theme, so they are plain when colours are disabled (for example when captured
#' into a file or stripped with [cli::ansi_strip()]). `print()` does not draw a
#' plot.
#'
#' @param x An object of class `efa_screen` (output from [efa_screen()]).
#' @param digits Integer. The number of decimal places the reported values are
#'   rounded to. Default is 3.
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines (styled to the active console theme;
#'   plain when colours are disabled).
#'
#' @family factor analysis suitability
#'
#' @export
#'
#' @method print efa_screen
#'
#' @examples
#' # From raw data
#' efa_screen(iris[, 1:4])
#'
#' # From a correlation matrix (supply N for Bartlett's test of sphericity)
#' efa_screen(test_models$baseline$cormat, N = 500)
#'
#' # format() returns the same lines as a character vector:
#' writeLines(format(efa_screen(test_models$baseline$cormat, N = 500)))
#'
print.efa_screen <- function(x, digits = 3, ...) {
  cat(format(x, digits = digits, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.efa_screen
#' @export
#' @method format efa_screen
format.efa_screen <- function(x, digits = 3, ...) {

  raw <- !is.null(x$per_item)

  # Overall KMO band (reused Kaiser & Rice, 1974 cut-offs) and the significance flags
  # that several sections and the recommendations share.
  kmo <- x$kmo$KMO
  bart_sig <- .screen_is_sig(x$bartlett$p_value)
  mvn_nonnormal <- raw && !inherits(x$normality, "efa_screen_no_mvn") &&
    (.screen_is_sig(x$normality$mardia$skewness_p) ||
       .screen_is_sig(x$normality$mardia$kurtosis_p) ||
       .screen_is_sig(x$normality$hz$p_value))

  cli::cli_format_method({

    # -- Sampling adequacy and sphericity --------------------------------------
    .print_efa_rule("Sampling adequacy and sphericity")

    if (!is.null(kmo) && !is.na(kmo)) {
      band <- .kmo_band(kmo)
      kval <- round(kmo, digits)
      band$alert("The overall KMO value for your data is {band$label} (Overall KMO = {kval}).",
                 wrap = TRUE)
      cli::cli_text("These data are {band$suitability} suitable for factor analysis.")
    } else {
      cli::cli_alert_warning("The overall KMO value for your data is not available.",
                             wrap = TRUE)
    }

    cli::cli_text("")

    bart <- x$bartlett
    if (is.null(bart)) {
      cli::cli_alert_warning(
        "Bartlett's test of sphericity was not computed; no sample size (N) was supplied.",
        wrap = TRUE)
    } else if (!is.null(bart$p_value) && !is.na(bart$p_value)) {
      # Bartlett significance wording reused verbatim from format.efa_bartlett.
      if (bart$p_value < .05) {
        cli::cli_alert_success(
          "The Bartlett's test of sphericity was significant at an alpha level of .05.",
          wrap = TRUE)
        cli::cli_text("These data are probably suitable for factor analysis.")
      } else {
        cli::cli_alert_danger(
          "The Bartlett's test of sphericity was not significant at an alpha level of .05.",
          wrap = TRUE)
        cli::cli_text("These data are probably not suitable for factor analysis.")
      }
      cli::cli_verbatim(paste0("\U03C7\U00B2(", bart$df, ") = ", round(bart$chisq, 2),
                               ", ", cli::style_italic("p"), .screen_p_str(bart$p_value)))
    } else {
      cli::cli_alert_warning("The Bartlett's test of sphericity did not render a result.",
                             wrap = TRUE)
    }

    # -- Multicollinearity -----------------------------------------------------
    .print_efa_rule("Multicollinearity")

    # Multicollinearity thresholds computed once here (a near-zero determinant, Field
    # 2018; a condition index above 30, Belsley-Kuh-Welsch 1980) drive both the section
    # verdicts and the consolidated recommendation, so the cut-offs live in one place.
    det_R <- x$determinant
    det_low <- !is.na(det_R) && det_R < 1e-5
    dstr <- format(det_R, digits = digits, scientific = FALSE)
    if (det_low) {
      cli::cli_alert_danger(
        "Determinant: {dstr}. Multicollinearity likely (a value near 0 signals it).",
        wrap = TRUE)
    } else {
      cli::cli_alert_success(
        "Determinant: {dstr}. No concern (a value near 0 signals multicollinearity).",
        wrap = TRUE)
    }

    cond <- x$condition
    cstr <- formatC(cond, format = "f", digits = digits)
    # The 10/30 cut-offs are Belsley-Kuh-Welsch thresholds for the condition *index*, the
    # square root of the condition number, so the index is shown alongside the number the
    # section reports: quoting the rule without the quantity it governs invites the value
    # to be read against the wrong scale.
    ci <- sqrt(cond)
    cistr <- formatC(ci, format = "f", digits = digits)
    ci_strong <- !is.na(ci) && ci > 30
    if (ci_strong) {
      cli::cli_alert_danger(
        "Condition number: {cstr} (condition index {cistr}). Strong multicollinearity (index above 30; Belsley, Kuh & Welsch, 1980).",
        wrap = TRUE)
    } else if (!is.na(ci) && ci > 10) {
      cli::cli_alert_warning(
        "Condition number: {cstr} (condition index {cistr}). Moderate multicollinearity (index 10 to 30; Belsley, Kuh & Welsch, 1980).",
        wrap = TRUE)
    } else {
      cli::cli_alert_success(
        "Condition number: {cstr} (condition index {cistr}). No concern (index below 10; Belsley, Kuh & Welsch, 1980).",
        wrap = TRUE)
    }

    # -- Per-variable diagnostics ----------------------------------------------
    .print_efa_rule("Per-variable diagnostics")
    # The standard display names are MSA (the per-variable measure of sampling adequacy,
    # stored as kmo_i) and SMC (the squared multiple correlation). The per-variable
    # values are laid out as base R's data-frame print, emitted verbatim so cli does not
    # reflow the aligned columns (as format.efa_kmo does for the KMO vector). This is the
    # one table not routed through .efa_num(), so it is also the one place a leading zero
    # is kept: the frame mixes bounded coefficients with variances and missing-data
    # percentages on the data's own scale, for which the leading zero is not redundant.
    disp <- .screen_per_item_display(x, digits)
    cli::cli_verbatim(utils::capture.output(print(disp)))

    # -- Multivariate normality (raw data only) --------------------------------
    if (raw) {
      .print_efa_rule("Multivariate normality")
      nm <- x$normality
      if (inherits(nm, "efa_screen_no_mvn")) {
        cli::cli_alert_warning(
          "Multivariate normality tests were skipped (singular complete-case covariance).",
          wrap = TRUE)
      } else {
        md <- nm$mardia
        .screen_mvn_line(.screen_is_sig(md$skewness_p),
                         paste0("Mardia's skewness: \U03C7\U00B2(", md$skewness_df, ") = ",
                                round(md$skewness, 2), ", p", .screen_p_str(md$skewness_p), "."))
        .screen_mvn_line(.screen_is_sig(md$kurtosis_p),
                         paste0("Mardia's kurtosis: z = ", round(md$kurtosis, 2),
                                ", p", .screen_p_str(md$kurtosis_p), "."))
        .screen_mvn_line(.screen_is_sig(nm$hz$p_value),
                         paste0("Henze-Zirkler: HZ = ", round(nm$hz$statistic, 2),
                                ", p", .screen_p_str(nm$hz$p_value), "."))
        if (isTRUE(md$small_sample)) {
          cli::cli_text("(Mardia's skewness uses the small-sample correction for n < 20.)")
        }
        cli::cli_text(if (isTRUE(mvn_nonnormal)) "These data depart from multivariate normality."
                      else "These data are consistent with multivariate normality.")
      }
    }

    # -- Outliers (raw data only) ----------------------------------------------
    if (raw) {
      .print_efa_rule("Outliers")
      o <- x$outliers
      if (inherits(o, "efa_screen_no_outliers")) {
        cli::cli_alert_warning("Outlier diagnostics were skipped.")
        # The recorded reason, as in the classical-fallback branch below: which of the
        # three ways the complete-case covariance failed decides what the user should do,
        # and "singular" alone points at collinearity even when the cause is missingness.
        if (!is.null(o$reason)) cli::cli_text(o$reason)
      } else {
        if (identical(o$method, "classical")) {
          cli::cli_alert_warning(paste(
            "A robust (MCD) covariance could not be computed; classical Mahalanobis",
            "distances were used."), wrap = TRUE)
          # The recorded reason (too few complete cases, or an exact fit) and the
          # consequence of the fallback: the classical covariance is computed from every
          # observation, outliers included, so the diagnostic no longer has the
          # high-breakdown property the MCD gives it.
          if (!is.null(o$fallback_reason)) cli::cli_text(o$fallback_reason)
          cli::cli_text(paste(
            "These distances come from a covariance the outliers themselves inflate, so",
            "the diagnostic is no longer high-breakdown and tends to under-flag."))
        }
        nf <- length(o$flagged)
        cstr_o <- format(o$cutoff, digits = digits, scientific = FALSE)
        # Name the distance for the estimate actually used (robust MCD, or classical
        # Mahalanobis when the robust covariance could not be formed).
        dist_label <- if (identical(o$method, "classical")) "Mahalanobis distance" else "robust distance"
        cli::cli_alert_info(paste0(
          nf, " of ", o$n_complete, " observations ",
          .screen_plural(nf, "was", "were"), " flagged as ",
          .screen_plural(nf, "a multivariate outlier", "multivariate outliers"),
          " (", dist_label, " > ", cstr_o, ")."), wrap = TRUE)
      }
    }

    # -- Recommendations -------------------------------------------------------
    .print_efa_rule("Recommendations")
    recs <- .screen_recommendations(x, raw, kmo, bart_sig, mvn_nonnormal,
                                    det_low || ci_strong, digits)
    # cli treats each bullet as a glue/cli template, so double any braces coming from
    # user variable names to render them literally (as .efa_emit_bullets does).
    msg <- gsub("}", "}}", gsub("{", "{{", recs$message, fixed = TRUE), fixed = TRUE)
    cli::cli_bullets(stats::setNames(msg, recs$symbol))

  })
}

# Kaiser & Rice (1974) verbal bands for the overall KMO value, and the single place they
# are written down: the highest band the value clears gives its label, the suitability it
# implies, and the severity (success >= .7, warning >= .6, danger below) as an alert
# function, the matching colour, and the equivalent cli_bullets() symbol, so that callers
# reporting the same value as an alert, a coloured label, or a bullet cannot drift apart.
.kmo_band <- function(kmo) {
  bands <- list(
    list(min = .9,   label = "marvellous",   colour = cli::col_green,  alert = cli::cli_alert_success, symbol = "v", suitability = "probably"),
    list(min = .8,   label = "meritorious",  colour = cli::col_green,  alert = cli::cli_alert_success, symbol = "v", suitability = "probably"),
    list(min = .7,   label = "middling",     colour = cli::col_green,  alert = cli::cli_alert_success, symbol = "v", suitability = "probably"),
    list(min = .6,   label = "mediocre",     colour = cli::col_yellow, alert = cli::cli_alert_warning, symbol = "!", suitability = "probably"),
    list(min = .5,   label = "miserable",    colour = cli::col_red,    alert = cli::cli_alert_danger,  symbol = "x", suitability = "hardly"),
    list(min = -Inf, label = "unacceptable", colour = cli::col_red,    alert = cli::cli_alert_danger,  symbol = "x", suitability = "not")
  )
  Find(function(b) kmo >= b$min, bands)
}

# TRUE when a p-value is present and below .05 (guards NA / NULL p-values).
.screen_is_sig <- function(p) !is.null(p) && !is.na(p) && p < .05

# p-value tail formatting shared by the Bartlett and MVN lines and by format.efa_bartlett.
.screen_p_str <- function(p) {
  if (is.null(p) || is.na(p)) " = NA" else if (p < .001) " < .001" else paste0(" = ", round(p, 3))
}

# Singular/plural picker for a scalar count (avoids cli's numeric pluralizer).
.screen_plural <- function(n, singular, plural) if (n == 1L) singular else plural

# One multivariate-normality test line: a danger alert when significant, a success alert
# otherwise. `text` already carries the test statistic and p-value (no cli braces). Left
# unwrapped, unlike the prose alerts around it: it is a fixed-format statistic line, and
# wrapping would let a break fall inside the statistic.
.screen_mvn_line <- function(sig, text) {
  if (sig) cli::cli_alert_danger(text) else cli::cli_alert_success(text)
}

# Oxford-comma list of names ("a", "a and b", "a, b, and c").
.screen_and_list <- function(x) {
  if (length(x) == 1L) return(x)
  if (length(x) == 2L) return(paste(x, collapse = " and "))
  paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
}

# Build the per-variable display data frame: the full per-item table for raw data
# (columns renamed to the standard MSA/SMC abbreviations), or MSA and SMC alone for a
# correlation-matrix input. Numeric columns are rounded to `digits`.
.screen_per_item_display <- function(x, digits) {
  if (!is.null(x$per_item)) {
    disp <- x$per_item
    num <- vapply(disp, is.numeric, logical(1))
    disp[num] <- lapply(disp[num], round, digits)
    names(disp)[names(disp) == "smc"] <- "SMC"
    names(disp)[names(disp) == "kmo_i"] <- "MSA"
    # `missing` is a percentage; the bare header reads as a count next to `variance`,
    # which is on the data's own scale.
    names(disp)[names(disp) == "missing"] <- "missing%"
    # A `flags` entry is NA for a variable treated as continuous, i.e. "no category
    # screening applies". Printed as <NA> that reads as a failed computation, so drop
    # the column when it says nothing at all and render the remaining NAs as a dash.
    if (all(is.na(disp$flags))) disp$flags <- NULL
    else disp$flags[is.na(disp$flags)] <- "-"
    disp
  } else {
    # make.unique() guards against duplicate variable names (a correlation matrix with
    # duplicate dimnames is a valid input), which base data.frame() row names reject; it
    # matches how the raw per-item table is keyed (efa_screen()).
    data.frame(MSA = round(x$kmo$KMO_i, digits),
               SMC = round(x$smc, digits),
               row.names = make.unique(names(x$kmo$KMO_i)))
  }
}

# Assemble the consolidated recommendations as parallel symbol/message vectors. Each
# check fires only when it applies; when no warning fires a single success line is
# returned, and a correlation-matrix input always closes with a raw-data note. Every
# recommendation is anchored in the literature cited in efa_screen()'s @source.
.screen_recommendations <- function(x, raw, kmo, bart_sig, mvn_nonnormal, multicollinear,
                                    digits) {
  sym <- character(0)
  msg <- character(0)
  push <- function(s, m) { sym[[length(sym) + 1L]] <<- s; msg[[length(msg) + 1L]] <<- m }

  # Few-category items (fewer than 5 response categories) are detected up front because
  # they steer the estimator recommendation away from the continuous-non-normal advice.
  few <- character(0)
  if (raw && !is.null(x$categories)) {
    n_cat <- vapply(x$categories, function(ct) {
      if (length(ct) == 1L && is.na(ct[1])) NA_integer_ else length(ct)
    }, integer(1))
    few <- names(n_cat)[!is.na(n_cat) & n_cat < 5L]
  }

  # Low overall sampling adequacy (KMO < .6): the set may not be factorable.
  if (!is.null(kmo) && !is.na(kmo) && kmo < .6) {
    push("!", paste0("Overall sampling adequacy is ", .kmo_band(kmo)$label,
                     " (KMO = ", round(kmo, digits),
                     "); the variables may not form a factorable set."))
  }

  # Variables with a low individual MSA (< .5) share little variance and are drop
  # candidates (Kaiser & Rice, 1974). make.unique() keys the names as the per-item table
  # and the sparse/empty recommendations do, so a variable is labelled consistently even
  # when the input has duplicate variable names.
  msa <- x$kmo$KMO_i
  low <- make.unique(names(msa))[!is.na(msa) & msa < .5]
  if (length(low)) {
    push("!", paste0(length(low), " ", .screen_plural(length(low), "variable has", "variables have"),
                     " a low individual MSA (< .5): ", .screen_and_list(low),
                     "; consider removing ", .screen_plural(length(low), "it", "them"),
                     " (little shared variance)."))
  }

  # Few-category ordinal data: prefer a categorical estimator (Rhemtulla et al., 2012).
  if (length(few)) {
    push("!", paste0(length(few), " ", .screen_plural(length(few), "item has", "items have"),
                     " fewer than 5 response categories; treating ",
                     .screen_plural(length(few), "it", "them"),
                     " as ordinal (polychoric correlations with a categorical estimator such as DWLS)",
                     " is less biased than normal-theory ML."))
  }

  # Continuous / 5+-category non-normal data: normal-theory standard errors are biased.
  # Suppressed when few-category items already routed the user to an ordinal estimator.
  if (isTRUE(mvn_nonnormal) && !length(few)) {
    push("!", paste("These data depart from multivariate normality; normal-theory standard errors",
                    "and fit statistics may be biased - prefer robust (sandwich) or bootstrapped",
                    "standard errors."))
  }

  # An over-powered Bartlett's test: significant but built on the normality it violates.
  if (bart_sig && isTRUE(mvn_nonnormal)) {
    push("!", paste("Bartlett's test is significant, but it assumes multivariate normality and grows",
                    "more sensitive as N increases; because these data are non-normal, treat it as",
                    "uninformative here and rely on the KMO."))
  }

  # Multicollinearity: a near-zero determinant (Field, 2018) or a high condition index
  # (Belsley, Kuh & Welsch, 1980), decided once by the caller from the same thresholds
  # used for the section verdicts.
  if (isTRUE(multicollinear)) {
    push("!", paste("A near-zero determinant or high condition number indicates multicollinearity;",
                    "look for redundant or linearly dependent items (r > .8) and consider removing them."))
  }

  if (raw) {
    flags <- x$per_item$flags
    items <- rownames(x$per_item)
    sparse <- items[!is.na(flags) & grepl("sparse", flags)]
    if (length(sparse)) {
      push("!", paste0(length(sparse), " ", .screen_plural(length(sparse), "variable has", "variables have"),
                       " a sparse response category (< 5 responses): ", .screen_and_list(sparse),
                       "; a low-frequency category can destabilise polychoric estimates",
                       " - consider collapsing it into an adjacent category."))
    }
    empty <- items[!is.na(flags) & grepl("empty", flags)]
    if (length(empty)) {
      push("!", paste0(length(empty), " ", .screen_plural(length(empty), "variable has", "variables have"),
                       " an unused interior response category (a gap in the scale): ",
                       .screen_and_list(empty), "; check the coding."))
    }
    if (!inherits(x$outliers, "efa_screen_no_outliers")) {
      nf <- length(x$outliers$flagged)
      nominal <- 1 - (x$settings$outlier_cutoff %||% 0.975)
      frac <- nf / x$outliers$n_complete
      # Far more flags than the cutoff nominally admits is not a contamination count to
      # work through: a high-breakdown estimate fitted to the most concentrated half of a
      # clustered or mixture sample legitimately calls the rest distant, so the excess is
      # evidence that the data are not elliptically distributed. The rate alone is not
      # enough to draw that conclusion - at a tight cutoff, or in a small sample, a couple
      # of genuine outliers clear it - so the flagged set must also be too long to work
      # through case by case, which is what makes the ordinary advice below unhelpful.
      if (nf >= 10L && frac > 5 * nominal) {
        push("!", paste0(nf, " of ", x$outliers$n_complete, " observations (",
                         round(100 * frac), "%) exceed the outlier cutoff, far above the ",
                         # keep a significant digit: a tight cutoff has a nominal rate well
                         # below 0.1%, which rounds to a self-defeating "0%"
                         format(signif(100 * nominal, 2), trim = TRUE, scientific = FALSE),
                         "% expected under multivariate normality; this usually means the",
                         " data are not elliptically distributed (subgroups or a mixture)",
                         " rather than that this many cases are contaminated."))
      } else if (nf > 0L) {
        push("!", paste0(nf, " ", .screen_plural(nf, "observation was", "observations were"),
                         " flagged as ",
                         .screen_plural(nf, "a potential multivariate outlier", "potential multivariate outliers"),
                         "; inspect ", .screen_plural(nf, "it", "them"),
                         " (see `$outliers$flagged`) before down-weighting or excluding."))
      }
    }
  }

  # A clean bill of health when nothing above fired; the raw-data note always closes a
  # correlation-matrix report (its raw-only sections were absent above).
  if (!length(msg)) {
    push("v", "The data appear suitable for factor analysis.")
  }
  if (!raw) {
    push("i", paste("Per-item variance, missing-data, category, normality, and outlier diagnostics",
                    "require raw data; only a correlation matrix was supplied."))
  }

  list(symbol = sym, message = msg)
}
