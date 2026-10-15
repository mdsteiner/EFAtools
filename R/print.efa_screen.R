#' Print and format an efa_screen object
#'
#' `print()` turns the diagnostics computed by [efa_screen()] into a sectioned report
#' with colour-coded verdicts: sampling adequacy and sphericity (the Kaiser-Meyer-Olkin
#' measure and Bartlett's test of sphericity), multicollinearity (the condition index),
#' the per-variable diagnostics, and, when raw data were supplied, multivariate normality
#' and multivariate outliers. It closes with a list of recommendations (for example,
#' which items to consider dropping, or whether to prefer an ordinal or a robust
#' estimator). `format()` assembles the same report and returns it as a character
#' vector; `print()` is `cat(format(x), sep = "\n")`. The lines follow the active console
#' theme, so they are plain when colours are disabled (for example when captured into a
#' file or stripped with [cli::ansi_strip()]). `print()` does not draw a plot.
#'
#' @param x An object of class `efa_screen` (output from [efa_screen()]).
#' @param digits Integer. The number of decimal places the reported values are
#'   rounded to. Default is 3.
#' @param ... Not used; for consistency with the generic.
#'
#' @returns `print()` returns its argument `x` invisibly. `format()` returns a
#'   character vector with the report lines.
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
  kmo <- x$kmo$KMO
  # The normality verdict counts the available tests: the kurtosis statistic is NA in a very
  # small sample.
  mvn_p <- if (raw && !inherits(x$normality, "efa_screen_no_mvn")) {
    c(x$normality$mardia$skewness_p, x$normality$mardia$kurtosis_p)
  }
  mvn_n_avail <- sum(!is.na(mvn_p))
  mvn_n_rej <- sum(vapply(mvn_p, .screen_is_sig, logical(1)))
  mvn_nonnormal <- mvn_n_rej > 0L
  ci_band <- .ci_band(sqrt(x$condition))

  cli::cli_format_method({

    # -- Sampling adequacy and sphericity --------------------------------------
    .print_efa_rule("Sampling adequacy and sphericity")
    .print_kmo_line(kmo, digits)
    .print_bartlett_line(x$bartlett)

    # -- Multicollinearity -----------------------------------------------------
    .print_efa_rule("Multicollinearity")

    cistr <- formatC(sqrt(x$condition), format = "f", digits = digits)
    if (ci_band$flag) {
      ci_band$alert(
        "Condition index: {cistr}, above 30: a {ci_band$strength} near linear dependency (Belsley, 1991).",
        wrap = TRUE)
    } else {
      ci_band$alert("Condition index: {cistr}. {ci_band$note}", wrap = TRUE)
    }

    # -- Per-variable diagnostics ----------------------------------------------
    .print_efa_rule("Per-variable diagnostics")
    # Emitted verbatim so cli does not reflow the aligned columns.
    disp <- .screen_per_item_display(x, digits)
    cli::cli_verbatim(utils::capture.output(print(disp)))

    # -- Multivariate normality (raw data only) --------------------------------
    if (raw) {
      .print_efa_rule("Multivariate normality")
      nm <- x$normality
      if (inherits(nm, "efa_screen_no_mvn")) {
        cli::cli_alert_info(nm$message, wrap = TRUE)
      } else {
        md <- nm$mardia
        .screen_mvn_line(.screen_is_sig(md$skewness_p),
                         paste0("Mardia's skewness: \U03C7\U00B2(", md$skewness_df, ") = ",
                                round(md$skewness, 2), ", p", .screen_p_str(md$skewness_p), "."))
        # An undefined kurtosis statistic is neither a pass nor a rejection.
        if (is.na(md$kurtosis)) {
          cli::cli_alert_info(
            "Mardia's kurtosis: not available at this number of complete cases.")
        } else {
          .screen_mvn_line(.screen_is_sig(md$kurtosis_p),
                           paste0("Mardia's kurtosis: z = ", round(md$kurtosis, 2),
                                  ", p", .screen_p_str(md$kurtosis_p), "."))
        }
        both <- mvn_n_avail == 2L
        cli::cli_text(if (mvn_nonnormal) {
          paste0("These data depart from multivariate normality: ",
                 if (!both) "the test rejects it."
                 else if (mvn_n_rej == 2L) "both tests reject it."
                 else "one of the two tests rejects it.")
        } else {
          paste0("These data are consistent with multivariate normality: ",
                 if (both) "neither test rejects it." else "the test does not reject it.")
        })
        # With missing data the tests describe the complete cases only, so name them.
        if (.screen_incomplete(nm$n_complete, x$settings$n_obs)) {
          cli::cli_text("(Computed from {nm$n_complete} complete cases of the
                         {x$settings$n_obs} rows supplied.)")
        }
      }
    }

    # -- Outliers (raw data only) ----------------------------------------------
    if (raw) {
      .print_efa_rule("Outliers")
      o <- x$outliers
      if (inherits(o, "efa_screen_no_outliers")) {
        cli::cli_alert_warning(o$message, wrap = TRUE)
      } else if (identical(o$method, "ranking")) {
        top <- names(sort(o$distances, decreasing = TRUE))[seq_len(min(5L, length(o$distances)))]
        cli::cli_alert_info(
          "No case is flagged: some variables have fewer than ten distinct values (categorical), and a chi-square cutoff is not valid for their Mahalanobis distances. The largest distances are in rows {top}.",
          wrap = TRUE)
        if (.screen_incomplete(o$n_complete, x$settings$n_obs)) {
          cli::cli_text("(Computed from {o$n_complete} complete cases of the
                         {x$settings$n_obs} rows supplied.)")
        }
      } else {
        classical <- identical(o$method, "classical")
        if (classical) {
          cli::cli_alert_warning(
            "No robust (MCD) estimate: {o$fallback_reason}. Classical Mahalanobis distances were used; they tend to under-flag.",
            wrap = TRUE)
        }
        nf <- length(o$flagged)
        cstr_o <- format(o$cutoff, digits = digits, scientific = FALSE)
        dist_label <- if (classical) "Mahalanobis distance" else "robust distance"
        cli::cli_alert_info(paste0(
          nf, " of ", .screen_case_base(o$n_complete, x$settings$n_obs), " ",
          .screen_plural(nf, "was", "were"), " flagged as ",
          .screen_plural(nf, "a multivariate outlier", "multivariate outliers"),
          " (", dist_label, " > ", cstr_o, ")."), wrap = TRUE)
      }
    }

    # -- Recommendations -------------------------------------------------------
    .print_efa_rule("Recommendations")
    recs <- .screen_recommendations(x, raw, kmo, mvn_nonnormal, ci_band$flag, digits)
    # cli reads each bullet as a template, so double the braces that variable names can
    # carry (as .efa_emit_bullets does).
    msg <- gsub("}", "}}", gsub("{", "{{", recs$message, fixed = TRUE), fixed = TRUE)
    cli::cli_bullets(stats::setNames(msg, recs$symbol))

  })
}

# The overall KMO line, shared by format.efa_screen() and format.efa_retain().
.print_kmo_line <- function(kmo, digits) {
  if (!is.null(kmo) && !is.na(kmo)) {
    band <- .kmo_band(kmo)
    kval <- round(kmo, digits)
    band$alert("Overall KMO = {kval}: {band$label} (Kaiser & Rice, 1974).", wrap = TRUE)
  } else {
    cli::cli_alert_warning("The overall KMO value is not available.", wrap = TRUE)
  }
  invisible(NULL)
}

# The Bartlett's test line, shared by format.efa_screen() and format.efa_retain(). `bart` holds
# `chisq`, `df`, and `p_value`, or is NULL when no N was supplied.
.print_bartlett_line <- function(bart) {
  if (is.null(bart)) {
    cli::cli_alert_warning(
      "Bartlett's test of sphericity was not computed; no sample size (N) was supplied.",
      wrap = TRUE)
  } else if (!is.na(bart$p_value)) {
    # A fixed-format statistic line, so it is not wrapped.
    sig <- bart$p_value < .05
    text <- paste0("Bartlett's test of sphericity: \U03C7\U00B2(", bart$df, ") = ",
                   round(bart$chisq, 2), ", ", cli::style_italic("p"),
                   .screen_p_str(bart$p_value),
                   if (sig) " (significant)." else " (not significant).")
    if (sig) cli::cli_alert_success(text) else cli::cli_alert_danger(text)
  } else {
    .print_bartlett_no_result()
  }
  invisible(NULL)
}

# Kaiser & Rice (1974) verbal bands for the overall KMO value, in one place: the highest band
# the value clears gives its label, the suitability it implies, and its severity as an alert
# function and a colour, so all reports band a value identically.
.kmo_band <- function(kmo) {
  bands <- list(
    list(min = .9,   label = "marvellous",   colour = cli::col_green,  alert = cli::cli_alert_success, suitability = "probably"),
    list(min = .8,   label = "meritorious",  colour = cli::col_green,  alert = cli::cli_alert_success, suitability = "probably"),
    list(min = .7,   label = "middling",     colour = cli::col_green,  alert = cli::cli_alert_success, suitability = "probably"),
    list(min = .6,   label = "mediocre",     colour = cli::col_yellow, alert = cli::cli_alert_warning, suitability = "probably"),
    list(min = .5,   label = "miserable",    colour = cli::col_red,    alert = cli::cli_alert_danger,  suitability = "hardly"),
    list(min = -Inf, label = "unacceptable", colour = cli::col_red,    alert = cli::cli_alert_danger,  suitability = "not")
  )
  Find(function(b) kmo >= b$min, bands)
}

# Belsley's (1991) reading of the condition index: the highest band the index clears gives
# the alert, the verdict text, and `flag`, which gates both the section verdict and the
# recommendation. The bands are steps on the progression 1, 3, 10, 30, 100, 300, 1000 that
# Belsley reads relative strength from (p. 38). An index of 10 or less is rarely of interest
# (pp. 43, 46); 30 is his example value for a near dependency, which an index must exceed
# (p. 38); 30 to 100 is moderate (p. 42), and his worked example grades indexes of 35, 153,
# and 455 as moderately strong, strong, and very strong (pp. 44-45).
.ci_band <- function(ci) {
  bands <- list(
    list(above = 300,  alert = cli::cli_alert_danger,  flag = TRUE, strength = "very strong"),
    list(above = 100,  alert = cli::cli_alert_danger,  flag = TRUE, strength = "strong"),
    list(above = 30,   alert = cli::cli_alert_danger,  flag = TRUE, strength = "moderate"),
    list(above = 10,   alert = cli::cli_alert_warning, flag = FALSE,
         note = "The index is above 10, but not above 30, the value that flags a near linear dependency (Belsley, 1991)."),
    list(above = -Inf, alert = cli::cli_alert_success, flag = FALSE,
         note = "An index of 10 or less is rarely of interest (Belsley, 1991).")
  )
  # An index that is not available takes the lowest band.
  if (is.na(ci)) return(bands[[length(bands)]])
  Find(function(b) ci > b$above, bands)
}

# TRUE when the complete cases used by the normality and outlier diagnostics are fewer than
# the rows supplied.
.screen_incomplete <- function(n_complete, n_obs) {
  !is.null(n_obs) && !is.na(n_obs) && !is.null(n_complete) && n_obs > n_complete
}

# The base of a complete-case count: with missing values, both the complete cases and the
# rows supplied are named, so a count cannot be read against the wrong base.
.screen_case_base <- function(n_complete, n_obs) {
  if (.screen_incomplete(n_complete, n_obs)) {
    paste0(n_complete, " complete cases of the ", n_obs, " rows supplied")
  } else {
    paste0(n_complete, " observations")
  }
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
# otherwise. `text` carries the statistic and p-value (no cli braces) and is not wrapped, so
# a break cannot fall inside the statistic.
.screen_mvn_line <- function(sig, text) {
  if (sig) cli::cli_alert_danger(text) else cli::cli_alert_success(text)
}

# Oxford-comma list of names ("a", "a and b", "a, b, and c").
.screen_and_list <- function(x) {
  if (length(x) == 1L) return(x)
  if (length(x) == 2L) return(paste(x, collapse = " and "))
  paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
}

# The per-variable display table, with numeric columns rounded to `digits`: the full
# per-item table for raw data (MSA/SMC headers; `missing%` only when a value is missing;
# `flags` only when a variable was screened for categories), or MSA and SMC alone for a
# correlation-matrix input.
.screen_per_item_display <- function(x, digits) {
  if (!is.null(x$per_item)) {
    disp <- x$per_item
    if (all(disp$missing == 0)) disp$missing <- NULL
    num <- vapply(disp, is.numeric, logical(1))
    disp[num] <- lapply(disp[num], round, digits)
    names(disp)[names(disp) == "smc"] <- "SMC"
    names(disp)[names(disp) == "kmo_i"] <- "MSA"
    names(disp)[names(disp) == "missing"] <- "missing%"
    # NA flags mean "treated as continuous"; <NA> would read as a failed computation.
    if (all(is.na(disp$flags))) disp$flags <- NULL
    else disp$flags[is.na(disp$flags)] <- "-"
    disp
  } else {
    # make.unique(): duplicate dimnames are valid in a correlation matrix, but not as
    # data-frame row names.
    data.frame(MSA = round(x$kmo$KMO_i, digits),
               SMC = round(x$smc, digits),
               row.names = make.unique(names(x$kmo$KMO_i)))
  }
}

# The recommendations as parallel symbol/message vectors. Each check fires only when it
# applies; a single success line is returned when none fires, and a correlation-matrix input
# always closes with the raw-data note.
.screen_recommendations <- function(x, raw, kmo, mvn_nonnormal, multicollinear, digits) {
  sym <- character(0)
  msg <- character(0)
  push <- function(s, m) { sym[[length(sym) + 1L]] <<- s; msg[[length(msg) + 1L]] <<- m }

  # Items with fewer than five response categories, for the ordinal advice. Not relevant when
  # the correlations are already polychoric or tetrachoric.
  few <- character(0)
  if (raw && !is.null(x$categories) && !.is_poly_cor(x$settings$cor_method)) {
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

  # Bartlett's test is almost always significant at usual N, so only its failure is news.
  bart <- x$bartlett
  if (!is.null(bart)) {
    if (is.na(bart$p_value)) {
      push("!", paste("Bartlett's test of sphericity has no result: N is too small relative",
                      "to the number of variables for a factor analysis."))
    } else if (bart$p_value >= .05) {
      push("x", paste0("Bartlett's test of sphericity is not significant (p",
                       .screen_p_str(bart$p_value),
                       "); the variables may not correlate enough for a factor analysis."))
    }
  }

  # Variables with a low individual MSA (< .5) share little variance and are drop
  # candidates (Kaiser & Rice, 1974). make.unique() labels them as the per-item table does.
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

  # Non-normal data with 5 or more categories: normal-theory standard errors are biased. Not
  # given when the few-category advice already points to an ordinal estimator.
  if (isTRUE(mvn_nonnormal) && !length(few)) {
    push("!", paste("These data depart from multivariate normality; normal-theory standard errors",
                    "and fit statistics may be biased - prefer robust (sandwich) or bootstrapped",
                    "standard errors."))
  }

  # A condition index above 30 (Belsley, 1991), named by the variables involved.
  if (isTRUE(multicollinear)) {
    cv <- x$collinear
    push("!", if (length(cv)) {
      paste0("A condition index above 30 points to a near linear dependency among ",
             .screen_and_list(cv), "; consider removing one or more of them.")
    } else {
      paste("A condition index above 30 points to a near linear dependency among the",
            "variables; look for redundant items and consider removing some of them.")
    })
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
      # Many more flags than the cutoff admits: a high-breakdown fit to the most concentrated
      # subset of a clustered or mixture sample calls the rest distant, so the excess shows
      # that the data are not elliptical. The count must also be too long to inspect case by
      # case (at a tight cutoff or in a small sample a few genuine outliers clear the rate).
      if (nf >= 10L && frac > 5 * nominal) {
        # With missing data both bases are named: against the rows supplied the rate would
        # look milder than it is.
        of_those <- if (.screen_incomplete(x$outliers$n_complete, x$settings$n_obs)) {
          "% of those cases"
        } else {
          "%"
        }
        push("!", paste0(nf, " of ",
                         .screen_case_base(x$outliers$n_complete, x$settings$n_obs),
                         " (", round(100 * frac), of_those,
                         ") exceed the outlier cutoff, far above the ",
                         # signif() keeps a tight cutoff's nominal rate from rounding to "0%"
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

  if (!length(msg)) {
    push("v", "The data appear suitable for factor analysis.")
  }
  if (!raw) {
    push("i", x$note$message)
  }

  list(symbol = sym, message = msg)
}
