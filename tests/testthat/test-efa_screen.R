scr_cor <- efa_screen(test_models$baseline$cormat, N = 500)
# GRiPS has an exact collinearity (enjoy == commonly for most respondents), so the robust
# MCD covariance is (near-)singular and the outlier diagnostic falls back to classical
# distances with a warning; suppress it here (it is asserted in its own test below). A seed
# makes the MCD random subsets - and hence the outlier diagnostics and the print snapshot -
# reproducible; it does not affect the seed-independent KMO/Bartlett/per-item/normality parts.
scr_raw <- suppressWarnings(efa_screen(GRiPS_raw, seed = 1))
scr_iris <- efa_screen(iris[, 1:4], seed = 1)

dat_nonames <- test_models$baseline$cormat
colnames(dat_nonames) <- NULL
rownames(dat_nonames) <- NULL
scr_nona <- efa_screen(dat_nonames, N = 500)

# correlation matrix without N: Bartlett's test skipped
scr_non <- suppressWarnings(efa_screen(test_models$baseline$cormat))

test_that("output class and structure are correct", {
  expect_s3_class(scr_cor, "efa_screen")
  expect_s3_class(scr_raw, "efa_screen")

  expect_named(scr_cor$kmo, c("KMO", "KMO_i"))
  expect_named(scr_cor$bartlett, c("chisq", "p_value", "df"))

  expect_length(scr_cor$kmo$KMO_i, ncol(test_models$baseline$cormat))
  expect_length(scr_cor$smc, ncol(test_models$baseline$cormat))

  # raw-only sections are NULL when a correlation matrix is analysed, with a note
  expect_null(scr_cor$per_item)
  expect_null(scr_cor$normality)
  expect_null(scr_cor$outliers)
  expect_null(scr_cor$categories)
  expect_s3_class(scr_cor$note, "efa_screen_no_raw")

  # the raw path fills the per-item table and clears the note
  expect_s3_class(scr_raw$per_item, "data.frame")
  expect_named(scr_raw$per_item, c("variance", "missing", "smc", "kmo_i", "flags"))
  expect_null(scr_raw$note)
})

test_that("KMO, Bartlett, determinant, condition, and SMC match the standalone tools", {
  R <- test_models$baseline$cormat

  # KMO consistent with the standalone KMO()
  k <- KMO(R)
  expect_equal(scr_cor$kmo$KMO, k$KMO, tolerance = 1e-10)
  expect_equal(unname(scr_cor$kmo$KMO_i), unname(k$KMO_i), tolerance = 1e-10)

  # Bartlett consistent with the standalone BARTLETT()
  b <- BARTLETT(R, N = 500)
  expect_equal(scr_cor$bartlett$chisq, b$chisq, tolerance = 1e-8)
  expect_equal(scr_cor$bartlett$df, b$df)
  expect_equal(scr_cor$bartlett$p_value, b$p_value, tolerance = 1e-8)

  # raw-data branch: N recovered from the data
  expect_equal(scr_raw$bartlett$chisq, suppressMessages(BARTLETT(GRiPS_raw))$chisq,
               tolerance = 1e-8)
  expect_equal(scr_raw$settings$N, nrow(GRiPS_raw))

  # determinant and condition number (eigenvalue-ratio convention)
  expect_equal(scr_cor$determinant, det(R), tolerance = 1e-8)
  expect_equal(scr_cor$condition, kappa(R, exact = TRUE), tolerance = 1e-8)
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(scr_cor$condition, max(ev) / min(ev), tolerance = 1e-6)

  # SMC == the closed form 1 - 1/diag(R^-1)
  expect_equal(unname(scr_cor$smc), unname(1 - 1 / diag(solve(R))), tolerance = 1e-10)
  expect_length(scr_cor$smc, ncol(R))
})

test_that("KMO, Bartlett, and SMC match psych", {
  skip_on_cran()
  skip_if_not_installed("psych")

  R <- test_models$baseline$cormat

  ps <- psych::KMO(R)
  expect_equal(scr_cor$kmo$KMO, unname(ps$MSA), tolerance = 1e-4)
  expect_equal(unname(scr_cor$kmo$KMO_i), unname(ps$MSAi), tolerance = 1e-4)

  expect_equal(scr_cor$bartlett$chisq, psych::cortest.bartlett(R, n = 500)$chisq,
               tolerance = 1e-4)

  expect_equal(unname(scr_cor$smc), unname(psych::smc(R)), tolerance = 1e-6)
})

test_that("the raw-data path is silent and matches the standalone tools", {
  # inform_from_data = FALSE: the "computing correlations from the raw data" note
  # is suppressed here (unlike the standalone suitability tools)
  expect_no_message(suppressWarnings(efa_screen(GRiPS_raw)), class = "efa_cor_from_data")

  # the correlation-side diagnostics on the raw path go through the same pipeline
  # as the standalone KMO(), including the per-variable naming
  km <- suppressMessages(KMO(GRiPS_raw))
  expect_equal(scr_raw$kmo$KMO, km$KMO, tolerance = 1e-10)
  expect_equal(unname(scr_raw$kmo$KMO_i), unname(km$KMO_i), tolerance = 1e-10)
  expect_equal(names(scr_raw$kmo$KMO_i), colnames(GRiPS_raw))
  expect_equal(names(scr_raw$smc), colnames(GRiPS_raw))
})

test_that("per-variable diagnostics get a V-fallback when unnamed", {
  expect_equal(names(scr_nona$kmo$KMO_i), paste0("V", seq_len(ncol(dat_nonames))))
  expect_equal(names(scr_nona$smc), paste0("V", seq_len(ncol(dat_nonames))))
})

test_that("Bartlett's test is skipped with a note when N is unavailable", {
  expect_warning(efa_screen(test_models$baseline$cormat),
                 class = "efa_suitability_no_n")

  expect_null(scr_non$bartlett)
  expect_true(is.na(scr_non$settings$N))
  # the remaining diagnostics are still returned
  expect_false(is.null(scr_non$kmo))
  expect_false(is.null(scr_non$smc))
  expect_true(is.finite(scr_non$determinant))
  expect_true(is.finite(scr_non$condition))
})

test_that("a too-small N is reported rather than left as a bare NA", {
  # N is known here, so the test is attempted and comes back undefined: the multiplier
  # N - 1 - (2p + 5)/6 is non-positive at N = 5 for p = 18. Same condition and class as
  # efa_bartlett(), so a user reads the same reason from either entry point.
  expect_warning(scr_small <- efa_screen(test_models$baseline$cormat, N = 5),
                 class = "efa_bartlett_n_too_small")
  expect_true(is.na(scr_small$bartlett$chisq))
  expect_true(is.na(scr_small$bartlett$p_value))
  # the remaining suitability diagnostics are unaffected
  expect_true(is.finite(scr_small$kmo$KMO))
  expect_true(is.finite(scr_small$determinant))

  # the printed report carries the same reason, so it is not only in the warning. Matched on
  # the collapsed report rather than line by line, so the assertion does not depend on where
  # cli happens to break the sentence.
  expect_match(paste(cli::ansi_strip(format(scr_small)), collapse = " "),
               "Bartlett multiplier", fixed = TRUE)

  expect_no_warning(efa_screen(test_models$baseline$cormat, N = 500),
                    class = "efa_bartlett_n_too_small")
})

test_that("settings are returned correctly", {
  expect_named(scr_cor$settings, c("N", "n_obs", "use", "cor_method", "mcd_alpha",
                                   "outlier_cutoff", "seed"))
  expect_equal(scr_cor$settings$N, 500)
  expect_equal(scr_cor$settings$use, "pairwise.complete.obs")
  expect_equal(scr_cor$settings$cor_method, "pearson")
  expect_equal(scr_cor$settings$outlier_cutoff, 0.975)

  # n_obs is the number of supplied raw rows (NA for a correlation-matrix input)
  expect_true(is.na(scr_cor$settings$n_obs))
  expect_equal(scr_raw$settings$n_obs, nrow(GRiPS_raw))

  # the MCD arguments that govern the outlier diagnostics are recorded, so a saved
  # object says how its robust distances were obtained
  expect_equal(scr_cor$settings$mcd_alpha, 0.5)
  expect_null(scr_cor$settings$seed)
  # on the raw path both are the values the outlier diagnostics actually ran with
  expect_equal(scr_raw$settings$seed, 1)
  expect_equal(scr_raw$settings$mcd_alpha, 0.5)
  expect_equal(
    efa_screen(test_models$baseline$cormat, N = 500, mcd_alpha = 0.75)$settings$mcd_alpha,
    0.75)
})

test_that("raw-data per-item diagnostics: variance, missing, and category flags (GRiPS)", {
  pi <- scr_raw$per_item

  # per-item variance, and (zero) missingness for the complete GRiPS data
  expect_equal(pi$variance, apply(GRiPS_raw, 2, var),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_true(all(pi$missing == 0))

  # smc / kmo_i in the table are the same vectors reported at the top level
  expect_equal(pi$smc, unname(scr_raw$smc), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(pi$kmo_i, unname(scr_raw$kmo$KMO_i), tolerance = 1e-12,
               ignore_attr = TRUE)

  # GRiPS items are 6-category (0-5) with no gaps; enjoy/hurt/attracted each have a
  # category with fewer than five responses (sparse), the other five are clean
  expect_equal(pi[c("enjoy", "hurt", "attracted"), "flags"], rep("sparse", 3L),
               ignore_attr = TRUE)
  expect_equal(pi[c("fun", "friends", "part", "commonly", "chances"), "flags"],
               rep("", 5L), ignore_attr = TRUE)

  # categories: one count table per item, summing to N, in category order
  expect_type(scr_raw$categories, "list")
  expect_length(scr_raw$categories, ncol(GRiPS_raw))
  expect_true(all(vapply(scr_raw$categories, sum, numeric(1)) == nrow(GRiPS_raw)))
  expect_equal(names(scr_raw$categories$fun), as.character(0:5))
})

test_that("raw-data flags separate sparse, empty, and continuous variables", {
  set.seed(123)
  n <- 60L
  cont1 <- rnorm(n)
  cont2 <- rnorm(n)
  miss <- rnorm(n); miss[sample.int(n, 6L)] <- NA           # 10% missing, continuous
  sparse <- sample(c(rep(1, 20), rep(2, 20), rep(3, 17), rep(4, 3)))  # 4 appears 3x
  empty <- sample(c(rep(1, 15), rep(2, 15), rep(4, 15), rep(5, 15)))  # 3 is skipped
  dat <- cbind(cont1, cont2, miss, sparse, empty)

  scr <- efa_screen(dat)
  pi <- scr$per_item

  # variance and percentage missing computed on all supplied rows, independent of use
  expect_equal(pi$variance, apply(dat, 2, var, na.rm = TRUE),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(pi["miss", "missing"], 10)
  expect_equal(pi[c("cont1", "cont2", "sparse", "empty"), "missing"], rep(0, 4L),
               ignore_attr = TRUE)
  expect_equal(scr$settings$n_obs, n)

  # a low-frequency category is flagged sparse, a skipped interior category empty,
  # and variables with many distinct values are treated as continuous (flags NA)
  expect_equal(pi["sparse", "flags"], "sparse")
  expect_equal(pi["empty", "flags"], "empty")
  expect_true(all(is.na(pi[c("cont1", "cont2", "miss"), "flags"])))

  # the offending categories are visible in the tabulation; continuous columns are NA
  expect_equal(unname(scr$categories$sparse["4"]), 3L)
  expect_false("3" %in% names(scr$categories$empty))
  expect_true(is.na(scr$categories$cont1))
})

test_that("factor columns are coded (not corrupted) and an unused ordinal category is empty", {
  set.seed(1)
  n <- 120L
  # ordered factors carry the response order the polychoric path needs; f1 never uses its
  # third level ("c"), leaving an interior gap
  f1 <- ordered(sample(c("a", "b", "d"), n, replace = TRUE),
                levels = c("a", "b", "c", "d"))
  f2 <- ordered(sample(c("lo", "mid", "hi"), n, replace = TRUE),
                levels = c("lo", "mid", "hi"))
  f3 <- ordered(sample(c("lo", "mid", "hi"), n, replace = TRUE),
                levels = c("lo", "mid", "hi"))
  dat <- data.frame(f1 = f1, f2 = f2, f3 = f3)

  scr <- suppressWarnings(efa_screen(dat, cor_method = "poly"))

  # data.matrix() coding (not as.matrix()): variance is finite and matches the codes
  expect_true(all(is.finite(scr$per_item$variance)))
  expect_equal(scr$per_item$variance, apply(data.matrix(dat), 2, var),
               tolerance = 1e-10, ignore_attr = TRUE)
  # f1's unused middle level (c -> code 3, between b = 2 and d = 4) is an empty category
  expect_match(scr$per_item["f1", "flags"], "empty")

  # the category counts are labelled by the original levels, not by the integer codes,
  # and only the observed levels appear (f1 never uses "c")
  expect_equal(names(scr$categories$f1), c("a", "b", "d"))
  expect_equal(names(scr$categories$f2), c("lo", "mid", "hi"))
  expect_equal(scr$categories$f1, table(dat$f1)[c("a", "b", "d")], ignore_attr = TRUE)
})

test_that("raw data with duplicate column names is handled positionally, not by name", {
  set.seed(7)
  x <- data.frame(a = rnorm(150), a = rnorm(150, 5, 3), b = rnorm(150),
                  check.names = FALSE)
  scr <- efa_screen(x)

  # no crash; per-item rows carry make.unique()'d variable names
  expect_equal(rownames(scr$per_item), c("a", "a.1", "b"))
  # each column reports its OWN diagnostics (the duplicate name is not aliased to
  # the first matching column)
  expect_equal(scr$per_item$variance, apply(data.matrix(x), 2, var),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_length(scr$categories, ncol(x))
  # per_item and categories are keyed the same way, so named lookup agrees
  expect_equal(names(scr$categories), rownames(scr$per_item))
})

test_that("a wide integer-code span does not overflow the empty-category test", {
  set.seed(9)
  # an integer MATRIX (data.matrix keeps integer storage) whose ordinal column spans
  # more than the 32-bit integer range: max - min must be formed in double, or it
  # overflows to NA and mislabels the flag "NA" instead of "empty"
  mi <- matrix(0L, 200L, 3L)
  mi[, 1] <- sample(c(-2000000000L, 0L, 1L, 2000000000L), 200L, replace = TRUE)
  mi[, 2] <- sample(1:5, 200L, replace = TRUE)
  mi[, 3] <- sample(1:5, 200L, replace = TRUE)
  colnames(mi) <- c("v1", "v2", "v3")
  expect_identical(typeof(mi), "integer")

  scr <- suppressWarnings(efa_screen(mi))
  expect_equal(scr$per_item["v1", "flags"], "empty")
})

test_that("large category codes are labelled as plain values, not scientific notation", {
  set.seed(10)
  dat <- data.frame(big = sample(c(1, 2, 100000000), 80L, replace = TRUE),
                    v2 = rnorm(80), v3 = rnorm(80))
  scr <- efa_screen(dat)
  expect_true("100000000" %in% names(scr$categories$big))
  expect_false("1e+08" %in% names(scr$categories$big))
})

test_that("a non-numeric (character) matrix is rejected on the polychoric path", {
  set.seed(11)
  cm <- matrix(sample(c("lo", "mid", "hi"), 300L, replace = TRUE), nrow = 100L,
               dimnames = list(NULL, c("v1", "v2", "v3")))
  # A character matrix carries no response order: coding it would rank the labels
  # alphabetically (hi < lo < mid), silently attenuating the polychoric correlations, so the
  # polychoric path refuses it rather than returning a misleading screening report.
  expect_error(suppressMessages(efa_screen(cm, cor_method = "poly")),
               class = "efa_cor_unordered_factor")
})

test_that("a wide-span integer code is flagged empty without enumerating its range", {
  set.seed(8)
  # a sentinel code far from the ordinal range: the empty-category test must not
  # enumerate the full 1..2e9 span (which would exhaust memory)
  dat <- data.frame(v1 = sample(c(1L, 2L, 3L, 2000000000L), 200, replace = TRUE),
                    v2 = rnorm(200), v3 = rnorm(200))
  scr <- efa_screen(dat)
  expect_equal(scr$per_item["v1", "flags"], "empty")
})

test_that("multivariate normality diagnostics are present on the raw path", {
  nm <- scr_raw$normality
  expect_named(nm, c("mardia", "hz", "n_complete"))
  expect_named(nm$mardia, c("skewness", "skewness_df", "skewness_p",
                            "kurtosis", "kurtosis_p", "b1p", "b2p"))
  expect_named(nm$hz, c("statistic", "p_value"))

  # GRiPS is complete, so every row feeds the tests
  expect_equal(nm$n_complete, nrow(GRiPS_raw))

  # statistics finite and p-values in [0, 1]
  expect_true(all(is.finite(c(nm$mardia$skewness, nm$mardia$kurtosis, nm$hz$statistic))))
  expect_true(nm$mardia$skewness_p >= 0 && nm$mardia$skewness_p <= 1)
  expect_true(nm$mardia$kurtosis_p >= 0 && nm$mardia$kurtosis_p <= 1)
  expect_true(nm$hz$p_value >= 0 && nm$hz$p_value <= 1)

  # NULL for a correlation-matrix input
  expect_null(scr_cor$normality)
})

test_that("the Henze-Zirkler p-value is withheld where its null approximation degenerates", {
  # The lognormal null variance si2 falls geometrically with the number of variables
  # while the null mean mu goes to 1, so si2 + mu^2 rounds to mu^2 and the z-score
  # divides by a zero psi. Without the guard the p-value is then exactly 0 or exactly 1,
  # decided by a rounding-level residual, on data from an exact multivariate normal.
  set.seed(4)
  X <- matrix(stats::rnorm(150 * 60), 150, 60)

  expect_warning(nm <- EFAtools:::.screen_normality(X), class = "efa_screen_no_hz")
  expect_s3_class(nm$hz, "efa_screen_no_hz")
  expect_true(is.na(nm$hz$p_value))
  # only the null approximation fails: the statistic and Mardia's tests are unaffected
  expect_true(is.finite(nm$hz$statistic))
  expect_true(is.finite(nm$mardia$skewness_p))
  expect_true(is.finite(nm$mardia$kurtosis_p))

  # and the report neither passes nor rejects on it, and counts the tests that remain
  local_reproducible_output()
  colnames(X) <- paste0("v", seq_len(60))
  expect_warning(scr <- efa_screen(X, seed = 1), class = "efa_screen_no_hz")
  txt <- cli::ansi_strip(paste(format(scr), collapse = " "))
  expect_match(txt, "p not available at this number of variables", fixed = TRUE)
  expect_match(txt, "of the 2 tests", fixed = TRUE)
})

test_that("the multivariate normality verdict counts the tests that reject", {
  local_reproducible_output()
  # GRiPS departs from normality on all three tests, so the verdict names all three; the
  # count is what separates a unanimous rejection from a single marginal one.
  expect_match(cli::ansi_strip(paste(format(scr_raw), collapse = " ")),
               "3 of the 3 tests reject it", fixed = TRUE)
  # a clean multivariate normal sample takes the other branch, which names the
  # denominator too (all three p-values here are above .96, so the fixture is not tied)
  set.seed(381)
  clean <- matrix(stats::rnorm(200 * 5), 200, 5)
  colnames(clean) <- paste0("v", seq_len(5))
  expect_match(cli::ansi_strip(paste(format(efa_screen(clean, seed = 1)), collapse = " ")),
               "consistent with multivariate normality: none of the 3 tests rejects it",
               fixed = TRUE)
})

test_that("Mardia's skewness/kurtosis match psych after the divisor rescaling", {
  skip_on_cran()
  skip_if_not_installed("psych")

  X <- as.matrix(iris[, 1:4])
  n <- nrow(X)
  md <- scr_iris$normality$mardia
  m <- psych::mardia(X, plot = FALSE)

  # efa_screen uses the biased (divisor-n) covariance; psych::mardia uses cov()
  # (divisor n - 1). The two are related exactly by (n / (n - 1))^j, so the COEFFICIENTS
  # agree after rescaling. The two statistics built on them do not, by design: psych
  # standardises b2p with the asymptotic null moments where efa_screen uses the exact
  # ones, and psych drops Mardia's skewness correction above 20 observations where
  # efa_screen keeps it. Both divergences have their own tests below.
  expect_equal(md$b1p, m$b1p * (n / (n - 1))^3, tolerance = 1e-6)
  expect_equal(md$b2p, m$b2p * (n / (n - 1))^2, tolerance = 1e-6)
})

test_that("Mardia's kurtosis is standardised with Mardia's (1974) exact null moments", {
  md <- scr_iris$normality$mardia
  n <- nrow(iris)
  p <- 4L

  kurt_mean <- p * (p + 2) * (n - 1) / (n + 1)
  kurt_var <- 8 * p * (p + 2) * (n - 3) * (n - p - 1) * (n - p + 1) /
    ((n + 1)^2 * (n + 3) * (n + 5))
  expect_equal(md$kurtosis, (md$b2p - kurt_mean) / sqrt(kurt_var), tolerance = 1e-10)
})

test_that("Mardia's kurtosis holds its nominal level on exact multivariate normal data", {
  # The exact null moments matter most when the number of variables is large relative to
  # the sample size. The asymptotic pair p(p + 2) and 8p(p + 2)/n that most other
  # implementations use overstates both: at n = 50 and p = 20 it puts the statistic at
  # about -2.1 on average with a standard deviation near 0.54, and rejects data from an
  # exact multivariate normal about 59% of the time at the .05 level. Correcting the mean
  # alone leaves the standard deviation at 0.54 and turns that into a test that almost
  # never rejects, so all three checks below are needed to pin both corrections.
  set.seed(20240817)
  n <- 50L
  p <- 20L
  md <- lapply(seq_len(300L), function(i) {
    EFAtools:::.screen_normality(matrix(stats::rnorm(n * p), n, p))$mardia
  })
  z <- vapply(md, `[[`, numeric(1), "kurtosis")
  pv <- vapply(md, `[[`, numeric(1), "kurtosis_p")

  expect_lt(abs(mean(z)), 0.25)
  expect_gt(stats::sd(z), 0.85)
  expect_lt(stats::sd(z), 1.20)
  # the reported p-value drives the level, so assert it rather than rebuilding it
  expect_lt(mean(pv < .05), 0.12)
  # and it is the two-sided normal tail of the reported statistic
  expect_equal(pv, 2 * stats::pnorm(abs(z), lower.tail = FALSE), tolerance = 1e-12)
})

test_that("Mardia's kurtosis is NA when its exact null variance is zero", {
  # At n = p + 1 the centred data span the whole space orthogonal to the mean, so b2p
  # equals (n - 1)^2 with probability 1 and carries no information; its exact null
  # variance is 0 and the standardised statistic does not exist.
  set.seed(1)
  n <- 8L
  p <- 7L
  nrm <- EFAtools:::.screen_normality(matrix(stats::rnorm(n * p), n, p))

  expect_equal(nrm$mardia$b2p, (n - 1)^2, tolerance = 1e-8)
  expect_true(is.na(nrm$mardia$kurtosis))
  expect_true(is.na(nrm$mardia$kurtosis_p))
  # the skewness test is unaffected
  expect_true(is.finite(nrm$mardia$skewness_p))

  # the printed report calls it unavailable rather than showing it as a pass, and drops
  # it from the denominator of the verdict count
  set.seed(4)
  Xn <- matrix(stats::rnorm(12 * 11), 12, 11)
  colnames(Xn) <- paste0("v", seq_len(11))
  txt <- cli::ansi_strip(paste(suppressWarnings(format(efa_screen(Xn, seed = 1))),
                               collapse = " "))
  expect_match(txt, "Mardia's kurtosis: not available", fixed = TRUE)
  expect_match(txt, "of the 2 tests", fixed = TRUE)
})

test_that("the Henze-Zirkler statistic and p-value match the reference value on iris", {
  # Henze-Zirkler (1990) statistic for the four numeric iris columns, computed once from
  # this implementation and cross-checked against MVN::hz(iris[, 1:4]) = 2.336.
  expect_equal(scr_iris$normality$hz$statistic, 2.3363942003, tolerance = 1e-6)
  # p-value pinned relatively (an absolute tolerance is meaningless at ~1e-19)
  expect_equal(scr_iris$normality$hz$p_value / 4.1413116299e-19, 1, tolerance = 1e-6)
})

test_that("the corrected Mardia skewness statistic matches psych below 20 observations", {
  skip_on_cran()
  skip_if_not_installed("psych")

  set.seed(1)
  Xs <- matrix(rnorm(15 * 3), 15, 3)
  nrm <- EFAtools:::.screen_normality(Xs)
  ms <- psych::mardia(Xs, plot = FALSE)

  # below 20 observations psych also applies the correction, as `small.skew`, so the two
  # agree after the divisor rescaling
  expect_equal(nrm$mardia$skewness, ms$small.skew * (15 / 14)^3, tolerance = 1e-6)
})

test_that("Mardia's skewness correction applies at every sample size", {
  # Mardia's (1974) correction k, eq. (5.5), is what makes the statistic's expectation
  # equal its degrees of freedom; it is not a small-sample device. psych and most other
  # implementations apply it only below 20 observations, so above that the two diverge by
  # exactly k. Here n = 60, far above any small-sample threshold.
  set.seed(4)
  n <- 60L
  p <- 6L
  X <- matrix(stats::rnorm(n * p), n, p)
  nrm <- EFAtools:::.screen_normality(X)
  m <- psych::mardia(X, plot = FALSE)

  k <- (p + 1) * (n + 1) * (n + 3) / (n * ((n + 1) * (p + 1) - 6))
  expect_gt(k, 1)
  # psych's uncorrected `skew` is our statistic divided by k, after the divisor rescaling
  expect_equal(nrm$mardia$skewness, k * m$skew * (n / (n - 1))^3, tolerance = 1e-6)
  # and matches psych's own corrected form, which psych computes but does not use here
  expect_equal(nrm$mardia$skewness, m$small.skew * (n / (n - 1))^3, tolerance = 1e-6)
})

test_that("the corrected Mardia skewness statistic is centred on its degrees of freedom", {
  # The correction makes E(statistic) = df exactly at every n. Without it the statistic is
  # biased low by the factor k: at n = 50 and p = 20 its expectation is df / k = 1444.7
  # against df = 1540, which is 1.72 chi-square standard deviations below, and the test
  # then almost never rejects.
  set.seed(20240818)
  n <- 50L
  p <- 20L
  df <- p * (p + 1) * (p + 2) / 6
  s <- vapply(seq_len(200L), function(i) {
    EFAtools:::.screen_normality(matrix(stats::rnorm(n * p), n, p))$mardia$skewness
  }, numeric(1))

  # mean within 0.5 chi-square sd of the df; the uncorrected statistic sits 1.72 sd below
  expect_lt(abs(mean(s) - df) / sqrt(2 * df), 0.5)
})

test_that("normality tests abort on a correlation matrix and on a singular covariance", {
  expect_error(EFAtools:::.screen_normality(test_models$baseline$cormat),
               class = "efa_screen_normality_cormat")

  set.seed(1)
  # n_complete <= p: the complete-case covariance is singular (centring leaves rank at
  # most n - 1). Both n < p and the n == p boundary abort deterministically - at n == p
  # a rounding-positive pivot can let chol() alone succeed on the singular covariance.
  expect_error(EFAtools:::.screen_normality(matrix(rnorm(3 * 5), 3, 5)),
               class = "efa_screen_singular")
  expect_error(EFAtools:::.screen_normality(matrix(rnorm(4 * 4), 4, 4)),
               class = "efa_screen_singular")
})

test_that("efa_screen skips the normality tests when the covariance is singular", {
  set.seed(3)
  Sig <- matrix(0.35, 4, 4)
  diag(Sig) <- 1
  # Drawn through the Cholesky factor: for a positive definite Sigma that factor is unique
  # (the positive-diagonal convention fixes it), so the seeded sample is reproducible to
  # rounding on any LAPACK build, which an eigen-based draw is not -- eigenvector signs and
  # the basis of a repeated eigenspace are left free by the decomposition.
  dat <- matrix(stats::rnorm(200 * 4), 200) %*% chol(Sig)
  colnames(dat) <- paste0("v", seq_len(4))
  # a rotating block of missing values leaves no complete case (so the MVN tests cannot
  # run) while every pair of variables stays well observed (so the pairwise correlation
  # matrix, and the other diagnostics, are unaffected)
  blk <- rep(1:4, each = 50)
  for (i in seq_len(200)) dat[i, blk[i]] <- NA

  # with no complete case neither the normality tests nor the robust outlier diagnostics
  # can run, so each degrades with its own classed warning and note
  expect_warning(
    expect_warning(scr <- efa_screen(dat), class = "efa_screen_no_mvn"),
    class = "efa_screen_no_outliers"
  )
  expect_s3_class(scr$normality, "efa_screen_no_mvn")
  expect_s3_class(scr$outliers, "efa_screen_no_outliers")
  # the remaining diagnostics are still returned
  expect_false(is.null(scr$kmo))
  expect_s3_class(scr$per_item, "data.frame")
})

# ---- robust-Mahalanobis outlier diagnostics (FAST-MCD) ----

# injected-outlier fixture: a clean multivariate-normal core with a block of far
# outliers; a robust covariance should recover the clean centre/scatter and flag the
# injected rows, where the classical covariance (inflated by the outliers) would not
set.seed(2024)
n_out <- 200L
p_out <- 4L
Sig_out <- matrix(0.4, p_out, p_out)
diag(Sig_out) <- 1
X_out <- matrix(rnorm(n_out * p_out), n_out, p_out) %*% chol(Sig_out)
colnames(X_out) <- paste0("v", seq_len(p_out))
inj <- 1:12
# a shifted, jittered block of outliers (jittered so they do not collapse onto a single
# point, which would make their own covariance degenerate)
X_out[inj, ] <- matrix(rep(c(7, -7, 7, -7), each = length(inj)), ncol = p_out) +
  matrix(rnorm(length(inj) * p_out), ncol = p_out)
scr_out <- efa_screen(X_out, seed = 1)

test_that("robust outlier diagnostics have the expected structure", {
  o <- scr_out$outliers
  expect_named(o, c("distances", "cutoff", "flagged", "center", "cov", "method",
                    "fallback_reason", "n_complete"))
  expect_equal(o$method, "mcd")
  # nothing to explain when the robust estimate was available
  expect_null(o$fallback_reason)
  expect_equal(o$n_complete, n_out)
  expect_length(o$distances, n_out)
  # cutoff is on the distance scale (comparable to `distances`), and `flagged` is exactly
  # reproducible from the two reported fields
  expect_equal(o$cutoff, sqrt(qchisq(0.975, p_out)))
  expect_equal(sort(as.integer(names(o$distances)[o$distances > o$cutoff])), o$flagged)
  expect_true(all(is.finite(o$distances)))
  # the cutoff probability is exposed
  expect_equal(efa_screen(X_out, seed = 1, outlier_cutoff = 0.99)$outliers$cutoff,
               sqrt(qchisq(0.99, p_out)))
  # NULL on a correlation-matrix input
  expect_null(scr_cor$outliers)
})

test_that("the MCD flags the injected outliers and recovers the clean centre", {
  o <- scr_out$outliers
  # every injected row is flagged, and the injected rows are cleanly separated from the
  # rest (their robust distances exceed every non-injected distance)
  expect_true(all(inj %in% o$flagged))
  expect_gt(min(o$distances[inj]), max(o$distances[-inj]))
  # few false positives (near the 2.5% nominal rate), not a swamped solution
  expect_lte(length(setdiff(o$flagged, inj)), 15L)
  # the robust centre recovers the clean mean (~0) despite the outliers
  expect_equal(unname(o$center), rep(0, p_out), tolerance = 0.2)
  # the consistency + small-sample scaling makes the robust covariance unbiased for the
  # clean Sigma (unit variances); a scale-invariant correlation comparison cannot see this,
  # so a dropped consistency factor (which halves the scatter) would be caught here
  expect_equal(diag(o$cov), rep(1, p_out), tolerance = 0.35, ignore_attr = TRUE)
  # flagged indices refer to rows of the supplied data
  expect_true(all(o$flagged %in% seq_len(n_out)))
})

test_that("robust outlier diagnostics are reproducible via a seed and leave the RNG alone", {
  a <- efa_screen(X_out, seed = 7)
  b <- efa_screen(X_out, seed = 7)
  expect_identical(a$outliers$flagged, b$outliers$flagged)
  expect_equal(a$outliers$distances, b$outliers$distances)
  # a different seed can give a different (still-valid) solution but still flags injected
  expect_true(all(inj %in% efa_screen(X_out, seed = 99)$outliers$flagged))
  # the caller's random-number-generator stream is preserved
  set.seed(4321)
  before <- .Random.seed
  invisible(efa_screen(X_out, seed = 3))
  expect_identical(before, .Random.seed)
})

test_that("outlier row indices map back to the supplied rows when data are incomplete", {
  # incomplete rows are placed *before* the outliers, so a complete-case position differs
  # from the original row index; the reported flags and distance names must refer to the
  # supplied data's rows, not to compacted positions in the complete-case submatrix
  set.seed(77)
  n2 <- 220L
  S2 <- matrix(0.4, p_out, p_out)
  diag(S2) <- 1
  X2 <- matrix(rnorm(n2 * p_out), n2, p_out) %*% chol(S2)
  colnames(X2) <- paste0("v", seq_len(p_out))
  out2 <- 200:211
  X2[out2, ] <- matrix(rep(c(7, -7, 7, -7), each = length(out2)), ncol = p_out) +
    matrix(rnorm(length(out2) * p_out), ncol = p_out)
  X2[1:20, 1] <- NA                          # 20 incomplete rows before the outliers
  o <- efa_screen(X2, seed = 1)$outliers

  expect_equal(o$method, "mcd")
  expect_equal(o$n_complete, 200L)
  expect_true(all(out2 %in% o$flagged))
  # the outlier rows keep their supplied-data indices (>= 200); a missing remap would
  # report compacted positions <= 200
  expect_gte(max(o$flagged), 200L)
  expect_true("200" %in% names(o$distances))
  expect_false("1" %in% names(o$distances))  # row 1 is incomplete, so it has no distance
})

test_that("a non-default mcd_alpha threads through the whole outlier pipeline", {
  # a larger coverage (alpha = 0.75) enlarges the MCD subset but still yields a robust
  # estimate that flags the injected outliers
  a <- efa_screen(X_out, seed = 1, mcd_alpha = 0.75)$outliers
  expect_equal(a$method, "mcd")
  expect_true(all(inj %in% a$flagged))
  # the subset size follows the coverage: larger alpha -> larger h
  expect_gt(EFAtools:::.mcd_hsize(n_out, p_out, 0.75),
            EFAtools:::.mcd_hsize(n_out, p_out, 0.5))
})

test_that(".fast_mcd returns a genuine concentration fixed point, never a stopped-early subset", {
  # binary items make many h-subsets (near-)degenerate; the FAST-MCD search must return a
  # stable concentration fixed point (the h observations closest to its centre reproduce its
  # covariance) or abort to the classical fallback. It must never return an unconcentrated
  # (p + 1)-point or stopped-early over-tight scatter, which would silently over-flag.
  ld <- function(S) as.numeric(determinant(S, logarithm = TRUE)$modulus)
  set.seed(4)
  Sig <- 0.5 * matrix(1, 8L, 8L) + diag(0.5, 8L)
  ok <- replicate(20L, {
    B <- (matrix(rnorm(50L * 8L), 50L, 8L) %*% chol(Sig) > 0) + 0    # binary items
    B <- B[, apply(B, 2L, function(col) length(unique(col)) > 1L), drop = FALSE]
    h <- EFAtools:::.mcd_hsize(nrow(B), ncol(B), 0.5)
    fit <- tryCatch(EFAtools:::.fast_mcd(B, h), error = function(e) NULL)
    if (is.null(fit)) return(TRUE)                      # aborted -> classical fallback: fine
    d2 <- stats::mahalanobis(B, fit$center, fit$cov)
    idx <- order(d2)[seq_len(h)]
    abs(ld(stats::cov(B[idx, , drop = FALSE])) - ld(fit$cov)) < 0.1
  })
  expect_true(all(ok))
})

test_that("the robust covariance cross-checks against MASS::cov.rob", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  h <- floor((n_out + p_out + 1) / 2)
  # MASS::cov.rob(method = "mcd") returns a *reweighted* MCD estimate (with its own
  # consistency scaling), so the comparable quantity is efa_screen's reweighted robust
  # location/scatter; they agree closely on this well-separated fixture
  set.seed(123)
  mm <- MASS::cov.rob(X_out, method = "mcd", quantile.used = h, nsamp = "best")
  o <- scr_out$outliers
  expect_equal(unname(o$center), unname(mm$center), tolerance = 0.1)
  expect_equal(cov2cor(o$cov), cov2cor(mm$cov), tolerance = 0.05, ignore_attr = TRUE)

  # the raw MCD subset location (before the finite-sample correction) is a robust centre
  # in the same region
  set.seed(1)
  fm <- EFAtools:::.fast_mcd(X_out, h, nsamp = 500L)
  expect_equal(unname(fm$center), unname(mm$center), tolerance = 0.25)
})

test_that("the MCD subset size, consistency, and small-sample factors match their forms", {
  # subset size: floor((n + p + 1) / 2) at alpha = 0.5 (max breakdown), rising to n at 1
  expect_equal(EFAtools:::.mcd_hsize(200, 4, 0.5), floor((200 + 4 + 1) / 2))
  expect_equal(EFAtools:::.mcd_hsize(200, 4, 1), 200)
  # consistency factor: alpha / F_{chi^2_{p + 2}}(qchisq(alpha, p))
  expect_equal(EFAtools:::.mcd_consistency(4, 0.5),
               0.5 / pchisq(qchisq(0.5, 4), 4 + 2), tolerance = 1e-12)
  # the small-sample correction inflates the scatter at small n and tends to 1 as n grows
  expect_gt(EFAtools:::.mcd_cnp2(4, 40, 0.5, reweighted = FALSE), 1)
  expect_equal(EFAtools:::.mcd_cnp2(4, 1e6, 0.5, reweighted = FALSE), 1, tolerance = 1e-3)
  # regression pins for both correction tables (raw and reweighted, all three p branches),
  # locking the Pison et al. (2002) coefficients used to scale the robust covariance
  expect_equal(EFAtools:::.mcd_cnp2(4, 40, 0.5, reweighted = FALSE), 1.26016971771, tolerance = 1e-8)
  expect_equal(EFAtools:::.mcd_cnp2(2, 50, 0.5, reweighted = FALSE), 1.15095829889, tolerance = 1e-8)
  expect_equal(EFAtools:::.mcd_cnp2(1, 50, 0.5, reweighted = FALSE), 1.13894209891, tolerance = 1e-8)
  expect_equal(EFAtools:::.mcd_cnp2(4, 40, 0.5, reweighted = TRUE),  1.06512515802, tolerance = 1e-8)
  expect_equal(EFAtools:::.mcd_cnp2(2, 50, 0.5, reweighted = TRUE),  1.01272751567, tolerance = 1e-8)
  expect_equal(EFAtools:::.mcd_cnp2(1, 50, 0.5, reweighted = TRUE),  1.00806218791, tolerance = 1e-8)
})

test_that("outlier diagnostics fall back to classical distances and degrade gracefully", {
  # n <= 2p: the robust covariance is undefined, so classical Mahalanobis is used
  set.seed(5)
  Xsmall <- matrix(rnorm(8 * 4), 8, 4)
  expect_warning(
    o <- EFAtools:::.screen_outliers(Xsmall, mcd_alpha = 0.5, outlier_cutoff = 0.975,
                                     seed = 1),
    class = "efa_screen_mcd_fallback"
  )
  expect_equal(o$method, "classical")
  expect_length(o$distances, 8L)
  # the recorded reason names the branch actually taken: n <= 2p is decided before
  # .fast_mcd() runs, so this case is "too few complete cases", never an exact fit
  expect_match(o$fallback_reason, "too few complete cases")

  # collinear variables with too few cases: even the classical covariance is singular,
  # so the diagnostic is skipped with a classed note
  set.seed(6)
  Xd <- matrix(rnorm(8 * 4), 8, 4)
  Xd[, 4] <- Xd[, 1] + Xd[, 2]
  expect_warning(
    o2 <- EFAtools:::.screen_outliers(Xd, mcd_alpha = 0.5, outlier_cutoff = 0.975,
                                      seed = 1),
    class = "efa_screen_no_outliers"
  )
  expect_s3_class(o2, "efa_screen_no_outliers")
  # each of the three ways the complete-case covariance can fail names itself, so the
  # report never sends a reader after collinearity when the cause is missingness or a
  # non-finite value
  expect_match(o2$reason, "linearly dependent")

  set.seed(6)
  Xnf <- matrix(rnorm(8 * 4), 8, 4)
  Xnf[3, 2] <- Inf
  expect_warning(
    o4 <- EFAtools:::.screen_outliers(Xnf, mcd_alpha = 0.5, outlier_cutoff = 0.975,
                                      seed = 1),
    class = "efa_screen_no_outliers"
  )
  expect_match(o4$reason, "not finite")

  X0 <- matrix(rnorm(5 * 4), 5, 4)
  X0[, 1] <- NA                              # no complete case at all
  expect_warning(
    o5 <- EFAtools:::.screen_outliers(X0, mcd_alpha = 0.5, outlier_cutoff = 0.975,
                                      seed = 1),
    class = "efa_screen_no_outliers"
  )
  expect_match(o5$reason, "too few complete cases")

  # near-collinear variables are skipped on the same terms as exactly collinear ones:
  # the covariance is numerically singular (rcond far below the double-precision epsilon)
  # even though it is positive definite, so it must not reach the distance step
  set.seed(6)
  Xn <- matrix(rnorm(8 * 4), 8, 4)
  Xn[, 4] <- Xn[, 1] + Xn[, 2] + 1e-9 * rnorm(8)
  expect_warning(
    o3 <- EFAtools:::.screen_outliers(Xn, mcd_alpha = 0.5, outlier_cutoff = 0.975,
                                      seed = 1),
    class = "efa_screen_no_outliers"
  )
  expect_s3_class(o3, "efa_screen_no_outliers")

  # a correlation matrix reaching the helper aborts
  expect_error(
    EFAtools:::.screen_outliers(test_models$baseline$cormat, mcd_alpha = 0.5,
                                outlier_cutoff = 0.975),
    class = "efa_screen_outliers_cormat"
  )
})

test_that("tied responses on coarse items fall back to classical distances (GRiPS)", {
  # GRiPS items are 6-point and heavily tied: more than h respondents give the same answer
  # on at least one item pair, so an h-subset sits exactly on the hyperplane x_i = x_j and
  # no robust MCD covariance exists. The outlier diagnostic then falls back to classical
  # Mahalanobis distances with a warning. A covering subset of determinant zero is the
  # smallest one there can be, so which random subsets the search happens to draw cannot
  # change that conclusion: the verdict, and the whole diagnostic behind it, must come out
  # the same at every seed rather than depending on one.
  seeds <- c(1L, 2L, 4L, 7L, 12L)
  scrs <- lapply(seeds, function(s) {
    expect_warning(scr <- efa_screen(GRiPS_raw, seed = s), class = "efa_screen_mcd_fallback")
    scr
  })
  outs <- lapply(scrs, `[[`, "outliers")
  expect_true(all(vapply(outs, function(o) o$method, character(1)) == "classical"))
  expect_true(all(vapply(outs, function(o) identical(o$flagged, outs[[1]]$flagged),
                         logical(1))))

  scr <- scrs[[1]]
  expect_length(scr$outliers$distances, nrow(GRiPS_raw))
  expect_equal(scr$outliers$cutoff, sqrt(qchisq(0.975, ncol(GRiPS_raw))))

  # with n = 810 well above 2p = 16 the cause cannot be too few complete cases, and the
  # correlation matrix is well conditioned, so the reported reason must be the exact fit
  expect_gt(nrow(GRiPS_raw), 2 * ncol(GRiPS_raw))
  expect_lt(sqrt(scr$condition), 10)
  expect_match(scr$outliers$fallback_reason, "exact fit")
})

test_that("the outlier diagnostic does not depend on the variables' units", {
  # The Mahalanobis distance, and the MCD location and scatter behind it, are affine
  # equivariant: recording a variable in cents rather than in euros cannot make anyone an
  # outlier who was not one before. Multiplying one column by 100 - a narrower spread of
  # units than an income variable beside a Likert item - must therefore leave the verdict,
  # the flagged rows and the distances where they were, and must carry the reported centre
  # and scatter into the new units exactly.
  set.seed(1809)
  X <- matrix(stats::rnorm(300 * 4), 300, 4)
  X[1:6, ] <- X[1:6, ] + 5                      # a block of genuine outliers
  colnames(X) <- paste0("V", seq_len(4))
  d <- c(100, 1, 1, 1)
  Xs <- sweep(X, 2L, d, "*", check.margin = FALSE)

  o <- efa_screen(X, seed = 1)$outliers
  os <- efa_screen(Xs, seed = 1)$outliers

  expect_equal(o$method, "mcd")
  expect_equal(os$method, "mcd")
  expect_identical(os$flagged, o$flagged)
  expect_equal(os$distances, o$distances, tolerance = 1e-8)
  expect_equal(os$center, o$center * d, tolerance = 1e-8)
  expect_equal(os$cov, o$cov * outer(d, d), tolerance = 1e-8)

  # a scale far enough apart to make the scatter singular once it is back in the supplied
  # units: the distances have to be taken where the conditioning was checked, or the
  # inversion inside mahalanobis() refuses a scatter that was passed as usable
  dw <- c(1e9, 1, 1, 1)
  ow <- efa_screen(sweep(X, 2L, dw, "*", check.margin = FALSE), seed = 1)$outliers
  expect_equal(ow$method, "mcd")
  expect_identical(ow$flagged, o$flagged)

  # the rescaling is taken from the complete cases, so an incomplete row cannot leave a
  # column unscaled - the median absolute deviation of a column holding a missing value is
  # itself missing - and quietly put the diagnostic back at the mercy of the units
  Xm <- Xs
  Xm[c(50, 120, 200), 1] <- NA
  expect_equal(efa_screen(Xm, seed = 1)$outliers$method, "mcd")
})

test_that("a near-collinear fallback is reported as collinearity, not as an exact fit", {
  # The robust scatter can be unusable for reasons that send a reader to quite different
  # places: near-collinear variables are dealt with by dropping a redundant item, while an
  # exact fit is a property of a covering subset of the respondents (tied answers on coarse
  # items) with nothing wrong with the variables at all. Here the fourth variable is a
  # near-exact sum of the first two, so the complete-case covariance is itself
  # ill-conditioned and the recorded reason has to say so.
  set.seed(12)
  Xc <- matrix(stats::rnorm(300 * 4), 300, 4)
  Xc[, 4] <- Xc[, 1] + Xc[, 2] + 0.01 * stats::rnorm(300)
  expect_warning(
    o <- EFAtools:::.screen_outliers(Xc, mcd_alpha = 0.5, outlier_cutoff = 0.975, seed = 1),
    class = "efa_screen_mcd_fallback"
  )
  expect_equal(o$method, "classical")
  expect_match(o$fallback_reason, "ill-conditioned")
  expect_false(grepl("exact fit", o$fallback_reason, fixed = TRUE))
})

test_that("outlier_cutoff is bounded to the range in which it defines a cutoff", {
  # 0 puts the threshold at zero (everything flagged) and 1 at infinity (nothing
  # flagged); neither is an outlier diagnostic. The argument check runs before anything
  # is computed, so the cheap correlation-matrix input reaches it - and pinning both
  # endpoints against the values just outside them fixes where the bound sits, which the
  # error class alone (a checkmate assertion carries no package class) cannot do
  R <- test_models$baseline$cormat
  expect_no_error(efa_screen(R, N = 500, outlier_cutoff = 0.5))
  expect_no_error(efa_screen(R, N = 500, outlier_cutoff = 0.9999))
  expect_error(efa_screen(R, N = 500, outlier_cutoff = 0.4999), class = "simpleError")
  expect_error(efa_screen(R, N = 500, outlier_cutoff = 1), class = "simpleError")
  expect_error(efa_screen(R, N = 500, outlier_cutoff = 0), class = "simpleError")
})

test_that("fewer than three variables is refused rather than reported on", {
  # At p = 1 the KMO is 0/0 and Bartlett's test has no degrees of freedom; at p = 2 the
  # KMO is algebraically .5 for every correlation. Both produced a confident report.
  expect_error(efa_screen(data.frame(a = rnorm(50))),
               class = "efa_screen_too_few_vars")
  expect_error(efa_screen(data.frame(a = rnorm(50), b = rnorm(50))),
               class = "efa_screen_too_few_vars")
  expect_error(efa_screen(matrix(1, 1, 1), N = 100),
               class = "efa_screen_too_few_vars")
  expect_error(efa_screen(diag(2), N = 100), class = "efa_screen_too_few_vars")
  # Three variables is the first usable width and still runs.
  expect_s3_class(efa_screen(test_models$baseline$cormat[1:3, 1:3], N = 500),
                  "efa_screen")
})

test_that("an exactly singular matrix names the redundancy that caused it", {
  # The abort is deliberate -- every measure here is built from R^-1 -- but the user
  # screening for redundant items is exactly the one who cannot get a report, so the
  # condition carries the culprit pairs as data (`$pairs`), which is what is asserted
  # here rather than the wording that reports them.
  set.seed(4)
  dup <- data.frame(a = rnorm(120), b = rnorm(120), c = rnorm(120))
  dup$copy_of_a <- dup$a
  expect_error(efa_screen(dup), class = "efa_cor_singular")
  cnd <- tryCatch(efa_screen(dup), efa_cor_singular = function(e) e)
  expect_identical(cnd$pairs, "a-copy_of_a")
  # The abort is attributed to efa_screen(), not to the internal handler that builds it.
  expect_identical(rlang::call_name(conditionCall(cnd)), "efa_screen")

  # A correlation-matrix input takes the same route, and an unnamed one falls back to
  # the V-names used elsewhere in the report.
  R_dup <- stats::cor(dup)
  expect_error(efa_screen(R_dup, N = 120), class = "efa_cor_singular")
  R_un <- unname(R_dup)
  cnd_un <- tryCatch(efa_screen(R_un, N = 120), efa_cor_singular = function(e) e)
  expect_identical(cnd_un$pairs, "V1-V4")

  # Rank deficiency without a perfect pair (a sum score) has no pair to name, so it
  # reports none rather than an empty pair label.
  cnd_sum <- tryCatch(efa_screen(sing_raw), efa_cor_singular = function(e) e)
  expect_identical(cnd_sum$pairs, character())
})

test_that("a near-singular matrix still produces the full report", {
  # The distinction is the point: perturbing a duplicate by 1e-7 is statistically the
  # same data, and it must come back as a report with multicollinearity flags rather
  # than as the singular abort.
  set.seed(5)
  near <- data.frame(a = rnorm(120), b = rnorm(120), c = rnorm(120))
  near$copy_of_a <- near$a + rnorm(120, sd = 1e-7)
  scr <- suppressWarnings(efa_screen(near, seed = 1))
  expect_s3_class(scr, "efa_screen")
  expect_lt(scr$determinant, 1e-10)
  expect_gt(scr$condition, 30)
})

test_that("the multicollinearity verdict follows the condition index, not the determinant", {
  local_reproducible_output()

  # A clean one-factor population with every loading 0.5: every off-diagonal correlation
  # is exactly 0.25 and the smallest eigenvalue exactly 0.75 at every p, so this pool is
  # no worse conditioned at 60 variables than at 10 (index 4.58). The determinant is the
  # product of the p eigenvalues, so it still falls below the 0.00001 that is commonly
  # quoted - a cut-off on it flags a pool with nothing wrong with it.
  p <- 60
  L <- matrix(0.5, p, 1)
  R60 <- L %*% t(L)
  diag(R60) <- 1
  dimnames(R60) <- list(paste0("v", seq_len(p)), paste0("v", seq_len(p)))
  scr_wide <- efa_screen(R60, N = 500)
  expect_lt(scr_wide$determinant, 1e-5)
  expect_lt(sqrt(scr_wide$condition), 10)

  txt <- cli::ansi_strip(paste(format(scr_wide), collapse = " "))
  # the determinant is still reported as a number, but carries no verdict of its own
  expect_match(txt, "Determinant: 0.0000006", fixed = TRUE)
  expect_match(txt, "so the condition index below carries the verdict", fixed = TRUE)
  expect_no_match(txt, "Multicollinearity likely", fixed = TRUE)
  # and no recommendation sends the reader after redundant items that do not exist
  expect_no_match(txt, "nearly linearly dependent", fixed = TRUE)

  # A genuinely ill-conditioned matrix still gets the verdict and the recommendation. A
  # duplicate perturbed by 1e-7 puts the index near 2e7, so far above the highest step of
  # the progression (300) that the graded label cannot move with the BLAS.
  set.seed(5)
  near2 <- data.frame(a = rnorm(120), b = rnorm(120), c = rnorm(120))
  near2$copy_of_a <- near2$a + rnorm(120, sd = 1e-7)
  scr_ill <- suppressWarnings(efa_screen(near2, seed = 1))
  expect_gt(sqrt(scr_ill$condition), 300)
  txt_ill <- cli::ansi_strip(paste(format(scr_ill), collapse = " "))
  expect_match(txt_ill, "which indicates a near linear dependency", fixed = TRUE)
  expect_match(txt_ill, "its relative strength is very strong", fixed = TRUE)
  expect_match(txt_ill, "A high condition index indicates multicollinearity", fixed = TRUE)
})

test_that("the condition-index verdict is graded on Belsley's progression", {
  local_reproducible_output()

  # An equicorrelated matrix of three variables has eigenvalues 1 + 2r and 1 - r (twice),
  # so its condition index is exactly sqrt((1 + 2r) / (1 - r)). Each r below puts the
  # index far from the band edges (10, 30, 100, 300), so the expected verdict cannot move
  # with the BLAS: the four indexes are 17.3, 44.7, 141.4, and 447.2.
  eqcor <- function(r) {
    R <- matrix(r, 3, 3)
    diag(R) <- 1
    dimnames(R) <- list(paste0("v", 1:3), paste0("v", 1:3))
    R
  }
  screened <- function(r) {
    scr <- efa_screen(eqcor(r), N = 500)
    expect_equal(sqrt(scr$condition), sqrt((1 + 2 * r) / (1 - r)), tolerance = 1e-6,
                 label = paste0("condition index at r = ", r))
    cli::ansi_strip(paste(format(scr), collapse = " "))
  }

  # Above 10 but not above the 30 that flags a near dependency: no flag, no recommendation.
  mid <- screened(0.99)
  expect_match(mid, "The index is above 10, but not above 30", fixed = TRUE)
  expect_no_match(mid, "A high condition index indicates", fixed = TRUE)

  # Flagged, and graded by the position on the progression the index sits at.
  for (case in list(list(r = 0.99850,  strength = "moderate"),
                    list(r = 0.99985,  strength = "strong"),
                    list(r = 0.999985, strength = "very strong"))) {
    txt <- screened(case$r)
    expect_match(txt, paste("its relative strength is", case$strength), fixed = TRUE)
    expect_match(txt, "A high condition index indicates", fixed = TRUE)
  }
})

# singular / non-positive-definite inputs
test_that("errors and warnings are thrown correctly", {
  # The rank-deficient raw data and its correlation matrix (`sing_raw` / `sing_cor`, with
  # `sing_N` rows) and the invertible-but-indefinite `cor_nposdef` come from
  # helper-singular.R, so both branches are exercised on the same inputs everywhere.
  expect_error(efa_screen(1:5), class = "efa_input_not_matrix")
  expect_error(efa_screen(sing_raw), class = "efa_cor_singular")
  expect_error(efa_screen(sing_cor, N = sing_N), class = "efa_cor_singular")
  expect_warning(efa_screen(cor_nposdef, N = 10), class = "efa_cor_smoothed")
})

test_that("a flag rate far above the nominal one is reported as a calibration issue", {
  # iris is three species clusters, so it is not elliptically distributed: a high-breakdown
  # estimate fitted to the most concentrated half legitimately calls a large share of the
  # rest distant. Far more flags than the 2.5% cutoff admits is evidence of that, not a list
  # of contaminated cases to work through.
  recs <- EFAtools:::.screen_recommendations(scr_iris, raw = TRUE,
                                             kmo = scr_iris$kmo$KMO, bart_sig = TRUE,
                                             mvn_nonnormal = TRUE, multicollinear = FALSE,
                                             digits = 3)
  nf <- length(scr_iris$outliers$flagged)
  expect_gt(nf / scr_iris$outliers$n_complete, 5 * (1 - scr_iris$settings$outlier_cutoff))
  expect_true(any(grepl("not elliptically distributed", recs$message, fixed = TRUE)))
  expect_false(any(grepl("before down-weighting or excluding", recs$message, fixed = TRUE)))

  # a rate in line with the cutoff keeps the ordinary inspect-them advice
  recs_ok <- EFAtools:::.screen_recommendations(scr_out, raw = TRUE,
                                                kmo = scr_out$kmo$KMO, bart_sig = TRUE,
                                                mvn_nonnormal = FALSE, multicollinear = FALSE,
                                                digits = 3)
  expect_true(any(grepl("before down-weighting or excluding", recs_ok$message, fixed = TRUE)))

  # the rate alone does not license the conclusion: a short flagged list is a list to
  # inspect however far above the nominal rate it sits (here 3 of 20 = 15% against 2.5%)
  few <- scr_iris
  few$outliers$flagged <- scr_iris$outliers$flagged[1:3]
  few$outliers$n_complete <- 20L
  recs_few <- EFAtools:::.screen_recommendations(few, raw = TRUE, kmo = few$kmo$KMO,
                                                 bart_sig = TRUE, mvn_nonnormal = TRUE,
                                                 multicollinear = FALSE, digits = 3)
  expect_gt(3 / 20, 5 * (1 - few$settings$outlier_cutoff))
  expect_true(any(grepl("before down-weighting or excluding", recs_few$message, fixed = TRUE)))
  expect_false(any(grepl("not elliptically distributed", recs_few$message, fixed = TRUE)))
})

test_that("the per-variable display marks the missing percentage and hides an empty flags column", {
  # `missing` is a percentage; the display says so
  expect_true("missing%" %in% names(EFAtools:::.screen_per_item_display(scr_raw, 3)))

  # all-continuous data: every flags entry is NA ("no category screening applies"), which
  # prints as a column of <NA>, so the column is dropped instead
  expect_true(all(is.na(scr_iris$per_item$flags)))
  expect_false("flags" %in% names(EFAtools:::.screen_per_item_display(scr_iris, 3)))

  # with some categorical variables the column stays and the NAs render as a dash
  set.seed(12)
  dat <- data.frame(cont = rnorm(80), cat = sample(1:5, 80, replace = TRUE),
                    v3 = rnorm(80))
  scr <- efa_screen(dat, seed = 1)
  disp <- EFAtools:::.screen_per_item_display(scr, 3)
  expect_true("flags" %in% names(disp))
  expect_equal(disp["cont", "flags"], "-")
})

test_that("the complete-case denominator is named when rows are incomplete", {
  # The normality tests and the outlier diagnostic are computed from the complete cases,
  # which under non-ignorable missingness are not a random subsample of the rows
  # supplied. Both bases therefore have to appear, so neither the normality verdict nor
  # the outlier rate can be read against the wrong one.
  local_reproducible_output(width = 200)
  set.seed(31)
  n3 <- 240L
  S3 <- matrix(0.4, p_out, p_out)
  diag(S3) <- 1
  X3 <- matrix(rnorm(n3 * p_out), n3, p_out) %*% chol(S3)
  colnames(X3) <- paste0("v", seq_len(p_out))
  X3[1:60, 1] <- NA
  scr_na <- efa_screen(X3, seed = 1)
  expect_equal(scr_na$settings$n_obs, 240L)
  expect_equal(scr_na$normality$n_complete, 180L)

  txt <- cli::ansi_strip(paste(format(scr_na), collapse = " "))
  expect_match(txt, "Computed from 180 complete cases of the 240 rows supplied",
               fixed = TRUE)
  expect_match(txt, "of 180 complete cases of the 240 rows supplied", fixed = TRUE)

  # With no missing values the two bases coincide, so the report stays as it was.
  scr_full <- efa_screen(X3[stats::complete.cases(X3), ], seed = 1)
  txt_full <- cli::ansi_strip(paste(format(scr_full), collapse = " "))
  expect_no_match(txt_full, "rows supplied", fixed = TRUE)
  expect_match(txt_full, "of 180 observations", fixed = TRUE)
})

test_that("print output is stable", {
  local_reproducible_output()

  # correlation-matrix input: KMO + Bartlett verdicts, multicollinearity, and the
  # per-variable MSA/SMC table, closing with the raw-data note (no raw-only sections)
  expect_snapshot(print(scr_cor), transform = scrub_num)

  # raw input: the full report incl. the per-item table, multivariate normality, the
  # classical-fallback outlier line, and the recommendations block (the non-normality,
  # over-powered-Bartlett caveat, sparse-category, and outlier recommendations)
  expect_snapshot(print(scr_raw), transform = scrub_num)

  # correlation-matrix input without N: Bartlett's test is reported as not computed
  expect_snapshot(print(scr_non), transform = scrub_num)
})

test_that("the report tracks a narrow console", {
  # Every other print snapshot in the package is recorded at the default width 80, so
  # nothing pins the behaviour of the report's advisory lines once the console is
  # narrower than they are. efa_screen() is the widest report in the package (banded
  # verdicts with a parenthetical rule attached to each), so it is the one that has to
  # reflow rather than overflow.
  local_reproducible_output(width = 60)

  expect_snapshot(print(scr_cor), transform = scrub_num)
})

rm(scr_cor, scr_raw, scr_iris, scr_nona, scr_non, dat_nonames, scr_out, X_out, Sig_out,
   n_out, p_out, inj)
