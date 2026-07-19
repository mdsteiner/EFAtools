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

test_that("settings are returned correctly", {
  expect_named(scr_cor$settings, c("N", "n_obs", "use", "cor_method"))
  expect_equal(scr_cor$settings$N, 500)
  expect_equal(scr_cor$settings$use, "pairwise.complete.obs")
  expect_equal(scr_cor$settings$cor_method, "pearson")

  # n_obs is the number of supplied raw rows (NA for a correlation-matrix input)
  expect_true(is.na(scr_cor$settings$n_obs))
  expect_equal(scr_raw$settings$n_obs, nrow(GRiPS_raw))
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
  # ordered factors; f1 never uses its third level ("c"), leaving an interior gap
  f1 <- factor(sample(c("a", "b", "d"), n, replace = TRUE),
               levels = c("a", "b", "c", "d"))
  f2 <- factor(sample(c("lo", "mid", "hi"), n, replace = TRUE),
               levels = c("lo", "mid", "hi"))
  f3 <- factor(sample(c("lo", "mid", "hi"), n, replace = TRUE),
               levels = c("lo", "mid", "hi"))
  dat <- data.frame(f1 = f1, f2 = f2, f3 = f3)

  scr <- suppressWarnings(efa_screen(dat, cor_method = "poly"))

  # data.matrix() coding (not as.matrix()): variance is finite and matches the codes
  expect_true(all(is.finite(scr$per_item$variance)))
  expect_equal(scr$per_item$variance, apply(data.matrix(dat), 2, var),
               tolerance = 1e-10, ignore_attr = TRUE)
  # f1's unused middle level (c -> code 3, between b = 2 and d = 4) is an empty category
  expect_match(scr$per_item["f1", "flags"], "empty")
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

test_that("a non-numeric (character) matrix is coded, not left as strings", {
  set.seed(11)
  cm <- matrix(sample(c("lo", "mid", "hi"), 300L, replace = TRUE), nrow = 100L,
               dimnames = list(NULL, c("v1", "v2", "v3")))
  scr <- suppressWarnings(efa_screen(cm, cor_method = "poly"))
  # data.matrix() alone leaves a character matrix as strings; efa_screen codes it
  # via as.data.frame() so the per-item diagnostics are computed on integer codes
  expect_true(all(is.finite(scr$per_item$variance)))
  expect_false(anyNA(scr$per_item$flags))
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
  expect_named(nm$mardia, c("skewness", "skewness_df", "skewness_p", "small_sample",
                            "kurtosis", "kurtosis_p", "b1p", "b2p"))
  expect_named(nm$hz, c("statistic", "p_value"))

  # GRiPS is complete, so every row feeds the tests and no small-sample correction applies
  expect_equal(nm$n_complete, nrow(GRiPS_raw))
  expect_false(nm$mardia$small_sample)

  # statistics finite and p-values in [0, 1]
  expect_true(all(is.finite(c(nm$mardia$skewness, nm$mardia$kurtosis, nm$hz$statistic))))
  expect_true(nm$mardia$skewness_p >= 0 && nm$mardia$skewness_p <= 1)
  expect_true(nm$mardia$kurtosis_p >= 0 && nm$mardia$kurtosis_p <= 1)
  expect_true(nm$hz$p_value >= 0 && nm$hz$p_value <= 1)

  # NULL for a correlation-matrix input
  expect_null(scr_cor$normality)
})

test_that("Mardia's skewness/kurtosis match psych after the divisor rescaling", {
  skip_on_cran()
  skip_if_not_installed("psych")

  X <- as.matrix(iris[, 1:4])
  n <- nrow(X)
  p <- ncol(X)
  md <- scr_iris$normality$mardia
  m <- psych::mardia(X, plot = FALSE)

  # efa_screen uses the biased (divisor-n) covariance; psych::mardia uses cov()
  # (divisor n - 1). The two are related exactly by (n / (n - 1))^k, so the
  # coefficients agree after rescaling.
  expect_equal(md$b1p, m$b1p * (n / (n - 1))^3, tolerance = 1e-6)
  expect_equal(md$b2p, m$b2p * (n / (n - 1))^2, tolerance = 1e-6)
  expect_equal(md$skewness, m$skew * (n / (n - 1))^3, tolerance = 1e-6)
  expect_false(md$small_sample)

  # kurtosis is the standardised (b2p - p(p + 2))
  expect_equal(md$kurtosis, (md$b2p - p * (p + 2)) * sqrt(n / (8 * p * (p + 2))),
               tolerance = 1e-10)
})

test_that("the Henze-Zirkler statistic and p-value match the reference value on iris", {
  # Henze-Zirkler (1990) statistic for the four numeric iris columns, computed once from
  # this implementation and cross-checked against MVN::hz(iris[, 1:4]) = 2.336.
  expect_equal(scr_iris$normality$hz$statistic, 2.3363942003, tolerance = 1e-6)
  # p-value pinned relatively (an absolute tolerance is meaningless at ~1e-19)
  expect_equal(scr_iris$normality$hz$p_value / 4.1413116299e-19, 1, tolerance = 1e-6)
})

test_that("the Mardia skewness statistic is small-sample corrected for n < 20", {
  skip_on_cran()
  skip_if_not_installed("psych")

  set.seed(1)
  Xs <- matrix(rnorm(15 * 3), 15, 3)
  nrm <- EFAtools:::.screen_normality(Xs)
  ms <- psych::mardia(Xs, plot = FALSE)

  expect_true(nrm$mardia$small_sample)
  # with n < 20 the reported statistic is Mardia's (1974) small-sample correction,
  # matching psych's small.skew after the divisor rescaling
  expect_equal(nrm$mardia$skewness, ms$small.skew * (15 / 14)^3, tolerance = 1e-6)
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
  skip_if_not_installed("MASS")

  set.seed(3)
  Sig <- matrix(0.35, 4, 4)
  diag(Sig) <- 1
  dat <- MASS::mvrnorm(200, rep(0, 4), Sig)
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
                    "n_complete"))
  expect_equal(o$method, "mcd")
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

test_that("collinear raw data falls back to classical distances (GRiPS)", {
  # GRiPS has an exact collinearity: 'enjoy' and 'commonly' are equal for most respondents,
  # so at least h observations lie on a hyperplane and the robust MCD covariance is singular.
  # The outlier diagnostic then falls back to classical Mahalanobis distances with a warning.
  expect_warning(scr <- efa_screen(GRiPS_raw, seed = 1), class = "efa_screen_mcd_fallback")
  expect_equal(scr$outliers$method, "classical")
  expect_length(scr$outliers$distances, nrow(GRiPS_raw))
  expect_equal(scr$outliers$cutoff, sqrt(qchisq(0.975, ncol(GRiPS_raw))))
})

# singular / non-positive-definite inputs
set.seed(42)
x <- rnorm(10)
y <- rnorm(10)
z <- x + y
dat_sing <- matrix(c(x, y, z), ncol = 3)
cor_sing <- stats::cor(dat_sing)
cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)

test_that("errors and warnings are thrown correctly", {
  expect_error(efa_screen(1:5), class = "efa_input_not_matrix")
  expect_error(efa_screen(dat_sing), class = "efa_cor_singular")
  expect_error(efa_screen(cor_sing, N = 10), class = "efa_cor_singular")
  expect_warning(efa_screen(cor_nposdef, N = 10), class = "efa_cor_smoothed")
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

rm(scr_cor, scr_raw, scr_iris, scr_nona, scr_non, dat_nonames, x, y, z, dat_sing,
   cor_sing, cor_nposdef, scr_out, X_out, Sig_out, n_out, p_out, inj)
