# Unit tests for the shared correlation-input preparer: correlation-matrix
# detection, the N-handling policy, and the singular / non-positive-definite
# branches.

set.seed(7)                                            # keep the drawn columns deterministic
x <- rnorm(10)                                         # two varying columns for the
y <- rnorm(10)                                         # constant-column fixture below
cormat <- test_models$baseline$cormat                  # well-behaved correlation matrix

test_that("a correlation matrix is detected and returned unchanged", {
  prep <- .prepare_cor_input(cormat, N = 500)

  expect_named(prep, c("R", "N", "is_cormat", "weights", "Gamma", "fiml"))
  expect_true(prep$is_cormat)
  expect_equal(prep$R, cormat)
  expect_equal(prep$N, 500)
  expect_null(prep$fiml)                                # only the FIML path populates it

  # a clean correlation matrix triggers no message, warning, or error
  expect_silent(.prepare_cor_input(cormat, N = 500))
})

test_that("an exact -1 correlation is within the accepted range", {
  R <- matrix(c(1, -1, -1, 1), 2, 2)
  expect_true(.is_cormat(R))
  expect_silent(.assert_cor_input(R))
})

test_that("a correlation matrix supplied as a data frame is detected and coerced", {
  # read.csv(file, row.names = 1) is how a published correlation table usually reaches R, and
  # it yields a data frame. It has to be classified as a correlation matrix -- diag() on a data
  # frame is a hard base error, so the classifier coerces first -- and the preparer has to hand
  # a real matrix downstream, since determinant(), the eigendecompositions and the compiled
  # estimators all reject a data frame.
  cormat_df <- as.data.frame(cormat)

  expect_true(.is_cormat(cormat_df))
  expect_silent(.assert_cor_input(cormat_df))

  prep_df <- .prepare_cor_input(cormat_df, N = 500)
  expect_true(prep_df$is_cormat)
  expect_true(is.matrix(prep_df$R))
  expect_equal(prep_df$R, cormat)
  expect_equal(prep_df$N, 500)
  expect_silent(.prepare_cor_input(cormat_df, N = 500))       # not mistaken for raw data

  # every entry point that accepts a correlation matrix returns what the matrix route returns
  expect_equal(efa_kmo(cormat_df), efa_kmo(cormat))
  expect_equal(efa_bartlett(cormat_df, N = 500), efa_bartlett(cormat, N = 500))
  expect_equal(efa_map(cormat_df), efa_map(cormat))
  expect_equal(efa_kgc(cormat_df), efa_kgc(cormat))
  expect_equal(efa_scree(cormat_df), efa_scree(cormat))
  expect_equal(efa_ekc(cormat_df, N = 500), efa_ekc(cormat, N = 500))
  expect_equal(efa_fit(cormat_df, n_factors = 3, N = 500),
               efa_fit(cormat, n_factors = 3, N = 500))
  expect_equal(efa_screen(cormat_df, N = 500), efa_screen(cormat, N = 500))
  expect_equal(efa_retain(cormat_df, N = 500, criteria = c("KGC", "MAP")),
               efa_retain(cormat, N = 500, criteria = c("KGC", "MAP")))
  # including the superseded wrappers, which forward to the same preparer
  expect_equal(EFA(cormat_df, n_factors = 3, N = 500), EFA(cormat, n_factors = 3, N = 500))
  expect_equal(KMO(cormat_df), KMO(cormat))

  # a non-numeric square data frame is still not a correlation matrix, and a covariance data
  # frame is still rejected as one
  expect_false(.is_cormat(data.frame(a = factor(c("x", "y")), b = c(1, 2))))
  expect_error(.assert_cor_input(as.data.frame(diag(c(2, 3, 4)))),
               class = "efa_input_is_covmat")
})

test_that("a 1x1 input is not classified as a correlation matrix", {
  # .is_covmat() has always required at least two columns; .is_cormat() applies the same guard,
  # so a 1x1 input takes the raw-data route and is rejected there rather than being accepted
  # and failing somewhere downstream.
  expect_false(.is_cormat(matrix(1)))
  expect_false(.is_covmat(matrix(1)))
})

test_that("a non-symmetric correlation matrix is named rather than reported as singular", {
  # Neither classifier accepts an asymmetric matrix, so without this branch it reaches
  # stats::cor() as if its p columns were p cases and is reported as singular -- sending the
  # user after collinearity that does not exist.
  Rl <- cormat
  Rl[upper.tri(Rl)] <- 0                       # a lower triangle transcribed on its own
  expect_error(.assert_cor_input(Rl), class = "efa_input_not_symmetric")

  Ra <- cormat
  Ra[1, 2] <- Ra[1, 2] + 1e-7                  # values that do not mirror to machine precision
  expect_error(.assert_cor_input(Ra), class = "efa_input_not_symmetric")

  # an empty (NA) upper triangle is the same transcription error, and must take the same route
  # rather than the missing-values one
  Rna <- cormat
  Rna[upper.tri(Rna)] <- NA
  expect_error(.assert_cor_input(Rna), class = "efa_input_not_symmetric")

  # a data frame is classified identically
  expect_error(.assert_cor_input(as.data.frame(Rl)), class = "efa_input_not_symmetric")

  # the rejection reaches every entry point, including the raw-data-only route, which points
  # at the observations instead of at mirroring
  expect_error(suppressMessages(efa_kmo(Rl)), class = "efa_input_not_symmetric")
  expect_error(suppressMessages(efa_fit(Rl, n_factors = 3, N = 500)),
               class = "efa_input_not_symmetric")
  expect_error(suppressMessages(efa_screen(Rl)), class = "efa_input_not_symmetric")
  expect_error(suppressMessages(efa_map(Rl)), class = "efa_input_not_symmetric")
  expect_error(suppressMessages(efa_cd(Rl)), class = "efa_input_not_symmetric")

  # the suggested repair is accepted and recovers the original matrix
  Rl[upper.tri(Rl)] <- t(Rl)[upper.tri(Rl)]
  expect_silent(.assert_cor_input(Rl))
  expect_equal(Rl, cormat)

  # The hint has to name the triangle that is actually empty. Mirroring the lower triangle
  # onto the upper one is the repair for a lower-triangle transcription and DESTROYS an
  # upper-triangle one -- it overwrites every entered correlation with the empty triangle and
  # leaves the identity matrix, which passes every later check and analyses silently.
  Ru <- cormat
  Ru[lower.tri(Ru)] <- 0                       # only the UPPER triangle was entered
  expect_error(.assert_cor_input(Ru), class = "efa_input_not_symmetric")
  Ru[lower.tri(Ru)] <- t(Ru)[lower.tri(Ru)]    # the repair this input needs
  expect_silent(.assert_cor_input(Ru))
  expect_equal(Ru, cormat)
  # the other direction really would have wiped it
  Rwipe <- cormat
  Rwipe[lower.tri(Rwipe)] <- 0
  Rwipe[upper.tri(Rwipe)] <- t(Rwipe)[upper.tri(Rwipe)]
  expect_equal(unname(Rwipe), diag(ncol(cormat)))

  # square raw data is untouched: it is asymmetric but not unit-diagonal
  set.seed(3)
  expect_silent(.assert_cor_input(matrix(stats::runif(9, -1, 1), 3)))
  expect_silent(.assert_cor_input(GRiPS_raw))
})

test_that("raw data is converted to a correlation matrix", {
  expect_message(prep <- .prepare_cor_input(GRiPS_raw), class = "efa_cor_from_data")

  expect_false(prep$is_cormat)
  expect_equal(prep$R,
               stats::cor(GRiPS_raw, use = "pairwise.complete.obs",
                          method = "pearson"))
  expect_equal(colnames(prep$R), colnames(GRiPS_raw))
  expect_equal(prep$N, nrow(GRiPS_raw))

  # inform_from_data = FALSE silences the "computed from raw data" message
  expect_no_message(.prepare_cor_input(GRiPS_raw, inform_from_data = FALSE))
})

test_that("N_policy governs how N is handled", {
  # required: a correlation matrix without N aborts, but a supplied N is kept
  expect_error(.prepare_cor_input(cormat, N_policy = "required"),
               class = "efa_n_required")
  prep_req <- .prepare_cor_input(cormat, N = 200, N_policy = "required")
  expect_equal(prep_req$N, 200)

  # none: N is passed through untouched and the both-N-and-raw warning is muted
  expect_no_warning(
    prep_none <- suppressMessages(
      .prepare_cor_input(GRiPS_raw, N = 999, N_policy = "none")))
  expect_equal(prep_none$N, 999)

  # optional (default): supplying N with raw data warns and N is taken from data
  expect_warning(
    prep_opt <- suppressMessages(.prepare_cor_input(GRiPS_raw, N = 999)),
    class = "efa_n_from_data")
  expect_equal(prep_opt$N, nrow(GRiPS_raw))
})

test_that("listwise-deletion modes set N to the number of complete cases", {
  # raw data with two incomplete rows: under "complete.obs"/"na.or.complete"
  # stats::cor() drops them listwise, so N must be the complete-case count.
  set.seed(11)
  dat_miss <- matrix(rnorm(30 * 3), ncol = 3)
  dat_miss[2, 1] <- NA
  dat_miss[7, 3] <- NA
  n_complete <- sum(stats::complete.cases(dat_miss))

  for (u in c("complete.obs", "na.or.complete")) {
    prep <- suppressMessages(.prepare_cor_input(dat_miss, use = u))
    expect_equal(prep$N, n_complete)
    expect_lt(prep$N, nrow(dat_miss))
  }

  # other modes keep the raw row count
  prep_pw <- suppressMessages(
    .prepare_cor_input(dat_miss, use = "pairwise.complete.obs"))
  expect_equal(prep_pw$N, nrow(dat_miss))
})

test_that("a non-positive-definite matrix is smoothed, or aborts under posdef_abort", {
  # default: the matrix is smoothed to positive-definiteness with a classed warning
  expect_warning(prep_sm <- .prepare_cor_input(cor_nposdef, N = 10),
                 class = "efa_cor_smoothed")
  expect_true(prep_sm$is_cormat)
  # Positive definite with a margin, not merely positive: the smallest eigenvalue of a smoothed
  # matrix sits on the projection's own floor (about 1e-10), so asserting a floor of 1e-12
  # states the guarantee the warning makes. Asserting only "> 0" would instead be satisfied by
  # an unsmoothed matrix whose smallest eigenvalue happens to be a positive rounding artefact.
  expect_gt(min(eigen(prep_sm$R, symmetric = TRUE, only.values = TRUE)$values), 1e-12)

  # posdef_abort = TRUE turns the same input into an error instead
  expect_error(.prepare_cor_input(cor_nposdef, N = 10, posdef_abort = TRUE),
               class = "efa_cor_not_posdef")

  # an already-smoothed matrix is not re-flagged on a subsequent pass (the
  # non-positive-definite check uses psych::cor.smooth()'s own trigger), and it
  # is returned unchanged
  smoothed <- prep_sm$R
  expect_no_warning(prep_again <- .prepare_cor_input(smoothed, N = 10))
  expect_equal(prep_again$R, smoothed)
})

test_that("the positive-definite projection floors the spectrum and keeps the correlation form", {
  # The fallback for a matrix psych::cor.smooth() declines to act on. It has to return a
  # correlation matrix -- unit diagonal, symmetric, dimnames intact -- whose smallest
  # eigenvalue is clear of the boundary rather than sitting at the rounding floor.
  R_pd <- .project_cor_pd(cor_nposdef)
  expect_gt(min(eigen(R_pd, symmetric = TRUE, only.values = TRUE)$values), 1e-12)
  expect_equal(diag(R_pd), rep(1, ncol(cor_nposdef)))
  expect_true(isSymmetric(R_pd))

  # the eigendecomposition drops dimnames, so they are restored explicitly
  named <- cor_nposdef
  dimnames(named) <- list(letters[1:3], letters[1:3])
  expect_equal(dimnames(.project_cor_pd(named)), dimnames(named))

  # an already positive definite matrix passes through: no eigenvalue is floored, so the
  # spectrum is unchanged and the rebuild returns the input up to rounding
  expect_equal(.project_cor_pd(cormat), cormat, tolerance = 1e-8)
})

test_that("a matrix psych::cor.smooth() leaves untouched is projected anyway", {
  # cor.smooth() re-decides whether to act on its own eigendecomposition, which asks for the
  # eigenvectors where the gate that calls it reads a values-only one; LAPACK runs a different
  # driver for each, so on a matrix that is singular to working precision the two verdicts
  # disagree and cor.smooth() returns the matrix untouched. Which side of that a build lands
  # on cannot be arranged from the data, so the branch is driven by standing in for
  # cor.smooth() with the no-op it becomes when it declines. The matrix must still come back
  # positive definite, under the same classed warning.
  local_mocked_bindings(cor.smooth = function(x, ...) x, .package = "psych")

  expect_warning(prep <- .prepare_cor_input(cor_nposdef, N = 10),
                 class = "efa_cor_smoothed")
  expect_gt(min(eigen(prep$R, symmetric = TRUE, only.values = TRUE)$values), 1e-12)
  expect_equal(diag(prep$R), rep(1, ncol(cor_nposdef)))
})

test_that("NA in the raw-data correlation matrix aborts with a classed error", {
  # under use = "everything" a missing value leaves NAs in the computed
  # correlations; this must be a clear classed error, not an opaque base crash
  dat_na <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8,
                     8, 6, 7, 5, 3, 4, 1, 2,
                     2, 5, 1, 8, 3, 7, 4, NA), ncol = 3)
  expect_error(suppressMessages(.prepare_cor_input(dat_na, use = "everything")),
               class = "efa_cor_na")
  # under use = "all.obs" stats::cor() throws a hard base error before any NA
  # check; this must route to the same classed abort, not an unclassed crash
  expect_error(suppressMessages(.prepare_cor_input(dat_na, use = "all.obs")),
               class = "efa_cor_na")
})

test_that("non-numeric raw data aborts with a distinct classed error", {
  # A character column makes stats::cor() fail for a reason other than NAs (the
  # data here has none); report that as efa_cor_uncomputable rather than blaming
  # missing values.
  dat_chr <- data.frame(a = c("x", "y", "z", "w"), b = c(1, 2, 3, 4),
                        d = c(4, 3, 2, 1))
  expect_error(suppressMessages(.prepare_cor_input(dat_chr)),
               class = "efa_cor_uncomputable")

  # A SQUARE data frame reaches .is_cormat(), where a factor column would otherwise
  # dispatch to Ops.factor and fail with an unclassed base error; it must take the same
  # classed route as any other non-numeric raw data.
  dat_fct <- data.frame(a = factor(c("x", "y", "z")), b = c(1, 2, 3), d = c(3, 1, 2))
  expect_false(.is_cormat(dat_fct))
  expect_error(suppressMessages(.prepare_cor_input(dat_fct)),
               class = "efa_cor_uncomputable")

  # A constant column takes the same classed route. stats::cor() warns about it in base R's
  # own words, which are neither classed nor translated; that warning must not escape ahead
  # of the abort that diagnoses the same column.
  dat_const <- cbind(x, y, const = 1)
  expect_no_warning(
    expect_error(suppressMessages(.prepare_cor_input(dat_const)),
                 class = "efa_cor_uncomputable"))
  # every correlation method routes through the same abort
  for (m in c("pearson", "spearman", "kendall")) {
    expect_no_warning(
      expect_error(suppressMessages(.prepare_cor_input(dat_const, cor_method = m)),
                   class = "efa_cor_uncomputable"))
  }

  # A column the constancy test cannot judge must not turn the classed abort into an unclassed
  # base error: sd() is NaN for a column carrying an infinity and NA for fewer than two
  # values, and a single NA in that test would make the branch selection itself fail.
  expect_error(
    suppressMessages(.prepare_cor_input(cbind(a = c(1, 2, 3, 4), b = c(Inf, 1, 2, 3),
                                              d = c(4, 3, 2, 1)))),
    class = "efa_cor_uncomputable")
  expect_error(suppressMessages(.prepare_cor_input(matrix(c(1, 2, 3), nrow = 1))),
               class = "efa_cor_uncomputable")
  expect_error(suppressMessages(.prepare_cor_input(matrix(numeric(0), ncol = 3))),
               class = "efa_cor_uncomputable")
  # the 1x1 input the ncol < 2 guard sends down this route is rejected here, as claimed
  expect_error(suppressMessages(.prepare_cor_input(matrix(1))),
               class = "efa_cor_uncomputable")

  # Logical columns are usable: stats::cor() accepts is.numeric() OR is.logical(), so
  # dichotomous TRUE/FALSE items must not be blamed as non-numeric while the column that
  # actually broke the correlation goes unnamed.
  # (the wording is pinned by the snapshots below, which show only the constant column named)
  dat_lgl <- data.frame(a = c(TRUE, FALSE, TRUE, FALSE), b = rep(1, 4), d = c(4, 1, 3, 2))
  expect_error(suppressMessages(.prepare_cor_input(dat_lgl)),
               class = "efa_cor_uncomputable")
  # and an all-logical data set correlates without reaching this branch at all
  dat_lgl2 <- data.frame(a = c(TRUE, FALSE, TRUE, FALSE), b = c(TRUE, TRUE, FALSE, FALSE))
  expect_no_error(suppressMessages(.prepare_cor_input(dat_lgl2)))
})

test_that("the uncomputable-correlation abort names the offending columns", {
  local_reproducible_output()

  # Which column is at fault, and what to do about it, is the whole content of this message:
  # a factor or character column is what an ordinal item read from a file looks like, and
  # cor_method = "poly" is the fix rather than a reason to drop it.
  dat_fct <- data.frame(item_1 = factor(c("x", "y", "z")), item_2 = c(1, 2, 3),
                        item_3 = c(3, 1, 2))
  expect_snapshot(error = TRUE,
                  .prepare_cor_input(dat_fct, inform_from_data = FALSE))

  # a constant column has a different remedy, and both causes can be present at once
  dat_const <- data.frame(item_1 = c(1, 2, 3, 4), item_2 = rep(2, 4), item_3 = c(4, 1, 3, 2))
  expect_snapshot(error = TRUE,
                  .prepare_cor_input(dat_const, inform_from_data = FALSE))

  dat_both <- data.frame(item_1 = letters[1:4], item_2 = rep(2, 4), item_3 = c(4, 1, 3, 2))
  expect_snapshot(error = TRUE,
                  .prepare_cor_input(dat_both, inform_from_data = FALSE))

  # an unnamed matrix falls back to positional labels
  expect_snapshot(error = TRUE,
                  .prepare_cor_input(cbind(c(1, 2, 3, 4), rep(2, 4), c(4, 1, 3, 2)),
                                     inform_from_data = FALSE))

  # a logical column is usable, so only the constant column is named and the polychoric hint
  # is not offered for a column that needs nothing done to it
  dat_lgl <- data.frame(a = c(TRUE, FALSE, TRUE, FALSE), b = rep(1, 4), d = c(4, 1, 3, 2))
  expect_snapshot(error = TRUE,
                  .prepare_cor_input(dat_lgl, inform_from_data = FALSE))

  # an infinite value is its own cause, with its own remedy
  dat_inf <- data.frame(item_1 = c(1, 2, 3, 4), item_2 = c(Inf, 1, 2, 3),
                        item_3 = c(4, 3, 2, 1))
  expect_snapshot(error = TRUE,
                  .prepare_cor_input(dat_inf, inform_from_data = FALSE))

  # fewer than two observations leaves no column judgeable, so the generic bullet stands
  expect_snapshot(error = TRUE,
                  .prepare_cor_input(matrix(c(1, 2, 3), nrow = 1), inform_from_data = FALSE))

  # An asymptotic covariance reduces the data to its listwise-complete rows before the
  # correlation is attempted, so a column can be constant there and vary in the data the user
  # supplied. The verdict has to say which rows it describes rather than advising a drop.
  set.seed(9)
  dat_lw <- cbind(a = stats::rnorm(60), b = stats::rnorm(60), c = stats::rnorm(60))
  keep <- c(3, 11, 25, 40)
  dat_lw[-keep, "a"] <- NA                    # only these rows are complete
  dat_lw[keep, "c"] <- 7                      # constant among the complete rows only
  expect_gt(stats::sd(dat_lw[, "c"]), 0)      # but not in the data as supplied
  expect_snapshot(error = TRUE,
                  .prepare_cor_input(dat_lw, acov = "full", inform_from_data = FALSE))

  # the plural forms of both bullets, and the five-name cap on the list
  dat_many <- data.frame(a = letters[1:4], b = LETTERS[1:4],
                         c1 = rep(2, 4), c2 = rep(3, 4), ok = 1:4)
  expect_snapshot(error = TRUE,
                  .prepare_cor_input(dat_many, inform_from_data = FALSE))

  dat_wide <- as.data.frame(matrix(rep(1, 35), nrow = 5))
  dat_wide$ok <- 1:5
  expect_snapshot(error = TRUE,
                  .prepare_cor_input(dat_wide, inform_from_data = FALSE))
})

test_that("the non-symmetric abort names the triangle that is actually empty", {
  local_reproducible_output()

  # Which triangle to mirror is read off the data: printing the lower-triangle repair for an
  # upper-triangle transcription would tell the user to overwrite everything they entered.
  R_lower <- cormat[1:5, 1:5]
  R_lower[upper.tri(R_lower)] <- 0
  expect_snapshot(error = TRUE, .assert_cor_input(R_lower))

  R_upper <- cormat[1:5, 1:5]
  R_upper[lower.tri(R_upper)] <- 0
  expect_snapshot(error = TRUE, .assert_cor_input(R_upper))

  # an empty triangle entered as NA rather than 0 is the same transcription, and the message
  # must not claim the missing cells were checked and found to be in range
  R_na <- cormat[1:5, 1:5]
  R_na[upper.tri(R_na)] <- NA
  expect_snapshot(error = TRUE, .assert_cor_input(R_na))

  # both triangles filled but disagreeing: neither is authoritative, so no mirroring is offered
  R_mismatch <- cormat[1:5, 1:5]
  R_mismatch[1, 2] <- R_mismatch[1, 2] + 1e-7
  expect_snapshot(error = TRUE, .assert_cor_input(R_mismatch))

  # neither triangle filled in (a table pre-allocated but never transcribed): mirroring an
  # empty triangle onto the other would silently produce the identity matrix, so it must not
  # be suggested
  R_none <- cormat[1:5, 1:5]
  R_none[upper.tri(R_none)] <- NA
  R_none[lower.tri(R_none)] <- 0
  expect_snapshot(error = TRUE, .assert_cor_input(R_none))

  # the raw-data-only route points at the observations instead
  expect_snapshot(error = TRUE, .assert_cor_input(R_lower, raw_only = TRUE))
})

test_that("a singular matrix aborts unless the check is disabled", {
  # LAPACK returns the eigenvalues of a rank-deficient matrix with absolute error of order
  # eps * ||R||, so the smallest computed eigenvalue IS that rounding error and a rank
  # tolerance of a bare eps decides by coin flip. The shared fixture is sized so the ratio
  # sits well below ncol(R) * eps; it is the archetypal EFA singularity -- a composite
  # entered alongside its components -- and must be refused rather than inverted downstream.
  expect_equal(qr(sing_cor, tol = 1e-10)$rank, 9L)         # rank deficient by one
  expect_error(.prepare_cor_input(sing_cor, N = sing_N), class = "efa_cor_singular")
  expect_error(suppressMessages(.prepare_cor_input(sing_raw)),
               class = "efa_cor_singular")
  expect_error(suppressMessages(efa_fit(sing_cor, n_factors = 2, N = sing_N)),
               class = "efa_cor_singular")

  # The two causes take different bullets. Here the rank deficiency is in the variables, so
  # the message hands the reader on to the function that names them.
  expect_snapshot(error = TRUE, .prepare_cor_input(sing_cor, N = sing_N))

  # A sample no larger than the number of variables is the other cause, and no variable is at
  # fault: the matrix is taken about the sample means, so its rank is at most N - 1 and the
  # numbers are the whole diagnosis. The boundary is N <= p, so the equality case takes this
  # branch too -- with N = p, efa_screen() would report an exact linear combination that is an
  # artefact of the sample size rather than a property of the variables.
  expect_snapshot(error = TRUE,
                  suppressMessages(.prepare_cor_input(GRiPS_raw[1:4, ])))
  expect_snapshot(error = TRUE,
                  suppressMessages(.prepare_cor_input(GRiPS_raw[1:8, ])))

  # check_singular = FALSE skips the abort. Whether the smoothing branch is entered at all is a
  # property of the build rather than of the data: for an exactly rank-deficient matrix the
  # smallest computed eigenvalue IS the eigensolver's rounding, so its sign -- and with it the
  # verdict of a gate at .Machine$double.eps -- differs between LAPACK implementations. What
  # holds on every build is that nothing is refused and the matrix that comes back is a
  # correlation matrix that is not materially indefinite; the positive definiteness that
  # smoothing guarantees is pinned above on `cor_nposdef`, whose smallest eigenvalue is -0.41
  # and which therefore takes the branch everywhere.
  # The two outcomes are asserted as the disjunction they are, rather than as a bound loose
  # enough to admit both: either the branch was not entered and the matrix comes back exactly
  # as it went in, or it was and the result is positive definite with the projection's margin.
  prep_ns <- suppressWarnings(
    .prepare_cor_input(sing_cor, N = sing_N, check_singular = FALSE))
  expect_true(
    identical(prep_ns$R, sing_cor) ||
      min(eigen(prep_ns$R, symmetric = TRUE, only.values = TRUE)$values) > 1e-12)

  # and a well-conditioned matrix is nowhere near the threshold
  ev_ok <- eigen(cormat, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(abs(ev_ok)) / max(abs(ev_ok)),
            1e6 * ncol(cormat) * .Machine$double.eps)
  expect_silent(.prepare_cor_input(cormat, N = 500))
})

test_that("a covariance matrix is rejected with a pointer to cov2cor()", {
  # A covariance matrix is square and symmetric but carries variances on its diagonal, so
  # it is neither a correlation matrix nor raw data. Passed on, it would be fed to
  # stats::cor() as if its columns were cases -- singular by construction and a wrong
  # diagnosis -- so the shared input assertion rejects it up front.
  sds <- seq(1, 5, length.out = ncol(cormat))
  S <- diag(sds) %*% cormat %*% diag(sds)               # a valid covariance matrix

  expect_error(.assert_cor_input(S), class = "efa_input_is_covmat")
  # a diagonal (uncorrelated) covariance is caught too: variances on the diagonal
  expect_error(.assert_cor_input(diag(c(2, 3, 4))), class = "efa_input_is_covmat")

  # an INDEFINITE covariance matrix -- as arises from pairwise deletion or a table
  # transcribed at few decimals -- is a covariance matrix too: PSD is deliberately not
  # required (an indefinite correlation matrix is likewise accepted and smoothed), so
  # exactly the problematic inputs do not fall through to the raw-data misdiagnosis
  S_ind <- S
  S_ind[1, 2] <- S_ind[2, 1] <- sqrt(S[1, 1] * S[2, 2]) * 1.05
  expect_lt(min(eigen(S_ind, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_error(.assert_cor_input(S_ind), class = "efa_input_is_covmat")

  # a non-finite entry must not abort inside the classifier (eigen() is no longer called);
  # the matrix is simply not classified as a covariance matrix
  S_inf <- S
  S_inf[1, 2] <- S_inf[2, 1] <- Inf
  expect_false(.is_covmat(S_inf))
  expect_silent(.assert_cor_input(S_inf))

  # a zero diagonal (a distance-like matrix) is not a covariance matrix
  expect_false(.is_covmat(as.matrix(stats::dist(t(cormat[1:4, 1:4])))))

  # the rejection reaches every entry point that accepts a correlation matrix or raw data
  expect_error(suppressMessages(efa_fit(S, n_factors = 2, N = 200)),
               class = "efa_input_is_covmat")
  expect_error(suppressMessages(efa_kmo(S)), class = "efa_input_is_covmat")
  expect_error(suppressMessages(efa_screen(S)), class = "efa_input_is_covmat")
  expect_error(suppressMessages(efa_bartlett(S, N = 200)), class = "efa_input_is_covmat")
  # including the retention criteria, which reach the shared preparer by their own routes
  expect_error(suppressMessages(efa_parallel(S, N = 200)), class = "efa_input_is_covmat")
  expect_error(suppressMessages(efa_map(S)), class = "efa_input_is_covmat")

  # a raw-data-only function points at the observations, not at cov2cor() (which would only
  # produce a correlation matrix it rejects in turn)
  expect_error(suppressMessages(efa_cd(S)), class = "efa_input_is_covmat")

  # standardising first is accepted; a genuine correlation matrix and raw data still pass
  expect_silent(.assert_cor_input(stats::cov2cor(S)))
  expect_silent(.assert_cor_input(cormat))
  expect_silent(.assert_cor_input(GRiPS_raw))
})
