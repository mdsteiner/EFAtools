# the uncomputable-correlation abort names the offending columns

    Code
      .prepare_cor_input(dat_fct, inform_from_data = FALSE)
    Condition
      Error:
      ! The correlation matrix could not be computed from the raw data.
      x Column "item_1" is not numeric.
      i Ordinal items stored as factors or character strings are correlated with `cor_method = "poly"` ("tetra" for binary items); anything else has to be converted or dropped.

---

    Code
      .prepare_cor_input(dat_const, inform_from_data = FALSE)
    Condition
      Error:
      ! The correlation matrix could not be computed from the raw data.
      x Column "item_2" has zero variance.
      i A constant variable correlates with nothing; drop it before the analysis.

---

    Code
      .prepare_cor_input(dat_both, inform_from_data = FALSE)
    Condition
      Error:
      ! The correlation matrix could not be computed from the raw data.
      x Column "item_1" is not numeric.
      i Ordinal items stored as factors or character strings are correlated with `cor_method = "poly"` ("tetra" for binary items); anything else has to be converted or dropped.
      x Column "item_2" has zero variance.
      i A constant variable correlates with nothing; drop it before the analysis.

---

    Code
      .prepare_cor_input(cbind(c(1, 2, 3, 4), rep(2, 4), c(4, 1, 3, 2)),
      inform_from_data = FALSE)
    Condition
      Error:
      ! The correlation matrix could not be computed from the raw data.
      x Column "V2" has zero variance.
      i A constant variable correlates with nothing; drop it before the analysis.

---

    Code
      .prepare_cor_input(dat_lgl, inform_from_data = FALSE)
    Condition
      Error:
      ! The correlation matrix could not be computed from the raw data.
      x Column "b" has zero variance.
      i A constant variable correlates with nothing; drop it before the analysis.

---

    Code
      .prepare_cor_input(dat_inf, inform_from_data = FALSE)
    Condition
      Error:
      ! The correlation matrix could not be computed from the raw data.
      x Column "item_2" contains infinite values.
      i An infinite value has no correlation with anything; check for a division by zero or an out-of-range missing-value code.

---

    Code
      .prepare_cor_input(matrix(c(1, 2, 3), nrow = 1), inform_from_data = FALSE)
    Condition
      Error:
      ! The correlation matrix could not be computed from the raw data.
      i Check that all columns are numeric and have non-zero variance. Correlations also need at least two observations.

---

    Code
      .prepare_cor_input(dat_lw, acov = "full", inform_from_data = FALSE)
    Message
      i An asymptotic covariance requires complete cases; incomplete rows were dropped (listwise), overriding `use = "pairwise.complete.obs"`.
    Condition
      Error:
      ! The correlation matrix could not be computed from the raw data.
      x Column "c" has zero variance.
      i The variance was computed on the listwise-complete rows an `acov` needs, which can be far fewer than the data supplied; either supply more complete cases or drop the constant column.

---

    Code
      .prepare_cor_input(dat_many, inform_from_data = FALSE)
    Condition
      Error:
      ! The correlation matrix could not be computed from the raw data.
      x Columns "a" and "b" are not numeric.
      i Ordinal items stored as factors or character strings are correlated with `cor_method = "poly"` ("tetra" for binary items); anything else has to be converted or dropped.
      x Columns "c1" and "c2" have zero variance.
      i A constant variable correlates with nothing; drop them before the analysis.

---

    Code
      .prepare_cor_input(dat_wide, inform_from_data = FALSE)
    Condition
      Error:
      ! The correlation matrix could not be computed from the raw data.
      x Columns "V1", "V2", "V3", "V4", "V5", and 2 more have zero variance.
      i A constant variable correlates with nothing; drop them before the analysis.

# the non-symmetric abort names the triangle that is actually empty

    Code
      .assert_cor_input(R_lower)
    Condition
      Error:
      ! `x` looks like a correlation matrix but is not symmetric.
      x Its diagonal is all ones and every entry lies in [-1, 1], but `x[i, j]` and `x[j, i]` differ.
      i Only the lower triangle carries entries; mirror it: `x[upper.tri(x)] <- t(x)[upper.tri(x)]`.

---

    Code
      .assert_cor_input(R_upper)
    Condition
      Error:
      ! `x` looks like a correlation matrix but is not symmetric.
      x Its diagonal is all ones and every entry lies in [-1, 1], but `x[i, j]` and `x[j, i]` differ.
      i Only the upper triangle carries entries; mirror it: `x[lower.tri(x)] <- t(x)[lower.tri(x)]`.

---

    Code
      .assert_cor_input(R_na)
    Condition
      Error:
      ! `x` looks like a correlation matrix but is not symmetric.
      x Every entry that is present lies in [-1, 1] and no diagonal entry departs from one, but part of the matrix is missing and `x[i, j]` and `x[j, i]` differ.
      i Only the lower triangle carries entries; mirror it: `x[upper.tri(x)] <- t(x)[upper.tri(x)]`.

---

    Code
      .assert_cor_input(R_mismatch)
    Condition
      Error:
      ! `x` looks like a correlation matrix but is not symmetric.
      x Its diagonal is all ones and every entry lies in [-1, 1], but `x[i, j]` and `x[j, i]` differ.
      i Both triangles carry entries but they disagree, so neither can be mirrored onto the other; check the entries of the pairs that differ.

---

    Code
      .assert_cor_input(R_none)
    Condition
      Error:
      ! `x` looks like a correlation matrix but is not symmetric.
      x Every entry that is present lies in [-1, 1] and no diagonal entry departs from one, but part of the matrix is missing and `x[i, j]` and `x[j, i]` differ.
      i Neither triangle carries a correlation, so there is nothing to mirror; enter the off-diagonal correlations.

---

    Code
      .assert_cor_input(R_lower, raw_only = TRUE)
    Condition
      Error:
      ! `x` looks like a (non-symmetric) correlation matrix, not a data frame/matrix of raw data.
      x Its diagonal is all ones and every entry lies in [-1, 1], but `x[i, j]` and `x[j, i]` differ.
      i Supply the raw observations the correlation matrix was computed from.

