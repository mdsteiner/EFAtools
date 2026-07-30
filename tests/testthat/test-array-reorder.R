# Unit tests for .array_reorder() (R/averaging.R), which aligns each efa_average()
# replicate's factors to a common column order and sign before they are averaged.
# Inputs are literal arrays with a known permutation, so both the reordered output
# and the optimality of the chosen permutation can be asserted exactly.

arr_re_NA <- .array_reorder(L = array(c(rep(.6, 6), rep(0, 12),
                                     rep(0, 6), rep(.6, 6), rep(0, 6),
                                     rep(0, 6), rep(0, 6), rep(.6, 6),
                                     rep(.6, 6), rep(0, 12),
                                     rep(0, 6), rep(.6, 6), rep(0, 6),
                                     rep(0, 6), rep(0, 6), rep(.6, 6),
                                     rep(0, 6), rep(.6, 6), rep(0, 6),
                                     rep(-.6, 6), rep(0, 12),
                                     rep(0, 6), rep(0, 6), rep(.6, 6)),
                                   c(18, 3, 3)),
                            vars_accounted = array(c(rep(.2, 3),
                                                     rep(.3, 3),
                                                     rep(.4, 3),
                                                     rep(.2, 3),
                                                     rep(.3, 3),
                                                     rep(.4, 3),
                                                     rep(.3, 3),
                                                     rep(.2, 3),
                                                     rep(.4, 3)),
                                                   c(3, 3, 3)),
                         L_corres = array(as.numeric(abs(c(rep(.6, 6), rep(0, 12),
                                            rep(0, 6), rep(.6, 6), rep(0, 6),
                                            rep(0, 6), rep(0, 6), rep(.6, 6),
                                            rep(.6, 6), rep(0, 12),
                                            rep(0, 6), rep(.6, 6), rep(0, 6),
                                            rep(0, 6), rep(0, 6), rep(.6, 6),
                                            rep(0, 6), rep(.6, 6), rep(0, 6),
                                            rep(.6, 6), rep(0, 12),
                                            rep(0, 6), rep(0, 6), rep(.6, 6)))> .3),
                                          c(18, 3, 3)),
                         phi = NA, extract_phi = FALSE, n_factors = 3)

arr_re <- .array_reorder(L = array(c(rep(.6, 6), rep(0, 12),
                                        rep(0, 6), rep(.6, 6), rep(0, 6),
                                        rep(0, 6), rep(0, 6), rep(.6, 6),
                                        rep(.6, 6), rep(0, 12),
                                        rep(0, 6), rep(.6, 6), rep(0, 6),
                                        rep(0, 6), rep(0, 6), rep(.6, 6),
                                        rep(0, 6), rep(.6, 6), rep(0, 6),
                                        rep(-.6, 6), rep(0, 12),
                                        rep(0, 6), rep(0, 6), rep(.6, 6)),
                                      c(18, 3, 3)),
                         vars_accounted = array(c(rep(.2, 3),
                                                  rep(.3, 3),
                                                  rep(.4, 3),
                                                  rep(.2, 3),
                                                  rep(.3, 3),
                                                  rep(.4, 3),
                                                  rep(.3, 3),
                                                  rep(.2, 3),
                                                  rep(.4, 3)),
                                                c(3, 3, 3)),
                            L_corres = array(as.numeric(abs(c(rep(.6, 6), rep(0, 12),
                                                              rep(0, 6), rep(.6, 6), rep(0, 6),
                                                              rep(0, 6), rep(0, 6), rep(.6, 6),
                                                              rep(.6, 6), rep(0, 12),
                                                              rep(0, 6), rep(.6, 6), rep(0, 6),
                                                              rep(0, 6), rep(0, 6), rep(.6, 6),
                                                              rep(0, 6), rep(.6, 6), rep(0, 6),
                                                              rep(.6, 6), rep(0, 12),
                                                              rep(0, 6), rep(0, 6), rep(.6, 6)))> .3),
                                             c(18, 3, 3)),
                            phi = array(rep(c(1, .3, .4, .3, 1, .2, .4, .2, 1), 3), c(3, 3, 3)),
                         extract_phi = TRUE, n_factors = 3)
test_that(".array_reorder works", {
  ### tests for arr_re_NA with phi = NA and extract_phi = FALSE
  expect_type(arr_re_NA, "list")
  expect_named(arr_re_NA, c("L", "L_corres", "phi", "vars_accounted"))
  checkmate::expect_array(arr_re_NA$L)
  expect_equal(dim(arr_re_NA$L), c(18, 3, 3))
  expect_equal(arr_re_NA$L,
               array(c(rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6),
                       rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6),
                       rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6)),
                      c(18, 3, 3)))
  checkmate::expect_array(arr_re_NA$L_corres)
  expect_equal(dim(arr_re_NA$L_corres), c(18, 3, 3))
  expect_equal(arr_re_NA$L_corres,
               array(as.numeric(c(rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6),
                       rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6),
                       rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6)) > .3),
                     c(18, 3, 3)))
  expect_equal(arr_re_NA$phi, NA)
  checkmate::expect_array(arr_re_NA$vars_accounted)
  expect_equal(arr_re_NA$vars_accounted,
               array(c(rep(.2, 3),
                       rep(.3, 3),
                       rep(.4, 3),
                       rep(.2, 3),
                       rep(.3, 3),
                       rep(.4, 3),
                       rep(.2, 3),
                       rep(.3, 3),
                       rep(.4, 3)),
                     c(3, 3, 3)))

  ### tests for arr_re with phi = array() and extract_phi = TRUE
  expect_type(arr_re, "list")
  expect_named(arr_re, c("L", "L_corres", "phi", "vars_accounted"))
  checkmate::expect_array(arr_re$L)
  expect_equal(dim(arr_re$L), c(18, 3, 3))
  expect_equal(arr_re$L,
               array(c(rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6),
                       rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6),
                       rep(.6, 6), rep(0, 12),
                       rep(0, 6), rep(.6, 6), rep(0, 6),
                       rep(0, 6), rep(0, 6), rep(.6, 6)),
                     c(18, 3, 3)))
  checkmate::expect_array(arr_re$L_corres)
  expect_equal(dim(arr_re$L_corres), c(18, 3, 3))
  expect_equal(arr_re$L_corres,
               array(as.numeric(c(rep(.6, 6), rep(0, 12),
                                  rep(0, 6), rep(.6, 6), rep(0, 6),
                                  rep(0, 6), rep(0, 6), rep(.6, 6),
                                  rep(.6, 6), rep(0, 12),
                                  rep(0, 6), rep(.6, 6), rep(0, 6),
                                  rep(0, 6), rep(0, 6), rep(.6, 6),
                                  rep(.6, 6), rep(0, 12),
                                  rep(0, 6), rep(.6, 6), rep(0, 6),
                                  rep(0, 6), rep(0, 6), rep(.6, 6)) > .3),
                     c(18, 3, 3)))
  checkmate::expect_array(arr_re$phi)
  expect_equal(dim(arr_re$phi), c(3, 3, 3))
  expect_equal(arr_re$phi,
               array(c(1, .3, .4, .3, 1, .2, .4, .2, 1,
                       1, .3, .4, .3, 1, .2, .4, .2, 1,
                       1, -.3, -.2, -.3, 1, .4, -.2, .4, 1), c(3, 3, 3)))
  checkmate::expect_array(arr_re$vars_accounted)
  expect_equal(arr_re$vars_accounted,
               array(c(rep(.2, 3),
                       rep(.3, 3),
                       rep(.4, 3),
                       rep(.2, 3),
                       rep(.3, 3),
                       rep(.4, 3),
                       rep(.2, 3),
                       rep(.3, 3),
                       rep(.4, 3)),
                     c(3, 3, 3)))

})

test_that(".array_reorder uses an optimal permutation alignment", {
  # Two solutions engineered so that a greedy per-row congruence match maps both
  # target factors 1 and 3 onto source column 1 (order 1, 2, 1): column 1 would
  # be duplicated and column 3 dropped, corrupting the average. The optimal
  # (linear-sum-assignment) alignment instead recovers the true permutation.
  L1 <- cbind(c(.9, .8, 0, 0, 0, 0),
              c(0, 0, .9, .8, 0, 0),
              c(0, 0, 0, 0, .9, .8))
  Ln <- cbind(c(.2, .1, 0, 0, .9, .8),
              c(0, 0, .9, .8, 0, 0),
              c(0, 0, 0, 0, .2, .1))

  # The greedy which.max alignment used previously collapses to a non-permutation.
  greedy_order <- apply(abs(.factor_congruence(L1, Ln, skip_checks = TRUE)), 1,
                        which.max)
  expect_equal(greedy_order, c(1, 2, 1))

  L_arr <- array(c(L1, Ln), dim = c(6, 3, 2))
  corres_arr <- array(as.numeric(abs(c(L1, Ln)) > .3), dim = c(6, 3, 2))
  va_arr <- array(1, dim = c(3, 3, 2))

  res <- .array_reorder(vars_accounted = va_arr, L = L_arr, L_corres = corres_arr,
                        phi = NA, extract_phi = FALSE, n_factors = 3)

  aligned <- res$L[, , 2]
  # All three source factors are retained (no duplicated/dropped column) ...
  expect_equal(nrow(unique(round(t(aligned), 8))), 3)
  # ... and the recovered order is the true permutation, leaving Ln in place.
  expect_equal(aligned, Ln)
})


rm(arr_re_NA, arr_re)
