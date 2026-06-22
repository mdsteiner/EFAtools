# Opt-in gate for slow tests. The heavy retention (HULL, PARALLEL, NEST, N_FACTORS),
# np-boot, polychoric, multimodal-rotation, lavaan-oracle, and grid-of-EFAs fixtures
# dominate the suite's runtime; they are skipped by default and run only when
# EFATOOLS_TEST_SLOW is set truthy, e.g. `Sys.setenv(EFATOOLS_TEST_SLOW = "true")` (or
# `EFATOOLS_TEST_SLOW=true` in the shell) before devtools::test() / devtools::check().
# `is_slow_test()` is the predicate the gate is built on, so file-top fixture blocks can
# guard with `if (is_slow_test()) { ... }` to avoid building heavy fixtures by default.
is_slow_test <- function() {
  isTRUE(as.logical(Sys.getenv("EFATOOLS_TEST_SLOW")))
}
skip_if_not_slow <- function() {
  testthat::skip_if_not(
    is_slow_test(),
    "slow test; set EFATOOLS_TEST_SLOW=true to run"
  )
}
