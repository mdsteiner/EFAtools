# Opt-in gate for slow tests. The heavy retention (HULL, PARALLEL, NEST, efa_retain),
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

# Counts the gated blocks in the test files themselves, so the number the notice below
# reports cannot drift away from the suite. Returns 0 when called from outside the test
# directory, where the files are not visible.
n_slow_gated <- function() {
  sum(vapply(
    list.files(pattern = "^test-.*\\.[Rr]$"),
    function(path) {
      sum(grepl("skip_if_not_slow()", readLines(path, warn = FALSE), fixed = TRUE))
    },
    integer(1)
  ))
}

# A skip is not a failure, so a default run can look complete when it is not. State once
# what the gate holds back -- a fact about the suite, not about the blocks this particular
# run happens to touch. local() keeps the count out of the environment the tests run in.
if (!is_slow_test()) local({
  n_gated <- n_slow_gated()
  if (n_gated > 0) {
    cli::cli_inform(c("i" = paste(
      "The suite has {n_gated} slow test block{?s} gated off; set",
      "{.envvar EFATOOLS_TEST_SLOW}=true (with {.envvar NOT_CRAN}=true) to run them."
    )))
  }
})
