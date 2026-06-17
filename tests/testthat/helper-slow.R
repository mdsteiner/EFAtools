# Opt-in gate for the slow rotation tests. The multimodal simplimax criterion fully optimizes
# every random start (the smooth criteria only screen the starts cheaply), so its high-restart
# parity checks dominate the suite's runtime. They are skipped by default and run only when
# EFATOOLS_TEST_SLOW is set truthy, e.g. `Sys.setenv(EFATOOLS_TEST_SLOW = "true")` before
# devtools::test() or devtools::check().
skip_if_not_slow <- function() {
  testthat::skip_if_not(
    isTRUE(as.logical(Sys.getenv("EFATOOLS_TEST_SLOW"))),
    "slow rotation test; set EFATOOLS_TEST_SLOW=true to run"
  )
}
