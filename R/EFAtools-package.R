#' @keywords internal
#' @useDynLib EFAtools, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils modifyList
## `lifecycle` is a run-time dependency, not a documentation-only one: `EFA()`,
## `efa_average()` and `efa_mi()` give their legacy arguments
## `lifecycle::deprecated()` defaults and resolve them with
## `lifecycle::is_present()` / `lifecycle::deprecate_warn()`. It is additionally
## used at documentation build time by the `lifecycle::badge()` calls that mark the
## superseded functions and the deprecated argument. Every call is qualified, so no
## import is declared.
"_PACKAGE"
