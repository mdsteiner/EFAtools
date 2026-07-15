#' @keywords internal
#' @useDynLib EFAtools, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils modifyList
## `lifecycle` is used only at documentation build time (the superseded badges in
## R/EFAtools-superseded.R); importing `deprecated` keeps it a used dependency so
## `R CMD check` does not flag it as an unused Import.
#' @importFrom lifecycle deprecated
"_PACKAGE"
