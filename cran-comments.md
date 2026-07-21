## Summary

This is a major release. The public interface of the package is now the
lowercase `efa_*` function family (`efa_fit()`, `efa_retain()`, `efa_screen()`,
and friends), configured through the new control objects `estimate_control()`
and `rotate_control()`. The release also adds a vignette on migrating to the new
interface and one on EFA with ordinal and missing data, and rewrites the
workflow vignette.

The uppercase function names of earlier versions (`EFA()`, `N_FACTORS()`,
`CD()`, `HULL()`, and the others) are superseded by their `efa_*` equivalents,
but they remain exported with unchanged arguments and emit no deprecation
warning. Existing code therefore needs no changes, and neither do the reverse
dependencies.

The submission also includes various bug-fixes, including the 
additional issues raised in the CRAN tests.

## Test environments

* local Windows 11 installation, R 4.6.0
* win-builder (release, devel, and oldrelease)
* mac-builder (release)
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release,
  and oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes
