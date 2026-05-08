## Resubmission
This is a resubmission. In this version we have:

* Fixed a test to accomodate changes in an upcoming release of the psych package, which caused reverse-dependency issues in psych.
* `EFA()`: Added `randomStarts` argument passed to GPArotation functions, as suggested by Coen Bernaards.
* `FACTOR_SCORES()`: Added `rho` argument, thanks to Andreas Soteriades.
* `EFA_POOLED()`: Fixed issue that could lead to averaged Phi not being symmetric.

## Test environments
* mac-builder (release, devel)
* local Windows 11 installation, R 4.6.0
* win-builder (release, devel, and oldrelease)

## R CMD check results

0 errors | 0 warnings | 1 note

Notes:

Namespace in Imports field not imported from: ‘progress’
-> package is needed in EFA_AVERAGE() (accessed indirectly via progressr package)

