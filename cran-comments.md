## Resubmission
This is a resubmission. In this version we have:

* Added new functions:
  * `MAP()` computes the Velicer MAP criterion (both TR2 and TR4).
  * `PROCRUSTES()` to perform orthogonal and oblique Procrustes / target rotation.
  * `CONSENSUS_PROCRUSTES()` to perform Procrustes on a list of targets to obtain a common target.
  * `residuals.EFA()` to extract and print residuals and, if computed, standardized residuals.
  * `EFA_POOLED()` to run EFA on a list of multiple imputated datasets and pool the results.
  * `print.EFA_POOLED()`, print method adapted from `print.EFA()`
* Added new features to existing Functions:
  * `N_FACTORS()` can also compute the MAP criterion.
  * `EFA()`:
    * Now returns and prints residuals and, if SEs are computed, standardized residuals.
    * Calculates and prints RMSR.
    * Can calculate bootstrap standard errors and CIs of parameters and fit indices.
  * `print.EFA()` now prints more information.
* Implemented the following bug fixes:
  * Fixed a bug in `OMEGA()` that led to incorrect omega, H, and ECV values for `lavaan` bifactor models. Tnanks to Christopher D. King for bug report and suggested fix.
  * Small fix in the documentation of `EFA_AVERAGE()` and added sign reflection in factor correlations if necessary.
  * Fixed incorrect calculation of RMSEA in `EFA()`.
  * Fixed a bug in the `EKC()` function.

## Test environments
* mac-builder (release, devel)
* local Windows 11 installation, R 4.6.0
* win-builder (release, devel, and oldrelease)

## R CMD check results

0 errors | 0 warnings | 1 note

Notes:

Namespace in Imports field not imported from: ‘progress’
-> package is needed in EFA_AVERAGE() (accessed indirectly via progressr package)

