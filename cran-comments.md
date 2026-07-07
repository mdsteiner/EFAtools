## Resubmission
This is a resubmission. In this version we have:

Added new functions:
* `efa_group()` performs multi-group EFA.
* `efa_power()` to conduct power analysis (analytic power based on RMSEA and simulation based power analysis).
* `efa_reliability()` to calculate various reliability indices.
* `efa_screen()` to screen data for multivariate normality and suitability for factor analysis.
* `efa_simulate()` to simulate data with various distributions and missing data mechanisms.

Changes to functions:
* computation of polychoric correlations
* raw-data missing value handling using FIML estimation
* added sandwich and information standard errors
* implemented all criterion-based rotations internally in C++ and moved GPArotation from imports to suggests (only used in tests)
* various bug-fixes

## Test environments
* mac-builder (release, devel)
* local Windows 11 installation, R 4.6.0
* win-builder (release, devel, and oldrelease)

## R CMD check results

0 errors | 0 warnings | 0 notes


