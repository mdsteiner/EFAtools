## Resubmission
This is a resubmission. In this version we have:

* Fixed a bug in the `EKC()` function.

## Test environments
* mac-builder m1, R 4.4.2
* local Windows 11 installation, R 4.5.1
* win-builder (release, devel, and oldrelease)

## R CMD check results

0 errors | 0 warnings | 1 note

Notes:

Namespace in Imports field not imported from: ‘progress’
-> package is needed in EFA_AVERAGE() (accessed indirectly via progressr package)

