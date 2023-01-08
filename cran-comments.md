## Resubmission
This is a resubmission. In this version we have:

* Added a small adaptation in a test of .gof() that should take care of the ATLAS issue when using R-devel on x86_64 Fedora 36 Linux with alternative BLAS/LAPACK.
* added `\dontrun{}` to examples of `EFA_AVERAGE()` and its print and plot methods as these were causing issues on R-oldrel which were not directly related to EFAtools and thus could not be fixed from within the package.
* `print.EFA()`: Added arguments `cutoff`, `digits` and `max_name_length` that are passed to `print.LOADINGS()`.
* `print.LOADINGS()`: New Argument `max_name_length` to control the maximum length of the displayed variable names (names longer than this will be cut on the right side). Previously, this was fixed to 10 (which is now the default).

## Test environments
* mac OS 11.7 (via GitHub Actions), R 4.2.1
* mac-builder m1, R 4.2.1
* ubuntu 20.04 (via GitHub Actions), R 4.2.1
* Microsoft Windows Server 2022 (via GitHub Actions), R 4.2.1
* local Windows 10 installation, R 4.2.2
* win-builder (release, devel, and oldrelease)

## R CMD check results

0 errors | 0 warnings | 2 notes

Notes:

Namespace in Imports field not imported from: ‘progress’
-> package needed in EFA_AVERAGE() (accessed indirectly via progressr package)

Found the following (possibly) invalid URLs:
  URL: https://psycnet.apa.org/record/2019-02883-001
    From: inst/doc/EFAtools.html
    Status: 403
    Message: Forbidden
-> I double checked the link and it works fine
