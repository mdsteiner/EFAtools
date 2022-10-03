## Resubmission
This is a resubmission. In this version we have:

* Added a small adaptation to .gof() that should take care of the MKL issue when using R-devel on x86_64 Fedora 34 Linux with alternative BLAS/LAPACK.

## Test environments
* mac OS 11.7 (via GitHub Actions), R 4.2.1
* mac-builder m1, R 4.2.1
* ubuntu 20.04 (via GitHub Actions), R 4.2.1
* Microsoft Windows Server 2022 (via GitHub Actions), R 4.2.1
* local Windows 10 installation, R 4.2.1
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
