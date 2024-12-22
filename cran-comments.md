## Resubmission
This is a resubmission. In this version we have:

* Fixed an error that occurred due to an update in the psych R-Package.
* Updated CITATION to include bibenty() rather than citeEntry()
* Removed C++11 standard

## Test environments
* mac-builder m1, R 4.4.0
* local Windows 10 installation, R 4.4.0
* win-builder (release, devel, and oldrelease)

## R CMD check results

0 errors | 0 warnings | 2 notes

Notes:

Namespace in Imports field not imported from: ‘progress’
-> package is needed in EFA_AVERAGE() (accessed indirectly via progressr package)

Found the following (possibly) invalid URLs:
  URL: https://psycnet.apa.org/record/2019-02883-001
    From: inst/doc/EFAtools.html
    Status: 403
    Message: Forbidden
-> I double checked the link and it works fine
