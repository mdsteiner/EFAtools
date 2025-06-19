## Resubmission
This is a resubmission. In this version we have:

* Added `NEST()` to perform the Next Eigenvalue Sufficiency Test (Achim, 2017).
* Updated `EKC()` to include both circulating implementation of the empirical Kaiser criterion.
* `N_FACTORS()`: 
  * Updated default settings to more sensible values.
  * Added NEST as additional factor retention method.
  * Added arguments to control NEST and to reflect changes in EKC.

## Test environments
* mac-builder m1, R 4.4.2
* local Windows 11 installation, R 4.5.0
* win-builder (release, devel, and oldrelease)

## R CMD check results

0 errors | 0 warnings | 2 notes

Notes:

Namespace in Imports field not imported from: ‘progress’
-> package is needed in EFA_AVERAGE() (accessed indirectly via progressr package)

Found the following (possibly) invalid URLs:
  URL: https://osf.io/x5cz2/?view_only=d03efba1fd0f4c849a87db82e6705668
    From: man/CD.Rd
    Status: Error
    Message: SSL connect error [osf.io]: schannel: failed to receive handshake, SSL/TLS connection failed
-> I double checked the link and it works fine
