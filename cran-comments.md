## Resubmission
This is a resubmission. In this version we have:

* Fixed the `print.EFA()` function to return explained variances and sum of squared loadings based on the rotated rather than unrotated factor solution, if a rotation was performed.
* Added communalities and uniquenesses to the print function `print.EFA()`.
* Calculate and return the model implied correlation matrix in `EFA()`.

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
