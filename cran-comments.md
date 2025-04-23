## Resubmission
This is a resubmission. In this version we have:

* Fixed a test of the .gof function, to take care of the ATLAS issue in version 0.4.5 of the package (the respective goodness of fit measure is set to the worst value and a warning is issued, if the relevant matrix cannot be inverted).
* Updated a url in the main vignette to a doi-url.

## Test environments
* mac-builder m1, R 4.4.2
* local Windows 10 installation, R 4.4.3
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
