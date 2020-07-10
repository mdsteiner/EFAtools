## Resubmission
This is a resubmission. In this version we have:

* Deleted the packages "parallel" and "plotly" from Imports as they are no
  longer used.  Moreover, we relocated the package "microbenchmark" from Imports
  to Suggest, as it is only used in the vignette.

* fixed a bug in parallel_summarise to fix the 'clang-ASAN', 'gcc-ASAN', 'noLD', and 
  'valgrind' issues.
  
* Added an error message in PARALLEL if no solution has been found after 25 tries.

* Fixed a bug in PARALLEL in method "EFA" 

* Updated the readme with CRAN installation instructions

* Updated the EFA-test that threw an error on some systems, fixing 'OpenBLAS' issue.

* Updated input to an N_FACTORS-test to ensure an exactly singular matrix is
  entered and the expected error is thrown, fixing 'noLD' issue.
  
* Updated an OMEGA-helper-test to avoid 'MKL' issue.



## Test environments
* local OS X R installation (Mojave 10.14.6), R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.0
* win-builder (release and devel)
* Windows Server 2008 R2 SP1 (r-hub), R-devel, 32/64 bit

## R CMD check results

0 errors | 0 warnings | 1 note

* Days since last update: 3 -> this is a patch fixing issues that were only
  discovered in the extended CRAN Package Checks after the initial release.
