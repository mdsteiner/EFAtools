# EFAtools 0.1.1.9000

## General

* Created new vignette *Replicate_SPSS_psych* to show replication of original `psych` and `SPSS` EFA solutions with `EFAtools`.

## Changes to Functions

* `EFA()`: Updated the EFAtools type in PAF and Pomax.
* `N_FACTORS()`: Added option to do a scree plot if "SCREE" is included in the `criteria` argument.

## Bug Fixes

* `PARALLEL()`: Fixed a bug that occurred when using `decision_rule = "percentile"`

## New Functions

* Added function `SCREE()` that does a scree plot.


# EFAtools 0.1.1

## Minor Changes

* Added an error message in PARALLEL if no solution has been found after 25 tries.

## Bug Fixes

* Updated different tests

* Deleted no longer used packages from Imports and Suggests in DESCRIPTION

* `PARALLEL()`: fixed a bug in indexing if method `"EFA"` was used.


# EFAtools 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* Initial CRAN submission
