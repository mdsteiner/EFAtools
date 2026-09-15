## Summary

This release fixes statistical and input-validation defects. NEWS.md lists
all relevant changes.

This release also fixes the test failures that CRAN currently reports for
macOS x86_64. The old tests checked implementation details of the
polychoric/tetrachoric ACOV code. The new tests check the statistical result
instead.

## Test environments

* Windows 11, R 4.6.0 (local)
* win-builder: oldrel, release, devel
* GitHub Actions: ubuntu (oldrel, release, devel, full test suite), windows,
  macOS (incl. macos-15-intel, vecLib)
* R-hub: linux, windows, macos, macos-arm64, atlas, mkl, nold, clang-asan,
  gcc-asan, clang-ubsan

## R CMD check results

0 errors | 0 warnings | 0 notes
