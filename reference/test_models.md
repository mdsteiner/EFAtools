# Four test models used in Grieder and Steiner (2022)

Correlation matrices created from simulated data from four of the
[`population_models`](https://mdsteiner.github.io/EFAtools/reference/population_models.md)
cases, each with strong factor intercorrelations. These are used in
Grieder & Steiner (2022) to compare the psych and SPSS implementations
in this package with the actual implementations of the programs. For
details on the cases, see
[`population_models`](https://mdsteiner.github.io/EFAtools/reference/population_models.md).

## Usage

``` r
test_models
```

## Format

A list of 4 lists "baseline", "case_1a", "case_6b", and "case_11b", each
with the following elements.

- cormat:

  (matrix) - The correlation matrix of the simulated data.

- n_factors:

  (numeric) - The true number of factors.

- N:

  (numeric) - The sample size of the generated data.

## Source

Grieder, S., & Steiner, M. D. (2022). Algorithmic jingle jungle: A
comparison of implementations of principal axis factoring and promax
rotation in R and SPSS. Behavior Research Methods, 54, 54–74. doi:
10.3758/s13428-021-01581-x
