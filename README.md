# EFAtools

**Note:** The EFAtools package is still in development.

The EFAtools package provides functions to perform exploratory factor analysis (EFA) procedures and compare their solutions. The goal is to provide state-of-the-art factor retention methods and a high degree of flexibility in the EFA procedures. This way, implementations from R psych and SPSS can be compared. Moreover, functions for Schmid-Leiman transformation, and computation of omegas are provided. To speed up the analyses, some of the iterative procedures like principal axis factoring (PAF) are implementded in C++.

## Installation

You can install EFAtools with:

``` r
install.packages("devtools")
devtools::install_github("mdsteiner/EFAtools")
```

## Example

Here are a few examples on how to perform the analyses with the different types and how to compare the results using the `compare` function.

``` r
# Run all possible factor retention methods
nfac_all <- N_FACTORS(test_models$baseline$cormat, N = 500, method = "ML")

# A type SPSS EFA to mimick the SPSS implementation with
# promax rotation
EFA_SPSS <- EFA(test_models$baseline$cormat, n_factors = 3, type = "SPSS",
                  rotation = "promax")

# A type psych EFA to mimick the psych::fa() implementation with
# promax rotation
EFA_psych <- EFA(test_models$baseline$cormat, n_factors = 3, type = "psych",
                  rotation = "promax")

# compare the type psych and type SPSS implementations
compare(EFA_SPSS$rot_loadings, EFA_psych$rot_loadings)

# Perform a Schmid-Leiman transformation
SL <- SL(EFA_psych)

# Compute omegas of the Schmid-Leiman solution
OMEGA(SL, factor_corres = rep(c(3, 2, 1), each = 6))
```

