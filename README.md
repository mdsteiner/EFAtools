# EFAdiff

The EFAdiff package provides functions to perform and compare exploratory factor analysis (EFA) implementations from R psych and SPSS, MacOrtho and Omega. It provides functions for principal axis factoring (PAF), promax or varimax rotation, Schmid-Leiman transformation, and computation of omegas.

## Installation

You can install EFAdiff with:

``` r
install.packages("devtools")
devtools::install_github("mdsteiner/EFAdiff")
```

## Example

Here are a few examples on how to perform the analyses with the different types and how to compare the results using the `compare` function.

``` r
# A type SG (as presented in Steiner and Grieder, 2019) EFA with
# promax rotation
EFA_SG_5 <- EFA(IDS2_R, n_factors = 5, type = "SG", rotation = "promax")

# A type SPSS EFA to mimick the SPSS implementation with
# promax rotation
EFA_SPSS_5 <- EFA(IDS2_R, n_factors = 5, type = "SPSS", rotation = "promax")

# A type psych EFA to mimick the psych::fa() implementation with
# promax rotation
EFA_psych_5 <- EFA(IDS2_R, n_factors = 5, type = "psych", rotation = "promax")

# compare the type psych and type SPSS implementations
compare(EFA_SPSS_5$unrot_loadings, EFA_psych_5$unrot_loadings)

# Perform a Schmid-Leiman transformation of the type SG EFA
SL_SG_5 <- SL(EFA_SG_5, type = "SG")

# Compute omegas of the Schmid-Leiman solution
OMEGA_SG_SG_5 <- OMEGA(SL_SG_SG_5,
                       factor_corres = c(4, 4, 5, 5, 2, 2, 3, 3, 4, 4, 1, 1,1, 1),
                       cormat = IDS2_R,
                       type = "SG")
```

