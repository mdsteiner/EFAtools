---
title: 'EFAtools: An R package with fast and flexible implementations of exploratory
  factor analysis tools'
tags:
- R
- exploratory factor analysis
- factor retention methods
- hierarchical factor analysis
- comparison of implementations
date: "19 July 2020"
output:
  html_document:
    df_print: paged
authors:
- name: Markus D. Steiner
  orcid: 0000-0002-8126-0757
  affiliation: 1
- name: Silvia Grieder
  orcid: 0000-0002-0118-7722
  affiliation: 2
bibliography: paper.bib
affiliations:
- name: Center for Cognitive and Decision Sciences, Department of Psychology, University
    of Basel, Switzerland
  index: 1
- name: Division of Developmental and Personality Psychology, Department of Psychology,
    University of Basel, Switzerland
  index: 2
---

# Summary

In the social sciences, factor analysis is a widely used tool to identify latent constructs underlying task performance or the answers to questionnaire items. Exploratory factor analysis (EFA) is a data-driven approach to factor analysis and is used to extract a smaller number of common factors that represent or explain the common variance of a larger set of manifest variables [see @costello2005best; and @watkins2018exploratory for an overview]. Several decisions have to be made in advance when performing an EFA, including the number of factors to extract, and the extraction and rotation method to be used. After a factor solution has been found, especially for data structures in the field of intelligence research where usually high, positive factor intercorrelations occur, it is useful to subject the resulting factor solution to an orthogonalization procedure to achieve a hierarchical factor solution with one general and several specific factors. From this orthogonalized, hierarchical solution, the variance can then be partitioned to estimate the relative importance of the general versus the specific factors using omega reliability coefficients [e.g., @mcdonald1999test].

*EFAtools* is an R package [@R2018R] that enables fast and flexible analyses in an EFA framework, from tests for suitability of the data for factor analysis and factor retention criteria to hierarchical factor analysis with Schmid-Leiman transformation [@schmid1957development] and McDonald's omegas [e.g., @mcdonald1999test]. The package's core functionalities are listed in Table 1. 

# Statement of Need

Compared to other R packages with which EFA can be performed, *EFAtools* has several advantages, including fast implementations using *Rcpp* [@eddelbuettel2017extending; @eddelbuettel2014rcpparmadillo], more flexibility in the adjustment of implementation features, the ability to reproduce the R *psych* [@revelle2019psych] and SPSS [@ibm_spss_2015] implementations of some analyses methods (see vignette *Replicate SPSS and R psych results with EFAtools*), as well as the inclusion of recommended implementations for these methods based on simulation analyses by the authors (pending publication). Finally, the package includes the implementation of the, as of yet, most comprehensive set of factor retention criteria in R, including recently developed criteria such as the Hull method [@lorenzoseva_hull_2011], comparison data [@ruscio_determining_2012], and the empirical Kaiser criterion [@braeken_empirical_2016]. As recommended by @auerswald_how_2019, multiple factor retention criteria should be examined simultaneously to check their convergence, which now is easily possible with a comprehensive function in *EFAtools* incorporating all implemented factor retention criteria for simultaneous application. Minor advantages over and above the existing implementations in R include that when intending to perform a Schmid-Leiman transformation, this can be done on an obliquely rotated solution obtained with functions from the *EFAtools* or the *psych* package instead of being forced to perform the whole EFA procedure again. Moreover, our implementation of McDonald's omegas calculations include the possibility of manual variable-to-factor correspondences (as are needed for variance partitioning for predetermined / theoretical composites) in addition to automatically determined variable-to-factor correspondences (as done, for example, in the *psych* package). Further, the *EFAtools* function to compute McDonald's omegas can easily be applied on *EFAtools* and *psych* Schmid-Leiman solutions as well as on *lavaan* [@rosseel_lavaan_2012] second-order, bifactor, and single factor solutions (including solutions from multiple group analyses).

# Development and Purpose

*EFAtools* was designed for use in the social sciences in general and is especially suitable for research on cognitive abilities or other hierarchically organized constructs as well as for more time-consuming applications such as in simulation analyses. Its development arose from the need for a tool for easy replication and comparison of EFA solutions from different programs, namely R and SPSS (pending publication), and has already been used in another publication [@grieder_exploratory_2019]. The package was then expanded for a broader, easy, fast, and flexible use of EFA tools such that it is now suitable for most projects within the EFA framework.


Table: Core functionalities of *EFAtools*.

| Topic                | Method                     | Function   
|----------------------|----------------------------|---------------|
|Suitability for factor analysis | Bartlett's test of sphericity | `BARTLETT()` |
|                                | Kaiser-Meyer-Olkin criterion | `KMO()` |
|Factor retention criteria | Comparison data                    | `CD()` |
|                          | Empirical Kaiser criterion         | `EKC()` |
|                          | Hull method                        | `HULL()` |
|                          | Kaiser-Guttman criterion           | `KGC()` |
|                          | Parallel analysis                  | `PARALLEL()` |
|                          | Scree plot                         | `SCREE()` |
|                          | Sequential model tests             | `SMT()` |
|                          | RMSEA lower bound criterion        | `SMT()` |
|                          | AIC criterion                      | `SMT()` |
|Factor extraction methods | Principal axis factoring           | `EFA()` |
|                          | Maximum likelihood                 | `EFA()` |
|                          | Unweighted least squares           | `EFA()` |
|Rotation methods | Orthogonal: Varimax, equamax, quartimax, geominT, bentlerT, bifactorT | `EFA()` |
|                 | Oblique: Promax, oblimin, quartimin, simplimax, bentlerQ, geominQ, bifactorQ | `EFA()` |
|Hierarchical factor analysis | Schmid-Leiman transformation   | `SL()` |
|                          | McDonald's omegas                 | `OMEGA()` |
*Note*. All functions for suitability for factor analysis and factor retention criteria can be called in any desired combination using the `N_FACTORS()` function.


# Installation

The *EFAtools* package can be installed from CRAN using `install.packages("EFAtools")`. Moreover, the development version can be installed from GitHub (https://github.com/mdsteiner/EFAtools) using `devtools::install_github("mdsteiner/EFAtools", build_vignettes = TRUE)`.

# Acknowledgements

We thank Dirk Wulff for helpful suggestions concerning the C++ implementations.

# References
