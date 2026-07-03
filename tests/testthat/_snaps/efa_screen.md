# print output is stable

    Code
      print(scr_cor)
    Output
      
      -- Sampling adequacy and sphericity --------------------------------------------
      
      v The overall KMO value for your data is marvellous (Overall KMO = <num>).
      These data are probably suitable for factor analysis.
      
      v The Bartlett's test of sphericity was significant at an alpha level of <num>.
      These data are probably suitable for factor analysis.
      𝜒²(153) = <num>, p < <num>
      
      -- Multicollinearity -----------------------------------------------------------
      
      v Determinant: <num>. No concern (a value near 0 signals multicollinearity).
      v Condition number: <num>. No concern (large values signal near-collinear variables).
      
      -- Per-variable diagnostics ----------------------------------------------------
      
            MSA   SMC
      V1 <num> <num>
      V2 <num> <num>
      V3 <num> <num>
      V4 <num> <num>
      V5 <num> <num>
      V6 <num> <num>
      V7 <num> <num>
      V8 <num> <num>
      V9 <num> <num>
      V10 <num> <num>
      V11 <num> <num>
      V12 <num> <num>
      V13 <num> <num>
      V14 <num> <num>
      V15 <num> <num>
      V16 <num> <num>
      V17 <num> <num>
      V18 <num> <num>
      
      -- Recommendations -------------------------------------------------------------
      
      v The data appear suitable for factor analysis.
      i Per-item variance, missing-data, category, normality, and outlier diagnostics
        require raw data; only a correlation matrix was supplied.

---

    Code
      print(scr_raw)
    Output
      
      -- Sampling adequacy and sphericity --------------------------------------------
      
      v The overall KMO value for your data is marvellous (Overall KMO = <num>).
      These data are probably suitable for factor analysis.
      
      v The Bartlett's test of sphericity was significant at an alpha level of <num>.
      These data are probably suitable for factor analysis.
      𝜒²(28) = <num>, p < <num>
      
      -- Multicollinearity -----------------------------------------------------------
      
      v Determinant: <num>. No concern (a value near 0 signals multicollinearity).
      v Condition number: <num>. No concern (large values signal near-collinear variables).
      
      -- Per-variable diagnostics ----------------------------------------------------
      
                variance missing   SMC   MSA  flags
      fun <num>       0 <num> <num>       
      friends <num>       0 <num> <num>       
      enjoy <num>       0 <num> <num> sparse
      hurt <num>       0 <num> <num> sparse
      part <num>       0 <num> <num>       
      commonly <num>       0 <num> <num>       
      chances <num>       0 <num> <num>       
      attracted <num>       0 <num> <num> sparse
      
      -- Multivariate normality ------------------------------------------------------
      
      x Mardia's skewness: 𝜒²(120) = <num>, p < <num>.
      x Mardia's kurtosis: z = <num>, p < <num>.
      x Henze-Zirkler: HZ = <num>, p < <num>.
      These data depart from multivariate normality.
      
      -- Outliers --------------------------------------------------------------------
      
      ! A robust (MCD) covariance could not be computed; classical Mahalanobis distances were used.
      i 71 of 810 observations were flagged as multivariate outliers (Mahalanobis distance > <num>).
      
      -- Recommendations -------------------------------------------------------------
      
      ! These data depart from multivariate normality; normal-theory standard errors
        and fit statistics may be biased - prefer robust (sandwich) or bootstrapped
        standard errors.
      ! Bartlett's test is significant, but it assumes multivariate normality and
        grows more sensitive as N increases; because these data are non-normal, treat
        it as uninformative here and rely on the KMO.
      ! 3 variables have a sparse response category (< 5 responses): enjoy, hurt, and
        attracted; a low-frequency category can destabilise polychoric estimates -
        consider collapsing it into an adjacent category.
      ! 71 observations were flagged as potential multivariate outliers; inspect them
        (see `$outliers$flagged`) before down-weighting or excluding.

---

    Code
      print(scr_non)
    Output
      
      -- Sampling adequacy and sphericity --------------------------------------------
      
      v The overall KMO value for your data is marvellous (Overall KMO = <num>).
      These data are probably suitable for factor analysis.
      
      ! Bartlett's test of sphericity was not computed; no sample size (N) was supplied.
      
      -- Multicollinearity -----------------------------------------------------------
      
      v Determinant: <num>. No concern (a value near 0 signals multicollinearity).
      v Condition number: <num>. No concern (large values signal near-collinear variables).
      
      -- Per-variable diagnostics ----------------------------------------------------
      
            MSA   SMC
      V1 <num> <num>
      V2 <num> <num>
      V3 <num> <num>
      V4 <num> <num>
      V5 <num> <num>
      V6 <num> <num>
      V7 <num> <num>
      V8 <num> <num>
      V9 <num> <num>
      V10 <num> <num>
      V11 <num> <num>
      V12 <num> <num>
      V13 <num> <num>
      V14 <num> <num>
      V15 <num> <num>
      V16 <num> <num>
      V17 <num> <num>
      V18 <num> <num>
      
      -- Recommendations -------------------------------------------------------------
      
      v The data appear suitable for factor analysis.
      i Per-item variance, missing-data, category, normality, and outlier diagnostics
        require raw data; only a correlation matrix was supplied.

