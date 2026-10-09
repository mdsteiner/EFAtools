# print output is stable

    Code
      print(scr_cor)
    Output
      
      -- Sampling adequacy and sphericity --------------------------------------------
      
      v Overall KMO = <num>: marvellous (Kaiser & Rice, 1974).
      v Bartlett's test of sphericity: χ²(153) = <num>, p < <num> (significant).
      
      -- Multicollinearity -----------------------------------------------------------
      
      v Condition index: <num>. An index of 10 or less is rarely of interest
      (Belsley, 1991).
      
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
      
      v Overall KMO = <num>: marvellous (Kaiser & Rice, 1974).
      v Bartlett's test of sphericity: χ²(28) = <num>, p < <num> (significant).
      
      -- Multicollinearity -----------------------------------------------------------
      
      v Condition index: <num>. An index of 10 or less is rarely of interest
      (Belsley, 1991).
      
      -- Per-variable diagnostics ----------------------------------------------------
      
                variance   SMC   MSA  flags
      fun <num> <num> <num>       
      friends <num> <num> <num>       
      enjoy <num> <num> <num> sparse
      hurt <num> <num> <num> sparse
      part <num> <num> <num>       
      commonly <num> <num> <num>       
      chances <num> <num> <num>       
      attracted <num> <num> <num> sparse
      
      -- Multivariate normality ------------------------------------------------------
      
      x Mardia's skewness: χ²(120) = <num>, p < <num>.
      x Mardia's kurtosis: z = <num>, p < <num>.
      These data depart from multivariate normality: both tests reject it.
      
      -- Outliers --------------------------------------------------------------------
      
      i No case is flagged: some variables have fewer than ten distinct values
      (categorical), and a chi-square cutoff is not valid for their Mahalanobis
      distances. The largest distances are in rows 561, 238, 622, 549, and 186.
      
      -- Recommendations -------------------------------------------------------------
      
      ! These data depart from multivariate normality; normal-theory standard errors
        and fit statistics may be biased - prefer robust (sandwich) or bootstrapped
        standard errors.
      ! 3 variables have a sparse response category (< 5 responses): enjoy, hurt, and
        attracted; a low-frequency category can destabilise polychoric estimates -
        consider collapsing it into an adjacent category.

---

    Code
      print(scr_non)
    Output
      
      -- Sampling adequacy and sphericity --------------------------------------------
      
      v Overall KMO = <num>: marvellous (Kaiser & Rice, 1974).
      ! Bartlett's test of sphericity was not computed; no sample size (N) was
      supplied.
      
      -- Multicollinearity -----------------------------------------------------------
      
      v Condition index: <num>. An index of 10 or less is rarely of interest
      (Belsley, 1991).
      
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

# the report tracks a narrow console

    Code
      print(scr_cor)
    Output
      
      -- Sampling adequacy and sphericity ------------------------
      
      v Overall KMO = <num>: marvellous (Kaiser & Rice, 1974).
      v Bartlett's test of sphericity: χ²(153) = <num>, p < <num> (significant).
      
      -- Multicollinearity ---------------------------------------
      
      v Condition index: <num>. An index of 10 or less is rarely
      of interest (Belsley, 1991).
      
      -- Per-variable diagnostics --------------------------------
      
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
      
      -- Recommendations -----------------------------------------
      
      v The data appear suitable for factor analysis.
      i Per-item variance, missing-data, category, normality, and
        outlier diagnostics require raw data; only a correlation
        matrix was supplied.

