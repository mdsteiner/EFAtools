# the mixed-se abort names the imputations that broke the shared se

    Code
      EFAtools:::.efa_pooled_mixed_se_abort(fits_mixed)
    Condition
      Error in `EFAtools:::.efa_pooled_mixed_se_abort()`:
      ! The component `efa_fit()` fits use different `se` methods, so their standard errors cannot be pooled.
      x Imputation 2 recorded `se` "sandwich"; the remaining 2 recorded "information".
      i Re-fit every imputation with the same `se` (all "none", "information", "sandwich", or "np-boot").

