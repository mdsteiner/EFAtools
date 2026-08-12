# the mixed-se abort names the imputations that broke the shared se

    Code
      EFAtools:::.efa_pooled_mixed_se_abort(fits_boot, cor_input)
    Condition
      Error in `EFAtools:::.efa_pooled_mixed_se_abort()`:
      ! The component `efa_fit()` fits use different `se` methods, so their standard errors cannot be pooled.
      x Imputation 2 recorded `se` "none"; the remaining 2 recorded "np-boot".
      i Imputation 2 is a correlation matrix, which cannot be bootstrapped, so `efa_fit()` fitted it with `se = "none"`.

---

    Code
      EFAtools:::.efa_pooled_mixed_se_abort(fits_mixed, raw)
    Condition
      Error in `EFAtools:::.efa_pooled_mixed_se_abort()`:
      ! The component `efa_fit()` fits use different `se` methods, so their standard errors cannot be pooled.
      x Imputation 2 recorded `se` "sandwich"; the remaining 2 recorded "information".
      i Re-fit every imputation with the same `se` (all "none", "information", "sandwich", or "np-boot").

