
# list of EFAs from similar datasets (10 bootstrap samples)
set.seed(42)
efa_list <- lapply(1:10,
                   function (x) {
                     EFA(DOSPERT_raw[sample(1:nrow(DOSPERT_raw), replace = TRUE),],
                         n_factors = 6, method = "ML", rotation = "promax")
                   })

# loadings lists
unrot_loadings <- lapply(efa_list, `[[`, "unrot_loadings")
rot_loadings <- lapply(efa_list, `[[`, "rot_loadings")


test_that("PROCRUSTES matches psych::Procrustes and GPArotation::targetQ outputs", {
  expect_equivalent(
    psych::Procrustes(unrot_loadings[[2]], rot_loadings[[1]])$loadings,
    PROCRUSTES(unrot_loadings[[2]], rot_loadings[[1]])$loadings
    )
  expect_equivalent(
    GPArotation::targetQ(unrot_loadings[[2]], Target = rot_loadings[[1]])$loadings,
    PROCRUSTES(unrot_loadings[[2]], rot_loadings[[1]], rotation = "oblique")$loadings,
    tolerance = 1e-3
  )
})

