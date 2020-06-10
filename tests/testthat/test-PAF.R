# MAYBE USE POPULATION MODEL INSTEAD OF TEST MODEL HERE (SUCH THAT WE KNOW
# TRUE LOADINGS)
paf_efatools <- .PAF(test_models$baseline$cormat, n_factors, N = 500,
                     max_iter = NULL, type = "EFAtools", init_comm = NULL,
                     criterion = NULL, criterion_type = NULL, abs_eigen = NULL)
paf_psych <- .PAF(test_models$baseline$cormat, n_factors, N = 500,
                     max_iter = NULL, type = "psych", init_comm = NULL,
                     criterion = NULL, criterion_type = NULL, abs_eigen = NULL)
paf_spss <- .PAF(test_models$baseline$cormat, n_factors, N = 500,
                  max_iter = NULL, type = "SPSS", init_comm = NULL,
                  criterion = NULL, criterion_type = NULL, abs_eigen = NULL)
paf_none <- .PAF(test_models$baseline$cormat, n_factors, N = 500,
                 max_iter = 500, type = "none", init_comm = "unity",
                 criterion = 1e-4, criterion_type = "sums", abs_eigen = TRUE)
