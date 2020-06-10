# SEE IF ALL THIS IS NEEDED GIVEN THAT PAF ETC FUNCTIONS ARE TESTED SEPARATELY
efa_cor <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500)
efa_raw <- EFA(GRiPS_raw, n_factors = 1)

# different types
efa_psych <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                 type = "psych", rotation = "promax")
efa_spss <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                type = "SPSS", rotation = "promax")

# different methods
efa_ml <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
              method = "ML")
efa_uls <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
               method = "ULS")

# different rotation methods from GPA rotation package (orthogonal and oblique)
efa_equa <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                rotation = "equamax")
efa_quart <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                 rotation = "quartimin")

# PAF with promax rotation without a specified type
efa_none <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500,
                type = "none", method = "PAF", rotation = "promax",
                max_iter = 500, init_comm = "mac", criterion = 1e-4,
                criterion_type = "sums", abs_eigen = FALSE, k = 3,
                P_type = "unnorm", precision= 1e-5, order_type = "eigen")
    #'
# with random data
set.seed(500)
efa_rand <- EFA(matrix(rnorm(100), ncol = 4), n_factors = 1, max_iter = 1000)
