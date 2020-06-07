efa_cor <- EFA(test_models$baseline$cormat, n_factors = 3, N = 500)
efa_raw <- EFA(GRiPS_raw, n_factors = 1)

set.seed(500)
# HIER STEHENGEBLIEBEN --> evtl. rumspielen mit 4, 3, 2 indicators
efa_rand <- EFA(matrix(rnorm(100), ncol = 4))
