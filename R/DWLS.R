## Diagonally-Weighted Least Squares Estimation of Factor Loadings
## Thin fitter: returns the raw DWLS results; shared post-processing happens in
## .finalize_fit() / .estimate_model(). Unlike ULS, DWLS optimises the loadings
## directly (weighted off-diagonal least squares) so the per-element weights shape the
## solution; see .fit_dwls() and the DwlsFunctor in src/estimate.cpp.
.DWLS <- function(x, n_factors, weights) {

  if (is.null(weights)) {
    cli::cli_abort(
      c("DWLS estimation requires a per-element weight matrix.",
        "i" = "Weights are produced by {.fn .prepare_cor_input} with {.code acov = \"diag\"} on raw ordinal data."),
      class = "efa_dwls_no_weights"
    )
  }

  # Get correlation matrix entered or created in EFA
  R <- x

  dwls <- .fit_dwls(R, n_factors, weights)

  L <- dwls$loadings
  orig_R <- R
  h2 <- rowSums(L^2)            # diag(L L'), without forming the full p x p product
  diag(R) <- h2

  # raw fit, finalized by .estimate_model(). `Fm` is the weighted off-diagonal objective the
  # C++ backend already minimised. DWLS optimises the loadings with no lower bound on the
  # uniquenesses, so the boundary-uniqueness Heywood heuristic in .finalize_fit() (which only
  # applies to the box-constrained ML/ULS optimisers) does not apply here; psi is left NULL
  # so a Heywood case is flagged only by a communality at or above 1, which the unconstrained
  # optimiser can still produce.
  list(
    L = L,
    h2 = h2,
    psi = NULL,
    Fm = dwls$Fm,
    iter = dwls$iter,
    convergence = dwls$convergence,
    orig_R = orig_R,
    R_final = R
  )
}

# Obtain the DWLS fit. The C++ backend warm-starts with the weighted-ULS solution (the
# eigen extraction that minimises the weighted objective over the uniquenesses) and then
# polishes the free loadings to the diagonally weighted least squares optimum.
.fit_dwls <- function(R, n_fac, weights) {

  dwls <- .fit_dwls_cpp(R, n_fac, weights)

  list(loadings = dwls$loadings, Fm = dwls$Fm, iter = dwls$iter,
       convergence = dwls$convergence)
}
