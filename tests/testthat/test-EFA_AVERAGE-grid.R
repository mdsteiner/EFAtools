# Guards .build_avg_grid() against the previous explicit, per-(method, type)
# construction of the EFA_AVERAGE implementation grid. `ref_build_avg_grid()` is
# the earlier construction kept verbatim; the test asserts the current builder
# produces a byte-identical grid (values, column types, and row names) for every
# combination of estimators, types, and rotations, including vector-valued
# user arguments. Both implementations call the same .type_grid() helper, so any
# difference can only come from the arguments assembled per cell.

ref_build_avg_grid <- function(method, type, rotation, init_comm, criterion,
                               criterion_type, abs_eigen, start_method, k_promax,
                               normalize, P_type, precision, varimax_type,
                               k_simplimax) {

  grid_list <- list()

  if ("PAF" %in% method) {

    if ("EFAtools" %in% type) {
      grid_list[["ftls_pf"]] <- .type_grid(method = "PAF", init_comm = "smc",
                                           criterion = 1e-3, criterion_type = "sum",
                                           abs_eigen = TRUE, start_method = NA,
                                           rotation = rotation, k_promax = 4,
                                           normalize = TRUE, P_type = "norm",
                                           precision = 1e-5,
                                           varimax_type = "kaiser",
                                           k_simplimax = k_simplimax)
    }

    if ("psych" %in% type) {
      grid_list[["psch_pf"]] <- .type_grid(method = "PAF", init_comm = "smc",
                                           criterion = 1e-3, criterion_type = "sum",
                                           abs_eigen = FALSE, start_method = NA,
                                           rotation = rotation, k_promax = 4,
                                           normalize = TRUE, P_type = "unnorm",
                                           precision = 1e-5,
                                           varimax_type = "svd",
                                           k_simplimax = k_simplimax)
    }

    if ("SPSS" %in% type) {
      grid_list[["spss_pf"]] <- .type_grid(method = "PAF", init_comm = "smc",
                                           criterion = 1e-3, criterion_type = "max_individual",
                                           abs_eigen = TRUE, start_method = NA,
                                           rotation = rotation, k_promax = 4,
                                           normalize = TRUE, P_type = "norm",
                                           precision = 1e-5,
                                           varimax_type = "kaiser",
                                           k_simplimax = k_simplimax)
    }

    if ("none" %in% type) {

        grid_list[["nn_pf"]] <- .type_grid(method = "PAF", init_comm = init_comm,
                                           criterion = criterion, criterion_type = criterion_type,
                                           abs_eigen = abs_eigen, start_method = NA,
                                           rotation = rotation, k_promax = k_promax,
                                           normalize = normalize, P_type = P_type,
                                           precision = precision,
                                           varimax_type = varimax_type,
                                           k_simplimax = k_simplimax)
    }
  }

    if ("ML" %in% method) {

      if ("EFAtools" %in% type) {
        grid_list[["ftls_ml"]] <- .type_grid(method = "ML", init_comm = NA,
                                             criterion = NA, criterion_type = NA,
                                             abs_eigen = NA, start_method = "psych",
                                             rotation = rotation, k_promax = 4,
                                             normalize = TRUE, P_type = "norm",
                                             precision = 1e-5,
                                             varimax_type = "kaiser",
                                             k_simplimax = k_simplimax)
      }

      if ("psych" %in% type) {
        grid_list[["psch_ml"]] <- .type_grid(method = "ML", init_comm = NA,
                                             criterion = NA, criterion_type = NA,
                                             abs_eigen = NA, start_method = "psych",
                                             rotation = rotation, k_promax = 4,
                                             normalize = TRUE, P_type = "unnorm",
                                             precision = 1e-5,
                                             varimax_type = "svd",
                                             k_simplimax = k_simplimax)
      }

      if ("SPSS" %in% type) {
        grid_list[["spss_ml"]] <- .type_grid(method = "ML", init_comm = NA,
                                             criterion = NA, criterion_type = NA,
                                             abs_eigen = NA, start_method = "psych",
                                             rotation = rotation, k_promax = 4,
                                             normalize = TRUE, P_type = "norm",
                                             precision = 1e-5,
                                             varimax_type = "kaiser",
                                             k_simplimax = k_simplimax)
      }

      if ("none" %in% type) {
          grid_list[["nn_ml"]] <- .type_grid(method = "ML", init_comm = NA,
                                             criterion = NA, criterion_type = NA,
                                             abs_eigen = NA, start_method = start_method,
                                             rotation = rotation, k_promax = k_promax,
                                             normalize = normalize, P_type = P_type,
                                             precision = precision,
                                             varimax_type = varimax_type,
                                             k_simplimax = k_simplimax)

      }
    }

      if ("ULS" %in% method) {

        if ("EFAtools" %in% type) {
          grid_list[["ftls_ls"]] <- .type_grid(method = "ULS", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, start_method = NA,
                                               rotation = rotation, k_promax = 4,
                                               normalize = TRUE, P_type = "norm",
                                               precision = 1e-5,
                                               varimax_type = "kaiser",
                                               k_simplimax = k_simplimax)
        }

        if ("psych" %in% type) {
          grid_list[["psch_ls"]] <- .type_grid(method = "ULS", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, start_method = NA,
                                               rotation = rotation, k_promax = 4,
                                               normalize = TRUE, P_type = "unnorm",
                                               precision = 1e-5,
                                               varimax_type = "svd",
                                               k_simplimax = k_simplimax)
        }

        if ("SPSS" %in% type) {
          grid_list[["spss_ls"]] <- .type_grid(method = "ULS", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, start_method = NA,
                                               rotation = rotation, k_promax = 4,
                                               normalize = TRUE, P_type = "norm",
                                               precision = 1e-5,
                                               varimax_type = "kaiser",
                                               k_simplimax = k_simplimax)
        }

        if ("none" %in% type) {

            grid_list[["nn_ls"]] <- .type_grid(method = "ULS", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, start_method = NA,
                                               rotation = rotation, k_promax = k_promax,
                                               normalize = normalize, P_type = P_type,
                                               precision = precision,
                                               varimax_type = varimax_type,
                                               k_simplimax = k_simplimax)

        }
      }

  unique(do.call(rbind, grid_list))
}

test_that(".build_avg_grid reproduces the original grid construction", {

  method_sets <- list("PAF", "ML", "ULS", c("PAF", "ML", "ULS"))
  type_sets   <- list("EFAtools", "psych", "SPSS", "none",
                      c("EFAtools", "psych", "SPSS", "none"))
  # Each set is a single rotation kind (none / oblique / orthogonal), as required.
  rot_sets    <- list("none", "oblique", "orthogonal", "promax", "varimax",
                      "simplimax", c("oblimin", "geominQ"),
                      c("varimax", "quartimax"))

  base_args <- list(init_comm = "smc", criterion = 1e-3, criterion_type = "sum",
                    abs_eigen = TRUE, start_method = "psych", k_promax = 4,
                    normalize = TRUE, P_type = "norm", precision = 1e-5,
                    varimax_type = "kaiser", k_simplimax = 1e-6)

  for (mm in method_sets) {
    for (tt in type_sets) {
      for (rr in rot_sets) {
        args <- c(list(method = mm, type = tt, rotation = rr), base_args)
        desc <- paste(paste(mm, collapse = "+"), "|", paste(tt, collapse = "+"),
                      "|", paste(rr, collapse = "+"))
        expect_equal(do.call(.build_avg_grid, args),
                     do.call(ref_build_avg_grid, args), info = desc)
      }
    }
  }

  # Vector-valued user arguments are only forwarded by type "none"; they must
  # expand to the same multi-row grid through both implementations.
  vec_args <- list(method = c("PAF", "ML", "ULS"), type = "none",
                   rotation = c("oblimin", "geominQ"),
                   init_comm = c("smc", "mac"), criterion = c(1e-3, 1e-4),
                   criterion_type = "sum", abs_eigen = TRUE,
                   start_method = c("psych", "factanal"), k_promax = c(3, 4),
                   normalize = TRUE, P_type = c("norm", "unnorm"),
                   precision = c(1e-4, 1e-5), varimax_type = "kaiser",
                   k_simplimax = c(1e-5, 1e-6))
  expect_equal(do.call(.build_avg_grid, vec_args),
               do.call(ref_build_avg_grid, vec_args))
})
