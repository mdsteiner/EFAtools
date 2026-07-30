# Guards .build_avg_grid() against the previous explicit, per-(estimator, type)
# construction of the efa_average implementation grid. `ref_build_avg_grid()` is
# the earlier construction kept verbatim; the test asserts the current builder
# produces a byte-identical grid (values, column types, and row names) for every
# combination of estimators, types, and rotations, including vector-valued
# user arguments. Both implementations call the same .type_grid() helper, so any
# difference can only come from the arguments assembled per cell.

ref_build_avg_grid <- function(estimator, type, rotation, init_comm, criterion,
                               criterion_type, abs_eigen, max_iter, start_method,
                               k_promax, normalize, P_type, precision, varimax_type,
                               k_simplimax) {

  grid_list <- list()

  if ("PAF" %in% estimator) {

    if ("EFAtools" %in% type) {
      grid_list[["ftls_pf"]] <- .type_grid(estimator = "PAF", init_comm = "smc",
                                           criterion = 1e-3, criterion_type = "sum",
                                           abs_eigen = TRUE, max_iter = 300,
                                           start_method = NA,
                                           rotation = rotation, k_promax = 4,
                                           normalize = TRUE, P_type = "norm",
                                           precision = 1e-5,
                                           varimax_type = "kaiser",
                                           k_simplimax = k_simplimax)
    }

    if ("psych" %in% type) {
      grid_list[["psch_pf"]] <- .type_grid(estimator = "PAF", init_comm = "smc",
                                           criterion = 1e-3, criterion_type = "sum",
                                           abs_eigen = FALSE, max_iter = 50,
                                           start_method = NA,
                                           rotation = rotation, k_promax = 4,
                                           normalize = TRUE, P_type = "unnorm",
                                           precision = 1e-5,
                                           varimax_type = "svd",
                                           k_simplimax = k_simplimax)
    }

    if ("SPSS" %in% type) {
      grid_list[["spss_pf"]] <- .type_grid(estimator = "PAF", init_comm = "smc",
                                           criterion = 1e-3, criterion_type = "max_individual",
                                           abs_eigen = TRUE, max_iter = 25,
                                           start_method = NA,
                                           rotation = rotation, k_promax = 4,
                                           normalize = TRUE, P_type = "norm",
                                           precision = 1e-5,
                                           varimax_type = "kaiser",
                                           k_simplimax = k_simplimax)
    }

    if ("none" %in% type) {

        grid_list[["nn_pf"]] <- .type_grid(estimator = "PAF", init_comm = init_comm,
                                           criterion = criterion, criterion_type = criterion_type,
                                           abs_eigen = abs_eigen, max_iter = max_iter,
                                           start_method = NA,
                                           rotation = rotation, k_promax = k_promax,
                                           normalize = normalize, P_type = P_type,
                                           precision = precision,
                                           varimax_type = varimax_type,
                                           k_simplimax = k_simplimax)
    }
  }

    if ("ML" %in% estimator) {

      if ("EFAtools" %in% type) {
        grid_list[["ftls_ml"]] <- .type_grid(estimator = "ML", init_comm = NA,
                                             criterion = NA, criterion_type = NA,
                                             abs_eigen = NA, max_iter = NA,
                                             start_method = "psych",
                                             rotation = rotation, k_promax = 4,
                                             normalize = TRUE, P_type = "norm",
                                             precision = 1e-5,
                                             varimax_type = "kaiser",
                                             k_simplimax = k_simplimax)
      }

      if ("psych" %in% type) {
        grid_list[["psch_ml"]] <- .type_grid(estimator = "ML", init_comm = NA,
                                             criterion = NA, criterion_type = NA,
                                             abs_eigen = NA, max_iter = NA,
                                             start_method = "psych",
                                             rotation = rotation, k_promax = 4,
                                             normalize = TRUE, P_type = "unnorm",
                                             precision = 1e-5,
                                             varimax_type = "svd",
                                             k_simplimax = k_simplimax)
      }

      if ("SPSS" %in% type) {
        grid_list[["spss_ml"]] <- .type_grid(estimator = "ML", init_comm = NA,
                                             criterion = NA, criterion_type = NA,
                                             abs_eigen = NA, max_iter = NA,
                                             start_method = "psych",
                                             rotation = rotation, k_promax = 4,
                                             normalize = TRUE, P_type = "norm",
                                             precision = 1e-5,
                                             varimax_type = "kaiser",
                                             k_simplimax = k_simplimax)
      }

      if ("none" %in% type) {
          grid_list[["nn_ml"]] <- .type_grid(estimator = "ML", init_comm = NA,
                                             criterion = NA, criterion_type = NA,
                                             abs_eigen = NA, max_iter = NA,
                                             start_method = start_method,
                                             rotation = rotation, k_promax = k_promax,
                                             normalize = normalize, P_type = P_type,
                                             precision = precision,
                                             varimax_type = varimax_type,
                                             k_simplimax = k_simplimax)

      }
    }

      if ("ULS" %in% estimator) {

        if ("EFAtools" %in% type) {
          grid_list[["ftls_ls"]] <- .type_grid(estimator = "ULS", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, max_iter = NA,
                                               start_method = NA,
                                               rotation = rotation, k_promax = 4,
                                               normalize = TRUE, P_type = "norm",
                                               precision = 1e-5,
                                               varimax_type = "kaiser",
                                               k_simplimax = k_simplimax)
        }

        if ("psych" %in% type) {
          grid_list[["psch_ls"]] <- .type_grid(estimator = "ULS", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, max_iter = NA,
                                               start_method = NA,
                                               rotation = rotation, k_promax = 4,
                                               normalize = TRUE, P_type = "unnorm",
                                               precision = 1e-5,
                                               varimax_type = "svd",
                                               k_simplimax = k_simplimax)
        }

        if ("SPSS" %in% type) {
          grid_list[["spss_ls"]] <- .type_grid(estimator = "ULS", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, max_iter = NA,
                                               start_method = NA,
                                               rotation = rotation, k_promax = 4,
                                               normalize = TRUE, P_type = "norm",
                                               precision = 1e-5,
                                               varimax_type = "kaiser",
                                               k_simplimax = k_simplimax)
        }

        if ("none" %in% type) {

            grid_list[["nn_ls"]] <- .type_grid(estimator = "ULS", init_comm = NA,
                                               criterion = NA, criterion_type = NA,
                                               abs_eigen = NA, max_iter = NA,
                                               start_method = NA,
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

  estimator_sets <- list("PAF", "ML", "ULS", c("PAF", "ML", "ULS"))
  type_sets   <- list("EFAtools", "psych", "SPSS", "none",
                      c("EFAtools", "psych", "SPSS", "none"))
  # Each set is a single rotation kind (none / oblique / orthogonal), as required.
  rot_sets    <- list("none", "oblique", "orthogonal", "promax", "varimax",
                      "simplimax", c("oblimin", "geominQ"),
                      c("varimax", "quartimax"))

  # A non-preset max_iter for the type "none" rows, so the assembled grid carries a
  # genuinely varying max_iter column on the PAF rows (25/50/300 for the named types,
  # 500 here); it stays NA on the ML and ULS rows, which never use the cap.
  base_args <- list(init_comm = "smc", criterion = 1e-3, criterion_type = "sum",
                    abs_eigen = TRUE, max_iter = 500, start_method = "psych",
                    k_promax = 4, normalize = TRUE, P_type = "norm", precision = 1e-5,
                    varimax_type = "kaiser", k_simplimax = 1e-6)

  for (mm in estimator_sets) {
    for (tt in type_sets) {
      for (rr in rot_sets) {
        args <- c(list(estimator = mm, type = tt, rotation = rr), base_args)
        desc <- paste(paste(mm, collapse = "+"), "|", paste(tt, collapse = "+"),
                      "|", paste(rr, collapse = "+"))
        expect_equal(do.call(.build_avg_grid, args),
                     do.call(ref_build_avg_grid, args), info = desc)
      }
    }
  }

  # Vector-valued user arguments are only forwarded by type "none"; they must
  # expand to the same multi-row grid through both implementations.
  vec_args <- list(estimator = c("PAF", "ML", "ULS"), type = "none",
                   rotation = c("oblimin", "geominQ"),
                   init_comm = c("smc", "mac"), criterion = c(1e-3, 1e-4),
                   criterion_type = "sum", abs_eigen = TRUE, max_iter = 500,
                   start_method = c("psych", "factanal"), k_promax = c(3, 4),
                   normalize = TRUE, P_type = c("norm", "unnorm"),
                   precision = c(1e-4, 1e-5), varimax_type = "kaiser",
                   k_simplimax = c(1e-5, 1e-6))
  expect_equal(do.call(.build_avg_grid, vec_args),
               do.call(ref_build_avg_grid, vec_args))
})


### test .oblq_grid

obl_grid_1 <- .oblq_grid(c("PAF"), c("smc", "mac"), .001,
                         c("sum", "max_individual"), c(FALSE, TRUE),
                         500, NA, c("promax", "simplimax", "oblimin"),
                         c(3, 4), TRUE, c("norm", "unnorm"), 1e-5, c("kaiser", "svd"),
                         30)
obl_grid_2 <- .oblq_grid(c("PAF"), c("smc", "mac"), .001,
                         c("sum", "max_individual"), c(FALSE, TRUE),
                         500, NA, c("simplimax", "oblimin"),
                         c(3, 4), TRUE, c("norm", "unnorm"), 1e-5, c("kaiser", "svd"),
                         30)
obl_grid_3 <- .oblq_grid(c("PAF"), c("smc", "mac"), .001,
                         c("sum", "max_individual"), c(FALSE, TRUE),
                         500, NA, c("simplimax", "oblimin"),
                         NA, TRUE, NA, 1e-5, NA, 30)
obl_grid_4 <- .oblq_grid("ML", NA, NA, NA, NA, 500, c("psych", "factanal"),
                         "oblimin", NA, TRUE, NA, 1e-5, NA, NA)

test_that(".oblq_grid works", {
  ### tests for obl_grid_1 with vector-valued PAF arguments
  expect_s3_class(obl_grid_1, "data.frame")
  expect_named(obl_grid_1, c("estimator", "init_comm", "criterion", "criterion_type",
                            "abs_eigen", "max_iter", "start_method", "rotation", "k_promax",
                            "normalize", "P_type", "precision", "varimax_type",
                            "k_simplimax"))
  expect_equal(nrow(obl_grid_1), 80)
  expect_equal(sum(is.na(obl_grid_1$k_simplimax)), 72)
  expect_equal(sum(is.na(obl_grid_1$k_promax)), 16)
  expect_equal(sum(is.na(obl_grid_1$varimax_type)), 16)
  expect_equal(unique(obl_grid_1$rotation), c("promax", "simplimax", "oblimin"))

  expect_s3_class(obl_grid_2, "data.frame")
  expect_named(obl_grid_2, c("estimator", "init_comm", "criterion", "criterion_type",
                             "abs_eigen", "max_iter", "start_method", "rotation", "k_promax",
                             "normalize", "P_type", "precision", "varimax_type",
                             "k_simplimax"))
  expect_equal(nrow(obl_grid_2), 16)
  expect_equal(sum(is.na(obl_grid_2$k_promax)), 16)
  expect_equal(sum(is.na(obl_grid_2$varimax_type)), 16)
  expect_equal(unique(obl_grid_2$rotation), c("simplimax", "oblimin"))

  expect_equal(obl_grid_2, obl_grid_3)

  expect_s3_class(obl_grid_4, "data.frame")
  expect_named(obl_grid_4, c("estimator", "init_comm", "criterion", "criterion_type",
                             "abs_eigen", "max_iter", "start_method", "rotation", "k_promax",
                             "normalize", "P_type", "precision", "varimax_type",
                             "k_simplimax"))
  expect_equal(nrow(obl_grid_4), 2)
  expect_equal(sum(is.na(obl_grid_4$k_simplimax)), 2)
  expect_equal(sum(is.na(obl_grid_4$k_promax)), 2)
  expect_equal(sum(is.na(obl_grid_4$varimax_type)), 2)
  expect_equal(unique(obl_grid_4$rotation), c("oblimin"))
  expect_equal(sum(is.na(obl_grid_4$init_comm)), 2)


})
### test .orth_grid

orth_grid_1 <- .orth_grid(c("PAF"), c("smc", "mac"), .001,
                         c("sum", "max_individual"), c(FALSE, TRUE),
                         500, NA, c("varimax", "quartimax"),
                         TRUE, 1e-5, c("kaiser", "svd"))
orth_grid_2 <- .orth_grid(c("PAF"), c("smc", "mac"), .001,
                          c("sum", "max_individual"), c(FALSE, TRUE),
                          500, NA, c("quartimax"), TRUE, 1e-5,
                          c("kaiser", "svd"))
orth_grid_3 <- .orth_grid(c("PAF"), c("smc", "mac"), .001,
                          c("sum", "max_individual"), c(FALSE, TRUE),
                          500, NA, c("quartimax"), TRUE, 1e-5, NA)
orth_grid_4 <- .orth_grid("ML", NA, NA, NA, NA, 500, c("psych", "factanal"),
                         "quartimax", TRUE, 1e-5, NA)

test_that(".orth_grid works", {
  ### tests for orth_grid_1 with vector-valued PAF arguments
  expect_s3_class(orth_grid_1, "data.frame")
  expect_named(orth_grid_1, c("estimator", "init_comm", "criterion", "criterion_type",
                             "abs_eigen", "max_iter", "start_method", "rotation", "k_promax",
                             "normalize", "P_type", "precision", "varimax_type",
                             "k_simplimax"))
  expect_equal(nrow(orth_grid_1), 24)
  expect_equal(sum(is.na(orth_grid_1$varimax_type)), 8)
  expect_equal(sum(is.na(orth_grid_1$k_promax)), 24)
  expect_equal(unique(orth_grid_1$rotation), c("varimax", "quartimax"))

  expect_s3_class(orth_grid_2, "data.frame")
  expect_named(orth_grid_2, c("estimator", "init_comm", "criterion", "criterion_type",
                             "abs_eigen", "max_iter", "start_method", "rotation", "k_promax",
                             "normalize", "P_type", "precision", "varimax_type",
                             "k_simplimax"))
  expect_equal(nrow(orth_grid_2), 8)
  expect_equal(sum(is.na(orth_grid_2$varimax_type)), 8)
  expect_equal(sum(is.na(orth_grid_2$k_promax)), 8)
  expect_equal(unique(orth_grid_2$rotation), c("quartimax"))

  expect_equal(orth_grid_2, orth_grid_3)

  expect_s3_class(orth_grid_4, "data.frame")
  expect_named(orth_grid_4, c("estimator", "init_comm", "criterion", "criterion_type",
                             "abs_eigen", "max_iter", "start_method", "rotation", "k_promax",
                             "normalize", "P_type", "precision", "varimax_type",
                             "k_simplimax"))
  expect_equal(nrow(orth_grid_4), 2)
  expect_equal(sum(is.na(orth_grid_4$k_simplimax)), 2)
  expect_equal(sum(is.na(orth_grid_4$k_promax)), 2)
  expect_equal(sum(is.na(orth_grid_4$varimax_type)), 2)
  expect_equal(unique(orth_grid_4$rotation), c("quartimax"))
  expect_equal(sum(is.na(orth_grid_4$init_comm)), 2)

})


### test .type_grid

tg_ob <- .type_grid("PAF", c("smc", "mac"), .001,
                   c("sum", "max_individual"), c(FALSE, TRUE),
                   500, NA, "oblique", c(3, 4), TRUE, c("norm", "unnorm"),
                   1e-5, c("kaiser", "svd"), 30)
tg_ob2 <- .type_grid("PAF", c("smc", "mac"), .001,
                    c("sum", "max_individual"), c(FALSE, TRUE),
                    500, NA, c("oblimin", "promax"), c(3, 4), TRUE, c("norm", "unnorm"),
                    1e-5, c("kaiser", "svd"), 30)
tg_orth <- .type_grid("PAF", c("smc", "mac"), .001,
                    c("sum", "max_individual"), c(FALSE, TRUE),
                    500, NA, "orthogonal", c(3, 4), TRUE, c("norm", "unnorm"),
                    1e-5, c("kaiser", "svd"), 30)
tg_orth2 <- .type_grid("PAF", c("smc", "mac"), .001,
                      c("sum", "max_individual"), c(FALSE, TRUE),
                      500, NA, c("varimax", "quartimax"), c(3, 4), TRUE,
                      c("norm", "unnorm"), 1e-5, c("kaiser", "svd"), 30)
tg_nn <- .type_grid("PAF", c("smc", "mac"), .001,
                      c("sum", "max_individual"), c(FALSE, TRUE),
                      500, NA, "none", c(3, 4), TRUE, c("norm", "unnorm"),
                      1e-5, c("kaiser", "svd"), 30)

test_that(".type_grid works", {
  ### test errors
  expect_error(.type_grid("PAF", NA, NA, NA, NA, NA, NA, c("oblique", "none"), NA,
                          NA, NA, NA, NA, NA),
               class = "efa_rotation_length")
  expect_error(.type_grid("PAF", NA, NA, NA, NA, NA, NA, c("oblique", "varimax"), NA,
                          NA, NA, NA, NA, NA),
               class = "efa_rotation_length")
  expect_error(.type_grid("PAF", NA, NA, NA, NA, NA, NA, c("orthogonal", "varimax"),
                          NA, NA, NA, NA, NA, NA),
               class = "efa_rotation_length")
  expect_error(.type_grid("PAF", NA, NA, NA, NA, NA, NA, c("promax", "varimax"),
                          NA, NA, NA, NA, NA, NA),
               class = "efa_rotation_mismatch")

  expect_s3_class(tg_ob, "data.frame")
  expect_named(tg_ob, c("estimator", "init_comm", "criterion", "criterion_type",
                              "abs_eigen", "max_iter", "start_method", "rotation", "k_promax",
                              "normalize", "P_type", "precision", "varimax_type",
                              "k_simplimax"))
  expect_equal(nrow(tg_ob), 112)
  expect_equal(sum(is.na(tg_ob$varimax_type)),
               nrow(tg_ob) - sum(tg_ob$rotation == "promax"))
  expect_equal(sum(is.na(tg_ob$k_promax)),
               nrow(tg_ob) - sum(tg_ob$rotation == "promax"))
  expect_equal(sum(is.na(tg_ob$P_type)),
               nrow(tg_ob) - sum(tg_ob$rotation == "promax"))
  expect_equal(sum(is.na(tg_ob$k_simplimax)),
               nrow(tg_ob) - sum(tg_ob$rotation == "simplimax"))
  expect_equal(sort(unique(tg_ob$rotation)),
               sort(c("promax", "oblimin", "quartimin", "simplimax",
                 "bentlerQ", "geominQ", "bifactorQ")))

  expect_s3_class(tg_ob2, "data.frame")
  expect_named(tg_ob2, c("estimator", "init_comm", "criterion", "criterion_type",
                        "abs_eigen", "max_iter", "start_method", "rotation", "k_promax",
                        "normalize", "P_type", "precision", "varimax_type",
                        "k_simplimax"))
  expect_equal(nrow(tg_ob2), 72)
  expect_equal(sum(is.na(tg_ob2$varimax_type)),
               nrow(tg_ob2) - sum(tg_ob2$rotation == "promax"))
  expect_equal(sum(is.na(tg_ob2$k_promax)),
               nrow(tg_ob2) - sum(tg_ob2$rotation == "promax"))
  expect_equal(sum(is.na(tg_ob2$P_type)),
               nrow(tg_ob2) - sum(tg_ob2$rotation == "promax"))
  expect_equal(sum(is.na(tg_ob2$k_simplimax)),
               nrow(tg_ob2) - sum(tg_ob2$rotation == "simplimax"))
  expect_equal(unique(tg_ob2$rotation),
               c("promax", "oblimin"))

  expect_s3_class(tg_orth, "data.frame")
  expect_named(tg_orth, c("estimator", "init_comm", "criterion", "criterion_type",
                        "abs_eigen", "max_iter", "start_method", "rotation", "k_promax",
                        "normalize", "P_type", "precision", "varimax_type",
                        "k_simplimax"))
  expect_equal(nrow(tg_orth), 56)
  expect_equal(sum(is.na(tg_orth$varimax_type)),
               nrow(tg_orth) - sum(tg_orth$rotation == "varimax"))
  expect_equal(sort(unique(tg_orth$rotation)),
               sort(c("varimax", "quartimax", "equamax",
                      "bentlerT", "geominT", "bifactorT")))

  expect_s3_class(tg_orth2, "data.frame")
  expect_named(tg_orth2, c("estimator", "init_comm", "criterion", "criterion_type",
                         "abs_eigen", "max_iter", "start_method", "rotation", "k_promax",
                         "normalize", "P_type", "precision", "varimax_type",
                         "k_simplimax"))
  expect_equal(nrow(tg_orth2), 24)
  expect_equal(sum(is.na(tg_orth2$varimax_type)),
               nrow(tg_orth2) - sum(tg_orth2$rotation == "varimax"))
  expect_equal(unique(tg_orth2$rotation),
               c("varimax", "quartimax"))


  expect_s3_class(tg_nn, "data.frame")
  expect_named(tg_nn, c("estimator", "init_comm", "criterion", "criterion_type",
                           "abs_eigen", "max_iter", "start_method", "rotation", "k_promax",
                           "normalize", "P_type", "precision", "varimax_type",
                           "k_simplimax"))
  expect_equal(nrow(tg_nn), 8)
  expect_true(all(is.na(tg_nn$varimax_type)))
  expect_true(all(tg_nn$rotation == "none"))
  expect_true(all(is.na(tg_nn$k_promax)))
  expect_true(all(is.na(tg_nn$normalize)))
  expect_true(all(is.na(tg_nn$P_type)))
  expect_true(all(is.na(tg_nn$precision)))
  expect_true(all(is.na(tg_nn$k_simplimax)))

})


rm(obl_grid_1, obl_grid_2, obl_grid_3, obl_grid_4, orth_grid_1, orth_grid_2,
   orth_grid_3, orth_grid_4, tg_ob, tg_ob2, tg_orth, tg_orth2, tg_nn)
