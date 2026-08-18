test_that(".new_efa_retention rejects unknown criterion ids", {
  rec <- list(list(name = "X", label = "X", n_factors = 1, plot_type = "none"))
  expect_error(
    .new_efa_retention("NOPE", results = rec, settings = list()),
    class = "efa_unknown_criterion"
  )
})

test_that(".new_efa_retention builds the documented shape for a known id", {
  rec <- list(
    list(name = "TR2", label = "Original", n_factors = 3, plot_type = "eigen"),
    list(name = "TR4", label = "Revised", n_factors = 2, plot_type = "eigen")
  )
  out <- .new_efa_retention("MAP", results = rec, settings = list(a = 1))

  expect_s3_class(out, "efa_retention")
  expect_equal(unname(out$criterion["label"]), "Minimum average partial")
  expect_named(out$n_factors, c("TR2", "TR4"))
  expect_equal(unname(out$n_factors), c(3, 2))
})

test_that(".n_factors_ctl() defaults mirror efa_retain()'s", {
  # efa_power() calls .n_factors_ctl() relying on its defaults, so a default drifting from
  # efa_retain()'s would silently retain factors under undocumented settings.
  ctl <- .n_factors_ctl()
  shared <- intersect(names(ctl), names(formals(efa_retain)))
  ret <- lapply(formals(efa_retain)[shared], eval)
  # efa_retain() resolves these with a single-choice match.arg(), so its effective default
  # is the first element; every other default (including the several.ok ones, whose whole
  # vector is the default) has to match .n_factors_ctl()'s in full
  single <- c("use", "cor_method", "estimator", "eigen_type_HULL", "decision_rule")
  ret[single] <- lapply(ret[single], `[`, 1L)
  expect_equal(ctl[shared], ret[shared])
})

test_that(".eigen_subtitle names the eigenvalue types the same way for every criterion", {
  # The criteria that report which eigenvalues they used share this one sentence, so
  # pin both forms here: efa_parallel()'s `detail` variant is otherwise only covered
  # by a print snapshot behind the slow gate.
  expect_equal(.eigen_subtitle("SMC"), "Eigenvalues found using SMC.")
  expect_equal(.eigen_subtitle(c("PCA", "SMC")),
               "Eigenvalues found using PCA and SMC.")
  expect_equal(.eigen_subtitle(c("PCA", "SMC", "EFA")),
               "Eigenvalues found using PCA, SMC, and EFA.")

  # a detail clause replaces the closing full stop rather than following it
  expect_equal(.eigen_subtitle("SMC", "1000 simulated datasets"),
               "Eigenvalues found using SMC; 1000 simulated datasets.")

  # and the criteria really route through it
  expect_equal(efa_scree(test_models$baseline$cormat, eigen_type = "PCA")$subtitle,
               .eigen_subtitle("PCA"))
  expect_equal(efa_kgc(test_models$baseline$cormat, eigen_type = "PCA")$subtitle,
               .eigen_subtitle("PCA"))
})

test_that("printed bullets name what the count belongs to, not the eigenvalues", {
  # Every bullet reads "<variant>: <suggested count>": the eigenvalue-based criteria
  # key theirs by eigenvalue type, the Hull method by goodness-of-fit index, and the
  # single-variant criteria name the quantity itself. The eigenvalues the counts were
  # derived from belong in the subtitle, so that a count is never labelled as if it
  # were the eigenvalues.
  kgc <- efa_kgc(test_models$baseline$cormat)
  expect_equal(vapply(kgc$results, function(r) r$label, character(1)),
               c("PCA", "SMC", "EFA"))
  expect_equal(kgc$subtitle, "Eigenvalues found using PCA, SMC, and EFA.")

  scree <- efa_scree(test_models$baseline$cormat, eigen_type = "SMC")
  expect_equal(vapply(scree$results, function(r) r$label, character(1)), "SMC")
  expect_equal(scree$subtitle, "Eigenvalues found using SMC.")
  # the one purely visual criterion names the call that shows the plot
  expect_match(scree$note, "plot(x)", fixed = TRUE)
})

test_that(".retention_count() renders a count independently of the number options", {
  # Counts are whole by construction, so they are rendered the same way whatever
  # options(scipen) asks for. The two properties the call sites rely on:
  skip_if_not_installed("withr")
  withr::local_options(scipen = -5, digits = 2, OutDec = ",")

  # fixed notation, for a double as much as for an integer
  expect_identical(.retention_count(3), "3")
  expect_identical(.retention_count(3L), "3")
  expect_identical(.retention_count(1000), "1000")

  # and no padding to a common width: the hull draws a vector of point labels, and
  # .retention_summary() lists tied modal values, so a padded element would be visible
  expect_identical(.retention_count(c(0, 1, 2, 10)), c("0", "1", "2", "10"))
})

test_that("counts are rendered the same way whatever options(scipen) is", {
  # The criteria differ in the storage mode they put in their record: a count derived with
  # cumprod() or which.min() comes out double, one taken from an index integer. Under a
  # negative options(scipen) a double is rendered in scientific notation ("3e+00" instead
  # of "3") by as.character() in the report, and also by grid, which coerces a numeric plot
  # label at draw time. Both surfaces are therefore checked for every criterion, because
  # each renders its count on its own path.
  #
  # The option is set BEFORE the fixtures are built, because a criterion can render a
  # count into a string when it constructs its object rather than when the object is
  # printed (efa_parallel() does this for its subtitle). Fixtures built first would only
  # exercise the format()-time paths.
  skip_if_not_installed("withr")
  withr::local_options(scipen = -5)
  cmat <- test_models$baseline$cormat
  set.seed(42)
  # PARALLEL keeps a two-digit `n_datasets` on purpose: its subtitle names that number, so
  # this is the fixture whose exponent is not "e+00", and it holds the pattern below to a
  # general one. The other simulation counts are as low as still yields a record.
  fits <- list(
    EKC = efa_ekc(cmat, N = 500),
    KGC = efa_kgc(cmat, eigen_type = c("PCA", "SMC")),
    MAP = efa_map(cmat),
    PARALLEL = efa_parallel(cmat, N = 500, n_datasets = 20, eigen_type = "SMC"),
    SCREE = efa_scree(cmat, eigen_type = "SMC"),
    SMT = efa_smt(cmat, N = 500),
    HULL = efa_hull(cmat, N = 500, gof = "CAF"),
    NEST = efa_nest(cmat, N = 500, n_datasets = 5),
    # CD is the remaining criterion that stores its count as a double; it needs raw data
    CD = suppressMessages(efa_cd(GRiPS_raw[1:300, ], N_pop = 1000, N_samples = 20))
  )

  # every text a plot draws, coerced the way grid coerces it
  plot_labels <- function(p) {
    unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) {
      if ("label" %in% names(d)) as.character(d$label) else NULL
    }))
  }

  # any exponent, not just "e+00": a count of ten or more carries a different one
  sci <- "e[+-][0-9]"
  for (id in names(fits)) {
    expect_false(any(grepl(sci, format(fits[[id]]))), info = id)
  }

  # the eigenvalue plots label the retained solution, the hull plot every point; the
  # criteria that draw neither (MAP, SMT, and SCREE, which retains nothing) have no label.
  # CD is checked on its report only: it is simulation-based, and a draw that suggested no
  # factors would leave its record with no point to label at all.
  for (id in c("EKC", "KGC", "PARALLEL", "NEST", "HULL")) {
    labels <- plot_labels(plot(fits[[id]]))
    expect_gt(length(labels), 0)
    expect_false(any(grepl(sci, labels)), info = id)
  }
})

test_that("the retention report fits a narrow console", {
  local_reproducible_output(width = 40)
  objs <- list(efa_kgc(test_models$baseline$cormat),
               efa_scree(test_models$baseline$cormat),
               efa_ekc(test_models$baseline$cormat, N = 500),
               efa_map(test_models$baseline$cormat))
  for (obj in objs) {
    expect_true(all(nchar(format(obj)) <= 40),
                info = obj$criterion[["id"]])
  }
})

test_that("format.efa_retention is the source of truth and honours the colour state", {
  # Two fixtures, because no one criterion covers both halves. The note is the only
  # styled element a report carries, and the scree plot is the criterion that always has
  # one -- but it makes no numeric suggestion, so its report has no bullets. EKC is the
  # other way round. Together they reach the rule, the note, and the bullet list.
  scree <- efa_scree(test_models$baseline$cormat, eigen_type = "SMC")
  ekc <- efa_ekc(test_models$baseline$cormat, N = 500)
  expect_true(all(is.na(scree$n_factors)))   # the bullets branch is skipped ...
  expect_true(any(!is.na(ekc$n_factors)))    # ... and taken

  # print() is exactly cat(format(x), sep = "\n"), so the two agree line for line.
  expect_identical(utils::capture.output(print(scree)), format(scree))
  expect_identical(utils::capture.output(print(ekc)), format(ekc))

  old <- options(cli.num_colors = 256)
  on.exit(options(old), add = TRUE)

  # With colours on the report embeds ANSI ...
  expect_true(any(grepl("\033", format(scree), fixed = TRUE)))

  # ... and with colours off it is plain.
  options(cli.num_colors = 1)
  expect_false(any(grepl("\033", format(scree), fixed = TRUE)))
  expect_false(any(grepl("\033", format(ekc), fixed = TRUE)))
})
