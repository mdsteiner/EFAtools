# Test helper: return the record with the given `name` from an efa_retention
# object (or NULL if no such record is present).
.retention_record <- function(obj, name) {
  rec <- Filter(function(r) r$name == name, obj$results)
  if (length(rec) == 0) NULL else rec[[1]]
}

# Test helper: collect the class strings of every warning raised while evaluating
# `expr`, so a single expected class can be asserted on a path where several fire
# without expect_warning() stopping at the first one.
.warn_classes <- function(expr) {
  classes <- character()
  withCallingHandlers(
    expr,
    warning = function(w) {
      classes <<- c(classes, class(w))
      invokeRestart("muffleWarning")
    })
  classes
}

# Test helper: an 11-variable correlation matrix that is not positive definite
# (smallest eigenvalue -0.0245). Used by the retention tests that need an input
# which has to be smoothed and whose factor solutions are ill-behaved.
.burt_cormat <- function() {
  matrix(c(1.00,  0.83,  0.81,  0.80,   0.71, 0.70, 0.54, 0.53,  0.59,  0.24, 0.13,
           0.83,  1.00,  0.87,  0.62,   0.59, 0.44, 0.58, 0.44,  0.23,  0.45,  0.21,
           0.81,  0.87,  1.00,  0.63,   0.37, 0.31, 0.30, 0.12,  0.33,  0.33,  0.36,
           0.80,  0.62,  0.63,  1.00,   0.49, 0.54, 0.30, 0.28,  0.42,  0.29, -0.06,
           0.71,  0.59,  0.37,  0.49,   1.00, 0.54, 0.34, 0.55,  0.40,  0.19, -0.10,
           0.70,  0.44,  0.31,  0.54,   0.54, 1.00, 0.50, 0.51,  0.31,  0.11,  0.10,
           0.54,  0.58,  0.30,  0.30,   0.34, 0.50, 1.00, 0.38,  0.29,  0.21,  0.08,
           0.53,  0.44,  0.12,  0.28,   0.55, 0.51, 0.38, 1.00,  0.53,  0.10, -0.16,
           0.59,  0.23,  0.33,  0.42,   0.40, 0.31, 0.29, 0.53,  1.00, -0.09, -0.10,
           0.24,  0.45,  0.33,  0.29,   0.19, 0.11, 0.21, 0.10, -0.09,  1.00,  0.41,
           0.13,  0.21,  0.36, -0.06,  -0.10, 0.10, 0.08, -0.16, -0.10, 0.41,  1.00),
         nrow = 11, ncol = 11)
}
