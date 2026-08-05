# The shared rank-deficient input fixture, used wherever a test pins the
# `efa_cor_singular` branch.
#
# A composite variable entered alongside its own components is the archetypal way an EFA
# correlation matrix comes out rank deficient, so that is what this builds: nine
# independent variables plus their tenth column, the sum of the first two.
#
# The dimension is not incidental. The singularity check compares
# min(abs(ev)) / max(abs(ev)) against ncol(R) * .Machine$double.eps, and for an exactly
# rank-deficient matrix the smallest computed eigenvalue IS LAPACK's rounding, of order
# eps * ||R||. A 3x3 fixture therefore leaves the ratio within a factor of about 1.5 of
# the threshold -- close enough that a LAPACK build with a slightly larger noise floor
# decides the other way and no error is raised. For this seeded p = 10 fixture the ratio
# is 2.7e-16 against a threshold of 2.2e-15, a factor of 8, so the verdict is the same
# whatever the build's noise floor. The seed keeps the realization fixed on top of that.
#
# set.seed() here is deliberate and its effect is global: helper files are sourced once,
# before any test file, so the ambient stream every unseeded fixture in the suite draws
# from starts from this point. Changing the draw below therefore re-rolls those fixtures.
set.seed(1)
sing_raw <- matrix(stats::rnorm(400 * 9), 400)
sing_raw <- cbind(sing_raw, sing_raw[, 1] + sing_raw[, 2])
sing_cor <- stats::cor(sing_raw)
sing_N <- nrow(sing_raw)

# Invertible but not positive definite: a valid correlation matrix cannot have a
# negative eigenvalue, and this one does, so it exercises the separate non-PD branch.
cor_nposdef <- matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), ncol = 3)
