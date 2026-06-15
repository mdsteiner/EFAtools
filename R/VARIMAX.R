# SPSS varimax implementation as described in the SPSS manual: ftp://public.dhe.ibm.com/software/analytics/spss/documentation/statistics/23.0/en/client/Manuals/IBM_SPSS_Statistics_Algorithms.pdf

.VARIMAX_SPSS <- function(x, normalize = TRUE, precision = 1e-5) {

  # get dimensions
  n_cols <- ncol(x)
  n_rows <- nrow(x)

  if (isTRUE(normalize)) {
    # Kaiser row normalisation. Floor the weights so a zero-communality row
    # (h2 = 0) does not divide by zero and inject NaN into the rotation; such a
    # row stays zero throughout and is restored unchanged below. The length-nrow
    # w recycles column-major, scaling row i by 1/w[i] (mirrors PROCRUSTES.R).
    h2 <- diag(x %*% t(x))
    w <- sqrt(h2)
    w[!is.finite(w) | w < 1e-15] <- 1
    x <- x / w
  }


  # initialize rotation matrix as identity matrix
  rotmat <- diag(1, n_cols)

  SV_now <- .SV(x)


  for (it in 1:1000) {

    for (col_j in 1:(n_cols - 1)) {
      for (col_k in (col_j + 1):n_cols) {

        u_p <- x[,col_j] ** 2 - x[,col_k] ** 2
        v_p <- 2 * x[,col_j] * x[,col_k]

        A <- sum(u_p)
        B <- sum(v_p)
        C <- sum(u_p ** 2 - v_p ** 2)
        D <- sum(2 * u_p * v_p)

        X <- D - 2 * A * B / n_rows
        Y <- C - (A ** 2 - B ** 2) / n_rows

        P <- 1/4 * atan2(X, Y)

        if (abs(sin(P)) > 1e-15) {
          sub_rot <- matrix(c(cos(P), -sin(P),
                              sin(P), cos(P)),
                            ncol = 2, byrow = TRUE)

          x[, c(col_j, col_k)] <- x[, c(col_j, col_k)] %*% sub_rot
          rotmat[, c(col_j, col_k)] <- rotmat[, c(col_j, col_k)] %*% sub_rot
        }

      }
    }

    SV_old <- SV_now
    SV_now <- .SV(x)

    if (abs(SV_now - SV_old) <= precision) {
      break
    }

  }

  if (isTRUE(normalize)) {
    x <- x * w
  }

  # Create output list
  list(rotmat = rotmat, loadings = x, iter = it)

}
