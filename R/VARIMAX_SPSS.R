#' Varimax implementation based on the SPSS algorithms
#'
#' This is an implementation exactly following the SPSS algorithms manual.
#'
#' @param loadings matrix. Unrotated loadings, e.g. from a PAF.
#' @param precision numeric. The criterion after which to end the iterative
#'  procedure if the difference from one svd to the other is smaller than the
#'  criterion. Default is 1e-5 (same as SPSS).
#' @param normalize logical. Whether Kaiser normalization should be performed.
#'  Default is TRUE.
#'
#' @return A list of two.
#' \item{rot_loadings}{The varimax rotated loadings (the pattern matrix).}
#' \item{rotmat}{The rotation matrix.}
#' @export
#' @source IBM Corp. (2014). IBM SPSS Statistics 23 algorithms [Computer software manual]. Armonk, NY
VARIMAX_SPSS <- function(loadings, precision = 1e-5, normalize = TRUE) {

  # get dimensions
  n_cols <- ncol(loadings)
  n_rows <- nrow(loadings)
  dim_names <- dimnames(loadings)

  if (n_cols == 1) {
    return(loadings)
  }

  if (normalize) {
    h2 <- diag(loadings %*% t(loadings))
    loadings <- loadings/sqrt(h2)
  }


  # initialize rotation matrix as identity matrix
  rotmat <- diag(1, n_cols)

  SV_now <- SV(loadings)


  for (it in 1:10000) {

    for (col_j in 1:(n_cols - 1)) {
      for (col_k in (col_j + 1):n_cols) {

        u_p <- loadings[,col_j] ** 2 - loadings[,col_k] ** 2
        v_p <- 2 * loadings[,col_j] * loadings[,col_k]

        A <- sum(u_p)
        B <- sum(v_p)
        C <- sum(u_p ** 2 - v_p ** 2)
        D <- sum(2 * u_p * v_p)

        X <- D - 2 * A * B / n_rows
        Y <- C - (A ** 2 - B ** 2) / n_rows

        P <- 1/4 *atan2(X, Y)

        sub_rot <- matrix(c(cos(P), -sin(P),
                            sin(P), cos(P)),
                          ncol = 2, byrow = TRUE)

        loadings[, c(col_j, col_k)] <- loadings[, c(col_j, col_k)] %*% sub_rot
        rotmat[, c(col_j, col_k)] <- rotmat[, c(col_j, col_k)] %*% sub_rot


      }
    }

    SV_old <- SV_now
    SV_now <- SV(loadings)

    if (abs(SV_now - SV_old) <= precision) {
      break
    }

  }

  if (it >= 10000) {
    warning("Maximum number of iterations (10'000) reached without convergence.")
  }

  if (normalize) {
    loadings <- loadings * sqrt(h2)
  }

  # Create output list
  list(rotmat = rotmat, loadings = loadings)

}
