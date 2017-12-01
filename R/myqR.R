#' QR Decomposition
#'
#' This function performs QR decomposition on a matrix  A into a product A = QR of an orthogonal matrix Q and an upper triangular matrix R
#' @param A an n by m  matrix whose QR decomposition is to be computed
#' @export
#' @return Q, an orthogonal matrix Q  and an upper triangular matrix R
#' @examples
#' A <- rbind(c(-13, 14, -2), c(0, 2, 2), c(3, 9, -7))
#' myQR(A)
myQR <- function(A){
  n <- nrow(A)
  m <- ncol(A)
  Q <- diag(n)
  R <- A

  for(k in 1:(m - 1)){
    x      <- rep(0, n)
    x[k:n] <- R[k:n, k]
    s      <- -1 * sign(x[k])
    v      <- x
    v[k]   <- x[k] - s * norm(x, type = "2")
    u      <- v / norm(v, type = "2")

    R <- R - 2 * u %*% t(u) %*% R
    Q <- Q - 2 * u %*% t(u) %*% Q

  }
  return(list("Q" = t(Q), "R" = R))

}
