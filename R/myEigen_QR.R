#' PCA using QR Decomposition
#'
#' This function performs PCA using QR decomposition.
#' @param A A square, symmetric matrix
#' @param numIter the number of iterations to run QR decomposition on, default is 1000
#' @export
#' @return Returns a list of eigenvalues and eigenvectors of A
#' @examples
#' set.seed(239)
#' n <- 100
#' p <- 3
#' X <- matrix(rnorm(n * p), nrow = 100)
#' A <- t(X) %*% X
#'
myEigen_QR <- function(A, numIter = 1000){
  r <- nrow(A)
  c <- ncol(A)
  V <- matrix(runif(r * r), r, r)

  for(i in 1:numIter){
    Q <- myQR(V)$Q
    V <- A %*% Q
  }

  eigen_out <- myQR(V)
  Q <- eigen_out$Q
  R <- eigen_out$R
  return(list("D" = diag(R), "V" = Q))
}
