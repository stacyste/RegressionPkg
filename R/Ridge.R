#' Ridge Regression
#'
#' This function performs ridge regression of Y on X.
#' @param X an n x p matrix of explanatory variables.
#' @param Y an n vector or matrix of dependent variables.
#' @param lambda regularization parameter for lambda >= 0
#' @export
#' @return Returns a vector of betas, the ridge regression solution with p + 1 elements
#' @examples
#' set.seed(12837)
#' X = runif(50)
#' X = matrix(sort(X), nrow = 50)
#' sigma = .1
#' Y =X^2 + rnorm(50)*sigma
#' lambda = 0
#' beta_r <- myRidge(X, Y, lambda)
#' Yhat <- cbind(rep(1, length(X)), X)%*%beta_r
#' plot(X, Y, col = 'indianred')
#' par(new =TRUE
#' plot(X, Yhat, type = 'l', col='slateblue')
#' abline(coef(lm(Y~X)))

myRidge <- function(X, Y, lambda){

  if(class(X) == "matrix"){
    n <- nrow(X)
    p <- ncol(X)
  }
  else{
    n <- length(X)
    p <- 1
  }

  Z <- cbind(rep(1,n), X, Y)
  A <- t(Z)%*%Z
  D <- diag(rep(lambda, p+2))
  D[1,1] <- 0
  D[p+2,p+2] <- 0
  A <- A + D
  S <- mySweep(A, p+1)
  beta_ridge <- S[1:(p+1), p+2]
  ## Function should output the vector beta_ridge, the
  ## solution to the ridge regression problem. beta_ridge
  ## should have p + 1 elements.
  return(beta_ridge)

}
