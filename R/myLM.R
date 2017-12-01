#' Linear Model Function
#'
#' This function calculates the regression coefficients for linear regression as well as standard errors
#' @param X an n x p matrix of explanatory variables
#' @param Y an n dimensional vector of responses
#' @param intercept Boolean for whether an intercept term should be included, default is true
#' @export
#' @return Returns the OLS betas along with their respective standard errors
#' @examples
#' set.seed(12083)
#' X <- matrix(rnorm(100), 50, 2)
#' Y <- matrix(3 * X[,1] + 2 * X[,2] + rnorm(50))
#' myLM(X, Y)
#'
myLM <- function(X, Y, intercept = TRUE){

  n <- nrow(X)
  p <- ncol(X)

  if(intercept){
    Z <- cbind(rep(1,n), X, Y)
    R <- myQR(Z)$R

    R1 <- R[1:(p + 1), 1:(p + 1)]
    Y1 <- R[1:(p + 1), p + 2]
    Y2 <- R[(p+2):n, p+2]
  }
  else{
    Z <- cbind(X, Y)
    R <- myQR(Z)$R

    R1 <- R[1:p, 1:p]
    Y1 <- R[1:p, p + 1]
    Y2 <- R[(p+1):n, p+1]
  }

  beta_ls <- solve(R1) %*% Y1
  #calculating the standard errors
  sigmasq <-  sum(Y2^2)/(n-p-1)
  se <- sqrt(diag(solve(t(R1)%*%R1))*sigmasq)
  return(data.frame("beta" = beta_ls, "error" = se))

}
