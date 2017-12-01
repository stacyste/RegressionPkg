#' myLasso
#'
#' This function finds the lasso solution path for various values of the regularization parameter
#' lambda using corrdinate descent. Returns a matrix containing the lasso solution vector
#' beta for each regularization parameter.
#'
#' @param X an n x p matrix of explanatory variables.
#' @param Y n dimensional response vector
#' @param lambda_all vector of regularization parameters
#' @export
#' @return Returns a p+1 x length of lambda_all matrix of regularized betas for each specified level of lambda
#' @examples
#' dim_r=50
#' dim_c=10
#' s=7
#' lam = (100:1)*10
#' x = matrix(rnorm(dim_r*dim_c), nrow = dim_r)
#' beta_true = matrix(rep(0,dim_c), nrow = dim_c)
#' beta_true[1:s] = 1:7
#' y = x%*%beta_true + rnorm(dim_r)
#' beta_true <- c(0, beta_true)
#' lass <- myLasso(x, y, lam)
#' matplot(t(matrix(rep(1, 11), nrow = 1)%*%abs(lass)), t(lass), type = "l")
#'
myLasso <- function(X, Y, lambda_all){

  lambda_all <- sort(lambda_all, decreasing = TRUE)
  X <- cbind(rep(1, nrow(X)), X)

  #initializing variables
  L <- length(lambda_all)
  runs <- 100
  p <- ncol(X)

  #initialize matrices
  db = matrix(rep(0, p), nrow = p)
  beta_all = matrix(rep(0, p*runs), nrow = p)
  beta <- matrix(rep(0,p), nrow = p)
  R <- Y

  SS <- rep(0,p)
  #save ||X||^2
  for(i in 1:p){
    SS[i] <- sum(X[,i]^2)
  }

  #for each lambda starting at a large lambda
  for(l in 1:L){
    lambda <- lambda_all[l]
    for(t in 1:runs){ #for each run of the iteration
      for(j in 1:p){ #for each beta_j
        beta_old <- beta
        db <- sum(R*X[,j])/SS[j]
        b <- beta[j] + db
        if(j != 1){
          b <- sign(b)*max(0, abs(b) - lambda/SS[j])
        }
        db <- b - beta[j]
        R <- R - X[,j]*db #update R
        beta[j] <- b
      }
    }
    beta_all[,l] <- beta
  }
  ## Function should output the matrix beta_all, the
  ## solution to the lasso regression problem for all
  ## the regularization parameters.
  ## beta_all is (p+1) x length(lambda_all)
  return(beta_all)
}


