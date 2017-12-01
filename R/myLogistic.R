#' Logistic Regression Fitting
#'
#' This function performs the logistic regression of Y on X and returns the coeffients. No intercept term is assumed
#' @param X an n x p matrix of explanatory variables
#' @param Y an n dimensional vector of binary responses
#' @export
#' @return A list of the logistic betas and their respective standard errors
#' @examples
#' set.seed(18274)
#' n <- 500
#' p <- 4
#' X    <- matrix(rnorm(n * p), nrow = n)
#' beta <- c(12, -2,-3, 4)
#' Y    <- 1 * (runif(n) <  1 / (1 + exp(-(X %*% beta))))
#' logistic_beta <- myLogistic(X, Y)
#' logistic_beta
#'
myLogistic <- function(X, Y){

  expit <- function(x){
    1 / (1 + exp(-x))
  }

  n = nrow(X)
  p = ncol(X)

  beta = rep(0,p)
  pr = rep(0,n)
  epsilon = 10^(-6)
  error = 10

  while(error > epsilon){
    eta = X%*%beta
    pr = expit(eta)
    weight = pr*(1-pr)
    z = eta + (Y - pr)/weight

    x_tilde = rep(sqrt(weight), p)*X
    y_tilde = sqrt(weight)*z

    lm = myLM(x_tilde, y_tilde, intercept = FALSE)
    beta_new = lm$beta
    error = sum(abs(beta_new-beta))
    beta = beta_new
  }
  se = sqrt(diag(solve(t(x_tilde)%*%x_tilde)))

  ## Function returns the logistic regression solution vector
  return(data.frame("beta" = beta,"error" = se ))
}
