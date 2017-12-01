#' Sweep Function
#'
#' This function performs sweep operations on a matrix in order to give the negative inverse
#' @param A a square matrix
#' @param m the number of pivot elements in A
#' @param lambda regularization parameter for lambda >= 0
#' @export
#' @return Returns a swept matrix
#' @examples
#' A <- rbind(c(-13, 14, -2), c(0, 2, 2), c(3, 9, -7))
#' mySweep(A, 3)
#'
mySweep <- function(A, m){

  n <- nrow(A)

  for(k in 1:m){
    for(i in 1:n)
      for(j in 1:n)
        if(i != k  & j != k)
          A[i,j] <- A[i,j] - A[i,k]*A[k,j]/A[k,k]

        for(i in 1:n)
          if(i != k)
            A[i,k] <- A[i,k]/A[k,k]

          for(j in 1:n)
            if(j != k)
              A[k,j] <- A[k,j]/A[k,k]
            A[k,k] <- - 1/A[k,k]
  }
  return(A)
}
