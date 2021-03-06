% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sweep.R
\name{mySweep}
\alias{mySweep}
\title{Sweep Function}
\usage{
mySweep(A, m)
}
\arguments{
\item{A}{a square matrix}

\item{m}{the number of pivot elements in A}

\item{lambda}{regularization parameter for lambda >= 0}
}
\value{
Returns a swept matrix
}
\description{
This function performs sweep operations on a matrix in order to give the negative inverse
}
\examples{
A <- rbind(c(-13, 14, -2), c(0, 2, 2), c(3, 9, -7))
mySweep(A, 3)

}
