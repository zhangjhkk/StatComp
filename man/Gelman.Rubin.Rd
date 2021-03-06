% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functionsH.R
\name{Gelman.Rubin}
\alias{Gelman.Rubin}
\title{Calculate the Gelman-Rubin statistic}
\usage{
Gelman.Rubin(psi)
}
\arguments{
\item{psi}{A matrix of scalar summary statistic with n rows and k columns}
}
\value{
A sample that obeys standard Laplace distribution and acceptance rates.
}
\description{
The scalar summary statistic psi_ij is the mean of the ith chain up to time j.
The convergence rate is monitored by calculating Gelman-Rubin statistics of 
multiple data chains at different times.
}
\examples{
\dontrun{
x0<-c(1,5,10)
for(i in 1: 3){
rw[i,] <- rw.Metropolis(1, x0[i], 1000)#generate random walk Metropolis sample
}
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))
}
}
