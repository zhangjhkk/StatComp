% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{metro}
\alias{metro}
\title{Random Walk Metropolis Sampler}
\usage{
metro(n, sigma, x0)
}
\arguments{
\item{n}{the size of sample}

\item{sigma}{the standard deviation of proposal Gaussian distribution}

\item{x0}{the initial value of random walk}
}
\value{
A sample that obeys standard Laplace distribution and acceptance rates.
}
\description{
Generate  a sample that obeys standard Laplace distribution by Random walk Metropolis, which uses proposal distribution normal (X_t,sigma^2). The data generation process records the acceptance rate.
}
\examples{
\dontrun{
a1<-as.numeric(metro(4000,1,10))
a1[1]
}
}
