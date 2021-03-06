% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functionsH.R
\name{loglike_obs}
\alias{loglike_obs}
\title{Observed data likelihood}
\usage{
loglike_obs(theta)
}
\arguments{
\item{theta}{a Gene frequency vector (p,q)}
}
\value{
Log of Likelihood of frequencies under observed data
}
\description{
As for the A-B-O blood type problem, observed data include nA, nB and nOO. 
Under these given conditions, the likelihood function of frequency p and q (and r=1-p-q) 
is shown.
}
\examples{
\dontrun{
loglike_obs(c(0.3,0.3))
}
}
