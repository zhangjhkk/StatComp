% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functionsH.R
\name{ranfe}
\alias{ranfe}
\title{Generate sample from the special distribution}
\usage{
ranfe(n)
}
\arguments{
\item{n}{the sample size}
}
\value{
a random sample of size \code{n} from fe
}
\description{
The pdf of distribution if fe=3/4*(1-x^2),|x|<=1
Devroye and Gyorfi give the following algorithm for simulation from this distribution.
}
\examples{
\dontrun{
rangfe(10)
}
}
