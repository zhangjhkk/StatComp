% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functionsH.R
\name{Tn}
\alias{Tn}
\title{Nearest Neighbor Statistic for Equal Distribution Test}
\usage{
Tn(z, ix, sizes, k)
}
\arguments{
\item{z}{A data set composed of two multidimensional samples}

\item{ix}{index used during bootstrap}

\item{sizes}{A vector of two sample sizes}

\item{k}{the number of chosen neighbors}
}
\value{
nearest neighbor test statistic
}
\description{
Find the nearest k neighbors by function "nn2".According to 
whether each sample point in the sample 1 and its corresponding k neighbors
are from the same sample, assign 1 or 0 to points of sample 1 and take the 
average as a statistic.
}
\examples{
\dontrun{
  n1<-10;n2<-20;n<-n1+n2
  sizes<-c(n1,n2);k=3
  x<-rnorm(10,0,1)
  y<-rnorm(20,0,2)
  z <- c(x, y)
  ix<-1:n
  tn<-Tn(z,ix,sizes,k)
}
}
