\name{approximatePoints}
\alias{approximatePoints}
\title{Approximate Data Points}
\description{Returns additional point estimates based on an existing estimate.}
\usage{approximatePoints(estimate, estimationPoints)}
\arguments{
\item{estimate}{the pdfe object returned from estimatePDF or convertToPDFe}
\item{estimationPoints}{a vector of additional points to estimate.}
}

\details{
This method approximates density estimates for the points specified by performing a linear interpolation on an existing probability density function. For a more precise point estimation, call estimatePDF with the estimationPoints argument.
}

\references{Farmer, J. and D. Jacobs (2018). "High throughput nonparametric probability density estimation." PLoS One 13(5): e0196937.}
\author{Jenny Farmer, Donald Jacobs}
\examples{
#Estimates a normal distribution with 1000 sample points using default
# parameters, then prints approximate probability density at points -3, 0, and 1
  
sampleSize = 1000
sample = rnorm(sampleSize, 0, 1)
dist = estimatePDF(sample)
approximatePoints(dist, c(-3, 0, 1))
}
