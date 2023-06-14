\name{estimatePDFmv}
\alias{estimatePDFmv}
\title{Multivariate Nonparametric Density Estimation}
\description{Estimates the multivariate probability density function for a data sample containing up to 3 variables.}
\usage{estimatePDFmv(sample, debug = 0, resolution = NULL)}
\arguments{
	\item{sample}{data sample from which to calculate the density estimate. Each column of data represents an independent variable.}
	\item{debug}{verbose output printed to console}
	\item{resolution}{grid length of data points for each independent variable.}
}
\details{
A multivariate nonparametric density estimator based on the maximum-entropy method.  Accurately predicts a probability density function (PDF) for random data for 1, 2, or 3 variables.
}
\value{
	\item{x}{estimated range of density data}
	\item{pdf}{estimated probability density function}
}
\references{Farmer, J. and D. Jacobs (2018). "High throughput nonparametric probability density estimation." PLoS One 13(5): e0196937.}
\author{Jenny Farmer, Donald Jacobs}
\examples{
#Estimates a 2-variable normal distribution with 10000 sample points

library(MultiRNG)
nSamples = 5000
cmat = matrix(c(1.0, 0.0, 0.0, 1.0), nrow = 2, ncol = 2)
meanvec = c(0, 0)
sample = draw.d.variate.normal(no.row = nSamples, d = 2,
                               mean.vec = meanvec, cov.mat = cmat)
mvPDF = estimatePDFmv(sample)


}
