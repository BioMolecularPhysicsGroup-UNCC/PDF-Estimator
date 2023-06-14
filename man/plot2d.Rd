\name{plot2d}
\alias{plot2d}
\title{Plot two-dimensional probability density estimate}
\description{The plot method for two-dimensional pdfEstimator objects.}
\usage{plot2d(x, xlab = "x", ylab = "y", zlab = "PDF")}
\arguments{
	\item{x}{an "estimatePDFmv" object}	
	\item{xlab}{x-axis label for pdf}
	\item{ylab}{y-axis label for pdf}
	\item{zlab}{z-axis label for pdf}
}
\references{Farmer, J. and D. Jacobs (2018). "High throughput nonparametric probability density estimation." PLoS One 13(5): e0196937.}
\author{Jenny Farmer, Donald Jacobs}
\examples{
library(MultiRNG)
nSamples = 10000
cmat = matrix(c(1.0, 0.0, 0.0, 1.0), nrow = 2, ncol = 2)
meanvec = c(0, 0)
sample = draw.d.variate.normal(no.row = nSamples, d = 2,
                               mean.vec = meanvec, cov.mat = cmat)
mvPDF = estimatePDFmv(sample, resolution = 50)

plot2d(mvPDF)

}
