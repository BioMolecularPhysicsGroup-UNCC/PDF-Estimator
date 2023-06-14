\name{getTarget}
\alias{getTarget}
\title{Define Target Outliers}
\description{calculates position-dependent threshold values about the mean according to a beta distribution with parameters k and (n + 1 - k), where k is the position and n is the total number of positions.  These beta distributions represent probability per position for sort order statistics for a uniform distribution.  This function returns a two-column matrix defining the 
upper and lower variances of the scaled quantile residual for the target threshold}
\usage{getTarget(Ns, target)}
\arguments{
  \item{Ns}{number of samples}
  \item{target}{target confidence threshold}
}
\details{plotTarget is intended for use with plot.PDFe density estimation objects for plotting scaled quantile residuals, but can be called as a stand-alone user method as well.
}
\value{
\item{bounds}{a two dimensional matrix defining the  upper and lower 
variance boundaries for the requested target.}
}
\references{Farmer, J. and D. Jacobs (2018). "High throughput nonparametric probability density estimation." PLoS One 13(5): e0196937.}
\author{Jenny Farmer, Donald Jacobs}
\seealso{plot.PDFe}
\examples{
#returns boundaries of position-dependent variance calculated for 100 data samples
#  for a threshold of 40\%
getTarget(100, 40)
}
