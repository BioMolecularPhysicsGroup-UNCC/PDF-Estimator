\name{lines.PDFe}
\alias{lines.PDFe}
\alias{lines}
\title{Plot Lines Method for Nonparametric Density Estimation}
\description{The lines method for pdfEstimator objects.}
\usage{
  \method{lines}{PDFe}(x, showOutlierPercent = 0, outlierColor = "red3",
  lwd = 2, ...)}
\arguments{
	\item{x}{an "estimatePDF" object}	
	\item{showOutlierPercent}{specify confidence threshold for outliers}
	\item{outlierColor}{color for outliers positions outside of threshold defined in showOutlierPercent}
	\item{lwd}{line width for pdf. If plotPDF = FALSE and plotSQR = TRUE, then the sqr plot uses this line width}
	\item{...}{further plotting parameters}
}
\references{Farmer, J. and D. Jacobs (2018). "High throughput nonparametric probability density estimation." PLoS One 13(5): e0196937.}
\author{Jenny Farmer, Donald Jacobs}
\examples{
plot(estimatePDF(rnorm(1000, 0, 1)))
lines(estimatePDF(rnorm(1000, 0, 1)), col = "gray")
}
