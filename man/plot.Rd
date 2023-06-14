\name{plot.PDFe}
\alias{plot.PDFe}
\alias{plot}
\title{Plot Method for Nonparametric Density Estimation}
\description{The plot method for pdfEstimator objects.}
\usage{
  \method{plot}{PDFe}(x, plotPDF = TRUE, plotSQR = FALSE,
  plotShading = FALSE, shadeResolution = 100, 
  showOutlierPercent = 0, outlierColor = "red3", sqrPlotThreshold = 2, 
  sqrColor = "steelblue4", type="l", lwd = 2, xlab = "x", ylab = "PDF", 
                      legendcex = 0.9, ...)}
\arguments{
	\item{x}{an "estimatePDF" object}	
	\item{plotPDF}{plot the probability density function}
	\item{plotSQR}{plot the scaled quantile residual of the estimate}
	\item{plotShading}{plot a gray background shading representing the probability density of the scaled quantile residuals}
	\item{shadeResolution}{the number of sample points plotted in the background if plotShading = TRUE. Increasing resolution will provide sharper contours and take longer to plot.}
	\item{showOutlierPercent}{specify confidence threshold for outliers}
	\item{outlierColor}{color for outliers positions outside of threshold defined in showOutlierPercent}
  \item{sqrPlotThreshold}{magnitude of ylim above and below zero for SQR plot}
	\item{sqrColor}{color for sqr plot for positions within the threshold defined in showOutlierPercentage}
	\item{type}{plot type for pdf. If plotPDF = FALSE and plotSQR = TRUE, then the sqr plot uses this type}
	\item{lwd}{line width for pdf. If plotPDF = FALSE and plotSQR = TRUE, then the sqr plot uses this line width}
	\item{xlab}{x-axis label for pdf. If plotPDF = FALSE and plotSQR = TRUE, then the sqr plot uses this label}
	\item{ylab}{y-axis label for pdf. If plotPDF = FALSE and plotSQR = TRUE, then the sqr plot uses this label}
	\item{legendcex}{expansion factor for legend point size with sqr plot type, for plotPDF = FALSE and plotSQR = TRUE}
	\item{...}{further plotting parameters}
}
\references{Farmer, J. and D. Jacobs (2018). "High throughput nonparametric probability density estimation." PLoS One 13(5): e0196937.}
\author{Jenny Farmer, Donald Jacobs}
\examples{

plot(estimatePDF(rnorm(1000, 0, 1)), plotSQR = TRUE, showOutlierPercent = 99)

}
