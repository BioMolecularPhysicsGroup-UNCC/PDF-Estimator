\name{estimatePDF}
\alias{estimatePDF}
\title{Nonparametric Density Estimation}
\description{Estimates the probability density function for a data sample.}
\usage{estimatePDF(sample, pdfLength = NULL, estimationPoints = NULL, 
lowerBound = NULL, upperBound = NULL, target = 70, lagrangeMin = 1, 
lagrangeMax = 200, debug = 0, outlierCutoff = 7, smooth = TRUE)}
\arguments{
	\item{sample}{the data sample from which to calculate the density estimate. If the sample has more than 1 column, the multivariate estimation function, estimatePDFmv(), is called instead.}
	\item{pdfLength}{the desired length of the estimate returned.  Default value is calculated based on sample length.  Overriding this calculation can increase or decrease the resolution of the estimate.}
	\item{estimationPoints}{a vector containing the points to estimate.  If not specified, this is calculated automatically to span the entire sample data.}
	\item{lowerBound}{the lower bound of the PDF, if known.  Default value is calculated based on the range of the data sample.}
	\item{upperBound}{the upper bound of the PDF, if known.  Default value is calculated based on the range of the data sample.}
	\item{target}{a value from 1 to 100 representing the desired confidence percentage for the estimate score.  The default of 70\% represents the most likely score based on empirical simulations.  A lower value may smooth estimates.  A higher value tends to overfit to the sample and is not recommended.}
	\item{lagrangeMin}{minimum number of lagrange multipliers}
	\item{lagrangeMax}{maximum number of lagrange multipliers}
	\item{debug}{verbose output printed to console}
	\item{outlierCutoff}{outliers are automatically detected and removed according to the formula: < Q1 - outlierCutoff * IQR; or > Q3 + outlierCutoff * IQR, where Q1, Q3, and IQR represent the first quartile, third quartile, and inter-quartile range, respectively.  Setting outlierCutoff = 0 turns off outlier detection.}
	\item{smooth}{minimizes noise in estimates, particularly in areas of low data density}
}
\details{
A nonparametric density estimator based on the maximum-entropy method.  Accurately predicts a probability density function (PDF) for random data using a novel iterative scoring function to determine the best fit without overfitting to the sample.
}
\value{
	\item{failedSolution}{returns true if the pdf calculated is not considered an acceptable estimate of the data according to the scoring function.}
	\item{threshold}{represents the quality of the solution returned. Values of 40 to 70  indicate high confidence in the estimate. Values less than 5 are considered to be of poor quality. For more information on scoring see the referenced publication.}
	\item{x}{estimated range of density data}
	\item{pdf}{estimated probability density function}
	\item{cdf}{estimated cummulative density function}
	\item{sqr}{scaled quantile residual.  Provides a sample-size invariant measure of the fluctuations in the estimate.}
	\item{sqrSize}{length of the returned scaled quantile residual.  In most cases, this is the size of the input sample. Exceptions are if outliers are detected and/or if the failedSolution flag is true.}
	\item{lagrange}{values of lagrange multipliers.  Can be used to reproduce the expansions for an analytical solution.}
	\item{r}{inverse of cdf for the sample.}
}
\references{Farmer, J. and D. Jacobs (2018). "High throughput nonparametric probability density estimation." PLoS One 13(5): e0196937.}
\author{Jenny Farmer, Donald Jacobs}
\examples{
#Estimates a normal distribution with 1000 sample points using default parameters

sampleSize = 1000
sample = rnorm(sampleSize, 0, 1)
dist = estimatePDF(sample)

}
