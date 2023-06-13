\name{PDFEstimator-package}
\alias{PDFEstimator-package}
\alias{PDFEstimator}
\title{Nonparametric Probability Density Estimation and Analysis}
\description{
   This package provides tools for nonparametric density estimation according to the maximum entropy method described in Farmer and Jacobs (2018).  PDFEstimator includes functionality for creating a robust data-driven estimate from a data sample requiring minimal user intervention, thus suitable for high-throughput applications.  
   
   Additionally, the package includes advanced plotting and visual diagnostics for confidence thresholding and identification of potentially poorly fitted regions of the estimate.  These diagnostics are made available to other density estimation methods through a custom conversion utility, allowing for equitable comparison between estimates.
}
\details{
\tabular{ll}{
Main function for estimating the density from a data sample:  \tab \code{\link{estimatePDF}}\cr\cr 
Customized plotting function for visual inspection and analysis:  \tab \code{\link{plot}}\cr\cr
Plotting function for densities with 2 variables:  \tab \code{\link{plot2d}}\cr\cr
Plotting function for densities with 3 variables:  \tab \code{\link{plot3d}}\cr\cr
Conversion utility for estimates obtained by other methods:  \tab \code{\link{convertToPDFe}}\cr\cr
Calculation of boundaries for user-defined confidence levels:  \tab \code{\link{getTarget}}\cr\cr
Optional background shading outlining expected variance by position:  \tab \code{\link{plotBeta}}\cr\cr
Utility for additional point approximation for an existing estimate:  \tab \code{\link{approximatePoints}}\cr\cr
}
}
\references{
Farmer, J. and D. Jacobs (2018). "High throughput nonparametric probability density estimation." PLoS One 13(5): e0196937. \doi{10.1371/journal.pone.0196937}.
}
\author{
Jenny Farmer, University of North Carolina at Charlotte.  \email{jfarmer@carolina.rr.com}.

Donald Jacobs, University of North Carolina at Charlotte.  \email{djacobs1@uncc.edu}.

}

