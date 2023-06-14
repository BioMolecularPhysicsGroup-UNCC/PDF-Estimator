\name{plot3d}
\alias{plot3d}
\title{Plot three-dimensional probability density estimate}
\description{The plot method for three-dimensional pdfEstimator objects.  Plots two-dimensional cross-sectional slices.}
\usage{plot3d(x, xs = c(0), ys = c(0), zs = NULL, xlab = "X1", ylab = "X2", zlab = "X3")}
\arguments{
	\item{x}{an "estimatePDFmv" object}	
	\item{xlab}{x-axis label for pdf}
	\item{ylab}{y-axis label for pdf}
	\item{zlab}{z-axis label for pdf}
	\item{xs, ys, zs}{Vectors or matrices. Vectors specify the positions in x, y or z where the slices (planes) are to be drawn.}
}
\references{Farmer, J. and D. Jacobs (2018). "High throughput nonparametric probability density estimation." PLoS One 13(5): e0196937.}
\author{Jenny Farmer, Donald Jacobs}

