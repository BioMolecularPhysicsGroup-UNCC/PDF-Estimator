plot.PDFe <- function(x, plotPDF = TRUE, plotSQR = FALSE, 
                      plotShading = FALSE, shadeResolution = 100, 
                      showOutlierPercent = 0, outlierColor = "red3",
                      sqrPlotThreshold = 2, sqrColor = "steelblue4", 
                      type="l", lwd = 2, xlab = "x", ylab = "PDF", 
                      legendcex = 0.9, ...){
  
  
  
  if (!is.numeric(lwd)) {
    stop("plot line width must be numeric")
  }
  if (plotPDF == FALSE && plotSQR == FALSE && plotShading == FALSE) {
    stop("no plot type specified: plotPDF, plotSQR, or plotShading")
  }
  
  if (!is.numeric(showOutlierPercent) || 
      showOutlierPercent < 0 || showOutlierPercent >= 100) {
    stop("showOutlierPercent must be greater than 0 and less than 100")
  }
  
  estimate = x
  x = estimate$x
  pdf = estimate$pdf
  sqr = estimate$sqr[1:estimate$sqrSize]
  sample = estimate$sample
  mean = seq(0, 1, 1/(length(sqr) - 1))
  sqrRange = c(-sqrPlotThreshold, sqrPlotThreshold)
  
  if (showOutlierPercent) {
    targets = suppressWarnings(getTarget(length(sqr), showOutlierPercent))
    above = sqr > targets[ , 1]
    below = sqr < targets[ , 2]
    outside = above | below
    iOutside = which(outside == TRUE)
    nOutside = length(iOutside)
  }
  
  if (plotShading) {
    if (plotPDF) {
      y = sort(sample)
      y = y[y >= min(x) & y <= max(x)]
      plotBeta(y, xPlotRange = range(x), resolution = shadeResolution, sqrPlotThreshold = sqrPlotThreshold)
    } else {
      plotBeta(mean, xPlotRange = range(mean), resolution = shadeResolution, sqrPlotThreshold = sqrPlotThreshold)
    }
    par(new = TRUE)
  }
  
  if (plotSQR && plotPDF) {
    xMin = max(min(x), min(sample))
    xMax = min(max(x), max(sample))
    y = sort(sample)
    if (length(sample) > length(sqr)) {
      inRange = which((y > min(x)) & (y < max(x)))
      y = y[inRange]
    }
    if (length(y) == length(sqr)) {
      plot(y, sqr, type = "l", ylim = sqrRange, xlim = range(x), axes = FALSE, col = sqrColor, xlab ="", ylab = "")
      axis(4, ylim = sqrRange, ylab = "SQR", col = sqrColor, col.ticks = sqrColor, col.axis = sqrColor)
      par(new = TRUE)
      if (showOutlierPercent) {
        if (nOutside > 0) {
          symbols(y[iOutside], sqr[iOutside], add = TRUE, circles = rep(0.01, length(iOutside)), bg = outlierColor, fg = outlierColor, inches = 0.01, xlab ="", ylab = "", ylim = sqrRange, xlim = c(0, 1))
          par(new = TRUE)
        }
      }
    }
  }
  
  if (plotSQR && !plotPDF) {
    if (missing(xlab)) {
      if (missing(ylab)) {
        plot(mean, sqr, xlab = "order fraction", ylab = "SQR", ylim = sqrRange, type = type, lwd = lwd, col = sqrColor, ...)
      } else {
        plot(mean, sqr, xlab = "order fraction", ylim = sqrRange, type = type, lwd = lwd, col = sqrColor, ...)
      }
    } else {
      if (missing(ylab)) {
        plot(mean, sqr, ylab = "SQR", ylim = sqrRange, type = type, lwd = lwd, col = sqrColor, ...)
      } else {
        plot(mean, sqr, ylim = sqrRange, type = type, lwd = lwd, col = sqrColor, ...)
      }
    }
    par(new = TRUE)
    if (showOutlierPercent) {
      plot(mean, targets[ , 1], type = "l", lty = 2, ylim = sqrRange, xlab = "", ylab = "")
      par(new = TRUE)
      plot(mean, targets[ , 2], type = "l", lty = 2, ylim = sqrRange, xlab = "", ylab = "")
      par(new = TRUE)
      percentActual = ceiling(nOutside / length(sqr) * 100)
      legend("bottomright", inset = 0.01, horiz = TRUE, legend = paste("outliers ", percentActual, "%", sep = ""), col = outlierColor, lty=1, cex=legendcex, box.lty = 0)
      legend("topright", inset = 0.01, horiz = TRUE, legend = paste("confidence ", showOutlierPercent, "%", sep = ""), col = sqrColor, lty=2, cex=legendcex, box.lty = 0)
      if (nOutside > 0) {
        symbols(mean[iOutside], sqr[iOutside], add = TRUE, circles = rep(0.01, length(iOutside)), bg = outlierColor, fg = outlierColor, inches = 0.01, xlab ="", ylab = "", ylim = sqrRange, xlim = c(0, 1))
        par(new = TRUE)
      }
    }
  }
  
  if (plotPDF) {
    plot.default(x, pdf, type = type, lwd = lwd, xlab = xlab, ylab = ylab, ...)
    if (!plotSQR && showOutlierPercent) {
      if (nOutside > 0) {
        xMin = max(min(x), min(sample))
        xMax = min(max(x), max(sample))
        d = seq(xMin, xMax, length.out = length(sqr))
        outliers = approx(x, pdf, xout = d[iOutside])
        symbols(outliers$x, outliers$y, add = TRUE, circles = rep(0.01, length(iOutside)), bg = outlierColor, fg = outlierColor, inches = 0.01, xlab ="", ylab = "", ylim = sqrRange, xlim = c(0, 1))
      }
    }
  }
  
  par(new = FALSE)
  
}