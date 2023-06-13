lines.PDFe <- function(x, showOutlierPercent = 0, 
                      outlierColor = "red3", lwd = 2, ...){
  
  if (!is.numeric(showOutlierPercent) || 
      showOutlierPercent < 0 || showOutlierPercent >= 100) {
    stop("showOutlierPercent must be greater than 0 and less than 100")
  }
  
  
  estimate = x
  x = estimate$x
  pdf = estimate$pdf
  sqr = estimate$sqr[1:estimate$sqrSize]
  sample = estimate$sample
  
  lines.default(x, pdf, lwd = lwd, ...)
  
  if (showOutlierPercent) {
    targets = suppressWarnings(getTarget(length(sqr), showOutlierPercent))
    above = sqr > targets[ , 1]
    below = sqr < targets[ , 2]
    outside = above | below
    iOutside = which(outside == TRUE)
    if (length(iOutside) > 0) {
      xMin = max(min(x), min(sample))
      xMax = min(max(x), max(sample))
      d = seq(xMin, xMax, length.out = length(sqr))
      outliers = approx(x, pdf, xout = d[iOutside])
      symbols(outliers$x, outliers$y, add = TRUE, circles = rep(0.01, length(iOutside)), bg = outlierColor, fg = outlierColor, inches = 0.01, xlab ="", ylab = "", ylim = c(-2, 2), xlim = c(0, 1))
    }
  }
}