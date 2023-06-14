approximatePoints <- function(estimate, estimationPoints) {
 
  if (is(estimate) != "PDFe") {
    stop("estimate must be a PDFe class")
  }
  if (!is.numeric(estimationPoints)) {
    stop("estimationPoints must numeric")
  }
  
  approximatePDF = approx(x = estimate$x, y = estimate$pdf, 
                          xout = estimationPoints, yleft = 0, yright = 0)
  estimate$estimationPoints = approximatePDF$x
  estimate$pointsPDF = approximatePDF$y
  estimate$isSpecifyPoints = TRUE
  estimate$estimationLength = length(approximatePDF$x)
  print(estimate)
  return(invisible(approximatePDF))
  
}