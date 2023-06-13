summary.PDFe <- function(object, ...){
  estimate = object
  
  lowBound = Inf
  highBound = Inf
  if (estimate$isLow) {
    lowBound = estimate$low
  }
  if (estimate$isHigh) {
    highBound = estimate$high
  }
  
  if (estimate$outlierCutoff > 0) {
    eliminateOutlier = "YES"
  } else {
    eliminateOutlier = "NO"
  }
  
  print(sprintf("PARAMETERS"), quote = FALSE)
  print(sprintf("   specified boundaries:         [%0.2f  %0.2f]", 
                lowBound, highBound), quote = FALSE)
  print(sprintf("   confidence threshold target:  %0.2f%%", 
                estimate$target), quote = FALSE)
  print(sprintf("   detect outliers:              %s", 
                eliminateOutlier), quote = FALSE)
  print("", quote = FALSE)
  
  print(sprintf("INPUT SAMPLE"), quote = FALSE)
  print(sprintf("   original sample size:        %d", 
                estimate$inputLength), quote = FALSE)
  print(sprintf("   sample range:                [%0.4f  %0.4f]", 
                min(estimate$sample), max(estimate$sample)), quote = FALSE)
  print(sprintf("   calculated range:            [%0.4f  %0.4f]", 
                min(estimate$x), max(estimate$x)), quote = FALSE)
  print(sprintf("   number of outliers removed:  %d", 
                estimate$inputLength - estimate$sqrSize), quote = FALSE)
  print("", quote = FALSE)
  
  print(sprintf("DENSITY ESTIMATE"), quote = FALSE)
  print(sprintf("   minimum density:                x = %0.4f;  pdf = %0.4f", 
                estimate$x[which.min(estimate$pdf)], min(estimate$pdf)), quote = FALSE)
  print(sprintf("   maximum density:                x = %0.4f;  pdf = %0.4f", 
                estimate$x[which.max(estimate$pdf)], max(estimate$pdf)), quote = FALSE)
  print(sprintf("   confidence threshold achieved:  %0.2f%%", 
                estimate$threshold), quote = FALSE)
  print(sprintf("   number of Lagrange multipliers: %d", 
                length(which(estimate$lagrange != 0))), quote = FALSE)
  
  return(invisible(object))
}


summary.mvPDFe <- function(object, ...){
  estimate = object 
  
  print(sprintf("Multivariate Probability Density"), quote = FALSE)
  print("", quote = FALSE)
  print(sprintf("   number of variables:  %d", 
                estimate$variableLength), quote = FALSE)
  print(sprintf("   sample length:        %d", 
                estimate$sampleLength), quote = FALSE)
  print(sprintf("   grid resolution:      %d", 
                estimate$resolution), quote = FALSE)
  
  print(sprintf("   minimum density:      %0.4f", 
                min(estimate$pdf)), quote = FALSE)
  print(sprintf("   maximum density:      %0.4f", 
                max(estimate$pdf)), quote = FALSE)
  print("", quote = FALSE)
  
  print(sprintf("Estimated Ranges"), quote = FALSE)
  print("", quote = FALSE)
  i = 1;
  for (v in 1:estimate$variableLength) {
    temp = estimate$x[i:(i + estimate$resolution - 1)];
    print(sprintf("  variable %d:  [%0.4f  %0.4f]", v, min(temp), max(temp)), quote = FALSE)
    i = i + estimate$resolution
  }
  
  return(invisible(object))
}
  