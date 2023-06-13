print.PDFe <- function(x, ...){
  estimate = x
  
  if (estimate$isSpecifyPoints) {
    print(sprintf("x        pdf"), quote = FALSE)
    print("", quote = FALSE)
    for (i in 1:estimate$estimationLength) {
      print(sprintf("%0.4f   %0.4f", 
          estimate$estimationPoints[i], estimate$pointsPDF[i]), quote = FALSE)
    }
  } else {
    print(sprintf("x        pdf       cdf"), quote = FALSE)
    print("", quote = FALSE)
    for (i in 1:estimate$outputLength) {
      print(sprintf("%0.4f   %0.4f    %0.4f", 
          estimate$x[i], estimate$pdf[i], estimate$cdf[i]), quote = FALSE)
    }
  }
  return(invisible(x))
}